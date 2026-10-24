########testing soundscape protocol#######
library(data.table) #load raw data For loading large CSV files
library(lunar) #install.packages("lunar"), extract moon phase
library(dplyr) #data manipulation 
library(tidyverse) # For data manipulation

library(purrr)
library(readr)
library(openmeteo) # collect weather data for case 2 # install.packages('openmeteo') only available for certain versions of R
library(unmarked)  # For occupancy modeling

library(AICcmodavg) # For model comparison
library(lubridate) # For handling dates (if needed)


library(scales) # label y axis in percentage
library(ggpubr) #combine figures

##############Multi-Season Occupancy Model (case 1 for NC NABat data)################

dfo <- fread("occu3.csv")
dfo$lunar <- lunar.illumination(as.Date(dfo$date, format = "%Y%m%d")) # extract moon phase
dfo$year <-as.numeric(dfo$year)

#reorder categorical variables
dfo$structure <- factor(dfo$structure, levels = c("interior", "corridor", "edge","field"))
dfo$canopy <- factor(dfo$canopy, levels = c("open", "half", "closed"))
dfo$water <- factor(dfo$water, levels = c("yes", "no"))

occu_data <- dfo

##################epfu####################
# Create detection history and turn into detection matrix#
detection_history <- occu_data %>%
  # Select relevant columns (site, year, visit, species)
  select(site, cell, year, night, EPTFUS) 
#replace EPTFUS in species
# epfu EPTFUS, labo LASBOR, laci LASCIN, lano LASNOC, mylu MYOLUC, #
# nyhu NYCHUM, pesu PERSUB, tabr TADBRA#
# mygr MYOGRI, myau MYOAUS, myse MYOSEP

# Convert to detection matrix - SAFE METHOD, not sure if efficient
y_matrix <- detection_history %>%
  pivot_wider(names_from = night, values_from = EPTFUS,
              values_fill = NA, names_sort = TRUE) %>%
  arrange(site, cell, year)

# Create 3D array with proper dimension checking
n_sites <- n_distinct(y_matrix$site)
n_seasons <- n_distinct(y_matrix$year)
n_visits <- 4  # Fixed number of visits

y_array <- array(NA, dim = c(n_sites, n_visits, n_seasons))

for(i in seq_len(n_sites)){
  current_site <- unique(y_matrix$site)[i]
  
  for(t in seq_len(n_seasons)){
    current_year <- unique(y_matrix$year)[t]
    
    site_year_data <- y_matrix %>%
      filter(site == current_site, year == current_year)
    
    if(nrow(site_year_data) == 1){
      # Safely extract visit data
      visit_data <- as.numeric(site_year_data[1, paste0("night", 1:4)])
      
      # Verify length matches expected visits
      if(length(visit_data) == n_visits){
        y_array[i, , t] <- visit_data
      } else {
        warning(paste("Unexpected number of visits for site", current_site, 
                      "year", current_year))
      }
    }
    # Otherwise leaves as NA (for missing surveys)
  }
}

# Verify the array structure
str(y_array)
table(apply(y_array, 3, function(x) sum(!is.na(x)))) # Surveys per year

y_matrix1 <- matrix(y_array, 
                    nrow = dim(y_array)[1], 
                    ncol = dim(y_array)[2] * dim(y_array)[3])


# Site-level covariates (constant across years)
site_covs <- occu_data %>%
  distinct(site, structure, canopy, water, mean_2019_250, sd_2019_250, region, wns_sus) %>%
  arrange(site)

# Observation or year level covariates (change across years and nights)
yearly_covs <- occu_data[,c(1,2,32:34,49,50, 51)]

# Prepare yearlySiteCovs (must be named lists of matrices)
yearly_covs_list <- list(
  tmax = matrix(yearly_covs$tmax, nrow = n_sites, ncol = n_seasons),
  tmin = matrix(yearly_covs$tmin, nrow = n_sites, ncol = n_seasons),
  prcp = matrix(yearly_covs$prcp, nrow = n_sites, ncol = n_seasons),
  gis_tree = matrix(yearly_covs$gis_tree, nrow = n_sites, ncol = n_seasons),
  gis_water = matrix(yearly_covs$gis_water, nrow = n_sites, ncol = n_seasons),
  lunar = matrix(yearly_covs$lunar, nrow = n_sites, ncol = n_seasons)
)

umf1 <- unmarkedMultFrame(
  y = y_matrix1,            # 2D detection matrix
  numPrimary = n_seasons,  # Number of primary periods (years)
  siteCovs = site_covs,    # Site-level covariates
  yearlySiteCovs = yearly_covs_list  # Yearly-varying covariates
)

summary(umf1)
# View detection history structure
image(umf1@y, main = "Detection History", 
      xlab = "Site", ylab = "Visit×Year")

#  Find Best Detection Model (p)
# Null model (constant detection)
m0 <- colext(psiformula = ~1,      # Initial occupancy
             gammaformula = ~1,    # Colonization
             epsilonformula = ~1,  # Extinction
             pformula = ~1,        # Detection (null)
             data = umf1)

# Weather models
m1 <- colext(~1, ~1, ~1, ~tmax, umf1)               # p ~ tmax
m2 <- colext(~1, ~1, ~1, ~tmin, umf1)        # p ~ tmin
m3 <- colext(~1, ~1, ~1, ~prcp, umf1)        # p ~ precipitation
m4 <- colext(~1, ~1, ~1, ~lunar, umf1)          # p ~ moon phase

# Site structure models
m5 <- colext(~1, ~1, ~1, ~structure, umf1)          # p ~ structure type
m6 <- colext(~1, ~1, ~1, ~canopy, umf1)             # p ~ canopy cover
m7 <- colext(~1, ~1, ~1, ~water, umf1)              # p ~ water presence

# GIS models
m8 <- colext(~1, ~1, ~1, ~gis_tree, umf1)   # p ~ gis tree 
m9 <- colext(~1, ~1, ~1, ~gis_water, umf1)      # p ~ gis water

det_models1 <- fitList(
  "p(.)" = m0,
  "p(tmax)" = m1,
  "p(tmin)" = m2,
  "p(prcp)" = m3,
  "p(lunar)" = m4,
  "p(structure)" = m5,
  "p(canopy)" = m6,
  "p(water)" = m7,
  "p(gis_tree)" = m8,
  "p(gis_water)" = m9
)

# Rank models by AIC
aic_table <- modSel(det_models1)
print(aic_table)

m6 #top model for EPFU
#graph

# Create prediction data
new_p <- expand.grid(
  canopy = factor(levels(umf1@siteCovs$canopy)[1:3],
                  levels = levels(umf1@siteCovs$canopy))
)
# Predict detection p
p_pred <- predict(m6, 
                  type = "det", 
                  newdata = new_p, 
                  se.fit = TRUE)

plot_p <- cbind(new_p, p_pred)


fig_epfu <- ggplot(plot_p, aes(x = canopy, y = Predicted)) +
  geom_point(size = 3, fill = "navy") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 1,
                alpha = 0.3,
                color = "steelblue") +
  scale_x_discrete(labels= c( "open" = "Open", 
                              "half" = "Intermediate",
                              "closed"= "Closed")) +
  labs(x = "Site canopy", y = "EPFU detection probabitliy") +
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )

##################labo####################
# Step 1: Create detection history and turn into detection matrix#
detection_history <- occu_data %>%
  # Select relevant columns (site, year, visit, species)
  select(site, cell, year, night, LASBOR) 
#replace EPTFUS in species
# epfu EPTFUS, labo LASBOR, laci LASCIN, lano LASNOC, mylu MYOLUC, #
# nyhu NYCHUM, pesu PERSUB, tabr TADBRA#
# mygr MYOGRI, myau MYOAUS, myse MYOSEP

# Convert to detection matrix - SAFE METHOD, not sure if efficient
y_matrix <- detection_history %>%
  pivot_wider(names_from = night, values_from = LASBOR,
              values_fill = NA, names_sort = TRUE) %>%
  arrange(site, cell, year)

# Create 3D array with proper dimension checking
n_sites <- n_distinct(y_matrix$site)
n_seasons <- n_distinct(y_matrix$year)
n_visits <- 4  # Fixed number of visits

y_array <- array(NA, dim = c(n_sites, n_visits, n_seasons))

for(i in seq_len(n_sites)){
  current_site <- unique(y_matrix$site)[i]
  
  for(t in seq_len(n_seasons)){
    current_year <- unique(y_matrix$year)[t]
    
    site_year_data <- y_matrix %>%
      filter(site == current_site, year == current_year)
    
    if(nrow(site_year_data) == 1){
      # Safely extract visit data
      visit_data <- as.numeric(site_year_data[1, paste0("night", 1:4)])
      
      # Verify length matches expected visits
      if(length(visit_data) == n_visits){
        y_array[i, , t] <- visit_data
      } else {
        warning(paste("Unexpected number of visits for site", current_site, 
                      "year", current_year))
      }
    }
    # Otherwise leaves as NA (for missing surveys)
  }
}

# Verify the array structure
str(y_array)
table(apply(y_array, 3, function(x) sum(!is.na(x)))) # Surveys per year

y_matrix1 <- matrix(y_array, 
                    nrow = dim(y_array)[1], 
                    ncol = dim(y_array)[2] * dim(y_array)[3])


# Site-level covariates (constant across years)
site_covs <- occu_data %>%
  distinct(site, structure, canopy, water, mean_2019_250, sd_2019_250, region, wns_sus) %>%
  arrange(site)

# Observation or year level covariates (change across years and nights)
yearly_covs <- occu_data[,c(1,2,32:34,49,50, 51)]

# Prepare yearlySiteCovs (must be named lists of matrices)
yearly_covs_list <- list(
  tmax = matrix(yearly_covs$tmax, nrow = n_sites, ncol = n_seasons),
  tmin = matrix(yearly_covs$tmin, nrow = n_sites, ncol = n_seasons),
  prcp = matrix(yearly_covs$prcp, nrow = n_sites, ncol = n_seasons),
  gis_tree = matrix(yearly_covs$gis_tree, nrow = n_sites, ncol = n_seasons),
  gis_water = matrix(yearly_covs$gis_water, nrow = n_sites, ncol = n_seasons),
  lunar = matrix(yearly_covs$lunar, nrow = n_sites, ncol = n_seasons)
)

umf1 <- unmarkedMultFrame(
  y = y_matrix1,            # 2D detection matrix
  numPrimary = n_seasons,  # Number of primary periods (years)
  siteCovs = site_covs,    # Site-level covariates
  yearlySiteCovs = yearly_covs_list  # Yearly-varying covariates
)

summary(umf1)
# View detection history structure
image(umf1@y, main = "Detection History", 
      xlab = "Site", ylab = "Visit×Year")


#  Find Best Detection Model (p)
# Null model (constant detection)
m0 <- colext(psiformula = ~1,      # Initial occupancy
             gammaformula = ~1,    # Colonization
             epsilonformula = ~1,  # Extinction
             pformula = ~1,        # Detection (null)
             data = umf1)

# Weather models
m1 <- colext(~1, ~1, ~1, ~tmax, umf1)               # p ~ tmax
m2 <- colext(~1, ~1, ~1, ~tmin, umf1)        # p ~ tmin
m3 <- colext(~1, ~1, ~1, ~prcp, umf1)        # p ~  precipitation
m4 <- colext(~1, ~1, ~1, ~lunar, umf1)          # p ~ moon phase

# Site structure models
m5 <- colext(~1, ~1, ~1, ~structure, umf1)          # p ~ structure type
m6 <- colext(~1, ~1, ~1, ~canopy, umf1)             # p ~ canopy cover
m7 <- colext(~1, ~1, ~1, ~water, umf1)              # p ~ water presence

# GIS models
m8 <- colext(~1, ~1, ~1, ~gis_tree, umf1)   # p ~ gis tree
m9 <- colext(~1, ~1, ~1, ~gis_water, umf1)      # p ~ gis water


det_models1 <- fitList(
  "p(.)" = m0,
  "p(tmax)" = m1,
  "p(tmin)" = m2,
  "p(prcp)" = m3,
  "p(lunar)" = m4,
  "p(structure)" = m5,
  "p(canopy)" = m6,
  "p(water)" = m7,
  "p(gis_tree)" = m8,
  "p(gis_water)" = m9
)

# Rank models by AIC
aic_table <- modSel(det_models1)
print(aic_table)

m6 #top model for labo
#graph

# Create prediction data
new_p <- expand.grid(
  canopy = factor(levels(umf1@siteCovs$canopy)[1:3],
                  levels = levels(umf1@siteCovs$canopy))
)
# Predict detection p
p_pred <- predict(m6, 
                  type = "det", 
                  newdata = new_p, 
                  se.fit = TRUE)

plot_p <- cbind(new_p, p_pred)


fig_labo <- ggplot(plot_p, aes(x = canopy, y = Predicted)) +
  geom_point(size = 3, fill = "navy") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 1,
                alpha = 0.3,
                color = "steelblue")+  
  scale_x_discrete(labels= c( "open" = "Open", 
                              "half" = "Intermediate",
                              "closed"= "Closed")) +
  labs(x = "Site canopy", y = "LABO detection probabitliy") +
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )

##################laci####################
# Step 1: Create detection history and turn into detection matrix#
detection_history <- occu_data %>%
  # Select relevant columns (site, year, visit, species)
  select(site, cell, year, night, LASCIN) 
#replace EPTFUS in species
# epfu EPTFUS, labo LASBOR, laci LASCIN, lano LASNOC, mylu MYOLUC, #
# nyhu NYCHUM, pesu PERSUB, tabr TADBRA#
# mygr MYOGRI, myau MYOAUS, myse MYOSEP

# Convert to detection matrix - SAFE METHOD, not sure if efficient
y_matrix <- detection_history %>%
  pivot_wider(names_from = night, values_from = LASCIN,
              values_fill = NA, names_sort = TRUE) %>%
  arrange(site, cell, year)

# Create 3D array with proper dimension checking
n_sites <- n_distinct(y_matrix$site)
n_seasons <- n_distinct(y_matrix$year)
n_visits <- 4  # Fixed number of visits

y_array <- array(NA, dim = c(n_sites, n_visits, n_seasons))

for(i in seq_len(n_sites)){
  current_site <- unique(y_matrix$site)[i]
  
  for(t in seq_len(n_seasons)){
    current_year <- unique(y_matrix$year)[t]
    
    site_year_data <- y_matrix %>%
      filter(site == current_site, year == current_year)
    
    if(nrow(site_year_data) == 1){
      # Safely extract visit data
      visit_data <- as.numeric(site_year_data[1, paste0("night", 1:4)])
      
      # Verify length matches expected visits
      if(length(visit_data) == n_visits){
        y_array[i, , t] <- visit_data
      } else {
        warning(paste("Unexpected number of visits for site", current_site, 
                      "year", current_year))
      }
    }
    # Otherwise leaves as NA (for missing surveys)
  }
}

# Verify the array structure
str(y_array)
table(apply(y_array, 3, function(x) sum(!is.na(x)))) # Surveys per year

y_matrix1 <- matrix(y_array, 
                    nrow = dim(y_array)[1], 
                    ncol = dim(y_array)[2] * dim(y_array)[3])


# Site-level covariates (constant across years)
site_covs <- occu_data %>%
  distinct(site, structure, canopy, water, mean_2019_250, sd_2019_250, region, wns_sus) %>%
  arrange(site)

# Observation or year level covariates (change across years and nights)
yearly_covs <- occu_data[,c(1,2,32:34,49,50, 51)]

# Prepare yearlySiteCovs (must be named lists of matrices)
yearly_covs_list <- list(
  tmax = matrix(yearly_covs$tmax, nrow = n_sites, ncol = n_seasons),
  tmin = matrix(yearly_covs$tmin, nrow = n_sites, ncol = n_seasons),
  prcp = matrix(yearly_covs$prcp, nrow = n_sites, ncol = n_seasons),
  gis_tree = matrix(yearly_covs$gis_tree, nrow = n_sites, ncol = n_seasons),
  gis_water = matrix(yearly_covs$gis_water, nrow = n_sites, ncol = n_seasons),
  lunar = matrix(yearly_covs$lunar, nrow = n_sites, ncol = n_seasons)
)

umf1 <- unmarkedMultFrame(
  y = y_matrix1,            # 2D detection matrix
  numPrimary = n_seasons,  # Number of primary periods (years)
  siteCovs = site_covs,    # Site-level covariates
  yearlySiteCovs = yearly_covs_list  # Yearly-varying covariates
)

summary(umf1)
# View detection history structure
image(umf1@y, main = "Detection History", 
      xlab = "Site", ylab = "Visit×Year")


# Find Best Detection Model (p)
# Null model (constant detection)
m0 <- colext(psiformula = ~1,      # Initial occupancy
             gammaformula = ~1,    # Colonization
             epsilonformula = ~1,  # Extinction
             pformula = ~1,        # Detection (null)
             data = umf1)

# Weather models
m1 <- colext(~1, ~1, ~1, ~tmax, umf1)               # p ~ tmax
m2 <- colext(~1, ~1, ~1, ~tmin, umf1)        # p ~  tmin
m3 <- colext(~1, ~1, ~1, ~prcp, umf1)        # p ~ precipitation
m4 <- colext(~1, ~1, ~1, ~lunar, umf1)          # p ~ lunar

# Site structure models
m5 <- colext(~1, ~1, ~1, ~structure, umf1)          # p ~ structure type
m6 <- colext(~1, ~1, ~1, ~canopy, umf1)             # p ~ canopy cover
m7 <- colext(~1, ~1, ~1, ~water, umf1)              # p ~ water presence

# GIS models
m8 <- colext(~1, ~1, ~1, ~gis_tree, umf1)   # p ~ gis tree
m9 <- colext(~1, ~1, ~1, ~gis_water, umf1)      # p ~ gis water

det_models1 <- fitList(
  "p(.)" = m0,
  "p(tmax)" = m1,
  "p(tmin)" = m2,
  "p(prcp)" = m3,
  "p(lunar)" = m4,
  "p(structure)" = m5,
  "p(canopy)" = m6,
  "p(water)" = m7,
  "p(gis_tree)" = m8,
  "p(gis_water)" = m9
)

# Rank models by AIC
aic_table <- modSel(det_models1)
print(aic_table)

m6 #top model for LACI
#graph

# Create prediction data
new_p <- expand.grid(
  canopy = factor(levels(umf1@siteCovs$canopy)[1:3],
                  levels = levels(umf1@siteCovs$canopy))
)
# Predict detection p
p_pred <- predict(m6, 
                  type = "det", 
                  newdata = new_p, 
                  se.fit = TRUE)

plot_p <- cbind(new_p, p_pred)

fig_laci <- ggplot(plot_p, aes(x = canopy, y = Predicted)) +
  geom_point(size = 3, fill = "navy") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 1,
                alpha = 0.3,
                color = "steelblue") +
  scale_x_discrete(labels= c( "open" = "Open", 
                              "half" = "Intermediate",
                              "closed"= "Closed")) +
  labs(x = "Site canopy", y = "LACI detection probabitliy") +  
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )

##################lano####################
# Create detection history and turn into detection matrix#
detection_history <- occu_data %>%
  # Select relevant columns (site, year, visit, species)
  select(site, cell, year, night, LASNOC) 
#replace EPTFUS in species
# epfu EPTFUS, labo LASBOR, laci LASCIN, lano LASNOC, mylu MYOLUC, #
# nyhu NYCHUM, pesu PERSUB, tabr TADBRA#
# mygr MYOGRI, myau MYOAUS, myse MYOSEP

# Convert to detection matrix - SAFE METHOD, not sure if efficient
y_matrix <- detection_history %>%
  pivot_wider(names_from = night, values_from = LASNOC,
              values_fill = NA, names_sort = TRUE) %>%
  arrange(site, cell, year)

# Create 3D array with proper dimension checking
n_sites <- n_distinct(y_matrix$site)
n_seasons <- n_distinct(y_matrix$year)
n_visits <- 4  # Fixed number of visits

y_array <- array(NA, dim = c(n_sites, n_visits, n_seasons))

for(i in seq_len(n_sites)){
  current_site <- unique(y_matrix$site)[i]
  
  for(t in seq_len(n_seasons)){
    current_year <- unique(y_matrix$year)[t]
    
    site_year_data <- y_matrix %>%
      filter(site == current_site, year == current_year)
    
    if(nrow(site_year_data) == 1){
      # Safely extract visit data
      visit_data <- as.numeric(site_year_data[1, paste0("night", 1:4)])
      
      # Verify length matches expected visits
      if(length(visit_data) == n_visits){
        y_array[i, , t] <- visit_data
      } else {
        warning(paste("Unexpected number of visits for site", current_site, 
                      "year", current_year))
      }
    }
    # Otherwise leaves as NA (for missing surveys)
  }
}

# Verify the array structure
str(y_array)
table(apply(y_array, 3, function(x) sum(!is.na(x)))) # Surveys per year

y_matrix1 <- matrix(y_array, 
                    nrow = dim(y_array)[1], 
                    ncol = dim(y_array)[2] * dim(y_array)[3])


# Site-level covariates (constant across years)
site_covs <- occu_data %>%
  distinct(site, structure, canopy, water, mean_2019_250, sd_2019_250, region, wns_sus) %>%
  arrange(site)

# Observation or year level covariates (change across years and nights)
yearly_covs <- occu_data[,c(1,2,32:34,49,50, 51)]

# Prepare yearlySiteCovs (must be named lists of matrices)
yearly_covs_list <- list(
  tmax = matrix(yearly_covs$tmax, nrow = n_sites, ncol = n_seasons),
  tmin = matrix(yearly_covs$tmin, nrow = n_sites, ncol = n_seasons),
  prcp = matrix(yearly_covs$prcp, nrow = n_sites, ncol = n_seasons),
  gis_tree = matrix(yearly_covs$gis_tree, nrow = n_sites, ncol = n_seasons),
  gis_water = matrix(yearly_covs$gis_water, nrow = n_sites, ncol = n_seasons),
  lunar = matrix(yearly_covs$lunar, nrow = n_sites, ncol = n_seasons)
)

umf1 <- unmarkedMultFrame(
  y = y_matrix1,            # 2D detection matrix
  numPrimary = n_seasons,  # Number of primary periods (years)
  siteCovs = site_covs,    # Site-level covariates
  yearlySiteCovs = yearly_covs_list  # Yearly-varying covariates
)

summary(umf1)
# View detection history structure
image(umf1@y, main = "Detection History", 
      xlab = "Site", ylab = "Visit×Year")


# Find Best Detection Model (p)
# Null model (constant detection)
m0 <- colext(psiformula = ~1,      # Initial occupancy
             gammaformula = ~1,    # Colonization
             epsilonformula = ~1,  # Extinction
             pformula = ~1,        # Detection (null)
             data = umf1)

# Weather models
m1 <- colext(~1, ~1, ~1, ~tmax, umf1)               # p ~ tmax
m2 <- colext(~1, ~1, ~1, ~tmin, umf1)        # p ~ tmin
m3 <- colext(~1, ~1, ~1, ~prcp, umf1)        # p ~  precipitation
m4 <- colext(~1, ~1, ~1, ~lunar, umf1)          # p ~ moon phase

# Site structure models
m5 <- colext(~1, ~1, ~1, ~structure, umf1)          # p ~ structure type
m6 <- colext(~1, ~1, ~1, ~canopy, umf1)             # p ~ canopy cover
m7 <- colext(~1, ~1, ~1, ~water, umf1)              # p ~ water presence

# GIS models
m8 <- colext(~1, ~1, ~1, ~gis_tree, umf1)   # p ~ gis tree
m9 <- colext(~1, ~1, ~1, ~gis_water, umf1)      # p ~ gis water

det_models1 <- fitList(
  "p(.)" = m0,
  "p(tmax)" = m1,
  "p(tmin)" = m2,
  "p(prcp)" = m3,
  "p(lunar)" = m4,
  "p(structure)" = m5,
  "p(canopy)" = m6,
  "p(water)" = m7,
  "p(gis_tree)" = m8,
  "p(gis_water)" = m9
)

# Rank models by AIC
aic_table <- modSel(det_models1)
print(aic_table)

m6 #top model for LANO
#graph
# Create prediction data
new_p <- expand.grid(
  canopy = factor(levels(umf1@siteCovs$canopy)[1:3],
                  levels = levels(umf1@siteCovs$canopy))
)
# Predict detection p
p_pred <- predict(m6, 
                  type = "det", 
                  newdata = new_p, 
                  se.fit = TRUE)

plot_p <- cbind(new_p, p_pred)


fig_lano <- ggplot(plot_p, aes(x = canopy, y = Predicted)) +
  geom_point(size = 3, fill = "navy") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 1,
                alpha = 0.3,
                color = "steelblue") +
  scale_x_discrete(labels= c( "open" = "Open", 
                              "half" = "Intermediate",
                              "closed"= "Closed")) +
  labs(x = "Site canopy", y = "LANO detection probabitliy") +  
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )

##################nyhu####################
# Create detection history and turn into detection matrix#
detection_history <- occu_data %>%
  # Select relevant columns (site, year, visit, species)
  select(site, cell, year, night, NYCHUM) 
#replace EPTFUS in species
# epfu EPTFUS, labo LASBOR, laci LASCIN, lano LASNOC, mylu MYOLUC, #
# nyhu NYCHUM, pesu PERSUB, tabr TADBRA#
# mygr MYOGRI, myau MYOAUS, myse MYOSEP

# Convert to detection matrix - SAFE METHOD, not sure if efficient
y_matrix <- detection_history %>%
  pivot_wider(names_from = night, values_from = NYCHUM,
              values_fill = NA, names_sort = TRUE) %>%
  arrange(site, cell, year)

# Create 3D array with proper dimension checking
n_sites <- n_distinct(y_matrix$site)
n_seasons <- n_distinct(y_matrix$year)
n_visits <- 4  # Fixed number of visits

y_array <- array(NA, dim = c(n_sites, n_visits, n_seasons))

for(i in seq_len(n_sites)){
  current_site <- unique(y_matrix$site)[i]
  
  for(t in seq_len(n_seasons)){
    current_year <- unique(y_matrix$year)[t]
    
    site_year_data <- y_matrix %>%
      filter(site == current_site, year == current_year)
    
    if(nrow(site_year_data) == 1){
      # Safely extract visit data
      visit_data <- as.numeric(site_year_data[1, paste0("night", 1:4)])
      
      # Verify length matches expected visits
      if(length(visit_data) == n_visits){
        y_array[i, , t] <- visit_data
      } else {
        warning(paste("Unexpected number of visits for site", current_site, 
                      "year", current_year))
      }
    }
    # Otherwise leaves as NA (for missing surveys)
  }
}

# Verify the array structure
str(y_array)
table(apply(y_array, 3, function(x) sum(!is.na(x)))) # Surveys per year

y_matrix1 <- matrix(y_array, 
                    nrow = dim(y_array)[1], 
                    ncol = dim(y_array)[2] * dim(y_array)[3])


# Site-level covariates (constant across years)
site_covs <- occu_data %>%
  distinct(site, structure, canopy, water, mean_2019_250, sd_2019_250, region, wns_sus) %>%
  arrange(site)

# Observation or year level covariates (change across years and nights)
yearly_covs <- occu_data[,c(1,2,32:34,49,50, 51)]

# Prepare yearlySiteCovs (must be named lists of matrices)
yearly_covs_list <- list(
  tmax = matrix(yearly_covs$tmax, nrow = n_sites, ncol = n_seasons),
  tmin = matrix(yearly_covs$tmin, nrow = n_sites, ncol = n_seasons),
  prcp = matrix(yearly_covs$prcp, nrow = n_sites, ncol = n_seasons),
  gis_tree = matrix(yearly_covs$gis_tree, nrow = n_sites, ncol = n_seasons),
  gis_water = matrix(yearly_covs$gis_water, nrow = n_sites, ncol = n_seasons),
  lunar = matrix(yearly_covs$lunar, nrow = n_sites, ncol = n_seasons)
)

umf1 <- unmarkedMultFrame(
  y = y_matrix1,            # 2D detection matrix
  numPrimary = n_seasons,  # Number of primary periods (years)
  siteCovs = site_covs,    # Site-level covariates
  yearlySiteCovs = yearly_covs_list  # Yearly-varying covariates
)

summary(umf1)
# View detection history structure
image(umf1@y, main = "Detection History", 
      xlab = "Site", ylab = "Visit×Year")

# Find Best Detection Model (p)
# Null model (constant detection)
m0 <- colext(psiformula = ~1,      # Initial occupancy
             gammaformula = ~1,    # Colonization
             epsilonformula = ~1,  # Extinction
             pformula = ~1,        # Detection (null)
             data = umf1)

# Weather models
m1 <- colext(~1, ~1, ~1, ~tmax, umf1)               # p ~ tmax
m2 <- colext(~1, ~1, ~1, ~tmin, umf1)        # p ~ tmin
m3 <- colext(~1, ~1, ~1, ~prcp, umf1)        # p ~ precipitation
m4 <- colext(~1, ~1, ~1, ~lunar, umf1)          # p ~ moon phase

# Site structure models
m5 <- colext(~1, ~1, ~1, ~structure, umf1)          # p ~ structure type
m6 <- colext(~1, ~1, ~1, ~canopy, umf1)             # p ~ canopy cover
m7 <- colext(~1, ~1, ~1, ~water, umf1)              # p ~ water presence

# GIS models
m8 <- colext(~1, ~1, ~1, ~gis_tree, umf1)   # p ~ gis tree
m9 <- colext(~1, ~1, ~1, ~gis_water, umf1)      # p ~ gis water

det_models1 <- fitList(
  "p(.)" = m0,
  "p(tmax)" = m1,
  "p(tmin)" = m2,
  "p(prcp)" = m3,
  "p(lunar)" = m4,
  "p(structure)" = m5,
  "p(canopy)" = m6,
  "p(water)" = m7,
  "p(gis_tree)" = m8,
  "p(gis_water)" = m9
)

# Rank models by AIC
aic_table <- modSel(det_models1)
print(aic_table)

m9 #top model for NYHU
#graph

# Create prediction data
new_p <- expand.grid(
  gis_water = seq(min(umf1@yearlySiteCovs$gis_water, na.rm = TRUE),
             max(umf1@yearlySiteCovs$gis_water, na.rm = TRUE),
             length.out=100)
)

# Predict detection p
p_pred <- predict(m9, 
                  type = "det", 
                  newdata = new_p, 
                  se.fit = TRUE)

plot_p <- cbind(new_p, p_pred)

fig_nyhu <- ggplot(plot_p, aes(x = gis_water, y = Predicted)) +
  geom_line(linewidth = 1, color = "navy") +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.3, 
              fill = "steelblue") +
  labs(x = "GIS derived water cover percentage", y = "NYHU detection probabitliy") +
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )

##################pesu####################
# Create detection history and turn into detection matrix#
detection_history <- occu_data %>%
  # Select relevant columns (site, year, visit, species)
  select(site, cell, year, night, PERSUB) 
#replace EPTFUS in species
# epfu EPTFUS, labo LASBOR, laci LASCIN, lano LASNOC, mylu MYOLUC, #
# nyhu NYCHUM, pesu PERSUB, tabr TADBRA#
# mygr MYOGRI, myau MYOAUS, myse MYOSEP

# Convert to detection matrix - SAFE METHOD, not sure if efficient
y_matrix <- detection_history %>%
  pivot_wider(names_from = night, values_from = PERSUB,
              values_fill = NA, names_sort = TRUE) %>%
  arrange(site, cell, year)

# Create 3D array with proper dimension checking
n_sites <- n_distinct(y_matrix$site)
n_seasons <- n_distinct(y_matrix$year)
n_visits <- 4  # Fixed number of visits

y_array <- array(NA, dim = c(n_sites, n_visits, n_seasons))

for(i in seq_len(n_sites)){
  current_site <- unique(y_matrix$site)[i]
  
  for(t in seq_len(n_seasons)){
    current_year <- unique(y_matrix$year)[t]
    
    site_year_data <- y_matrix %>%
      filter(site == current_site, year == current_year)
    
    if(nrow(site_year_data) == 1){
      # Safely extract visit data
      visit_data <- as.numeric(site_year_data[1, paste0("night", 1:4)])
      
      # Verify length matches expected visits
      if(length(visit_data) == n_visits){
        y_array[i, , t] <- visit_data
      } else {
        warning(paste("Unexpected number of visits for site", current_site, 
                      "year", current_year))
      }
    }
    # Otherwise leaves as NA (for missing surveys)
  }
}

# Verify the array structure
str(y_array)
table(apply(y_array, 3, function(x) sum(!is.na(x)))) # Surveys per year

y_matrix1 <- matrix(y_array, 
                    nrow = dim(y_array)[1], 
                    ncol = dim(y_array)[2] * dim(y_array)[3])


# Site-level covariates (constant across years)
site_covs <- occu_data %>%
  distinct(site, structure, canopy, water, mean_2019_250, sd_2019_250, region, wns_sus) %>%
  arrange(site)

# Observation or year level covariates (change across years and nights)
yearly_covs <- occu_data[,c(1,2,32:34,49,50, 51)]

# Prepare yearlySiteCovs (must be named lists of matrices)
yearly_covs_list <- list(
  tmax = matrix(yearly_covs$tmax, nrow = n_sites, ncol = n_seasons),
  tmin = matrix(yearly_covs$tmin, nrow = n_sites, ncol = n_seasons),
  prcp = matrix(yearly_covs$prcp, nrow = n_sites, ncol = n_seasons),
  gis_tree = matrix(yearly_covs$gis_tree, nrow = n_sites, ncol = n_seasons),
  gis_water = matrix(yearly_covs$gis_water, nrow = n_sites, ncol = n_seasons),
  lunar = matrix(yearly_covs$lunar, nrow = n_sites, ncol = n_seasons)
)

umf1 <- unmarkedMultFrame(
  y = y_matrix1,            # 2D detection matrix
  numPrimary = n_seasons,  # Number of primary periods (years)
  siteCovs = site_covs,    # Site-level covariates
  yearlySiteCovs = yearly_covs_list  # Yearly-varying covariates
)

summary(umf1)
# View detection history structure
image(umf1@y, main = "Detection History", 
      xlab = "Site", ylab = "Visit×Year")

# Find Best Detection Model (p)
# Null model (constant detection)
m0 <- colext(psiformula = ~1,      # Initial occupancy
             gammaformula = ~1,    # Colonization
             epsilonformula = ~1,  # Extinction
             pformula = ~1,        # Detection (null)
             data = umf1)

# Weather models
m1 <- colext(~1, ~1, ~1, ~tmax, umf1)               # p ~ tmax
m2 <- colext(~1, ~1, ~1, ~tmin, umf1)        # p ~ tmin
m3 <- colext(~1, ~1, ~1, ~prcp, umf1)        # p ~ precipitation
m4 <- colext(~1, ~1, ~1, ~lunar, umf1)          # p ~ moon phase

# Site structure models
m5 <- colext(~1, ~1, ~1, ~structure, umf1)          # p ~ structure type
m6 <- colext(~1, ~1, ~1, ~canopy, umf1)             # p ~ canopy cover
m7 <- colext(~1, ~1, ~1, ~water, umf1)              # p ~ water presence

# GIS models
m8 <- colext(~1, ~1, ~1, ~gis_tree, umf1)   # p ~ gis tree
m9 <- colext(~1, ~1, ~1, ~gis_water, umf1)      # p ~ gis water

det_models1 <- fitList(
  "p(.)" = m0,
  "p(tmax)" = m1,
  "p(tmin)" = m2,
  "p(prcp)" = m3,
  "p(lunar)" = m4,
  "p(structure)" = m5,
  "p(canopy)" = m6,
  "p(water)" = m7,
  "p(gis_tree)" = m8,
  "p(gis_water)" = m9
)

# Rank models by AIC
aic_table <- modSel(det_models1)
print(aic_table)

m7 #top model for PESU
#graph
# Create prediction data
new_p <- expand.grid(
  water = factor(levels(umf1@siteCovs$water)[1:2],
                  levels = levels(umf1@siteCovs$water))
)
# Predict detection p
p_pred <- predict(m7, 
                  type = "det", 
                  newdata = new_p, 
                  se.fit = TRUE)

plot_p <- cbind(new_p, p_pred)


fig_pesu <- ggplot(plot_p, aes(x = water, y = Predicted)) +
  geom_point(size = 3, fill = "navy") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 1,
                alpha = 0.3,
                color = "steelblue") +
  scale_x_discrete(labels= c( "yes" = "Facing water", 
                              "no" = "Not facing water")) +
  labs(x = "Water facing", y = "PESU detection probabitliy") +  
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )

##################tabr####################
#Create detection history and turn into detection matrix#
detection_history <- occu_data %>%
  # Select relevant columns (site, year, visit, species)
  select(site, cell, year, night, TADBRA) 
#replace EPTFUS in species
# epfu EPTFUS, labo LASBOR, laci LASCIN, lano LASNOC, mylu MYOLUC, #
# nyhu NYCHUM, pesu PERSUB, tabr TADBRA#
# mygr MYOGRI, myau MYOAUS, myse MYOSEP

# Convert to detection matrix - SAFE METHOD, not sure if efficient
y_matrix <- detection_history %>%
  pivot_wider(names_from = night, values_from = TADBRA,
              values_fill = NA, names_sort = TRUE) %>%
  arrange(site, cell, year)

# Create 3D array with proper dimension checking
n_sites <- n_distinct(y_matrix$site)
n_seasons <- n_distinct(y_matrix$year)
n_visits <- 4  # Fixed number of visits

y_array <- array(NA, dim = c(n_sites, n_visits, n_seasons))

for(i in seq_len(n_sites)){
  current_site <- unique(y_matrix$site)[i]
  
  for(t in seq_len(n_seasons)){
    current_year <- unique(y_matrix$year)[t]
    
    site_year_data <- y_matrix %>%
      filter(site == current_site, year == current_year)
    
    if(nrow(site_year_data) == 1){
      # Safely extract visit data
      visit_data <- as.numeric(site_year_data[1, paste0("night", 1:4)])
      
      # Verify length matches expected visits
      if(length(visit_data) == n_visits){
        y_array[i, , t] <- visit_data
      } else {
        warning(paste("Unexpected number of visits for site", current_site, 
                      "year", current_year))
      }
    }
    # Otherwise leaves as NA (for missing surveys)
  }
}

# Verify the array structure
str(y_array)
table(apply(y_array, 3, function(x) sum(!is.na(x)))) # Surveys per year

y_matrix1 <- matrix(y_array, 
                    nrow = dim(y_array)[1], 
                    ncol = dim(y_array)[2] * dim(y_array)[3])


# Site-level covariates (constant across years)
site_covs <- occu_data %>%
  distinct(site, structure, canopy, water, mean_2019_250, sd_2019_250, region, wns_sus) %>%
  arrange(site)

# Observation or year level covariates (change across years and nights)
yearly_covs <- occu_data[,c(1,2,32:34,49,50, 51)]

# Prepare yearlySiteCovs (must be named lists of matrices)
yearly_covs_list <- list(
  tmax = matrix(yearly_covs$tmax, nrow = n_sites, ncol = n_seasons),
  tmin = matrix(yearly_covs$tmin, nrow = n_sites, ncol = n_seasons),
  prcp = matrix(yearly_covs$prcp, nrow = n_sites, ncol = n_seasons),
  gis_tree = matrix(yearly_covs$gis_tree, nrow = n_sites, ncol = n_seasons),
  gis_water = matrix(yearly_covs$gis_water, nrow = n_sites, ncol = n_seasons),
  lunar = matrix(yearly_covs$lunar, nrow = n_sites, ncol = n_seasons)
)

umf1 <- unmarkedMultFrame(
  y = y_matrix1,            # 2D detection matrix
  numPrimary = n_seasons,  # Number of primary periods (years)
  siteCovs = site_covs,    # Site-level covariates
  yearlySiteCovs = yearly_covs_list  # Yearly-varying covariates
)

summary(umf1)
# View detection history structure
image(umf1@y, main = "Detection History", 
      xlab = "Site", ylab = "Visit×Year")


# Find Best Detection Model (p)
# Null model (constant detection)
m0 <- colext(psiformula = ~1,      # Initial occupancy
             gammaformula = ~1,    # Colonization
             epsilonformula = ~1,  # Extinction
             pformula = ~1,        # Detection (null)
             data = umf1)

# Weather models
m1 <- colext(~1, ~1, ~1, ~tmax, umf1)               # p ~ tmax
m2 <- colext(~1, ~1, ~1, ~tmin, umf1)        # p ~ tmin
m3 <- colext(~1, ~1, ~1, ~prcp, umf1)        # p ~ precipitation
m4 <- colext(~1, ~1, ~1, ~lunar, umf1)          # p ~ moon phase

# Site structure models
m5 <- colext(~1, ~1, ~1, ~structure, umf1)          # p ~ structure type
m6 <- colext(~1, ~1, ~1, ~canopy, umf1)             # p ~ canopy cover
m7 <- colext(~1, ~1, ~1, ~water, umf1)              # p ~ water presence

# GIS models
m8 <- colext(~1, ~1, ~1, ~gis_tree, umf1)   # p ~ gis tree
m9 <- colext(~1, ~1, ~1, ~gis_water, umf1)      # p ~ gis water

det_models1 <- fitList(
  "p(.)" = m0,
  "p(tmax)" = m1,
  "p(tmin)" = m2,
  "p(prcp)" = m3,
  "p(lunar)" = m4,
  "p(structure)" = m5,
  "p(canopy)" = m6,
  "p(water)" = m7,
  "p(gis_tree)" = m8,
  "p(gis_water)" = m9
)

# Rank models by AIC
aic_table <- modSel(det_models1)
print(aic_table)

m6
#graph
new_p <- expand.grid(
  canopy = factor(levels(umf1@siteCovs$canopy)[1:3],
                  levels = levels(umf1@siteCovs$canopy))
)
# Predict detection p
p_pred <- predict(m6, 
                  type = "det", 
                  newdata = new_p, 
                  se.fit = TRUE)

plot_p <- cbind(new_p, p_pred)


fig_tabr <- ggplot(plot_p, aes(x = canopy, y = Predicted)) +
  geom_point(size = 3, fill = "navy") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 1,
                alpha = 0.3,
                color = "steelblue") +
  scale_x_discrete(labels= c( "open" = "Open", 
                              "half" = "Intermediate",
                              "closed"= "Closed")) +
  labs(x = "Site canopy", y = "TABR detection probabitliy") +
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )
  
##################mylu####################
# Create detection history and turn into detection matrix#
detection_history <- occu_data %>%
  # Select relevant columns (site, year, visit, species)
  select(site, cell, year, night, MYOLUC) 
#replace EPTFUS in species
# epfu EPTFUS, labo LASBOR, laci LASCIN, lano LASNOC, mylu MYOLUC, #
# nyhu NYCHUM, pesu PERSUB, tabr TADBRA#
# mygr MYOGRI, myau MYOAUS, myse MYOSEP

# Convert to detection matrix - SAFE METHOD, not sure if efficient
y_matrix <- detection_history %>%
  pivot_wider(names_from = night, values_from = MYOLUC,
              values_fill = NA, names_sort = TRUE) %>%
  arrange(site, cell, year)

# Create 3D array with proper dimension checking
n_sites <- n_distinct(y_matrix$site)
n_seasons <- n_distinct(y_matrix$year)
n_visits <- 4  # Fixed number of visits

y_array <- array(NA, dim = c(n_sites, n_visits, n_seasons))

for(i in seq_len(n_sites)){
  current_site <- unique(y_matrix$site)[i]
  
  for(t in seq_len(n_seasons)){
    current_year <- unique(y_matrix$year)[t]
    
    site_year_data <- y_matrix %>%
      filter(site == current_site, year == current_year)
    
    if(nrow(site_year_data) == 1){
      # Safely extract visit data
      visit_data <- as.numeric(site_year_data[1, paste0("night", 1:4)])
      
      # Verify length matches expected visits
      if(length(visit_data) == n_visits){
        y_array[i, , t] <- visit_data
      } else {
        warning(paste("Unexpected number of visits for site", current_site, 
                      "year", current_year))
      }
    }
    # Otherwise leaves as NA (for missing surveys)
  }
}

# Verify the array structure
str(y_array)
table(apply(y_array, 3, function(x) sum(!is.na(x)))) # Surveys per year

y_matrix1 <- matrix(y_array, 
                    nrow = dim(y_array)[1], 
                    ncol = dim(y_array)[2] * dim(y_array)[3])


# Site-level covariates (constant across years)
site_covs <- occu_data %>%
  distinct(site, structure, canopy, water, mean_2019_250, sd_2019_250, region, wns_sus) %>%
  arrange(site)

# Observation or year level covariates (change across years and nights)
yearly_covs <- occu_data[,c(1,2,32:34,49,50, 51)]

# Prepare yearlySiteCovs (must be named lists of matrices)
yearly_covs_list <- list(
  tmax = matrix(yearly_covs$tmax, nrow = n_sites, ncol = n_seasons),
  tmin = matrix(yearly_covs$tmin, nrow = n_sites, ncol = n_seasons),
  prcp = matrix(yearly_covs$prcp, nrow = n_sites, ncol = n_seasons),
  gis_tree = matrix(yearly_covs$gis_tree, nrow = n_sites, ncol = n_seasons),
  gis_water = matrix(yearly_covs$gis_water, nrow = n_sites, ncol = n_seasons),
  lunar = matrix(yearly_covs$lunar, nrow = n_sites, ncol = n_seasons)
)

umf1 <- unmarkedMultFrame(
  y = y_matrix1,            # 2D detection matrix
  numPrimary = n_seasons,  # Number of primary periods (years)
  siteCovs = site_covs,    # Site-level covariates
  yearlySiteCovs = yearly_covs_list  # Yearly-varying covariates
)

summary(umf1)
# View detection history structure
image(umf1@y, main = "Detection History", 
      xlab = "Site", ylab = "Visit×Year")

# Find Best Detection Model (p)
# Null model (constant detection)
m0 <- colext(psiformula = ~1,      # Initial occupancy
             gammaformula = ~1,    # Colonization
             epsilonformula = ~1,  # Extinction
             pformula = ~1,        # Detection (null)
             data = umf1)

# Weather models
m1 <- colext(~1, ~1, ~1, ~tmax, umf1)               # p ~ tmax
m2 <- colext(~1, ~1, ~1, ~tmin, umf1)        # p ~ tmin
m3 <- colext(~1, ~1, ~1, ~prcp, umf1)        # p ~ precipitation
m4 <- colext(~1, ~1, ~1, ~lunar, umf1)          # p ~ moon phase

# Site structure models
m5 <- colext(~1, ~1, ~1, ~structure, umf1)          # p ~ structure type
m6 <- colext(~1, ~1, ~1, ~canopy, umf1)             # p ~ canopy cover
m7 <- colext(~1, ~1, ~1, ~water, umf1)              # p ~ water presence

# GIS models
m8 <- colext(~1, ~1, ~1, ~gis_tree, umf1)   # p ~ gis tree
m9 <- colext(~1, ~1, ~1, ~gis_water, umf1)      # p ~ gis water

det_models1 <- fitList(
  "p(.)" = m0,
  "p(tmax)" = m1,
  "p(tmin)" = m2,
  "p(prcp)" = m3,
  "p(lunar)" = m4,
  "p(structure)" = m5,
  "p(canopy)" = m6,
  "p(water)" = m7,
  "p(gis_tree)" = m8,
  "p(gis_water)" = m9
)

# Rank models by AIC
aic_table <- modSel(det_models1)
print(aic_table)

m5 #top model for mylu
#graph
new_p <- expand.grid(
  structure = factor(levels(umf1@siteCovs$structure)[1:4],
                  levels = levels(umf1@siteCovs$structure))
)
# Predict detection p
p_pred <- predict(m5, 
                  type = "det", 
                  newdata = new_p, 
                  se.fit = TRUE)

plot_p <- cbind(new_p, p_pred)


fig_mylu <- ggplot(plot_p, aes(x = structure, y = Predicted)) +
  geom_point(size = 3, fill = "navy") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 1,
                alpha = 0.3,
                color = "steelblue") +
  scale_x_discrete(labels= c( "interior" = "Interior", 
                              "corridor" = "Corridor", 
                              "edge" = "Edge",
                              "field" = "Field"), limits = rev) +
  labs(x = "Site structure", y = "MYLU detection probabitliy") +
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )

##################myse####################
#only present in mountain and coastal regions#
occu_data1 <- subset(dfo, dfo$region=="66"|dfo$l3_key=="63  Middle Atlantic Coastal Plain")
unique(occu_data1$region_name)

#Create detection history and turn into detection matrix#
detection_history <- occu_data1 %>%
  # Select relevant columns (site, year, visit, species)
  select(site, cell, year, night, MYOSEP) 
#replace EPTFUS in species
# epfu EPTFUS, labo LASBOR, laci LASCIN, lano LASNOC, mylu MYOLUC, #
# nyhu NYCHUM, pesu PERSUB, tabr TADBRA#
# mygr MYOGRI, myau MYOAUS, myse MYOSEP

# Convert to detection matrix - SAFE METHOD, not sure if efficient
y_matrix <- detection_history %>%
  pivot_wider(names_from = night, values_from = MYOSEP,
              values_fill = NA, names_sort = TRUE) %>%
  arrange(site, cell, year)

# Create 3D array with proper dimension checking
n_sites <- n_distinct(y_matrix$site)
n_seasons <- n_distinct(y_matrix$year)
n_visits <- 4  # Fixed number of visits

y_array <- array(NA, dim = c(n_sites, n_visits, n_seasons))

for(i in seq_len(n_sites)){
  current_site <- unique(y_matrix$site)[i]
  
  for(t in seq_len(n_seasons)){
    current_year <- unique(y_matrix$year)[t]
    
    site_year_data <- y_matrix %>%
      filter(site == current_site, year == current_year)
    
    if(nrow(site_year_data) == 1){
      # Safely extract visit data
      visit_data <- as.numeric(site_year_data[1, paste0("night", 1:4)])
      
      # Verify length matches expected visits
      if(length(visit_data) == n_visits){
        y_array[i, , t] <- visit_data
      } else {
        warning(paste("Unexpected number of visits for site", current_site, 
                      "year", current_year))
      }
    }
    # Otherwise leaves as NA (for missing surveys)
  }
}

# Verify the array structure
str(y_array)
table(apply(y_array, 3, function(x) sum(!is.na(x)))) # Surveys per year

y_matrix1 <- matrix(y_array, 
                    nrow = dim(y_array)[1], 
                    ncol = dim(y_array)[2] * dim(y_array)[3])


# Site-level covariates (constant across years)
site_covs <- occu_data1 %>%
  distinct(site, structure, canopy, water, mean_2019_250, sd_2019_250, region, wns_sus) %>%
  arrange(site)

# Observation or year level covariates (change across years and nights)
yearly_covs <- occu_data1[,c(1,2,32:34,49,50, 51)]

# Prepare yearlySiteCovs (must be named lists of matrices)
yearly_covs_list <- list(
  tmax = matrix(yearly_covs$tmax, nrow = n_sites, ncol = n_seasons),
  tmin = matrix(yearly_covs$tmin, nrow = n_sites, ncol = n_seasons),
  prcp = matrix(yearly_covs$prcp, nrow = n_sites, ncol = n_seasons),
  gis_tree = matrix(yearly_covs$gis_tree, nrow = n_sites, ncol = n_seasons),
  gis_water = matrix(yearly_covs$gis_water, nrow = n_sites, ncol = n_seasons),
  lunar = matrix(yearly_covs$lunar, nrow = n_sites, ncol = n_seasons)
)

umf1 <- unmarkedMultFrame(
  y = y_matrix1,            # 2D detection matrix
  numPrimary = n_seasons,  # Number of primary periods (years)
  siteCovs = site_covs,    # Site-level covariates
  yearlySiteCovs = yearly_covs_list  # Yearly-varying covariates
)

summary(umf1)
# View detection history structure
image(umf1@y, main = "Detection History", 
      xlab = "Site", ylab = "Visit×Year")

# Find Best Detection Model (p)
# Null model (constant detection)
m0 <- colext(psiformula = ~1,      # Initial occupancy
             gammaformula = ~1,    # Colonization
             epsilonformula = ~1,  # Extinction
             pformula = ~1,        # Detection (null)
             data = umf1)

# Weather models
m1 <- colext(~1, ~1, ~1, ~tmax, umf1)               # p ~ tmax
m2 <- colext(~1, ~1, ~1, ~tmin, umf1)        # p ~ tmin
m3 <- colext(~1, ~1, ~1, ~prcp, umf1)        # p ~ precipitation
m4 <- colext(~1, ~1, ~1, ~lunar, umf1)          # p ~ moon phase

# Site structure models
m5 <- colext(~1, ~1, ~1, ~structure, umf1)          # p ~ structure type
m6 <- colext(~1, ~1, ~1, ~canopy, umf1)             # p ~ canopy cover
m7 <- colext(~1, ~1, ~1, ~water, umf1)              # p ~ water presence

# GIS models
m8 <- colext(~1, ~1, ~1, ~gis_tree, umf1)   # p ~ gis tree
m9 <- colext(~1, ~1, ~1, ~gis_water, umf1)      # p ~ gis water 

det_models1 <- fitList(
  "p(.)" = m0,
  "p(tmax)" = m1,
  "p(tmin)" = m2,
  "p(prcp)" = m3,
  "p(lunar)" = m4,
  "p(structure)" = m5,
  "p(canopy)" = m6,
  "p(water)" = m7,
  "p(gis_tree)" = m8,
  "p(gis_water)" = m9#,
#  "p(mean_2019_250)" = m10, 
#  "p(sd_2019_250)" = m11
)

# Rank models by AIC
aic_table <- modSel(det_models1)
print(aic_table)

m6 #top model for myse
#graph
new_p <- expand.grid(
  canopy = factor(levels(umf1@siteCovs$canopy)[1:3],
                  levels = levels(umf1@siteCovs$canopy))
)
# Predict detection p
p_pred <- predict(m6, 
                  type = "det", 
                  newdata = new_p, 
                  se.fit = TRUE)

plot_p <- cbind(new_p, p_pred)


fig_myse <- ggplot(plot_p, aes(x = canopy, y = Predicted)) +
  geom_point(size = 3, fill = "navy") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 1,
                alpha = 0.3,
                color = "steelblue") +
  scale_x_discrete(labels= c( "open" = "Open", 
                              "half" = "Intermediate",
                              "closed"= "Closed")) +
  labs(x = "Site canopy", y = "MYSE detection probabitliy") +
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )

##################myau####################
##only coastal region presence
occu_data1 <- subset(dfo, dfo$region_name=="Middle Atlantic Coastal Plain"|
                       dfo$region_name=="Southeastern Plains")
unique(occu_data$region_name)

#Create detection history and turn into detection matrix#
detection_history <- occu_data1 %>%
  # Select relevant columns (site, year, visit, species)
  select(site, cell, year, night, MYOAUS) 
#replace EPTFUS in species
# epfu EPTFUS, labo LASBOR, laci LASCIN, lano LASNOC, mylu MYOLUC, #
# nyhu NYCHUM, pesu PERSUB, tabr TADBRA#
# mygr MYOGRI, myau MYOAUS, myse MYOSEP

# Convert to detection matrix - SAFE METHOD, not sure if efficient
y_matrix <- detection_history %>%
  pivot_wider(names_from = night, values_from = MYOAUS,
              values_fill = NA, names_sort = TRUE) %>%
  arrange(site, cell, year)

# Create 3D array with proper dimension checking
n_sites <- n_distinct(y_matrix$site)
n_seasons <- n_distinct(y_matrix$year)
n_visits <- 4  # Fixed number of visits

y_array <- array(NA, dim = c(n_sites, n_visits, n_seasons))

for(i in seq_len(n_sites)){
  current_site <- unique(y_matrix$site)[i]
  
  for(t in seq_len(n_seasons)){
    current_year <- unique(y_matrix$year)[t]
    
    site_year_data <- y_matrix %>%
      filter(site == current_site, year == current_year)
    
    if(nrow(site_year_data) == 1){
      # Safely extract visit data
      visit_data <- as.numeric(site_year_data[1, paste0("night", 1:4)])
      
      # Verify length matches expected visits
      if(length(visit_data) == n_visits){
        y_array[i, , t] <- visit_data
      } else {
        warning(paste("Unexpected number of visits for site", current_site, 
                      "year", current_year))
      }
    }
    # Otherwise leaves as NA (for missing surveys)
  }
}

# Verify the array structure
str(y_array)
table(apply(y_array, 3, function(x) sum(!is.na(x)))) # Surveys per year

y_matrix1 <- matrix(y_array, 
                    nrow = dim(y_array)[1], 
                    ncol = dim(y_array)[2] * dim(y_array)[3])


# Site-level covariates (constant across years)
site_covs <- occu_data1 %>%
  distinct(site, structure, canopy, water, mean_2019_250, sd_2019_250, region, wns_sus) %>%
  arrange(site)

# Observation or year level covariates (change across years and nights)
yearly_covs <- occu_data1[,c(1,2,32:34,49,50, 51)]

# Prepare yearlySiteCovs (must be named lists of matrices)
yearly_covs_list <- list(
  tmax = matrix(yearly_covs$tmax, nrow = n_sites, ncol = n_seasons),
  tmin = matrix(yearly_covs$tmin, nrow = n_sites, ncol = n_seasons),
  prcp = matrix(yearly_covs$prcp, nrow = n_sites, ncol = n_seasons),
  gis_tree = matrix(yearly_covs$gis_tree, nrow = n_sites, ncol = n_seasons),
  gis_water = matrix(yearly_covs$gis_water, nrow = n_sites, ncol = n_seasons),
  lunar = matrix(yearly_covs$lunar, nrow = n_sites, ncol = n_seasons)
)

umf1 <- unmarkedMultFrame(
  y = y_matrix1,            # 2D detection matrix
  numPrimary = n_seasons,  # Number of primary periods (years)
  siteCovs = site_covs,    # Site-level covariates
  yearlySiteCovs = yearly_covs_list  # Yearly-varying covariates
)

summary(umf1)
# View detection history structure
image(umf1@y, main = "Detection History", 
      xlab = "Site", ylab = "Visit×Year")

#Find Best Detection Model (p)
# Null model (constant detection)
m0 <- colext(psiformula = ~1,      # Initial occupancy
             gammaformula = ~1,    # Colonization
             epsilonformula = ~1,  # Extinction
             pformula = ~1,        # Detection (null)
             data = umf1)

# Weather models
m1 <- colext(~1, ~1, ~1, ~tmax, umf1)               # p ~ tmax
m2 <- colext(~1, ~1, ~1, ~tmin, umf1)        # p ~ tmin
m3 <- colext(~1, ~1, ~1, ~prcp, umf1)        # p ~ precipitation
m4 <- colext(~1, ~1, ~1, ~lunar, umf1)          # p ~ moon phase

# Site structure models
m5 <- colext(~1, ~1, ~1, ~structure, umf1)          # p ~ structure type
m6 <- colext(~1, ~1, ~1, ~canopy, umf1)             # p ~ canopy cover
m7 <- colext(~1, ~1, ~1, ~water, umf1)              # p ~ water presence

# GIS models
m8 <- colext(~1, ~1, ~1, ~gis_tree, umf1)   # p ~ gis tree
m9 <- colext(~1, ~1, ~1, ~gis_water, umf1)      # p ~ gis water

det_models1 <- fitList(
  "p(.)" = m0,
  "p(tmax)" = m1,
  "p(tmin)" = m2,
  "p(prcp)" = m3,
  "p(lunar)" = m4,
  "p(structure)" = m5,
  "p(canopy)" = m6,
  "p(water)" = m7,
  "p(gis_tree)" = m8,
  "p(gis_water)" = m9
)

# Rank models by AIC
aic_table <- modSel(det_models1)
print(aic_table)

m6 #top model for myau
#graph

new_p <- expand.grid(
  canopy = factor(levels(umf1@siteCovs$canopy)[1:3],
                  levels = levels(umf1@siteCovs$canopy))
)
# Predict detection p
p_pred <- predict(m6, 
                  type = "det", 
                  newdata = new_p, 
                  se.fit = TRUE)

plot_p <- cbind(new_p, p_pred)


fig_myau <- ggplot(plot_p, aes(x = canopy, y = Predicted)) +
  geom_point(size = 3, fill = "navy") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 1,
                alpha = 0.3,
                color = "steelblue") +
  scale_x_discrete(labels= c( "open" = "Open", 
                              "half" = "Intermediate",
                              "closed"= "Closed")) +
  labs(x = "Site canopy", y = "MYAU detection probabitliy") +
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )

##################mygr####################
##only mountain region presence
occu_data1 <- subset(dfo, dfo$region=="66")
unique(occu_data1$region_name)

# Create detection history and turn into detection matrix#
detection_history <- occu_data1 %>%
  # Select relevant columns (site, year, visit, species)
  select(site, cell, year, night, MYOGRI) 
#replace EPTFUS in species
# epfu EPTFUS, labo LASBOR, laci LASCIN, lano LASNOC, mylu MYOLUC, #
# nyhu NYCHUM, pesu PERSUB, tabr TADBRA#
# mygr MYOGRI, myau MYOAUS, myse MYOSEP

# Convert to detection matrix - SAFE METHOD, not sure if efficient
y_matrix <- detection_history %>%
  pivot_wider(names_from = night, values_from = MYOGRI,
              values_fill = NA, names_sort = TRUE) %>%
  arrange(site, cell, year)

# Create 3D array with proper dimension checking
n_sites <- n_distinct(y_matrix$site)
n_seasons <- n_distinct(y_matrix$year)
n_visits <- 4  # Fixed number of visits

y_array <- array(NA, dim = c(n_sites, n_visits, n_seasons))

for(i in seq_len(n_sites)){
  current_site <- unique(y_matrix$site)[i]
  
  for(t in seq_len(n_seasons)){
    current_year <- unique(y_matrix$year)[t]
    
    site_year_data <- y_matrix %>%
      filter(site == current_site, year == current_year)
    
    if(nrow(site_year_data) == 1){
      # Safely extract visit data
      visit_data <- as.numeric(site_year_data[1, paste0("night", 1:4)])
      
      # Verify length matches expected visits
      if(length(visit_data) == n_visits){
        y_array[i, , t] <- visit_data
      } else {
        warning(paste("Unexpected number of visits for site", current_site, 
                      "year", current_year))
      }
    }
    # Otherwise leaves as NA (for missing surveys)
  }
}

# Verify the array structure
str(y_array)
table(apply(y_array, 3, function(x) sum(!is.na(x)))) # Surveys per year

y_matrix1 <- matrix(y_array, 
                    nrow = dim(y_array)[1], 
                    ncol = dim(y_array)[2] * dim(y_array)[3])


# Site-level covariates (constant across years)
site_covs <- occu_data1 %>%
  distinct(site, structure, canopy, water, mean_2019_250, sd_2019_250, region, wns_sus) %>%
  arrange(site)

# Observation or year level covariates (change across years and nights)
yearly_covs <- occu_data1[,c(1,2,32:34,49,50, 51)]

# Prepare yearlySiteCovs (must be named lists of matrices)
yearly_covs_list <- list(
  tmax = matrix(yearly_covs$tmax, nrow = n_sites, ncol = n_seasons),
  tmin = matrix(yearly_covs$tmin, nrow = n_sites, ncol = n_seasons),
  prcp = matrix(yearly_covs$prcp, nrow = n_sites, ncol = n_seasons),
  gis_tree = matrix(yearly_covs$gis_tree, nrow = n_sites, ncol = n_seasons),
  gis_water = matrix(yearly_covs$gis_water, nrow = n_sites, ncol = n_seasons),
  lunar = matrix(yearly_covs$lunar, nrow = n_sites, ncol = n_seasons)
)

umf1 <- unmarkedMultFrame(
  y = y_matrix1,            # 2D detection matrix
  numPrimary = n_seasons,  # Number of primary periods (years)
  siteCovs = site_covs,    # Site-level covariates
  yearlySiteCovs = yearly_covs_list  # Yearly-varying covariates
)

summary(umf1)
# View detection history structure
image(umf1@y, main = "Detection History", 
      xlab = "Site", ylab = "Visit×Year")


# Find Best Detection Model (p)
# Null model (constant detection)
m0 <- colext(psiformula = ~1,      # Initial occupancy
             gammaformula = ~1,    # Colonization
             epsilonformula = ~1,  # Extinction
             pformula = ~1,        # Detection (null)
             data = umf1)

# Weather models
m1 <- colext(~1, ~1, ~1, ~tmax, umf1)               # p ~ tmax
m2 <- colext(~1, ~1, ~1, ~tmin, umf1)        # p ~ tmin
m3 <- colext(~1, ~1, ~1, ~prcp, umf1)        # p ~ precipitation
m4 <- colext(~1, ~1, ~1, ~lunar, umf1)          # p ~ moon phase

# Site structure models
m5 <- colext(~1, ~1, ~1, ~structure, umf1)          # p ~ structure type
m6 <- colext(~1, ~1, ~1, ~canopy, umf1)             # p ~ canopy cover
m7 <- colext(~1, ~1, ~1, ~water, umf1)              # p ~ water presence

# GIS models
m8 <- colext(~1, ~1, ~1, ~gis_tree, umf1)   # p ~ gis tree
m9 <- colext(~1, ~1, ~1, ~gis_water, umf1)      # p ~ gis water

det_models1 <- fitList(
  "p(.)" = m0,
  "p(tmax)" = m1,
  "p(tmin)" = m2,
  "p(prcp)" = m3,
  "p(lunar)" = m4,
  "p(structure)" = m5,
  "p(canopy)" = m6,
  "p(water)" = m7,
  "p(gis_tree)" = m8,
  "p(gis_water)" = m9#,
  #  "p(mean_2019_250)" = m10, 
  #  "p(sd_2019_250)" = m11
)

# Rank models by AIC
aic_table <- modSel(det_models1)
print(aic_table)

m5 #top model for mygr
#graph
new_p <- expand.grid(
  structure = factor(levels(umf1@siteCovs$structure)[1:4],
                     levels = levels(umf1@siteCovs$structure))
)
# Predict detection p
p_pred <- predict(m5, 
                  type = "det", 
                  newdata = new_p, 
                  se.fit = TRUE)

plot_p <- cbind(new_p, p_pred)


fig_mygr <- ggplot(plot_p, aes(x = structure, y = Predicted)) +
  geom_point(size = 3, fill = "navy") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 1,
                alpha = 0.3,
                color = "steelblue") +
  scale_x_discrete(labels= c( "interior" = "Interior", 
                              "corridor" = "Corridor", 
                              "edge" = "Edge",
                              "field" = "Field"), limits = rev) +
  labs(x = "Site structure", y = "MYGR detection probabitliy") +
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )


##################case1 wrap up graph##########
case1 <- ggarrange(fig_epfu, fig_labo, fig_laci, fig_lano, fig_myau, fig_mygr,
                   fig_mylu, fig_myse, fig_nyhu, fig_pesu, fig_tabr,
                   ncol = 4,nrow=3,
                   #        labels = c( "LACI","NYHU", "PESU"),
                   label.x = 0.1,
                   label.y = 0.98,
                   font.label = list(size = 9, color = "black",   family = "serif"),
                   align = "hv",
                   common.legend = FALSE)

ggsave("Figure_case1.tiff",plot=case1, 
       dpi=500, dev='tiff', height=20, width=35, units="cm")

##################Single-season Occupancy model (case 2 for NE)###########
####################load data ###########
Survey_Covariates_Wide <- read.csv("Survey_Covariates_Wide.csv")

Survey_Results <- read.csv("ThreeBat_SurveyResults.csv")

Survey_Results <- Survey_Results %>%
  pivot_wider(
    names_from = Night_Number,
    values_from = c(MYLU, PESU, LACI, NYHU, LANO)
  )

Survey_Covariates_Wide_Standardized <- Survey_Covariates_Wide %>% #Wide format for use in models
  mutate(
    MinTemp_1_Scaled = scale(MinTemp_1),
    MinTemp_2_Scaled = scale(MinTemp_2),
    MinTemp_3_Scaled = scale(MinTemp_3),
    MinTemp_4_Scaled = scale(MinTemp_4),
    Rain_1_Scaled = scale(Rain_1),
    Rain_2_Scaled = scale(Rain_2),
    Rain_3_Scaled = scale(Rain_3),
    Rain_4_Scaled = scale(Rain_4),
    LI_1_Scaled = scale(LI_1),
    LI_2_Scaled = scale(LI_2),
    LI_3_Scaled = scale(LI_3),
    LI_4_Scaled = scale(LI_4), 
    maxTemp_1_Scaled = scale(maxTemp_1),
    maxTemp_2_Scaled = scale(maxTemp_2),
    maxTemp_3_Scaled = scale(maxTemp_3),
    maxTemp_4_Scaled = scale(maxTemp_4)
  )

#fwrite(Data, "NE_2024.csv")

Data <- cbind(Survey_Results, Survey_Covariates_Wide_Standardized)
Data$SiteStructure_1 <- factor(Data$SiteStructure_1,
                               levels=c("Open", "Edge", "Corridor", "Interior"))
Data$SiteStructure_2 <- factor(Data$SiteStructure_2,
                               levels=c("Open", "Edge", "Corridor", "Interior"))
Data$SiteStructure_3 <- factor(Data$SiteStructure_3,
                               levels=c("Open", "Edge", "Corridor", "Interior"))
Data$SiteStructure_4 <- factor(Data$SiteStructure_4,
                               levels=c("Open", "Edge", "Corridor", "Interior"))

####NYHU####
#Creating unmarked frame

NYHU_Frame <-unmarkedFrameOccu(
  y=Data[,c("NYHU_1","NYHU_2","NYHU_3","NYHU_4")],
  obsCovs=list(MinTemp=Data[,c("MinTemp_1","MinTemp_2","MinTemp_3","MinTemp_4")], 
               maxTemp=Data[,c("maxTemp_1","maxTemp_2","maxTemp_3","maxTemp_4")],
               Rain=Data[,c("Rain_1","Rain_2","Rain_3","Rain_4")],
               LI=Data[,c("LI_1","LI_2","LI_3","LI_4")], 
               WaterFacing=Data[,c("WaterFacing_1","WaterFacing_2","WaterFacing_3","WaterFacing_4")], 
               SiteCanopy=Data[,c("SiteCanopy_1", "SiteCanopy_2", "SiteCanopy_3", "SiteCanopy_4")], 
               SiteStructure=Data[,c("SiteStructure_1", "SiteStructure_2","SiteStructure_3", "SiteStructure_4")])
)

#Creating Single Covariate Models

#Null model

NYHU_Null <- occu(~1 ~1, data = NYHU_Frame)

#Site covariate models

NYHU_MinTemp <- occu(~MinTemp ~1, data = NYHU_Frame)
NYHU_maxTemp <- occu(~maxTemp ~1, data = NYHU_Frame)
NYHU_Rain <- occu(~Rain ~1, data = NYHU_Frame)
NYHU_LI <- occu(~LI ~1, data = NYHU_Frame)
NYHU_WaterFacing <- occu(~WaterFacing ~1, data = NYHU_Frame)
NYHU_SiteCanopy <- occu(~SiteCanopy ~1, data = NYHU_Frame)
NYHU_SiteStructure <- occu(~SiteStructure ~1, data = NYHU_Frame)


NYHU_Single_SC_Model_List <- fitList(
  
  "psi, p" = NYHU_Null,
  "psi, p(MinTemp)" = NYHU_MinTemp,
  "psi, p(maxTemp)" = NYHU_maxTemp,
  "psi, p(Rain)" = NYHU_Rain,
  "psi, p(LI)" = NYHU_LI,
  "psi, p(WaterFacing)" = NYHU_WaterFacing,
  "psi, p(SiteCanopy)" = NYHU_SiteCanopy,
  "psi, p(SiteStructure)" = NYHU_SiteStructure
  
)

modSel(NYHU_Single_SC_Model_List)

NYHU_maxTemp 
#graph
new_p <- expand.grid(
  maxTemp = seq(min(NYHU_Frame@obsCovs$maxTemp, na.rm = TRUE),
                  max(NYHU_Frame@obsCovs$maxTemp, na.rm = TRUE),
                  length.out=100)
)
# Predict detection p
p_pred <- predict(NYHU_maxTemp, 
                  type = "det", 
                  newdata = new_p, 
                  se.fit = TRUE)

plot_p <- cbind(new_p, p_pred)

c2_nyhu <- ggplot(plot_p, aes(x = maxTemp, y = Predicted)) +
  geom_line(linewidth = 1, color = "navy") +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.3, 
              fill = "steelblue") +
  labs(x = expression("Daily temperature maximum (" * degree * "C)"), 
       y = "NYHU detection probabitliy") +
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )

####PESU####
#Creating unmarked frame

PESU_Frame <-unmarkedFrameOccu(
  y=Data[,c("PESU_1","PESU_2","PESU_3","PESU_4")],
  obsCovs=list(MinTemp=Data[,c("MinTemp_1","MinTemp_2","MinTemp_3","MinTemp_4")], 
               maxTemp=Data[,c("maxTemp_1","maxTemp_2","maxTemp_3","maxTemp_4")],
               Rain=Data[,c("Rain_1","Rain_2","Rain_3","Rain_4")],
               LI=Data[,c("LI_1","LI_2","LI_3","LI_4")], 
               WaterFacing=Data[,c("WaterFacing_1","WaterFacing_2","WaterFacing_3","WaterFacing_4")], 
               SiteCanopy=Data[,c("SiteCanopy_1", "SiteCanopy_2", "SiteCanopy_3", "SiteCanopy_4")], 
               SiteStructure=Data[,c("SiteStructure_1", "SiteStructure_2","SiteStructure_3", "SiteStructure_4")])
)
PESU_Frame@obsCovs$SiteStructure <- factor(PESU_Frame@obsCovs$SiteStructure,
                                           levels=c("Open", "Edge", "Corridor", "Interior"))

#Creating Single Covariate Models

#Null model

PESU_Null <- occu(~1 ~1, data = PESU_Frame)

#Site covariate models

PESU_MinTemp <- occu(~MinTemp ~1, data = PESU_Frame)
PESU_maxTemp <- occu(~maxTemp ~1, data = PESU_Frame)
PESU_Rain <- occu(~Rain ~1, data = PESU_Frame)
PESU_LI <- occu(~LI ~1, data = PESU_Frame)
PESU_WaterFacing <- occu(~WaterFacing ~1, data = PESU_Frame)
PESU_SiteCanopy <- occu(~SiteCanopy ~1, data = PESU_Frame)
PESU_SiteStructure <- occu(~SiteStructure ~1, data = PESU_Frame)
PESU_m1 <- occu(~SiteStructure + maxTemp ~1, data = PESU_Frame)
PESU_m2 <- occu(~SiteStructure * maxTemp ~1, data = PESU_Frame)

PESU_Single_SC_Model_List <- fitList(
  
  "psi, p" = PESU_Null,
  "psi, p(MinTemp)" = PESU_MinTemp,
  "psi, p(maxTemp)" = PESU_maxTemp,
  "psi, p(Rain)" = PESU_Rain,
  "psi, p(LI)" = PESU_LI,
  "psi, p(WaterFacing)" = PESU_WaterFacing,
  "psi, p(SiteCanopy)" = PESU_SiteCanopy,
  "psi, p(SiteStructure)" = PESU_SiteStructure
)

modSel(PESU_Single_SC_Model_List)

PESU_SiteStructure

#graph
new_data_psi <- expand.grid(
  SiteStructure = unique(Survey_Covariates_Wide$SiteStructure_1)
)

prediction <- predict(PESU_SiteStructure, 
                      type = "det", 
                      newdata = new_data_psi, 
                      se.fit = TRUE)

PESU_Plot_Data_Psi <- cbind(new_data_psi, prediction)

#PESU_Plot_Data_Psi$Species <- "PESU"
PESU_Plot_Data_Psi$SiteStructure <- factor(PESU_Plot_Data_Psi$SiteStructure,
                                           levels=c("Open", "Edge", "Corridor", "Interior"))

c2_PESU <- ggplot(PESU_Plot_Data_Psi, aes(x = SiteStructure, y = Predicted)) +
  geom_point(size = 3, color = "navy")+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 1,
                alpha = 0.3,
                color = "steelblue") +
  xlab("Site Structure") +  
  scale_x_discrete(labels= c( "Open" = "Field",
                              "Edge" = "Edge",
                              "Corridor" = "Corridor", 
                              "Interior" = "Interior" )) +
  ylab( "PESU detection probabitliy")+
  scale_y_continuous(labels = percent_format(accuracy = 1))+
 # ylim(0,1)+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )

####LACI####
#Creating unmarked frame

LACI_Frame <-unmarkedFrameOccu(
  y=Data[,c("LACI_1","LACI_2","LACI_3","LACI_4")],
  obsCovs=list(MinTemp=Data[,c("MinTemp_1","MinTemp_2","MinTemp_3","MinTemp_4")], 
               maxTemp=Data[,c("maxTemp_1","maxTemp_2","maxTemp_3","maxTemp_4")],
               Rain=Data[,c("Rain_1","Rain_2","Rain_3","Rain_4")],
               LI=Data[,c("LI_1","LI_2","LI_3","LI_4")], 
               WaterFacing=Data[,c("WaterFacing_1","WaterFacing_2","WaterFacing_3","WaterFacing_4")], 
               SiteCanopy=Data[,c("SiteCanopy_1", "SiteCanopy_2", "SiteCanopy_3", "SiteCanopy_4")], 
               SiteStructure=Data[,c("SiteStructure_1", "SiteStructure_2","SiteStructure_3", "SiteStructure_4")])
)

LACI_Frame@obsCovs$SiteStructure <- factor(LACI_Frame@obsCovs$SiteStructure,
                                           levels=c("Open", "Edge", "Corridor", "Interior"))

#Creating Single Covariate Models

#Null model

LACI_Null <- occu(~1 ~1, data = LACI_Frame)

#Site covariate models

LACI_MinTemp <- occu(~MinTemp ~1, data = LACI_Frame)
LACI_maxTemp <- occu(~maxTemp ~1, data = LACI_Frame)
LACI_Rain <- occu(~Rain ~1, data = LACI_Frame)
LACI_LI <- occu(~LI~1, data = LACI_Frame)
LACI_WaterFacing <- occu(~WaterFacing ~1, data = LACI_Frame)
LACI_SiteCanopy <- occu(~SiteCanopy ~1, data = LACI_Frame)
LACI_SiteStructure <- occu(~SiteStructure ~1, data = LACI_Frame)

#LACI_m1 <- occu(~MinTemp + SiteStructure ~ 1, data = LACI_Frame)
#LACI_m2 <- occu(~MinTemp * SiteStructure ~ 1, data = LACI_Frame)

LACI_Single_SC_Model_List <- fitList(
  
  "psi, p" = LACI_Null,
  "psi, p(MinTemp)" = LACI_MinTemp,
  "psi, p(maxTemp)" = LACI_maxTemp,
  "psi, p(Rain)" = LACI_Rain,
  "psi, p(LI)" = LACI_LI,
  "psi, p(WaterFacing)" = LACI_WaterFacing,
  "psi, p(SiteCanopy)" = LACI_SiteCanopy,
  "psi, p(SiteStructure)" = LACI_SiteStructure
)

modSel(LACI_Single_SC_Model_List)
LACI_SiteStructure
#graph
new_data_psi <- expand.grid(
  SiteStructure = unique(Survey_Covariates_Wide$SiteStructure_1), 
  MinTemp = median(LACI_Frame@obsCovs$MinTemp)
)

prediction <- predict(LACI_m1, 
                      type = "det", 
                      newdata = new_data_psi, 
                      se.fit = TRUE)

LACI_Plot_Data_Psi <- cbind(new_data_psi, prediction)


LACI_Plot_Data_Psi$SiteStructure <- factor(LACI_Plot_Data_Psi$SiteStructure,
                                           levels=c("Open", "Edge", "Corridor", "Interior"))

c2_LACI <- ggplot(LACI_Plot_Data_Psi, aes(x = SiteStructure, y = Predicted)) +
  geom_point(size = 3, color = "navy")+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 1,
                alpha = 0.3,
                color = "steelblue") +
  xlab("Site structure") +  
  scale_x_discrete(labels= c( "Open" = "Field",
                              "Edge" = "Edge",
                              "Corridor" = "Corridor", 
                              "Interior" = "Interior" )) +
  ylab( "LACI detection probabitliy")+
  scale_y_continuous(labels = percent_format(accuracy = 1))+
  # ylim(0,1)+
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "right",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =12), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks.y = element_blank()
  )

####LANO####
#Creating unmarked frame
LANO_Frame <-unmarkedFrameOccu(
  y=Data[,c("LANO_1","LANO_2","LANO_3","LANO_4")],
  obsCovs=list(MinTemp=Data[,c("MinTemp_1","MinTemp_2","MinTemp_3","MinTemp_4")], 
               maxTemp=Data[,c("maxTemp_1","maxTemp_2","maxTemp_3","maxTemp_4")],
               Rain=Data[,c("Rain_1","Rain_2","Rain_3","Rain_4")],
               LI=Data[,c("LI_1","LI_2","LI_3","LI_4")], 
               WaterFacing=Data[,c("WaterFacing_1","WaterFacing_2","WaterFacing_3","WaterFacing_4")], 
               SiteCanopy=Data[,c("SiteCanopy_1", "SiteCanopy_2", "SiteCanopy_3", "SiteCanopy_4")], 
               SiteStructure=Data[,c("SiteStructure_1", "SiteStructure_2","SiteStructure_3", "SiteStructure_4")])
)

#Creating Single Covariate Models

#Null model

LANO_Null <- occu(~1 ~1, data = LANO_Frame)

#Site covariate models

LANO_MinTemp <- occu(~MinTemp ~1, data = LANO_Frame)
LANO_maxTemp <- occu(~maxTemp ~1, data = LANO_Frame)
LANO_Rain <- occu(~Rain ~1, data = LANO_Frame)
LANO_LI <- occu(~LI~1, data = LANO_Frame)
LANO_WaterFacing <- occu(~WaterFacing ~1, data = LANO_Frame)
LANO_SiteCanopy <- occu(~SiteCanopy ~1, data = LANO_Frame)
LANO_SiteStructure <- occu(~SiteStructure ~1, data = LANO_Frame)


LANO_Single_SC_Model_List <- fitList(
  
  "psi, p" = LANO_Null,
  "psi, p(MinTemp)" = LANO_MinTemp,
  "psi, p(maxTemp)" = LANO_maxTemp,
  "psi, p(Rain)" = LANO_Rain,
  "psi, p(LI)" = LANO_LI,
  "psi, p(WaterFacing)" = LANO_WaterFacing,
  "psi, p(SiteCanopy)" = LANO_SiteCanopy,
  "psi, p(SiteStructure)" = LANO_SiteStructure
  
)

modSel(LANO_Single_SC_Model_List)

####case2 wrap up graph##########
case2 <- ggarrange(c2_LACI, c2_nyhu, c2_PESU,
          ncol = 1,nrow=3,
  #        labels = c( "LACI","NYHU", "PESU"),
          label.x = 0.1,
          label.y = 0.98,
          font.label = list(size = 9, color = "black",   family = "serif"),
          align = "hv",
          common.legend = FALSE)

ggsave("Figure_case2_v2.tiff",plot=case2, 
       dpi=500, dev='tiff', height=25, width=14, units="cm")

##################DeltaAIC graph###############
df <- fread("heatmap.csv")

unique(df$model)
df$model = factor(df$model, 
                  levels = c( "p(.)",  "p(water)", "p(canopy)","p(structure)",
                              "p(tmax)", "p(tmin)", "p(prcp)", "p(lunar)",
                              "p(gis_tree)", "p(gis_water)"))


tcase1 <- ggplot(subset(df, df$case=="nc"), aes(model, species, fill = code)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("null" = "#DEF2F1", 
                               "top" = "#2B7A78", 
                               "yes" = "#3AAFA9",
                               "no" = "#DEF2F1")) +
  scale_x_discrete(labels = c("p(.)" = "Null",
                              "p(canopy)" = "Site canopy", 
                              "p(gis_tree)" = "GIS vegetation",
                              "p(gis_water)" = "GIS water",
                              "p(lunar)" = "Moon phase",
                              "p(prcp)" = "Precipitation",    
                              "p(structure)" = "Site structure", 
                              "p(tmax)" = "Max temperature",   
                              "p(tmin)" = "Min temperature",  
                              "p(water)" = "Water facing" )) +
  labs (x="Model", y="Species") +
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "none",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "white"),
    axis.line.y = element_line(colour = "white"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =10), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks = element_blank()
  )

tcase2 <- ggplot(subset(df, df$case=="ne"), aes(model, species, fill = code)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("null" = "#DEF2F1", 
                               "top" = "#2B7A78", 
                               "yes" = "#3AAFA9",
                               "no" = "#DEF2F1")) +
  scale_x_discrete(labels = c("p(.)" = "Null",
                              "p(canopy)" = "Site canopy", 
                              "p(gis_tree)" = "GIS vegetation",
                              "p(gis_water)" = "GIS water",
                              "p(lunar)" = "Moon phase",
                              "p(prcp)" = "Precipitation",    
                              "p(structure)" = "Site structure", 
                              "p(tmax)" = "Max temperature",   
                              "p(tmin)" = "Min temperature",  
                              "p(water)" = "Water facing" )) +
  labs (x="Model", y="Species") +
  theme(
    #panel.spacing = unit(0, "lines"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12 ), #,face = "bold"
    legend.position = "none",
    legend.title = element_text(family="serif",size=12),
    legend.text = element_text(family="serif",size=12),
    axis.line.x  = element_line(colour = "white"),
    axis.line.y = element_line(colour = "white"),
    axis.title= element_text(family="serif", size=12),
    axis.text = element_text(family="serif", size =12),
    axis.text.x = element_text(angle = 0,family="serif", size =10), #,  hjust=1vjust = 0.5,,
    axis.text.y = element_text(hjust = 0),
    axis.ticks = element_blank()
  )

heatmap <- ggarrange(tcase1, tcase2,
                   ncol = 1,nrow=2,
                   heights = c(10, 4),
                   labels = c( "a", "b"),
                   label.x = 0.005,
                   label.y = 1.01,
                   font.label = list(size = 12, color = "black",   family = "serif"),
                   align = "hv",
                   common.legend = FALSE)
ggsave("Fig_heat.tiff",plot=heatmap, 
       dpi=500, dev='tiff', height=20, width=30, units="cm")      