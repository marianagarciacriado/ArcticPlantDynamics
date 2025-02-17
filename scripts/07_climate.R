## Plant diversity dynamics over space and time in a warming Arctic
## Mariana Garcia Criado (mariana.garcia.criado@gmail.com)
## Script 7. Climate


# Climatic data can be downloaded from https://chelsa-climate.org/


## LIBRARIES ----
library(tidyverse)
library(raster)
library(brms)
library(broom)
library(cowplot)


## THEME ----
bio.theme <- theme(legend.position = "right",
                   axis.title.x = element_text(face="bold", size=20),
                   axis.text.x  = element_text(vjust=0.5, size=18, colour = "black"), 
                   axis.title.y = element_text(face="bold", size=20),
                   axis.text.y  = element_text(vjust=0.5, size=18, colour = "black"),
                   panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
                   panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
                   panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                   plot.title = element_text(color = "black", size = 18, face = "bold", hjust = 0.5),
                   plot.margin = unit(c(1,1,1,1), units = , "cm"))


## FUNCTIONS ----
`%notin%` <- Negate(`%in%`)


  
## CLIMATOLOGIES ----

# CHELSA CLIMATOLOGIES:
# mean annual temp (bio1), temperature of warmest quarter (bio10), 
# temperature of coldest quarter (bio11) and annual precipitation amount (bio12).

# Load climatic data
mat <- raster("C:/wgetdown/CHELSA_bio10_01.tif")
warmq <- raster("C:/wgetdown/CHELSA_bio10_10.tif")
coldq <- raster("C:/wgetdown/CHELSA_bio10_11.tif")
prec <- raster("C:/wgetdown/CHELSA_bio10_12.tif")

# Check that everything makes sense
mat # values: -32768, 32767

# Coordinates
load("data/itex_dec22.RData")
coords0 <- itex.dec22

# Extract climatic data for each pair of coords
coords0$mat <- raster::extract(mat, cbind(coords0$LONG, coords0$LAT))
coords0$mat <- coords0$mat/10
hist(coords0$mat)

coords0$warmq <- raster::extract(warmq, cbind(coords0$LONG, coords0$LAT))
coords0$warmq <- coords0$warmq/10
hist(coords0$warmq)

coords0$coldq <- raster::extract(coldq, cbind(coords0$LONG, coords0$LAT))
coords0$coldq <- coords0$coldq/10
hist(coords0$coldq)

coords0$prec <- raster::extract(prec, cbind(coords0$LONG, coords0$LAT))
hist(coords0$prec)

# Save file with all rows and records
save(coords0, file = "data/22clim_coords.RData")

#load("data/22clim_coords.RData")

# Extract climate data only
clim.only <- coords0 %>% dplyr::select(SiteSubsitePlot, mat, warmq, coldq, prec) %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>% filter(!is.na(mat))

write.csv(clim.only, file = "data/22clim_only.csv")

#clim.only <- read.csv("data/22clim_only.csv")

clim.lat <- coords0 %>% dplyr::select(SiteSubsitePlot, LAT, LONG, mat, warmq, coldq, prec) %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE)



## PLOTTING ----
coords <- coords0 %>% 
  tidyr::pivot_longer(mat:prec, names_to = "clim_var", values_to = "clim_value") %>%
  distinct(LAT, LONG, clim_var, .keep_all = TRUE) %>% 
  dplyr::select(SiteSubsitePlot, LAT, LONG, clim_var, clim_value) %>%
  tidyr::pivot_wider(names_from = clim_var, values_from = clim_value)


# Climatic space
(clim.space <- ggplot(coords) + 
    geom_errorbar(aes(xmin = coldq, xmax = warmq, y = prec), width = 0, size = 0.6, colour = "lightgrey") +
    geom_point(aes(x = mat, y = prec), fill = "#FAA7A1", colour = "black", size = 1.5, shape = 21) + 
    geom_point(aes(x = coldq, y = prec), size = 1.5, fill = "#FDD2CF", colour = "black", shape = 21) + 
    geom_point(aes(x = warmq, y = prec), size = 1.5, fill = "#F8766D", colour = "black", shape = 21) + 
    scale_fill_manual(values = c("#FAA7A1" = "MAT", "#FDD2CF" = "CQT", "#F8766D" = "WQT")) +
    ylab("Annual Precipitation (mm)\n") + xlab("\nTemperature (°C)") + 
   theme(axis.title.x = element_text(size=8),
         axis.text.x  = element_text(vjust=0.5, size=8, colour = "black"), 
         axis.title.y = element_text(size=8),
         axis.text.y  = element_text(vjust=0.5, size=8, colour = "black"),
         panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
         panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
         panel.background = element_blank(), axis.line = element_line(colour = "black"), 
         legend.background = element_rect(fill = NA),
         plot.margin = unit(c(0,0,0,0), units ="cm")))


## CORRELATIONS ----
clim.only <- read.csv("data/22clim_only.csv")

mat.warmq <- brm(mat ~ warmq, data = clim.only, 
                 file = "models/22mat_warmq_mod")
summary(mat.warmq) # positive significant
bayes_R2(mat.warmq) #35.7%

mat.coldq <- brm(mat ~ coldq, 
                 data = clim.only, file = "models/22mat_coldq_mod")
summary(mat.coldq) # positive significant
bayes_R2(mat.coldq) #93.7%

warmq.coldq <- brm(warmq ~ coldq, 
                   data = clim.only, file = "models/22warmq_coldq_mod")
summary(warmq.coldq) # positive significant
bayes_R2(warmq.coldq) #15.7%

warmq.lat <- brm(LAT ~ warmq, 
                 data = clim.lat, file = "models/22warmq_lat_mod")
summary(warmq.lat) # negative significant




## TIMESERIES ----

# CHELSA TIMESERIES:
# mean daily mean air temperatures of the warmest quarter (bio10), 
# annual precipitation amount (bio12)


## Extracting Warm Quarter data per plot

# defining the filepath
folderpath_mat <- "C:/wgetdown/bio10_ts"
filenames_mat <- list.files(folderpath_mat, pattern = "*.tif")
filepath_mat = paste0(folderpath_mat, "/", filenames_mat)

# create raster stack
mat_stack <- stack(filepath_mat)

# add CRS
crs(mat_stack) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# extract the coordinates as a spatial points object 
latlong <- itex.dec22 %>% dplyr::select(LONG, LAT) %>% drop_na() %>% distinct() %>% SpatialPoints()
crs(latlong) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# extract MAT values for each pair of coordinates
master <- mat_stack %>% raster::extract(., latlong, df = TRUE)

# convert the SpatialPoints object into a dataframe 
latlong2 <- as.data.frame(latlong)

# reassign the ID to the latlong and assign it to the mastersheet too
latlong2$ID <- row.names(latlong2)

# merge the two dataframes
master_mat <- merge(master, latlong2, by = c("ID"))

# reshape from wide to long format and convert to celsius
master_mat2 <- master_mat %>% 
  tidyr::pivot_longer(names_to = "warmq_year", values_to = "value", cols = 2:36) %>%
  separate(warmq_year, c("CHELSA", "variable", "year")) %>% dplyr::select(., -CHELSA) %>%
  mutate(value = value/10 - 273.15) 

hist(master_mat2$value)

# save dataframe as a file
write.csv(master_mat2, "data/22warmq_timeseries.csv")





## Extracting Annual Precipitation data per plot

# defining the filepath
folderpath_map <- "C:/wgetdown/bio12_ts"
filenames_map <- list.files(folderpath_map, pattern = "*.tif")
filepath_map = paste0(folderpath_map, "/", filenames_map)

# create raster stack
map_stack <- stack(filepath_map)

# add CRS
crs(map_stack) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# extract MAP values for each pair of coordinates
master.map <- map_stack %>% raster::extract(., latlong, df = TRUE)

# merge the two dataframes
master_map <- merge(master.map, latlong2, by = c("ID"))

# reshape from wide to long format 
master_map2 <- master_map %>% 
  tidyr::pivot_longer(names_to = "prec_year", values_to = "value", cols = 2:36) %>%
  separate(prec_year, c("CHELSA", "variable", "year")) %>% dplyr::select(., -CHELSA)

hist(master_map2$value)

# save dataframe as a file
write.csv(master_map2, "data/22prec_timeseries.csv")



## HIGH LEVEL MODEL ----
prec.overtime <- read.csv("data/22prec_timeseries.csv")
temp.overtime <- read.csv("data/22warmq_timeseries.csv")

# Need subsite info
load("data/itex_dec22.RData")

# keep data at the plot level only
subsite.info <- itex.dec22 %>% dplyr::select(SiteSubsitePlot, SiteSubsite, LAT, LONG) %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE)

# combine with prec and temp data through coordinates
prec.alldata <- left_join(prec.overtime, subsite.info, by = c("LAT", "LONG"))
temp.alldata <- left_join(temp.overtime, subsite.info, by = c("LAT", "LONG"))

hist(prec.alldata$value)
hist(temp.alldata$value)

prec.alldata2 <- prec.alldata %>% mutate(YearCentred = year - 1996)
temp.alldata2 <- temp.alldata %>% mutate(YearCentred = year - 1996)

#mean(prec.alldata$year)
#mean(temp.alldata$year)





## CHANGE OVER TIME ----
temp.overtime$year <- as.numeric(temp.overtime$year)
prec.overtime$year <- as.numeric(prec.overtime$year)

# Calculate slopes (annual change in warmest quarter) for each plot 
warmq.slopes <- temp.overtime %>%
  nest_by(ID) %>%  # replace group_by() with nest_by() to convert your model data to a vector of lists
  mutate(mod = list(lm(value ~ year, data = data))) %>% # change do() to mutate(), then add list() before your model make sure to change data = .  to data = data
  summarise(tidy(mod)) %>%
  filter(term == "year") %>%
  dplyr::select(ID, estimate) %>% 
  rename(WarmQSlope = estimate) %>% unnest()


# Calculate slopes (annual change in precipitation) for each plot 
prec.slopes <- prec.overtime %>%
  nest_by(ID) %>%  # replace group_by() with nest_by() to convert your model data to a vector of lists
  mutate(mod = list(lm(value ~ year, data = data))) %>% # change do() to mutate(), then add list() before your model make sure to change data = .  to data = data
  summarise(tidy(mod)) %>%
  filter(term == "year") %>%
  dplyr::select(ID, estimate) %>% 
  rename(PrecSlope = estimate) %>% unnest()


# Join with plot data
latlong2$ID <- as.numeric(latlong2$ID)

# First with ID and lat/long
id.slopes <- left_join(latlong2, warmq.slopes, by = "ID")
id.slopes2 <- left_join(id.slopes, prec.slopes, by = "ID")

# Now with plot name data
clim.all <- left_join(clim.lat, id.slopes2, by = c("LAT", "LONG"))

# save dataframe as a file
write.csv(clim.all, "data/22clim_slopes.csv")



## DENSITY PLOTS ----
clim.all <- read.csv("data/22clim_slopes.csv")

clim.allX <- clim.all %>% drop_na(WarmQSlope)

# Temperature change over time
(tempch.dens <- ggplot(clim.allX, aes(x=WarmQSlope)) + 
    geom_density(fill="#F8766D", colour="#F8766D", alpha = 0.5, size = 0.7) + 
    geom_vline(aes(xintercept = mean(WarmQSlope)), color="#F8766D", linetype="dotted", size=0.7) + 
    geom_vline(aes(xintercept = 0), color = "black", linetype = "solid", size = 0.1) +
    annotate("text", x = 0.065, y = 30, label = "All plots experienced \ntemperature increases", size = 2) +
    xlab("\nMTWQ change (°C per year)") + ylab("") + bio.theme +
    theme(axis.title.x = element_text(face = "plain", size = 8.5),
          axis.text.x  = element_text(vjust=0.5, size=8, colour = "black"), 
          axis.title.y = element_text(size = 8),
          axis.text.y  = element_text(vjust = 0.5, size=8, colour = "black"),
          plot.margin = unit(c(0,0,0,0), units = "cm")))


precinc <- filter(clim.allX, PrecSlope > 0) # 1905/2174 = 87.6%
precdec <- filter(clim.allX, PrecSlope < 0) # 269/2174 = 12.4%


# Precipitation change over time
(precch.dens <- ggplot(clim.allX, aes(x=PrecSlope)) + 
    geom_density(fill="#00BFC4", colour="#00BFC4", alpha = 0.5, size = 0.7) + 
    geom_vline(aes(xintercept = mean(PrecSlope)), color="#00BFC4", linetype="dotted", size=0.7) +
    geom_vline(aes(xintercept = 0), color = "black", linetype = "solid", size = 0.1) +
    annotate("text", x = 8, y = 0.36, label = "87.6% plots had \nprecipitation \nincreases", size = 2) +
    xlab("\nPrecipitation change (mm per year)") + ylab("") + bio.theme +
    theme(axis.title.x = element_text(face = "plain", size = 8),
          axis.text.x  = element_text(vjust=0.5, size=8, colour = "black"), 
          axis.title.y = element_text(size = 8),
          axis.text.y  = element_text(vjust = 0.5, size=8, colour = "black"),
          plot.margin = unit(c(0,0,0,0), units = , "cm")))


hist(clim.allX$PrecSlope)
hist(clim.allX$WarmQSlope)


# Panel
(climch.panel <- plot_grid(clim.space, tempch.dens, precch.dens, ncol = 3, 
                         align = "hv", labels = c("a", "b", "c"), label_size = 14))
ggplot2::ggsave(climch.panel, filename = "figures/figure_s4.jpg", 
                width = 18, height = 6.5, units = "cm")





## CLIMATIC SPACE ----

# selection of random points
arctic.points <- read.csv("data/GIMMS_Whittaker_data_tundra_10000_random_points_JTKerby.csv")

arctic.points2 <- arctic.points %>% filter(Tundra == "Tundra") %>% 
  filter(FID %notin% c(19, 234, 272)) %>% mutate(prec_mm = MA_prp*10)

# Extract coordinates at subsite level
coords.rich <- coords0 %>% 
  tidyr::pivot_longer(mat:prec, names_to = "clim_var", values_to = "clim_value") %>%
  distinct(SiteSubsite, clim_var, .keep_all = TRUE) %>% 
  dplyr::select(SiteSubsite, LAT, LONG, clim_var, clim_value) %>%
  tidyr::pivot_wider(names_from = clim_var, values_from = clim_value)

# calculate mean plot-level richness at the subsite level
mean.rich <- coords0 %>% dplyr::select(SiteSubsite, SiteSubsitePlot, MeanRichness) %>%
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>% group_by(SiteSubsite) %>% 
  mutate(MeanPlotRichSubs = mean(MeanRichness)) %>% distinct(SiteSubsite, .keep_all = TRUE) %>% 
  dplyr::select(SiteSubsite, MeanPlotRichSubs)

# join them both
coords.mean.rich <- left_join(coords.rich, mean.rich, by = "SiteSubsite")

# plot climate space
(clim.space.all <- ggplot() + 
    geom_point(data = arctic.points2, aes(x = MA_tmp, y = prec_mm), colour = "lightgrey", size = 3) +
    geom_point(data = coords.mean.rich, aes(x = mat, y = prec, fill = MeanPlotRichSubs), 
               colour = "black", size = 7, shape = 21) + 
    ylab("Mean Annual Precipitation (mm)\n") + xlab("\nMean Annual Temperature (°C)") + 
    scale_fill_viridis(option = "plasma", name = "Mean plot-level\nrichness per subsite") +
    theme(axis.title.x = element_text(size=24),
          axis.text.x  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.title.y = element_text(size=24),
          axis.text.y  = element_text(vjust=0.5, size=22, colour = "black"),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.title = element_text(size = 18),
          legend.background = element_rect(fill = NA),
          legend.position = c(0.25, 0.85),
          legend.text = element_text(size=18),
          plot.margin = unit(c(1,1,1,1), units = , "cm")))

ggplot2::ggsave(clim.space.all, filename = "figures/Figure_1b.png", 
                width = 20, height = 20, units = "cm")
