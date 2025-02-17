## Plant diversity dynamics over space and time in a warming Arctic
## Mariana Garcia Criado (mariana.garcia.criado@gmail.com)
## Snow analysis

# Thanks to Noemie Boulanger-Lapointe for data extraction code! 

# Snow data from ARCLIM is available at: 
# https://springernature.figshare.com/collections/ARCLIM_bioclimatic_indices_for_the_terrestrial_Arctic/6216368 



## LIBRARIES ----
library(terra)
library(data.table)
library(tidyverse)
library(brms)


## FUNCTIONS ----
`%notin%` <- Negate(`%in%`)



## ARCLIM DATA ----

# Open points csv file
sites <- read.csv("data/itex_coords_only.csv")

# Convert point file to vector file (projection is the same for the points and raster file)
v <- vect(sites, geom = c("LONG", "LAT"), crs = "epsg:4326")

# read file
fl_terra <- terra::rast('C:/Users/mgarcia5/Desktop/arclim_trends.tif')

fl_terra
levels(fl_terra)

# Extract value at the points for each nc files, note that nc files are opened as raster
site_value <- terra::extract(fl_terra, v) %>%
      mutate(site = sites$SiteSubsite)

# keep only relevant values
site_value_short <- site_value %>%
  dplyr::select(site, arclim_trends_11, arclim_trends_12, arclim_trends_13) %>%
  rename(SnowSeasonLength = arclim_trends_11,
         OnsetSnowSeason = arclim_trends_12,
         EndSnowSeason = arclim_trends_13,
         SiteSubsite = site)

# histogram - make pretty and include in response to reviewers
hist(site_value_short$SnowSeasonLength) # snow season getting shorter overall
hist(site_value_short$OnsetSnowSeason) # snow falling later in the year
hist(site_value_short$EndSnowSeason) # snow ending both ealrier and later

# histogram version
(snow1 <- ggplot(site_value_short, aes(x=SnowSeasonLength)) + 
    geom_histogram(fill="darkgrey", bins = 10) + 
    geom_vline(xintercept = -0.14, linetype = "dashed") +
    xlab("\nChange in snow season length") + ylab("Count\n") +
    bio.theme)


# histogram version
(snow2 <- ggplot(site_value_short, aes(x=OnsetSnowSeason)) + 
    geom_histogram(fill="darkgrey", bins = 10) + 
    geom_vline(xintercept = 0.09, linetype = "dashed") +
    xlab("\nChange in snow season onset") + ylab("Count\n") +
    bio.theme)

# histogram version
(snow3 <- ggplot(site_value_short, aes(x=EndSnowSeason)) + 
    geom_histogram(fill="darkgrey", bins = 10) + 
    geom_vline(xintercept = -0.04, linetype = "dashed") +
    xlab("\nChange in snow season end") + ylab("Count\n") +
    bio.theme)



# mean change
mean(site_value_short$SnowSeasonLength, na.rm = TRUE) # on average getting 0.14 days shorter every year
mean(site_value_short$OnsetSnowSeason, na.rm = TRUE) # snow falling 0.09 days later in the year
mean(site_value_short$EndSnowSeason, na.rm = TRUE) # snow ending 0.04 days earlier in the year


## MODELS ----
rich.time.change <- read.csv("data/rich_time_change.csv")

# join richness change data with snow data
rich.time.change.snow <- rich.time.change %>% left_join(site_value_short, by = "SiteSubsite")

# model with snow variables
rich.time.snow.plot.mod <- brm(Slope ~ SnowSeasonLength + OnsetSnowSeason + EndSnowSeason + 
                                 LogPlotSize + Duration + (1|SiteSubsite), 
                                data = rich.time.change.snow, family = 'gaussian',
                                iter = 3000, chains = 4, warmup = 500, 
                                file = "models/24_03_rich_time_plot_snow_mod")

summary(rich.time.snow.plot.mod) #ns


# turnover data

turnover.change <- read.csv("data/turnover_change.csv")

turnover.change.snow <- turnover.change %>% left_join(site_value_short, by = "SiteSubsite")


# Jaccard snow  model
jac.snow.mod <- brm(bf(Jaccard ~ SnowSeasonLength + OnsetSnowSeason + EndSnowSeason + 
                           LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                      family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                      data = turnover.change.snow, iter = 2000, chains = 4, warmup = 500, 
                      file = "models/24_04_jaccard_snow_mod")

summary(jac.snow.mod) #ns


# Bray-Curtis snow model
bray.snow.mod <- brm(bf(BrayCurtis ~ SnowSeasonLength + OnsetSnowSeason + EndSnowSeason + 
                             LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~1),
                        family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                        data = turnover.change.snow, iter = 2000, chains = 4, warmup = 500, 
                        file = "models/24_04_bray_snow_mod")
summary(bray.snow.mod) #ns


# trajectories data
ext.tur.change <- read.csv("data/ext_tur_change.csv")
col.tur.change <- read.csv("data/col_tur_change.csv")
pers.tur.change <- read.csv("data/pers_tur_change.csv")


ext.snow <- ext.tur.change %>% left_join(site_value_short, by = "SiteSubsite")
col.snow <- col.tur.change %>% left_join(site_value_short, by = "SiteSubsite")
pers.snow <- pers.tur.change %>% left_join(site_value_short, by = "SiteSubsite")



# persisting
pers.snow.mod <- brm(bf(ProportionSpecies ~ SnowSeasonLength + OnsetSnowSeason + EndSnowSeason
                        + LogPlotSize + MeanRichness + Duration + (1|SiteSubsite),
                               zoi ~ 1, coi ~ 1),
                            family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                            data = pers.snow, iter = 2000, chains = 4, warmup = 500, 
                            file = "models/22_04_pers_snow_mod")
summary(pers.snow.mod) #ns


# extinct
ext.snow.mod <- brm(bf(ProportionSpecies ~  SnowSeasonLength + OnsetSnowSeason + EndSnowSeason + 
                         LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                         family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                         data = ext.snow, iter = 2000, chains = 4, warmup = 500, 
                         file = "models/22_04_ext_snow_mod")
summary(ext.snow.mod) #ns


#gains
col.snow.mod <- brm(bf(ProportionSpecies ~ SnowSeasonLength + OnsetSnowSeason + EndSnowSeason + 
                         LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                         family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                         data = col.snow, iter = 2000, chains = 4, warmup = 500, 
                         file = "models/22_04_col_snow_mod")
summary(col.snow.mod) #ns


# evenness
even.time.change2 <- read.csv("data/even_time_change2.csv")

even.snow <- even.time.change2 %>% left_join(site_value_short, by = "SiteSubsite")

even.time.snow.mod <- brm(Slope ~ SnowSeasonLength + OnsetSnowSeason + EndSnowSeason + 
                                  MeanRichness + Duration + (1|SiteSubsite), 
                                data = even.snow, family = 'gaussian',
                                iter = 2000, chains = 4, warmup = 400, 
                                file = "models/24_05_even_snow_mod")
summary(even.time.snow.mod) #ns
