## Plant diversity dynamics over space and time in a warming Arctic
## Mariana Garcia Criado (mariana.garcia.criado@gmail.com)
## Script 2. Richness analyses


## PACKAGES ----
library(tidyverse)
library(brms)
library(fitdistrplus)
library(viridis)
library(ggOceanMaps)
library(broom)
library(ggrepel)
library(cowplot)
library(ggeffects)
library(ggnewscale)
library(ggdist)


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



## DATA ----
load("data/itex_dec22.RData")
richness.end.plot <- read.csv("data/end_itex_rich22.csv")




## SPATIAL MODS ----

# What's the mean richness per plot? (Intercept-only)
rich.intonly.mod <- brm(CurrentRichness ~ 1, data = richness.end.plot, 
                        iter = 2000, chains = 4, warmup = 400, 
                        file = "models/22_00_rich_int_mod")
summary(rich.intonly.mod) # mean richness per plot is 7.21


# Centre latitude and divide cover by 100
richness.end.plot2 <- richness.end.plot %>% filter(!is.na(LAT)) %>% 
  filter(!is.na(SurveyedArea)) %>%
  mutate(LatitudeCentred = LAT - mean(LAT)) %>% 
  mutate(LogPlotSize = log(SurveyedArea)) %>%
  mutate(ShrubCentred = PlotShrubMean - mean(PlotShrubMean), 
         ForbCentred = PlotForbMean - mean(PlotForbMean), 
         GramCentred = PlotGraminoidMean - mean(PlotGraminoidMean), 
         ShrubDivided = PlotShrubMean/100, 
         ForbDivided = PlotForbMean/100, 
         GramDivided = PlotGraminoidMean/100)



#### 1) GEOGRAPHICAL MODEL ####
rich.geo.mod <- brm(CurrentRichness ~ LAT + LogPlotSize + Region + (1|SiteSubsite),
                    data = richness.end.plot2, family = 'negbinomial', 
                    prior = prior(gamma(0.01, 0.01), class = shape), 
                    iter = 2500, chains = 4, warmup = 400, 
                    file = "models/22_01_rich_geo_mod")

print(summary(rich.geo.mod), digits = 5) # latitude negative significant, plot size ns, region ns
conditional_effects(rich.geo.mod)
bayes_R2(rich.geo.mod) #0.67 [conditional]
bayes_R2(rich.geo.mod, re.form = NA) #0.09 [marginal]

summary(rich.geo.mod, prob = 0.975) # latitude negative significant, plot size ns, region ns
conditional_effects(rich.geo.mod, prob = 0.975)

# Model predictions - get number of species per degree from here
rich.geo.pred <- ggpredict(rich.geo.mod, terms = "LAT [60:82, sample = 30]")
colnames(rich.geo.pred) = c('Latitude', 'fit', 'lwr', 'upr', 'dunno')


## version of this model with 1m2 plots only
richness.end.1m <- richness.end.plot2 %>% filter(SurveyedArea == 1)

rich.geo.mod.1m <- brm(CurrentRichness ~ LAT + Region + (1|SiteSubsite),
                    data = richness.end.1m, family = 'negbinomial', 
                    prior = prior(gamma(0.01, 0.01), class = shape), 
                    iter = 2500, chains = 4, warmup = 400, 
                    file = "models/22_01b_rich_geo_mod_1m2")

summary(rich.geo.mod.1m) # lat positive ns, region ns
conditional_effects(rich.geo.mod.1m)


#### 2) CLIMATE MODEL ####

# file created from climate script
clim.only <- read.csv("data/22clim_only.csv")

# Combine data
richness.end.plot3 <- left_join(richness.end.plot2, clim.only, by = "SiteSubsitePlot")

# remove empty values
richness.end.plot4 <- richness.end.plot3 %>% 
  filter(!is.na(warmq)) %>% 
  filter(MOISTURE != "")


# Fit model
rich.clim.mod <- brm(CurrentRichness ~ warmq + prec + MOISTURE + LogPlotSize + (1|SiteSubsite), 
                     data = richness.end.plot4, family = 'negbinomial', 
                     prior = prior(gamma(0.01, 0.01), class = shape), 
                     iter = 3000, chains = 4, warmup = 500,
                     file = "models/22_02_rich_clim_mod")
summary(rich.clim.mod) # positive sig warm q, prec ns, moisture ns, plot size ns
print(summary(rich.clim.mod), digits = 5)
plot(rich.clim.mod)
conditional_effects(rich.clim.mod)
bayes_R2(rich.clim.mod) #0.63
bayes_R2(rich.clim.mod, re.form = NA) #0.16

summary(rich.clim.mod, prob = 0.975) # temp positive sign, prec ns, plot size ns, moisture ns
conditional_effects(rich.clim.mod, prob = 0.975)

# get number of species per degree temp - 1sp gain every 2deg
rich.clim.pred <- ggpredict(rich.clim.mod, terms = "warmq [all]")
colnames(rich.clim.pred) = c('temp', 'fit', 'lwr', 'upr', 'dunno')



#### 3) FG-DOM MODEL ####

# Model using mean FG proportion with all FGs doesn't converge given that they are all negatively correlated

# breaking it up into 3 models
rich.fg.shb <- brm(CurrentRichness ~ ShrubDivided + LogPlotSize + (1|SiteSubsite),
                 data = richness.end.plot2, family = 'negbinomial', 
                 prior = prior(gamma(0.01, 0.01), class = shape), 
                 iter = 2000, chains = 4, warmup = 500, 
                 file = "models/22_03a_rich_fg_mod")

summary(rich.fg.shb) # shrub ns, plot size positive significant
summary(rich.fg.shb, prob = 0.975) # shrub ns, plot size positive significant

# forb
rich.fg.frb <- brm(CurrentRichness ~ ForbDivided + LogPlotSize + (1|SiteSubsite),
                   data = richness.end.plot2, family = 'negbinomial', 
                   prior = prior(gamma(0.01, 0.01), class = shape), 
                   iter = 2000, chains = 4, warmup = 500, 
                   file = "models/22_03b_rich_fg_mod")

summary(rich.fg.frb) # forb positive significant, plot size positive significant
summary(rich.fg.frb, prob = 0.975) # forb positive significant, plot size ns

# gram
rich.fg.gram <- brm(CurrentRichness ~ GramDivided + LogPlotSize + (1|SiteSubsite),
                   data = richness.end.plot2, family = 'negbinomial', 
                   prior = prior(gamma(0.01, 0.01), class = shape), 
                   iter = 2000, chains = 4, warmup = 500, 
                   file = "models/22_03c_rich_fg_mod")

summary(rich.fg.gram) # gram negative significant, plot size positive significant
summary(rich.fg.gram, prob = 0.975) # gram negative significant, plot size positive significant




## SPATIAL MAP ----

# if these packages are not installed, they can be installed directly from github like so:
#devtools::install_github("MikkoVihtakari/ggOceanMapsData")
#devtools::install_github("MikkoVihtakari/ggOceanMaps")


# Richness at the site level
itex.end <- read.csv("data/end_only_itex22.csv")

# Calculate mean richness per plot
rich.site <- itex.end %>% distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>%
  group_by(SITE) %>% mutate(MeanSiteRichness = mean(AnnualRichness)) %>% ungroup() %>% 
  distinct(SITE, .keep_all = TRUE) %>% dplyr::select(SITE, LAT, LONG, MeanSiteRichness)

# Convert coords to utm
richsite.coords <- transform_coord(x = rich.site, lon = "LONG", lat = "LAT", 
                                   new.names = c("lon.utm", "lat.utm"), 
                                   verbose = FALSE, bind = TRUE)


# Figure 1: map with richness per site (with labels)
(richness.site.map <- basemap(limits = 57, land.col = "#c9c7c1") +
    geom_point(data = richsite.coords, aes(x = lon.utm, y = lat.utm, colour = MeanSiteRichness),
               size = 15, alpha = 0.7) + scale_colour_viridis(option = "plasma") +
    labs(colour='Mean plot-level\nrichness per study area') +
    geom_label_repel(data = subset(richsite.coords, SITE %in% c("ALEXFIORD", "ABISKO", "DISKO",
                                                                "THINGVELLIR", "TOOLIK", "ENDALEN", "DARING")),
                     aes(lon.utm, lat.utm, label = SITE), color = "black", box.padding = 2,
                     segment.color = "black", segment.size = 1.5, fill = "white", label.size = 1,
                     size=10, max.iter = 50000) +
    theme(legend.title = element_text(size = 30), 
          legend.text=element_text(size = 30),
          legend.position="top", 
          legend.margin=margin(1,1,1,1),
          legend.box.margin=margin(10,10,10,10), 
          plot.margin=grid::unit(c(0,0,0,0), "mm")))

ggsave(richness.site.map, filename = "figures/Figure_1a.png",
       width = 30, height = 30, units = "cm")



# ED Figure 6b: Plot map at site level - unique species per site (without labels)
(richness.site.map2 <- basemap(limits = 57, land.col = "#c9c7c1") + 
    geom_point(data = richsite.coords, aes(x = lon.utm, y = lat.utm, colour = MeanSiteRichness), 
               size = 5, alpha = 0.7) + 
    scale_colour_viridis(option = "plasma") +
    labs(colour='Mean plot-level\nrichness per study area') +
    theme(legend.title = element_text(size = 8), 
          legend.text=element_text(size = 8), 
          legend.position="top", 
          legend.margin=margin(1,1,1,1),
          legend.key.height = unit(0.2, 'cm'),
          legend.box.margin=margin(1,1,1,11), 
          plot.margin=grid::unit(c(0,0,0,0), "mm")))


# ED figure 6a: Scatterplot with richness across latitude (each point is a plot)

# Colour by richness at site level
rich.site.only <- rich.site %>% dplyr::select(SITE, MeanSiteRichness)
rich.all.vars <- left_join(richness.end.plot, rich.site.only, by = "SITE")

(rich.lat.plot <- ggplot() + 
    geom_point(data = rich.all.vars, aes(x = LAT, y = CurrentRichness, colour = MeanSiteRichness), alpha = 0.7, size = 2) + 
    scale_colour_viridis(option = "plasma") +
    xlab("\nLatitude (°)") + ylab("Current richness\n") +
    geom_line(data = rich.geo.pred, aes(x = Latitude, y = fit), colour = "black", size = 1) + 
    geom_ribbon(data = rich.geo.pred, aes(x = Latitude, ymax = upr, ymin = lwr), fill = "black", alpha = 0.2) + 
    theme(axis.title.x = element_text(face="plain", size=8),
          axis.text.x  = element_text(vjust=0.5, size=8, colour = "black"), 
          axis.title.y = element_text(face="plain", size=8),
          axis.text.y  = element_text(vjust=0.5, size=8, colour = "black"),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", legend.key = element_blank()))

# ED figure 6 panel
(C2 <- plot_grid(rich.lat.plot, NULL, richness.site.map2, 
                 ncol = 3, rel_widths = c(0.5, 0.01, 1),  
                 align = "h", labels = c("a", "", "b"), label_size = 14))

ggsave(C2, filename = "figures/figure_s6.jpg", 
       width = 18, height = 9.5, units = "cm")







## TEMPORAL MODS ----

# We need annual values of richness to fit the models
load("data/itex_dec22.RData")
itex.dec22$SurveyedArea <- as.numeric(itex.dec22$SurveyedArea)

# Keep a row per plotXyear only for those plots that were surveyed more than once
rich.time <- itex.dec22 %>% distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  group_by(SiteSubsitePlot) %>% mutate(NumYears = length(unique(YEAR))) %>% 
  ungroup() %>% filter(NumYears > 1) %>% group_by(SiteSubsitePlot) %>%
  mutate(Duration = max(YEAR) - min(YEAR)) %>% filter(Duration > 4) %>% ungroup() %>% 
  mutate(LogPlotSize = log(SurveyedArea)) %>% 
  mutate(YearCentred = YEAR - 2003)

length(unique(rich.time$SiteSubsitePlot))
length(unique(rich.time$SITE))



#### 0) HIGH-LEVEL MODEL WITH RANDOM SLOPES ####
hist(rich.time$AnnualRichness) # count data starting at 1

# variation is greater than the mean, so neg-binomial distribution
mean(rich.time$AnnualRichness) # 7.24
var(rich.time$AnnualRichness) # 8.71

# Negative binomial distribution and nested random effect
big.rich.time.nest.mod <- brm(AnnualRichness ~ YearCentred + (YearCentred|SiteSubsite/SiteSubsitePlot), 
                         family = "negbinomial", data = rich.time, init = 0,
                         iter = 3000, chains = 4, warmup = 400, 
                         file = "models/24_00_rich_year_ranslopnest_mod")
print(summary(big.rich.time.nest.mod), digits = 4) # 0.0021, CI  = -0.0002 to 0.0043
plot(big.rich.time.nest.mod)
conditional_effects(big.rich.time.nest.mod)
bayes_R2(big.rich.time.nest.mod) # 0.63
bayes_R2(big.rich.time.nest.mod, re.form = NA) # 0.003

# NB. this model is not uploaded to GitHub as its size is greater than the permitted size limit.

# predictions
big.rich.pred <- ggpredict(big.rich.time.nest.mod, terms = "YearCentred [sample = 30]",
                           type = "random", condition = c(SiteSubsite = 1266), 
                           allow_new_levels = TRUE)
colnames(big.rich.pred) = c('year', 'fit', 'lwr', 'upr', 'group')
# 0.01 species increase per year, 0.16 species per decade (from predictions)



  ## SUBSITE LEVEL SLOPES ##

# Extract subsite-level slopes to fit univariate models
site.richch.slopes <- as.data.frame(ranef(big.rich.time.nest.mod)$SiteSubsite[,,2], optional = F)

site.richch.slopes$SiteSubsite <- row.names(site.richch.slopes)


# Put together latitude and climate info
change_clim <- read.csv("data/22clim_slopes.csv")

# merge by plot
lat.warm <- left_join(rich.time, change_clim, by = "SiteSubsitePlot")

# keep subsite level data
lat.warm2 <- lat.warm %>% distinct(SiteSubsite, .keep_all = TRUE) %>% 
  dplyr::select(SITE, SiteSubsite, WarmQSlope, LAT.x, warmq)

# join richness chane slopes with temp change slopes
site.richch.slopes.latwarm <- left_join(site.richch.slopes, lat.warm2, by = "SiteSubsite")

# Gaussian distribution
hist(site.richch.slopes.latwarm$Estimate)

# Model richness change at the subsite level
site.rich.mod <- brm(Estimate ~ WarmQSlope, data = site.richch.slopes.latwarm, 
                     iter = 3000, chains = 4, warmup = 500, 
                     file = "models/24_richch_tempch_subsite_mod")

print(summary(site.rich.mod), digits = 4) # positive ns
conditional_effects(site.rich.mod)

fix.site.rich.mod <- as.data.frame(fixef(site.rich.mod))

# predictions 
site.rich.mod.pred <- ggpredict(site.rich.mod, terms = "WarmQSlope [all]")
colnames(site.rich.mod.pred) = c('WarmQSlope', 'fit', 'lwr', 'upr')

# rich change vs temp change (at subsite level)
(site.slopes.warm <- ggplot() +
    geom_point(data = site.richch.slopes.latwarm, aes(x = WarmQSlope, y = Estimate, fill = warmq, colour = warmq), size = 6, shape = 21, alpha = 0.8) + 
    geom_line(data = site.rich.mod.pred, aes(x = WarmQSlope, y = fit), size = 2, linetype = "dashed", colour = "black") +
    scale_colour_gradient(low = "#66ccff", high = "#ff6666", guide = "none") +
    scale_fill_gradient(low = "#66ccff", high = "#ff6666", name = "MTWQ\n(°C)") +
    geom_ribbon(data = site.rich.mod.pred, aes(x = WarmQSlope, ymin = lwr, ymax = upr), fill = "black", alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    xlab("\nChange in Warmest Quarter Temperature\n(°C per year)") +
    ylab("Richness change at subsite level\n(log(species) per year)\n") + 
    bio.theme +
    theme(axis.title.x = element_text(face="plain", size=24), 
          axis.title.y = element_text(face="plain", size=24),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20)))

ggsave(site.slopes.warm, filename = "figures/Figure_2d.png", 
       width = 24, height = 18, units = "cm")




#### 1) RICHNESS CHANGE OVER TIME ####

# Calculate slopes (annual richness change) for each plot 
rich.time.slopes <- rich.time %>%
  nest_by(SiteSubsitePlot) %>%  # replace group_by() with nest_by() to convert your model data to a vector of lists
  mutate(mod = list(lm(AnnualRichness ~ YEAR, data = data))) %>% # change do() to mutate(), then add list() before your model make sure to change data = .  to data = data
  summarise(tidy(mod)) %>%
  filter(term == "YEAR") %>%
  dplyr::select(SiteSubsitePlot, estimate) %>% 
  rename(Slope = estimate)

write.csv(rich.time.slopes, "data/22rich_time_slopes.csv")
rich.time.slopes <- read.csv("data/22rich_time_slopes.csv")

# Richness change over time
hist(rich.time.slopes$Slope)

# What's the mean richness change across plots? (Intercept-only)
richtime.intonly.mod <- brm(Slope ~ 1, data = rich.time.slopes, 
                            iter = 2000, chains = 4, warmup = 400, 
                            file = "models/24_00_richtime_int_mod")
print(summary(richtime.intonly.mod), digits = 5) # mean slope of richness change is 0.002 [-0.007 to 0.001]
plot(richtime.intonly.mod)

# extract model results for plot below
mean <- 0.002
low.ci <- -0.0002
high.ci <- 0.0043
df.eb <- data.frame(mean, low.ci, high.ci)

# Plot slopes per plot
(richslopes <- ggplot() + 
    geom_histogram(data = rich.time.slopes, aes(x=Slope, y = ..count..), breaks = seq(-1, 0.7, by = 0.1), fill="lightgray", alpha = 0.7) +
    geom_density(data = rich.time.slopes, aes(x=Slope, y = ..density..*(1259*0.1)), fill="#75C2D0", colour="#75C2D0", alpha = 0.4, size = 1.2) + 
    geom_vline(aes(xintercept=0.002), color="#75C2D0", linetype="dashed", size=1.5) + 
    annotate("text", x = -0.58, y = 420, label = "Mean richness change:\n0.002 log(species) per year\n", size = 7.5) +
    xlab("\nRichness change at plot level\n(log(species) per year)") + ylab("Count\n") + bio.theme +
    theme(axis.title.x = element_text(face="plain", size=24), 
          axis.title.y = element_text(face="plain", size=24)))

ggplot2::ggsave(richslopes, filename = "figures/figure_2c.png", 
                width = 25, height = 20, units = "cm")





#### 2) SUBSITE LEVEL MODEL ####

# Now we need one row per plot only
rich.time.linexplot0 <- rich.time %>% distinct(SiteSubsitePlot, .keep_all = TRUE)

# Join with slope data
rich.time.linexplot <- left_join(rich.time.linexplot0, rich.time.slopes, by = "SiteSubsitePlot")

write.csv(rich.time.linexplot, "data/22rich_time_all.csv")
rich.time.linexplot <- read.csv("data/22rich_time_all.csv")

# Post-hoc subsite-level model 
rich.time.subsite.mod <- brm(Slope ~ LAT + Region + LogPlotSize + Duration + (1|SiteSubsite), 
                             data = rich.time.linexplot, family = 'gaussian',
                             iter = 2000, chains = 4, warmup = 400, 
                             file = "models/24_04_rich_time_subs_mod")
print(summary(rich.time.subsite.mod), digits = 5) # latitude ns, plot size ns, duration ns
conditional_effects(rich.time.subsite.mod) # region ns

print(summary(rich.time.subsite.mod, prob = 0.975), digits = 5) # latitude ns, plot size ns, duration ns
conditional_effects(rich.time.subsite.mod, prob = 0.975) # region ns


# Make predictions for figure below
rich.time.pred <- ggpredict(rich.time.subsite.mod, terms = "LAT")
colnames(rich.time.pred) = c('Latitude', 'fit', 'lwr', 'upr', 'dunno')
glimpse(rich.time.pred)




#### 3) PLOT LEVEL MODEL WITH CHANGE (CLIMATE + FUNCTIONAL GROUP) ####
clim.only <- read.csv("data/22clim_only.csv")

# Combine data
rich.time.linexplot2 <- left_join(rich.time.linexplot, clim.only, by = "SiteSubsitePlot")

# Remove empty values
rich.time.linexplot3 <- rich.time.linexplot2 %>% filter(!is.na(warmq)) %>% 
  filter(!is.na(SurveyedArea)) %>%
  filter(MOISTURE != "") 

# Add change data
change_clim <- read.csv("data/22clim_slopes.csv")
change_fg <- read.csv("data/22fg_slopes.csv")

rich.time.change0 <- left_join(rich.time.linexplot3, change_clim, by = "SiteSubsitePlot")
rich.time.change <- left_join(rich.time.change0, change_fg, by = "SiteSubsitePlot")

# Save dataset with only change 
only.chg <- rich.time.change %>% 
  dplyr::select(Slope, MOISTURE, WarmQSlope, PrecSlope,
                ShrubSlope, LogPlotSize, Duration, SiteSubsite, SiteSubsitePlot)

write.csv(only.chg, "data/select_covs.csv")

# plot-level model - one group at a time to facilitate convergence
rich.time.shbch.plot.mod <- brm(Slope ~ MOISTURE + WarmQSlope + PrecSlope +
                               ShrubSlope + LogPlotSize + Duration + (1|SiteSubsite), 
                             data = rich.time.change, family = 'gaussian',
                             iter = 3000, chains = 4, warmup = 500, 
                             file = "models/24_03_rich_time_plot_shbch_mod")

print(summary(rich.time.shbch.plot.mod), digits = 5) # temp ns, prec ns, shrub negative, plot size ns, duration positive
summary(rich.time.shbch.plot.mod, prob = 0.975) # temp ns, prec ns, shrub negative, plot size ns, duration ns
conditional_effects(rich.time.shbch.plot.mod) # moisture ns
ranef(rich.time.shbch.plot.mod)
bayes_R2(rich.time.shbch.plot.mod) #0.16

bayes_R2(rich.time.shbch.plot.mod, re.form = NULL) #same as above, this is indeed the default [conditional, includes REfs]
bayes_R2(rich.time.shbch.plot.mod, re.form = NA) #marginal, only random effects = 0.05



## version of the model without outliers
rich.time.chg.outliers <- rich.time.change %>% drop_na(ShrubSlope) %>% 
  mutate(ShrubChangeSD = sd(ShrubSlope), ShrubChangeSD2 = ShrubChangeSD*2, 
         ShrubChangeSD3 = ShrubChangeSD*3, ShrubChangeSD3Neg = -(ShrubChangeSD3)) %>% 
  filter(ShrubSlope < ShrubChangeSD3) %>% # this gets rid of the positive outliers but not the negative ones
  filter(ShrubSlope > ShrubChangeSD3Neg)


# model without outliers
rich.time.shbch.plot.outliers.mod <- brm(Slope ~ MOISTURE + WarmQSlope + PrecSlope +
                                           ShrubSlope + LogPlotSize + Duration + (1|SiteSubsite), 
                                         data = rich.time.chg.outliers, family = 'gaussian',
                                         iter = 3500, chains = 2, warmup = 500, 
                                         file = "models/24_03_rich_time_plot_shbch_outliers_mod")

summary(rich.time.shbch.plot.outliers.mod) # temp ns, prec ns, shrub negative, plot size ns, duration positive
summary(rich.time.shbch.plot.outliers.mod, prob = 0.975) # temp ns, prec ns, shrub negative, plot size ns, duration ns
conditional_effects(rich.time.shbch.plot.outliers.mod) # moisture ns
bayes_R2(rich.time.shbch.plot.outliers.mod)
bayes_R2(rich.time.shbch.plot.outliers.mod, re.form = NA)

# model predictions
shrubchg.outliers.pred <- ggpredict(rich.time.shbch.plot.outliers.mod, terms = "ShrubSlope", 
                                    type = "random", 
                                    condition = c(SiteSubsite = 1266),
                                    allow_new_levels = TRUE)
colnames(shrubchg.outliers.pred) = c('ShrubSlope', 'fit', 'lwr', 'upr', 'group')


# plot without outliers
(shrubchg.out.plot <- ggplot(rich.time.chg.outliers) + 
    geom_point(aes(x = ShrubSlope, y = Slope), size = 2, colour = "#5fbf33", alpha = 0.5) +
    geom_line(data = shrubchg.outliers.pred, aes(x = ShrubSlope, y = fit), colour = "#215028", linewidth = 1) + 
    geom_ribbon(data = shrubchg.outliers.pred, aes(x = ShrubSlope, ymin = lwr, ymax = upr), fill = "#215028", alpha = 0.2) +
    xlab("\nShrub cover change\n(% per year)") +
    ylab("Richness change\n(species/year)\n") + 
    bio.theme + 
    theme(axis.title.x = element_text(face="plain", size=9), 
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y  = element_text(vjust=0.5, size=9, colour = "black"), 
          legend.key.size = unit(1, 'cm'), 
          plot.margin = unit(c(0.2,0.2,0.2,0.2), units = , "cm")))

## now only with plot size = 1m2
rich.time.change.1m <- rich.time.change %>% filter(SurveyedArea == 1)

# Post-hoc plot-level model - one group at a time to facilitate convergence
rich.time.shbch.plot.1m.mod <- brm(Slope ~ MOISTURE + WarmQSlope + PrecSlope +
                                  ShrubSlope + Duration + (1|SiteSubsite), 
                                data = rich.time.change.1m, family = 'gaussian',
                                iter = 3000, chains = 4, warmup = 500, 
                                file = "models/24_03c_rich_time_plot_shbch_1m_mod")
print(summary(rich.time.shbch.plot.1m.mod), digits = 5) # temp ns, prec ns, shrub negative, duration ns
summary(rich.time.shbch.plot.mod, prob = 0.975) # temp ns, prec ns, shrub negative, duration ns
conditional_effects(rich.time.shbch.plot.1m.mod) # moisture ns

# Forb change model
rich.time.forbch.plot.mod <- brm(Slope ~ MOISTURE + WarmQSlope + PrecSlope +
                                  ForbSlope + LogPlotSize + Duration + (1|SiteSubsite), 
                                data = rich.time.change, family = 'gaussian',
                                iter = 3500, chains = 2, warmup = 500, 
                                file = "models/24_03_rich_time_plot_frbch_mod")

summary(rich.time.forbch.plot.mod) # temp ns, prec ns, forb positive, plot size ns, duration positve
summary(rich.time.forbch.plot.mod, prob = 0.975) # temp ns, prec ns, forb positive, plot size ns, duration ns
conditional_effects(rich.time.forbch.plot.mod) # moisture ns
bayes_R2(rich.time.forbch.plot.mod) #0.18
bayes_R2(rich.time.forbch.plot.mod, re.form = NA) #0.07

# predictions 
forbchg.pred <- ggpredict(rich.time.forbch.plot.mod, terms = "ForbSlope", 
                          type = "random", condition = c(SiteSubsite = 1266), 
                          allow_new_levels = TRUE)
colnames(forbchg.pred) = c('ForbSlope', 'fit', 'lwr', 'upr', 'group')


# plot
(forbchg.plot <- ggplot(rich.time.change) + 
    geom_point(aes(x = ForbSlope, y = Slope, colour = PlotForbMean, fill = PlotForbMean), size = 6, alpha = 0.8) +
    geom_line(data = forbchg.pred, aes(x = ForbSlope, y = fit), colour = "#490439", linewidth = 1.2) + 
    geom_ribbon(data = forbchg.pred, aes(x = ForbSlope, ymin = lwr, ymax = upr), fill = "#490439", alpha = 0.2) +
    scale_colour_gradient(low = "#DDAFD3", high = "#490439", guide = "none") +
    scale_fill_gradient(low = "#DDAFD3",  high = "#490439", name = "Mean Forb\nCover (%)") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    xlab("\nForb cover change\n(% per year)") +
    ylab("Richness change\n(species/year)\n") + 
    bio.theme + 
    theme(axis.title.x = element_text(face="plain", size=24), 
          axis.text.x  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.title.y = element_text(face="plain", size=24),
          axis.text.y  = element_text(vjust=0.5, size=22, colour = "black"),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          legend.key.size = unit(1, 'cm')))

ggsave(forbchg.plot, filename = "scripts/mgarciacriado/figures/Figure_2f.png", 
       width = 24, height = 18, units = "cm")



## now without outliers 
rich.time.chg.outliers2 <- rich.time.change %>% drop_na(ForbSlope) %>% 
  mutate(ForbChangeSD = sd(ForbSlope), ForbChangeSD2 = ForbChangeSD*2, 
         ForbChangeSD3 = ForbChangeSD*3, ForbChangeSD3Neg = -(ForbChangeSD3)) %>% 
  filter(ForbSlope < ForbChangeSD3) %>% # this gets rid of the positive outliers but not the negative ones
  filter(ForbSlope > ForbChangeSD3Neg)


# Forb change model without outliers
rich.time.forbch.plot.out.mod <- brm(Slope ~ MOISTURE + WarmQSlope + PrecSlope +
                                   ForbSlope + LogPlotSize + Duration + (1|SiteSubsite), 
                                 data = rich.time.chg.outliers2, family = 'gaussian',
                                 iter = 3500, chains = 2, warmup = 500, 
                                 file = "models/24_03_rich_time_plot_frbch_outliers_mod")

summary(rich.time.forbch.plot.out.mod) # temp ns, prec ns, forb positive, plot size ns, duration positve
summary(rich.time.forbch.plot.out.mod, prob = 0.975) # temp ns, prec ns, forb positive, plot size ns, duration positive
conditional_effects(rich.time.forbch.plot.out.mod) # moisture ns
bayes_R2(rich.time.forbch.plot.out.mod)
bayes_R2(rich.time.forbch.plot.out.mod, re.form = NA)


# predictions 
forbchg.out.pred <- ggpredict(rich.time.forbch.plot.out.mod, terms = "ForbSlope", 
                              type = "random", condition = c(SiteSubsite = 1266), 
                              allow_new_levels = TRUE)
colnames(forbchg.out.pred) = c('ForbSlope', 'fit', 'lwr', 'upr', 'group')


# plot
(forbchg.out.plot <- ggplot(rich.time.chg.outliers2) + 
    geom_point(aes(x = ForbSlope, y = Slope), size = 2, colour = "#95558b", alpha = 0.5) +
    geom_line(data = forbchg.out.pred, aes(x = ForbSlope, y = fit), colour = "#490439", linewidth = 1) + 
    geom_ribbon(data = forbchg.out.pred, aes(x = ForbSlope, ymin = lwr, ymax = upr), fill = "#490439", alpha = 0.2) +
    xlab("\nForb cover change\n(% per year)") +
    ylab("Richness change\n(species/year)\n") + 
    bio.theme + 
    theme(axis.title.x = element_text(face="plain", size=9), 
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y  = element_text(vjust=0.5, size=9, colour = "black"), 
          legend.key.size = unit(1, 'cm'),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), units = , "cm")))



# Graminoid change model
rich.time.gramch.plot.mod <- brm(Slope ~ MOISTURE + WarmQSlope + PrecSlope +
                                   GraminoidSlope + LogPlotSize + Duration + (1|SiteSubsite), 
                                 data = rich.time.change, family = 'gaussian',
                                 iter = 2000, chains = 4, warmup = 400, 
                                 file = "models/24_03_rich_time_plot_gramch_mod")

summary(rich.time.gramch.plot.mod) # prec ns, temp ns, gram ns, plot size ns, duration positive
summary(rich.time.gramch.plot.mod, prob = 0.975) # prec ns, temp ns, gram ns, plot size ns, duration positive
conditional_effects(rich.time.gramch.plot.mod) # moisture ns



# predictions 
gramchg.pred <- ggpredict(rich.time.gramch.plot.mod, terms = "GraminoidSlope")
colnames(gramchg.pred) = c('GraminoidSlope', 'fit', 'lwr', 'upr')


# plot
(gramchg.plot <- ggplot(rich.time.change) + 
    geom_point(aes(x = GraminoidSlope, y = Slope), size = 5, colour = "#dc9537", alpha = 0.5) +
    geom_line(data = gramchg.pred, aes(x = GraminoidSlope, y = fit), colour = "#ca6115", linewidth = 1.2) + 
    geom_ribbon(data = gramchg.pred, aes(x = GraminoidSlope, ymin = lwr, ymax = upr), fill = "#ca6115", alpha = 0.2) +
    xlab("\nGraminoid cover change\n(% per year, slopes)") +
    ylab("Richness change\n(species/year, slopes)\n") + 
    bio.theme + 
    theme(axis.title.x = element_text(face="bold", size=20), 
          axis.text.x  = element_text(vjust=0.5, size=18, colour = "black"), 
          axis.title.y = element_text(face="bold", size=20),
          axis.text.y  = element_text(vjust=0.5, size=18, colour = "black"), 
          legend.key.size = unit(1, 'cm')))





## now without outliers 
rich.time.chg.outliers3 <- rich.time.change %>% drop_na(GraminoidSlope) %>% 
  mutate(GramChangeSD = sd(GraminoidSlope), GramChangeSD2 = GramChangeSD*2, 
         GramChangeSD3 = GramChangeSD*3, GramChangeSD3Neg = -(GramChangeSD3)) %>% 
  filter(GraminoidSlope < GramChangeSD3) %>% # this gets rid of the positive outliers but not the negative ones
  filter(GraminoidSlope > GramChangeSD3Neg)


# Graminoid change model
rich.time.gramch.plot.out.mod <- brm(Slope ~ MOISTURE + WarmQSlope + PrecSlope +
                                   GraminoidSlope + LogPlotSize + Duration + (1|SiteSubsite), 
                                 data = rich.time.chg.outliers3, family = 'gaussian',
                                 iter = 2000, chains = 4, warmup = 400, 
                                 file = "models/24_03_rich_time_plot_gramch_out_mod")

print(summary(rich.time.gramch.plot.out.mod), digits = 5) # prec ns, temp ns, gram ns, plot size ns, duration positive
summary(rich.time.gramch.plot.out.mod, prob = 0.975) # prec ns, temp ns, gram ns, plot size ns, duration positive
conditional_effects(rich.time.gramch.plot.out.mod) # moisture ns
bayes_R2(rich.time.gramch.plot.out.mod)
bayes_R2(rich.time.gramch.plot.out.mod, re.form = NA)

# predictions 
gramchg.out.pred <- ggpredict(rich.time.gramch.plot.out.mod, terms = "GraminoidSlope", 
                              type = "random", condition = c(SiteSubsite = 1266), 
                              allow_new_levels = TRUE)
colnames(gramchg.out.pred) = c('GraminoidSlope', 'fit', 'lwr', 'upr', 'group')

# plot
(gramchg.out.plot <- ggplot(rich.time.chg.outliers3) + 
    geom_point(aes(x = GraminoidSlope, y = Slope), size = 2, colour = "#dc9537", alpha = 0.5) +
    geom_line(data = gramchg.out.pred, aes(x = GraminoidSlope, y = fit), colour = "#ca6115", linewidth = 1, linetype = "dashed") + 
    geom_ribbon(data = gramchg.out.pred, aes(x = GraminoidSlope, ymin = lwr, ymax = upr), fill = "#ca6115", alpha = 0.2) +
    xlab("\nGraminoid cover change\n(% per year)") +
    ylab("Richness change\n(species/year)\n") + 
    bio.theme + 
    theme(axis.title.x = element_text(face="plain", size=9), 
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y  = element_text(vjust=0.5, size=9, colour = "black"), 
          legend.key.size = unit(1, 'cm'),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), units = , "cm")))


# ED Fig 7b-d: Panel without outliers
(rchg.out.panel <- plot_grid(shrubchg.out.plot, forbchg.out.plot, gramchg.out.plot, ncol = 3, 
                         align = "hv", labels = c("b", "c", "d"), label_size = 14))


## Plot-level model with temp x temp change interaction
rich.temp.int.mod <- brm(Slope ~ WarmQSlope*warmq.x + LogPlotSize + Duration + (1|SiteSubsite), 
                                data = rich.time.change, family = 'gaussian',
                                iter = 3500, chains = 2, warmup = 500, 
                                file = "models/24_rich_temp_int_mod")
summary(rich.temp.int.mod)
conditional_effects(rich.temp.int.mod) # interaction ns
bayes_R2(rich.temp.int.mod)
bayes_R2(rich.temp.int.mod, re.form = NA)





## SHRUB ARROW PLOT ----

# Arrow plot for shrub change
rcc <- rich.time %>% group_by(SiteSubsitePlot) %>% 
  filter(YEAR %in% c(min(YEAR), max(YEAR))) %>% 
  arrange(SiteSubsitePlotYear) %>% mutate(StartEnd = ifelse(YEAR == min(YEAR), "Start", "End")) %>%
  mutate(StartShrubCover = case_when(StartEnd == "Start" ~ ShrubCover)) %>% 
  mutate(EndShrubCover = case_when(StartEnd == "End" ~ ShrubCover)) %>% 
  mutate(StartRichness = case_when(StartEnd == "Start" ~ AnnualRichness)) %>%
  mutate(EndRichness = case_when(StartEnd == "End" ~ AnnualRichness)) %>% 
  summarise_each(funs(first(.[!is.na(.)]))) %>% ungroup()

# Add types of relationships
rcc2 <- rcc %>% mutate(RelType = case_when(EndShrubCover > StartShrubCover & EndRichness > StartRichness ~ "Rel1", 
                                           EndShrubCover > StartShrubCover & EndRichness < StartRichness ~ "Rel2", 
                                           EndShrubCover < StartShrubCover ~ "Rel34", TRUE ~ "Rel56"))

rel1 <- rcc2 %>% filter(RelType == "Rel1") # 186/1266 = 14.7%
rel2 <- rcc2 %>% filter(RelType == "Rel2") # 220/1266 = 17.3%
rel3 <- rcc2 %>% filter(RelType == "Rel34") # 556/1266 = 43.9%
rel5 <- rcc2 %>% filter(RelType == "Rel56") # 304/1266 = 24%

# Keep only the increases in shrub as we're interested in this question
relxxx <- rcc2 %>% filter(RelType %in% c("Rel1", "Rel2"))

# Calculate shrub cover change per plot
rich.shb.slopes <- rich.time %>%
  nest_by(SiteSubsitePlot) %>%  
  mutate(mod = list(lm(ShrubCover ~ YEAR, data = data))) %>% 
  summarise(tidy(mod)) %>%
  filter(term == "YEAR") %>%
  dplyr::select(SiteSubsitePlot, estimate, p.value) %>% 
  rename(ShrubSlope = estimate)

## Model with shrub increases (including both significant and non-significant)

# filter only those plots where shrub has increased 
rich.shb.slopes.filterXX <- rich.shb.slopes %>% filter(ShrubSlope > 0)

# extract name of plots as a vector
plots.vectorXX <- rich.shb.slopes.filterXX$SiteSubsitePlot

# extract the info for these sites from the main database
selected.plotsXX <- rich.time %>% filter(SiteSubsitePlot %in% plots.vectorXX) %>% 
  group_by(SiteSubsitePlot) %>% filter(YEAR == min(YEAR)) %>% mutate(StartShrubCover = ShrubCover) %>% 
  dplyr::select(SiteSubsite, SiteSubsitePlot, StartShrubCover, Duration)

# extract the richness slopes
rich.time.slopes <- read.csv("data/22rich_time_slopes.csv")

rich.slopes.filterXX <- rich.time.slopes %>% filter(SiteSubsitePlot %in% plots.vectorXX) %>% 
  dplyr::select(SiteSubsitePlot, Slope) %>% rename(RichnessSlope = Slope)

# join datasets together
arrow.final.dbXX <- left_join(rich.shb.slopes.filterXX, selected.plotsXX, by = "SiteSubsitePlot")
arrow.final.db2XX <- left_join(arrow.final.dbXX, rich.slopes.filterXX, by = "SiteSubsitePlot")

# Model with start shrub cover
arrow.new.modXX <- brm(RichnessSlope ~ ShrubSlope*StartShrubCover + (1|SiteSubsite), 
                     data = arrow.final.db2XX, iter = 3500, chains = 2, warmup = 500, 
                     file = "models/24_arrow_new_modXX")

print(summary(arrow.new.modXX), digits = 4) # all ns
plot(conditional_effects(arrow.new.modXX), points = TRUE)

# Model with duration interaction
arrow.new.mod <- brm(RichnessSlope ~ ShrubSlope*Duration + (1|SiteSubsite), 
                     data = arrow.final.db2XX, iter = 3500, chains = 2, warmup = 500, 
                     file = "models/24_arrow_new_mod_zzzz")
print(summary(arrow.new.mod), digits = 3) # all ns
plot(conditional_effects(arrow.new.mod), points = TRUE)


# Arrow plot (shrub increases only)
orangeandblue <- left_join(relxxx, rich.shb.slopes, by = "SiteSubsitePlot")

(arrow.plot.onlyincreases <- ggplot(orangeandblue) +
    geom_segment(aes(x = StartShrubCover, y = StartRichness, xend = EndShrubCover, yend = EndRichness, 
                     size = ShrubSlope, colour = RelType), arrow = arrow(length = unit(0.10, "inches"))) + 
    scale_colour_manual(values = c("#5F90E2", "#E59434"), name = "Trend", 
                        labels = c("Positive", "Negative")) +
    scale_size_continuous(range = c(0.01, 1.2), name = "Shrub change\n(%/year)") + 
    xlab("\nShrub cover (%)") +
    ylab("Plot richness (number)\n") + bio.theme +
    guides(size = guide_legend(order = 2), colour = guide_legend(order = 1)) +
    theme(legend.title = element_text(face = "bold", size = 9),
          legend.text = element_text(size = 9), 
          legend.key = element_blank(),
          axis.title.x = element_text(face="plain", size=9),
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y = element_text(vjust=0.5, size=9),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), units = "cm")))

fig.s7a <- plot_grid(arrow.plot.onlyincreases, labels = c("a"), label_size = 14)
            
# ED Fig 7: together with arrow plot (7a, created below)
(fig.s7 <- plot_grid(fig.s7a, rchg.out.panel, ncol = 1, nrow = 2, rel_heights = c(1.5, 1),
                     align = "hv"))
ggplot2::ggsave(fig.s7, filename = "figures/figure_s7.jpg", 
                width = 18, height = 14, units = "cm")


## RICHNESS VS RICHNESS CHANGE ----
rich.richch.mod <- brm(Slope ~ MeanRichness, data = rich.time.linexplot, 
                       iter = 2000, chains = 4, warmup = 400, 
                       file = "models/24_richchange_vs_rich_mod")
print(summary(rich.richch.mod), digits = 5) # ns
bayes_R2(rich.richch.mod)



## RICHNESS CHANGE MAP ----

# Calculate mean richness change per site
rich.ch.site <- rich.time.linexplot %>% group_by(SITE) %>% 
  mutate(MeanRichnessChange = mean(Slope)) %>%
  ungroup() %>% distinct(SITE, .keep_all = TRUE) %>% 
  mutate(AbsoluteRichnessChange = abs(MeanRichnessChange)) %>% 
  arrange(-AbsoluteRichnessChange)

# Convert coords to utm
rich.ch.site.coords <- transform_coord(x = rich.ch.site, lon = "LONG", lat = "LAT", 
                                   new.names = c("lon.utm", "lat.utm"), 
                                   verbose = FALSE, bind = TRUE)

# Figure 2a - map of richness change at site level 
(rich.ch.site.map <- basemap(limits = 57, land.col = "#c9c7c1") + 
    geom_point(data = rich.ch.site.coords, 
               aes(x = lon.utm, y = lat.utm, 
                   colour = MeanRichnessChange, size = AbsoluteRichnessChange), 
               alpha = 0.8) + 
    scale_color_viridis(option = "plasma") +
    scale_radius(range=c(10, 50)) + 
    labs(colour = "Mean plot richness\nchange per study area\n(species/year)", size = "") +
    guides(colour = guide_colorbar(direction = "horizontal", barwidth = 25, 
                                   title.position = "top", title.hjust = 0.5)) +
    theme(legend.title = element_text(size = 32, hjust = 0.5), 
          legend.text = element_text(size = 32, hjust = 0.5)))

ggsave(rich.ch.site.map, filename = "figures/Figure_2a.png", 
       width = 50, height = 50, units = "cm")

# Figure 2b
(plot.richslopes.time <- ggplot() +
    geom_point(data = rich.time, aes(x = YEAR, y = AnnualRichness, fill = LAT, colour = LAT), size = 4, shape = 21, alpha = 0.2) + 
    scale_colour_viridis(option = "plasma", guide = "none", direction = -1) +
    scale_fill_viridis(option = "plasma", name = "Latitude (°)", direction = -1) +
    geom_line(data = big.rich.pred, aes(x = year + 2003, y = fit), colour = "black", linewidth = 2, linetype = "dashed") +
    geom_ribbon(data = big.rich.pred, aes(x = year + 2003, ymin = lwr, ymax = upr), fill = "black", alpha = 0.1) +
    xlab("\nYear") +
    ylab("Richness\n(number of species)\n") + 
    bio.theme +
    theme(axis.title.x = element_text(face = "plain", size = 24), 
          axis.title.y = element_text(face = "plain", size = 24),
          legend.position = "top", 
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20)))

ggplot2::ggsave(plot.richslopes.time, filename = "figures/Figure_2b.png", width = 22, height = 18, units = "cm")





## POSTHOC ANALYSES ----


#### SHRUB INCREASE ~ LATITUDE ####
shrubslopeinc <- rich.time.change %>% filter(ShrubSlope > 0) %>% 
  mutate(ShrubSlopeDivided = ShrubSlope/100)

# model
shrubinc.mod <- brm(ShrubSlopeDivided ~ LAT.y + (1|SiteSubsite), data = shrubslopeinc, 
                    family = "beta", iter = 3500, chains = 2, warmup = 500, 
                    file = "models/24_shbinc_lat_mod")

summary(shrubinc.mod)

# predictions
shrubinc.pred <- ggpredict(shrubinc.mod, terms = "LAT.y", type = "random",
                           condition = c(SiteSubsite = 1266), allow_new_levels = TRUE)
colnames(shrubinc.pred) = c('Latitude', 'fit', 'lwr', 'upr', 'group')

# figure
(shrubchg.plot.col <- ggplot(shrubslopeinc) + 
    geom_point(aes(y = ShrubSlopeDivided*100, x = LAT.y, colour = WarmQSlope), size = 1.5, alpha = 0.5) +
    scale_colour_viridis(option = "inferno", name = "Temperature\nchange (°C)") +
    geom_line(data = shrubinc.pred, aes(x = Latitude, y = fit*100), colour = "black", linewidth = 1, linetype = "dashed") + 
    geom_ribbon(data = shrubinc.pred, aes(x = Latitude, ymin = lwr*100, ymax = upr*100), fill = "black", alpha = 0.2) +
    ylab("Shrub cover change\n(% per year)\n") +
    xlab("\nLatitude (°)") + 
    bio.theme + 
    theme(axis.title.x = element_text(face="plain", size=9), 
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y  = element_text(vjust=0.5, size=9, colour = "black"), 
          legend.background=element_blank(),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 5),
          legend.position = c(0.2,0.84),
          legend.key.size = unit(0.2, 'cm'),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "lines")))




#### SHRUB CATEGORIES ####

# Add in shrub functional group category from Garcia Criado et al. (2023) traits paper
load("data/2022_master_new_superfinal.RData")

# extract categories 
fg.only <- master.new.superfinal %>% dplyr::select(sp, gf) %>% 
  distinct(sp, .keep_all = T) %>% rename(SPECIES_NAME = sp)

# add column in ITEX data
itex.fg <- left_join(itex.dec22, fg.only, by = "SPECIES_NAME")

# Check which ones don't have a FG
itex.fg2 <- itex.fg %>% filter(FuncGroup == "Shrub") %>% filter(is.na(gf))
unique(itex.fg2$SPECIES_NAME) # species without func group data

itex.fg$gf[itex.fg$SPECIES_NAME == "Ledum palustre"] <- "Low shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Salix ovalifolia"] <- "Dwarf shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Arctostaphylos rubra"] <- "Low shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Salix alaxensis"] <- "Tall shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Empetrum hermaphroditum"] <- "Low shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Juniperus communis"] <- "Tall shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Arctostaphylos alpina"] <- "Low shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Cassiope hypnoides"] <- "Dwarf shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Betula pubescens"] <- "Tall shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Salix arctophila"] <- "Low shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Vaccinium microcarpum"] <- "Dwarf shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Ledum decumbens"] <- "Low shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Salix chamissonis"] <- "Dwarf shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Salix fuscescens"] <- "Tall shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Oxycoccus microcarpus"] <- "Dwarf shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Salix planifolia"] <- "Tall shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Salix arbusculoides"] <- "Tall shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Salix lanata x reticulata"] <- "Low shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Dryas integrifolia x octopetala"] <- "Dwarf shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "XXXDryas:ZACKENBERG"] <- "Dwarf shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "XXXDryas:TOOLIK"] <- "Dwarf shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "XXXWOODYDRYAS:SADVENT"] <- "Dwarf shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "Salix arctica/arctophila"] <- "Low shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "XXXSALIX:JAMESONLAND"] <- "Low shrub"
itex.fg$gf[itex.fg$SPECIES_NAME == "XXXWOODYSALIX:SADVENT"] <- "Dwarf shrub"
# Not sure about "XXXSalix:TOOLIK" as it could be a tall (most likely) or dwarf shrub species
# Not sure about "XXXSALIX:AUYUITTUQ" as it could be a dwarf (most likely) or tall shrub species


# check that all spps have a FG now
itex.fg3 <- itex.fg %>% filter(FuncGroup == "Shrub") %>% filter(is.na(gf))
unique(itex.fg3$SPECIES_NAME)

# transform into dwarf and not-dwarf
itex.fg.shb <- itex.fg %>% 
  filter(FuncGroup == "Shrub") %>%
  mutate(gf = case_when(gf %in% c("Low shrub", "Tall shrub") ~ "Non-dwarf",
                        gf == "Dwarf shrub" ~ "Dwarf"))

check <- itex.fg.shb %>% filter(is.na(gf)) 
unique(check$SPECIES_NAME)

# Calculate cover values per plotXyear for dwarf vs non-dwarf cover
fg.values <- itex.fg.shb %>% group_by(SiteSubsitePlot) %>%
  filter(length(unique(YEAR)) > 2 & (max(YEAR) - min(YEAR)) > 4) %>% ungroup %>%
  group_by(SiteSubsitePlotYear, gf) %>% 
  mutate(DwarfCover = case_when(gf == "Dwarf" ~ sum(RelCover)),
         NonDwarfCover = case_when(gf == "Non-dwarf" ~ sum(RelCover))) %>%
  ungroup() %>% 
  distinct(SiteSubsitePlotYear, DwarfCover, NonDwarfCover, .keep_all = TRUE) %>% 
  dplyr::select(SiteSubsite, SiteSubsitePlot, SiteSubsitePlotYear, YEAR, DwarfCover, NonDwarfCover) %>%
  group_by(SiteSubsitePlotYear) %>% 
  summarise_each(funs(first(.[!is.na(.)]))) %>% 
  ungroup() 

# Add values of NA for those functional groups that were never present
fg.values2 <- fg.values %>% replace(is.na(.), 0) %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(DwarfShrubTotal = sum(DwarfCover)) %>%
  mutate(NonDwarfShrubTotal = sum(NonDwarfCover)) %>%
  ungroup() %>%
  mutate(DwarfShrubCoverNew = ifelse(DwarfShrubTotal == 0, NA, DwarfCover)) %>% 
  mutate(NonDwarfShrubCoverNew = ifelse(NonDwarfShrubTotal == 0, NA, NonDwarfCover)) %>% 
  dplyr::select(SiteSubsite, SiteSubsitePlot, SiteSubsitePlotYear, YEAR, 
                DwarfShrubCoverNew, NonDwarfShrubCoverNew)

# Convert to long format (with NA)
fg.values.long <- fg.values2 %>% 
  pivot_longer(cols = c(DwarfShrubCoverNew, NonDwarfShrubCoverNew), 
               names_to = "FuncGroup", values_to = "CoverFG")

# Long format without NA
fg.values.long.nona <- fg.values.long %>% filter(!is.na(CoverFG))


# Calculate slopes (annual change in functional group) for each plot 
fg.dwarf.slopes <- fg.values.long.nona %>%
  nest_by(SiteSubsitePlot, FuncGroup) %>%  # replace group_by() with nest_by() to convert your model data to a vector of lists
  mutate(mod = list(lm(CoverFG ~ YEAR, data = data))) %>% # change do() to mutate(), then add list() before your model make sure to change data = .  to data = data
  summarise(tidy(mod)) %>%
  filter(term == "YEAR") %>%
  dplyr::select(SiteSubsitePlot, FuncGroup, estimate) %>% 
  rename(FGSlope = estimate) %>% unnest()

write.csv(fg.dwarf.slopes, "data/dwarf_slopes.csv")

#fg.dwarf.slopes <- read.csv("data/dwarf_slopes.csv")

# Convert to wide format, keep those which weren't there at the start and at the end as NA
fg.dwarf.slopes.wide <- fg.dwarf.slopes %>% 
  tidyr::pivot_wider(names_from = FuncGroup, values_from = FGSlope) %>% 
  rename(DwarfShrubSlope = DwarfShrubCoverNew, NonDwarfShrubSlope = NonDwarfShrubCoverNew)

# And subsite
subs <- itex.dec22 %>% dplyr::select(SiteSubsite, SiteSubsitePlot) %>% 
  distinct(SiteSubsitePlot, .keep_all = T)

fg.dwarf.slopes.wide2 <- left_join(fg.dwarf.slopes, subs, by = "SiteSubsitePlot")



#### RICHNESS CHANGE VS SHRUB CHANGE WITH CATEGORIES ####

# Does the relationship between shrub cover change and richness change depend on FG?

# add in richness change slopes
rich.time.slopes <- read.csv("data/22rich_time_slopes.csv")
all.slopes.fg <- left_join(fg.dwarf.slopes.wide2, rich.time.slopes, by = "SiteSubsitePlot")


# model
rich.shb.fg <- brm(Slope ~ FGSlope*FuncGroup + (1|SiteSubsite), data = all.slopes.fg, 
                   family = 'gaussian', iter = 3500, chains = 2, warmup = 500, 
                   file = "models/24_rich_shb_fg_mod")

print(summary(rich.shb.fg), digits = 4)
conditional_effects(rich.shb.fg) # significant interaction
bayes_R2(rich.shb.fg) #0.08
bayes_R2(rich.shb.fg, re.form = NA)

# model predcitions (FG model)
dwarf.pred <- ggpredict(rich.shb.fg, 
                        terms = c("FGSlope [-15:10, sample = 30]", "FuncGroup"),
                        type = "random", condition = c(SiteSubsite = 1266),
                        allow_new_levels = TRUE)
colnames(dwarf.pred) = c('ShrubSlope', 'fit', 'lwr', 'upr', "FuncGroup")


# plot
(shrubchg.plot <- ggplot() + 
    geom_point(data = rich.time.change, aes(x = ShrubSlope, y = Slope, colour = PlotShrubMean, fill = PlotShrubMean), size = 5, alpha = 0.5) +
    scale_colour_gradient(low = "#F0E68C", high = "#28a428", guide = "none") +
    scale_fill_gradient(low = "#F0E68C", high = "#28a428", name = "Mean Shrub\nCover (%)") +
    geom_line(aes(x = ShrubSlope, y = fit), colour = "#8b4539", linewidth = 2, linetype = "dashed",
              data = dwarf.pred %>% filter(FuncGroup == "DwarfShrubCoverNew")) + 
    geom_ribbon(aes(x = ShrubSlope, ymin = lwr, ymax = upr), fill = "#8b4539", alpha = 0.2,
                data = dwarf.pred %>%
                  filter(FuncGroup == "DwarfShrubCoverNew")) +
    geom_line(aes(x = ShrubSlope, y = fit), colour = "#215028", linewidth = 2, linetype = "solid",
              data = dwarf.pred %>% filter(FuncGroup == "NonDwarfShrubCoverNew")) + 
    geom_ribbon(aes(x = ShrubSlope, ymin = lwr, ymax = upr), fill = "#215028", alpha = 0.2,
                data = dwarf.pred %>%
                  filter(FuncGroup == "NonDwarfShrubCoverNew")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    xlab("\nShrub cover change\n(% per year)") +
    ylab("Richness change\n(species/year)\n") + 
    bio.theme + 
    theme(axis.title.x = element_text(face="plain", size=24), 
          axis.text.x  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.title.y = element_text(face="plain", size=24),
          axis.text.y  = element_text(vjust=0.5, size=22, colour = "black"), 
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          legend.key.size = unit(1, 'cm')))

ggsave("figures/Figure_2e.png", shrubchg.plot, width = 24, height = 18, units = "cm")





#### SHRUB CHANGE ~ TEMPERATURE CHANGE ####

# Are rates of temperature change related to rates of shrub cover change?
fg.dwarf.slopes.temp <- left_join(fg.dwarf.slopes, change_clim, by = "SiteSubsitePlot")


# Model with functional group interaction
dwarf.mod <- brm(FGSlope ~ FuncGroup*WarmQSlope + (1|SiteSubsite), 
                 data = fg.dwarf.slopes.temp2, 
                family = 'gaussian', iter = 3500, chains = 2, warmup = 500, 
                file = "models/23_dwarf_shbchg_temp_mod")
summary(dwarf.mod) #ns
conditional_effects(dwarf.mod)


# model predictions 
dwarf.temp.pred <- ggpredict(dwarf.mod, terms = c("WarmQSlope", "FuncGroup"), type = "random",
                             condition = c(SiteSubsite = 1266), allow_new_levels = TRUE)
colnames(dwarf.temp.pred) = c('WarmQSlope', 'fit', 'lwr', 'upr', "FuncGroup")


# plot
(dwarf.temp.plot <- ggplot() + 
    geom_point(data = rich.time.change, aes(y = ShrubSlope, x = WarmQSlope, colour = PlotShrubMean, fill = PlotShrubMean), size = 1.5, alpha = 0.5) +
    scale_colour_gradient(low = "#F0E68C", high = "#28a428", guide = "none") +
    scale_fill_gradient(low = "#F0E68C", high = "#28a428", name = "Mean Shrub\nCover (%)") +
    geom_line(aes(x = WarmQSlope, y = fit), colour = "#8b4539", linewidth = 1, linetype = "dashed",
              data = dwarf.temp.pred %>% filter(FuncGroup == "DwarfShrubCoverNew")) + 
    geom_ribbon(aes(x = WarmQSlope, ymin = lwr, ymax = upr), fill = "#8b4539", alpha = 0.2,
                data = dwarf.temp.pred %>%
                  filter(FuncGroup == "DwarfShrubCoverNew")) +
    geom_line(aes(x = WarmQSlope, y = fit), colour = "#215028", linewidth = 1, linetype = "dashed",
              data = dwarf.temp.pred %>% filter(FuncGroup == "NonDwarfShrubCoverNew")) + 
    geom_ribbon(aes(x = WarmQSlope, ymin = lwr, ymax = upr), fill = "#215028", alpha = 0.2,
                data = dwarf.temp.pred %>%
                  filter(FuncGroup == "NonDwarfShrubCoverNew")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    xlab("\nTemperature change\n(°C/year)") +
    ylab("Shrub cover change\n(%/year)\n") + 
    bio.theme + 
    theme(axis.title.x = element_text(face="plain", size=9), 
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y  = element_text(vjust=0.5, size=9, colour = "black"), 
          legend.background=element_blank(),
          legend.key.size = unit(0.2, 'cm'),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 6),
          legend.direction = "vertical",
          legend.position = c(0.2, 0.2),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "lines")))


#### SENSITIVITY ANALYSIS ####

# temperature change data
temp.overtime <- read.csv("data/22warmq_timeseries.csv") %>% 
  rename(YEAR = year)

# keep temp data from 5 years prior to each survey
subsite.info.year <- itex.dec22 %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>%
  dplyr::select(SiteSubsitePlot, SiteSubsite, LAT, LONG, YEAR) %>% 
  mutate(MeasureYear = YEAR) %>%
  mutate(YEAR1 = YEAR - 1, YEAR2 = YEAR - 2, YEAR3 = YEAR - 3, YEAR4 = YEAR - 4) %>%
  pivot_longer(
    cols = starts_with("YEAR"),
    names_to = "YEAR_prior",
    names_prefix = "YEAR",
    values_to = "YEAR",
    values_drop_na = T
  )

# combine with prec and temp data through coordinates
temp.alldata <- left_join(temp.overtime, subsite.info.year, by = c("LAT", "LONG", "YEAR"))

yr5_mean <- temp.alldata %>% group_by(SiteSubsitePlot, SiteSubsite, LAT, LONG, MeasureYear) %>% 
  summarise(value_mean5yr = mean(value)) %>% rename(YEAR = MeasureYear)

write.csv(temp.alldata.mean, file = "data/clim_time_series_5yr.csv")


# merge with shrub data
shrub.temp <- left_join(fg.values.long.nona, yr5_mean, by = c("SiteSubsitePlot", "YEAR")) %>% 
  filter(!is.na(value_mean5yr)) %>%
  mutate(ShrubCoverProp = (CoverFG/100)) %>% 
  group_by(SiteSubsite.y) %>%
  dplyr::mutate(value_mean5yr_cent = scale(value_mean5yr, scale = F)) %>%
  na.omit() %>%
  ungroup()


# Shrub change with warming change and FG interaction
shrub.temp.mod <- brm(bf(ShrubCoverProp ~ value_mean5yr_cent*FuncGroup + 
                        (value_mean5yr_cent|SiteSubsite.y/SiteSubsitePlot), zoi ~ 1, coi ~ 1),
                      family = "zero_one_inflated_beta",
                      data = shrub.temp, init = 0, 
                      prior = prior(gamma(0.01, 0.01), class = phi),
                      iter = 3500, chains = 4, warmup = 500,
                      file = "scripts/mgarciacriado/models/shrub_temp_mod")

summary(shrub.temp.mod) # negative significant
conditional_effects(shrub.temp.mod)


# model predictions
shrub.temp.pred <- ggpredict(shrub.temp.mod, terms = c("value_mean5yr_cent", "FuncGroup"), type = "random",
                             condition = c(SiteSubsite.y = 1266), allow_new_levels = TRUE)
colnames(shrub.temp.pred) = c('MeanFiveYears', 'fit', 'lwr', 'upr', "FuncGroup")


# plot
(shrub.temp.plot <- ggplot() +
  geom_point(data = shrub.temp, aes(x = value_mean5yr_cent, y = ShrubCoverProp*100, colour = FuncGroup), size = 1.5, alpha = 0.1) +
  scale_colour_manual(values = c("#8b4539", "#215028")) +
  geom_line(aes(x = MeanFiveYears, y = fit*100), colour = "#8b4539", linewidth = 1, linetype = "dashed",
            data = shrub.temp.pred %>% filter(FuncGroup == "DwarfShrubCoverNew")) + 
  geom_ribbon(aes(x = MeanFiveYears, ymin = lwr*100, ymax = upr*100), fill = "#8b4539", alpha = 0.2,
              data = shrub.temp.pred %>%
                filter(FuncGroup == "DwarfShrubCoverNew")) +
  geom_line(aes(x = MeanFiveYears, y = fit*100), colour = "#215028", linewidth = 1, linetype = "dashed",
            data = shrub.temp.pred %>% filter(FuncGroup == "NonDwarfShrubCoverNew")) + 
  geom_ribbon(aes(x = MeanFiveYears, ymin = lwr*100, ymax = upr*100), fill = "#215028", alpha = 0.2,
              data = shrub.temp.pred %>%
                filter(FuncGroup == "NonDwarfShrubCoverNew")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
    ylim(0, 100) +
    xlab("\nMTWQ average from previous\n5 years (scaled, °C)") + ylab("Shrub cover (%)\n") + 
    bio.theme + 
    theme(legend.position = "none", 
          axis.title.x = element_text(face = "plain", size = 9),
          axis.text.x  = element_text(vjust = 0.5, size = 9, colour = "black"), 
          axis.title.y = element_text(face = "plain", size = 9),
          axis.text.y  = element_text(vjust = 0.5, size = 9, colour = "black"),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "lines")))



## ED Fig 9
panel.s9 <- plot_grid(shrubchg.plot.col, shrub.temp.plot, dwarf.temp.plot, 
                       labels = c("a", "b", "c"), label_size = 14, nrow = 1, ncol = 3)
ggplot2::ggsave(panel.s9, filename = "figures/figure_s9.jpg", 
                width = 18, height = 6.4, units = "cm")
