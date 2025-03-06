## Plant diversity dynamics over space and time in a warming Arctic
## Mariana Garcia Criado (mariana.garcia.criado@gmail.com)
## Script 5. Evenness


## LIBRARIES ----
library(tidyverse)
library(vegan)
library(viridis)
library(brms)
library(cowplot)
library(broom)
library(ggeffects)


## THEME ----
theme_fancy <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica"),
          axis.text = element_text(size = 16, color="black"), 
          axis.title = element_text(size = 18, color="black"),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0, color="black"),
          legend.text = element_text(size = 16, color="black"),          
          legend.title = element_text(size = 18, color="black"),                              
          legend.position = c(0.9, 0.9), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))}


## LOAD DATA ----

# Composition data
load("data/itex_dec22.RData")
itex.dec22$SurveyedArea <- as.numeric(itex.dec22$SurveyedArea)

# Climate data
clim.only <- read.csv("data/22clim_only.csv")

# Turnover data
beta_jaccard_final <- read.csv("data/22beta_jaccard_final.csv")




## SPATIAL EVENNESS ----

## END YEAR COMPARISON OF EVENNESS ACROSS SITES

# Pielou’s evenness: J = H/ log(S), where H is Shannon's diversity index and S is the total number of species in a sample.
# We can calculate Shannon easily with vegan package: H <- diversity(BCI)
# Data structure for diversity() is a column for sites and then a column per species with the abundance

# We will only keep last year here, but below when we compare start/end, we will filter for those with 2 timepoints
last.year <- itex.dec22  %>% group_by(SiteSubsitePlot) %>% 
  filter(YEAR == max(YEAR)) %>% ungroup() %>%
  dplyr::select(SiteSubsite, SiteSubsitePlot, SiteSubsitePlotYear, YEAR, SPECIES_NAME, RelCover, MOISTURE, LAT, LONG,
                SurveyedArea, Region, PlotShrubMean, PlotGraminoidMean, PlotForbMean, AnnualRichness, MeanRichness)

# Long-format for vegan, replace NAs with 0s
shannon.plots <- last.year %>% 
  dplyr::select(SiteSubsitePlot, SPECIES_NAME, RelCover) %>% 
  pivot_wider(names_from = SPECIES_NAME, values_from = RelCover) %>%
  replace(is.na(.), 0) %>% arrange(SiteSubsitePlot)

# Calculate Shannon index with vegan package
shannon.ind <- as.data.frame(diversity(shannon.plots[-1], index = "shannon"))
colnames(shannon.ind) <- "Shannon"

# Extract plot data - one row per site, A-Z
plot.data <- last.year %>% distinct(SiteSubsitePlot, .keep_all = TRUE) %>% 
  arrange(SiteSubsitePlot)

# Merge the two dataframes
even.df <- cbind(plot.data, shannon.ind)

# Calculate Pilou's evenness (J = H/ log(S))
evenness <- even.df %>% mutate(PilouEvenness = Shannon / log(AnnualRichness)) # 1327 obs
hist(evenness$PilouEvenness) 
# values between 0.04 and 1 - beta distribution (but we bring down the 1s with a constant below)

# Save this for the temporal richness models
write.csv(evenness, "data/22evenness_spatial.csv")

# What's the mean evenness across plots? (Intercept-only)
even.intonly.mod <- brm(PilouEvenness ~ 1, data = evenness, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/22_00_even_int_mod")
summary(even.intonly.mod) # mean evenness is 0.7, plots are pretty even




#### 1) GEOGRAPHICAL MODEL ####
evenness <- read.csv("data/22evenness_spatial.csv")

# Prepare variables for model, bring down Pilou so beta distribution is acceptable
evenness2 <- evenness %>% mutate(PilouEvenCorrected = PilouEvenness - 0.0001) %>%
  mutate(LogPlotSize = log(SurveyedArea)) %>% 
  mutate(ShrubDivided = PlotShrubMean/100, ForbDivided = PlotForbMean/100, GramDivided = PlotGraminoidMean/100)


# Fit geo model with random effect and no plot size
even.geo.noplot.mod <- brm(PilouEvenCorrected ~ LAT + Region + MeanRichness + (1|SiteSubsite), 
                    data = evenness2, iter = 2000, chains = 4, warmup = 500, 
                    family = "beta", file = "models/22_01_even_geo_noplot_mod")

print(summary(even.geo.noplot.mod), digits = 4) # latitude ns, plot size ns, richness positive
conditional_effects(even.geo.noplot.mod) # NAM-wast greater than GreenIceland and Eurasia

summary(even.geo.noplot.mod, prob = 0.975) # latitude ns, plot size ns, richness positive
conditional_effects(even.geo.noplot.mod, prob = 0.975) # NAM-west greater than GreenIceland and Eurasia

posterior_summary(even.geo.noplot.mod)



#### 2) CLIMATE MODEL ####

# Join evenness with climate data
even_clim <- left_join(evenness2, clim.only, by = "SiteSubsitePlot")

# Clean up variables
even_clim2 <- even_clim %>% filter(MOISTURE != "")


# Climate model without plot size
even.clim.nolog.mod <- brm(PilouEvenCorrected ~ MOISTURE + warmq + prec + 
                       MeanRichness + (1|SiteSubsite),
                     data = even_clim2, iter = 2000, chains = 4, warmup = 500, 
                     family = "beta", file = "models/22_02_even_clim_nolog_mod")
print(summary(even.clim.nolog.mod), digits = 5) # temp ns, precip ns, richness positive, moisture ns
conditional_effects(even.clim.nolog.mod)

print(summary(even.clim.nolog.mod, prob = 0.975), digits = 5) # temp ns, precip ns, richness positive, moisture ns
conditional_effects(even.clim.nolog.mod, prob = 0.975)

posterior_summary(even.clim.nolog.mod)



#### 3) FG-DOM MODEL ####


# Model using mean FG proportion over time rather than final cover [no plot size]
even.shb.nolog.mod <- brm(PilouEvenCorrected ~ ShrubDivided + MeanRichness + (1|SiteSubsite),
                    data = even_clim2, family = "beta",  
                    iter = 2000, chains = 4, warmup = 500, 
                    file = "models/22_03_even_shb_nolog_mod")
print(summary(even.shb.nolog.mod), digits = 3) # shrub negative, richness positive  
summary(even.shb.nolog.mod, prob = 0.975) # shrub negative, richness positive

posterior_summary(even.shb.nolog.mod)


# forb
even.forb.nolog.mod <- brm(PilouEvenCorrected ~ ForbDivided + MeanRichness + (1|SiteSubsite),
                     data = even_clim2, family = "beta",  
                     iter = 2000, chains = 4, warmup = 500, 
                     file = "models/22_03_even_forb_nolog_mod")

print(summary(even.forb.nolog.mod), digits = 3) # forb positive, richness positive
summary(even.forb.nolog.mod, prob = 0.975) # forb positive, richness positive

posterior_summary(even.forb.nolog.mod)


# gram
even.gram.nolog.mod <- brm(PilouEvenCorrected ~ GramDivided + MeanRichness + (1|SiteSubsite),
                     data = even_clim2, family = "beta",  
                     iter = 2000, chains = 4, warmup = 500, 
                     file = "models/22_03_even_gram_nolog_mod")

print(summary(even.gram.nolog.mod), digits = 3) # gram ns, richness positive
summary(even.gram.nolog.mod, prob = 0.975) # gram ns, richness positive

posterior_summary(even.gram.nolog.mod)




## EVENNESS OVER TIME ----

## COMPARISON OF EVENNESS WITHIN A COMMUNITY AND OVER TIME

# Keep all years, long-format for vegan, replace NAs with 0s
shannon.plots.time <- itex.dec22 %>% group_by(SiteSubsitePlot) %>%
  filter(length(unique(YEAR)) > 1) %>% ungroup() %>% 
  group_by(SiteSubsitePlot) %>%
  mutate(Duration = max(YEAR) - min(YEAR)) %>% filter(Duration > 4) %>% ungroup() %>% 
  dplyr::select(SiteSubsitePlotYear, SPECIES_NAME, RelCover) %>% 
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = SPECIES_NAME, values_from = RelCover)  %>%
  group_by(SiteSubsitePlotYear) %>% summarise_each(funs(first(.[!is.na(.)]))) %>%
  replace(is.na(.), 0) %>% dplyr::select(-row) %>% arrange(SiteSubsitePlotYear)

# Calculate Shannon index with vegan package
shannon.ind.time <- as.data.frame(diversity(shannon.plots.time[-1], index = "shannon"))
colnames(shannon.ind.time) <- "Shannon"

# Merge back with original data - one line per plotXyear
plot.data.time <- itex.dec22 %>% group_by(SiteSubsitePlot) %>% 
  filter(length(unique(YEAR)) > 1) %>% mutate(Duration = max(YEAR) - min(YEAR)) %>% 
  filter(Duration > 4) %>% ungroup() %>%
  dplyr::select(SiteSubsitePlotYear, SiteSubsitePlot, YEAR, LAT, LONG,
                SiteSubsite, RelCover, MOISTURE, SurveyedArea, Region, PlotShrubMean, 
                PlotGraminoidMean, PlotForbMean, AnnualRichness, MeanRichness) %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% arrange(SiteSubsitePlotYear) %>%
  group_by(SiteSubsitePlot) %>% mutate(Duration = max(YEAR) - min(YEAR)) %>% ungroup() %>% 
  mutate(LogPlotSize = log(SurveyedArea))

# Merge the two dataframes
even.df.time <- cbind(plot.data.time, shannon.ind.time)

# Calculate Pilou's evenness (J = H/ log(S))
evenness.time <- even.df.time %>% 
  mutate(PilouEvenness = Shannon / log(AnnualRichness)) %>% 
  tidyr::drop_na(PilouEvenness) # 4822 obs
# removing NaN values of J because Shannon = 1, this is only where there is one single species that is 100% cover
# so it would be misleading to replace it by a 0, so just removing all together.

hist(evenness.time$PilouEvenness) # various 1, no 0. 

# Save 
write.csv(evenness.time, "data/22even_time.csv")



#### 00) HIGH-LEVEL MODEL ####
evenness.time <- read.csv("data/22even_time.csv")

# substracting a tiny constant so that beta distribution fits
evenness.timeX <- evenness.time %>% 
  mutate(PilouBeta = PilouEvenness - 0.001) %>% 
  mutate(YearCentred = YEAR - mean(YEAR))

# High-level model with nested plot random effect
big.even.time.nest.mod <- brm(PilouBeta ~ YearCentred + (YearCentred|SiteSubsite/SiteSubsitePlot), 
                              family = "beta", data = evenness.timeX, init = 0,
                              iter = 2000, chains = 4, warmup = 400, 
                              file = "models/24_00_even_year_ranslopnest_mod")
print(summary(big.even.time.nest.mod), digits = 4) # ns
plot(big.even.time.nest.mod)



#### 1) EVENNESS CHANGE OVER TIME ####
#evenness.time <- read.csv("data/22even_time.csv")

# Calculate slopes (annual evenness change) for each plot 
even.time.slopes <- evenness.time %>%
  nest_by(SiteSubsitePlot) %>%  # replace group_by() with nest_by() to convert your model data to a vector of lists
  dplyr::mutate(mod = list(lm(PilouEvenness ~ YEAR, data = data))) %>% # change do() to mutate(), then add list() before your model make sure to change data = .  to data = data
  summarise(tidy(mod)) %>%
  filter(term == "YEAR") %>%
  dplyr::select(SiteSubsitePlot, estimate) %>% 
  rename(Slope = estimate) %>%
  tidyr::drop_na(Slope)

write.csv(even.time.slopes, "data/22_even_time_slopes.csv")

even.time.slopes <- read.csv("data/22_even_time_slopes.csv")

# Evenness change over time
hist(even.time.slopes$Slope) # centred on 0

# What's the mean evenness change across plots? (Intercept-only)
eventime.intonly.mod <- brm(Slope ~ 1, data = even.time.slopes, iter = 2000, chains = 4, warmup = 400, 
                            file = "models/24_00_eventime_int_mod")
print(summary(eventime.intonly.mod), digits = 5) # mean slope of evenness change is 0.00018 (overlapping 0)




#### 2) SUBSITE LEVEL MODEL ####

# Now we need one row per plot only
even.time.linexplot0 <- evenness.time %>% distinct(SiteSubsitePlot, .keep_all = TRUE) %>% 
  dplyr::select(SiteSubsite, SiteSubsitePlot, MOISTURE, LAT, Region, LogPlotSize, 
                MeanRichness, Duration, PlotForbMean, PlotShrubMean, PlotGraminoidMean)

# Join with slope data
even.time.linexplot <- left_join(even.time.linexplot0, even.time.slopes, by = "SiteSubsitePlot")


# model without plot size
even.time.subsite.nolog.mod <- brm(Slope ~ LAT + Region + 
                               MeanRichness + Duration + (1|SiteSubsite), 
                             data = even.time.linexplot, family = 'gaussian',
                             iter = 2000, chains = 4, warmup = 400, 
                             file = "scripts/mgarciacriado/models/24_04_even_time_subs_nolog_mod")
print(summary(even.time.subsite.nolog.mod), digits = 5) # lat ns, richness ns, duration ns
conditional_effects(even.time.subsite.nolog.mod) # region ns. 

print(summary(even.time.subsite.nolog.mod, prob = 0.975), digits = 5) # lat ns, richness ns, duration ns
conditional_effects(even.time.subsite.nolog.mod, prob = 0.975) # region ns.




#### 3) PLOT LEVEL MODEL WITH CHANGE (CLIMATE + FUNCTIONAL GROUP) ####
change_clim <- read.csv("data/22clim_slopes.csv")
change_fg <- read.csv("data/22fg_slopes.csv")

even.time.change0 <- left_join(even.time.linexplot, change_clim, by = "SiteSubsitePlot")
even.time.change <- left_join(even.time.change0, change_fg, by = "SiteSubsitePlot")
even.time.change2 <- even.time.change %>% filter(MOISTURE != "")

hist(even.time.change2$Slope)



# PCHG model - one by one FG [no plot size]

# shrubs
even.time.shbch.nolog.plot.mod <- brm(Slope ~ MOISTURE + WarmQSlope + PrecSlope + ShrubSlope + 
                                  MeanRichness + Duration + (1|SiteSubsite), 
                                data = even.time.change2, family = 'gaussian',
                                iter = 2000, chains = 4, warmup = 400, 
                                file = "models/24_05_even_time_plot_shbch_noplot_mod")

print(summary(even.time.shbch.nolog.plot.mod), digits = 5) # temp ns, prec ns, shrub negative, richness ns, duration ns
conditional_effects(even.time.shbch.nolog.plot.mod) # wet more even than dry

print(summary(even.time.shbch.nolog.plot.mod, prob = 0.975), digits = 5) # temp ns, prec ns, shrub negative, richness ns, duration ns
conditional_effects(even.time.shbch.nolog.plot.mod, prob = 0.975) # wet more even than dry



# forbs - one by one FG
even.time.forbch.nolog.plot.mod <- brm(Slope ~ MOISTURE + WarmQSlope + PrecSlope + ForbSlope + 
                                   MeanRichness + Duration + (1|SiteSubsite), 
                                 data = even.time.change2, family = 'gaussian',
                                 iter = 2000, chains = 4, warmup = 400, 
                                 file = "models/24_06_even_time_plot_forbch_noplot_mod")

print(summary(even.time.forbch.nolog.plot.mod), digits = 5) # temp ns, prec ns, forb positive, richness positive, duration ns
conditional_effects(even.time.forbch.nolog.plot.mod) # wet became more even than dry

print(summary(even.time.forbch.nolog.plot.mod, prob = 0.975), digits = 5) # temp ns, prec ns, forb positive, richness ns, duration ns
conditional_effects(even.time.forbch.nolog.plot.mod, prob = 0.975) # moisture ns

posterior_summary(even.time.forbch.nolog.plot.mod)


# grams - one by one FG
even.time.gramch.plot.nolog.mod <- brm(Slope ~ MOISTURE + WarmQSlope + PrecSlope + GraminoidSlope + 
                                   MeanRichness + Duration + (1|SiteSubsite), 
                                 data = even.time.change2, family = 'gaussian',
                                 iter = 2000, chains = 4, warmup = 400, 
                                 file = "models/24_07_even_time_plot_gramch_noplot_mod")


print(summary(even.time.gramch.plot.nolog.mod), digits = 5) # temp ns, prec ns, gram ns, richness ns, duration ns
conditional_effects(even.time.gramch.plot.nolog.mod) # wet became more even than dry

print(summary(even.time.gramch.plot.nolog.mod, prob = 0.975), digits = 5) # temp ns, prec ns, gram ns, richness ns, duration ns
conditional_effects(even.time.gramch.plot.nolog.mod, prob = 0.975) # moisture ns
posterior_summary(even.time.gramch.plot.nolog.mod)





## SUMMARY FIGS ----

## START VS END EVENNESS

# Comparing start and end year
evenness.time.start <- evenness.time %>% group_by(SiteSubsitePlot) %>% filter(YEAR == min(YEAR)) %>% 
  ungroup() %>% filter(!is.na(PilouEvenness))
evenness.time.end <- evenness.time %>% group_by(SiteSubsitePlot) %>% filter(YEAR == max(YEAR)) %>% 
  ungroup() %>% filter(!is.na(PilouEvenness))

# What's the mean start evenness across plots? (Intercept-only)
even.start.intonly.mod <- brm(PilouEvenness ~ 1, data = evenness.time.start, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/24even_start_intonly_mod")
summary(even.start.intonly.mod) # 0.72 (0.71-0.73)

# What's the mean end evenness across plots? (Intercept-only), should be same as spatial value (it is!)
even.end.intonly.mod <- brm(PilouEvenness ~ 1, data = evenness.time.end, iter = 2000, chains = 4, warmup = 400, 
                              file = "models/24even_end_intonly_mod")
summary(even.end.intonly.mod) # 0.72 (0.71-0.73)

# Calculate mean evenness at start and end time (to cross-check with intercept-only mod results)
start.even <- evenness.time.start %>% dplyr::pull(PilouEvenness) %>% mean() # 0.72
end.even <- evenness.time.end %>% dplyr::pull(PilouEvenness) %>% mean() #0.72






## EVENNESS VS TRAJECTORIES ----

## Trajectory and turnover data
all.trends.wide <- read.csv("data/22all_trends_wide.csv")
beta_turnover <- read.csv("data/22beta_turnover.csv")

# Calculate mean evenness per plot
evenness.mean <- evenness.time %>% group_by(SiteSubsitePlot) %>% 
  mutate(EvennessMean = mean(PilouEvenness)) %>% ungroup %>%
  distinct(SiteSubsitePlot, .keep_all = TRUE)


# Merge with trajectories
all.tomerge <- all.trends.wide %>% dplyr::select(SiteSubsitePlot, Persisting, Extinct, Coloniser)
even.mean.trends <- left_join(evenness.mean, all.tomerge, by = "SiteSubsitePlot")

brayonly <- beta_turnover %>% dplyr::select(SiteSubsitePlot, BrayCurtis)
even.mean.trends2 <- left_join(even.mean.trends, brayonly, by = "SiteSubsitePlot")



#### 01) EXTINCT #### 
even.ext.mod <- brm(bf(Extinct ~ EvennessMean + (1|SiteSubsite), zi ~ 1), 
                     data = even.mean.trends2, 
                    iter = 2000, chains = 4, warmup = 400, 
                    family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi),
                    file = "models/24even_ext_mod")
summary(even.ext.mod) # negative significant
print(summary(even.ext.mod), prob = 0.975) # negative significant


# predictions
even.ext.predXX <- ggpredict(even.ext.mod, term = "EvennessMean",  
                             type = "random", condition = c(SiteSubsite = 1266), allow_new_levels = TRUE)
colnames(even.ext.predXX) = c('EvennessMean', 'fit', 'lwr', 'upr', 'group')


# plot
(even.ext.plotXX <- ggplot() + 
    geom_point(data = even.mean.trends2, aes(x = EvennessMean, y = Extinct*100, colour = BrayCurtis), size = 2) +
    geom_line(data = even.ext.predXX, aes(x = EvennessMean, y = fit*100), colour = "black", size = 1) + 
    geom_ribbon(data = even.ext.predXX, aes(x = EvennessMean, ymin = lwr*100, ymax = upr*100), fill = "black", alpha = 0.2) +
    xlab("\nEvenness (plot mean)") +
    ylab("Species losses (%)\n") + ylim(0, 100) +
    scale_colour_viridis("Turnover (B-C)") + 
    theme_fancy() + 
    theme(axis.title.x = element_text(face="plain", size=9), 
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y  = element_text(vjust=0.5, size=9, colour = "black"), 
          legend.key.size = unit(0.3, 'cm'),
          legend.position = c(0.3, 0.8),
          legend.text = element_text(size = 9, color="black"),          
          legend.title = element_text(size = 9, color="black"),           
          plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")))





#### 02) COLONISERS ####
even.col.mod <- brm(bf(Coloniser ~ EvennessMean + (1|SiteSubsite), zi ~ 1), 
                     data = even.mean.trends2, iter = 2000, chains = 4, warmup = 400, 
                    family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi),
                    file = "models/24even_col_mod")
summary(even.col.mod) # negative significant
print(summary(even.col.mod), prob = 0.975) # negative significant

# predictions
even.col.predXX <- ggpredict(even.col.mod, term = "EvennessMean", type = "random", 
                             condition = c(SiteSubsite = 1266), allow_new_levels = TRUE)
colnames(even.col.predXX) = c('EvennessMean', 'fit', 'lwr', 'upr', 'group')


# plot
(even.col.plotXX <- ggplot() + 
    geom_point(data = even.mean.trends2, aes(x = EvennessMean, y = Coloniser*100, colour = BrayCurtis), size = 2) +
    geom_line(data = even.col.predXX, aes(x = EvennessMean, y = fit*100), colour = "black",  size = 1) + 
    geom_ribbon(data = even.col.predXX, aes(x = EvennessMean, ymin = lwr*100, ymax = upr*100), fill = "black", alpha = 0.2) +
    xlab("\nEvenness (plot mean)") +
    ylab("Species gains (%)\n") +
    scale_colour_viridis("Turnover (B-C)") + ylim(0, 100) +
    theme_fancy() +
    theme(axis.title.x = element_text(face="plain", size=9), 
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y  = element_text(vjust=0.5, size=9, colour = "black"), 
          legend.position='none',
          plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")))






#### 03) PERSISTERS ####
even.pers.mod <- brm(bf(Persisting ~ EvennessMean + (1|SiteSubsite), zoi ~ 1, coi ~ 1), 
                      data = even.mean.trends2, iter = 2000, chains = 4, warmup = 400, 
                     family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi),
                     file = "models/24even_pers_mod")
summary(even.pers.mod) # positive significant
print(summary(even.pers.mod), prob = 0.975) # positive significant

# predictions
even.pers.predXX <- ggpredict(even.pers.mod, term = "EvennessMean", type = "random", 
                              condition = c(SiteSubsite = 1266), allow_new_levels = TRUE)
colnames(even.pers.predXX) = c('EvennessMean', 'fit', 'lwr', 'upr', 'group')

# plot
(even.pers.plotXX <- ggplot() + 
    geom_point(data = even.mean.trends2, aes(x = EvennessMean, y = Persisting*100, colour = BrayCurtis), size = 2) +
    geom_line(data = even.pers.predXX, aes(x = EvennessMean, y = fit*100), colour = "black",  size = 1) + 
    geom_ribbon(data = even.pers.predXX, aes(x = EvennessMean, ymin = lwr*100, ymax = upr*100), fill = "black", alpha = 0.2) +
    xlab("\nEvenness (plot mean)") +
    ylab("Persisting species (%)\n") +
    scale_colour_viridis("Turnover (B-C)") + 
    theme_fancy() +
    theme(axis.title.x = element_text(face="plain", size=9), 
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y  = element_text(vjust=0.5, size=9, colour = "black"), 
          legend.position='none',
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "lines")))


# Panel
(even.extcol.panel <- plot_grid(even.ext.plotXX, even.pers.plotXX, even.col.plotXX, 
                                ncol = 3, nrow = 1, align = "hv", 
                                font.label = list(size = 14, face = "bold"), 
                                labels = c("d", "e", "f")))


## ED Fig 8
fig.s8 <- plot_grid(prop.pool.panel, even.extcol.panel, duration.panel,
                    ncol = 1, nrow = 3, align = "h")
ggplot2::ggsave(fig.s8, filename = "figures/figure_s8.jpg", 
                width = 18, height = 24, units = "cm")
