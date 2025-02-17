## Plant diversity dynamics over space and time in a warming Arctic
## Mariana Garcia Criado (mariana.garcia.criado@gmail.com)
## Script 4. Extinct/Coloniser dynamics


## LIBRARIES ----
library(tidyverse)
library(viridis)
library(brms)
library(cowplot)
library(ggthemes)
library(ggeffects)
library(modelr)
library(corrplot)



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
                                           size = 2, linetype = "blank"))
}


## LOAD DATA ----
load("data/itex_dec22.RData")
rich.time <- read.csv("data/22rich_time.csv")
clim.only <- read.csv("data/22clim_only.csv")
beta_jaccard_final <- read.csv("data/22beta_jaccard_final.csv")

turnover <- beta_jaccard_final %>% 
  dplyr::select(SiteSubsitePlot, richness_change, richness_change_abs, Jaccard, BrayCurtis)


## DATA PREP ----

# Remove plots with one recorded year only, keep baseline and resurvey only, add species trajectory
# Reminder: abundance = 0 rows should not be here but let's make sure to remove them just in case
col.ext <- itex.dec22 %>% group_by(SiteSubsitePlot) %>% 
  mutate(MinYear = min(YEAR)) %>% mutate(MaxYear = max(YEAR)) %>%
  filter(MinYear != MaxYear) %>% filter(YEAR == min(YEAR) | YEAR == max(YEAR)) %>% ungroup() %>% 
  mutate(Duration = MaxYear - MinYear) %>% filter(Duration > 4) %>% filter(RelCover > 0) %>%
  group_by(SiteSubsitePlot, SPECIES_NAME) %>% mutate(SpeciesTimePoint = n()) %>% 
  mutate(trend = case_when(SpeciesTimePoint > 1 ~ "Persisting",
                           SpeciesTimePoint = 1 & YEAR == MinYear ~ "Extinct",
                           SpeciesTimePoint = 1 & YEAR == MaxYear ~ "Coloniser",
                           TRUE ~ "none")) %>% ungroup()

# Filter per trend, keep only one row per species and plot
persisting <- col.ext %>% filter(trend == "Persisting") %>% distinct(SiteSubsitePlot, SPECIES_NAME, .keep_all = TRUE) %>%
  group_by(SiteSubsitePlot) %>% mutate(NbSpeciesPerTrendPlot = length(unique(SPECIES_NAME))) %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>% dplyr::select(SiteSubsitePlot, NbSpeciesPerTrendPlot, trend, LAT)

extinct <- col.ext %>% filter(trend == "Extinct") %>% distinct(SiteSubsitePlot, SPECIES_NAME, .keep_all = TRUE) %>%
  group_by(SiteSubsitePlot) %>% mutate(NbSpeciesPerTrendPlot = length(unique(SPECIES_NAME))) %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>% dplyr::select(SiteSubsitePlot, NbSpeciesPerTrendPlot, trend, LAT)

colonisations <- col.ext %>% filter(trend == "Coloniser") %>% distinct(SiteSubsitePlot, SPECIES_NAME, .keep_all = TRUE) %>%
  group_by(SiteSubsitePlot) %>% mutate(NbSpeciesPerTrendPlot = length(unique(SPECIES_NAME))) %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>% dplyr::select(SiteSubsitePlot, NbSpeciesPerTrendPlot, trend, LAT)

# One row per plot/trend
all.trends <- rbind(persisting, extinct, colonisations)


# Note that at this point every plot doesn't have 3 rows as there are particular trends that won't be represented in all plots.
# They will be introduced below with 0s so the models are representative.



## SPECIES & FG ----

# Which are the most common species per trend?

# Extract key of species/FG
name.gf <- col.ext %>% dplyr::select(SPECIES_NAME, FuncGroup) %>% distinct(SPECIES_NAME, .keep_all = TRUE)

# Proportion of a given trend for each species - TOP 10 LIST
xxx <- col.ext %>% distinct(SiteSubsitePlot, SPECIES_NAME, .keep_all = TRUE) %>% 
  group_by(SPECIES_NAME) %>% mutate(SpeciesTotalOccurrence = n())

# Note that the species are still repeated in 'xxx', it's not a single line per species.
# We need this to calculate the proportions below.

pers.prop <- xxx %>% filter(trend == "Persisting") %>% group_by(SPECIES_NAME) %>% 
  mutate(TimesInTrend = n()) %>% mutate(PropInTrend = (TimesInTrend/SpeciesTotalOccurrence)*100) %>% 
  distinct(SPECIES_NAME, .keep_all = TRUE) %>% dplyr::select(SPECIES_NAME, TimesInTrend, PropInTrend) %>% 
  arrange(-TimesInTrend)

pers.prop.gf <- left_join(pers.prop, name.gf, by = "SPECIES_NAME")
write.csv(pers.prop.gf, "data/22top_persisters.csv")

ext.prop <- xxx %>% filter(trend == "Extinct") %>% group_by(SPECIES_NAME) %>% 
  mutate(TimesInTrend = n()) %>% mutate(PropInTrend = (TimesInTrend/SpeciesTotalOccurrence)*100) %>% 
  distinct(SPECIES_NAME, .keep_all = TRUE) %>% dplyr::select(SPECIES_NAME, TimesInTrend, PropInTrend) %>% 
  arrange(-TimesInTrend)

ext.prop.gf <- left_join(ext.prop, name.gf, by = "SPECIES_NAME")
write.csv(ext.prop.gf, "data/22top_extinct.csv")

col.prop <- xxx %>% filter(trend == "Coloniser") %>% group_by(SPECIES_NAME) %>% 
  mutate(TimesInTrend = n()) %>% mutate(PropInTrend = (TimesInTrend/SpeciesTotalOccurrence)*100) %>% 
  distinct(SPECIES_NAME, .keep_all = TRUE) %>% dplyr::select(SPECIES_NAME, TimesInTrend, PropInTrend) %>%
  arrange(-TimesInTrend)

col.prop.gf <- left_join(col.prop, name.gf, by = "SPECIES_NAME")
write.csv(col.prop.gf, "data/22top_coloniser.csv")





## DONUGHT CHARTS

# Calculate proportion of persisters per functional group
persisting.pc <- pers.prop.gf %>% group_by(FuncGroup) %>% summarise(n = n()) %>% 
  ungroup() %>% mutate(freq = n / sum(n))

# Calculate proportion of extincts per functional group
extinct.pc <- ext.prop.gf %>% group_by(FuncGroup) %>% summarise(n = n()) %>% 
  ungroup() %>% mutate(freq = n / sum(n))

# Calculate proportion of extincts per functional group
coloniser.pc <- col.prop.gf %>% group_by(FuncGroup) %>% 
  summarise(n = n()) %>% ungroup() %>% mutate(freq = n / sum(n))


## Donut graph - persisting
persisting.pc$percent <- round((persisting.pc$freq*100), 1)

# Compute the cumulative percentages (top of each rectangle)
persisting.pc$ymax = cumsum(persisting.pc$freq)

# Compute the bottom of each rectangle
persisting.pc$ymin = c(0, head(persisting.pc$ymax, n=-1))

# Compute label position
persisting.pc$labelPosition <- (persisting.pc$ymax + persisting.pc$ymin) / 2

# Compute a good label
persisting.pc$label <- paste0(persisting.pc$FuncGroup, " ", persisting.pc$percent, "%")

# Total donught
(persis.don <- ggplot(persisting.pc, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=FuncGroup)) +
    geom_rect() +
    geom_text(x=2, aes(y=labelPosition, label=label), color="black", size=3) + # x here controls label position (inner / outer)
    scale_fill_manual(values = c("#95558b", "#dc9537", "#5fbf33")) +
    coord_polar(theta="y") +
    xlim(c(-1, 4)) + ggtitle("Persisters") +
    theme_void() +
    theme(legend.position = "none", 
          plot.title = element_text(size=10, face="plain", hjust = 0.5, colour = "black"),
          plot.margin = unit(c(0, 0, 0, 0), units = , "cm")))




## Donut graph - extinct
extinct.pc$percent <- round((extinct.pc$freq*100), 1)

# Compute the cumulative percentages (top of each rectangle)
extinct.pc$ymax = cumsum(extinct.pc$freq)

# Compute the bottom of each rectangle
extinct.pc$ymin = c(0, head(extinct.pc$ymax, n=-1))

# Compute label position
extinct.pc$labelPosition <- (extinct.pc$ymax + extinct.pc$ymin) / 2

# Compute a good label
extinct.pc$label <- paste0(extinct.pc$FuncGroup, " ", extinct.pc$percent, "%")

# Total donught
(extinct.don <- ggplot(extinct.pc, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=FuncGroup)) +
    geom_rect() +
    geom_text(x=2, aes(y=labelPosition, label=label), color="black", size=3) + 
    scale_fill_manual(values = c("#95558b", "#dc9537", "#5fbf33")) +
    coord_polar(theta="y") +
    xlim(c(-1, 4)) + ggtitle("Losses") +
    theme_void() +
    theme(legend.position = "none", 
          plot.title = element_text(size=10, face="plain", hjust = 0.5, colour = "black"),
          plot.margin = unit(c(0, 0, 0, 0), units = , "cm")))



## Donut graph - coloniser
coloniser.pc$percent <- round((coloniser.pc$freq*100), 1)

# Compute the cumulative percentages (top of each rectangle)
coloniser.pc$ymax = cumsum(coloniser.pc$freq)

# Compute the bottom of each rectangle
coloniser.pc$ymin = c(0, head(coloniser.pc$ymax, n=-1))

# Compute label position
coloniser.pc$labelPosition <- (coloniser.pc$ymax + coloniser.pc$ymin) / 2

# Compute a good label
coloniser.pc$label <- paste0(coloniser.pc$FuncGroup, " ", coloniser.pc$percent, "%")

# Total donught
(coloniser.don <- ggplot(coloniser.pc, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=FuncGroup)) +
    geom_rect() +
    geom_text(x=2, aes(y=labelPosition, label=label), color="black", size=3) + 
    scale_fill_manual(values = c("#95558b", "#dc9537", "#5fbf33")) +
    coord_polar(theta="y") +
    xlim(c(-1, 4)) + ggtitle("Gains") + 
    theme_void() +
    theme(legend.position = "none", 
          plot.title = element_text(size=10, face="plain", hjust = 0.5, colour = "black"),
          plot.margin = unit(c(0, 0, 0, 0), units = , "cm")))



## Overall proportions of FGs in the database
fg.prop <- col.ext %>% distinct(SiteSubsitePlot, SPECIES_NAME, .keep_all = TRUE) %>% 
  group_by(SPECIES_NAME) %>% summarise(n()) %>% ungroup()

# Add in functional group info
fg.gf <- left_join(fg.prop, name.gf, by = "SPECIES_NAME")

# Calculate overall proportion of functional group
fg.pc <- fg.gf %>% drop_na(FuncGroup) %>% group_by(FuncGroup) %>% summarise(n = n()) %>% 
  ungroup() %>% mutate(freq = n / sum(n))


## Donut graph - total in database
fg.pc$percent <- round((fg.pc$freq*100), 1)

# Compute the cumulative percentages (top of each rectangle)
fg.pc$ymax = cumsum(fg.pc$freq)

# Compute the bottom of each rectangle
fg.pc$ymin = c(0, head(fg.pc$ymax, n=-1))

# Compute label position
fg.pc$labelPosition <- (fg.pc$ymax + fg.pc$ymin) / 2

# Compute a good label
fg.pc$label <- paste0(fg.pc$FuncGroup, " ", fg.pc$percent, "%")

# Total donught
(fg.don <- ggplot(fg.pc, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=FuncGroup)) +
    geom_rect() +
    geom_text(x=2, aes(y=labelPosition, label=label), color="black", size=3) + 
    scale_fill_manual(values = c("#95558b", "#dc9537", "#5fbf33")) +
    coord_polar(theta="y") +
    xlim(c(-1, 4)) + ggtitle("Total") +
    theme_void() +
    theme(legend.position = "none", 
          plot.title = element_text(size=10, face="plain", hjust = 0.5, colour = "black"),
          plot.margin = unit(c(0, 0, 0, 0), units = , "cm")))


# Panel
(S5eh <- plot_grid(fg.don, extinct.don, persis.don, coloniser.don, 
                           ncol=4, nrow = 1, align="hv", 
                           labels = c("e", "f", "g", "h")))


## Chi-square test of independence
## To test if the proportions per FG and trend reflect that of the main database

# First retain count columns only
ext.chi <- extinct.pc %>% select(FuncGroup, n) %>% rename(ExtinctCount = n)
col.chi <- coloniser.pc %>% select(FuncGroup, n) %>% rename(ColoniserCount = n)
pers.chi <- persisting.pc %>% select(FuncGroup, n) %>% rename(PersistingCount = n)
total.chi <- fg.pc %>% select(FuncGroup, n) %>% rename(OverallCount = n)

# Merge into one dataframe
trends.chi <- left_join(ext.chi, col.chi, by = "FuncGroup") %>% 
  left_join(., pers.chi, by = "FuncGroup") %>% 
  left_join(., total.chi, by = "FuncGroup") %>% 
  column_to_rownames(var = "FuncGroup")

# Convert to contingency table
trends.matrix <- as.table(as.matrix(trends.chi))

chisq <- chisq.test(trends.matrix, correct=FALSE)
chisq # p = 0.4337, the rows and column are not statistically associated

chisq$observed
chisq$expected
round(chisq$residuals, 3)

corrplot(chisq$residuals, is.cor = FALSE)




# Test this with only two trends
two.chi <- left_join(ext.chi, total.chi, by = "FuncGroup") %>% column_to_rownames(var = "FuncGroup")
two.matrix <- as.table(as.matrix(two.chi))

two.chisq <- chisq.test(two.matrix)
two.chisq # p = 0.82


# Test this with only two trends
two.chi2 <- left_join(col.chi, total.chi, by = "FuncGroup") %>% column_to_rownames(var = "FuncGroup")
two.matrix2 <- as.table(as.matrix(two.chi2))

two.chisq2 <- chisq.test(two.matrix2)
two.chisq2 # p = 0.75


# Test this with only two trends
two.chi3 <- left_join(pers.chi, total.chi, by = "FuncGroup") %>% column_to_rownames(var = "FuncGroup")
two.matrix3 <- as.table(as.matrix(two.chi3))

two.chisq3 <- chisq.test(two.matrix3)
two.chisq3 # p = 0.08


# z-test to compare proportions when one is a subset of the other
prop.test(two.matrix3) # p = 0.08
prop.test(two.matrix2) # p = 0.75
prop.test(two.matrix) # p = 0.82

# p > 0.05 so the proportions of the two groups are similar between the global database and the groups



## PLOTS ----

# Calculate means
mu2 <- plyr::ddply(all.trends3, "trend", summarise, grp.mean = mean(ProportionSpecies))

(trend.dens2 <- ggplot(all.trends3, aes(x=ProportionSpecies, fill=trend, color=trend))+
    geom_density(alpha = 0.2, size = 1.2)+
    geom_vline(data = mu2, aes(xintercept = grp.mean, color = trend), linetype="dashed", size=1, show.legend = FALSE) +
    theme_bw()+ bio.theme + 
    annotate("text", x = 85, y = 0.025, label = "105 (8.3%) plots had\n only persisters", size = 4) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_blank() ,
          legend.title = element_blank(),
          legend.position = c(0.88, 0.87),
          legend.text = element_text(size = 9),
          axis.title.x = element_text(face="plain", size=9),
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y = element_text(vjust=0.5, face="plain", size=9),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), units = , "cm")) + 
    scale_fill_manual(values= c("#8738C7", "#d68317", "#38C787")) +
    scale_color_manual(values= c("#8738C7", "#d68317", "#38C787")) +
    xlab("\nSpecies per trajectory (%)")+
    ylab(""))

only.pers <- all.trends2 %>% filter(ProportionSpecies == 100) # 105 plots have only persisters.
length(unique(all.trends2$SiteSubsitePlot)) #1266 -- 8.3% plots have only persisters.

## ED Fig 5
S5d <- plot_grid(trend.dens2, labels = c("d"), label_size = 14)

# panel - top row
(S5.panel <- plot_grid(S5A, S5d, S5eh, ncol = 1, nrow = 3,  
                       rel_heights = c(1, 2, 1),
                  align = "hv"))

# full panel
ggsave(S5.panel, filename = "figures/figure_s5.jpg", 
       width = 18, height = 24, units = "cm")


## Slice into different datasets

# First we need to add 0s to those plots where there are no species of that trend
plot.names <- col.ext %>% distinct(SiteSubsitePlot)

# Join with trend data
extinct.full0 <- left_join(plot.names, extinct, by = "SiteSubsitePlot") 
coloniser.full0 <- left_join(plot.names, colonisations, by = "SiteSubsitePlot") 
persisters.full0 <- left_join(plot.names, persisting, by = "SiteSubsitePlot")

# Replace NA with 0s
extinct.full <- extinct.full0 %>% tidyr::replace_na(list(NbSpeciesPerTrendPlot = 0)) %>% 
  tidyr::replace_na(list(trend = "Extinct"))

coloniser.full <- coloniser.full0 %>% tidyr::replace_na(list(NbSpeciesPerTrendPlot = 0)) %>% 
  tidyr::replace_na(list(trend = "Coloniser"))

persisters.full <- persisters.full0 %>% tidyr::replace_na(list(NbSpeciesPerTrendPlot = 0)) %>% 
  tidyr::replace_na(list(trend = "Persisting"))

# Join all trends in one dataframe
alltrends.full <- rbind(extinct.full, coloniser.full, persisters.full)



# Calculate proportions
alltrends.full.prop <- alltrends.full %>% group_by(SiteSubsitePlot) %>% 
  mutate(TotalSpecies = sum(NbSpeciesPerTrendPlot)) %>% ungroup() %>% 
  mutate(ProportionSpecies = NbSpeciesPerTrendPlot/TotalSpecies) %>%
  dplyr::select(., -LAT)

# Get additional covariates
covs <- rich.time %>% distinct(SiteSubsitePlot, .keep_all = TRUE) %>%
  dplyr::select(SiteSubsite, SiteSubsitePlot, LAT, SurveyedArea, Region, MOISTURE, Duration,
                PlotGraminoidMean, PlotShrubMean, PlotForbMean, MeanRichness) 

# Add in extra columns
alltrends.full.prop2 <- left_join(alltrends.full.prop, covs, by = "SiteSubsitePlot")

# Convert to wide form and replace NA by 0
all.trends.wide <- alltrends.full.prop2 %>% 
  pivot_wider(names_from = trend, values_from = ProportionSpecies) %>% 
  group_by(SiteSubsitePlot) %>% summarise_each(funs(first(.[!is.na(.)]))) %>% 
  replace(is.na(.), 0) %>% ungroup()

write.csv(all.trends.wide, "data/22all_trends_wide.csv")
#all.trends.wide <- read.csv("data/22all_trends_wide.csv")

# Join dataframes
all.metrics <- left_join(all.trends.wide, turnover, by = "SiteSubsitePlot")





## MODELS ----

## What's the mean number of species per trend?

# Extinct: intercept-only
ext.intonly.mod <- brm(NbSpeciesPerTrendPlot ~ 1, data = extinct.full, iter = 2000, chains = 4, warmup = 400, 
                       file = "models/22ext_intonly_mod")
summary(ext.intonly.mod) # mean is 1.67 extinct species across sites

# Colonisers: intercept-only
col.intonly.mod <- brm(NbSpeciesPerTrendPlot ~ 1, data = coloniser.full, iter = 2000, chains = 4, warmup = 400, 
                       file = "models/22col_intonly_mod")
summary(col.intonly.mod) # mean is 1.84 coloniser species across sites

# Persisters: intercept-only
pers.intonly.mod <- brm(NbSpeciesPerTrendPlot ~ 1, data = persisters.full, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/22pers_intonly_mod")
summary(pers.intonly.mod) # mean is 5.49 persister species across sites




## Do the numbers per trend differ from each other?

# Model
trends.mod <- brm(NbSpeciesPerTrendPlot ~ trend, data = alltrends.full, iter = 2000, chains = 4, warmup = 400,
                  file = "models/22nbs_trend_mod")
summary(trends.mod) # way more persisters than both
conditional_effects(trends.mod)


# Model with proportions
prop.trends.mod <- brm(ProportionSpecies ~ trend, data = alltrends.full.prop2, iter = 2000, chains = 4, warmup = 400,
                  file = "models/22prop_trend_mod")
summary(prop.trends.mod) # way more persisters than both, more colonisers than extinct
conditional_effects(prop.trends.mod)


# Check Y distributions for the different trends
ext <- alltrends.full.prop2 %>% filter(trend == "Extinct") %>% mutate(LogPlotSize = log(SurveyedArea)) %>%
  mutate(GramDivided = PlotGraminoidMean/100, ShrubDivided = PlotShrubMean/100, ForbDivided = PlotForbMean/100)
pers <- alltrends.full.prop2 %>% filter(trend == "Persisting") %>% mutate(LogPlotSize = log(SurveyedArea)) %>%
  mutate(GramDivided = PlotGraminoidMean/100, ShrubDivided = PlotShrubMean/100, ForbDivided = PlotForbMean/100)
col <- alltrends.full.prop2 %>% filter(trend == "Coloniser") %>% mutate(LogPlotSize = log(SurveyedArea)) %>%
  mutate(GramDivided = PlotGraminoidMean/100, ShrubDivided = PlotShrubMean/100, ForbDivided = PlotForbMean/100)

hist(alltrends.full.prop2$ProportionSpecies)
hist(ext$ProportionSpecies) # (0 - 0.75) [zero-inflated beta]
hist(col$ProportionSpecies) # (0 - 0.75) [zero-inflated beta]
hist(pers$ProportionSpecies) # (0 - 1) [zero-one inflated beta]

# What is the proportion of 1s in the persister data?
prop.one <- pers %>% filter(ProportionSpecies == 1 | ProportionSpecies == 0) # 106/1266 = 0.08
# zoi = 0.08 (8% observations are 0s or 1s)
# coi = 0.98 (of those 8%, most are 1s[106])


# Extinct proportions: intercept-only
extprop.intonly.mod <- brm(ProportionSpecies ~ 1, data = ext, iter = 2000, chains = 4, warmup = 400, 
                           file = "models/22extprop_intonly_mod")
summary(extprop.intonly.mod) # mean proportion is 17% extinct species 

# Colonisers: intercept-only
colprop.intonly.mod <- brm(ProportionSpecies ~ 1, data = col, iter = 2000, chains = 4, warmup = 400, 
                           file = "models/22colprop_intonly_mod")
summary(colprop.intonly.mod) # mean proportion is 19% coloniser species

# Persisters: intercept-only
persprop.intonly.mod <- brm(ProportionSpecies ~ 1, data = pers, iter = 2000, chains = 4, warmup = 400, 
                            file = "models/22persprop_intonly_mod")
summary(persprop.intonly.mod) # mean proportion is 64% persisting species





## PERSISTERS ----


#### 01) GEOGRAPHICAL MODEL ####
pers.geo.mod <- brm(bf(ProportionSpecies ~ LAT + Region + LogPlotSize + 
                          MeanRichness + Duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1), 
                    data = pers, iter = 2000, chains = 4, warmup = 500, 
                    family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi),
                    file = "models/22_01_pers_geo_mod")

summary(pers.geo.mod) # latitude ns, plot size ns, richness positive, duration ns
conditional_effects(pers.geo.mod) # less persisters in GreenlandIceland than NAM-west

print(summary(pers.geo.mod, prob = 0.975)) # latitude ns, plot size ns, richness positive, duration ns
conditional_effects(pers.geo.mod, prob = 0.975) # less persisters in GreenlandIceland than NAM-west



#### 02) CLIMATIC MODEL ####

# Join turnover with climate data
pers_clim <- left_join(pers, clim.only, by = "SiteSubsitePlot")

# Clean up variables
pers_clim2 <- pers_clim %>% filter(MOISTURE != "")

# Persisters model 
pers.clim.mod <- brm(bf(ProportionSpecies ~ MOISTURE + warmq + prec + Duration + 
                        LogPlotSize + MeanRichness + (1|SiteSubsite), zoi ~ 1, coi ~ 1), init = 0,
                      family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                      data = pers_clim2, iter = 2000, chains = 4, warmup = 500, 
                      file = "models/22_02_pers_clim_mod")

print(summary(pers.clim.mod), digits = 5) # temp positive, prec negative, duration ns, plot size ns, richness positive
conditional_effects(pers.clim.mod) # moisture ns

print(summary(pers.clim.mod, prob = 0.975)) # temp positive, prec ns, duration ns, plot size ns, richness positive
conditional_effects(pers.clim.mod, prob = 0.975) # moisture ns

# predictions for climate figure
pers.clim.predXX <- ggpredict(pers.clim.mod, terms = "warmq")
colnames(pers.clim.predXX) = c('warmq', 'fit', 'lwr', 'upr')

pers.clim.predXXX <- ggpredict(pers.clim.mod, terms = "prec")
colnames(pers.clim.predXXX) = c('prec', 'fit', 'lwr', 'upr')




#### 03) FG COMPOSITION MODEL ####

# Model using mean FG proportion over time - one per FG

# shrub
pers.fg.shb.mod <- brm(bf(ProportionSpecies ~ ShrubDivided + LogPlotSize + MeanRichness + 
                        Duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                   family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                   data = pers, iter = 2000, chains = 4, warmup = 500, init = 0,
                   file = "models/22_03_pers_fg_shb_mod")

summary(pers.fg.shb.mod) # shrub ns, plot size ns, richness positive, duration ns
print(summary(pers.fg.shb.mod, prob = 0.975)) # shrub ns, plot size ns, richness positive, duration ns 


# forb
pers.fg.forb.mod <- brm(bf(ProportionSpecies ~ ForbDivided + LogPlotSize + MeanRichness + 
                             Duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                   family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                   data = pers, iter = 3000, chains = 4, warmup = 500, init = 0,
                   file = "models/22_03_pers_fg_forb_mod")

summary(pers.fg.forb.mod) # forb ns, plot size ns, richness positive, duration ns
print(summary(pers.fg.forb.mod, prob = 0.975)) # forb ns, plot size ns, richness positive, duration ns
 


# gram
pers.fg.gram.mod <- brm(bf(ProportionSpecies ~ GramDivided + 
                        LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                   family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                   data = pers, iter = 3000, chains = 4, warmup = 500, init = 0,
                   file = "models/22_03_pers_fg_gram_mod")

summary(pers.fg.gram.mod) # gram ns, plot size ns, richness positive, duration ns
print(summary(pers.fg.gram.mod, prob = 0.975)) # gram ns, plot size ns, richness positive, duration ns




#### 04) CHANGE OCCURRED MODEL ####


# Load up change in climate and FG
change_clim <- read.csv("data/22clim_slopes.csv")
change_fg <- read.csv("data/22fg_slopes.csv")

pers.tur.change0 <- left_join(pers, change_clim, by = "SiteSubsitePlot")
pers.tur.change <- left_join(pers.tur.change0, change_fg, by = "SiteSubsitePlot")



# One by one - shrub
pers.shrubchange.mod <- brm(bf(ProportionSpecies ~ WarmQSlope + PrecSlope + 
                              ShrubSlope + LogPlotSize + MeanRichness + Duration + (1|SiteSubsite),
                              zoi ~ 1, coi ~ 1),
                            family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                            data = pers.tur.change, iter = 2000, chains = 4, warmup = 500, 
                            file = "models/22_04_pers_shrubchange_mod")
summary(pers.shrubchange.mod) # temp negative, prec ns, shrub NS, plot size ns, richness positive, duration ns
print(summary(pers.shrubchange.mod, prob = 0.975)) # temp negative, prec ns, shrub ns, plot size ns, richness positive, duration ns

# One by one - forb
pers.forbchange.mod <- brm(bf(ProportionSpecies ~ WarmQSlope + PrecSlope + 
                             ForbSlope + LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                           family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                           data = pers.tur.change, iter = 2000, chains = 4, warmup = 500, 
                           file = "models/22_04_pers_forbchange_mod")
summary(pers.forbchange.mod) # temp negative, prec negative, forb negative, plot size ns, richness positive, duration ns
print(summary(pers.forbchange.mod, prob = 0.975)) # temp negative, prec ns, forb negative, plot size ns, richness positive, duration ns



# One by one - gram
pers.gramchange.mod <- brm(bf(ProportionSpecies ~ WarmQSlope + PrecSlope + 
                             GraminoidSlope + LogPlotSize + MeanRichness + 
                               Duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                           family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                           data = pers.tur.change, iter = 2000, chains = 4, warmup = 500, 
                           file = "models/22_04_pers_gramchange_mod")
summary(pers.gramchange.mod) # temp negative, prec ns, gram positive, plot size ns, richness positive, duration ns
print(summary(pers.gramchange.mod, prob = 0.975)) # temp negative, prec ns, gram positive, plot size ns, richness positive, duration ns





## EXTINCTS ----


#### 01) GEOGRAPHICAL MODEL ####
ext.geo.mod <- brm(bf(ProportionSpecies ~ LAT + Region + LogPlotSize + 
                      MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                    data = ext, iter = 2000, chains = 4, warmup = 500, 
                    family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi),
                    file = "models/22_01_ext_geo_mod")


summary(ext.geo.mod) # latitude ns, plot size ns, richness negative, duration ns
conditional_effects(ext.geo.mod) # region ns

print(summary(ext.geo.mod, prob = 0.975)) # latitude ns, plot size ns, richness negative, duration ns
conditional_effects(ext.geo.mod, prob = 0.975) # region ns



#### 02) CLIMATIC MODEL ####

# Join turnover with climate data
ext_clim <- left_join(ext, clim.only, by = "SiteSubsitePlot")

# Clean up variables
ext_clim2 <- ext_clim %>% filter(MOISTURE != "")

# Extincts model 
ext.clim.mod <- brm(bf(ProportionSpecies ~ MOISTURE + warmq + prec + Duration + 
                       LogPlotSize + MeanRichness + (1|SiteSubsite), zi ~ 1),
                     family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                     data = ext_clim2, iter = 2000, chains = 4, warmup = 500, init = 0,
                     file = "models/22_02_ext_clim_mod")

print(summary(ext.clim.mod), digits = 5) # temp negative, prec positive, duration ns, plot size ns, richness negative
conditional_effects(ext.clim.mod) # moisture ns

print(summary(ext.clim.mod, prob = 0.975), digits = 5) # temp negative, prec positive, duration ns, plot size ns, richness negative
conditional_effects(ext.clim.mod, prob = 0.975) # moisture ns



# predictions for climate figure
ext.clim.predX <- ggpredict(ext.clim.mod, terms = "warmq")
colnames(ext.clim.predX) = c('warmq', 'fit', 'lwr', 'upr')

ext.clim.predXX <- ggpredict(ext.clim.mod, terms = "prec")
colnames(ext.clim.predXX) = c('prec', 'fit', 'lwr', 'upr')



## Univariate Extinct ~ temp model
ext.clim.uni.mod <- brm(bf(ProportionSpecies ~ warmq + (1|SiteSubsite), zi ~ 1), init = 0,
                        family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                        data = ext_clim2, iter = 3500, chains = 2, warmup = 500, 
                        file = "models/22_02_ext_clim_univ_mod")

summary(ext.clim.uni.mod) # negative significant
conditional_effects(ext.clim.uni.mod)

# predictions
ext.clim.uni.pred <- ggpredict(ext.clim.uni.mod, terms = "warmq", type = "random",
                               condition = c(SiteSubsite = 1266), allow_new_levels = TRUE)
colnames(ext.clim.uni.pred) = c('warmq', 'fit', 'lwr', 'upr', 'group')


## Univariate Extinct ~ precip model
ext.clim.uni.mod2 <- brm(bf(ProportionSpecies ~ prec + (1|SiteSubsite), zi ~ 1), init = 0,
                        family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                        data = ext_clim2, iter = 3500, chains = 2, warmup = 500, 
                        file = "models/22_02_ext_clim_univ_mod2")

print(summary(ext.clim.uni.mod2), digits = 4) # positive ns
conditional_effects(ext.clim.uni.mod2)





#### 03) FG COMPOSITION MODEL ####

# shrub
ext.fg.shb.mod <- brm(bf(ProportionSpecies ~ ShrubDivided + 
                    LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                  family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                  data = ext, iter = 2000, chains = 4, warmup = 500, init = 0,
                  file = "models/22_03_ext_fg_shb_mod")

summary(ext.fg.shb.mod) # shrub ns, plot size ns, richness negative, duration ns
print(summary(ext.fg.shb.mod, prob = 0.975)) # shrub ns, plot size ns, richness negative, duration ns


# forb
ext.fg.forb.mod <- brm(bf(ProportionSpecies ~ ForbDivided + 
                           LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                      family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                      data = ext, iter = 2000, chains = 4, warmup = 500, init = 0,
                      file = "models/22_03_ext_fg_forb_mod")

summary(ext.fg.forb.mod) # forb ns, plot size ns, richness negative, duration ns
print(summary(ext.fg.forb.mod, prob = 0.975)) # forb ns, plot size ns, richness negative, duration ns


# gram
ext.fg.gram.mod <- brm(bf(ProportionSpecies ~ GramDivided + 
                            LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                       family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                       data = ext, iter = 2000, chains = 4, warmup = 500, init = 0,
                       file = "models/22_03_ext_fg_gram_mod")

summary(ext.fg.gram.mod) # gram ns, plot size ns, richness negative, duration ns
print(summary(ext.fg.gram.mod, prob = 0.975)) # gram ns, plot size ns, richness negative, duration ns




#### 04) CHANGE OCCURRED MODEL ####

ext.tur.change0 <- left_join(ext, change_clim, by = "SiteSubsitePlot")
ext.tur.change <- left_join(ext.tur.change0, change_fg, by = "SiteSubsitePlot")


# One by one - shrub
ext.shbchange.mod <- brm(bf(ProportionSpecies ~ WarmQSlope + PrecSlope + 
                           ShrubSlope + LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                         family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                         data = ext.tur.change, iter = 2000, chains = 4, warmup = 500, 
                         file = "models/22_04_ext_shbchange_mod")
summary(ext.shbchange.mod) # temp positive, prec positive, shrub positive, plot size ns, richness neg, duration neg
print(summary(ext.shbchange.mod, prob = 0.975)) # temp positive, prec ns, shrub positive, plot size ns, richness neg, duration ns


# univariate model
ext.shbchange.uni.mod <- brm(bf(ProportionSpecies ~ ShrubSlope + (1|SiteSubsite), zi ~ 1),
                             family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                             data = ext.tur.change, iter = 3500, chains = 2, warmup = 500, 
                             file = "models/22_04_ext_shbchange_univ_mod")

summary(ext.shbchange.uni.mod) # shrub positive significant
conditional_effects(ext.shbchange.uni.mod)


# model predictions
ext.shbchange.uni.mod.pred <- ggpredict(ext.shbchange.uni.mod, terms = "ShrubSlope", type = "random",
                                        allow_new_levels = TRUE, condition = c(SiteSubsite = 1266))
colnames(ext.shbchange.uni.mod.pred) = c('ShrubSlope', 'fit', 'lwr', 'upr', 'group')



## now without outliers 
ext.tur.change.outliers <- ext.tur.change %>% 
  drop_na(ShrubSlope) %>% 
  mutate(ShrubChangeSD = sd(ShrubSlope), ShrubChangeSD2 = ShrubChangeSD*2, 
         ShrubChangeSD3 = ShrubChangeSD*3, ShrubChangeSD3Neg = -(ShrubChangeSD3)) %>% 
  filter(ShrubSlope < ShrubChangeSD3) %>% # this gets rid of the positive outliers but not the negative ones
  filter(ShrubSlope > ShrubChangeSD3Neg)


## model without outliers
ext.shbchange.uni.out.mod <- brm(bf(ProportionSpecies ~ ShrubSlope + (1|SiteSubsite), zi ~ 1),
                             family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                             data = ext.tur.change.outliers, iter = 3500, chains = 2, warmup = 500, 
                             file = "models/22_04_ext_shbchange_univ_out_mod")

summary(ext.shbchange.uni.out.mod) # shrub positive significant
conditional_effects(ext.shbchange.uni.out.mod)


# model predictions
ext.shbchange.uni.out.mod.pred <- ggpredict(ext.shbchange.uni.out.mod, terms = "ShrubSlope", 
                                            type = "random",
                                            allow_new_levels = TRUE, condition = c(SiteSubsite = 1266))
colnames(ext.shbchange.uni.out.mod.pred) = c('ShrubSlope', 'fit', 'lwr', 'upr', 'group')


# plot
(ext.shrubchg.out.plot <- ggplot(ext.tur.change.outliers) + 
    geom_point(aes(x = ShrubSlope, y = ProportionSpecies*100), size = 5, colour = "#d68317", alpha = 0.5) +
    geom_line(data = ext.shbchange.uni.out.mod.pred, aes(x = ShrubSlope, y = fit*100), colour = "#9c5b06", linewidth = 1.2) + 
    geom_ribbon(data = ext.shbchange.uni.out.mod.pred, aes(x = ShrubSlope, ymin = lwr*100, ymax = upr*100), fill = "#9c5b06", alpha = 0.2) +
    xlab("\nShrub cover change\n(% per year, slopes)") +
    ylab("Species losses (%)\n") + 
    bio.theme + 
    theme(axis.title.x = element_text(face="plain", size=20), 
          axis.text.x  = element_text(vjust=0.5, size=18, colour = "black"), 
          axis.title.y = element_text(face="plain", size=20),
          axis.text.y  = element_text(vjust=0.5, size=18, colour = "black"), 
          legend.key.size = unit(1, 'cm')))



# One by one - forb
ext.forbchange.mod <- brm(bf(ProportionSpecies ~ WarmQSlope + PrecSlope + 
                            ForbSlope + LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                          family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                          data = ext.tur.change, iter = 2000, chains = 4, warmup = 500, 
                          file = "models/22_04_ext_forbchange_mod")
print(summary(ext.forbchange.mod), digits = 5) # temp positive, prec positive, forb ns, plot size ns, richness negative, duration ns
print(summary(ext.forbchange.mod, prob = 0.975)) # temp positive, prec positive, forb ns, plot size ns, richness negative, duration ns


# One by one - gram
ext.gramchange.mod <- brm(bf(ProportionSpecies ~ WarmQSlope + PrecSlope + 
                            GraminoidSlope + LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                          family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                          data = ext.tur.change, iter = 2000, chains = 4, warmup = 500, 
                          file = "models/22_04_ext_gramchange_mod")
summary(ext.gramchange.mod) # temp positive, prec positive, gram negative, plot size ns, richness negative, duration ns
print(summary(ext.gramchange.mod, prob = 0.975)) # temp positive, prec ns, gram ns, plot size ns, richness negative, duration ns



## univariate with temp change

# temp change mod only
ext.tempchange.uni.mod <- brm(bf(ProportionSpecies ~ WarmQSlope + (1|SiteSubsite), zi ~ 1),
                          family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                          data = ext.tur.change, iter = 2000, chains = 4, warmup = 500, 
                          file = "models/22_04_ext_tempchange_univ_mod")

summary(ext.tempchange.uni.mod) # positive significant


# predictions
ext.tempchange.uni.pred <- ggpredict(ext.tempchange.uni.mod, terms = "WarmQSlope",
                                     type = "random", condition = c(SiteSubsite = 1266),
                                     allow_new_levels = TRUE)
colnames(ext.tempchange.uni.pred) = c('WarmQSlope', 'fit', 'lwr', 'upr', 'group')





## COLONISERS ----


#### 01) GEOGRAPHICAL MODEL ####
col.geo.mod <- brm(bf(ProportionSpecies ~ LAT + Region + LogPlotSize + 
                      MeanRichness + Duration + (1|SiteSubsite), zi ~ 1), 
                    data = col, iter = 2000, chains = 4, warmup = 500, 
                    family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi),
                    file = "models/22_01_col_geo_mod")

summary(col.geo.mod) # lat ns, plot size ns, richness negative, duration ns
conditional_effects(col.geo.mod) # Greenland greater then Eurasia and NAm-west

print(summary(col.geo.mod, prob = 0.975)) # lat ns, plot size ns, richness negative, duration ns
conditional_effects(col.geo.mod, prob = 0.975) # region ns




#### 02) CLIMATIC MODEL ####

# Join turnover with climate data
col_clim <- left_join(col, clim.only, by = "SiteSubsitePlot")

# Clean up variables
col_clim2 <- col_clim %>% filter(MOISTURE != "")

# Model
col.clim.mod <- brm(bf(ProportionSpecies ~ MOISTURE + warmq + prec + Duration + 
                       LogPlotSize + MeanRichness + (1|SiteSubsite), zi ~ 1), 
                     family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                     data = col_clim2, iter = 2000, chains = 4, warmup = 500, 
                     file = "models/22_02_col_clim_mod")

print(summary(col.clim.mod), digits = 5) # temp negative, prec ns, duration ns, plot size ns, richness negative
conditional_effects(col.clim.mod) # moisture ns

print(summary(col.clim.mod, prob = 0.975)) # temp ns, prec ns, duration ns, plot size ns, richness negative
conditional_effects(col.clim.mod, prob = 0.975) # moisture ns



# predictions for climate figure
col.clim.predXX <- ggpredict(col.clim.mod, terms = "warmq")
colnames(col.clim.predXX) = c('warmq', 'fit', 'lwr', 'upr')

col.clim.predXXX <- ggpredict(col.clim.mod, terms = "prec")
colnames(col.clim.predXXX) = c('prec', 'fit', 'lwr', 'upr')




## Univariate colonisers ~ temp model
col.clim.uni.mod <- brm(bf(ProportionSpecies ~ warmq + (1|SiteSubsite), zi ~ 1), init = 0,
                        family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                        data = col_clim2, iter = 3500, chains = 2, warmup = 500, 
                        file = "models/22_02_col_clim_univ_mod")

summary(col.clim.uni.mod) # negative significant
conditional_effects(col.clim.uni.mod)

# predictions
col.clim.uni.pred <- ggpredict(col.clim.uni.mod, terms = "warmq", type = "random",
                               condition = c(SiteSubsite = 1266), allow_new_levels = TRUE)
colnames(col.clim.uni.pred) = c('warmq', 'fit', 'lwr', 'upr', 'group')


# Plot both trends in one panel
col.clim.uni.pred.clean <- col.clim.uni.pred %>% dplyr::select(!5) %>% mutate(Trend = "Colonisers")
ext.clim.uni.pred.clean <- ext.clim.uni.pred %>% dplyr::select(!5) %>% mutate(Trend = "Extinct")


# Temp graph - with points
(extcol.temp.points <- ggplot() +
    geom_point(data = col_clim2, aes(x=warmq, y = ProportionSpecies*100), size = 4, alpha = 0.2, shape = 16, colour = "#8738C7") + 
    geom_point(data = ext_clim2, aes(x=warmq, y = ProportionSpecies*100), size = 4, alpha = 0.2, shape = 16, colour = "#d68317") + 
    geom_line(data = ext.clim.uni.pred.clean, aes(y = fit*100, x = warmq), colour = "#d68317", size = 2) +
    geom_ribbon(data = ext.clim.uni.pred.clean, aes(x = warmq, ymin = lwr*100, ymax = upr*100), fill = "#d68317", alpha = 0.2) +
    geom_line(data = col.clim.uni.pred.clean, aes(y = fit*100, x = warmq), colour = "#8738C7",size = 2) + 
    geom_ribbon(data = col.clim.uni.pred.clean, aes(x = warmq, ymin = lwr*100, ymax = upr*100), fill = "#8738C7", alpha = 0.2) +
    xlab("\nWarmest Quarter Temperature (°C)") + ylab ("Species per trajectory (%)\n") + bio.theme +
    ylim(0, 40) +
    theme(axis.title.x = element_text(face="plain", size=24),
          axis.text.x  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.title.y = element_text(face="plain", size=24),
          axis.text.y  = element_text(vjust=0.5, size=22, colour = "black")))

ggsave(extcol.temp.points, filename = "figures/Figure_3c.png", 
       width = 20, height = 20, units = "cm")




#### 03) FG COMPOSITION MODEL ####

# shrub
col.fg.shb.mod <- brm(bf(ProportionSpecies ~ ShrubDivided + 
                    LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                  family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                  data = col, iter = 2000, chains = 4, warmup = 500, 
                  file = "models/22_03_col_fg_shb_mod")

print(summary(col.fg.shb.mod), digits = 5) # shrub ns, plot size ns, richness negative, duration positive
print(summary(col.fg.shb.mod, prob = 0.975)) # shrub ns, plot size ns, richness negative, duration ns



# forb
col.fg.forb.mod <- brm(bf(ProportionSpecies ~ ForbDivided + 
                        LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                   family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                   data = col, iter = 2000, chains = 4, warmup = 500, 
                   file = "models/22_03_col_fg_forb_mod")

summary(col.fg.forb.mod) # forb ns, plot size ns, richness negative, duration ns
print(summary(col.fg.forb.mod, prob = 0.975)) # forb ns, plot size ns, richness negative, duration ns


# gram
col.fg.gram.mod <- brm(bf(ProportionSpecies ~ GramDivided + 
                        LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                   family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                   data = col, iter = 2000, chains = 4, warmup = 500, 
                   file = "models/22_03_col_fg_gram_mod")

print(summary(col.fg.gram.mod), digits = 5) # graminoid ns, plot size ns, richness negative, duration ns
print(summary(col.fg.gram.mod, prob = 0.975)) # graminoid ns, plot size ns, richness negative, duration ns




#### 04) CHANGE OCCURRED MODEL ####

col.tur.change0 <- left_join(col, change_clim, by = "SiteSubsitePlot")
col.tur.change <- left_join(col.tur.change0, change_fg, by = "SiteSubsitePlot")


# Model using change in covariates
col.shbchange.mod <- brm(bf(ProportionSpecies ~ WarmQSlope + PrecSlope + 
                           ShrubSlope + LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                         family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                         data = col.tur.change, iter = 2000, chains = 4, warmup = 500, 
                         file = "models/22_04_col_shbchange_mod")
summary(col.shbchange.mod) # temp positive, prec ns, shrub ns, plot size ns, richness negative, duration ns
print(summary(col.shbchange.mod, prob = 0.975), digits = 5) # temp positive, prec ns, shrub ns, plot size ns, richness negative, duration ns


## univariate model
col.shbchange.uni.mod <- brm(bf(ProportionSpecies ~ ShrubSlope + (1|SiteSubsite), zi ~ 1),
                             family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                             data = col.tur.change, iter = 3500, chains = 2, warmup = 500, 
                             file = "models/22_04_col_shbchange_univ_mod")

summary(col.shbchange.uni.mod) # shrub negative ns


# model predictions
col.shbchange.uni.mod.pred <- ggpredict(col.shbchange.uni.mod, terms = "ShrubSlope", type = "random",
                                        condition = c(SiteSubsite = 1266), allow_new_levels = TRUE)
colnames(col.shbchange.uni.mod.pred) = c('ShrubSlope', 'fit', 'lwr', 'upr', 'group')


# Combine col/ext into one
(col.shrubchg.plot <- ggplot() + 
    geom_point(data = col.tur.change, aes(x = ShrubSlope, y = ProportionSpecies*100), size = 4, colour = "#8738C7", alpha = 0.2, shape = 16) +
    geom_point(data = ext.tur.change, aes(x = ShrubSlope, y = ProportionSpecies*100), size = 4, colour = "#d68317", alpha = 0.2, shape = 16) +
    geom_line(data = col.shbchange.uni.mod.pred, aes(x = ShrubSlope, y = fit*100), colour = "#410373", linetype = "dashed", linewidth = 2) + 
    geom_ribbon(data = col.shbchange.uni.mod.pred, aes(x = ShrubSlope, ymin = lwr*100, ymax = upr*100), fill = "#410373", alpha = 0.2) +
    geom_line(data = ext.shbchange.uni.mod.pred, aes(x = ShrubSlope, y = fit*100), colour = "#9c5b06", linewidth = 2) + 
    geom_ribbon(data = ext.shbchange.uni.mod.pred, aes(x = ShrubSlope, ymin = lwr*100, ymax = upr*100), fill = "#9c5b06", alpha = 0.2) +
    xlab("\nShrub cover change (% per year)") +
    ylab("Species gains and losses (%)\n") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    bio.theme + 
    theme(axis.title.x = element_text(face="plain", size=24), 
          axis.text.x  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.title.y = element_text(face="plain", size=24),
          axis.text.y  = element_text(vjust=0.5, size=22, colour = "black"), 
          legend.key.size = unit(1, 'cm')))


ggplot2::ggsave("figures/Figure_3e.png", col.shrubchg.plot, width = 20, height = 20, units = "cm")




# Model with plot size == 1 only
col.tur.change.1m <- col.tur.change %>% filter(SurveyedArea == 1)

col.shbchange.1m.mod <- brm(bf(ProportionSpecies ~ WarmQSlope + PrecSlope + 
                              ShrubSlope + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                         family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                         data = col.tur.change.1m, iter = 2000, chains = 4, warmup = 500, 
                         file = "models/22_04b_col_shbchange_1m_mod")
summary(col.shbchange.1m.mod) # temp positive, prec ns, shrub ns, plot size ns, richness negative, duration ns



# forb
col.forbchange.mod <- brm(bf(ProportionSpecies ~ WarmQSlope + PrecSlope + 
                            ForbSlope + LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                          family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                          data = col.tur.change, iter = 2000, chains = 4, warmup = 500, 
                          file = "models/22_04_col_forbchange_mod")
summary(col.forbchange.mod) # temp positive, prec ns, forb positive, plot size ns, richness negative, duration ns
print(summary(col.forbchange.mod, prob = 0.975)) # temp positive, prec ns, forb positive, plot size ns, richness negative, duration ns



# gram
col.gramchange.mod <- brm(bf(ProportionSpecies ~ WarmQSlope + PrecSlope + 
                            GraminoidSlope + LogPlotSize + MeanRichness + Duration + (1|SiteSubsite), zi ~ 1),
                          family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                          data = col.tur.change, iter = 2000, chains = 4, warmup = 500, 
                          file = "models/22_04_col_gramchange_mod")
summary(col.gramchange.mod) # temp positive, prec ns, gram ns, plot size ns, richness negative, duration ns
print(summary(col.gramchange.mod, prob = 0.975)) # temp positive, prec ns, gram ns, plot size ns, richness negative, duration ns





## univariate with temp change

# temp change mod only
col.tempchange.uni.mod <- brm(bf(ProportionSpecies ~ WarmQSlope + (1|SiteSubsite), zi ~ 1),
                              family = "zero_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                              data = col.tur.change, iter = 2000, chains = 4, warmup = 500, 
                              file = "models/22_04_col_tempchange_univ_mod")

summary(col.tempchange.uni.mod) # positive significant


# predictions
col.tempchange.uni.pred <- ggpredict(col.tempchange.uni.mod, terms = "WarmQSlope",
                                     type = "random", condition = c(SiteSubsite = 1266),
                                     allow_new_levels = TRUE)
colnames(col.tempchange.uni.pred) = c('WarmQSlope', 'fit', 'lwr', 'upr', 'group')


# Plot both trends in one panel
col.tempchange.uni.pred.clean <- col.tempchange.uni.pred %>% dplyr::select(!5) %>% mutate(Trend = "Colonisers")
ext.tempchange.uni.pred.clean <- ext.tempchange.uni.pred %>% dplyr::select(!5) %>% mutate(Trend = "Extinct")


# Temp graph
(extcol.temp.change <- ggplot() +
    geom_point(data = col.tur.change, aes(x=WarmQSlope, y = ProportionSpecies*100), size = 4, alpha = 0.2, shape = 16, colour = "#8738C7") + 
    geom_point(data = ext.tur.change, aes(x=WarmQSlope, y = ProportionSpecies*100), size = 4, alpha = 0.2, shape = 16, colour = "#d68317") + 
    geom_line(data = ext.tempchange.uni.pred.clean, aes(y = fit*100, x = WarmQSlope), colour = "#d68317", size = 2) +
    geom_ribbon(data = ext.tempchange.uni.pred.clean, aes(x = WarmQSlope, ymin = lwr*100, ymax = upr*100), fill = "#d68317", alpha = 0.2) +
    geom_line(data = col.tempchange.uni.pred.clean, aes(y = fit*100, x = WarmQSlope), colour = "#8738C7",size = 2) + 
    geom_ribbon(data = col.tempchange.uni.pred.clean, aes(x = WarmQSlope, ymin = lwr*100, ymax = upr*100), fill = "#8738C7", alpha = 0.2) +
    ylim(0, 40) +
    xlab("\nTemperature change (°C per year)") + ylab ("Species per trajectory (%)\n") + bio.theme +
    theme(axis.title.x = element_text(face="plain", size=24),
          axis.text.x  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.title.y = element_text(face="plain", size=24),
          axis.text.y  = element_text(vjust=0.5, size=22, colour = "black")))

ggsave(extcol.temp.change, filename = "figures/Figure_3d.png", 
       width = 20, height = 20, units = "cm")




## SPECIES POOL PANEL ----

## Persisters model predictions
pers.geo.pred <- ggpredict(pers.geo.mod, term = "MeanRichness", type = "random", 
                           condition = c(SiteSubsite = 1266),
                           allow_new_levels = TRUE)

colnames(pers.geo.pred) = c('MeanRichness', 'fit', 'lwr', 'upr', 'group')

# Plot predictions with x = species pool
(pers.geo.pool.plot <- ggplot() + 
    geom_point(data = pers, aes(x = MeanRichness, y = ProportionSpecies*100), size = 2, alpha = 0.5, colour = "#38c787") +
    xlab("\nSpecies richness (plot mean)") + ylab ("Persisting species (%)\n") + bio.theme + 
    geom_line(data = pers.geo.pred, aes(x = MeanRichness, y = fit*100), colour = "black",  size = 1) + 
    geom_ribbon(data = pers.geo.pred, aes(x = MeanRichness, ymin = lwr*100, ymax = upr*100), fill = "black", alpha = 0.2) +
    theme(legend.key = element_blank(), legend.background = element_blank(), 
          axis.title.x = element_text(face="plain", size=9),
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y  = element_text(vjust=0.5, size=9, colour = "black"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines")))



# Extinct model predictions
ext.geo.pred <- ggpredict(ext.geo.mod, term = "MeanRichness", type = "random", 
                          condition = c(SiteSubsite = 1266),
                          allow_new_levels = TRUE)
colnames(ext.geo.pred) = c('MeanRichness', 'fit', 'lwr', 'upr', 'group')

# Plot predictions with x = species pool
(ext.geo.pool.plot <- ggplot() + 
    geom_point(data = ext, aes(x = MeanRichness, y = ProportionSpecies*100), size = 2, alpha = 0.5, colour = "#d68317") +
    xlab("\nSpecies richness (plot mean)") + ylab ("Species losses (%)\n") + bio.theme + 
    ylim(0, 100) +
    geom_line(data = ext.geo.pred, aes(x = MeanRichness, y = fit*100), colour = "black",  size = 1) + 
    geom_ribbon(data = ext.geo.pred, aes(x = MeanRichness, ymin = lwr*100, ymax = upr*100), fill = "black", alpha = 0.2) +
    theme(legend.key = element_blank(), legend.background = element_blank(), 
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),
          axis.title.x = element_text(face="plain", size=9),
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y  = element_text(vjust=0.5, size=9, colour = "black")))




# Coloniser model predictions
col.geo.pred <- ggpredict(col.geo.mod, term = "MeanRichness", type = "random", 
                          condition = c(SiteSubsite = 1266),
                          allow_new_levels = TRUE)
colnames(col.geo.pred) = c('MeanRichness', 'fit', 'lwr', 'upr', 'group')

# Plot predictions with x = species pool
(col.geo.pool.plot <- ggplot() + 
    geom_point(data = col, aes(x = MeanRichness, y = ProportionSpecies*100), size = 2, alpha = 0.5, colour = "#8738C7") +
    xlab("\nSpecies richness (plot mean)") + ylab ("Species gains (%)\n") + bio.theme + 
    ylim(0, 100) +
    geom_line(data = col.geo.pred, aes(x = MeanRichness, y = fit*100), colour = "black",  size = 1) + 
    geom_ribbon(data = col.geo.pred, aes(x = MeanRichness, ymin = lwr*100, ymax = upr*100), fill = "black", alpha = 0.2) +
    theme(legend.key = element_blank(), legend.background = element_blank(), 
          axis.title.x = element_text(face="plain", size=9),
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y  = element_text(vjust=0.5, size=9, colour = "black"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines")))


## Panel: proportion of species vs species pool + prop species vs evenness
(prop.pool.panel <- plot_grid(ext.geo.pool.plot, pers.geo.pool.plot, col.geo.pool.plot, 
                              ncol=3, nrow = 1, align="hv", label_size = 14, 
                              labels = c("a", "b", "c")))




## RARE SPECIES ANALYSIS ----

# species losses
losses <- read.csv("data/22top_extinct.csv") 

# compared with how many times are they found across ITEX (and across how mnay sites)
load("data/itex_dec22.RData")

losses.sp <- itex.dec22 %>% group_by(SPECIES_NAME) %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>% 
  mutate(SpeciesInPlots = length(unique(SiteSubsitePlot))) %>%
  ungroup() %>%
  distinct(SPECIES_NAME, .keep_all = TRUE) %>%
  select(SPECIES_NAME, SpeciesInPlots)

losses.sp.sites <- itex.dec22 %>% group_by(SPECIES_NAME) %>% 
  distinct(SITE, .keep_all = TRUE) %>% 
  mutate(SpeciesInSites = length(unique(SITE))) %>%
  ungroup() %>%
  distinct(SPECIES_NAME, .keep_all = TRUE) %>%
  select(SPECIES_NAME, SpeciesInSites)


# join with losses database
losses.rare <- left_join(losses, losses.sp.sites, by = "SPECIES_NAME") %>% 
  mutate(PropInTrendDiv = PropInTrend/100,
         PropInTrendDivOffset = PropInTrendDiv - 0.00001)

hist(losses.rare$TimesInTrend) #negbin


# Are species that have been lost multiple times those that are more rare (i.e., present at fewer sites)?
str(losses.rare)

rare.mod.prop <- brm(bf(PropInTrendDivOffset ~ SpeciesInSites),
                     data = losses.rare, 
                     family = "beta", 
                     iter = 2500, chains = 4, warmup = 400, 
                     file = "models/losses_rarity_prop_mod")

summary(rare.mod.prop)
conditional_effects(rare.mod.prop)
bayes_R2(rare.mod.prop)
bayes_R2(rare.mod.prop, re.form = NA)


