## Plant diversity dynamics over space and time in a warming Arctic
## Mariana Garcia Criado (mariana.garcia.criado@gmail.com)
## Script 3. Species turnover 

# Thank you to Isla Myers-Smith and Dani Gargya for the initial loop code! 



## LIBRARIES ----
library(tidyverse)
library(vegan)
library(betapart)
library(viridis)
library(brms)
library(cowplot)
library(ggeffects)




## THEMES ----
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


## DATA ----
load("data/itex_dec22.RData")
str(itex.dec22$RelCover)


# Filter for those plots with at least two time points
bio_turnover <- itex.dec22 %>%
  dplyr::select(SiteSubsitePlot, YEAR, SPECIES_NAME, RelCover) %>%
  group_by(SiteSubsitePlot) %>%
  filter(length(unique(YEAR)) > 1) %>% 
  mutate(Duration = max(YEAR) - min(YEAR)) %>% filter(Duration > 4)
  filter(RelCover >= 0) %>%
  ungroup()

write.csv(bio_turnover, "data/22bio_turnover.csv")




## TURNOVER ----

# Create empty dataframe
beta_Jaccard <- data.frame(matrix(ncol = 10, nrow = length(unique(bio_turnover$SiteSubsitePlot)))) 
names(beta_Jaccard) <- c("SiteSubsitePlot", "duration", "richness", "richness_change", 
                         "Jbeta", "Jtu", "Jne", "Bbal", "Bgra", "Bbray") 
i = 1

# Loop to calculate turnover with betapart
for (i in 1:length(unique(bio_turnover$SiteSubsitePlot))) {
  SiteSubsitePlotName <- as.character(unique(bio_turnover$SiteSubsitePlot)[i])
  sub_bio_abundance <- filter(bio_turnover, SiteSubsitePlot == SiteSubsitePlotName)
  YearMin <- min(sub_bio_abundance$YEAR)
  YearMax <- max(sub_bio_abundance$YEAR)
  duration <- YearMax - YearMin
  # filters the dataframe for just the first and last observations per plot
  sub_bio_abundance_min <- filter(sub_bio_abundance, YEAR == YearMin)
  sub_bio_abundance_max <- filter(sub_bio_abundance, YEAR == YearMax)
  sub_bio_abundance <- rbind(sub_bio_abundance_min, sub_bio_abundance_max)
  richness <- length(unique(sub_bio_abundance$SPECIES_NAME))
  # averages any species that have multiple records per plot and time point 
  sub_bio_abundance <- sub_bio_abundance %>% group_by(SiteSubsitePlot, SPECIES_NAME, YEAR) %>% summarise(MeanRelCover = mean(RelCover)) %>% ungroup()
  # reshape to wide form
  sub_bio_abundance_wider <- tidyr::pivot_wider(sub_bio_abundance, names_from = SPECIES_NAME, 
                                         values_from = MeanRelCover, 
                                         values_fill = list(MeanRelCover = 0))
  # removes columns for beta.pair() function
  sub_bio_abundance_matrix <- dplyr::select(sub_bio_abundance_wider, -SiteSubsitePlot, -YEAR) 
  # creates presence-absence matrix
  sub_bio_presence_matrix <- with(sub_bio_abundance_matrix, ifelse(sub_bio_abundance_matrix > 0,1,0))
  # calculates Jaccard overall, turnover and nestedness
  J_components <- beta.pair(sub_bio_presence_matrix, index.family='jaccard')
  B_components <- beta.pair.abund(sub_bio_abundance_matrix, index.family='bray')
  # saves biodiversity metrics
  richness_change <- rowSums(sub_bio_presence_matrix)[2] - rowSums(sub_bio_presence_matrix)[1]
  Jbeta <- J_components$beta.jac
  Jtu <- J_components$beta.jtu
  Jne <- J_components$beta.jne
  Bbal <- B_components$beta.bray.bal
  Bgra <- B_components$beta.bray.gra
  Bbray <- B_components$beta.bray
  beta_Jaccard[i,] <- c(SiteSubsitePlotName, duration, richness, richness_change, Jbeta, Jtu, Jne, Bbal, Bgra, Bbray)
  i = i+1
}


# Convert in adequate format 
beta_Jaccard_final <- beta_Jaccard %>% 
  mutate(duration = as.numeric(duration), richness = as.numeric(richness), richness_change = as.numeric(richness_change),
         richness_change_abs = abs(as.numeric(richness_change)), Jtu = as.numeric(Jtu), Bbray = as.numeric(Bbray)) %>%
  dplyr::select(., -c(Jbeta, Jne, Bbal, Bgra)) %>% rename(Jaccard = Jtu, BrayCurtis = Bbray)

# All makes sense
glimpse(beta_Jaccard_final)

# Save dataframe
write.csv(beta_Jaccard_final, "data/22beta_jaccard_final.csv")


# beta.pair function computes turnover, nestedness and dissimilarity (Jaccard/Sorensen)
# index.family can be sorensen or jaccard. For index.family="jaccard" the three matrices are:
# beta.jtu: dissimilarity matrix accounting for spatial turnover, measured as the turnover-fraction of Jaccard pair-wise dissimilarity
# beta.jne: dissimilarity matrix accounting for nestedness-resultant dissimilarity, measured as the nestedness-fraction of Jaccard pair-wise dissimilarity
# beta.jac: dissimilarity matrix accounting for beta diversity, measured as Jaccard pair-wise dissimilarity (a monotonic transformation of beta diversity)



## JACCARD MODELS ----
beta_Jaccard_final <- read.csv("data/22beta_jaccard_final.csv")

# Examine data distributions
hist(beta_Jaccard_final$Jaccard) #[0 - 1], lots of 0s, just 1 value of 1

# Extract covariates for the models
add.covs <- itex.dec22 %>% distinct(SiteSubsitePlot, .keep_all = TRUE) %>%
  dplyr::select(SiteSubsite, SiteSubsitePlot, MOISTURE, LAT, LONG, SurveyedArea, 
                Region, PlotShrubMean, PlotGraminoidMean, PlotForbMean, MeanRichness)

# Merge the 2 dataframes
beta_turnover <- left_join(beta_Jaccard_final, add.covs, by = "SiteSubsitePlot")

# Plot size is not numerical
beta_turnover$SurveyedArea <- as.numeric(beta_turnover$SurveyedArea)

write.csv(beta_turnover, "data/22beta_turnover.csv")



#### 01) GEOGRAPHICAL MODEL ####
beta_turnover <- read.csv("data/22beta_turnover.csv")


# Prepare variables for model
beta_turnover2 <- beta_turnover %>% 
  mutate(LogPlotSize = log(SurveyedArea)) %>% 
  mutate(ShrubDivided = PlotShrubMean/100, 
         ForbDivided = PlotForbMean/100, 
         GramDivided = PlotGraminoidMean/100)


# GEO mod with fixed zoi/coi
jac.geo.mod <- brm(bf(Jaccard ~ LAT + Region + LogPlotSize + 
                      MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                    data = beta_turnover2, iter = 2000, chains = 4, warmup = 500, 
                    family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi),
                    file = "models/24_01_jaccard_geo_mod")

summary(jac.geo.mod) # latitude ns, plot size ns, richness negative, duration ns, greater turnover in GreenIceland than in NAM-west
conditional_effects(jac.geo.mod)

summary(jac.geo.mod, prob = 0.975) # latitude ns, plot size ns, richness negative, duration ns, region ns
conditional_effects(jac.geo.mod, prob = 0.975)


# ED Fig 8g: Jaccard ~ duration plot 

# predictions
jac.geo.pred <- ggpredict(jac.geo.mod, terms = "duration", type = "random", 
                          condition = c(SiteSubsite = 1266),
                          allow_new_levels = TRUE)
colnames(jac.geo.pred) = c('duration', 'fit', 'lwr', 'upr', 'group')


# plot relationship
(jac.geo.pred.plot <- ggplot() + 
    geom_point(data = beta_turnover2, aes(x=duration, y = Jaccard), 
               size = 2, pch = 21, colour = "black", fill = "#54C9AD70") + 
    xlab("\nDuration (years)") + ylab("Jaccard\n") +
    geom_line(data = jac.geo.pred, aes(x = duration, y = fit), linewidth = 1, colour = "#54C9AD70", linetype = "dashed") +
    geom_ribbon(data = jac.geo.pred, aes(x = duration, ymin = lwr, ymax = upr), 
                fill = "#54C9AD70", alpha = 0.3) + 
    bio.theme +
    theme(axis.title.x = element_text(face="plain", size=9),
  axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
  axis.title.y = element_text(face="plain", size=9),
  axis.text.y  = element_text(vjust=0.5, size=9, colour = "black"),
  plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines")))



#### 02) CLIMATE MODEL ####

# Bring in climate data
clim.only <- read.csv("data/22clim_only.csv")

# Join turnover with climate data
turn_clim <- left_join(beta_turnover2, clim.only, by = "SiteSubsitePlot")

# Prepare variables
turn_clim2 <- turn_clim %>% filter(MOISTURE != "")

# CLIM mod with fixed zoi and coi
jac.clim.mod <- brm(bf(Jaccard ~ MOISTURE + warmq + prec + duration + 
                       LogPlotSize + MeanRichness + (1|SiteSubsite), zoi ~ 1, coi ~ 1), init = 0,
                     family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                     data = turn_clim2, iter = 2000, chains = 4, warmup = 500, 
                     file = "models/24_02_jaccard_clim_mod")

summary(jac.clim.mod) # temp negative, prec positive, duration ns, plot size ns, richness negative, moisture ns
print(summary(jac.clim.mod), digits = 5)
conditional_effects(jac.clim.mod)

summary(jac.clim.mod, prob = 0.975) # temp negative, prec ns, duration ns, plot size ns, richness negative, moisture ns
conditional_effects(jac.clim.mod, prob = 0.975)


## Univariate Jaccard ~ temp model
jac.clim.uni.mod <- brm(bf(Jaccard ~ warmq + (1|SiteSubsite), zoi ~ 1, coi ~ 1), init = 0,
                    family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                    data = turn_clim2, iter = 2000, chains = 4, warmup = 500, 
                    file = "models/24_02_jaccard_clim_univ_mod")

summary(jac.clim.uni.mod) # negative significant (same as multivariate)

# predictions
jac.clim.uni.pred <- ggpredict(jac.clim.uni.mod, terms = "warmq", 
                               type = "random", condition = c(SiteSubsite = 1266),
                               allow_new_levels = TRUE)
colnames(jac.clim.uni.pred) = c('warmq', 'fit', 'lwr', 'upr', 'group')






#### 03) FG COMPOSITION MODEL ####


# One model per func group
jac.shb.mod <- brm(bf(Jaccard ~ ShrubDivided + LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                  family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                  data = beta_turnover2, iter = 2000, chains = 4, warmup = 500, 
                  file = "models/24_03_jaccard_shb_mod")

summary(jac.shb.mod) # shrub ns, plot size ns, richness negative, duration ns
summary(jac.shb.mod, prob = 0.975) # shrub ns, plot size ns, richness negative, duration ns

jac.forb.mod <- brm(bf(Jaccard ~ ForbDivided + LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                   family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                   data = beta_turnover2, iter = 2000, chains = 4, warmup = 500, 
                   file = "models/24_03_jaccard_forb_mod")

summary(jac.forb.mod) # forb ns, plot size ns, richness negative, duration ns
summary(jac.forb.mod, prob = 0.975) # forb ns, plot size ns, richness negative, duration ns


jac.gram.mod <- brm(bf(Jaccard ~ GramDivided + LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                    family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                    data = beta_turnover2, iter = 2000, chains = 4, warmup = 500, 
                    file = "models/24_03_jaccard_gram_mod")

summary(jac.gram.mod) # gram ns, plot size ns, richness negative, duration ns
summary(jac.gram.mod, prob = 0.975) # gram ns, plot size ns, richness negative, duration ns



#### 04) PLOT CHANGE MODEL ####
change_clim <- read.csv("data/22clim_slopes.csv")
change_fg <- read.csv("data/22fg_slopes.csv")

turnover.change0 <- left_join(beta_turnover2, change_clim, by = "SiteSubsitePlot")
turnover.change <- left_join(turnover.change0, change_fg, by = "SiteSubsitePlot")


# Split into 3 models for convergence


# PCHG shrub
jac.shbchg.mod <- brm(bf(Jaccard ~ WarmQSlope + PrecSlope + ShrubSlope + 
                        LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                      family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                      data = turnover.change, iter = 2000, chains = 4, warmup = 500, 
                      file = "models/24_04_jaccard_shbchg_mod")
summary(jac.shbchg.mod) # temp positive, prec ns, shrub ns, plot size ns, richness negative, duration ns
summary(jac.shbchg.mod, prob = 0.975) # temp positive, prec ns, shrub ns, plot size ns, richness negative, duration ns

# predictions
jac.shbchg.pred <- ggpredict(jac.shbchg.mod, terms = "ShrubSlope", type = "random", 
                          condition = c(SiteSubsite = 1266), allow_new_levels = TRUE)
colnames(jac.shbchg.pred) = c('ShrubSlope', 'fit', 'lwr', 'upr', 'group')




# PCHG forb
jac.forbchg.mod <- brm(bf(Jaccard ~ WarmQSlope + PrecSlope + ForbSlope + 
                         LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                       family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                       data = turnover.change, iter = 2000, chains = 4, warmup = 500, 
                       file = "models/24_04_jaccard_forbchg_mod")
summary(jac.forbchg.mod) # temp positive, prec ns, forb positive, plot size ns, richness negative, duration ns
summary(jac.forbchg.mod, prob = 0.975) #temp positive, prec ns, forb ns, plot size ns, richness negative, duration ns



# PCHG graminoid
jac.gramchg.mod <- brm(bf(Jaccard ~ WarmQSlope + PrecSlope + GraminoidSlope + 
                         LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                       family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                       data = turnover.change, iter = 2000, chains = 4, warmup = 500, 
                       file = "models/24_04_jaccard_gramchg_mod")
summary(jac.gramchg.mod) # temp positive, prec ns, gram negative, plot size ns, richness negative, duration ns
summary(jac.gramchg.mod, prob = 0.975) # temp positive, prec ns, gram ns, plot size ns, richness negative, duration ns


# univariate model
jac.uni.mod <- brm(bf(Jaccard ~ WarmQSlope + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                   family = "zero_one_inflated_beta", 
                   prior = prior(gamma(0.01, 0.01), class = phi), data = turnover.change, 
                   iter = 3500, chains = 2, warmup = 500, 
                   file = "models/24_04_jaccard_uni_mod")

summary(jac.uni.mod) # positive significant (same as univariate model)


# predictions
jac.uni.pred <- ggpredict(jac.uni.mod, terms = "WarmQSlope", type = "random", 
                          condition = c(SiteSubsite = 1266), allow_new_levels = TRUE)
colnames(jac.uni.pred) = c('WarmQSlope', 'fit', 'lwr', 'upr', 'group')


# plot individually
(bjac.uni.mod.plot <- ggplot() + 
    geom_point(data = turnover.change, aes(x=WarmQSlope, y = Jaccard), size = 2, alpha = 0.1, colour = "#307351") + 
    geom_line(data = jac.uni.pred, aes(x = WarmQSlope, y = fit)) +
    geom_ribbon(data = jac.uni.pred, aes(x = WarmQSlope, ymin = lwr, ymax = upr), fill = "#307351", alpha = 0.2) + 
    bio.theme)




## BRAY-CURTIS MODELS ----
hist(beta_turnover2$BrayCurtis) # a few 0s and a 1

#### 01) GEOGRAPHICAL MODEL ####
bray.geo.mod <- brm(bf(BrayCurtis ~ LAT + Region + LogPlotSize + 
                       MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                     data = beta_turnover2, iter = 2000, chains = 4, warmup = 500, 
                     family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi),
                     file = "models/24_01_bray_geo_mod")

summary(bray.geo.mod) # latitude negative, plot size ns, richness positive, duration positive
conditional_effects(bray.geo.mod) # Eurasia more turnover than NAm-west

summary(bray.geo.mod, prob = 0.975) # latitude negative, plot size ns, richness positive, duration positive, region ns
conditional_effects(bray.geo.mod, prob = 0.975) # Eurasia more turnover than NAm-west


# BC ~ duration plot

# predictions
bray.geo.pred <- ggpredict(bray.geo.mod, terms = "duration", type = "random", 
                           condition = c(SiteSubsite = 1266),
                           allow_new_levels = TRUE)
colnames(bray.geo.pred) = c('duration', 'fit', 'lwr', 'upr', 'group')


# plot relationship
(bray.geo.pred.plot <- ggplot() + 
    geom_point(data = beta_turnover2, aes(x=duration, y = BrayCurtis), 
               size = 2, pch = 21, colour = "black", fill = "#40498E70") + 
    xlab("\nDuration (years)") + ylab("Bray-Curtis\n") +
    geom_line(data = bray.geo.pred, aes(x = duration, y = fit), linewidth = 1, colour = "#40498E70") +
    geom_ribbon(data = bray.geo.pred, aes(x = duration, ymin = lwr, ymax = upr), 
                fill = "#40498E70", alpha = 0.3) + 
    bio.theme +
    theme(axis.title.x = element_text(face="plain", size=9),
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y  = element_text(vjust=0.5, size=9, colour = "black"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines")))

(duration.panel <- plot_grid(jac.geo.pred.plot, bray.geo.pred.plot, 
                            ncol = 2, align = "hv", labels = c("g", "h"), label_size = 14))


#### 02) CLIMATE MODEL ####
bray.clim.mod <- brm(bf(BrayCurtis ~ MOISTURE + warmq + prec + duration + 
                        LogPlotSize + MeanRichness + (1|SiteSubsite), zoi ~ 1, coi ~ 1), init = 0,
                      family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                      data = turn_clim2, iter = 2000, chains = 4, warmup = 500, 
                      file = "models/24_02_bray_clim_mod")

print(summary(bray.clim.mod), digits = 5) # temp ns, prec positive, duration ns, plot size ns, richness positive, moisture ns
conditional_effects(bray.clim.mod)

print(summary(bray.clim.mod, prob = 0.975), digits = 5) # temp ns, prec positive, duration ns, plot size ns, richness positive, moisture ns
conditional_effects(bray.clim.mod, prob = 0.975)


# univariate model
bray.clim.uni.mod <- brm(bf(BrayCurtis ~ warmq + (1|SiteSubsite), zoi ~ 1, coi ~ 1), init = 0,
                     family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                     data = turn_clim2, iter = 3500, chains = 2, warmup = 500, 
                     file = "models/24_02_bray_clim_uni_mod")

summary(bray.clim.uni.mod) # temp positive significant (multivariate was ns)

# predictions
bray.clim.uni.pred <- ggpredict(bray.clim.uni.mod, terms = "warmq", type = "random",
                                condition = c(SiteSubsite = 1266), allow_new_levels = TRUE)
colnames(bray.clim.uni.pred) = c('warmq', 'fit', 'lwr', 'upr', 'group')

# plot turnover ~ temp for both metrics (univariate mods)
(turn.temp.uni.mod.plot <- ggplot() + 
    geom_point(data = turn_clim2, aes(x = warmq, y = BrayCurtis), size = 4, shape = 16, colour = "#40498E70", alpha = 0.7) + 
    geom_point(data = turn_clim2, aes(x = warmq, y = Jaccard), size = 4, shape = 16, colour = "#54C9AD70", alpha = 0.7) + 
    geom_line(data = bray.clim.uni.pred, aes(x = warmq, y = fit), colour = "#40498E70", size = 2) +
    geom_line(data = jac.clim.uni.pred, aes(x = warmq, y = fit), colour = "#54C9AD70", size = 2) +
    geom_ribbon(data = bray.clim.uni.pred, aes(x = warmq, ymin = lwr, ymax = upr), fill = "#40498E70", alpha = 0.2) + 
    geom_ribbon(data = jac.clim.uni.pred, aes(x = warmq, ymin = lwr, ymax = upr), fill = "#54C9AD70", alpha = 0.3) + 
    ylim(0, 1) +
    xlab("\nWarmest Quarter Temperature (°C)") + ylab("Turnover\n") + bio.theme + 
    theme(axis.title.x = element_text(face="plain", size=24),
          axis.text.x  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.title.y = element_text(face="plain", size=24),
          axis.text.y  = element_text(vjust=0.5, size=22, colour = "black")))

ggplot2::ggsave(turn.temp.uni.mod.plot, filename = "scripts/mgarciacriado/figures/Figure_3a.png", 
                width = 20, height = 20, units = "cm")


#### 03) FG COMPOSITION MODEL ####

# One model per func group as a single model doesn't converge
bray.shb.mod <- brm(bf(BrayCurtis ~ ShrubDivided + LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                   family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                   data = beta_turnover2, iter = 2000, chains = 4, warmup = 500, 
                   file = "models/24_03_bray_shb_mod")

summary(bray.shb.mod) # shrub ns, plot size ns, richness positive, duration positive
summary(bray.shb.mod, prob = 0.975) # shrub ns, plot size ns, richness positive, duration positive


bray.forb.mod <- brm(bf(BrayCurtis ~ ForbDivided + LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                    family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                    data = beta_turnover2, iter = 2000, chains = 4, warmup = 500, 
                    file = "models/24_03_bray_forb_mod")

print(summary(bray.forb.mod), digits = 5) # forb positive, plot size ns, richness positive, duration positive
summary(bray.forb.mod, prob = 0.975) # forb ns, plot size ns, richness ns, duration positive


bray.gram.mod <- brm(bf(BrayCurtis ~ GramDivided + LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                    family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                    data = beta_turnover2, iter = 2000, chains = 4, warmup = 500, 
                    file = "models/24_03_bray_gram_mod")

summary(bray.gram.mod) # gram ns, plot size ns, richness positive, duration positive
summary(bray.gram.mod, prob = 0.975) # gram ns, plot size ns, richness positive, duration positive





#### 04) PLOT CHANGE MODEL ####

# Model with climate and FG change over time
bray.shb.chg.mod <- brm(bf(BrayCurtis ~ WarmQSlope + PrecSlope + ShrubSlope + 
                          LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~1),
                        family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                        data = turnover.change, iter = 2000, chains = 4, warmup = 500, 
                        file = "models/24_04_bray_shbchg_mod")
summary(bray.shb.chg.mod) # temp negative, prec ns, shrub ns, plot size ns, richness ns, duration positive
summary(bray.shb.chg.mod, prob = 0.975) # temp ns, prec ns, shrub ns, plot size ns, richness ns, duration positive


# predictions
bray.shbchg.pred <- ggpredict(bray.shb.chg.mod, terms = "ShrubSlope", type = "random", 
                             condition = c(SiteSubsite = 1266), allow_new_levels = TRUE)
colnames(bray.shbchg.pred) = c('ShrubSlope', 'fit', 'lwr', 'upr', 'group')

## Figure with both turnover metrics and shrub change
(turnover.shrub.plot <- ggplot() + 
    geom_point(data = turnover.change, aes(x=ShrubSlope, y = BrayCurtis), size = 4, alpha = 0.7, shape = 16, colour = "#40498E70") + 
    geom_point(data = turnover.change, aes(x=ShrubSlope, y = Jaccard), size = 4, alpha = 0.7, shape = 16, colour = "#54C9AD70") + 
    geom_line(data = bray.shbchg.pred, aes(x = ShrubSlope, y = fit), colour = "#40498E70", size = 2, linetype = "dashed") +
    geom_ribbon(data = bray.shbchg.pred, aes(x = ShrubSlope, ymin = lwr, ymax = upr), fill = "#40498E70", alpha = 0.2) + 
    geom_line(data = jac.shbchg.pred, aes(x = ShrubSlope, y = fit), colour = "#54C9AD70", size = 2, linetype = "dashed") +
    geom_ribbon(data = jac.shbchg.pred, aes(x = ShrubSlope, ymin = lwr, ymax = upr), fill = "#54C9AD70", alpha = 0.3) + 
    ylim(0, 1) +
    xlab("\nShrub cover change (% per year)") + ylab("Turnover\n") + bio.theme +
    theme(axis.title.x = element_text(face="plain", size=24),
          axis.text.x  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.title.y = element_text(face="plain", size=24),
          axis.text.y  = element_text(vjust=0.5, size=22, colour = "black")))

ggplot2::ggsave(turnover.shrub.plot, filename = "scripts/mgarciacriado/figures/Figure_3c.png", 
                width = 20, height = 20, units = "cm")



bray.forb.chg.mod <- brm(bf(BrayCurtis ~ WarmQSlope + PrecSlope + ForbSlope + 
                           LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1), init = 0,
                         family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                         data = turnover.change, iter = 2000, chains = 4, warmup = 500, 
                         control = list(adapt_delta = 0.9),
                         file = "models/24_04_bray_forbchg_mod")
print(summary(bray.forb.chg.mod), digits = 5) # temp ns, prec ns, forb ns, plot size ns, richness positive, duration positive
summary(bray.forb.chg.mod, prob = 0.975) # temp ns, prec ns, forb ns, plot size ns, richness ns, duration positive


bray.gram.chg.mod <- brm(bf(BrayCurtis ~ WarmQSlope + PrecSlope + GraminoidSlope + 
                           LogPlotSize + MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                         family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi), 
                         data = turnover.change, iter = 2000, chains = 4, warmup = 500, 
                         file = "models/24_04_bray_gramchg_mod")
summary(bray.gram.chg.mod) # temp negative, prec ns, gram ns, plot size ns, richness positive, duration positive
summary(bray.gram.chg.mod, prob = 0.975) # temp ns, prec ns, gram ns, plot size ns, richness positive, duration positive


# univariate model 
bray.uni.mod <- brm(bf(BrayCurtis ~ WarmQSlope + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                   family = "zero_one_inflated_beta", 
                   prior = prior(gamma(0.01, 0.01), class = phi), data = turnover.change, 
                   iter = 3500, chains = 2, warmup = 500, 
                   file = "models/24_04_bray_uni_mod")

summary(bray.uni.mod) # negative ns

# bray predictions
bray.uni.pred <- ggpredict(bray.uni.mod, terms = "WarmQSlope", 
                           type = "random", condition = c(SiteSubsite = 1266), 
                           allow_new_levels = TRUE)
colnames(bray.uni.pred) = c('WarmQSlope', 'fit', 'lwr', 'upr', 'group')


# plot univariate models (both Jac/BC) with points
(bjac.uni.mod.points.plot <- ggplot() + 
    geom_point(data = turnover.change, aes(x=WarmQSlope, y = BrayCurtis), size = 4, alpha = 0.7, shape = 16, colour = "#40498E70") + 
    geom_point(data = turnover.change, aes(x=WarmQSlope, y = Jaccard), size = 4, alpha = 0.7, shape = 16, colour = "#54C9AD70") + 
    geom_line(data = bray.uni.pred, aes(x = WarmQSlope, y = fit), colour = "#40498E70", size = 2, linetype = "dashed") +
    geom_ribbon(data = bray.uni.pred, aes(x = WarmQSlope, ymin = lwr, ymax = upr), fill = "#40498E70", alpha = 0.2) + 
    geom_line(data = jac.uni.pred, aes(x = WarmQSlope, y = fit), colour = "#54C9AD70", size = 2) +
    geom_ribbon(data = jac.uni.pred, aes(x = WarmQSlope, ymin = lwr, ymax = upr), fill = "#54C9AD70", alpha = 0.3) + 
    ylim(0, 1) +
    xlab("\nTemperature change (°C per year)") + ylab("Turnover\n") + bio.theme +
    theme(axis.title.x = element_text(face="plain", size=24),
          axis.text.x  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.title.y = element_text(face="plain", size=24),
          axis.text.y  = element_text(vjust=0.5, size=22, colour = "black")))

ggplot2::ggsave(bjac.uni.mod.points.plot, filename = "scripts/mgarciacriado/figures/Figure_3b.png", 
                width = 20, height = 20, units = "cm")






## MEANS & DISTRIBUTIONS ----
nochange.jac <- beta_turnover %>% filter(Jaccard == 0) #526/1266 = 41.5% plots haven't changed at all
nochange.bc <- beta_turnover %>% filter(BrayCurtis == 0) #6/1266 = 0.4% plots haven't changed at all

# Overlay histograms and Density plots of turnover and richness change at the PLOT level
(jaccard.dens <- ggplot(beta_turnover, aes(x=Jaccard)) + 
    geom_histogram(aes(y = ..count..), breaks = seq(0, 1, by = 0.05), fill="lightgray", alpha = 0.7) +
    geom_density(aes(y = ..density..*(1259*0.05)), fill="#307351", colour="#307351", alpha = 0.5, size = 1.2) + 
    geom_vline(aes(xintercept=mean(Jaccard)), color="#307351", linetype="dotted", size=1) + 
    annotate("text", x = 0.6, y = 500, label = "526 (41.5%) plots\ndid not change at all", size = 7.5) +
    xlab("\nJaccard Turnover\n(presence/absence)") + ylab("Count\n") + bio.theme)


(braycurtis.dens <- ggplot(beta_turnover, aes(x=BrayCurtis)) + 
    geom_histogram(aes(y = ..count..), breaks = seq(0, 1, by = 0.05), fill="lightgray", alpha = 0.7) +
    geom_density(aes(y = ..density..*(1259*0.05)), fill="#70163C", colour="#70163C", alpha = 0.5, size = 1.2) + 
    geom_vline(aes(xintercept=mean(BrayCurtis)), color="#70163C", linetype="dotted", size=1) + 
    annotate("text", x = 0.63, y = 168, label = "6 (0.4%) plots\ndid not change at all", size = 7.5) +
    xlab("\nBray-Curtis Turnover\n(presence/absence + abundance)") + ylab("Count\n") + bio.theme)



# What is the mean Jaccard turnover across plots? Intercept-only model
turn.intonly.mod <- brm(Jaccard ~ 1, data = beta_turnover, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/24turn_intonly_mod")
summary(turn.intonly.mod) # mean turnover is 0.22 

# What is the mean BC turnover across plots? Intercept-only model
bray.intonly.mod <- brm(BrayCurtis ~ 1, data = beta_turnover, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/24bray_intonly_mod")
summary(bray.intonly.mod) # mean turnover is 0.36 




## RELATIONSHIPS ---- 





## CHECKS ----

# Turnover value of 1 - the one species has been completely replaced by another
alex.check <- filter(itex.dec22, SiteSubsitePlot == "ALEXFIORD:DOMEDOLOMITE:Dom.c.4d")

# Turnover value of 0 - all species remain the same (just an example)
zack.check <- filter(itex.dec22, SiteSubsitePlot == "ZACKENBERG:VACCINIUM ULIGONOSUM HEATH:V5")

# Calculate richness per year
rich.timexx <- itex.dec22 %>% filter(RelCover > 0) %>% 
  group_by(SiteSubsitePlotYear) %>% 
  mutate(Richness = length(unique(SPECIES_NAME))) %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  ungroup() %>% group_by(SiteSubsitePlot) %>% 
  mutate(Duration = max(YEAR) - min(YEAR)) %>% 
  filter(Duration > 1) %>% ungroup()

# This file might come in handy
write.csv(rich.timexx, "data/22rich_time.csv")






## NULL MODEL ----
null.start <- bio_turnover %>% group_by(SiteSubsitePlot) %>% filter(YEAR == min(YEAR))

# only for the last year we add and remove species
null.end <- bio_turnover %>% group_by(SiteSubsitePlot) %>% filter(YEAR == max(YEAR)) %>% 
  mutate(Prop20 = round(0.2*unique(length(SPECIES_NAME)), digits = 0)) %>%
  # remove 20% of species in each plot, prop refers to the values that are kept
  slice_sample(prop = 0.8) %>% ungroup()

unique(null.end$Prop20) #0, 1, 2, 3, 4, 5
  
# Add in the extra colonising species
# Okay this is the ugliest code ever but am not inspired for a more elegant solution. 
# At least it works right?

null.end0 <- null.end %>% filter(Prop20 == 0)

null.end1 <- null.end %>% filter(Prop20 == 1) %>% group_by(SiteSubsitePlot, YEAR) %>% 
  group_modify(~ add_row(., .before=0)) %>% replace_na(list(SPECIES_NAME = "Species1"))
  
null.end2 <- null.end %>% filter(Prop20 == 2) %>% group_by(SiteSubsitePlot, YEAR) %>% 
  group_modify(~ add_row(., SPECIES_NAME = c("Species1", "Species2"), .before=0)) 

null.end3 <- null.end %>% filter(Prop20 == 3) %>% group_by(SiteSubsitePlot, YEAR) %>% 
  group_modify(~ add_row(., SPECIES_NAME = c("Species1", "Species2", "Species3"), .before=0)) 

null.end4 <- null.end %>% filter(Prop20 == 4) %>% group_by(SiteSubsitePlot, YEAR) %>% 
  group_modify(~ add_row(., SPECIES_NAME = c("Species1", "Species2", "Species3", "Species4"), .before=0)) 

null.end5 <- null.end %>% filter(Prop20 == 5) %>% group_by(SiteSubsitePlot, YEAR) %>% 
  group_modify(~ add_row(., SPECIES_NAME = paste0("Species", 1:5), .before=0)) 

null.end.full <- rbind(null.end0, null.end1, null.end2, null.end3, null.end4, null.end5)

# Fill in relative cover
null.end.cover <- null.end.full %>% group_by(SiteSubsitePlot) %>% 
  mutate(CoverSum = sum(RelCover, na.rm=TRUE), RemainingCover = 100 - CoverSum) %>%
  mutate(NumberNA = sum(is.na(RelCover))) %>% mutate(NewCover = RemainingCover/NumberNA) %>%
  mutate(CompleteCover = case_when(is.na(RelCover) ~ NewCover, TRUE ~ RelCover)) %>% ungroup()

null.end.final <- null.end.cover %>% dplyr::select(SiteSubsitePlot, YEAR, SPECIES_NAME, CompleteCover) %>% 
  rename(RelCover = CompleteCover)
        
# Combine start and end timepoints
null.full <- rbind(null.start, null.end.final)
null.full2 <- null.full %>% group_by(SiteSubsitePlot) %>% 
  filter(length(unique(YEAR)) > 1) %>% filter(RelCover >= 0) %>%
  ungroup()



# Create empty dataframe
beta_null <- data.frame(matrix(ncol = 10, nrow = length(unique(null.full2$SiteSubsitePlot)))) 
names(beta_null) <- c("SiteSubsitePlot", "duration", "richness", "richness_change", 
                         "Jbeta", "Jtu", "Jne", "Bbal", "Bgra", "Bbray") 
i = 1



# Loop to calculate turnover with betapart
for (i in 1:length(unique(null.full2$SiteSubsitePlot))) {
  SiteSubsitePlotName <- as.character(unique(null.full2$SiteSubsitePlot)[i])
  sub_bio_abundance <- filter(null.full2, SiteSubsitePlot == SiteSubsitePlotName)
  YearMin <- min(sub_bio_abundance$YEAR)
  YearMax <- max(sub_bio_abundance$YEAR)
  duration <- YearMax - YearMin
  # filters the dataframe for just the first and last observations per plot
  sub_bio_abundance_min <- filter(sub_bio_abundance, YEAR == YearMin)
  sub_bio_abundance_max <- filter(sub_bio_abundance, YEAR == YearMax)
  sub_bio_abundance <- rbind(sub_bio_abundance_min, sub_bio_abundance_max)
  richness <- length(unique(sub_bio_abundance$SPECIES_NAME))
  # averages any species that have multiple records per plot and time point 
  sub_bio_abundance <- sub_bio_abundance %>% group_by(SiteSubsitePlot, SPECIES_NAME, YEAR) %>% summarise(MeanRelCover = mean(RelCover)) %>% ungroup()
  # reshape to wide form
  sub_bio_abundance_wider <- tidyr::pivot_wider(sub_bio_abundance, names_from = SPECIES_NAME, 
                                                values_from = MeanRelCover, 
                                                values_fill = list(MeanRelCover = 0))
  # removes columns for beta.pair() function
  sub_bio_abundance_matrix <- dplyr::select(sub_bio_abundance_wider, -SiteSubsitePlot, -YEAR) 
  # creates presence-absence matrix
  sub_bio_presence_matrix <- with(sub_bio_abundance_matrix, ifelse(sub_bio_abundance_matrix > 0,1,0))
  # calculates Jaccard overall, turnover and nestedness
  J_components <- beta.pair(sub_bio_presence_matrix, index.family='jaccard')
  B_components <- beta.pair.abund(sub_bio_abundance_matrix, index.family='bray')
  # saves biodiversity metrics
  richness_change <- rowSums(sub_bio_presence_matrix)[2] - rowSums(sub_bio_presence_matrix)[1]
  Jbeta <- J_components$beta.jac
  Jtu <- J_components$beta.jtu
  Jne <- J_components$beta.jne
  Bbal <- B_components$beta.bray.bal
  Bgra <- B_components$beta.bray.gra
  Bbray <- B_components$beta.bray
  beta_null[i,] <- c(SiteSubsitePlotName, duration, richness, richness_change, Jbeta, Jtu, Jne, Bbal, Bgra, Bbray)
  i = i+1
}



# Convert in adequate format 
beta_null_final <- beta_null %>% 
  mutate(duration = as.numeric(duration), richness = as.numeric(richness), richness_change = as.numeric(richness_change),
         richness_change_abs = abs(as.numeric(richness_change)), Jtu = as.numeric(Jtu), Bbray = as.numeric(Bbray)) %>%
  dplyr::select(., -c(Jbeta, Jne, Bbal, Bgra)) %>% rename(Jaccard = Jtu, BrayCurtis = Bbray)

# All makes sense
glimpse(beta_null_final)

itex.dec22$SurveyedArea <- as.numeric(itex.dec22$SurveyedArea)



# What is the mean turnover across plots? Intercept-only model
jac.null.mod <- brm(Jaccard ~ 1, data = beta_null_final, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/24jac_null_mod")
summary(jac.null.mod) # mean turnover is 0.48

# What is the mean turnover across plots? Intercept-only model
bray.null.mod <- brm(BrayCurtis ~ 1, data = beta_null_final, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/24bray_null_mod")
summary(bray.null.mod) # mean turnover is 0.52




# Extract covariates for the models
add.covs <- itex.dec22 %>% distinct(SiteSubsitePlot, .keep_all = TRUE) %>%
  dplyr::select(SiteSubsite, SiteSubsitePlot, MOISTURE, LAT, LONG, SurveyedArea, 
                Region, PlotShrubMean, PlotGraminoidMean, PlotForbMean, MeanRichness) %>%
  mutate(LogPlotSize = log(SurveyedArea))

# Merge the 2 dataframes
beta_null_turnover <- left_join(beta_null_final, add.covs, by = "SiteSubsitePlot")


# Fit geo model - Jaccard
jac.null.geo.mod <- brm(bf(Jaccard ~ LAT + Region + LogPlotSize + 
                        MeanRichness + duration + (1|SiteSubsite), zoi ~ 1, coi ~ 1),
                   data = beta_null_turnover, iter = 2000, chains = 4, warmup = 500, 
                   family = "zero_one_inflated_beta", prior = prior(gamma(0.01, 0.01), class = phi),
                   file = "models/24_01_null_jaccard_geo_mod")


summary(jac.null.geo.mod) # latitude ns, plot size positive, richness ns, duration ns, greater turnover in GreenIceland
conditional_effects(jac.null.geo.mod)

# there are 1s but not 0s - scale so beta fits
beta_null_bray <- beta_null_turnover %>% mutate(BrayScaled = BrayCurtis - 0.0001)

# Fit geo model - BC
bray.null.geo.mod <- brm(BrayScaled ~ LAT + Region + LogPlotSize + 
                         MeanRichness + duration + (1|SiteSubsite),
                    data = beta_null_bray, iter = 2000, chains = 4, warmup = 500, 
                    family = "beta", file = "models/24_01_null_bray_geo_mod")

summary(bray.null.geo.mod) # latitude ns, plot size ns, richness negative, duration ns, greater turnover in Eurasia than NAm-west
conditional_effects(bray.null.geo.mod)


