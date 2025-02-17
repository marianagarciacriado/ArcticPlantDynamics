## Plant diversity dynamics over space and time in a warming Arctic
## Mariana Garcia Criado (mariana.garcia.criado@gmail.com)
## Script 6. Homogeneization/heterogeneization 

# Thanks to Tora Finderup Nielsen for initial cluster code!


## LIBRARIES ----
library(vegan)
library(ape) 
library(viridis)
library(cowplot)
library(tidyverse)
library(broom)
library(brms)
library(ggrepel)

#library(githubinstall)
#githubinstall("ape")

#library(devtools)
#devtools::install_github("emmanuelparadis/ape")

#library(remotes)
#remotes::install_github("emmanuelparadis/ape")




## DATA ----
load("data/itex_dec22.RData")

# Create dataset where composition was measured at at least two time points
abundance <- itex.dec22 %>%
  dplyr::select(SiteSubsitePlotYear, SiteSubsitePlot, SiteSubsite, SITE, YEAR, SPECIES_NAME, RelCover, 
                LAT, FuncGroup, DominatingFG, PlotYearDominatingFG) %>%
  group_by(SiteSubsitePlot) %>%
  filter(length(unique(YEAR)) > 1, RelCover >= 0) %>%
  mutate(Duration = max(YEAR) - min(YEAR)) %>% filter(Duration > 4) %>%
  ungroup()


## Calculate functional group change

# Abundance of func groups for start year
fg.start <- abundance %>% group_by(SiteSubsitePlot) %>% 
  filter(YEAR == min(YEAR)) %>% ungroup() %>%
  group_by(SiteSubsitePlotYear, FuncGroup) %>% mutate(StartCoverFG = sum(RelCover)) %>%
  ungroup() %>% distinct(SiteSubsitePlotYear, FuncGroup, .keep_all = TRUE) %>% 
  dplyr::select(SiteSubsitePlot, YEAR, FuncGroup, StartCoverFG, PlotYearDominatingFG) %>% 
  rename(StartYear = YEAR) %>% rename(StartDominatingFG = PlotYearDominatingFG)

# Abundance of func groups for end year
fg.end <- abundance %>% group_by(SiteSubsitePlot) %>% 
  filter(YEAR == max(YEAR)) %>% ungroup() %>% 
  group_by(SiteSubsitePlotYear, FuncGroup) %>% mutate(EndCoverFG = sum(RelCover)) %>%
  ungroup() %>% distinct(SiteSubsitePlotYear, FuncGroup, .keep_all = TRUE) %>% 
  dplyr::select(SiteSubsitePlot, YEAR, FuncGroup, EndCoverFG, PlotYearDominatingFG) %>% 
  rename(EndYear = YEAR) %>% rename(EndDominatingFG = PlotYearDominatingFG)

# it's okay that they have different number of rows because there are some func groups that are at the start and not at the end
# and viceversa. they will be filled in below with NA/0.




#### DATASET 1 (ABSOLUTE FG CHANGE) ####

# Merge the two dataframes - full join includes those FG that are only present at start or end.
fg.meta0 <- full_join(fg.end, fg.start, by = c("SiteSubsitePlot", "FuncGroup"))

# Fill in NA cover as 0s (they weren't present at the start or end)
fg.meta <- fg.meta0 %>% mutate(StartCoverFG = ifelse(is.na(StartCoverFG), 0, StartCoverFG)) %>% 
  mutate(EndCoverFG = ifelse(is.na(EndCoverFG), 0, EndCoverFG))
  
# Extract original missing start/end year and dominating group info
start.info <- fg.start %>% distinct(SiteSubsitePlot, .keep_all = TRUE) %>% 
  dplyr::select(SiteSubsitePlot, StartYear, StartDominatingFG)

end.info <- fg.end %>% distinct(SiteSubsitePlot, .keep_all = TRUE) %>% 
  dplyr::select(SiteSubsitePlot, EndYear, EndDominatingFG)

start.end.info <- left_join(start.info, end.info, by = "SiteSubsitePlot")

# Fill this back into the main db 
fg.meta2 <- left_join(fg.meta, start.end.info, by = "SiteSubsitePlot")

# Remove extra columns and calculate absolute change
fg.meta3 <- fg.meta2 %>% dplyr::select(., -c(StartDominatingFG.x, StartYear.x, EndYear.x, EndDominatingFG.x)) %>%
  rename(StartDominatingFG = StartDominatingFG.y, StartYear = StartYear.y, 
         EndYear = EndYear.y, EndDominatingFG = EndDominatingFG.y) %>%
  mutate(FuncGroupChange = EndCoverFG - StartCoverFG) %>%
  mutate(FGtrend = case_when(FuncGroupChange > 0 ~ "Increase",
                             FuncGroupChange == 0 ~ "Stable",
                             FuncGroupChange < 0 ~ "Decrease"))


# Extract change per functional group per site
fg.meta.wide <- fg.meta3 %>% 
  tidyr::pivot_wider(names_from = FuncGroup, values_from = FuncGroupChange) %>%
  dplyr::select(SiteSubsitePlot, Shrub, Graminoid, Forb) %>%
  group_by(SiteSubsitePlot) %>% 
  summarise_each(funs(first(.[!is.na(.)]))) %>%
  rename(ShrubChange = Shrub, GramChange = Graminoid, ForbChange = Forb) %>%
  mutate(ShrubTrend = case_when(ShrubChange > 0 ~ "Increase",
                                ShrubChange == 0 ~ "Stable",
                                ShrubChange < 0 ~ "Decrease")) %>%
  mutate(GramTrend = case_when(GramChange > 0 ~ "Increase",
                                GramChange == 0 ~ "Stable",
                                GramChange < 0 ~ "Decrease")) %>%
  mutate(ForbTrend = case_when(ForbChange > 0 ~ "Increase",
                                ForbChange == 0 ~ "Stable",
                                ForbChange < 0 ~ "Decrease"))

# those with start cover = 0 and end cover = 0 remain as NA
write.csv(fg.meta.wide, "data/22fg_change.csv")

#fg.meta.wide <- read.csv("scripts/mgarciacriado/data/22fg_change.csv")




#### DATASET 2 (FG CHANGE SLOPES) ####

# Abundance of func groups for all years
fg.values <- itex.dec22 %>% group_by(SiteSubsitePlot) %>%
  filter(length(unique(YEAR)) > 1) %>% 
  mutate(Duration = max(YEAR) - min(YEAR)) %>% filter(Duration > 4) %>% ungroup() %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  dplyr::select(SiteSubsitePlot, SiteSubsitePlotYear, YEAR, ShrubCover, GraminoidCover, ForbCover)

# Add values of NA for those functional groups that were never present
fg.values2 <- fg.values %>% group_by(SiteSubsitePlot) %>% 
  mutate(ShrubTotal = sum(ShrubCover), ForbTotal = sum(ForbCover), GramTotal = sum(GraminoidCover)) %>%
  mutate(ShrubCoverNew = ifelse(ShrubTotal == 0, NA, ShrubCover), 
         ForbCoverNew = ifelse(ForbTotal == 0, NA, ForbCover),
         GramCoverNew = ifelse(GramTotal == 0, NA, GraminoidCover)) %>% ungroup() %>%
  dplyr::select(SiteSubsitePlot, SiteSubsitePlotYear, YEAR, ShrubCoverNew, ForbCoverNew, GramCoverNew)
  

# Convert to long format (with NA)
fg.values.long <- fg.values2 %>% 
  pivot_longer(cols = c(ShrubCoverNew, GramCoverNew, ForbCoverNew), 
               names_to = "FuncGroup", values_to = "CoverFG")

# Long format without NA
fg.values.long.nona <- fg.values.long %>% filter(!is.na(CoverFG))


# Calculate slopes (annual change in functional group) for each plot 
fg.slopes <- fg.values.long.nona %>%
  nest_by(SiteSubsitePlot, FuncGroup) %>%  # replace group_by() with nest_by() to convert your model data to a vector of lists
  mutate(mod = list(lm(CoverFG ~ YEAR, data = data))) %>% # change do() to mutate(), then add list() before your model make sure to change data = .  to data = data
  summarise(tidy(mod)) %>%
  filter(term == "YEAR") %>%
  dplyr::select(SiteSubsitePlot, FuncGroup, estimate) %>% 
  rename(FGSlope = estimate) %>% unnest()

# Convert to wide format, keep those which weren't there at the start and at the end as NA
fg.slopes.wide <- fg.slopes %>% tidyr::pivot_wider(names_from = FuncGroup, values_from = FGSlope) %>% 
  rename(ForbSlope = ForbCoverNew, GraminoidSlope = GramCoverNew, ShrubSlope = ShrubCoverNew)

write.csv(fg.slopes.wide, "data/22fg_slopes.csv")



#### FG CHANGE OVERALL
#fg.slopes.wide <- read.csv("data/22fg_slopes.csv")

subsite.info <- abundance %>% distinct(SiteSubsitePlot, .keep_all = TRUE) %>% 
  dplyr::select(SiteSubsitePlot, SiteSubsite)

fg.slopes.full <- left_join(fg.slopes.wide, subsite.info, by = "SiteSubsitePlot")

# shrub change
shrub.change.mod <- brm(ShrubSlope ~ 1 + (1|SiteSubsite), data = fg.slopes.full, 
                        warmup = 500, iter = 2000, chains = 4, 
                        file = "models/24_shrub_change_mod")
summary(shrub.change.mod) # positive ns


# forb change
forb.change.mod <- brm(ForbSlope ~ 1 + (1|SiteSubsite), data = fg.slopes.full, 
                        warmup = 500, iter = 2000, chains = 4, 
                        file = "models/24_forb_change_mod")
summary(forb.change.mod) # positive ns


# graminoid change
gram.change.mod <- brm(GraminoidSlope ~ 1 + (1|SiteSubsite), data = fg.slopes.full, 
                       warmup = 500, iter = 2000, chains = 4, 
                       file = "models/24_gram_change_mod")
summary(gram.change.mod) # negative ns





#### DATASET 3 (COLUMN PER SPECIES AT SUBSITE LEVEL) ####

# To run the PCOA analyses at the subsite level 
# (average abundance of each species across all plots in a subsite)

# The issue here is that there are different starting and ending times for the different plots in the subsite
# which means that if we extract the min and the max they might be from different plots and not all species will be considered
# A solution is to take the start and end points for all plots in a subsite, and then do the average of all those.
abundance.subs <- abundance %>% group_by(SiteSubsitePlot) %>% 
  mutate(MinYear = min(YEAR)) %>% 
  mutate(MaxYear = max(YEAR)) %>% filter(MinYear != MaxYear) %>% 
  filter(YEAR == min(YEAR) | YEAR == max(YEAR)) %>% ungroup() %>% 
  mutate(Timepoint = case_when(YEAR == MinYear ~ "Begin", 
                               YEAR == MaxYear ~ "End")) %>%
  group_by(SiteSubsite, Timepoint, SPECIES_NAME) %>% 
  mutate(MeanSubsiteCover = mean(RelCover)) %>% ungroup() %>%
# now I need one value only of each species at the subsite level  
  group_by(SiteSubsite, Timepoint) %>% distinct(SPECIES_NAME, .keep_all = TRUE) %>%
  ungroup() %>% dplyr::select(SiteSubsite, Timepoint, SPECIES_NAME, MeanSubsiteCover) 
  

#kangertest <- abundance.subs %>% filter(SiteSubsite == "KANGER:DOPEY") # all means make sense!


# Convert to wide format, replace NAs by 0 - two rows per site (start/end)
abundance.subs.wide <- abundance.subs %>% 
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = SPECIES_NAME, values_from = MeanSubsiteCover) %>%
  mutate(SiteSubsiteTime = paste(SiteSubsite, Timepoint, sep = ":")) %>% 
  group_by(SiteSubsiteTime) %>% summarise_each(funs(first(.[!is.na(.)]))) %>%
  replace(is.na(.), 0) %>% dplyr::select(-c(row, Timepoint, SiteSubsite)) %>% 
  arrange(SiteSubsiteTime) %>%
  tibble::column_to_rownames("SiteSubsiteTime")

# Cross-check that all subsites have two rows: OK!
glimpse(abundance.subs.wide) # all looks good
which(is.na(abundance.subs.wide)) # no NAs

# Save file to start straight down below - 180 obs so 90 subsites
write.csv(abundance.subs.wide, "data/22cover_subs_wide.csv")




## FG STATS (SUMMARY TABLE) ----

# Load dataset 1
fg.meta.wide <- read.csv("scripts/mgarciacriado/data/22fg_change.csv")

# How many plots have increased/decreased per functional group?
fg.stats <- fg.meta.wide %>% group_by(ShrubTrend) %>% summarise(total = n()) %>% 
  mutate(percent = (total/sum(total)*100))

# How many plots have increased/decreased per functional group?
fg.stats2 <- fg.meta.wide %>% group_by(GramTrend) %>% summarise(total = n()) %>% 
  mutate(percent = (total/sum(total)*100))

# How many plots have increased/decreased per functional group?
fg.stats3 <- fg.meta.wide %>% group_by(ForbTrend) %>% summarise(total = n()) %>% 
  mutate(percent = (total/sum(total)*100))





## SUBSITE ANALYSES ----

# Use dataset 4 created above
abundance.subs.wide <- read.csv("data/22cover_subs_wide.csv")


# Make disimilarity-matrix for several methods
dist_subs_jacbi <- vegdist(abundance.subs.wide, method = "jaccard", binary = T)
dist_subs_S <- vegdist(abundance.subs.wide, method = "bray", binary = T)
dist_subs_BC <- vegdist(abundance.subs.wide, method = "bray")
dist_subs_aGw <- vegdist(abundance.subs.wide, method = "altGower")
dist_subs_man <- vegdist(abundance.subs.wide, method = "manhattan")
dist_subs_euc <- vegdist(abundance.subs.wide, method = "euclidian")


# Multivariate analyses as a function of time groups

# Make time vector to be grouped by: each plot is old/new time, total of 180 obs so we repeat the vector 90 times
time_vector <- c("Begin", "End")
Group_Time_Subs <- rep(time_vector, 90)

# Calculate mean distance to centroid for groups: homogeneity of variance
beta_subs_jacbi <- betadisper(dist_subs_jacbi, group = Group_Time_Subs)
beta_subs_jacbi
beta_subs_S <- betadisper(dist_subs_S, group = Group_Time_Subs)
beta_subs_S
beta_subs_BC <- betadisper(dist_subs_BC, group = Group_Time_Subs)
beta_subs_BC
beta_subs_aGw <- betadisper(dist_subs_aGw, group = Group_Time_Subs)
beta_subs_aGw
beta_subs_man <- betadisper(dist_subs_man, group = Group_Time_Subs)
beta_subs_man
beta_subs_euc <- betadisper(dist_subs_euc, group = Group_Time_Subs)
beta_subs_euc

## ANOVA to test if distances to centroid are different between groups 
# ("Null hypothesis of no difference in dispersion between groups")
anova(beta_subs_jacbi) #p = 0.6
anova(beta_subs_S) #p = 0.6
anova(beta_subs_BC) #p = 0.8
anova(beta_subs_aGw) #p = 0.95
anova(beta_subs_man) #p = 0.89
anova(beta_subs_euc) #p = 0.62

# No major differences here though between timepoints


## FIGURE S10

# Extract PCOA values for axes
labs.jac <- paste0("Dimension ", 1:4, " (", round(100*beta_subs_jacbi$eig / sum(beta_subs_jacbi$eig), 2), "%)")
labs.sor <- paste0("Dimension ", 1:4, " (", round(100*beta_subs_S$eig / sum(beta_subs_S$eig), 2), "%)")
labs.bray <- paste0("Dimension ", 1:4, " (", round(100*beta_subs_BC$eig / sum(beta_subs_BC$eig), 2), "%)")
labs.gow <- paste0("Dimension ", 1:4, " (", round(100*beta_subs_aGw$eig / sum(beta_subs_aGw$eig), 2), "%)")
labs.man <- paste0("Dimension ", 1:4, " (", round(100*beta_subs_man$eig / sum(beta_subs_man$eig), 2), "%)")
labs.euc <- paste0("Dimension ", 1:4, " (", round(100*beta_subs_euc$eig / sum(beta_subs_euc$eig), 2), "%)")


# PCOA plots (base R, for ggplots there's conversion issues)
# save as Export-png size 2000x1000

par(mfrow=c(2,6))

plot(beta_subs_jacbi, cex = 1.5, label=F, hull = TRUE, pch = c(17, 16), col = c("#F4B043", "#219EAE"), 
     seg.col = "lightgrey", seg.lty = "dotted", seg.lwd = 1, lwd = 3, lty = "solid", 
     main = "Jaccard", cex.lab=1.5, xlab=labs.jac[1], ylab=labs.jac[2], cex.main = 1.5, cex.axis = 1.45)
plot(beta_subs_S, cex = 1.5, label=F, hull = TRUE, pch = c(17, 16), col = c("#F4B043", "#219EAE"), 
     seg.col = "lightgrey", seg.lty = "dotted", seg.lwd = 1, lwd = 3, lty = "solid", 
     main = "Sørensen", cex.lab=1.5, xlab=labs.sor[1], ylab=labs.sor[2], cex.main = 1.5, cex.axis = 1.45)
plot(beta_subs_BC, cex = 1.5, label=F, hull = TRUE, pch = c(17, 16), col = c("#F4B043", "#219EAE"), 
     seg.col = "lightgrey", seg.lty = "dotted", seg.lwd = 1, lwd = 3, lty = "solid", 
     main = "Bray-Curtis", cex.lab=1.5, xlab=labs.bray[1], ylab=labs.bray[2], cex.main = 1.5, cex.axis = 1.45)
plot(beta_subs_aGw, cex = 1.5, label=F, hull = TRUE, pch = c(17, 16), col = c("#F4B043", "#219EAE"), 
     seg.col = "lightgrey", seg.lty = "dotted", seg.lwd = 1, lwd = 3, lty = "solid", 
     main = "Modified Gower", cex.lab=1.5, xlab=labs.gow[1], ylab=labs.gow[2], cex.main = 1.5, cex.axis = 1.45)
plot(beta_subs_man, cex = 1.5, label=F, hull = TRUE, pch = c(17, 16), col = c("#F4B043", "#219EAE"), 
     seg.col = "lightgrey", seg.lty = "dotted", seg.lwd = 1, lwd = 3, lty = "solid", 
     main = "Manhattan", cex.lab=1.5, xlab=labs.man[1], ylab=labs.man[2], cex.main = 1.5, cex.axis = 1.45)
plot(beta_subs_euc, cex = 1.5, label=F, hull = TRUE, pch = c(17, 16), col = c("#F4B043", "#219EAE"), 
     seg.col = "lightgrey", seg.lty = "dotted", seg.lwd = 1, lwd = 3, lty = "solid", 
     main = "Euclidian", cex.lab=1.5, xlab=labs.euc[1], ylab=labs.euc[2], cex.main = 1.5, cex.axis = 1.45)


# boxplots
boxplot(beta_subs_jacbi, las=2, notch = TRUE, col = c("#F4B043", "#219EAE"), 
        cex.lab=1.5, cex.axis = 1.45, xlab = "")

boxplot(beta_subs_S, las=2, notch = TRUE, col = c("#F4B043", "#219EAE"), 
        cex.lab=1.5, cex.axis = 1.45, xlab = "")

boxplot(beta_subs_BC, las=2, notch = TRUE, col = c("#F4B043", "#219EAE"), 
        cex.lab=1.5, cex.axis = 1.45, xlab = "")

boxplot(beta_subs_aGw, las=2, notch = TRUE, col = c("#F4B043", "#219EAE"), 
        cex.lab=1.5, cex.axis = 1.45, xlab = "")

boxplot(beta_subs_man, las=2, notch = TRUE, col = c("#F4B043", "#219EAE"), 
        cex.lab=1.5, cex.axis = 1.45, xlab = "")

boxplot(beta_subs_euc,las=2, notch = TRUE, col = c("#F4B043", "#219EAE"), 
        cex.lab=1.5, cex.axis = 1.45, xlab = "")

dev.off()



## Adonis

# PERMANOVA test: are the centroids significantly different to each other?
# R2 shows the percentage of variance explained by the groups

adon.subs.jac <- adonis2(dist_subs_jacbi ~ Group_Time_Subs, perm = 999)
adon.subs.jac # F = 0.14, R2 = 0.0008

adon.subs.sor <- adonis2(dist_subs_S ~ Group_Time_Subs, perm = 999) 
adon.subs.sor # F = 0.0567, R2 = 0.00032

adon.subs.bc <- adonis2(dist_subs_BC ~ Group_Time_Subs, perm = 999) 
adon.subs.bc # F = 0.1379, R2 = 0.00077

adon.subs.gow <- adonis2(dist_subs_aGw ~ Group_Time_Subs, perm = 999)
adon.subs.gow # F = 0.753, R2 = 0.00365

adon.subs.man <- adonis2(dist_subs_man ~ Group_Time_Subs, perm = 999)
adon.subs.man # F = 0.13, R2 = 0.00078

adon.subs.euc <- adonis2(dist_subs_euc ~ Group_Time_Subs, perm = 999) 
adon.subs.euc # F = 0.14, R2 = 0.00083

# p> 0.05 for all. The two groups are very similar


#### Greatest distance models

# Extract vector information
pcoa.subs.jac <- pcoa(dist_subs_jacbi)
pcoa.subs.bray <- pcoa(dist_subs_BC)



# Extract %
Axis1.pcoa.subs.jac <- pcoa.subs.jac$values$Relative_eig[[1]] * 100 # Dimension (i.e., Axis 1 (PCOA1)) - 12.52%
Axis2.pcoa.subs.jac <- pcoa.subs.jac$values$Relative_eig[[2]] * 100 # Dimension (i.e., Axis 2 (PCOA2)) - 7%

Axis1.pcoa.subs.bc <- pcoa.subs.bray$values$Relative_eig[[1]] * 100 # Dimension (i.e., Axis 1 (PCOA1)) - 13.83%
Axis2.pcoa.subs.bc <- pcoa.subs.bray$values$Relative_eig[[2]] * 100 # Dimension (i.e., Axis 2 (PCOA2)) - 10.38%


# Convert to dataframe
jac.subs.df <- data.frame(SiteSubsiteYear = as.character(rownames(pcoa.subs.jac$vectors)),
                     X = pcoa.subs.jac$vectors[, 1], Y = pcoa.subs.jac$vectors[, 2])

bc.subs.df <- data.frame(SiteSubsiteYear = as.character(rownames(pcoa.subs.bray$vectors)),
                          X = pcoa.subs.bray$vectors[, 1], Y = pcoa.subs.bray$vectors[, 2])

# Convert to wide format
jac.subs.df.wide <- jac.subs.df %>% mutate(TimePoint = rep(c("Start", "End"), 90)) %>% 
  mutate(SiteSubsite = case_when(TimePoint == "Start" ~ str_sub(SiteSubsiteYear, end = -7),
                                 TimePoint == "End" ~ str_sub(SiteSubsiteYear, end = -5)))
rownames(jac.subs.df.wide) <- NULL

write.csv(jac.subs.df.wide, "data/pcoa_jac_subs_scores.csv")

bc.subs.df.wide <- bc.subs.df %>% mutate(TimePoint = rep(c("Start", "End"), 90)) %>% 
  mutate(SiteSubsite = case_when(TimePoint == "Start" ~ str_sub(SiteSubsiteYear, end = -7),
                                 TimePoint == "End" ~ str_sub(SiteSubsiteYear, end = -5)))
rownames(bc.subs.df.wide) <- NULL

write.csv(bc.subs.df.wide, "data/pcoa_bray_subs_scores.csv")

# Calculate Cartesian distances
jac.subs.start <- jac.subs.df.wide %>% filter(TimePoint == "Start") %>% rename(X1 = X, Y1 = Y) %>%
  select(-c(SiteSubsiteYear, TimePoint))

jac.subs.end <- jac.subs.df.wide %>% filter(TimePoint == "End") %>% rename(X2 = X, Y2 = Y) %>%
  select(-c(SiteSubsiteYear, TimePoint))

jac.subs.all <- left_join(jac.subs.start, jac.subs.end, by = "SiteSubsite")

jacxy.subs <- jac.subs.all %>% mutate(CoordDistance = sqrt((X2 - X1)^2 + (Y2 - Y1)^2))
write.csv(jacxy.subs, "data/22jacxy_subs.csv")
# Which are the plots with the greatest distance? 
# BYLOT:MESPOLYGON, ZACKENBERG:ABRASION PLATEAU, AUYUITTUQ:OWL RIVER, 



bc.subs.start <- bc.subs.df.wide %>% filter(TimePoint == "Start") %>% rename(X1 = X, Y1 = Y) %>%
  select(-c(SiteSubsiteYear, TimePoint))

bc.subs.end <- bc.subs.df.wide %>% filter(TimePoint == "End") %>% rename(X2 = X, Y2 = Y) %>%
  select(-c(SiteSubsiteYear, TimePoint))

bc.subs.all <- left_join(bc.subs.start, bc.subs.end, by = "SiteSubsite")

bcxy.subs <- bc.subs.all %>% mutate(CoordDistance = sqrt((X2 - X1)^2 + (Y2 - Y1)^2))
write.csv(bcxy.subs, "data/22bcxy_subs.csv")
# Which are the plots with the greatest distance? 
# ABISKO:PEATLAND, LATNJA:TUSSOCK_TUNDRA, ZACKENBERG:FEN



## Add in climate data: need to calculate mean temp/prec change at the subsite level
clim.all <- read.csv("data/22clim_slopes.csv")

clim.all.subs <- clim.all %>% separate(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), ":") %>%
  unite("SiteSubsite", SITE:SUBSITE, sep = ":") %>% group_by(SiteSubsite) %>% 
  mutate(SubsiteTempChange = mean(WarmQSlope), SubsitePrecChange = mean(PrecSlope)) %>%
  ungroup() %>% distinct(SiteSubsite, .keep_all = TRUE) %>% 
  select(SiteSubsite, SubsiteTempChange, SubsitePrecChange)
  
## Add in FG data: need to calculate mean FG change at the subsite level
fg.slopes.wide <- read.csv("data/22fg_slopes.csv")

fg.slopes.subs <- fg.slopes.wide %>% separate(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), ":") %>%
  unite("SiteSubsite", SITE:SUBSITE, sep = ":") %>% group_by(SiteSubsite) %>% 
  mutate(SubsiteForbChange = mean(ForbSlope), SubsiteShrubChange = mean(ShrubSlope), 
         SubsiteGraminoidChange = mean(GraminoidSlope)) %>%
  ungroup() %>% distinct(SiteSubsite, .keep_all = TRUE) %>% 
  select(SiteSubsite, SubsiteForbChange, SubsiteShrubChange, SubsiteGraminoidChange)



## Extract rest of metadata - calculate mean for all explanatory variables
itex.dec22$SurveyedArea <- as.numeric(itex.dec22$SurveyedArea)

metadata.subs <- itex.dec22 %>% group_by(SiteSubsitePlot) %>% 
  mutate(Duration = max(YEAR) - min(YEAR)) %>% ungroup() %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>% group_by(SiteSubsite) %>%
  mutate(SubsiteLatitude = mean(LAT), SubsiteSurveyedArea = mean(SurveyedArea), 
         SubsiteRichness = mean(MeanRichness), SubsiteDuration = mean(Duration)) %>%
  distinct(SiteSubsite, .keep_all = TRUE) %>% 
  select(SiteSubsite, SubsiteLatitude, MOISTURE, Region, SubsiteSurveyedArea, SubsiteRichness, SubsiteDuration) %>%
  mutate(LogPlotSize = log(SubsiteSurveyedArea)) %>% filter(MOISTURE !="")

# Check that all plots within a subsite have the same moisture value - all OK
moisture.check <- metadata.subs %>% group_by(SiteSubsite) %>% 
  summarise(MOISTURE) %>% distinct()
  
# Combine data for modeling
jac.subs.meta <- left_join(jacxy.subs, fg.slopes.subs, by = "SiteSubsite")
jac.subs.meta2 <- left_join(jac.subs.meta, clim.all.subs, by = "SiteSubsite")
jac.subs.full <- left_join(jac.subs.meta2, metadata.subs, by = "SiteSubsite")

hist(jac.subs.full$CoordDistance) # beta
str(jac.subs.full) 


bc.subs.meta <- left_join(bcxy.subs, fg.slopes.subs, by = "SiteSubsite")
bc.subs.meta2 <- left_join(bc.subs.meta, clim.all.subs, by = "SiteSubsite")
bc.subs.full <- left_join(bc.subs.meta2, metadata.subs, by = "SiteSubsite")

hist(bc.subs.full$CoordDistance) # beta
str(bc.subs.full) 



#### 01) JACCARD ####
jac.subs.mod <- brm(bf(CoordDistance ~ Region + MOISTURE + SubsiteLatitude + SubsiteTempChange + SubsitePrecChange + 
                         LogPlotSize + SubsiteRichness + SubsiteDuration + SubsiteShrubChange),
                       data = jac.subs.full, family = "beta",
                       iter = 2000, chains = 4, warmup = 400, 
                       file = "models/24_01_jac_subs_shb_mod")
summary(jac.subs.mod) # lat ns, temp positive, prec ns, plot size ns, richness ns, duration ns, shrub ns
conditional_effects(jac.subs.mod) # moisture ns, NAm-east less change than Eurasia

print(summary(jac.subs.mod, prob = 0.975)) # lat ns, temp positive, prec ns, plot size ns, richness ns, duration ns, shrub ns
conditional_effects(jac.subs.mod, prob = 0.975) # moisture ns, region ns

# forb
jac.subs.forb.mod <- brm(bf(CoordDistance ~ Region + MOISTURE + SubsiteLatitude + SubsiteTempChange + SubsitePrecChange + 
                         LogPlotSize + SubsiteRichness + SubsiteDuration + SubsiteForbChange),
                    data = jac.subs.full, family = "beta",
                    iter = 2000, chains = 4, warmup = 400, 
                    file = "models/24_01_jac_subs_frb_mod")
summary(jac.subs.forb.mod) # lat ns, temp ns, prec ns, plot size positive, richness ns, duration ns, forb ns
conditional_effects(jac.subs.forb.mod) # region ns, moisture ns 

print(summary(jac.subs.forb.mod, prob = 0.975)) # lat ns, temp ns, prec ns, plot size positive, richness ns, duration ns, forb ns
conditional_effects(jac.subs.forb.mod, prob = 0.975) # moisture ns, region ns


# gram
jac.subs.gram.mod <- brm(bf(CoordDistance ~ Region + MOISTURE + SubsiteLatitude + SubsiteTempChange + SubsitePrecChange + 
                              LogPlotSize + SubsiteRichness + SubsiteDuration + SubsiteGraminoidChange),
                         data = jac.subs.full, family = "beta",
                         iter = 2000, chains = 4, warmup = 400, 
                         file = "models/24_01_jac_subs_gram_mod")
summary(jac.subs.gram.mod) # lat ns, temp ns, prec ns, plot size ns, richness positive, duration negative, gram ns
conditional_effects(jac.subs.gram.mod) # region ns, moisture ns 

print(summary(jac.subs.gram.mod, prob = 0.975)) # lat ns, temp ns, prec ns, plot size ns, richness positive, duration negative, gram ns
conditional_effects(jac.subs.gram.mod, prob = 0.975) # moisture ns, region ns


#### 02) BRAY-CURTIS ####

# shrub
bray.subs.mod <- brm(CoordDistance ~ Region + MOISTURE + SubsiteLatitude + SubsiteTempChange + SubsitePrecChange + 
                         LogPlotSize + SubsiteRichness + SubsiteDuration + SubsiteShrubChange,
                    data = bc.subs.full, family = "beta",
                    iter = 2000, chains = 4, warmup = 400, 
                    file = "models/24_01_bray_subs_shb_mod")
summary(bray.subs.mod) # lat ns, temp ns, prec ns, plot size ns, richness ns, duration ns, shrub ns
conditional_effects(bray.subs.mod) # moisture ns, region ns

print(summary(bray.subs.mod, prob = 0.975)) # lat ns, temp ns, prec ns, plot size ns, richness ns, duration ns, shrub ns
conditional_effects(bray.subs.mod, prob = 0.975) # moisture ns, region ns


# forb
bray.subs.forb.mod <- brm(CoordDistance ~ Region + MOISTURE + SubsiteLatitude + SubsiteTempChange + SubsitePrecChange + 
                       LogPlotSize + SubsiteRichness + SubsiteDuration + SubsiteForbChange,
                     data = bc.subs.full, family = "beta",
                     iter = 2000, chains = 4, warmup = 400, 
                     file = "models/24_01_bray_subs_forb_mod")
summary(bray.subs.forb.mod) # lat ns, temp ns, prec ns, plot size ns, richness ns, duration ns, forb ns
conditional_effects(bray.subs.forb.mod) # moisture ns, region ns

print(summary(bray.subs.forb.mod, prob = 0.975)) # lat ns, temp ns, prec ns, plot size ns, richness ns, duration ns, forb ns
conditional_effects(bray.subs.forb.mod, prob = 0.975) # moisture ns, region ns


# gram
bray.subs.gram.mod <- brm(CoordDistance ~ Region + MOISTURE + SubsiteLatitude + SubsiteTempChange + SubsitePrecChange + 
                            LogPlotSize + SubsiteRichness + SubsiteDuration + SubsiteGraminoidChange,
                          data = bc.subs.full, family = "beta",
                          iter = 2000, chains = 4, warmup = 400, 
                          file = "models/24_01_bray_subs_gram_mod")
print(summary(bray.subs.gram.mod), digits = 5) # lat ns, temp ns, prec ns, plot size ns, richness ns, duration ns, gram ns
conditional_effects(bray.subs.gram.mod) # moisture ns, region ns

print(summary(bray.subs.gram.mod, prob = 0.975)) # lat ns, temp ns, prec ns, plot size ns, richness ns, duration ns, gram ns
conditional_effects(bray.subs.gram.mod, prob = 0.975) # moisture ns, region ns




#### CENTROID ANALYSES
beta_subs_BC

# For the same plot, different distances to centroid for start/end plot
JacDisSubs <- as.data.frame(beta_subs_jacbi$distances)
SorDisSubs <- as.data.frame(beta_subs_S$distances)
BrayDisSubs <- as.data.frame(beta_subs_BC$distances)
AgwDisSubs <- as.data.frame(beta_subs_aGw$distances)
ManDisSubs <- as.data.frame(beta_subs_man$distances)
EucDisSubs <- as.data.frame(beta_subs_euc$distances)


# Calculate differences in distance to centroid for each plot (compare end - start distance).

## Jaccard
JacDisSubs2 <- JacDisSubs %>% rownames_to_column() %>% 
  tidyr::separate(rowname, c("Site", "Subsite", "Timepoint"), sep = ":") %>% 
  rename("DistanceToCentroid" = "beta_subs_jacbi$distances") %>%
  mutate(SiteSubsite = paste(Site, Subsite, sep = ":")) %>%
  group_by(SiteSubsite) %>%
  mutate(DistanceDifference = diff(DistanceToCentroid)) %>%
  mutate(AbsDistanceDifference = abs(DistanceDifference)) %>%
  distinct(SiteSubsite, .keep_all = TRUE) %>%
  select(SiteSubsite, AbsDistanceDifference)
# Greatest changes in presence-absence composition: AUYUITTUQ:OWL RIVER, DISKO:DISTURBANCE, ABISKO:ABISKOWET_MULTE


# Combine with other data for models
jac.dis.subs.meta <- left_join(JacDisSubs2, fg.slopes.subs, by = "SiteSubsite")
jac.dis.subs.meta2 <- left_join(jac.dis.subs.meta, clim.all.subs, by = "SiteSubsite")
jac.dis.subs.meta3 <- left_join(jac.dis.subs.meta2, metadata.subs, by = "SiteSubsite")

# Bit of admin
hist(jac.dis.subs.meta3$AbsDistanceDifference) # beta dist
str(jac.dis.subs.meta3)



## Bray-Curtis
BrayDisSubs2 <- BrayDisSubs %>% rownames_to_column() %>% 
  tidyr::separate(rowname, c("Site", "Subsite", "Timepoint"), sep = ":") %>% 
  rename("DistanceToCentroid" = "beta_subs_BC$distances") %>%
  mutate(SiteSubsite = paste(Site, Subsite, sep = ":")) %>%
  group_by(SiteSubsite) %>%
  mutate(DistanceDifference = diff(DistanceToCentroid)) %>%
  mutate(AbsDistanceDifference = abs(DistanceDifference)) %>%
  distinct(SiteSubsite, .keep_all = TRUE) %>%
  select(SiteSubsite, AbsDistanceDifference)



# Combine with other data for models
bc.dis.subs.meta <- left_join(BrayDisSubs2, fg.slopes.subs, by = "SiteSubsite")
bc.dis.subs.meta2 <- left_join(bc.dis.subs.meta, clim.all.subs, by = "SiteSubsite")
bc.dis.subs.meta3 <- left_join(bc.dis.subs.meta2, metadata.subs, by = "SiteSubsite")

# Bit of admin
hist(bc.dis.subs.meta3$AbsDistanceDifference) # beta dist
str(bc.dis.subs.meta3)


#### 01) CENTROID JACCARD ####
jac.cent.subs.mod <- brm(AbsDistanceDifference ~ Region + MOISTURE + SubsiteLatitude + SubsiteTempChange + SubsitePrecChange + 
                              SubsiteShrubChange + LogPlotSize + SubsiteRichness + SubsiteDuration,
                    data = jac.dis.subs.meta3, family = "beta",
                    iter = 2000, chains = 4, warmup = 400, 
                    file = "models/24_01_jac_cent_subs_shb_mod")
print(summary(jac.cent.subs.mod), digits = 5) # lat ns, temp ns, prec ns, shrub ns, plot size ns, richness ns, duration ns
conditional_effects(jac.cent.subs.mod) # moisture ns, region ns

print(summary(jac.cent.subs.mod, prob = 0.975)) # lat ns, temp ns, prec ns, shrub ns, plot size ns, richness ns, duration ns
conditional_effects(jac.cent.subs.mod, prob = 0.975) # moisture ns, region ns


# forb
jac.cent.subs.forb.mod <- brm(AbsDistanceDifference ~ Region + MOISTURE + SubsiteLatitude + SubsiteTempChange + SubsitePrecChange + 
                           SubsiteForbChange + LogPlotSize + SubsiteRichness + SubsiteDuration,
                         data = jac.dis.subs.meta3, family = "beta",
                         iter = 2000, chains = 4, warmup = 400, 
                         file = "models/24_01_jac_cent_subs_forb_mod")
summary(jac.cent.subs.forb.mod) # lat ns, temp ns, prec ns, forb ns, plot size ns, richness ns, duration ns
conditional_effects(jac.cent.subs.forb.mod) # moisture ns, region ns

print(summary(jac.cent.subs.forb.mod, prob = 0.975)) # lat ns, temp ns, prec ns, forb ns, plot size ns, richness ns, duration ns
conditional_effects(jac.cent.subs.forb.mod, prob = 0.975) # moisture ns, region ns



# gram
jac.cent.subs.gram.mod <- brm(AbsDistanceDifference ~ Region + MOISTURE + SubsiteLatitude + SubsiteTempChange + SubsitePrecChange + 
                                SubsiteGraminoidChange + LogPlotSize + SubsiteRichness + SubsiteDuration,
                              data = jac.dis.subs.meta3, family = "beta",
                              iter = 2000, chains = 4, warmup = 400, 
                              file = "models/24_01_jac_cent_subs_gram_mod")
summary(jac.cent.subs.gram.mod) # lat ns, temp ns, prec ns, gram ns, plot size ns, richness ns, duration negative
conditional_effects(jac.cent.subs.gram.mod) # moisture ns, region ns

print(summary(jac.cent.subs.gram.mod, prob = 0.975)) # lat ns, temp ns, prec ns, gram ns, plot size ns, richness ns, duration negative
conditional_effects(jac.cent.subs.gram.mod, prob = 0.975) # moisture ns, region ns




#### 02) CENTROID BRAY ####
bc.cent.subs.mod <- brm(AbsDistanceDifference ~ Region + MOISTURE + SubsiteLatitude + SubsiteTempChange + SubsitePrecChange + 
                           SubsiteShrubChange + LogPlotSize + SubsiteRichness + SubsiteDuration,
                         data = bc.dis.subs.meta3, family = "beta",
                         iter = 2000, chains = 4, warmup = 400, 
                         file = "models/24_02_bc_cent_subs_shb_mod")
summary(bc.cent.subs.mod) # lat ns, temp ns, prec ns, shrub ns, plot size ns, richness ns, duration ns
conditional_effects(bc.cent.subs.mod) # moisture ns, region ns

print(summary(bc.cent.subs.mod, prob = 0.975)) # lat ns, temp ns, prec ns, shrub ns, plot size ns, richness ns, duration ns
conditional_effects(bc.cent.subs.mod, prob = 0.975) # moisture ns, region ns


# forb
bc.cent.subs.forb.mod <- brm(AbsDistanceDifference ~ Region + MOISTURE + SubsiteLatitude + SubsiteTempChange + SubsitePrecChange + 
                                SubsiteForbChange + LogPlotSize + SubsiteRichness + SubsiteDuration,
                              data = bc.dis.subs.meta3, family = "beta",
                              iter = 2000, chains = 4, warmup = 400, 
                              file = "models/24_02_bc_cent_subs_forb_mod")
summary(bc.cent.subs.forb.mod) # lat ns, temp ns, prec ns, forb ns, plot size ns, richness ns, duration ns
conditional_effects(bc.cent.subs.forb.mod) # moisture ns, region ns

print(summary(bc.cent.subs.forb.mod, prob = 0.975)) # lat ns, temp ns, prec ns, forb ns, plot size ns, richness ns, duration ns
conditional_effects(bc.cent.subs.forb.mod, prob = 0.975) # moisture ns, region ns



# gram
bc.cent.subs.gram.mod <- brm(AbsDistanceDifference ~ Region + MOISTURE + SubsiteLatitude + SubsiteTempChange + SubsitePrecChange + 
                                SubsiteGraminoidChange + LogPlotSize + SubsiteRichness + SubsiteDuration,
                              data = bc.dis.subs.meta3, family = "beta",
                              iter = 2000, chains = 4, warmup = 400, 
                              file = "models/24_02_bc_cent_subs_gram_mod")
summary(bc.cent.subs.gram.mod) # lat ns, temp ns, prec ns, gram ns, plot size ns, richness ns, duration ns
conditional_effects(bc.cent.subs.gram.mod) # moisture ns, region ns

print(summary(bc.cent.subs.gram.mod, prob = 0.975)) # lat ns, temp ns, prec ns, gram ns, plot size ns, richness ns, duration ns
conditional_effects(bc.cent.subs.gram.mod, prob = 0.975) # moisture ns, region ns




## FIGURE 4 ----

#### Jaccard

# Make disimilarity-matrix
dist_OldNew_J <- vegdist(abundance.subs.wide, method = "jaccard")

# Jaccard: extract variance explained by PCOA tests
pcoa.j.subsite <- pcoa(dist_OldNew_J)
barplot(pcoa.j.subsite$values$Relative_eig[1:10])

# Now properly calculated
Axis1.pcoa.j.subsite <- pcoa.j.subsite$values$Relative_eig[[1]] * 100 # Dimension (i.e., Axis 1 (PCOA1)) - 14.04%
Axis2.pcoa.j.subsite <- pcoa.j.subsite$values$Relative_eig[[2]] * 100 # Dimension (i.e., Axis 2 (PCOA2)) - 10.68%

# Make time vector to be grouped by: each plot is old/new time
time_vector <- c("Begin", "End")
Group_Time <- rep(time_vector, length(unique(abundance.subs$SiteSubsite)))

# Calculate mean distance to centroid for groups: homogeneity of variance
beta_OldNew_J <- betadisper(dist_OldNew_J, group = Group_Time)

labs.jacc.subsite <- paste0("Dimension ", 1:4, " (", round(100*beta_OldNew_J$eig / sum(beta_OldNew_J$eig), 2), "%)")

# List of additional parameters
str(beta_OldNew_J)

# values of PCoA1 and PCoA2 per site/year and also the centroid values for PCoA1 and PCoA2
centroids_J <- as.data.frame(scores(beta_OldNew_J)$centroids) %>%
  tibble::rownames_to_column("Timepoint") %>%
  rename(PCoA1_centre = PCoA1, PCoA2_centre = PCoA2)

centroids_J_wide <- centroids_J %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = c(PCoA1_centre, PCoA2_centre)
  )

JaccDisTime <- as.data.frame(beta_OldNew_J$distances)

JaccScoresTime <- scores(beta_OldNew_J)

# extract coordinate points
pcoa.jacc.subsite.coords <- as.data.frame(JaccScoresTime$sites) %>%
  tibble::rownames_to_column(var = "SiteSubsite_row") %>%
  rename(SiteSubsite = SiteSubsite_row) %>%
  separate(SiteSubsite, into = c("Site", "Subsite", "Timepoint"), sep = ":", fill = "right") %>%
  unite(SiteSubsite, Site, Subsite, sep = ":", remove = TRUE) %>%
  unite(SiteSubsiteTime, SiteSubsite, Timepoint, sep = ":", remove = FALSE)

# merge with climate data
pcoa.jacc.subsite.coords.warm <- left_join(pcoa.jacc.subsite.coords, clim.change.subsite, by = "SiteSubsite", relationship = "many-to-many")

# reformat so we have a row per subsite
jacc.pcoa.subsite.arrow <- pcoa.jacc.subsite.coords.warm %>%
  pivot_wider(names_from = Timepoint, values_from = c(PCoA1, PCoA2)) 

# collapse so we have one row per subsite
jacc.pcoa.subsite.df <- jacc.pcoa.subsite.arrow %>% group_by(SiteSubsite) %>% 
  summarise_each(funs(first(.[!is.na(.)]))) %>% ungroup() %>% 
  dplyr::select(-SiteSubsiteTime)

# create convex hull
hull.jacc.pcoa.subsite <- pcoa.jacc.subsite.coords.warm %>%
  group_by(Timepoint) %>%
  slice(chull(PCoA1, PCoA2))

# Add latitude as a categorical variable
jacc.pcoa.subsite.df <- jacc.pcoa.subsite.df %>% 
  mutate(LatCat = case_when(LAT < 68.35765 ~ "Low",
                            LAT > 68.35765 & LAT < 71.31187 ~ "Mid",
                            LAT > 71.31187 ~ "High")) %>%
  cross_join(centroids_J_wide)


hull.jacc.pcoa.subsite <- hull.jacc.pcoa.subsite %>%
  mutate(LatCat = case_when(LAT < 68.35765 ~ "Low",
                            LAT > 68.35765 & LAT < 71.31187 ~ "Mid",
                            LAT > 71.31187 ~ "High")) %>%
  inner_join(centroids_J)


# long format so shape can also appear in the legend
jacc.pcoa.subsite.df.long.pcoa1 <- jacc.pcoa.subsite.df %>% 
  tidyr::pivot_longer(cols = c(PCoA1_Begin, PCoA1_End),
                      names_to = "Timepoint",
                      values_to = "PCoA1") %>% 
  mutate(Timepoint = sub('PCoA1_', '', Timepoint)) 

# PCOA2
jacc.pcoa.subsite.df.long.pcoa2 <- jacc.pcoa.subsite.df %>%
  tidyr::pivot_longer(cols = c(PCoA2_Begin, PCoA2_End),
                      names_to = "Timepoint",
                      values_to = "PCoA2") %>%
  mutate(Timepoint = sub('PCoA2_', '', Timepoint)) %>%
  dplyr::select(SiteSubsite, Timepoint, PCoA2)

# join in one df
jac.pcoa.long <- left_join(jacc.pcoa.subsite.df.long.pcoa1, jacc.pcoa.subsite.df.long.pcoa2, 
                           by = c("SiteSubsite", "Timepoint"))

# plot PCoA w ggplot and arrows between the points
(jacc.pcoa.subsite.plot <- ggplot() + 
   geom_polygon(data = hull.jacc.pcoa.subsite, aes(PCoA2, PCoA1, colour = Timepoint, fill = Timepoint), alpha = 0.05, size = 1) +
   scale_color_manual(values = c("#F4B043", "#219EAE"), name = "Time point") +
   scale_fill_manual(values = c("#F4B043", "#219EAE"), name = "Time point") +
   new_scale_colour() +
   new_scale_fill() +
   geom_point(data = jac.pcoa.long, aes(x = PCoA2, y = PCoA1, colour = LAT, shape = Timepoint), alpha = 0.9, size = 5) +
   scale_shape_manual(values = c(17, 16), name = "Time point") +
   geom_segment(data = jacc.pcoa.subsite.df, aes(x = PCoA2_Begin, y = PCoA1_Begin, xend = PCoA2_End, yend = PCoA1_End), colour = "black", linewidth = 0.8, arrow = arrow(length = unit(0.3, "cm"))) +
   scale_colour_gradient(low = "#E85C90", high = "#3CD1CF", name = "Latitude (°)") +
   ylab("PCoA Axis 1<br>") +
   xlab("<br>PCoA Axis 2") +
   xlim(-0.5, 0.5) +
   ylim(-0.5, 0.5) +
   bio.theme +
   theme(legend.position = "right", text = element_text(size = 24), 
         panel.background = element_blank(), 
         panel.border = element_rect(colour = "black", fill = NA), 
         panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
         panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
         axis.title.x = element_markdown(size = 24), 
         axis.title.y = element_markdown(size = 24),
         legend.title = element_markdown(size = 20),
         legend.key = element_rect(fill = NA, colour = NA))) +
  guides(colour = guide_legend(order = 1), fill = guide_legend(order = 2))

ggplot2::ggsave(jacc.pcoa.subsite.plot, filename = "figures/Figure_4a.png", width = 24, height = 18, units = "cm")


#### Bray-Curtis

# Make disimilarity-matrix - method can be changed (for options see: ?vegdist)
# jaccard =  _jacbi; bray bi = Sorensen;  altGower =  _aGw: manhattan =  _man; euclidian =  _euc
dist_OldNew_BC <- vegdist(abundance.subs.wide, method = "bray")

# Bray-Curtis: extract variance explained by PCOA tests
pcoa.bc.subsite <- pcoa(dist_OldNew_BC)
barplot(pcoa.bc.subsite$values$Relative_eig[1:10])

# Now properly calculated
Axis1.pcoa.bc.subsite <- pcoa.bc.subsite$values$Relative_eig[[1]] * 100 # Dimension (i.e., Axis 1 (PCOA1)) - 14.04%
Axis2.pcoa.bc.subsite <- pcoa.bc.subsite$values$Relative_eig[[2]] * 100 # Dimension (i.e., Axis 2 (PCOA2)) - 10.68%

# Make time vector to be grouped by: each plot is old/new time
time_vector <- c("Begin", "End")
Group_Time <- rep(time_vector, length(unique(abundance.subs$SiteSubsite)))

# Calculate mean distance to centroid for groups: homogeneity of variance - see ?betadisper
beta_OldNew_BC <- betadisper(dist_OldNew_BC, group = Group_Time)

labs.bray.subsite <- paste0("Dimension ", 1:4, " (", round(100*beta_OldNew_BC$eig / sum(beta_OldNew_BC$eig), 2), "%)")

# List of additional parameters
str(beta_OldNew_BC)

# values of PCoA1 and PCoA2 per site/year and also the centroid values for PCoA1 and PCoA2
centroids_BC <- as.data.frame(scores(beta_OldNew_BC)$centroids) %>%
  tibble::rownames_to_column("Timepoint") %>%
  rename(PCoA1_centre = PCoA1, PCoA2_centre = PCoA2)

centroids_BC_wide <- centroids_BC %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = c(PCoA1_centre, PCoA2_centre)
  )

BrayDisTime <- as.data.frame(beta_OldNew_BC$distances)

BrayScoresTime <- scores(beta_OldNew_BC)

# extract coordinate points
pcoa.bray.subsite.coords <- as.data.frame(BrayScoresTime$sites) %>%
  tibble::rownames_to_column(var = "SiteSubsite_row") %>%
  rename(SiteSubsite = SiteSubsite_row) %>%
  separate(SiteSubsite, into = c("Site", "Subsite", "Timepoint"), sep = ":", fill = "right") %>%
  unite(SiteSubsite, Site, Subsite, sep = ":", remove = TRUE) %>%
  unite(SiteSubsiteTime, SiteSubsite, Timepoint, sep = ":", remove = FALSE)

# merge with climate data
pcoa.bray.subsite.coords.warm <- left_join(pcoa.bray.subsite.coords, clim.change.subsite, by = "SiteSubsite", relationship = "many-to-many")

# reformat so we have a row per subsite
bray.pcoa.subsite.arrow <- pcoa.bray.subsite.coords.warm %>%
  pivot_wider(names_from = Timepoint, values_from = c(PCoA1, PCoA2)) 

# collapse so we have one row per subsite
bray.pcoa.subsite.df <- bray.pcoa.subsite.arrow %>% group_by(SiteSubsite) %>% 
  summarise_each(funs(first(.[!is.na(.)]))) %>% ungroup() %>% 
  dplyr::select(-SiteSubsiteTime)

# create convex hull
hull.bray.pcoa.subsite <- pcoa.bray.subsite.coords.warm %>%
  group_by(Timepoint) %>%
  slice(chull(PCoA1, PCoA2))

# Add latitude as a categorical variable
bray.pcoa.subsite.df <- bray.pcoa.subsite.df %>% 
  mutate(LatCat = case_when(LAT < 68.35765 ~ "Low",
                            LAT > 68.35765 & LAT < 71.31187 ~ "Mid",
                            LAT > 71.31187 ~ "High")) %>%
  cross_join(centroids_BC_wide)


hull.bray.pcoa.subsite <- hull.bray.pcoa.subsite %>%
  mutate(LatCat = case_when(LAT < 68.35765 ~ "Low",
                            LAT > 68.35765 & LAT < 71.31187 ~ "Mid",
                            LAT > 71.31187 ~ "High")) %>%
  inner_join(centroids_BC)


# long format so shape can also appear in the legend
bc.pcoa.subsite.df.long.pcoa1 <- bray.pcoa.subsite.df %>% 
  tidyr::pivot_longer(cols = c(PCoA1_Begin, PCoA1_End),
                      names_to = "Timepoint",
                      values_to = "PCoA1") %>% 
  mutate(Timepoint = sub('PCoA1_', '', Timepoint)) 


# PCOA 2
bc.pcoa.subsite.df.long.pcoa2 <- bray.pcoa.subsite.df %>%
  tidyr::pivot_longer(cols = c(PCoA2_Begin, PCoA2_End),
                      names_to = "Timepoint",
                      values_to = "PCoA2") %>%
  mutate(Timepoint = sub('PCoA2_', '', Timepoint)) %>%
  dplyr::select(SiteSubsite, Timepoint, PCoA2)

# join in one df
bc.pcoa.long <- left_join(bc.pcoa.subsite.df.long.pcoa1, bc.pcoa.subsite.df.long.pcoa2, 
                          by = c("SiteSubsite", "Timepoint"))



# plot PCoA w ggplot and arrows between the points
(bray.pcoa.subsite.plot <- ggplot() + 
   geom_polygon(data = hull.bray.pcoa.subsite, aes(PCoA2, PCoA1, colour = Timepoint, fill = Timepoint), alpha = 0.05, size = 1) +
   scale_color_manual(values = c("#F4B043", "#219EAE"), name = "Time point") +
   scale_fill_manual(values = c("#F4B043", "#219EAE"), name = "Time point") +
   #scale_y_reverse(limits = c(0.5, -0.5)) +
   new_scale_colour() +
   new_scale_fill() +
   geom_point(data = bc.pcoa.long, aes(x = PCoA2, y = PCoA1, colour = LAT, shape = Timepoint), alpha = 0.9, size = 5) +
   scale_shape_manual(values = c(17, 16), name = "Time point") +
   geom_segment(data = bray.pcoa.subsite.df, aes(x = PCoA2_Begin, y = PCoA1_Begin, xend = PCoA2_End, yend = PCoA1_End, colour = LAT), colour = "black", linewidth = 0.8, arrow = arrow(length = unit(0.3, "cm"))) +
   scale_colour_gradient(low = "#E85C90", high = "#3CD1CF", name = "Latitude (°)") +
   ylab("PCoA Axis 1<br>") +
   xlab("<br>PCoA Axis 2") +
   xlim(-0.5, 0.5) +
   ylim(-0.5, 0.5) +
   bio.theme +
   theme(legend.position = "right", text = element_text(size = 24), 
         panel.background = element_blank(), 
         panel.border = element_rect(colour = "black", fill = NA), 
         panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
         panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
         axis.title.x = element_markdown(size = 24), 
         axis.title.y = element_markdown(size = 24),
         legend.title = element_markdown(size = 20),
         legend.key = element_rect(fill = NA, colour = NA))) +
  guides(colour = guide_legend(order = 1), fill = guide_legend(order = 2))


ggplot2::ggsave(bray.pcoa.subsite.plot, filename = "figures/Figure_4b.png", width = 24, height = 18, units = "cm")





## Distance to centroid boxplots
jac.dist.cent <- read.csv("data/pcoa_jac_dist_centroid.csv")

# Change to start and order so that boxplot appears first
jac.dist.cent2 <- jac.dist.cent %>% mutate(Timepoint = case_when(Timepoint == "Begin" ~ "Start", 
                                                                 TRUE ~ "End"))

jac.dist.cent2$Timepoint <- factor(jac.dist.cent2$Timepoint, levels=c('Start', 'End'), ordered = TRUE)

(jac.dist.box <- ggplot(jac.dist.cent2, aes(Timepoint, DistanceToCentroid)) + 
    stat_boxplot(geom = "errorbar", width = 0.2) + 
    geom_boxplot(aes(fill = Timepoint), notch = TRUE) +
    scale_fill_manual(values = c("#F4B043", "#219EAE"), labels = c("Start", "End")) +
    theme_bw() + 
    xlab("\nTime point") + ylab("Distance to centroid\n") +
    ylim(0.5, 0.8) +
    theme(axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 24),
          axis.text = element_text(size = 20, colour = "black"),
          legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()))

ggsave(jac.dist.box, filename = "figures/jac_dist_box.png", 
       width = 16, height = 18, units = "cm")

# Mean distance to centroid
meanjacst <- jac.dist.cent2 %>% group_by(Timepoint) %>% 
  summarise(mean(DistanceToCentroid), sd(DistanceToCentroid))

mean(jac.dist.cent$DistanceToCentroid) # 0.657 + 0.03
sd(jac.dist.cent$DistanceToCentroid)


## Distance to centroid boxplots
bc.dist.cent <- read.csv("data/pcoa_bray_dist_centroid.csv")

# Change to start and order so that boxplot appears first
bc.dist.cent2 <- bc.dist.cent %>% 
  mutate(Timepoint = case_when(Timepoint == "Begin" ~ "Start", TRUE ~ "End"))

bc.dist.cent2$Timepoint <- factor(bc.dist.cent2$Timepoint, levels=c('Start', 'End'), ordered = TRUE)

(bc.dist.box <- ggplot(bc.dist.cent2, aes(Timepoint, DistanceToCentroid)) + 
    stat_boxplot(geom = "errorbar", width = 0.2) + 
    geom_boxplot(aes(fill = Timepoint), notch = TRUE) +
    scale_fill_manual(values = c("#F4B043", "#219EAE"), labels = c("Start", "End")) +
    theme_bw() + 
    #ggtitle("(d) Bray-Curtis PCoA scores") + 
    xlab("\nTime point") + ylab("Distance to centroid\n") +
    ylim(0.5, 0.8) +
    theme(plot.title = element_text(size = 20), 
          axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 24),
          axis.text = element_text(size = 20, colour = "black"),
          legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()))

ggsave(bc.dist.box, filename = "figures/Figure_4d.png", 
       width = 16, height = 18, units = "cm")


# Mean distance to centroid
meanbc <- bc.dist.cent2 %>% group_by(Timepoint) %>% 
  summarise(mean(DistanceToCentroid), sd(DistanceToCentroid))


## Distance between start and end time points
jacxy.subs <- read.csv("data/22jacxy_subs.csv")
bcxy.subs <- read.csv("data/22bcxy_subs.csv")

jacxy.subs2 <- jacxy.subs %>% dplyr::select(SiteSubsite, CoordDistance) %>% 
  mutate(Metric = "Jaccard")

bcxy.subs2 <- bcxy.subs %>% dplyr::select(SiteSubsite, CoordDistance) %>% 
  mutate(Metric = "Bray-Curtis")

mean(jacxy.subs2$CoordDistance) #0.0352155
mean(bcxy.subs2$CoordDistance) #0.03998489

sd(jacxy.subs2$CoordDistance) #0.02639436
sd(bcxy.subs2$CoordDistance) #0.03214765


distance.db <- bind_rows(bcxy.subs2, jacxy.subs2)
distance.db$Metric <- factor(distance.db$Metric, levels=c('Jaccard', 'Bray-Curtis'), ordered = TRUE)

(dist.coord.box <- ggplot(distance.db, aes(Metric, CoordDistance, fill = Metric)) + 
    stat_boxplot(geom = "errorbar", width = 0.2) + 
    geom_boxplot(notch = TRUE) +
    scale_fill_manual(values = c("#54C9AD", "#40498E"), labels = c('Jaccard', 'Bray-Curtis'), name = "Metric") +
    theme_bw() + 
    ylim(0, 0.2) +
    #ggtitle("(e) Distance per subsite") +
    xlab("\nDissimilarity metric") + ylab("Distance between timepoints\n") +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 24),
          axis.text = element_text(size = 20, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()))

ggsave(dist.coord.box, filename = "figures/Figure_4e.png", 
       width = 16, height = 18, units = "cm")


