## Plant diversity dynamics over space and time in a warming Arctic
## Mariana Garcia Criado (mariana.garcia.criado@gmail.com)
## Script 1. Data checks

# Summary statistics & support figures


## PACKAGES ----
library(tidyverse)
library(cowplot)
library(brms)
library(viridis)
library(ggtern)
library(ggtext)
library(paletteer)
library(randomcoloR)
library(ggpubr)


## FUNCTIONS ----
`%notin%` <- Negate(`%in%`)


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

coords.only <- itex.dec22 %>% distinct(SiteSubsite, .keep_all = T) %>% 
  select(SiteSubsite, LAT, LONG)
write.csv(coords.only, "data/itex_coords_only.csv")

min(itex.dec22$YEAR)

# Get general stats
Xplots <- length(unique(itex.dec22$SiteSubsitePlot))
Xsites <- length(unique(itex.dec22$SITE))
Xspecies <- length(unique(itex.dec22$SPECIES_NAME))
Xsubsites <- length(unique(itex.dec22$SiteSubsite))


## LATITUDE VS ELEVATION ----
(ele.lat <- ggplot(richness.end.plot) + geom_point(aes(x = LAT, y = ELEV)) +  
    geom_smooth(aes(x = LAT, y = ELEV), method = "lm", colour = "black") + bio.theme)

# Elevation and latitude are correlated so elevation won't go into the models.
ele.lat.mod <- brm(ELEV ~ LAT, data = richness.end.plot, iter = 2000, chains = 4, warmup = 400,
              file = "scripts/mgarciacriado/models/lat_elev_mod")
summary(ele.lat.mod) # negative significant





## NUMBER OF PLOTS ----

## How many plots? Including those that have been monitored at least once
nb.plots <- itex.dec22 %>% distinct(SiteSubsitePlot, .keep_all = TRUE) #2,174

# How many plots were surveyed more than once?
rich.time <- itex.dec22 %>% distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  group_by(SiteSubsitePlot) %>% mutate(NumYears = length(unique(YEAR))) %>% 
  ungroup() %>% filter(NumYears > 1) %>% distinct(SiteSubsitePlot, .keep_all = TRUE) # 1387 plots

# How many plots were surveyed more than once and for a minimum of 5 years?
rich.time.5y <- itex.dec22 %>% distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  group_by(SiteSubsitePlot) %>% mutate(duration = max(YEAR) - min(YEAR)) %>% 
  filter(duration > 4) %>% distinct(SiteSubsitePlot, .keep_all = TRUE) #1266 plots

# Which plots were surveyed more than once and but for fewer than 5 years?
rich.time.lessthan5y <- itex.dec22 %>% distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  group_by(SiteSubsitePlot) %>% mutate(NumYears = length(unique(YEAR))) %>% 
  mutate(duration = max(YEAR) - min(YEAR)) %>% 
  filter(NumYears > 1 & duration < 5) %>% distinct(SiteSubsitePlot, .keep_all = TRUE) #121 plots

# Study areas that contain plots monitored more than once but <5 years
unique(rich.time.lessthan5y$SITE)
# "LOGH"       "THUFUVER"   "ZACKENBERG" "ADVENT"     "DISKO"      "SADVENT"    "KANGER"     "IGLOOLIK"   "NARSARSUAQ"
# "ATQASUK"    "BARROW"  



## TIMELINE (ED FIG 3) ----
timeline.db <- itex.dec22 %>% distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>%
  select(SiteSubsitePlotYear, SiteSubsitePlot, SiteSubsite, SITE, YEAR) %>% 
  mutate(SITE = case_when(SITE == "RIRI" ~ "RITSEM",
                          SITE == "LORI" ~ "LANGFJALLET (LORI)",
                          SITE == "LOGH" ~ "LANGFJALLET (LOGH)",
                          SITE == "FURI" ~ "FULUFJALL",
                          SITE == "LATNJA" ~ "LATNJAJAURE",
                          SITE == "QHI" ~ "QIKIQTARUK",
                          SITE == "ADVENT" ~ "ADVENTDALEN",
                          TRUE ~ SITE)) %>% 
  mutate(SiteNoDup = SITE)


# to remove duplicates study areas in y axis
timeline.db$SiteNoDup[duplicated(timeline.db$SiteNoDup)] <- "  " 

advent <- filter(timeline.db, SITE == "ADVENTDALEN")
sadv <- filter(timeline.db, SITE == "SADVENT")

# create colour palette with many classes
n <- 45
palette <- distinctColorPalette(n)
pie(rep(1, n), col=palette)

# plot timeline
(timeline.plot.legend <- ggplot(timeline.db, aes(x = YEAR, y = SiteSubsitePlot, col = SITE)) + 
    geom_line(linewidth = 0.5, alpha = 0.3) + 
    geom_point(aes(YEAR, SiteSubsitePlot, col = SITE), size = 3) +
    scale_x_continuous(expand = c(0,0)) + # remove extra white space 
    scale_y_discrete(breaks = timeline.db$SiteSubsitePlot, labels = timeline.db$SiteNoDup, limits = rev) + # order alphabetically
    xlab("") + ylab("") + 
    labs(colour="Study\narea") +
    scale_color_manual(values = palette) +
    theme(legend.position = "bottom", 
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.key = element_rect(fill = NA), 
          axis.text.x  = element_text(vjust = 0.5, size = 8, colour = "black"), 
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), 
          legend.box.spacing = unit(0, "pt"),
          axis.line = element_line(colour = "black")))

ggsave(timeline.plot.legend, filename = "figures/figure_s3.jpg", 
       width = 18, height = 24, units = "cm")





## SUMMARY FIGURES ----
itex.end.rich <- read.csv("data/end_itex_rich22.csv")

# Compare richness values across methods (ED Fig 1a)
(met <- ggplot(itex.end.rich, aes(x = ValueType, y = CurrentRichness, fill = Method)) + geom_boxplot() + 
    scale_fill_manual(values = c("#35B0CA", "#CA35B0", "#B0CA35")) +
    xlab("\nMethod") + ylab("Plot richness\n") + theme_bw() + 
    theme(axis.title.x = element_text(face="plain", size=8.5),
          axis.text.x  = element_text(vjust=0.5, size=8.5, colour = "black"), 
          axis.title.y = element_text(face="plain", size=8.5),
          axis.text.y  = element_text(vjust=0.5, size=8.5, colour = "black"),
          legend.title = element_text(size=8.5, face="bold"),
          legend.text = element_text(size = 8.5),
          legend.position = "top",
          legend.box.spacing = unit(0, "pt"),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")))

## Method proportion for Appendix conceptual diagram (ED Fig 2)
methods <- itex.dec22 %>% distinct(SiteSubsitePlot, .keep_all = TRUE) %>% 
  group_by(ValueType) %>% summarise(PlotsPerMethod = n()) %>% ungroup() %>%
  mutate(PropPlotsMethod = (PlotsPerMethod / sum(PlotsPerMethod))*100)

# 2174 plots in total
sum(methods$PlotsPerMethod)

## Remove empty rows
plot.size <- itex.dec22 %>% distinct(SiteSubsitePlot, .keep_all = TRUE) %>% 
  filter(!is.na(SurveyedArea)) %>% filter(SurveyedArea != "")

plot.size$SurveyedArea <- as.numeric(plot.size$SurveyedArea)

(plot_size <- ggplot(plot.size, aes(x=SurveyedArea)) + 
    geom_histogram(fill="darkgrey", bins = 70) +
    xlab(expression(paste("\nPlot size (m"^bold("2"), ")"))) + 
    ylab("Count\n") + 
    bio.theme +
    theme(axis.title.y = element_text(face="plain", size=8.5),
          axis.title.x = element_text(face="plain", size=8.5),
          axis.text.x  = element_text(vjust=0.5, size=8.5, colour = "black"),
          axis.text.y  = element_text(vjust=0.5, size=8.5, colour = "black"),
          plot.margin = unit(c(0,0,0,0), "cm")))

# Plot size values 
mean(plot.size$SurveyedArea) # 0.5774178
min(plot.size$SurveyedArea) # 0.0484
max(plot.size$SurveyedArea) # 1

## Plots per site
plots.site <- itex.dec22 %>% distinct(SiteSubsitePlot, .keep_all = TRUE) %>%
  group_by(SITE) %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>% 
  mutate(PlotsPerSite = length(unique(SiteSubsitePlot))) %>% 
  distinct(SITE, .keep_all = TRUE)

# Plot per site values 
mean(plots.site$PlotsPerSite) # 48.31
min(plots.site$PlotsPerSite) # 5
max(plots.site$PlotsPerSite) # 276


## Plots per subsite
plots <- itex.dec22 %>% distinct(SiteSubsitePlot, .keep_all = TRUE) %>%
  group_by(SiteSubsite) %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>%
  mutate(PlotsPerSubsite = length(unique(SiteSubsitePlot))) %>% 
  distinct(SiteSubsite, .keep_all = TRUE)

# histogram version
(plots.plot2 <- ggplot(plots, aes(x=PlotsPerSubsite)) + 
    geom_histogram(fill="darkgrey", bins = 80) + 
    xlab("Plots per subsite") + ylab("Count\n") +
    bio.theme +
    theme(axis.title.y = element_text(face="plain", size=8.5),
          axis.title.x = element_text(face="plain", size=8.5),
          axis.text.x  = element_text(vjust=0.5, size=8.5, colour = "black"),
          axis.text.y  = element_text(vjust=0.5, size=8.5, colour = "black"),
          plot.margin = unit(c(0,0,0,0), "cm")))

# Plot per subsite values 
mean(plots$PlotsPerSubsite) # 14
min(plots$PlotsPerSubsite) # 1
max(plots$PlotsPerSubsite) # 87


## Total surveyed area
plots$SurveyedArea <- as.numeric(plots$SurveyedArea)

surv.area <- plots %>% mutate(TotalSurvArea = PlotsPerSubsite * SurveyedArea) %>% 
  dplyr::select(SITE, SiteSubsite, TotalSurvArea)

# save for #script 6
write.csv(surv.area, "data/total_surveyed_area.csv")

# histogram
(total.area <- ggplot(surv.area, aes(x=TotalSurvArea)) + 
    geom_histogram(fill="darkgrey", binwidth = 1) + 
    xlab(expression(paste("\nTotal surveyed area (m"^bold("2"), ")"))) + 
    ylab("Count\n") + bio.theme +
    theme(
      axis.title.y = element_text(face="plain", size=8.5),
      axis.title.x = element_text(face="plain", size=8.5),
      axis.text.x  = element_text(vjust=0.5, size=8.5, colour = "black"),
      axis.text.y  = element_text(vjust=0.5, size=8.5, colour = "black"),
      plot.margin = unit(c(0,0,0,0), "cm")))


# ED Fig 1b-d: total surveyed area panel
(surv.panel <- plot_grid(plot_size, plots.plot2, total.area, ncol = 3, nrow = 1,
                         align = "hv", labels = c("b", "c", "d"), 
                         label_size = 14))

# ED Fig 1 (full panel)
(fig.s1.panel <- plot_grid(met, surv.panel, nrow = 2, ncol = 1, 
                           labels = c("a", NULL, NULL, NULL), label_size = 14, 
                           rel_heights = c(1.5, 1)))

ggsave(fig.s1.panel, filename = "figures/figure_s1.jpg", 
       width = 18, height = 14, units = "cm")




## Subsites per site
subsites <- itex.dec22 %>% distinct(SiteSubsite, .keep_all = TRUE) %>% group_by(SITE) %>% 
  mutate(SubsitesPerSite = length(unique(SiteSubsite))) %>% distinct(SITE, .keep_all = TRUE)

# Subsites per Site values 
mean(subsites$SubsitesPerSite) # 3.4
min(subsites$SubsitesPerSite) # 1
max(subsites$SubsitesPerSite) # 11


## Number of study areas
studyareas <- itex.dec22 %>% distinct(SITE, .keep_all = TRUE) #45

## Duration of studies
duration <- itex.dec22 %>% group_by(SiteSubsitePlot) %>% 
  mutate(duration = max(YEAR) - min(YEAR)) %>% distinct(SiteSubsitePlot, .keep_all = TRUE)


# Figure 1d
duration0 <- duration %>% filter(duration < 5)
duration5 <- duration %>% filter(duration >= 5)

(duration.hist <- ggplot(duration, aes(duration)) +
    geom_histogram(data = duration0, alpha = 0.5, fill = "grey", binwidth = 1, boundary = 0.5) +
    geom_histogram(data = duration5, alpha = 1, fill = "darkgrey", binwidth = 1, boundary = 0.5) +
    geom_vline(xintercept = 4.5, linetype = "dashed") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    ylab("Count of plots\n") +
    xlab("\nDuration of study (years)") +
    bio.theme +
    theme(axis.title.x = element_text(face = "plain", size = 24), 
          axis.title.y = element_text(face = "plain", size = 24),
          legend.title = element_text(size = 20)))

ggplot2::ggsave(duration.hist, filename = "figures/Figure_1d.png", width = 60, height = 12, units = "cm")


# Duration of study period
mean(duration$duration) # 8.7
min(duration$duration) # 0
max(duration$duration) # 28


## Time points per plot
tp <- itex.dec22 %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(TimePointsPerPlot = length(unique(YEAR))) %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE)

tp.twice <- tp %>% filter(TimePointsPerPlot == 2) #490 plots
tp.three <- tp %>% filter(TimePointsPerPlot == 3) #299 plots (21.5%)
tp.four <- tp %>% filter(TimePointsPerPlot == 4) #274 plots (19.7)
tp.more <- tp %>% filter(TimePointsPerPlot > 1) #1387 plots

#490/1387 = 35.3% of plots surveyed twice

tp.five <- tp %>% filter(TimePointsPerPlot > 4) #324 plots
#324/1387 = 23.3%

tp.ten <- tp %>% filter(TimePointsPerPlot > 9) #7 plots (0.5%)

# Time points of monitoring 
mean(tp$TimePointsPerPlot) # 2.71
min(tp$TimePointsPerPlot) # 1
max(tp$TimePointsPerPlot) # 11


## Time between surveys
time <- itex.dec22 %>% distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  group_by(SiteSubsitePlot) %>% mutate(NbYears = n()) %>% filter(NbYears > 1) %>%
  arrange(YEAR, .by_group = TRUE) %>%
  mutate(diff = YEAR - lag(YEAR, default = first(YEAR))) %>%
  filter(diff > 0) # remove the lines with 0 since we are substracting from previous row and for the first row there is no previous year value

# Time between monitoring surveys
mean(time$diff) # 5.08
min(time$diff) # 1
max(time$diff) # 26


## Are all the plot sizes the same per subsite?
sub.plotsize <- itex.dec22 %>% group_by(SiteSubsite) %>% distinct(SurveyedArea, .keep_all = TRUE)
# ALL of them except for JAMESONLAND:TYSKIT which has 0.1 and 0.21




## PF COMPARISON ----

# Comparing different methods of point-framing in the dataset

## SUMMED
# (this is a raw file and thus not available in this script)
load("scripts/mgarciacriado/data/itex_june2022/pfplot_all.RData")

morpho10.vector0 <- read.csv("data/morpho10_vector.csv")
morpho10.vector <- morpho10.vector0$x

# Create vectors
live <- c("LIVE", "Live", "Alive", NA)
fg <- c("FORB", "SEVER", "SDECI", "SHRUBU", "SHRUB", "GRAMINOIDU", "GRASS", "SEDGE", "RUSH", "GRAMU")
bad.sites <- c("SADVENT:WET_PHOTO", "ALEXFIORD:LEVDOLOMITE", "ALEXFIORD:LEVGRANITE", "SVERDRUP:SVERDRUP", 
               "SADVENT:MES_PHOTO", "ALEXFIORD:LEVDOLOMITE", "LATNJA:MESIC_MEADOW", "LATNJA:CAREX")

# Make necessary filtering
pfplot_all_comp <- pfplot_all %>% 
  unite(SiteSubsitePlotYear, c("SITE", "SUBSITE", "PLOT", "YEAR"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>%
  filter(STATUS %in% live) %>% 
  filter(GFNARROWwalker %in% fg) %>%
  filter(TREATMENT %in% c("CTL", "CONTROL", NA)) %>%
  filter(SiteSubsite %notin% bad.sites) %>%
  filter(ABUNDANCE > 0) %>%
  filter(SiteSubsitePlotYear %notin% morpho10.vector) %>% 
  group_by(SiteSubsitePlot) %>% 
  filter(YEAR == max(YEAR))

# With this filter, no plots at all that have C/U data. So this issue doesn't apply in the original database
sppall <- sort(unique(pfplot_all_comp$SPECIES_NAME))
plotnames <- sort(unique(pfplot_all_comp$SiteSubsitePlot)) # compare with those with massive plot sizes. all NH
unique(pfplot_all_comp$ValueType)


# More relaxed filter just to compare with XY data (not representative of clean database)
cu <- pfplot_all %>%
  unite(SiteSubsitePlotYear, c("SITE", "SUBSITE", "PLOT", "YEAR"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>%
  filter(COVER_UNDERSTORY %in% c("C", "U")) %>% 
  filter(TREATMENT %in% c("CTL", "CONTROL")) %>%
  filter(GFNARROWwalker %in% fg) %>%
  filter(SiteSubsite %notin% bad.sites) %>%
  filter(ABUNDANCE > 0)

# 570 from 45054 = 1.2% records have C/U info
# 16 from 893 plots = 1.79% plots
length(unique(cu$SiteSubsitePlot)) # 16 plots 
length(unique(pfplot_all$PLOT)) # 893 plots 

# the only all plots that have C/U info are pf_topbot_plot in Lantja
# so we can only compare top-bottom to top and bottom, no middle hits to check

# Calculate richness
cu.rich <- cu %>% group_by(SiteSubsitePlot) %>% 
  filter(YEAR == max(YEAR))

# top hits only
cu.top <- cu.rich %>% filter(COVER_UNDERSTORY == "C") %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(Richness = length(unique(SPECIES_NAME))) %>% ungroup() %>%  
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>% 
  select(SiteSubsitePlot, Richness, COVER_UNDERSTORY) %>%
  mutate(COVER_UNDERSTORY = "Top hit only")

# bottom hits only majority are richness = 1
cu.bot <- cu.rich %>% filter(COVER_UNDERSTORY == "U") %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(Richness = length(unique(SPECIES_NAME))) %>% ungroup() %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>%
  select(SiteSubsitePlot, Richness, COVER_UNDERSTORY) %>%
  mutate(COVER_UNDERSTORY = "Bottom hit only")

# all hits
cu.both <- cu.rich %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(Richness = length(unique(SPECIES_NAME))) %>% ungroup() %>%  
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>%
  select(SiteSubsitePlot, Richness, COVER_UNDERSTORY) %>%
  mutate(COVER_UNDERSTORY = "Top-bottom hits")

# bind in dataframe
cu.all <- rbind(cu.top, cu.bot, cu.both)

# Show in graph
(cu.plot <- ggplot(cu.all, aes(x = COVER_UNDERSTORY, y = Richness)) + geom_boxplot() + 
    xlab("\nPoint-framing (sum)") + ylab("Plot Richness") + theme_bw())




## XY data
# (this is a raw file and thus not available in this script)
load("scripts/mgarciacriado/data/itex_june2022/pfxy_all.RData")

# If we apply filters including live specimens again this doesn't yield any results, so this no longer applies
pfxy_allhits <- pfxy_all %>% 
  unite(SiteSubsitePlotYear, c("SITE", "SUBSITE", "PLOT", "YEAR"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>%
  filter(ValueType %in% c("pf_all_XY", "pf_all_xy")) %>%
  filter(HIT %in% c("C", "U", "M", "bottom", "top", "middle1", "middle2", "middle3")) %>%
  filter(TREATMENT %in% c("CTL", "CONTROL")) %>%
  filter(GFNARROWwalker %in% fg) %>%
  filter(STATUS %in% live) %>% 
  filter(SiteSubsite %notin% bad.sites) %>%
  filter(ABUNDANCE > 0)


# Calculate richness
pfxy_allhits_rich <- pfxy_allhits %>% group_by(SiteSubsitePlot) %>% 
  filter(YEAR == max(YEAR))

# 195109 from 1563331 = 12.4% XY records have C/U info

# 308/999 = 30.83% plots have C/U info
unique(pfxy_all$PLOT) # 999 plots
unique(pfxy_allhits$SiteSubsitePlot) # 308 plots



# top hits only
xy.top <- pfxy_allhits_rich %>% filter(HIT %in% c("C", "top")) %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(Richness = length(unique(SPECIES_NAME))) %>% ungroup() %>%  
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>% 
  select(SiteSubsitePlot, Richness, HIT) %>%
  mutate(HIT = "Top hit only")

# bottom hits only 
xy.bot <- pfxy_allhits_rich %>% filter(HIT %in% c("U", "bottom")) %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(Richness = length(unique(SPECIES_NAME))) %>% ungroup() %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>%
  select(SiteSubsitePlot, Richness, HIT) %>%
  mutate(HIT = "Bottom hit only")

# top-bottom hits
xy.topbot <- pfxy_allhits_rich %>% filter(HIT %in% c("U", "bottom", "C", "top")) %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(Richness = length(unique(SPECIES_NAME))) %>% ungroup() %>%  
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>%
  select(SiteSubsitePlot, Richness, HIT) %>%
  mutate(HIT = "Top-bottom hits")

# all hits
xy.all <- pfxy_allhits_rich %>% 
  filter(HIT %in% c("U", "bottom", "C", "top", "M", "middle1", "middle2", "middle3")) %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(Richness = length(unique(SPECIES_NAME))) %>% ungroup() %>%  
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>%
  select(SiteSubsitePlot, Richness, HIT) %>%
  mutate(HIT = "All hits")

# bind in dataframe
xy.all.df <- rbind(xy.top, xy.bot, xy.topbot, xy.all)

# Show in graph
(xy.plot <- ggplot(xy.all.df, aes(x = HIT, y = Richness)) + geom_boxplot() + 
    xlab("\nPoint-framing (XY)") + ylab("Plot Richness") + theme_bw())






## FG CORRELATIONS ----

# Keep one line with all covers values per functional group 
dom.fg2 <- itex.dec22 %>% distinct(SiteSubsitePlot, .keep_all = TRUE)

# Shrub vs Graminoid
shrub.gram.mod <- brm(PlotGraminoidMean ~ PlotShrubMean, 
                      data = dom.fg2, iter = 2000, chains = 4, warmup = 400, 
                      file = "smodels/shrub_gram_mod")
summary(shrub.gram.mod) # negative significant: more shrubs, less grams. More grams, less shrubs

# Shrub vs Forb
shrub.forb.mod <- brm(PlotForbMean ~ PlotShrubMean, 
                      data = dom.fg2, iter = 2000, chains = 4, warmup = 400, 
                      file = "models/shrub_forb_mod")
summary(shrub.forb.mod) # negative significant: more shrubs, less forbs. More forbs, less shrubs.

# Forb vs Graminoid
gram.forb.mod <- brm(PlotGraminoidMean ~ PlotForbMean, 
                     data = dom.fg2, iter = 2000, chains = 4, warmup = 400, 
                     file = "models/gram_forb_mod")
summary(gram.forb.mod) # negative significant: more grams, less forbs - more forbs, less grams.


## TERNARY PLOT (FIG 1C)
(tern.plot <- ggtern(data = dom.fg2, aes(x = PlotForbMean, y = PlotGraminoidMean, z = PlotShrubMean)) +
    geom_point(size = 3, alpha = 0.5, aes(colour = MeanRichness)) +
    labs(x="Forb Cover", y="Graminoid Cover", z="Shrub Cover", colour = "Plot Richness") +
    scale_colour_viridis(option = "plasma", begin = 0, end = 0.95) + theme_bw() + 
   theme(legend.title = element_text(size = 26), legend.text=element_text(size = 20), 
         tern.axis.text.T = element_text(size =22), tern.axis.text.L = element_text(size =22), 
         tern.axis.text.R = element_text(size =22), tern.axis.title.T = element_text(size =28), 
         tern.axis.title.L = element_text(size =28), tern.axis.title.R = element_text(size =28)))

ggsave(tern.plot, filename = "figures/Figure_1c.png", 
       width = 25, height = 25, units = "cm")

# mean plot cover per functional group
shrub.mean <- dom.fg2 %>% summarise(mean(PlotShrubMean)) # 50%
gram.mean <- dom.fg2 %>% summarise(mean(PlotGraminoidMean)) # 37.4%
forb.mean <- dom.fg2 %>% summarise(mean(PlotForbMean)) # 12.6%


## Histograms
(shrub.hist <- ggplot(dom.fg2, aes(x=PlotShrubMean)) + 
    geom_histogram(binwidth = 5, fill = "#5fbf33") + 
    geom_vline(aes(xintercept = mean(PlotShrubMean)), colour = "black", linetype = "dotted", size = 1) +
    xlab("\nShrub Cover (%)") + ylab("Number of plots\n") + bio.theme +
    theme(axis.title.x = element_text(face="plain", size=9),
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y = element_text(vjust=0.5, face="plain", size=9),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), units = , "cm")))

(gram.hist <- ggplot(dom.fg2, aes(x=PlotGraminoidMean)) + 
    geom_histogram(binwidth = 5, fill = "#dc9537") + 
    geom_vline(aes(xintercept = mean(PlotGraminoidMean)), colour = "black", linetype = "dotted", size = 1) +
    xlab("\nGraminoid Cover (%)") + ylab("Number of plots\n") + bio.theme+
    theme(axis.title.x = element_text(face="plain", size=9),
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y = element_text(vjust=0.5, face="plain", size=9),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), units = , "cm")))

(forb.hist <- ggplot(dom.fg2, aes(x=PlotForbMean)) + 
    geom_histogram(binwidth = 5, fill = "#95558b") + 
    geom_vline(aes(xintercept = mean(PlotForbMean)), colour = "black", linetype = "dotted", size = 1) +
    xlab("\nForb Cover (%)") + ylab("Number of plots\n") + bio.theme  +
    theme(axis.title.x = element_text(face="plain", size=9),
          axis.text.x  = element_text(vjust=0.5, size=9, colour = "black"), 
          axis.title.y = element_text(face="plain", size=9),
          axis.text.y = element_text(vjust=0.5, face="plain", size=9),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), units = , "cm")))


# ES Fig 5-ac: 
(S5A <- plot_grid(shrub.hist, gram.hist, forb.hist, ncol = 3, nrow = 1,   
                 align = "h", labels = c("a", "b", "c"), label_size = 14))



## How many plots (%) are dominated by each functional group?
domfg.plots <- dom.fg2 %>% group_by(DominatingFG) %>% 
  summarise(total = n()) %>% 
  mutate(percent = (total/sum(total)*100))

## Model to see if increased richness is related to functional groups 
rich.forb.mod <- brm(PlotForbMean ~ MeanRichness, data = dom.fg2, iter = 2000, chains = 4, warmup = 400, 
                     file = "scripts/mgarciacriado/models/forb_rich_mod")
summary(rich.forb.mod) # positive significant


rich.gram.mod <- brm(PlotGraminoidMean ~ MeanRichness, data = dom.fg2, iter = 2000, chains = 4, warmup = 400, 
                     file = "scripts/mgarciacriado/models/gram_rich_mod")
summary(rich.gram.mod) # positive significant


rich.shb.mod <- brm(PlotShrubMean ~ MeanRichness, data = dom.fg2, iter = 2000, chains = 4, warmup = 400, 
                     file = "scripts/mgarciacriado/models/shb_rich_mod")
summary(rich.shb.mod) # negative significant

