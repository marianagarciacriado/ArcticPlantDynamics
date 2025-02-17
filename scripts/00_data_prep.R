## Plant diversity dynamics over space and time in a warming Arctic
## Mariana Garcia Criado (mariana.garcia.criado@gmail.com)
## Script 0. ITEX data cleaning 

# This script is provided to show the process of generating the main input file of plant composition (itex_dec22.RData), 
# but the original raw files to produce this file are not provided here. 
# The original raw files from the ITEX+ database will be provided when the data paper led by Dr Anne Bjorkman
# is published.


## PACKAGES ----
library(tidyverse)
library(brms)
library(Taxonstand)


## FUNCTIONS ----
`%notin%` <- Negate(`%in%`)


## DATA ----
load("scripts/mgarciacriado/data/itex_june2022/pfxy_all.RData")
load("scripts/mgarciacriado/data/itex_june2022/pfplot_all.RData")
load("scripts/mgarciacriado/data/itex_june2022/perccov_all.RData")
qhi <- read.csv("scripts/mgarciacriado/data/qhi_nov2022/qhi-1999-2022-clean-nov22.csv") #43064

unique(pfxy_all$SITE)
unique(pfplot_all$SITE)
unique(perccov_all$SITE)

# Check variable structure
str(pfxy_all)
str(qhi)

# Consistent QHI variable structure
qhi.good <- qhi %>% select(-1)
qhi.good$PLOT <- as.character(qhi.good$PLOT)
qhi.good$X <- as.character(qhi.good$X)
qhi.good$Y <- as.character(qhi.good$Y)
qhi.good$HIT <- as.character(qhi.good$HIT)
qhi.good$ABUNDANCE <- as.numeric(qhi.good$ABUNDANCE)

# Removing QHI from ITEX dataset and include the most recent cleaned one
pfxy_all.qhi <- pfxy_all %>% filter(SITE != "QHI") %>% bind_rows(., qhi.good) #1572755


# Raw figures of plots and species number
length(unique(pfxy_all$PLOT)) #999
length(unique(pfplot_all$PLOT)) #893
length(unique(perccov_all$PLOT)) #6757
length(unique(qhi.good$PLOT)) #6
#total = 8655


length(unique(pfxy_all$SPECIES_NAME)) #1033
length(unique(pfplot_all$SPECIES_NAME)) #615
length(unique(perccov_all$SPECIES_NAME)) #564
length(unique(qhi.good$SPECIES_NAME)) #82
# total = 2294


## STRUCTURE FIXES ----

# First of all fix Latnjajaure's name to Latnja to standardise
pfplot_all0 <- pfplot_all %>% mutate(SITE = case_when(SITE == "LATNJAJAURE" ~ "LATNJA", TRUE ~ SITE))


# Now fix plot names that include year - important because later on we will group by plot, not by plotXyear
pfplot_allX <- pfplot_all0 %>% 
  mutate(YearInPlot = ifelse(str_detect(PLOT, fixed(as.character(YEAR))), "YES", "NO")) %>%
  mutate(PLOT = ifelse(YearInPlot == "YES", str_remove(PLOT, pattern = as.character(YEAR)), PLOT)) %>%
  mutate(PLOT = str_remove(PLOT, pattern = "^_+"), # Removes leading "_"
         PLOT = str_remove(PLOT, pattern = "_$")) %>% # Removes trailing "_"  
  select(., -YearInPlot) 


# This is not the case in xy dataframe
# pfxy_all.qhiX <- pfxy_all.qhi %>% mutate(YearInPlot = ifelse(str_detect(PLOT, fixed(as.character(YEAR))), "YES", "NO"))

# cover dataframe
perccov_allX <- perccov_all %>% 
  mutate(YearInPlot = ifelse(str_detect(PLOT, fixed(as.character(YEAR))), "YES", "NO")) %>%
  mutate(PLOT = ifelse(YearInPlot == "YES", str_remove(PLOT, pattern = as.character(YEAR)), PLOT)) %>%
  mutate(PLOT = str_remove(PLOT, pattern = "^_+"), # Removes leading "_"
         PLOT = str_remove(PLOT, pattern = "_$")) %>% # Removes trailing "_"  
  select(., -YearInPlot) 

unique(perccov_allX$SITE)

# Not permanently marked plots are identified in the manner "651a, 651b and 651c" (as plot names) 
# where the letter indicates a different sampling year (also identified in the year column).
# I am only keeping the last year of each plot to include in the spatial analyses, but won't keep the previous years because
# it will definitely include uncertainty, increased turnover and local extinctions and colonisations that might be artifacts.

# These are all plots from the same year so no cleaning needed
# pfplot_allXX <- pfplot_allX %>% 
#  mutate(LetterInPlot = ifelse(str_detect(PLOT, "[abc]$"), "YES", "NO")) 

# These are named with some letters but it's always the same plot name, and the year is specified correctly
# pfxy_all.qhiXX <- pfxy_all.qhi %>% 
#  mutate(LetterInPlot = ifelse(str_detect(PLOT, "[abc]$"), "YES", "NO"))

# These are all plots where a = start year and b, c = subsequent years
perccov_allXX <- perccov_allX %>% 
  mutate(LetterInPlot = ifelse(str_detect(PLOT, "[abc]$"), "YES", "NO")) %>%
  mutate(PLOT = ifelse(LetterInPlot == "YES", str_remove(PLOT, "[abc]$"), PLOT)) 
  
  # database without the letter plots
perccov_noletter <- perccov_allXX %>% filter(LetterInPlot == "NO")
  
  # database with the letter plots, filter for only the last year
perccov_letter <- perccov_allXX %>% filter(LetterInPlot == "YES") %>% group_by(PLOT) %>% 
  filter(YEAR == max(YEAR)) %>% ungroup()

  # merge together
perccov_fixed <- bind_rows(perccov_noletter, perccov_letter) %>% select(., -LetterInPlot)


# Sites from where the early years were removed: 
# AKUREYRI, BLONDUOS, DALSMYNNI, HJARDARLAND, HOLTAVORDUHEIDI, MODRUVELLIR, OXNADALSHEIDI and THYKKVIBAER

# covletter <- perccov_allX %>% mutate(LetterInPlot = ifelse(str_detect(PLOT, "[abc]$"), "YES", "NO")) %>%
#  filter(LetterInPlot == "YES")

unique(perccov_fixed$SITE)



## FILTERING ----

# Keep control plots only
unique(perccov_fixed$TREATMENT) # "CTL"     "CTLNG"   "OTC"     "CONTROL"
unique(pfplot_allX$TREATMENT) # "OTC"     "CTL"     "OTCW"    "OTCWS"   "CONTROL" "WARMING" "DAMAGE" 
unique(pfxy_all.qhi$TREATMENT) # "OTC"      "CTL"      "CTLNG"    "CONTROL"  "WARMING"  NA         "CLT_GRA"  "CTL_NGRA"

trtmt.vector <- c("CTL", "CONTROL", "CLT_GRA", NA)

perccov_all2 <- perccov_fixed %>% filter(TREATMENT %in% trtmt.vector)
pfxy_all2 <- pfxy_all.qhi %>% filter(TREATMENT %in% trtmt.vector)
pfplot_all2 <-pfplot_allX  %>% filter(TREATMENT %in% trtmt.vector)


# only keeping alive hits and NA
unique(perccov_all2$STATUS) # "LIVE"         "LITTER"       "OTHER"        "DEAD"         "UNK"          "STANDINGDEAD" NA  
unique(pfplot_all2$STATUS) # "LIVE"   "LITTER"  "OTHER"  "DEAD"    "STANDINGDEAD" NA       "UNK"  "SOIL"   "ROCK"         
unique(pfxy_all2$STATUS) # "LIVE" "OTHER" "LITTER" "DEAD" "STANDINGDEAD"  "UNK" "STANDING DEAD" NA "ROCK_SOIL" "Alive" "Dead"  
# "UNKNOWN"       "SOIL"          "ROCK"      "Standing dead" "Live"          "Litter"      "N/A"

# Check the empty space - it's 2 records from QHI so we assume they are the same as NA
# emptyspace <- pfxy_all2 %>% filter(STATUS == "")

status.vector <- c("LIVE", NA, "Alive", "Live", "N/A", "")

perccov_all3 <- perccov_all2 %>% filter(STATUS %in% status.vector)
pfxy_all3 <- pfxy_all2 %>% filter(STATUS %in% status.vector)
pfplot_all3 <- pfplot_all2 %>% filter(STATUS %in% status.vector)


# Functional groups
unique(perccov_all3$GFNARROWwalker)
unique(pfplot_all3$GFNARROWwalker)
unique(pfxy_all3$GFNARROWwalker)

# Check FG = unknown
unk.cov <- filter(perccov_all3, GFNARROWwalker == "UNK") #	Cerastium arcticum
unique(unk.cov$SPECIES_NAME)
#[1] "XXXOTHERSP"         "XXXVASPLA"          "XXXVASC"            "XXXNONVASC"        
# "XXXUNKNOWN"  "XXXGRAEQU"  "Cerastium arcticum"

unk.pfxy <- filter(pfxy_all3, GFNARROWwalker %in% c("UNK", "UNKNOWN"))
unique(unk.pfxy$SPECIES_NAME)
#"XXXUNK"              "XXXERIST"            "XXXUNKPLA"              "XXXUNKNOWN"   " "        

unk.pfall <- filter(pfplot_all3, GFNARROWwalker == "UNK") # "XXXFORB"
unique(unk.pfall$SPECIES_NAME)
# XXXUNKPLA" " NA"         "XXXFORB" 


# Check FG = NA
na.cov <- perccov_all3 %>% filter(is.na(GFNARROWwalker))
unique(na.cov$SPECIES_NAME) # "Carex maritima"    "Bupleurum americanum" "XXXCRUST"    

na.xy <- pfxy_all3 %>% filter(is.na(GFNARROWwalker))
unique(na.xy$SPECIES_NAME) # NA             "XXXPOHSPP"    "XXXEquisetum" "XXXALGAE"     

na.all <- pfplot_all3 %>% filter(is.na(GFNARROWwalker)) #nothing


# Add in functional group to the clear vascular plants
perccov_all4 <- perccov_all3 %>% 
  mutate(GFNARROWwalker = case_when(GFNARROWwalker == "UNK" & SPECIES_NAME == "Cerastium arcticum" ~ "FORB", 
                                    is.na(GFNARROWwalker) & SPECIES_NAME == "Carex maritima" ~ "SEDGE",
                                    is.na(GFNARROWwalker) & SPECIES_NAME == "Bupleurum americanum" ~ "FORB", TRUE ~ GFNARROWwalker))

pfplot_all4 <- pfplot_all3 %>% 
  mutate(GFNARROWwalker = case_when(GFNARROWwalker == "UNK" & SPECIES_NAME == "XXXFORB" ~ "FORB", 
                                    TRUE ~ GFNARROWwalker))

fg.vector <- c("FORB", "SEVER", "SDECI", "SHRUBU", "SHRUB", "GRAMINOIDU", "GRASS", "SEDGE", "RUSH", 
               "GRAMU", "WOODYU")

# Note: did NOT retain NA (as they're fixed below), and "XXXNA" as status = other and will be removed regardless because of dead/experimental.

perccov_all5 <- perccov_all4 %>% filter(GFNARROWwalker %in% fg.vector)
pfplot_all5 <- pfplot_all4 %>% filter(GFNARROWwalker %in% fg.vector)
pfxy_all4 <- pfxy_all3 %>% filter(GFNARROWwalker %in% fg.vector)

# Create unique plotXyear ID
perccov_all6 <- perccov_all5 %>% 
  unite(SiteSubsitePlotYear, c("SITE", "SUBSITE", "PLOT", "YEAR"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>%
  tidyr::replace_na(list(ABUNDANCE = 0)) %>% # STEPSTONES:TUNDRA1 and STEPSTONES:TUNDRA2 which are NA should be 0
  filter(SiteSubsite != "SADVENT:WET_PHOTO") # incomplete subsite, removing it to avoid issues in the cover calc below

pfplot_all6 <- pfplot_all5 %>% 
  unite(SiteSubsitePlotYear, c("SITE", "SUBSITE", "PLOT", "YEAR"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>% 
  select(., -COVER_UNDERSTORY) %>% tidyr::replace_na(list(ABUNDANCE = 0)) # KANGERS should be 0 (investigated below)

# Investigate those with Abundance = NA
kanger.na <- pfplot_all5 %>% filter(SUBSITE %in% c("BASHFUL", "DOPEY", "SNEEZY"))
# When it's a 1 the info is in there, but there are no 0s. It's always the same species over the years so these should be Abundance = 0

pfxy_all5 <- pfxy_all4 %>% 
  unite(SiteSubsitePlotYear, c("SITE", "SUBSITE", "PLOT", "YEAR"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>%
  tidyr::replace_na(list(ABUNDANCE = 1)) #1 NA abundance value for QHI:HE only, should be 1 instead as no 0s are recorded


# "SUM" and "sum" don't belong there
unique(pfxy_all5$HIT) 

# "SUM" is Gavia which will be removed from the database (<60deg latitude) so I'm dropping it here
# "sum" is Zackenberg so moving it to pfplot_all
summed <- pfxy_all5 %>% filter(HIT == "sum") %>% select(., -c(X, Y, HIT)) #66 obs

# Remove from XY database
pfxy_all6 <- pfxy_all5 %>% filter(HIT != "sum" | is.na(HIT))

# add in to sum database
pfplot_all7 <- rbind(pfplot_all6, summed)




## COVER CONVERSION ----

# Do the cover values add up to 100?
cov <- perccov_all6 %>% filter(ValueType == "percent_cover") %>% 
  group_by(SiteSubsitePlotYear) %>% summarise(sum = sum(ABUNDANCE))
# Quite a lot of values over 100 so they need to be made proportional too so all values are comparable


#### Cover-equivalent ####

# Confirm that 1 row = 1 species
cov.test <- perccov_all6 %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(NumberRows = n()) %>% mutate(NumberSpecies = length(unique(SPECIES_NAME))) %>%
  mutate(SameOrNot = ifelse(NumberRows == NumberSpecies, "Same", "Different")) %>% ungroup()

# Check if this is because of the missing species names or actually there are repeated species names
dif <- cov.test %>% filter(SameOrNot == "Different")

# Add up values per species so we end up with only one row per species
dif2 <- dif %>% 
  filter(SiteSubsitePlotYear != "BARROW:CAREX_MOIST_MEADOW_MICROTOPO:BC02.5:1999") %>%
  group_by(SiteSubsitePlotYear, SPECIES_NAME) %>% mutate(AbundanceFixed = sum(ABUNDANCE)) %>% ungroup() %>%
  group_by(SiteSubsitePlotYear) %>% distinct(SPECIES_NAME, .keep_all = TRUE) %>% ungroup()

# Confirm that this has worked and 1 row = 1 species
dif3 <- dif2 %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(NumberRows = n()) %>% mutate(NumberSpecies = length(unique(SPECIES_NAME))) %>%
  mutate(SameOrNot = ifelse(NumberRows == NumberSpecies, "Same", "Different")) %>% ungroup()

# One site has duplicate values: all records have exactly the same abundance values twice
barrow.dup <- dif %>% filter(SiteSubsitePlotYear == "BARROW:CAREX_MOIST_MEADOW_MICROTOPO:BC02.5:1999") %>% 
  distinct(SPECIES_NAME, .keep_all = TRUE)


# 1) Remove inconsistent plots from original dataset
perccov_all7 <- cov.test %>% filter(SameOrNot != "Different")

# 2) Dataframe with summed up values
sum.fixed <- dif2 %>% mutate(ABUNDANCE = AbundanceFixed) %>% select(., -AbundanceFixed)

# 3) Barrow with no duplicates (barrow.dup)

# Bind all three into one fixed cover dataset
perccov_fixed0 <- rbind(perccov_all7, sum.fixed, barrow.dup)

# Keep relevant columns only
perccov_fixed <- perccov_fixed0 %>% filter(ABUNDANCE > 0) %>% select(., -c(NumberRows, NumberSpecies, SameOrNot))



# Convert all values to relative cover
itex.cov <- perccov_fixed %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalAbundance = sum(ABUNDANCE)) %>%
  mutate(RelCover = (ABUNDANCE/TotalAbundance)*100) %>% ungroup() # 14041 obs

# Confirm that total cover values add up to 100 in every plotXyear
cov.check <- itex.cov %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalCover = sum(RelCover)) %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE)




#### Point-framing (summed) ####

# Confirm that 1 row = 1 species
pfsum.test <- pfplot_all7 %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(NumberRows = n()) %>% mutate(NumberSpecies = length(unique(SPECIES_NAME))) %>%
  mutate(SameOrNot = ifelse(NumberRows == NumberSpecies, "Same", "Different")) %>% ungroup()

# Check if this is because of the missing species names or actually there are repeated species names
dif.sum <- pfsum.test %>% filter(SameOrNot == "Different")

# Duplicates: exactly the same species and values on repeat
ab.vector <- c("ABISKO:PEATLAND:AA1:2000", "ABISKO:PEATLAND:AA1:2002", "ABISKO:PEATLAND:AA1:2004", "ABISKO:PEATLAND:AA1:2006", "ABISKO:PEATLAND:AA1:2008",
               "ABISKO:PEATLAND:AA2:2000", "ABISKO:PEATLAND:AA2:2002", "ABISKO:PEATLAND:AA2:2004", "ABISKO:PEATLAND:AA2:2006", "ABISKO:PEATLAND:AA2:2008",
               "ABISKO:PEATLAND:AA3:2000", "ABISKO:PEATLAND:AA3:2002", "ABISKO:PEATLAND:AA3:2004", "ABISKO:PEATLAND:AA3:2006", "ABISKO:PEATLAND:AA3:2008",
               "ABISKO:PEATLAND:AA4:2000", "ABISKO:PEATLAND:AA4:2002", "ABISKO:PEATLAND:AA4:2004", "ABISKO:PEATLAND:AA4:2006", "ABISKO:PEATLAND:AA4:2008",
               "ABISKO:PEATLAND:AA5:2000", "ABISKO:PEATLAND:AA5:2002", "ABISKO:PEATLAND:AA5:2004", "ABISKO:PEATLAND:AA5:2006", "ABISKO:PEATLAND:AA5:2008")

# Multiple plots have duplicate values: all records have exactly the same values twice
ab.dup <- dif.sum %>% filter(SiteSubsitePlotYear %in% ab.vector) %>% group_by(SiteSubsitePlotYear) %>%
  distinct(SPECIES_NAME, .keep_all = TRUE) %>% ungroup()

# Add up values per species so we end up with only one row per species
dif.sum2 <- dif.sum %>% filter(SiteSubsitePlotYear %notin% ab.vector) %>% 
  group_by(SiteSubsitePlotYear, SPECIES_NAME) %>% mutate(AbundanceFixed = sum(ABUNDANCE)) %>% ungroup() %>%
  group_by(SiteSubsitePlotYear) %>% distinct(SPECIES_NAME, .keep_all = TRUE) %>% ungroup()

# Confirm that this has worked and 1 row = 1 species
dif.sum.test <- dif.sum2 %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(NumberRows = n()) %>% mutate(NumberSpecies = length(unique(SPECIES_NAME))) %>%
  mutate(SameOrNot = ifelse(NumberRows == NumberSpecies, "Same", "Different")) %>% ungroup()


# Dataframes to merge:

# 1) Remove inconsistent plots from original dataset
pfplot_all8 <- pfsum.test %>% filter(SameOrNot != "Different")

# 2) Dataframe with added values
dif.sum3 <- dif.sum2 %>% mutate(ABUNDANCE = AbundanceFixed) %>% select(., -AbundanceFixed)

# 3) Abisko with no duplicates (ab.dup)

# Bind all three into one fixed cover dataset
pfsum_fixed0 <- rbind(pfplot_all8, dif.sum3, ab.dup)

# Keep relevant columns only
pfsum_fixed <- pfsum_fixed0 %>% filter(ABUNDANCE > 0) %>% select(., -c(NumberRows, NumberSpecies, SameOrNot))




# Convert all values to relative cover
itex.pfsum <- pfsum_fixed %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalAbundance = sum(ABUNDANCE)) %>%
  mutate(RelCover = (ABUNDANCE/TotalAbundance)*100) %>% ungroup() # 15427 obs

# Confirm that total cover values add up to 100 in every plotXyear
pfsum.check <- itex.pfsum %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalCover = sum(RelCover)) %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE)





#### Point-framing (XY) ####

# There are XY data that only have one coordinate - I'm filling in the NA cell so we avoid problems with cover calculation
na.x <- pfxy_all6 %>% filter(is.na(X)) # all filled in
na.y <- pfxy_all6 %>% filter(is.na(Y)) # 139256 obs

# Replace Y coords that are NA by 0s, create a unique coordinate
pfxy_all7 <- pfxy_all6 %>% tidyr::replace_na(list(Y = "0")) %>% 
  unite(XY, c("X", "Y"), sep = "_", remove = FALSE)

# No point on checking that 1 row = 1 species because multiple entries of the same species per plotXyear


# STEP 1: Convert species abundance to presence/absence 
# (2D, not considering multiple hits of the same species at each xy coord, just 1)
pfxy_all_pa <- pfxy_all7 %>% mutate(ABUNDANCE = ifelse(ABUNDANCE > 1, 1, ABUNDANCE)) %>%
  filter(ABUNDANCE > 0) %>%
  group_by(SiteSubsitePlotYear, XY) %>% distinct(SPECIES_NAME, .keep_all = TRUE) %>% ungroup()
  

# STEP 2: Calculate unique species hits per plot and total unique species hits per plot
pfxy_all_pa2 <- pfxy_all_pa %>% group_by(SiteSubsitePlotYear, SPECIES_NAME) %>% 
  mutate(UniqueSpHitsPlot = n()) %>% 
  distinct(SiteSubsitePlotYear, SPECIES_NAME, .keep_all = TRUE) %>% 
  ungroup() %>% select(., -c(X, Y, XY, HIT)) %>% group_by(SiteSubsitePlotYear) %>%
  mutate(TotalUniqueSpHitsPlot = sum(UniqueSpHitsPlot)) %>% ungroup()

# STEP 3: Calculate cover per species
pfxy_all_cov <- pfxy_all_pa2 %>% mutate(RelCover = (UniqueSpHitsPlot/TotalUniqueSpHitsPlot)*100) #58542


# Confirm that total cover values add up to 100 in every plotXyear
pfxy.check <- pfxy_all_cov %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalCover = sum(RelCover)) %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE)



## BINDING ----

# Keep same number of relevant columns
itex.cov.f <- itex.cov %>% select(., -c(TotalAbundance, ABUNDANCE))
itex.pfsum.f <- itex.pfsum %>% select(., -c(TotalAbundance, ABUNDANCE))
pfxy_all_cov.f <- pfxy_all_cov %>% select(., -c(ABUNDANCE, UniqueSpHitsPlot, TotalUniqueSpHitsPlot))

# Bind all methods in one database
itex.all <- rbind(itex.cov.f, itex.pfsum.f, pfxy_all_cov.f) #88010



## SITE CHECKS ----

# Plots that had inconsistent surveyed areas/methods or didn't identify species in the whole subsite
out.subsites <- c("SADVENT:WET_PHOTO", "SADVENT:MES_PHOTO",
                  "ALEXFIORD:LEVDOLOMITE", "ALEXFIORD:LEVGRANITE", "SVERDRUP:SVERDRUP",
                  "LATNJA:MESIC_MEADOW", "LATNJA:CAREX")

# Remove inconsistent sites, Tibet and Mongolia out because they are more meadow than tundra
itex.all2 <- itex.all %>% filter(SiteSubsite %notin% out.subsites) %>% filter(SITE != "TIBET") # 86964



## METADATA ----

# Updated metadata file November 2022
metadata0 <- read.csv("scripts/mgarciacriado/data/itex_june2022/TVC_SITE_SUBSITE_UPDATED2022.csv")

# Keep only relevant columns
metadata <- metadata0 %>% unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>% 
  select(SiteSubsite, MOISTURE, LAT, LONG, ELEV, AZONE, SurveyedArea)

# Merge with composition data
itex.full <- left_join(itex.all2, metadata, by = "SiteSubsite")

# Keep Arctic and small plots only 
itex.full2 <- itex.full %>% filter(LAT > 60| is.na(LAT)) %>% 
  filter(SurveyedArea < 1.1 | is.na(SurveyedArea)) # 47943


# Check for missing NA
unique(itex.full2$MOISTURE)
unique(itex.full2$SurveyedArea)
unique(itex.full2$LAT)
unique(itex.full2$LONG)

# Which subsites didn't get metadata?
nomoist <- itex.full2 %>% filter(MOISTURE == "" | is.na(MOISTURE))
sort(unique(nomoist$SiteSubsite)) #29 subsites

nolat <- itex.full2 %>% filter(is.na(LAT))
sort(unique(nolat$SiteSubsite)) #3 subsites without latitude data

nolong <- itex.full2 %>% filter(is.na(LONG))
sort(unique(nolat$SiteSubsite)) #3 subsites without longitude data

nosize <- itex.full2 %>% filter(SurveyedArea == "" |is.na(SurveyedArea))
missing.size <- sort(unique(nosize$SiteSubsite)) # 1 subsite


# Fill in missing surveyed area - thanks Joe for Jameson code!
itex.full3 <- itex.full2 %>% 
  mutate(SurveyedArea = ifelse(SiteSubsite == "AUYUITTUQ:OWL RIVER", 1, SurveyedArea), 
         SurveyedArea = ifelse(SiteSubsite %in% c("JAMESONLAND:TYSKIT"), 0.2178, SurveyedArea),
         SurveyedArea = ifelse(SiteSubsite %in% c("JAMESONLAND:TYSKIT") & str_detect(PLOT, paste(c("13", "CASEMP", "BBN"),collapse = '|')), 0.1089, SurveyedArea))


# Identify sites known to have LAT and LONG that fall in the ocean (identified during climate extraction)
water.LAT.LONG <- itex.full3 %>% 
  dplyr::select(SiteSubsite, YEAR, LAT, LONG) %>% 
  filter(LAT %in% c(74.28, 74.29)) %>% 
  distinct(SiteSubsite, .keep_all = TRUE)

water.LAT.LONG <- unique(water.LAT.LONG$SiteSubsite)

# Incorrect LAT and LONG: 8 x "ZACKENGERG" SITES

# Manually input new LAT and LONG info for those missing or in water (all above 60 so not removed for being 'non-Arctic')
itex.full4 <- itex.full3 %>% 
  mutate(LAT = ifelse(SiteSubsite == "DISKO:DRYHEATH_FLUX", 69.27, LAT),
         LAT = ifelse(SiteSubsite == "DISKO:WETFEN_FLUX", 69.43, LAT),
         LAT = ifelse(SiteSubsite %in% water.LAT.LONG, 74.47427, LAT),
         LONG = ifelse(SiteSubsite == "DISKO:DRYHEATH_FLUX", -53.45, LONG),
         LONG = ifelse(SiteSubsite == "DISKO:WETFEN_FLUX", -53.78, LONG),
         LONG = ifelse(SiteSubsite %in% water.LAT.LONG,-20.52895, LONG),
         LAT = ifelse(SiteSubsite == "NARSARSUAQ:HIGH_ELEVATION", 61.16, LAT),
         LONG = ifelse(SiteSubsite == "NARSARSUAQ:HIGH_ELEVATION", -45.40, LONG))
# coordinate info taken from paper
# https://cdnsciencepub.com/doi/10.1139/AS-2020-0041

# These have all been inputted now
nolat2 <- itex.full4 %>% filter(is.na(LAT))
nolong2 <- itex.full4 %>% filter(is.na(LONG))
nosize <- itex.full4 %>% filter(is.na(SurveyedArea))



## SPOT CHECKS ----
hist(itex.full4$LAT) # makes sense
hist(itex.full4$LONG) # ok
unique(itex.full4$SiteSubsitePlot)
unique(itex.full4$SITE)
unique(itex.full4$RelCover)

# Calculate number of subsites and plots 
subsite.summary <- itex.full4 %>% group_by(SiteSubsite) %>% 
  mutate(NumberPlots = length(unique(SiteSubsitePlot))) %>%
  ungroup() %>% distinct(SiteSubsite, .keep_all = TRUE) %>% 
  arrange(SiteSubsite) %>% mutate(SubsiteNumber = 1:n()) %>% 
  select(SubsiteNumber, SiteSubsite, NumberPlots)




## SPECIES NAMES ----
spp <- unique(itex.full4$SPECIES_NAME)

# Identify empty species names
empty <- itex.full4 %>% filter(SPECIES_NAME == " ") # no empty cells

# Remove trailing white space so the same species are comparable
itex.full5 <- itex.full4 %>% mutate(SPECIES_NAME = str_trim(SPECIES_NAME))

# Check full species name again
unique(itex.full5$SPECIES_NAME)

# Standardise subspecies to species level
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Ledum palustre subsp. groenlandicum"] <- "Ledum palustre"
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Tephroseris integrifolia subsp. atropurpurea"] <- "Tephroseris integrifolia"
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Silene uralensis subsp. apetala"] <- "Silene uralensis"
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Cardamine bellidifolia subsp. alpina"] <- "Cardamine bellidifolia"
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Carex aquatilis var. minor"] <- "Carex aquatilis"
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Eriophorum angustifolium subsp. triste"] <- "Eriophorum angustifolium"
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Empetrum nigrum subsp. hermaphroditum"] <- "Empetrum nigrum"
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Anthoxanthum odoratum subsp. nipponicum"] <- "Anthoxanthum odoratum"
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Silene uralensis subsp. apetala"] <- "Silene uralensis"  
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Salix lanata subsp. richardsonii"] <- "Salix lanata"  
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Rumex alpestris subsp. lapponicus"] <- "Rumex alpestris"   
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Empetrum nigrum/Empetrum nigrum subsp. hermaphroditum"] <- "Empetrum nigrum"   
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Salix arctica/Salix arctica"] <- "Salix arctica"


# Check that species name are correct
sp.names <- TPL(unique(itex.full5$SPECIES_NAME))
head(sp.names)

sp.names[sp.names$New.Species != sp.names$Species,]
sp.names[sp.names$Typo == "TRUE",]
unique(sp.names$New.Taxonomic.status)
tochange <- sp.names[sp.names$Plant.Name.Index==FALSE,] # these are all good

# Fix spelling errors
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Aconitum delphiniifolium"] <- "Aconitum delphinifolium"
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Chamaeorchis alpina"] <- "Chamorchis alpina"
itex.full5$GENUS[itex.full5$SPECIES_NAME == "Chamaeorchis alpina"] <- "Chamorchis"
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Valeriana stichensis"] <- "Valeriana sitchensis"
itex.full5$SPECIES_NAME[itex.full5$SPECIES_NAME == "Cardamine digitalis"] <- "Cardamine digitata"


# Check if genus is correct so we can convert into morphospecies
genus.vector <- c("Luzula spicata/confusa", "Eriophorum scheuchzeri/chamissonis", 
                  "Empetrum nigrum/Phyllodoce caerulea", "Carex microchaeta&rupestris", "Deschampsia flexuosa/Juncus trifidus")

# These will be fixed below
genus.check <- itex.full5 %>% filter(SPECIES_NAME %in% genus.vector)


# Convert genus/family species names to morphospecies
morp.vector <- scan(text = "Alchemilla
                 Anemone 
                 Antennaria 	
                 Arnica
                 Astragalus
                 Calamagrostis
                 Cardamine
                 Carex
                 Cyperaceae
                 Deschampsia
                 Draba
                 Dryas
                 Epilobium
                 Eriophorum
                 Festuca
                 Galium
                 Gentiana
                 Hedysarum
                 Hepatica
                 Luzula
                 Minuartia
                 Oxytropis
                 Pedicularis
                 Petasites
                 Poa	
                 Poa sp.
                 Poaceae
                 Polemonium
                 Polygonum
                 Salix
                 Saxifraga
                 Senecio
                 Stellaria
                 Tofieldia
                 Viola", what="")

# This function is AMAZING! Never adding quotes and commas manually again!
morp.vector2 <- c(morp.vector, "Juncus NA", "Minuartia NA", "Sagina NA", 
                  "Luzula spicata/confusa", "Eriophorum scheuchzeri/chamissonis", "Carex microchaeta&rupestris", 
                  "Empetrum nigrum/Phyllodoce caerulea", "Deschampsia flexuosa/Juncus trifidus", "Graminoid unknown")

# Transform into morphospecies
itex.full6 <- itex.full5 %>% 
  mutate(SPECIES_NAME = ifelse(SPECIES_NAME %in% morp.vector2, paste0("XXX", SPECIES_NAME, ":", SITE), SPECIES_NAME))

# Double-check species names
spppp <- unique(itex.full6$SPECIES_NAME)


# Vector with the morphospecies that don't include site name
nosite.vector <- c("XXXCOMP", "XXXRUSH", "XXXUNKGRA", "XXXPOA", "XXXkobresia", "XXXotherherb", "XXXotheraster", "XXXGRASSUNK", 
                   "XXXasteraceae", "XXXoxytropis", "XXXotherforb", "XXXcarex", "XXXluzula", "XXXSaxifraga", "XXXDraba", 
                   "XXXGRASS", "XXXUNKDIC", "XXXSEDGE", "XXXFORB", "XXXLEGSPP", "XXXUNKFOR", "XXXGRASSUNK", "XXXGRAMINOID", 
                   "XXXLuzula", "XXXothergram", "XXXped", "XXXunkgram", "XXXUNIFOR", "XXXRanunculaceae", "XXXunkOxytropis")

# A few manual fixes
itex.full7 <- itex.full6 %>%
  mutate(SPECIES_NAME = ifelse(SPECIES_NAME == "XXXJuncus NA:THUFUVER", "XXXJuncus:THUFUVER", 
                               ifelse(SPECIES_NAME == "XXXMinuartia NA:THUFUVER", "XXXMinuartia:THUFUVER", 
                                      ifelse(SPECIES_NAME == "XXXSagina NA:THUFUVER", "XXXSagina:THUFUVER", 
                                             ifelse(SPECIES_NAME == "XXXGraminoid unknown:THINGVELLIR", "XXXGRAMINOID:THINGVELLIR", SPECIES_NAME))))) %>%
  mutate(SPECIES_NAME = ifelse(SPECIES_NAME %in% nosite.vector, paste0(SPECIES_NAME, ":", SITE), SPECIES_NAME))

sort(unique(itex.full7$SPECIES_NAME))
unique(itex.full7$GFNARROWwalker)


# Check that 1 func group per species and per genus

# How many func groups per genus?
gf.check <- itex.full7 %>% 
  group_by(GENUS) %>% 
  mutate(count_fg = length(unique(GFNARROWwalker))) %>% 
  ungroup() %>% 
  filter(count_fg > 1, !is.na(GENUS)) %>% 
  dplyr::select(SPECIES_NAME, GENUS, GFNARROWwalker, count_fg) %>% 
  unique() %>% 
  arrange(GENUS)


# Manually correct incorrect records
itex.full70 <- itex.full7 %>% 
  mutate(GFNARROWwalker = case_when(SPECIES_NAME == "Thymus paecox" ~ "SEVER", 
                                    SPECIES_NAME == "Vaccinium vitis-idaea" ~ "SEVER", 
                                    SPECIES_NAME == "Linnaea borealis" ~ "SEVER", 
                                    SPECIES_NAME == "Kobresia myosuroides" ~ "SEDGE", 
                                    SPECIES_NAME == "Dryas integrifolia" ~ "SEVER",
                                    SPECIES_NAME == "Carex rariflora" ~ "SEDGE", 
                                    SPECIES_NAME == "Potentilla palustris" ~ "FORB", 
                                    SPECIES_NAME == "Potentilla fruticosa" ~ "FORB",
                                    SPECIES_NAME == "XXXLuzula:ZACKENBERG" ~ "RUSH",
                                    SPECIES_NAME == "XXXLuzula:WOLFCREEK" ~ "RUSH",
                                    SPECIES_NAME == "XXXCarex microchaeta&rupestris:WOLFCREEK" ~ "SEDGE", 
                                    TRUE ~ GFNARROWwalker))


## GEO DATA ----
unique(itex.full70$SITE)

# Specify biogeographic region
eurasia <- c("ABISKO", "BILLEFJORDEN", "DOVRE", "ENDALEN", "FINSE", "FURI", "KILPISJARVI", "LATNJA", "LOGH", "LORI", 
             "RIRI", "SADVENT", "JOATKA", "NAKKALA", "NYALESUND", "MALAYA", "TAISETSU", "GAVIA", "STILLBERG", "VALBERCLA",
             "PYRAMIDEN", "BILLEFJORDEN", "ADVENT")
greenice <- c("KANGER", "ZACKENBERG", "DISKO", "AKUREYRI", "AUDKULUHEIDI", "BLONDUOS", "FAROE", "HJARDARLAND", "NARSARSUAQ",
              "HOLTAVORDUHEIDI", "MODRUVELLIR", "OXNADALSHEIDI", "THINGVELLIR", "THUFUVER", "THYKKVIBAER", "NUUK", "JAMESONLAND")
na.east <- c("ALEXFIORD", "TORNGATS", "BYLOT", "DALSMYNNI", "AUYUITTUQ", "IGLOOLIK", "QUTTINIRPAAQ")
na.west <- c("ANWR", "ATQASUK", "BARROW", "BROOKS", "DARING", "KLUANE", "NIWOT", "QHI", "TOOLIK", "WOLFCREEK")


# Add in categories
itex.full8 <- itex.full70 %>% 
  mutate(AlpArc = ifelse(LAT >= 66.4 & ELEV >= 1000, "Arctic-Alpine",
                         ifelse(LAT >= 66.4, "Arctic",
                                ifelse(ELEV >= 1000, "Alpine", "Subarctic")))) %>% 
  mutate(lat_grid = plyr::round_any(LAT, 0.5, f = floor)) %>% 
  mutate(lon_grid = ifelse(LONG >0, plyr::round_any(LONG, 0.5, f = floor), 
                           plyr::round_any(LONG, 0.5, f = ceiling))) %>%
  mutate(gridcell = paste0("_", lat_grid, "_", lon_grid)) %>%
  mutate(Region = ifelse(SITE %in% eurasia, "Eurasia", 
                         ifelse(SITE %in% greenice, "GreenIceLand", 
                                ifelse(SITE %in% na.east, "North America-East",
                                       ifelse(SITE %in% na.west, "North America-West", NA)))))



## MORPHOSPECIES ----
length(unique(itex.full8$SiteSubsitePlotYear)) # 6821 plotXyear in total

# Check the unidentified species names
xxxtest <- itex.full8 %>% 
  filter(str_detect(SPECIES_NAME, 'XXX|xxx')) %>% 
  group_by(SiteSubsitePlotYear) %>% 
  mutate(MorphoCover = sum(RelCover)) %>%
  ungroup() %>% select(SITE, SiteSubsite, SiteSubsitePlot, SiteSubsitePlotYear, MorphoCover) %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE)

# 2028 unique plotXyear that contain morphospecies
length(unique(xxxtest$SiteSubsitePlotYear))

# Check how many plots would be removed with different cut-offs
morpho10 <- xxxtest %>% filter(MorphoCover > 10) # 10% morphospecies cutoff would remove 920 plotXyear (13.4% of database)
morpho15 <- xxxtest %>% filter(MorphoCover > 15) # 15% morphospecies cutoff would remove 682 plotXyear (9.9% of database)
morpho20 <- xxxtest %>% filter(MorphoCover > 20) # 20% morphospecies cutoff would remove 538 plotXyear (7.8% of database)
morpho25 <- xxxtest %>% filter(MorphoCover > 25) # 25% morphospecies cutoff would remove 415 plotXyear (6% of database)

# Mean morphospecies cover is 15.9%
mean(xxxtest$MorphoCover)

# Visualize morphospecies cover across plots
(xxxmorph <- ggplot(xxxtest, aes(x=MorphoCover)) + 
    geom_histogram(binwidth = 0.5) + 
    geom_vline(aes(xintercept = mean(MorphoCover)), colour = "red", linetype = "dashed", size = 1) +
    xlab("Cover of morphospecies in plots"))

# Extract as vector from the main dataset those plotXyear with >10% morphospecies cover
morpho10.vec <- unique(morpho10$SiteSubsitePlotYear)

write.csv(morpho10.vec, "data/morpho10_vector.csv")

# Remove those plots
itex.full9 <- itex.full8 %>% filter(SiteSubsitePlotYear %notin% morpho10.vec)




## FG-DOMINATION ----
unique(itex.full9$ValueType)

# Standardise method and FG categories, add % cover per functional group & specify dominant FG
itex.full10 <- itex.full9 %>% 
  mutate(ValueType = case_when(ValueType == "pf_all_xy" ~ "pf_all_XY", 
                               ValueType == "pf_topbot_xy" ~ "pf_topbot_XY", 
                               TRUE ~ ValueType)) %>%
  mutate(Method = case_when(ValueType %in% c("BraunBlanquet", "percent_cover") ~ "Cover",
                            ValueType %in% c("pf_all_plot", "pf_top_plot", "pf_topbot_plot") ~ "Point-framing (sum)",
                            ValueType %in% c("pf_all_XY", "pf_top_XY", "pf_topbot_XY") ~ "Point-framing (XY)")) %>%
  mutate(FuncGroup = case_when(GFNARROWwalker %in% c("SEVER", "SDECI") ~ "Shrub",
                               GFNARROWwalker == "FORB" ~ "Forb",
                               GFNARROWwalker %in% c("RUSH", "SEDGE", "GRASS", "GRAMINOIDU", "GRAMU") ~ "Graminoid")) %>%
  group_by(SiteSubsitePlotYear, FuncGroup) %>% 
  mutate(PlotYearFuncCover = sum(RelCover)) %>% ungroup() %>%
  mutate(PlotYearDominatingFG = case_when(PlotYearFuncCover > 50 ~ paste0(FuncGroup, "-Dominated"), TRUE ~ "None")) # we need this column for homogenization

# This is where we specify the dominant FG for that particular plotXyear: long format including all cover values per FG
dom.fg <- itex.full10 %>% distinct(SiteSubsitePlotYear, FuncGroup, .keep_all = TRUE) %>% 
  select(SiteSubsitePlotYear, FuncGroup, PlotYearFuncCover) %>%
  pivot_wider(names_from = FuncGroup, values_from = PlotYearFuncCover, values_fill = list(PlotYearFuncCover = 0)) %>%
  mutate(DominatingFG = case_when(Shrub > 50 ~ "Shrub-Dominated",
                                  Graminoid > 50 ~ "Graminoid-Dominated",
                                  Forb > 50 ~ "Forb-Dominated", TRUE ~ "None")) %>%
  rename(ShrubCover = Shrub) %>% rename(GraminoidCover = Graminoid) %>% rename(ForbCover = Forb)

# Join with full dataset
itex.full11 <- left_join(itex.full10, dom.fg, by = "SiteSubsitePlotYear")



## MEAN COVER ----

# Calculate mean cover per plot over time (more representative of FG cover over time)
avg.fg <- itex.full11 %>% distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  group_by(SiteSubsitePlot) %>% mutate(PlotShrubMean = mean(ShrubCover)) %>% 
  mutate(PlotForbMean = mean(ForbCover)) %>% mutate(PlotGraminoidMean = mean(GraminoidCover)) %>%
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>% ungroup() %>%
  select(SiteSubsitePlot, PlotShrubMean, PlotForbMean, PlotGraminoidMean)

# Join with main dataset
itex.full12 <- left_join(itex.full11, avg.fg, by = "SiteSubsitePlot")

# Not removing the rows with cover = 0 fully because there might be a total cover of 0 in a plot and we would be removing a whole plotXyear survey
# So just filtering below to avoid any issues with richness calculation

# Calculate average richness over time and round to nearest integer
mean.rich <- itex.full12 %>% filter(RelCover > 0) %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(AnnualRichness = length(unique(SPECIES_NAME))) %>% ungroup() %>%
  distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(MeanRichness = round(mean(AnnualRichness))) %>% 
  ungroup() %>% distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  dplyr::select(SiteSubsitePlotYear, AnnualRichness, MeanRichness)

# Join with main dataset
itex.dec22 <- left_join(itex.full12, mean.rich, by = "SiteSubsitePlotYear") #42234

# this is the main input file for the rest of scripts
save(itex.dec22, file = "data/itex_dec22.RData")




## FINAL CHECKS ----

# Calculate number of subsites and plots
subsite.summary.nomp <- itex.dec22 %>% group_by(SiteSubsite) %>% 
  mutate(NumberPlots = length(unique(SiteSubsitePlot))) %>%
  ungroup() %>% distinct(SiteSubsite, .keep_all = TRUE) %>% 
  arrange(SiteSubsite) %>% mutate(SubsiteNumber = 1:n()) %>% 
  select(SubsiteNumber, SiteSubsite, NumberPlots)




## SPATIAL DATASET ----

# Retain last year of monitoring only
itex.end <- itex.dec22 %>% filter(RelCover > 0) %>% 
  group_by(SiteSubsitePlot) %>% 
  filter(YEAR == max(YEAR)) %>% ungroup()

write.csv(itex.end, "data/end_only_itex22.csv") #15672

# Calculating end richness
itex.end.rich <- itex.end %>% group_by(SiteSubsitePlot) %>% 
  mutate(CurrentRichness = length(unique(SPECIES_NAME))) %>% ungroup() %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>% 
  dplyr::select(., -c(STATUS, TISSUE, ORIGINAL_NAME, GENUS, GFNARROWwalker, 
                      SPECIES_NAME, TREATMENT, RelCover, FuncGroup, PlotYearFuncCover, PlotYearDominatingFG))

write.csv(itex.end.rich, "data/end_itex_rich22.csv") #2174





