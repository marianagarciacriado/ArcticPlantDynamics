## Plant diversity dynamics over space and time in a warming Arctic
## Mariana Garcia Criado (mariana.garcia.criado@gmail.com)
## Script 8. Model summary table


## LIBRARIES ----
library(dplyr)
library(broom)
library(stargazer)

# Let's break this down since otherwise R crashes
# Load a few models, run the table, and then close R and clean the environment
# Rinse and repeat with the rest of models



## FUNCTION ----

# Jonathan Chang's function
p_summarize <- function(model) {
  brms::posterior_summary(model) %>% 
    as_tibble(rownames = "parameter")
}



## SCRIPTS 2 & 3 MODEL TABLE ----

# Add model objects to list
models.list23 <- list(rich.geo.mod, rich.clim.mod, 
                     rich.fg.shb, rich.fg.frb, rich.fg.gram, 
                     rich.time.subsite.mod, 
                     rich.time.shbch.plot.mod, rich.time.forbch.plot.mod, rich.time.gramch.plot.mod, 
                     jac.geo.mod, jac.clim.mod, 
                     jac.shb.mod, jac.forb.mod, jac.gram.mod, 
                     jac.shbchg.mod, jac.forbchg.mod, jac.gramchg.mod, 
                     bray.geo.mod, bray.clim.mod,
                     bray.shb.mod, bray.forb.mod, bray.gram.mod,
                     bray.shb.chg.mod, bray.forb.chg.mod, bray.gram.chg.mod)


# Create dataframe with model names
models.name23 <- c("Richness GEO", "Richness CLIM", 
                  "Richness FG - Shrub", "Richness FG - Forb", "Richness FG - Graminoid",
                  "Richness change SUBS", 
                  "Richness change PCHG - Shrub", "Richness change PCHG - Forb", "Richness change PCHG - Graminoid",
                  "Jaccard GEO", "Jaccard CLIM",
                  "Jaccard FG - Shrub", "Jaccard FG - Forb", "Jaccard FG - Graminoid",
                  "Jaccard PCHG - Shrub", "Jaccard PCHG - Forb", "Jaccard PCHG - Graminoid",
                  "Bray-Curtis GEO", "Bray-Curtis CLIM",
                  "Bray-Curtis FG - Shrub", "Bray-Curtis FG - Forb", "Bray-Curtis FG - Graminoid",
                  "Bray-Curtis PCHG - Shrub", "Bray-Curtis PCHG - Forb", "Bray-Curtis PCHG - Graminoid")


# Bind them together
model_number23 <- 1:25
mod.df23 <- data.frame(model_number23, models.name23)

# Extract parameters
mod.table23 <- lapply(models.list23, p_summarize) %>% 
  bind_rows(.id = "model_number23") 

# Add model name to table
mod.table23$model_number23 <- as.integer(mod.table23$model_number23)
mod.table.final23 <- left_join(mod.table23, mod.df23, by = "model_number23")

# Clean model parameters
mod.table.final23.2 <- mod.table.final23 %>% dplyr::filter(parameter != "lp__") %>% 
  dplyr::filter(!str_detect(parameter, "^r_"))

mod.table.final23.2$models.name23[duplicated(mod.table.final23.2$models.name23)] <- "  "
mod.table.final23.2$models.name23 <- as.character(mod.table.final23.2$models.name23)
mod.table.final23.2$model_number23[duplicated(mod.table.final23.2$model_number23)] <- "  "

colnames(mod.table.final23.2) <- c("Model number", "Term", "Estimate", "Std. error", "Lower 95% CI", "Upper 95% CI", "Model name")
mod.table.final23.3 <- mod.table.final23.2[, c(1, 7, 2, 3, 4, 5, 6)]

#  Round to 3 decimals only because not working on stargazer function
mod.table.final23.4 <- mod.table.final23.3 %>% mutate_if(is.numeric, round, digits = 3)

# Save in csv
write.csv(mod.table.final23.4, "models/table_s8/table_s8_23.csv")

# Convert to table
# stargazer(mod.table.final23.4, type = "html", summary = FALSE) stargazer is not working due to underscores in names




## SCRIPTS 4 MODEL TABLE ----

# Add model objects to list
models.list4 <- list(pers.geo.mod, pers.clim.mod,
                     pers.fg.shb.mod, pers.fg.forb.mod, pers.fg.gram.mod,
                     pers.shrubchange.mod, pers.forbchange.mod, pers.gramchange.mod,
                     ext.geo.mod, ext.clim.mod,
                     ext.fg.shb.mod, ext.fg.forb.mod, ext.fg.gram.mod,
                     ext.shbchange.mod, ext.forbchange.mod, ext.gramchange.mod,
                     col.geo.mod, col.clim.mod,
                     col.fg.shb.mod, col.fg.forb.mod, col.fg.gram.mod,
                     col.shbchange.mod, col.forbchange.mod, col.gramchange.mod)


# Create dataframe with model names
models.name4 <- c("Persisters GEO", "Persisters CLIM",
                   "Persisters FG - Shrub", "Persisters FG - Forb", "Persisters FG - Graminoid",
                   "Persisters PCHG - Shrub", "Persisters PCHG - Forb", "Persisters PCHG - Graminoid",
                   "Losses GEO", "Losses CLIM",
                   "Losses FG - Shrub", "Losses FG - Forb", "Losses FG - Graminoid",
                   "Losses PCHG - Shrub", "Losses PCHG - Forb", "Losses PCHG - Graminoid",
                   "Gains GEO", "Gains CLIM",
                   "Gains FG - Shrub", "Gains FG - Forb", "Gains FG - Graminoid",
                   "Gains PCHG - Shrub", "Gains PCHG - Forb", "Gains PCHG - Graminoid")


# Bind them together
model_number4 <- 26:49
mod.df4 <- data.frame(model_number4, models.name4)

# Extract parameters
mod.table4 <- lapply(models.list4, p_summarize) %>% 
  bind_rows(.id = "model_number4") 

# Add model name to table
mod.table4$model_number4 <- as.integer(mod.table4$model_number4)
mod.table4$model_number4 <- (mod.table4$model_number4 + 25)
mod.table.final4 <- left_join(mod.table4, mod.df4, by = "model_number4")

# Clean model parameters
mod.table.final4.2 <- mod.table.final4 %>% filter(parameter != "lp__") %>% 
  filter(!str_detect(parameter, "^r_"))

mod.table.final4.2$models.name4[duplicated(mod.table.final4.2$models.name4)] <- "  "
mod.table.final4.2$models.name4 <- as.character(mod.table.final4.2$models.name4)
mod.table.final4.2$model_number4[duplicated(mod.table.final4.2$model_number4)] <- "  "

colnames(mod.table.final4.2) <- c("Model number", "Term", "Estimate", "Std. error", "Lower 95% CI", "Upper 95% CI", "Model name")
mod.table.final4.3 <- mod.table.final4.2[, c(1, 7, 2, 3, 4, 5, 6)]

#  Round to 3 decimals only because not working on stargazer function
mod.table.final4.4 <- mod.table.final4.3 %>% mutate_if(is.numeric, round, digits = 3)

# Save in csv
write.csv(mod.table.final4.4, "models/table_s8/table_s8_4.csv")




## SCRIPT 5 & 6 MODEL TABLE ----

# Add model objects to list
models.list56 <- list(even.geo.mod, even.clim.mod, 
                     even.shb.mod, even.forb.mod, even.gram.mod, 
                     even.time.subsite.mod,
                     even.time.shbch.plot.mod, even.time.forbch.plot.mod, even.time.gramch.plot.mod, 
                     jac.subs.mod, jac.subs.forb.mod, jac.subs.gram.mod, 
                     bray.subs.mod, bray.subs.forb.mod, bray.subs.gram.mod,
                     jac.cent.subs.mod, jac.cent.subs.forb.mod, jac.cent.subs.gram.mod,
                     bc.cent.subs.mod, bc.cent.subs.forb.mod, bc.cent.subs.gram.mod)


# Create dataframe with model names
models.name56 <- c("Evenness GEO", "Evenness CLIM",
                  "Evenness FG - Shrub", "Evenness FG - Forb", "Evenness FG - Graminoid",
                  "Evenness SUBS",
                  "Evenness PCHG - Shrub", "Evenness PCHG - Forb", "Evenness PCHG - Graminoid", 
                  "PCoA Jaccard Cartesian SUBS - Shrub", "PCoA Jaccard Cartesian SUBS - Forb", "PCoA Jaccard Cartesian SUBS - Graminoid",
                  "PCoA Bray-Curtis Cartesian SUBS - Shrub", "PCoA Bray-Curtis Cartesian SUBS - Forb", "PCoA Bray-Curtis Cartesian SUBS - Graminoid",
                  "PCoA Jaccard Centroid SUBS - Shrub", "PCoA Jaccard Centroid SUBS - Forb", "PCoA Jaccard Centroid SUBS - Graminoid",
                  "PCoA Bray-Curtis Centroid SUBS - Shrub", "PCoA Bray-Curtis Centroid SUBS - Forb", "PCoA Bray-Curtis Centroid SUBS - Graminoid")




# Bind them together
model_number56 <- 50:70
mod.df56 <- data.frame(model_number56, models.name56)

# Extract parameters
mod.table56 <- lapply(models.list56, p_summarize) %>% 
  bind_rows(.id = "model_number56") 

# Add model name to table
mod.table56$model_number56 <- as.integer(mod.table56$model_number56)
mod.table56$model_number56 <- (mod.table56$model_number56 + 49)
mod.table.final56 <- left_join(mod.table56, mod.df56, by = "model_number56")


# Clean model parameters
mod.table.final56.2 <- mod.table.final56 %>% filter(parameter != "lp__") %>% 
  filter(!stringr::str_detect(parameter, "^r_"))

mod.table.final56.2$models.name56[duplicated(mod.table.final56.2$models.name56)] <- "  "
mod.table.final56.2$models.name56 <- as.character(mod.table.final56.2$models.name56)
mod.table.final56.2$model_number56[duplicated(mod.table.final56.2$model_number56)] <- "  "

colnames(mod.table.final56.2) <- c("Model number", "Term", "Estimate", "Std. error", "Lower 95% CI", "Upper 95% CI", "Model name")
mod.table.final56.3 <- mod.table.final56.2[, c(1, 7, 2, 3, 4, 5, 6)]

#  Round to 3 decimals only because not working on stargazer function
mod.table.final56.4 <- mod.table.final56.3 %>% mutate_if(is.numeric, round, digits = 3)

# Save in csv
write.csv(mod.table.final56.4, "models/table_s8/table_s8_56.csv")




## Extra analysis (richness change ~ shrub categories interaction)
mod.table.extra <- lapply(rich.shb.fg, p_summarize)

# Add model objects to list
models.list.extra <- list(rich.shb.fg)

# Create dataframe with model names
models.name.extra <- c("extra mod")

# Bind them together
model_number_extra <- 1
mod.df_extra <- data.frame(model_number_extra, models.name.extra)

# Extract model outputs
mod.table_extra <- lapply(models.list.extra, p_summarize) %>% 
  bind_rows(.id = "model_extra")


## Models are put together in Excel and then in Word for Table S3.