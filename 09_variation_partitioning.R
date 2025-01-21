####
#### R script for Ohigashi et al (2024)
#### Variation partitioning of community composition and functional composition
#### 2024.11.03 written by Ohigashi
#### R 4.3.3
#### ref https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html


### load packages and functions
source("Function/F1_HelperFunctions.R")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(ade4); packageVersion("ade4")
library(adespatial); packageVersion("adespatial")
library(adegraphics); packageVersion("adegraphics")
library(spdep); packageVersion("spdep")
library(sp); packageVersion("sp")
library(vegan); packageVersion("vegan")


### load data
# sample location
gpsdata <- read.csv("Data/soil_gpsdata.csv", header = T)

# metadata
sample_data <- read.table("Data/soil_metadata.txt", header = T)

# ASV tables
# prokaryotes
b_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T)
# fungi
f_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T)

# function tables
# prokaryotic function
prok_func_all <- read.table("02_Function_analysis_out/PICRUSt2_full_function.txt", header = T)
# fungal lifestyle
fungilife <- read.table("02_Function_analysis_out/FungalTraits_primarilifestyle_percent.txt", header = T)


### format XY data
## add XY coordinates
gpsdata <- gpsdata |>
  mutate(xy_lon=lapply(Coord_long, dms_to_deg),
         xy_lat=lapply(Coord_lat, dms_to_deg)
  )

# add plot column to sample_data to merge with gpsdata
sample_data <- sample_data |>
  mutate(Plot = rep(rep(c("1", "2", "3"), each = 3), times = 10))

# merge
sample_data <- merge(sample_data, gpsdata, by=c("Country", "Site", "Landuse", "Plot"), all.x = T, sort=F)

# add Reps column
sample_data <- sample_data |> mutate(Reps = rep(rep(c("1", "2", "3")), times = 30))
env_data <- sample_data |> select(Sample, Country, Site, Landuse, Plot, Reps, 
                                  xy_lon, xy_lat)
env_data$xy_lon <- as.numeric(env_data$xy_lon)
env_data$xy_lat <- as.numeric(env_data$xy_lat)

# adjust longitude and latitude for the replications
env_data <- env_data %>%
  group_by(Country, Site, Landuse, Plot) %>%
  mutate(xy_lon = case_when(
    Reps == 1 ~ xy_lon + 0,
    Reps == 2 ~ xy_lon + 0.000003,
    Reps == 3 ~ xy_lon + 0.000006,
    TRUE ~ xy_lat),
    xy_lat = case_when(
      Reps == 1 ~ xy_lat + 0,
      Reps == 2 ~ xy_lat + 0.000003,
      Reps == 3 ~ xy_lat + 0.000006,
      TRUE ~ xy_lat)
  ) |>
  ungroup()

# get XY coordinates
mxy <- cbind(as.numeric(env_data$xy_lon), as.numeric(env_data$xy_lat))
colnames(mxy) <- c("x", "y")

# mxy for Natural and Farm
mxy.nat <- env_data |>
  filter(Sample %in% nat_samples) |>
  select(x = xy_lon, y = xy_lat) |>
  as.matrix()
mxy.farm <- env_data |>
  filter(Sample %in% farm_samples) |>
  select(x = xy_lon, y = xy_lat) |>
  as.matrix()


### format ASV and function tables
# get Sample names for Natural and Farm samples
nat_samples <- sample_data %>%
  filter(Landuse == "Natural") %>%
  pull(Sample)
farm_samples <- sample_data %>%
  filter(Landuse == "Farm") %>%
  pull(Sample)

# create data frame for natural and farm communities
b_ASV.nat <- b_ASV.table |> select(all_of(nat_samples))
b_ASV.farm <- b_ASV.table |> select(all_of(farm_samples))
b_ASV.nat.t <- t(b_ASV.nat) # transpose
b_ASV.farm.t <- t(b_ASV.farm) # transpose

f_ASV.nat <- f_ASV.table |> select(all_of(nat_samples))
f_ASV.farm <- f_ASV.table |> select(all_of(farm_samples))
f_ASV.nat.t <- t(f_ASV.nat) # transpose
f_ASV.farm.t <- t(f_ASV.farm) # transpose

# create data frame for natural and farm functions
prok_func_all <- prok_func_all |> rownames_to_column(var = "Sample")
b_func_all.nat <- prok_func_all |>
  filter(Sample %in% nat_samples) |>
  column_to_rownames(var = "Sample") |>
  as.matrix()
b_func_all.farm <- prok_func_all |>
  filter(Sample %in% farm_samples) |>
  column_to_rownames(var = "Sample") |>
  as.matrix()

fungilife <- fungilife |> filter(primary_lifestyle != "unassigned") # 2024.09.11
fungilife <- fungilife |> column_to_rownames(var = "primary_lifestyle")
fungilife.nat <- fungilife |> select(all_of(nat_samples))
fungilife.farm <- fungilife |> select(all_of(farm_samples))

fungilife.nat.t <- t(fungilife.nat)
fungilife.farm.t <- t(fungilife.farm)


### format environmental variables
env_vars <- sample_data |>
  select(Landuse, Gravimetric.water.content, pH, Carbon, Nitrogen, CN_ratio)

env.nat <- env_vars |>
  filter(Landuse == "Natural") |>
  select(-Landuse) |>
  as.matrix()

env.farm <- env_vars |>
  filter(Landuse == "Farm") |>
  select(-Landuse) |>
  as.matrix()


### spatial analysis
# # define spatial weighting matrices
# nbgab <- graph2nb(gabrielneigh(mxy), sym = TRUE) # consider the Gabriel graph
# nb2listw(nbgab)
# 
# # get distance between neighbors
# distgab <- nbdists(nbgab, mxy)
# # nbgab[[1]]
# # distgab[[1]]
# 
# # define a function to do weight by distance
# fdist <- lapply(distgab, function(x) 1 - x/max(dist(mxy)))
# 
# # spatial weighting matrix 
# listwgab <- nb2listw(nbgab, glist = fdist)
# listwgab
# 
# 
# ### MEM
# mem.gab <- mem(listwgab)
# mem.gab
# 
# # set rng.seed
# set.seed(123)
# 
# # calculate Maran's I
# moranI <- moran.randtest(mem.gab, listwgab, 99)
# significant_indices <- order(moranI$pvalue)[moranI$pvalue < 0.05] # subset by significance
# 
# # draw the result
# dir.create("09_variation_partitioning_out")
# png(filename = "09_variation_partitioning_out/mems_significant.png", width = 1000, height = 1000, res = 300)
# plot(mem.gab[,significant_indices], SpORcoords = mxy)
# dev.off()


######### variation partitioning ########
set.seed(123)
## 1. Prokaryotic communities
# 1-1. Prokaryotes (Natural)
pca.hell.prok.nat <- dudi.pca(b_ASV.nat.t, scale = FALSE, scannf = FALSE, nf = 2)
cand.lw.nat <- listw.candidates(mxy.nat, nb = c("gab", "rel"), weights = c("bin", "flin"))
sel.lw.prok.nat <- listw.select(pca.hell.prok.nat$tab, candidates = cand.lw.nat, nperm = 99)
vp_prok.nat <- varpart(pca.hell.prok.nat$tab, env.nat, sel.lw.prok.nat$best$MEM.select)
vp_prok.nat

# 1-2. Prokaryotes (Farm)
pca.hell.prok.farm <- dudi.pca(b_ASV.farm.t, scale = FALSE, scannf = FALSE, nf = 2)
cand.lw.farm <- listw.candidates(mxy.farm, nb = c("gab", "rel"), weights = c("bin", "flin"))
sel.lw.prok.farm <- listw.select(pca.hell.prok.farm$tab, candidates = cand.lw.farm, nperm = 99)
vp_prok.farm <- varpart(pca.hell.prok.farm$tab, env.farm, sel.lw.prok.farm$best$MEM.select)
vp_prok.farm


## 2. Fungal communities
# 2-1. Fungi (Natural)
pca.hell.fungi.nat <- dudi.pca(f_ASV.nat.t, scale = FALSE, scannf = FALSE, nf = 2)
cand.lw.nat <- listw.candidates(mxy.nat, nb = c("gab", "rel"), weights = c("bin", "flin"))
sel.lw.fungi.nat <- listw.select(pca.hell.fungi.nat$tab, candidates = cand.lw.nat, nperm = 99)
vp_fungi.nat <- varpart(pca.hell.fungi.nat$tab, env.nat, sel.lw.fungi.nat$best$MEM.select)
vp_fungi.nat

# 2-2. Fungi (Farm)
pca.hell.fungi.farm <- dudi.pca(f_ASV.farm.t, scale = FALSE, scannf = FALSE, nf = 2)
cand.lw.farm <- listw.candidates(mxy.farm, nb = c("gab", "rel"), weights = c("bin", "flin"))
sel.lw.fungi.farm <- listw.select(pca.hell.fungi.farm$tab, candidates = cand.lw.farm, nperm = 99)
vp_fungi.farm <- varpart(pca.hell.fungi.farm$tab, env.farm, sel.lw.fungi.farm$best$MEM.select)
vp_fungi.farm


## 3. Prokaryotic functions
# 3-1. Prokaryotic functions (Natural)
pca.hell.prokfunc.nat <- dudi.pca(b_func_all.nat, scale = FALSE, scannf = FALSE, nf = 2)
cand.lw.nat <- listw.candidates(mxy.nat, nb = c("gab", "rel"), weights = c("bin", "flin"))
sel.lw.prokfunc.nat <- listw.select(pca.hell.prokfunc.nat$tab, candidates = cand.lw.nat, nperm = 99)
vp_prokfunc.nat <- varpart(pca.hell.prokfunc.nat$tab, env.nat, sel.lw.prokfunc.nat$best$MEM.select)
vp_prokfunc.nat

# 3-2. Prokaryotic functions (Farm)
pca.hell.prokfunc.farm <- dudi.pca(b_func_all.farm, scale = FALSE, scannf = FALSE, nf = 2)
cand.lw.farm <- listw.candidates(mxy.farm, nb = c("gab", "rel"), weights = c("bin", "flin"))
sel.lw.prokfunc.farm <- listw.select(pca.hell.prokfunc.farm$tab, candidates = cand.lw.farm, nperm = 99)
vp_prokfunc.farm <- varpart(pca.hell.prokfunc.farm$tab, env.farm, sel.lw.prokfunc.farm$best$MEM.select)
vp_prokfunc.farm


## 4. Fungal communities
# 4-1. Fungi (Natural)
pca.hell.fungilife.nat <- dudi.pca(fungilife.nat.t, scale = FALSE, scannf = FALSE, nf = 2)
cand.lw.nat <- listw.candidates(mxy.nat, nb = c("gab", "rel"), weights = c("bin", "flin"))
sel.lw.fungilife.nat <- listw.select(pca.hell.fungilife.nat$tab, candidates = cand.lw.nat, nperm = 99)
vp_fungilife.nat <- varpart(pca.hell.fungilife.nat$tab, env.nat, sel.lw.fungilife.nat$best$MEM.select)
vp_fungilife.nat

# 4-2. Fungi (Farm)
pca.hell.fungilife.farm <- dudi.pca(fungilife.farm.t, scale = FALSE, scannf = FALSE, nf = 2)
cand.lw.farm <- listw.candidates(mxy.farm, nb = c("gab", "rel"), weights = c("bin", "flin"))
sel.lw.fungilife.farm <- listw.select(pca.hell.fungilife.farm$tab, candidates = cand.lw.farm, nperm = 99)
vp_fungilife.farm <- varpart(pca.hell.fungilife.farm$tab, env.farm, sel.lw.fungilife.farm$best$MEM.select)
vp_fungilife.farm


### sumarize and save data
dir.create("09_variation_partitioning_out")
# list of the results
pcalist <- list("prokcom.nat" = pca.hell.prok.nat,
                "prokcom.farm" = pca.hell.prokfunc.farm,
                "fungicom.nat" = pca.hell.fungi.nat,
                "fungicom.farm" = pca.hell.fungi.farm,
                "prokfunc.nat" = pca.hell.prokfunc.nat,
                "prokfunc.farm" = pca.hell.prokfunc.farm,
                "fungifunc.nat" = pca.hell.fungilife.nat,
                "fungifunc.farm" = pca.hell.fungilife.farm)
sel.lw.list <- list(sel.lw.prok.nat, sel.lw.prok.farm,
                    sel.lw.fungi.nat, sel.lw.fungi.farm,
                    sel.lw.prokfunc.nat, sel.lw.prokfunc.farm,
                    sel.lw.fungilife.nat, sel.lw.fungilife.farm)
vp_list <- list(vp_prok.nat, vp_prok.farm,
                vp_fungi.nat, vp_fungi.farm,
                vp_prokfunc.nat, vp_prokfunc.farm,
                vp_fungilife.nat, vp_fungilife.farm)

for (i in 1:length(pcalist)) {
  # get results from the lists
  pcares <- pcalist[[i]]
  sel.lw.res <- sel.lw.list[[i]]
  
  # set environmental variables
  if (grepl("nat", names(pcalist)[i])) {
    env <- env.nat
  } else {
    env <- env.farm
  }
  
  # perform CCA
  cca_env <- anova.cca(rda(pcares$tab, env))
  rownames(cca_env)[1] <- "Env(a+c)"
  # (space w/o controlling env)
  cca_space <- anova.cca(rda(pcares$tab, sel.lw.res$best$MEM.select))
  rownames(cca_space)[1] <- "Space(b+c)"
  # (env alone)
  cca_env_alone <- anova.cca(rda(pcares$tab, env, sel.lw.res$best$MEM.select))
  rownames(cca_env_alone)[1] <- "EnvAlone(a)"
  # (space alone)
  cca_space_alone <- anova.cca(rda(pcares$tab, sel.lw.res$best$MEM.select, env))
  rownames(cca_space_alone)[1] <- "SpaceAlone(b)"
  
  cca_res <- rbind(cca_env[1,],
                   cca_space[1,],
                   cca_env_alone[1,],
                   cca_space_alone[1,])
  write.csv(cca_res, file = sprintf("09_variation_partitioning_out/vp_stats_%s.csv", names(pcalist)[i]))
  
  # result of variation partitioning
  vp <- vp_list[[i]]
  R2df <- vp$part$indfract
  rownames(R2df) <- c("EnvAlone", "SpaceAlone", "Overlap", "Residuals")
  write.csv(R2df, file = sprintf("09_variation_partitioning_out/vp_R2_%s.csv", names(pcalist)[i]), quote = F, row.names = T)
}


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/09_SessionInfo_%s.txt", substr(Sys.time(), 1, 10))) 

