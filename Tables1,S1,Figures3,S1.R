rm(list = ls())

setwd("~/Documents/lab/spp_ages")

library(climateStability)
library(phytools)
library(geiger)
library(ggplot2)
library(maptools)
library(tidyverse)
library(stringi)
library(raster)
library(rgdal)
library(sf)
library(fasterize)
library(data.table)
library(caper)
library(patchwork)
library(terra)
library(largeList)
library(dplyr)
library(rgeos)

options(scipen=999)

ages <- read.csv("tables/ages.csv")

# Subsetting the geographical data
names_am <- ages$Species[ages$Taxon %in% c("Anuta", "Caudata", "Gymnophiona")]
names_sq <- ages$Species[ages$Taxon %in% c("Anguimorpha", "Gekkota", "Iguania",
                                           "Lacertoidea", "Scincoidea", "Serpentes")]
names_av <- ages$Species[ages$Taxon %in% c("Columbiformes", "Passeriformes", 
                                           "Piciformes", "Psittaciformes")]
names_ma <- ages$Species[ages$Taxon %in% c("Carnivora", "Cetartiodactyla", 
                                           "Chiroptera", "Diprotodontia", "Primates")]

map_am <- readOGR("~/Documents/lab/data/spatial/AMPHIBIANS/AMPHIBIANS.shp",
                  dropNULLGeometries = TRUE)
map_sq1 <- readOGR("~/Documents/lab/data/spatial/SQUAMATA/SCALED_REPTILES_PART1.shp", 
                   dropNULLGeometries = TRUE)
map_sq2 <- readOGR("~/Documents/lab/data/spatial/SQUAMATA/SCALED_REPTILES_PART2.shp",
                   dropNULLGeometries = TRUE)
map_av <- st_read(dsn = "~/Documents/lab/data/spatial/BOTW/BOTW.gdb", 
                  layer = "All_Species")
map_ma <- readOGR("~/Documents/lab/data/spatial/MAMMALS/MAMMALS.shp",
                  dropNULLGeometries = TRUE)

map_av <- map_av %>% filter(st_geometry_type(Shape) != "MULTISURFACE")

map_am$sci_name <- stri_replace_all_fixed(map_am$sci_name, " ", "_")
map_sq1$sci_name <- stri_replace_all_fixed(map_sq1$sci_name, " ", "_")
map_sq2$sci_name <- stri_replace_all_fixed(map_sq2$sci_name, " ", "_")
map_av$sci_name <- stri_replace_all_fixed(map_av$sci_name, " ", "_")
map_ma$sci_name <- stri_replace_all_fixed(map_ma$sci_name, " ", "_")

map2_am <- raster::subset(map_am, map_am$sci_name %in% names_am)
map2_sq1 <- raster::subset(map_sq1, map_sq1$sci_name %in% names_sq)
map2_sq2 <- raster::subset(map_sq2, map_sq2$sci_name %in% names_sq)
map2_av <- raster::subset(map_av, map_av$sci_name %in% names_av)
map2_ma <- raster::subset(map_ma, map_ma$sci_name %in% names_ma)

ages_sq <- ages_sq[ages_sq$Species %in% map_sq$sci_name, ]
ages_sq1 <- ages_sq[ages_sq$Species %in% map_sq1$sci_name, ]
ages_sq2 <- ages_sq[ages_sq$Species %in% map_sq2$sci_name, ]
ages_av <- ages_av[ages_av$Species %in% map_av$sci_name, ]
ages_ma <- ages_ma[ages_ma$Species %in% map_ma$sci_name, ]

#save(map2_am, file = "data/map2_am.RData")
#save(map2_sq1, file = "data/map2_sq1.RData")
#save(map2_sq2, file = "data/map2_sq2.RData")
#save(map2_av, file = "data/map2_av.RData")
#save(map2_ma, file = "data/map2_ma.RData")

# Loading phylogenies
tr_sq <- read.nexus("~/Documents/lab/data/trees/amphibia_VertLife_27JUL20.nex")
for (i in 1:length(tr_sq)) tr_sq[[i]] <- drop.tip(tr_sq[[i]], "Homo_sapiens")
tr_sq <- read.nexus("~/Documents/lab/data/trees/squamata_VertLife_27JUL20.nex")
tr_av <- read.nexus("~/Documents/lab/data/trees/aves_Ericson_VertLife_27JUL20.nex")
tr_ma <- read.nexus("~/Documents/lab/data/trees/mammalia_node_dated_VertLife_27JUL20.nex")

# Set the directory where the unziped climatic files are
climFiles <- c("data/PALEO-PGEM-Series/monthly_mean/monthly_precipitation01M.txt", 
               "data/PALEO-PGEM-Series/monthly_mean/monthly_precipitation02M.txt", 
               "data/PALEO-PGEM-Series/monthly_mean/monthly_precipitation03M.txt", 
               "data/PALEO-PGEM-Series/monthly_mean/monthly_precipitation04M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_precipitation05M.txt", 
               "data/PALEO-PGEM-Series/monthly_mean/monthly_precipitation06M.txt", 
               "data/PALEO-PGEM-Series/monthly_mean/monthly_precipitation07M.txt", 
               "data/PALEO-PGEM-Series/monthly_mean/monthly_precipitation08M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_precipitation09M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_precipitation10M.txt", 
               "data/PALEO-PGEM-Series/monthly_mean/monthly_precipitation11M.txt", 
               "data/PALEO-PGEM-Series/monthly_mean/monthly_precipitation12M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_temperature01M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_temperature02M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_temperature03M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_temperature04M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_temperature05M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_temperature06M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_temperature07M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_temperature08M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_temperature09M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_temperature10M.txt",
               "data/PALEO-PGEM-Series/monthly_mean/monthly_temperature11M.txt",   
               "data/PALEO-PGEM-Series/monthly_mean/monthly_temperature12M.txt")  

# Function to generate bioclimatic variables for the PALEO-PGEM-Series dataset
generate_biovars <- function(time = NULL) {
  
  p_list = list()
  files_prec = climFiles[c(1:12)]
  files_temp = climFiles[c(13:24)]  
  column = 5003 - time
  temp_list = prec_list = list()
  
  for (i in 1:length(files_temp)) {
    if(i == 1) {
      coord = fread(files_temp[i], select= c (1:2))
    }  
    temp_list[[i]] = as.matrix(fread(files_temp[i], select = column))
    prec_list[[i]] = as.matrix(fread(files_prec[i], select = column))
  }
  
  for (j in 1:length(column)) {
    temp = do.call(cbind, lapply(temp_list, "[",,j))
    prec = do.call(cbind, lapply(prec_list, "[",,j))
    
    p <- matrix(nrow = nrow(prec), ncol = 2)
    colnames(p) = c("bio1", "bio12")
    window <- function(x)  { 
      lng <- length(x)
      x <- c(x,  x[1:3])
      m <- matrix(ncol = 3, nrow = lng)
      for (i in 1:3) { m[, i] <- x[i:(lng + i - 1)] }
      apply(m, MARGIN = 1, FUN = sum)
    }
    
    # Annual Mean Temperature 
    p[, 1] <- apply(temp, 1, mean)
    # Annual Precipitation 
    p[, 2] <- apply(prec, 1, sum)

    p = round(p, 3)
    
    p_list[[paste("T", time[j], sep = "_")]] = as.matrix(cbind(coord, p)) 
    print(paste("T", time[j], sep = "_"))
  }
  return(p_list)
}

times_of_interest <- c(1:21) 

biovars <- generate_biovars(time = times_of_interest)

for (i in 1:length(biovars)) {
  
  if (i < 10) {
    num <- paste0("0", strsplit(names(biovars)[i], "_")[[1]][2])
  } else {
    num <- strsplit(names(biovars)[i], "_")[[1]][2]
  }
  
  vars_T <- rasterFromXYZ(biovars[[i]][, 1:3])
  writeRaster(vars_T, 
              paste0("data/PALEO-PGEM-Series/tempfiles/", num, "_temp_T", num, ".asc"), 
              format = "ascii")
  vars_P <- rasterFromXYZ(biovars[[i]][, c(1, 2, 4)])
  writeRaster(vars_P, 
              paste0("data/PALEO-PGEM-Series/precfiles/", num, "_prec_T", num, ".asc"), 
              format = "ascii")
}

# Calculating temperature and precipitation stability
tempDev1 <- deviationThroughTime(variableDirectory = "data/PaleoView/tempfiles",
                                 timeSlicePeriod = 1000, fileExtension = "asc")
tempStability1 <- stabilityCalc(tempDev1)

tempDev2 <- deviationThroughTime(variableDirectory = "data/PALEO-PGEM-Series/tempfiles",
                                 timeSlicePeriod = 1000, fileExtension = "asc")
tempStability2 <- stabilityCalc(tempDev2)

precDev1 <- deviationThroughTime(variableDirectory = "data/PaleoView/precfiles",
                                 timeSlicePeriod = 1000, fileExtension = "asc")
precStability1 <- stabilityCalc(precDev1)

precDev2 <- deviationThroughTime(variableDirectory = "data/PALEO-PGEM-Series/precfiles",
                                 timeSlicePeriod = 1000, fileExtension = "asc")
precStability2 <- stabilityCalc(precDev2)

# Getting averaged climatic variables in the ranges of each species and 
# creating rasters with values indicating species ages

r1 <- raster(ncols = 144, nrows = 72, ymn = -90, vals = 1:(144*72))
r2 <- raster(ncols = 360, nrows = 149, ymn = -65, vals = 1:(360*149))

temp_vals1 <- raster::extract(tempStability1, 1:length(r1[]))
temp_vals2 <- raster::extract(tempStability2, 1:length(r2[]))

prec_vals1 <- raster::extract(precStability1, 1:length(r1[]))
prec_vals2 <- raster::extract(precStability2, 1:length(r2[]))

getVarSpp <- function(maps, ages, path) {
  
  amp <- c("Anura", "Caudata", "Gymnophiona")
  squ <- c("Anguimorpha", "Gekkota", "Iguania", "Lacertoidea", "Scincoidea", 
           "Serpentes")
  ave <- c("Columbiformes", "Passeriformes", "Piciformes", "Psittaciformes")
  mam <- c("Carnivora", "Cetartiodactyla", "Chiroptera", "Diprotodontia", 
           "Primates")
  
  taxon <- ifelse(all(names(table(ages$Taxon)) %in% amp), "Amphibia",
                  ifelse(all(names(table(ages$Taxon)) %in% squ), "Squamata",
                         ifelse(all(names(table(ages$Taxon)) %in% ave),
                                "Aves", "Mammalia")))
  
  date <- paste(strsplit(format(Sys.time(), "%e %b %Y %H-%M"), " ")[[1]][c(2:5)], 
                collapse = "_")
  
  for (i in 1:nrow(ages)) {
    
    print(paste0(i, "/", nrow(ages)))
    
    species_i <- ages$Species[i]
    
    maps_i <- subset(maps, maps$sci_name == species_i) 
    
    if (taxon == "Aves") {
      for (j in 1:length(maps_i$Shape)) { 
        try(maps_i$Shape[[j]] <- st_cast(maps_i$Shape[[j]], 'MULTIPOLYGON'))
      }
      try(maps_i$Shape <- st_cast(maps_i$Shape, 'MULTIPOLYGON'))
      
      raster1_i <- fasterize(st_as_sf(maps_i$Shape), r1)
      raster2_i <- fasterize(st_as_sf(maps_i$Shape), r2)
    } else {
      raster1_i <- rasterize(maps_i, r1)
      raster2_i <- rasterize(maps_i, r2)
    }
  
    rastercells1 <- which(getValues(raster1_i) > 0)
    rastercells2 <- which(getValues(raster2_i) > 0)
    
    values_temp1 <- temp_vals1[rastercells1, ]
    values_temp2 <- temp_vals2[rastercells2, ]

    values_prec1 <- prec_vals1[rastercells1, ]
    values_prec2 <- prec_vals2[rastercells2, ]
    
    if (length(values_temp1) == 1) {
      ifelse(is.na(values_temp1[1]),
             temp_spp1 <- NA,
             temp_spp1 <- values_temp1[1])
    } else { 
      ifelse(all(is.na(values_temp1)),
             temp_spp1 <- NA,
             temp_spp1 <- mean(values_temp1, na.rm = TRUE))
    }
    
    if (length(values_temp2) == 1) {
      ifelse(is.na(values_temp2[1]),
             temp_spp2 <- NA,
             temp_spp2 <- values_temp2[1])
    } else { 
      ifelse(all(is.na(values_temp2)),
             temp_spp2 <- NA,
             temp_spp2 <- mean(values_temp2, na.rm = TRUE))
    }
    
    if (length(values_prec1) == 1) {
      ifelse(is.na(values_prec1[1]),
             prec_spp1 <- NA,
             prec_spp1 <- values_prec1[1])
    } else { 
      ifelse(all(is.na(values_prec1)),
             prec_spp1 <- NA,
             prec_spp1 <- mean(values_prec1, na.rm = TRUE))
    }
    
    if (length(values_prec2) == 1) {
      ifelse(is.na(values_prec2[1]),
             prec_spp2 <- NA,
             prec_spp2 <- values_prec2[1])
    } else { 
      ifelse(all(is.na(values_prec2)),
             prec_spp2 <- NA,
             prec_spp2 <- mean(values_prec2, na.rm = TRUE))
    }
    
    study_results <- data.frame(species = maps_i$sci_name[1],
                                ages$Ages[ages$Species == species_i][[1]],
                                temp_spp1,
                                temp_spp2,
                                prec_spp1,
                                prec_spp2)
    
    write.table(study_results, paste0("data/ClimSpp_", taxon, "_", date, ".csv"),
                sep = ",", col.names = FALSE, append = TRUE, 
                row.names = F)
    
    rastervals1 <- which(getValues(!is.na(raster1_i)))
    rastervals2 <- which(getValues(!is.na(raster2_i)))
    
    raster1_i[rastervals1] <- ages$Ages[i]
    raster2_i[rastervals2] <- ages$Ages[i]
    
    writeRaster(raster1_i, paste0(path, "/raster1_", species_i, ".tif"), format = "GTiff")
    writeRaster(raster2_i, paste0(path, "/raster2_", species_i, ".tif"), format = "GTiff")
  }
  
  res <- read.csv(paste0("data/ClimSpp_", taxon, "_", date, ".csv"),
                  row.names = 1, header = F)
  colnames(res) <- c("Ages", "Temp_PaleoView", "Temp_PALEOPGEM", 
                     "Prec_PaleoView", "Prec_PALEOPGEM")
  write.csv(res, paste0("data/ClimSpp_", taxon, "_", date, ".csv"))
  return(res)
  
}

env_sq <- getVarSpp(map_sq, ages_sq, "/Volumes/Personal/lab/spp_ages/amphibia")
env_sq1 <- getVarSpp(map_sq1, ages_sq1, "/Volumes/Personal/lab/spp_ages/squamata")
env_sq2 <- getVarSpp(map_sq2, ages_sq2, "/Volumes/Personal/lab/spp_ages/squamata")
env_av <- getVarSpp(map_av, ages_av, "/Volumes/Personal/lab/spp_ages/aves")
env_ma <- getVarSpp(map_ma, ages_ma, "/Volumes/Personal/lab/spp_ages/mammalia")

#env_sq <- read.csv("data/ClimSpp_Amphibia_6_Apr_2023_15-27.csv", row.names = 1)
#env_sq1 <- read.csv("data/ClimSpp_Squamata_6_Apr_2023_13-21.csv", row.names = 1)
#env_sq2 <- read.csv("data/ClimSpp_Squamata_6_Apr_2023_14-17.csv", row.names = 1)
#env_av <- read.csv("data/ClimSpp_Aves_Apr_2023_19-46_NA.csv", row.names = 1)
#env_ma <- read.csv("data/ClimSpp_Mammalia_7_Apr_2023_10-42.csv", row.names = 1)

env_sq <- rbind(env_sq1, env_sq2)

comp_sq <- lapply(tr_sq, comparative.data, data = cbind(env_sq, rownames(env_sq)), names.col = "rownames(env_sq)")
comp_sq <- lapply(tr_sq, comparative.data, data = cbind(env_sq, rownames(env_sq)), names.col = "rownames(env_sq)") 
comp_av <- lapply(tr_av, comparative.data, data = cbind(env_av, rownames(env_av)), names.col = "rownames(env_av)") 
comp_ma <- lapply(tr_ma, comparative.data, data = cbind(env_ma, rownames(env_ma)), names.col = "rownames(env_ma)") 

# Calculating PGLS

n <- 100

for (i in 1:n) {
  print(i)
  pgls1_int_am <- pgls(Ages ~ Temp_PaleoView + Prec_PaleoView + Temp_PaleoView*Prec_PaleoView, data = comp_am[[i]])
  pgls2_int_am <- pgls(Ages ~ Temp_PALEOPGEM + Prec_PALEOPGEM + Temp_PALEOPGEM*Prec_PALEOPGEM, data = comp_am[[i]])

  summ1_am <- summary(pgls1_int_am)
  summ2_am <- summary(pgls2_int_am)
  res_int_am <- data.frame(Intercept1 = round(summ1_am$coefficients[1, 1], 3),
                           Int1_SD = round(summ1_am$coefficients[1, 2], 3),
                           Int1_t = round(summ1_am$coefficients[1, 3], 3),
                           Int1_p = round(summ1_am$coefficients[1, 4], 3),
                           Temp1 = round(summ1_am$coefficients[2, 1], 3),
                           Temp1_SD = round(summ1_am$coefficients[2, 2], 3),
                           Temp1_t = round(summ1_am$coefficients[2, 3], 3), 
                           Temp1_p = round(summ1_am$coefficients[2, 4], 3),
                           Prec1 = round(summ1_am$coefficients[3, 1], 3),
                           Prec1_SD = round(summ1_am$coefficients[3, 2], 3),
                           Prec1_t = round(summ1_am$coefficients[3, 3], 3), 
                           Prec1_p = round(summ1_am$coefficients[3, 4], 3),
                           TempPrec1 = round(summ1_am$coefficients[4, 1], 3),
                           TempPrec1_SD = round(summ1_am$coefficients[4, 2], 3),
                           TempPrec1_t = round(summ1_am$coefficients[4, 3], 3),
                           TempPrec1_p = round(summ1_am$coefficients[4, 4], 3),
                           rSquared1 = round(summ1_am$r.squared, 3),
                           Intercept2 = round(summ2_am$coefficients[1, 1], 3),
                           Int2_SD = round(summ2_am$coefficients[1, 2], 3),
                           Int2_t = round(summ2_am$coefficients[1, 3], 3),
                           Int2_p = round(summ2_am$coefficients[1, 4], 3),
                           Temp2 = round(summ2_am$coefficients[2, 1], 3),
                           Temp2_SD = round(summ2_am$coefficients[2, 2], 3),
                           Temp2_t = round(summ2_am$coefficients[2, 3], 3), 
                           Temp2_p = round(summ2_am$coefficients[2, 4], 3),
                           Prec2 = round(summ2_am$coefficients[3, 1], 3),
                           Prec2_SD = round(summ2_am$coefficients[3, 2], 3),
                           Prec2_t = round(summ2_am$coefficients[3, 3], 3), 
                           Prec2_p = round(summ2_am$coefficients[3, 4], 3),
                           TempPrec2 = round(summ2_am$coefficients[4, 1], 3),
                           TempPrec2_SD = round(summ2_am$coefficients[4, 2], 3),
                           TempPrec2_t = round(summ2_am$coefficients[4, 3], 3),
                           TempPrec2_p = round(summ2_am$coefficients[4, 4], 3),
                           rSquared2 = round(summ2_am$r.squared, 3))
  
  if (i == 1) {
    write.table(res_int_am, "tables/Tables1,S1_unformatted_am.csv", sep = ",",
                col.names = T, row.names = F)
  } else {
    write.table(res_int_am, "tables/Tables1,S1_unformatted_am.csv", sep = ",",
                col.names = F, append = TRUE, row.names = F)
  }
}


for (i in 1:n) {
  print(i)
  pgls1_int_sq <- pgls(Ages ~ Temp_PaleoView + Prec_PaleoView + Temp_PaleoView*Prec_PaleoView, data = comp_sq[[i]])
  pgls2_int_sq <- pgls(Ages ~ Temp_PALEOPGEM + Prec_PALEOPGEM + Temp_PALEOPGEM*Prec_PALEOPGEM, data = comp_sq[[i]])

  summ1_sq <- summary(pgls1_int_sq)
  summ2_sq <- summary(pgls2_int_sq)
  res_int_sq <- data.frame(Intercept1 = round(summ1_sq$coefficients[1, 1], 3),
                           Int1_SD = round(summ1_sq$coefficients[1, 2], 3),
                           Int1_t = round(summ1_sq$coefficients[1, 3], 3),
                           Int1_p = round(summ1_sq$coefficients[1, 4], 3),
                           Temp1 = round(summ1_sq$coefficients[2, 1], 3),
                           Temp1_SD = round(summ1_sq$coefficients[2, 2], 3),
                           Temp1_t = round(summ1_sq$coefficients[2, 3], 3), 
                           Temp1_p = round(summ1_sq$coefficients[2, 4], 3),
                           Prec1 = round(summ1_sq$coefficients[3, 1], 3),
                           Prec1_SD = round(summ1_sq$coefficients[3, 2], 3),
                           Prec1_t = round(summ1_sq$coefficients[3, 3], 3), 
                           Prec1_p = round(summ1_sq$coefficients[3, 4], 3),
                           TempPrec1 = round(summ1_sq$coefficients[4, 1], 3),
                           TempPrec1_SD = round(summ1_sq$coefficients[4, 2], 3),
                           TempPrec1_t = round(summ1_sq$coefficients[4, 3], 3),
                           TempPrec1_p = round(summ1_sq$coefficients[4, 4], 3),
                           rSquared1 = round(summ1_sq$r.squared, 3),
                           Intercept2 = round(summ2_sq$coefficients[1, 1], 3),
                           Int2_SD = round(summ2_sq$coefficients[1, 2], 3),
                           Int2_t = round(summ2_sq$coefficients[1, 3], 3),
                           Int2_p = round(summ2_sq$coefficients[1, 4], 3),
                           Temp2 = round(summ2_sq$coefficients[2, 1], 3),
                           Temp2_SD = round(summ2_sq$coefficients[2, 2], 3),
                           Temp2_t = round(summ2_sq$coefficients[2, 3], 3), 
                           Temp2_p = round(summ2_sq$coefficients[2, 4], 3),
                           Prec2 = round(summ2_sq$coefficients[3, 1], 3),
                           Prec2_SD = round(summ2_sq$coefficients[3, 2], 3),
                           Prec2_t = round(summ2_sq$coefficients[3, 3], 3), 
                           Prec2_p = round(summ2_sq$coefficients[3, 4], 3),
                           TempPrec2 = round(summ2_sq$coefficients[4, 1], 3),
                           TempPrec2_SD = round(summ2_sq$coefficients[4, 2], 3),
                           TempPrec2_t = round(summ2_sq$coefficients[4, 3], 3),
                           TempPrec2_p = round(summ2_sq$coefficients[4, 4], 3),
                           rSquared2 = round(summ2_sq$r.squared, 3))
  
  if (i == 1) {
    write.table(res_int_sq, "tables/Tables1,S1_unformatted_sq.csv", sep = ",",
                col.names = T, row.names = F)
  } else {
    write.table(res_int_sq, "tables/Tables1,S1_unformatted_sq.csv", sep = ",",
                col.names = F, append = TRUE, row.names = F)
  }
}

for (i in 1:n) {
  print(i)
  pgls1_int_av <- pgls(Ages ~ Temp_PaleoView + Prec_PaleoView + Temp_PaleoView*Prec_PaleoView, data = comp_av[[i]])
  pgls2_int_av <- pgls(Ages ~ Temp_PALEOPGEM + Prec_PALEOPGEM + Temp_PALEOPGEM*Prec_PALEOPGEM, data = comp_av[[i]])

  summ1_av <- summary(pgls1_int_av)
  summ2_av <- summary(pgls2_int_av)
  res_int_av <- data.frame(Intercept1 = round(summ1_av$coefficients[1, 1], 3),
                           Int1_SD = round(summ1_av$coefficients[1, 2], 3),
                           Int1_t = round(summ1_av$coefficients[1, 3], 3),
                           Int1_p = round(summ1_av$coefficients[1, 4], 3),
                           Temp1 = round(summ1_av$coefficients[2, 1], 3),
                           Temp1_SD = round(summ1_av$coefficients[2, 2], 3),
                           Temp1_t = round(summ1_av$coefficients[2, 3], 3), 
                           Temp1_p = round(summ1_av$coefficients[2, 4], 3),
                           Prec1 = round(summ1_av$coefficients[3, 1], 3),
                           Prec1_SD = round(summ1_av$coefficients[3, 2], 3),
                           Prec1_t = round(summ1_av$coefficients[3, 3], 3), 
                           Prec1_p = round(summ1_av$coefficients[3, 4], 3),
                           TempPrec1 = round(summ1_av$coefficients[4, 1], 3),
                           TempPrec1_SD = round(summ1_av$coefficients[4, 2], 3),
                           TempPrec1_t = round(summ1_av$coefficients[4, 3], 3),
                           TempPrec1_p = round(summ1_av$coefficients[4, 4], 3),
                           rSquared1 = round(summ1_av$r.squared, 3),
                           Intercept2 = round(summ2_av$coefficients[1, 1], 3),
                           Int2_SD = round(summ2_av$coefficients[1, 2], 3),
                           Int2_t = round(summ2_av$coefficients[1, 3], 3),
                           Int2_p = round(summ2_av$coefficients[1, 4], 3),
                           Temp2 = round(summ2_av$coefficients[2, 1], 3),
                           Temp2_SD = round(summ2_av$coefficients[2, 2], 3),
                           Temp2_t = round(summ2_av$coefficients[2, 3], 3), 
                           Temp2_p = round(summ2_av$coefficients[2, 4], 3),
                           Prec2 = round(summ2_av$coefficients[3, 1], 3),
                           Prec2_SD = round(summ2_av$coefficients[3, 2], 3),
                           Prec2_t = round(summ2_av$coefficients[3, 3], 3), 
                           Prec2_p = round(summ2_av$coefficients[3, 4], 3),
                           TempPrec2 = round(summ2_av$coefficients[4, 1], 3),
                           TempPrec2_SD = round(summ2_av$coefficients[4, 2], 3),
                           TempPrec2_t = round(summ2_av$coefficients[4, 3], 3),
                           TempPrec2_p = round(summ2_av$coefficients[4, 4], 3),
                           rSquared2 = round(summ2_av$r.squared, 3))
  
  if (i == 1) {
    write.table(res_int_av, "tables/Tables1,S1_unformatted_av.csv", sep = ",",
                col.names = T, row.names = F)
  } else {
    write.table(res_int_av, "tables/Tables1,S1_unformatted_av.csv", sep = ",",
                col.names = F, append = TRUE, row.names = F)
  }
}

for (i in 1:n) {
  print(i)
  pgls1_int_ma <- pgls(Ages ~ Temp_PaleoView + Prec_PaleoView + Temp_PaleoView*Prec_PaleoView, data = comp_ma[[i]])
  pgls2_int_ma <- pgls(Ages ~ Temp_PALEOPGEM + Prec_PALEOPGEM + Temp_PALEOPGEM*Prec_PALEOPGEM, data = comp_ma[[i]])

  summ1_ma <- summary(pgls1_int_ma)
  summ2_ma <- summary(pgls2_int_ma)
  res_int_ma <- data.frame(Intercept1 = round(summ1_ma$coefficients[1, 1], 3),
                           Int1_SD = round(summ1_ma$coefficients[1, 2], 3),
                           Int1_t = round(summ1_ma$coefficients[1, 3], 3),
                           Int1_p = round(summ1_ma$coefficients[1, 4], 3),
                           Temp1 = round(summ1_ma$coefficients[2, 1], 3),
                           Temp1_SD = round(summ1_ma$coefficients[2, 2], 3),
                           Temp1_t = round(summ1_ma$coefficients[2, 3], 3), 
                           Temp1_p = round(summ1_ma$coefficients[2, 4], 3),
                           Prec1 = round(summ1_ma$coefficients[3, 1], 3),
                           Prec1_SD = round(summ1_ma$coefficients[3, 2], 3),
                           Prec1_t = round(summ1_ma$coefficients[3, 3], 3), 
                           Prec1_p = round(summ1_ma$coefficients[3, 4], 3),
                           TempPrec1 = round(summ1_ma$coefficients[4, 1], 3),
                           TempPrec1_SD = round(summ1_ma$coefficients[4, 2], 3),
                           TempPrec1_t = round(summ1_ma$coefficients[4, 3], 3),
                           TempPrec1_p = round(summ1_ma$coefficients[4, 4], 3),
                           rSquared1 = round(summ1_ma$r.squared, 3),
                           Intercept2 = round(summ2_ma$coefficients[1, 1], 3),
                           Int2_SD = round(summ2_ma$coefficients[1, 2], 3),
                           Int2_t = round(summ2_ma$coefficients[1, 3], 3),
                           Int2_p = round(summ2_ma$coefficients[1, 4], 3),
                           Temp2 = round(summ2_ma$coefficients[2, 1], 3),
                           Temp2_SD = round(summ2_ma$coefficients[2, 2], 3),
                           Temp2_t = round(summ2_ma$coefficients[2, 3], 3), 
                           Temp2_p = round(summ2_ma$coefficients[2, 4], 3),
                           Prec2 = round(summ2_ma$coefficients[3, 1], 3),
                           Prec2_SD = round(summ2_ma$coefficients[3, 2], 3),
                           Prec2_t = round(summ2_ma$coefficients[3, 3], 3), 
                           Prec2_p = round(summ2_ma$coefficients[3, 4], 3),
                           TempPrec2 = round(summ2_ma$coefficients[4, 1], 3),
                           TempPrec2_SD = round(summ2_ma$coefficients[4, 2], 3),
                           TempPrec2_t = round(summ2_ma$coefficients[4, 3], 3),
                           TempPrec2_p = round(summ2_ma$coefficients[4, 4], 3),
                           rSquared2 = round(summ2_ma$r.squared, 3))
  
  if (i == 1) {
    write.table(res_int_ma, "tables/Tables1,S1_unformatted_ma.csv", sep = ",",
                col.names = T, row.names = F)
  } else {
    write.table(res_int_ma, "tables/Tables1,S1_unformatted_ma.csv", sep = ",",
                col.names = F, append = TRUE, row.names = F)
  }
}

stats <- function(x) {
    mean <- round(mean(x), 3)
    min <- round(range(x)[1], 3)
    max <- round(range(x)[2], 3)

    res <- paste0(mean, " (", min, "-", max, ")")
    res
}

res_int_am <- read.csv("tables/Tables1,S1_unformatted_am.csv")[1:100,]
res_int_sq <- read.csv("tables/Tables1,S1_unformatted_sq.csv")[1:100,]
res_int_av <- read.csv("tables/Tables1,S1_unformatted_av.csv")[1:100,]
res_int_av[, 34] <- as.numeric(res_int_av[, 34])
res_int_ma <- read.csv("tables/Tables1,S1_unformatted_ma.csv")[1:100,]

res_pgls <- as.data.frame(matrix(ncol = 34, nrow = 4))
colnames(res_pgls) <- colnames(res_int_sq)
rownames(res_pgls) <- c("Amphibia", "Squamata", "Aves", "Mammalia")

res_pgls[1, ] <- apply(res_int_am, 2, stats)
res_pgls[2, ] <- apply(res_int_sq, 2, stats)
res_pgls[3, ] <- apply(res_int_av, 2, stats)
res_pgls[4, ] <- apply(res_int_ma, 2, stats)

write.csv(res_pgls, "tables/Tables1,S1_unformatted.csv")

###### Processing the rasters on QGIS to get average values of ages per cell

# Function to convert raster to be plotted by ggplot2
get_values <- function(full) {

  data(wrld_simpl, package = "maptools")

  rem1 <- extract(full, wrld_simpl, cellnumbers = T, weights = T, 
                  small = T)
  rem2 <- do.call(rbind.data.frame, rem1)[, 1]
  values(full)[-rem2] <- NA

  full_p <- rasterToPoints(full, spatial = TRUE)
  full_df  <- data.frame(full_p)
  colnames(full_df)[1] <- "index_1"
  #full_df <- full_df %>% mutate(index_2 = cut(index_1, breaks = brks))

  return(full_df)
} 

raster1_am <- get_values(raster("data/raster1_am.tif"))
raster2_am <- get_values(raster("data/raster2_am.tif"))

raster1_sq <- get_values(raster("data/raster1_sq.tif"))
raster2_sq <- get_values(raster("data/raster2_sq.tif"))

raster1_av <- get_values(raster("data/raster1_av.tif"))
raster2_av <- get_values(raster("data/raster2_av.tif"))

raster1_ma <- get_values(raster("data/raster1_ma.tif"))
raster2_ma <- get_values(raster("data/raster2_ma.tif"))

raster_temp1 <- get_values(raster(tempStability1))
raster_temp2 <- get_values(raster(tempStability2))

raster_prec1 <- get_values(raster(precStability1))
raster_prec2 <- get_values(raster(precStability2))

paltemp <- c("#FFC971", "#FFB627", "#FF9505", "#E2711D", "#CC5803")
palprec <- c("#A9D6E5", "#89C2D9", "#468FAF", "#2A6F97", "#01497C", "#012A4A")
palages <- c("#CAD2C5", "#84A98C", "#52796F", "#354F52", "#2F3E46")

r <- raster(ncols = 144, nrows = 72, ymn = -90)

ras_wrld1 <- rasterToPoints(rasterize(wrld_simpl, r), spatial = TRUE)
ras_wrld1  <- data.frame(ras_wrld1)

temp1 <- ggplot() + 
  geom_raster(data = ras_wrld1, aes(x = x, y = y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_temp1, aes(x = x, y = y, fill = index_1)) +
  scale_fill_gradientn(colours = paltemp, name = "Temperature\nstability") +
  theme_void() 

prec1 <- ggplot() + 
  geom_raster(data = ras_wrld1, aes(x = x, y = y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_prec1, aes(x = x, y = y, fill = index_1)) +
  scale_fill_gradientn(colours = palprec, name = "Precipitation\nstability") +
  theme_void()

plot1_am <- ggplot() +
  geom_raster(data = ras_wrld1, aes(x = x, y = y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster1_am, aes(x = x, y = y, fill = index_1)) +
  scale_fill_gradientn(colours = palages, name = "Ages") +
  theme_void()

plot1_sq <- ggplot() +
  geom_raster(data = ras_wrld1, aes(x = x, y = y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster1_sq, aes(x = x, y = y, fill = index_1)) +
  scale_fill_gradientn(colours = palages, name = "Ages") +
  theme_void() 


plot1_av <- ggplot() +
  geom_raster(data = ras_wrld1, aes(x = x, y = y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster1_av, aes(x = x, y = y, fill = index_1)) +
  scale_fill_gradientn(colours = palages, name = "Ages") +
  theme_void() 


plot1_ma <- ggplot() +
  geom_raster(data = ras_wrld1, aes(x = x, y = y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster1_ma, aes(x = x, y = y, fill = index_1)) +
  scale_fill_gradientn(colours = palages, name = "Ages") +
  theme_void() 

r_2 <- raster(ncols = 360, nrows = 149, ymn = -65)

ras_wrld2 <- rasterToPoints(rasterize(wrld_simpl, r_2), spatial = TRUE)
ras_wrld2  <- data.frame(ras_wrld2)

temp2 <- ggplot() + 
  geom_raster(data = ras_wrld2, aes(x = x, y = y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_temp2, aes(x = x, y = y, fill = index_1)) +
  scale_fill_gradientn(colours = paltemp, name = "Temperature\nstability") +
  theme_void() 

prec2 <- ggplot() + 
  geom_raster(data = ras_wrld2, aes(x = x, y = y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster_prec2, aes(x = x, y = y, fill = index_1)) +
  scale_fill_gradientn(colours = palprec, name = "Precipitation\nstability") +
  theme_void()

plot2_am <- ggplot() +
  geom_raster(data = ras_wrld2, aes(x = x, y = y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster2_am, aes(x = x, y = y, fill = index_1)) +
  scale_fill_gradientn(colours = palages, name = "Ages") +
  theme_void() 

plot2_sq <- ggplot() +
  geom_raster(data = ras_wrld2, aes(x = x, y = y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster2_sq, aes(x = x, y = y, fill = index_1)) +
  scale_fill_gradientn(colours = palages, name = "Ages") +
  theme_void() 


plot2_av <- ggplot() +
  geom_raster(data = ras_wrld2, aes(x = x, y = y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster2_av, aes(x = x, y = y, fill = index_1)) +
  scale_fill_gradientn(colours = palages, name = "Ages") +
  theme_void() 


plot2_ma <- ggplot() +
  geom_raster(data = ras_wrld2, aes(x = x, y = y), fill = "gray", show.legend = FALSE) +
  geom_raster(data = raster2_ma, aes(x = x, y = y, fill = index_1)) +
  scale_fill_gradientn(colours = palages, name = "Ages") +
  theme_void() 

pdf("figures/Figure3.pdf", width = 11)
(temp1 + prec1) / (plot1_am + plot1_sq) / (plot1_av + plot1_ma)
dev.off()

pdf("figures/FigureS1.pdf", width = 11)
(temp2 + prec2) / (plot2_am + plot2_sq) / (plot2_av + plot2_ma)
dev.off()
