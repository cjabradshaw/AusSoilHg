# soil Hg analysis
# determinants of Hg in soil for Australia
# Corey Bradshaw
# July 2025

# rm(list = ls()) # remove all objects from R environment

# increase memory to max
mem.maxVSize(v = Inf)

# libraries
library(ade4) # analyse ecological/environmental data
library(adegraphics) # visualisation of ade4 objects
library(adespatial) 
library(blockCV) # blocking of spatially structured data
library(boot) # bootstrapping
library(caret) # for model training and prediction
library(data.table)
library(dismo) # distribution modelling
library(dplyr) # data manipulation
library(fields) # distance matrix utilities
library(gbm) # boosted regression trees
library(geodata) # downloading spatial data
library(geosphere) # geographic distance calculations
library(ggplot2) # plotting
library(gstat) # for geostatistical modelling
library(gridExtra) # for grid graphics
library(kableExtra) # tabulation
library(ncdf4) # opening NetCDF files
library(patchwork) # multiplot layouts
library(performance) # model diagnostics
library(pdp) # partial dependence plots
library(randomForest) # random forest modelling
library(randomForestExplainer)
library(raster) # raster data handling
library(rnaturalearth) # natural Earth map data
library(rnaturalearthdata) # vector map data
library(sf) # simple features objects
library(sjPlot) # model diagnostics
library(sp) # methods for spatial data
library(spatialRF) # spatial regression with random forest
library(spatstat) # analysis of spatial point patterns
library(spdep) # package version & build info
library(terra) # spatial data manipulation & visualisation
library(tidyverse) # data management
library(usdm) # variable inflation
library(viridis) # for colour palettes

# import datasets
## geochemical dataset
geochem.dat <- read.csv("geochem.csv", header = T)
head(geochem.dat)
dim(geochem.dat)

geochem.TOS <- subset(geochem.dat, DEPTH == "TOS")
dim(geochem.TOS)
geochem.BOS <- subset(geochem.dat, DEPTH == "BOS")
dim(geochem.BOS)

# sort
geochem.sort <- geochem.dat[order(geochem.dat$SITEID),]
head(geochem.sort)

## field dataset
field.dat <- read.csv("field.csv", header = T)
head(field.dat)
dim(field.dat)

# sort
field.sort <- field.dat[order(field.dat$SITEID), ]
head(field.sort)

## Hg dataset
hg.dat <- read.csv("hg.csv", header = T)
head(hg.dat)
dim(hg.dat)

## add TARGSITEID from field.dat to hg.dat
head(field.dat)
head(hg.dat)
fieldIDs <- field.dat[, c("SITEID", "TARGSITEID")]
head(fieldIDs)
hg.dat2 <- merge(hg.dat, fieldIDs, by="SITEID", all.y=F, all.x=F)
head(hg.dat2)
hg.dat3 <- hg.dat2 %>% select(SITEID, TARGSITEID, SAMPLEID, HgCOMP)
head(hg.dat3)
dim(hg.dat3)

## add AREA from area.csv to hg3.dat & export
area.dat <- read.csv("areas.csv", header = T)
head(area.dat)
hg.dat4 <- merge(hg.dat3, area.dat, by="TARGSITEID", all.y=F, all.x=T)
head(hg.dat4)
dim(hg.dat4)
hg.dat4$SAMPLEID.MRG <- substr(hg.dat4$SAMPLEID, 1, 13)  # create proper merge column
head(hg.dat4)

which(is.na(hg.dat4$AREA_SQKM))
write.csv(hg.dat4, "hgTSID.csv", row.names=F)

areaXTARGSITID.stats <- hg.dat4 %>%
  group_by(TARGSITEID) %>%
  summarise(
    mean = mean(AREA_SQKM, na.rm = TRUE),
    n = n()
  )
areaXTARGSITID.stats

sum(areaXTARGSITID.stats$mean, na.rm = T)
area.aus <- 7688287
# proportion of catchments in Australia covered by Hg samples
sum(areaXTARGSITID.stats$mean, na.rm = T)/area.aus

## grain size dataset
gs.dat <- read.csv("gs.csv", header = T)
head(gs.dat)
dim(gs.dat)

# subset bulk rows only
gs.bulk <- subset(gs.dat, GS == "bulk")
dim(gs.bulk)

# grain size categories
gs.cat <- read.csv("gscats.csv", header = T)
head(gs.cat)

# create weighted mean grain size column
gs.bulk$gs.wmn <- NA
for (i in 1:dim(gs.bulk)[1]) {
  gs.bulk$gs.wmn[i] <- sum(gs.cat$mn * as.numeric(gs.bulk[i, 10:(dim(gs.bulk)[2]-1)])/100, na.rm=T)
}
gs.bulk$gs.wmn
head(gs.bulk[1:3,])

# merge datasets
head(geochem.dat)
geochem.sort$SAMPLEID.MRG <- substr(geochem.sort$SAMPLEID, 1, 13)  # create proper merge column
head(geochem.sort)
str(geochem.sort)

head(field.dat)
str(field.dat)
head(field.dat[,-c(20:25)])
BOSfield.dat <- field.dat[,-c(20:25)]
head(BOSfield.dat)
BOSfield.dat$SAMPLEID.MRG <- as.character(BOSfield.dat$BOSID)
head(BOSfield.dat)
dim(BOSfield.dat)

head(field.dat[,-c(26:dim(field.dat)[2])])
TOSfield.dat <- field.dat[,-c(26:dim(field.dat)[2])]
head(TOSfield.dat)
TOSfield.dat$SAMPLEID.MRG <- as.character(TOSfield.dat$TOSID)
head(TOSfield.dat)
dim(TOSfield.dat)

# merge geochem and field datasets
datmrg1 <- merge(geochem.sort, TOSfield.dat, by="SAMPLEID.MRG", no.dups=T, sort=T, all.x=F, all.y=T)
head(datmrg1)
dim(datmrg1)
colnames(datmrg1)[95:100] <- c("ID", "TYPE", "TOPD", "BD", "RAD", "pH")
head(datmrg1)
dim(datmrg1)
table(datmrg1$DEPTH)

datmrg2 <- merge(geochem.sort, BOSfield.dat, by="SAMPLEID.MRG", no.dups=T, sort=T, all.x=F, all.y=T)
head(datmrg2)
colnames(datmrg2)[95:100] <- c("ID", "TYPE", "TOPD", "BD", "RAD", "pH")
head(datmrg2)
dim(datmrg2)
table(datmrg2$DEPTH)
dim(datmrg1)[1] + dim(datmrg2)[1]

datmrg3 <- rbind(datmrg1, datmrg2)
head(datmrg3)
dim(datmrg3)

datmrg3.sort <- datmrg3[order(datmrg3$SITEID.x), ]
head(datmrg3.sort)
dim(datmrg3.sort)

head(hg.dat4)
datmrg4 <- merge(hg.dat4, datmrg3.sort, by="SAMPLEID.MRG", no.dups=T, sort=T, all.x=F, all.y=T)
head(datmrg4)
datmrg4.sort <- datmrg4[order(datmrg4$SITEID), ]
head(datmrg4.sort)
dim(datmrg4.sort)

head(gs.bulk)
gs.bulk$SAMPLEID.MRG <- substr(gs.bulk$SAMPLEID, 1, 13)
head(gs.bulk)
dim(gs.bulk)

datmrg5 <- merge(gs.bulk, datmrg4.sort, by="SAMPLEID.MRG", no.dups=T, sort=T, all.x=F, all.y=T)
head(datmrg5)
dim(datmrg5)

# remove unnecessary columns
colnames(datmrg5)
datmrg6 <- datmrg5[, !names(datmrg5) %in% c("SITEID.y", "SAMPLEID.x", "LON.y", "TARGSITEID.y",
                                            "DEPTH.y", "DB.y", "DUPL.y", "LAT.y", "DATE.y", "STATE.y",
                                            "SAMPLEID.y","GS", "TIME", "ID", "RAD", "pH")]
which(colnames(datmrg6) == "SITEID.x.1")
datmrg6 <- datmrg6[, -(which(colnames(datmrg6) == "SITEID.x.1"))]
colnames(datmrg6)
head(datmrg6)

setnames(datmrg6, old=c("SITEID.x", "DEPTH.x", "DUPL.x", "LAT.x", "DATE.x", "LON.x", "STATE.x", "TARGSITEID.x"),
         new = c("SITEID", "DEPTH", "DUPL", "LAT", "DATE", "LON", "STATE", "SAMPLEID"))
head(datmrg6)
datmrg7.sort <- datmrg6[order(datmrg6$SITEID), ]
head(datmrg7.sort)
dim(datmrg7.sort)

# remove duplicate columns
duplicate_names <- duplicated(names(datmrg7.sort))
names(datmrg7.sort)[duplicate_names]  # duplicate column names
datmrg7.clean <- datmrg7.sort[, !duplicated(t(datmrg7.sort))]
dim(datmrg7.clean)
head(datmrg7.clean)

# how many HgCOMP == NA?
length(which(is.na(datmrg7.clean$HgCOMP) == T))
datmrg8.sort <- datmrg7.clean[-which(is.na(datmrg7.clean$HgCOMP) == T), ]
head(datmrg8.sort)
dim(datmrg8.sort)

write.csv(datmrg8.sort, "datmrg.csv", row.names=F)

datmrg <- datmrg8.sort
str(datmrg)

# redo landuse recategorisation
datmrg$lucat <- ifelse(datmrg$LANDUSE1 == "1 0 0 Conservation and Natural Environments" |
                                       datmrg$LANDUSE1 == "2 1 0 Grazing natural vegetation" |
                                       datmrg$LANDUSE1 == "2 2 0 Production forestry of natural vegetation",
                                     "nat", "mod/agr")
datmrg$lucat <- ifelse(datmrg$LANDUSE1 == "5 5 0 Services" | 
                                       datmrg$LANDUSE1 == "5 4 0 Residential" |
                                       datmrg$LANDUSE1 == "5 6 0 Utilities" | 
                                       datmrg$lucat == "5 8 0 Mining",
                                     "built", datmrg$lucat)
table(datmrg$lucat)

# Hg sample location coordinates
Hgpts <- vect(cbind(datmrg$LON, datmrg$LAT), crs="+proj=longlat")
terra::plot(Hgpts)
Hgpts

## land use https://www.agriculture.gov.au/abares/aclump/land-use/data-download
lu <- rast('NLUM_v7_250_ALUMV8_2020_21_alb.tif')
terra::plot(lu)

# extract
lu.Hgpts <- terra::extract(lu, Hgpts)
head(lu.Hgpts)
table(lu.Hgpts$TERTV8)
lu.Hgpts$plu <- as.numeric(substr(lu.Hgpts$TERTV8, 1, 1))

lu.key <- data.frame(plu = 1:6, landuse = c("conservation/natural", "production-relatively natural", "production-dryland agr",
                                            "production-irrigated agr", "intensive", "water"))
# add code names for lu categories
lu.inf <- lu.Hgpts %>%
  left_join(lu.key, by="plu")
head(lu.inf)
table(lu.inf$landuse)
table(lu.Hgpts$plu)

# add to data
datmrg$landuse <- lu.inf$landuse

# date recognition
datmrg$posixdate <- as.POSIXct(datmrg$DATE, format="%d.%m.%Y")


## overlay data onto WWF ecoregions
# WWF ecoregion polygon # download shapefile from: https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
WWFecoregions <- vect("wwf_terr_ecos.shp") 
WWFecoregions
head(WWFecoregions['BIOME'])
table(WWFecoregions$BIOME)

# biome name key (https://databasin.org/datasets/68635d7c77f1475f9b6c1d1dbe0a4c4c/)
biome.names <- c('Tropical & Subtropical Moist Broadleaf Forests',
                 'Tropical & Subtropical Dry Broadleaf Forests',
                 'Tropical & Subtropical Coniferous Forests',
                 'Temperate Broadleaf & Mixed Forests',
                 'Temperate Conifer Forests',
                 'Boreal Forests/Taiga',
                 'Tropical & Subtropical Grasslands, Savannas & Shrublands',
                 'Temperate Grasslands, Savannas & Shrublands',
                 'Flooded Grasslands & Savannas',
                 'Montane Grasslands & Shrublands',
                 'Tundra',
                 'Mediterranean Forests, Woodlands & Scrub',
                 'Deserts & Xeric Shrublands',
                 'Mangroves')
biome.abbr <- c('TSMBF', "TSDBF", "TSCF", "TBMF", "TCF", "BFT", "TSGSS", "TGSS", "FGS", "MGS", "TNDRA", "MFW", "DXS", "MNGRV")
biome.mtype <- c('trop for', 'trop for', 'trop for', 'temp for', 'temp for', 'bor for', 'trop grass/sav',
                 'temp grass/sav', 'flood grass/sav', 'mont grass/shrub', 'tundra', 'Med for', 'des/xer', 'mangr')
biome.key <- data.frame('BIOME'=seq(1,14,1), 'name'=biome.names, "abbr"=biome.abbr, "type"=biome.mtype)
biome.key

# add Australia boundary
world <- geodata::world(path = tempdir())
aus <- world[world$NAME_0 == 'Australia',]
terra::plot(aus, add=T)

# extract
ecoreg.Hgpts <- terra::extract(WWFecoregions, Hgpts)

# add biome abbr & type
biome.inf <- ecoreg.Hgpts %>%
  left_join(biome.key, by="BIOME")
biome.pts <- biome.inf$abbr
btype.pts <- biome.inf$type

# add to data
datmrg$biome <- biome.pts
table(datmrg$biome)
datmrg$btype <- btype.pts
table(datmrg$btype)


## geology (1:1 M scale)
# https://ecat.ga.gov.au/geonetwork/js/api/records/c8856c41-0d5b-2b1d-e044-00144fdd4fa6
geol <- vect("GeologicUnitPolygons1M.shp")
geol
names(geol)
attr(table(geol$LITHOLOGY), 'names')

# import lithology reclassification key
lith.key <- read.csv("lithreclass.csv", header = T)
lith.key

lith.Hgpts <- terra::extract(geol, Hgpts)

lith.inf <- lith.Hgpts %>%
  left_join(lith.key, by="LITHOLOGY")
lith.pts <- lith.inf$LITHGRP

# add to data
datmrg$lithgrp <- lith.pts
table(datmrg$lithgrp)


## leaf area index
# https://ausenv.tern.org.au/aer/how-to-use-australias-environment-data-explorer/aer/australias-environment/index.html
# https://thredds.nci.org.au/thredds/catalog/ub8/au/OzWALD/annual/catalog.html
lai <- nc_open("OzWALD.LAI.AnnualMeans.nc")
print(lai)
lai.dat <- ncvar_get(lai, "AnnualMeans")
lon.lai <- ncvar_get(lai, "longitude")
lat.lai <- ncvar_get(lai, "latitude")

lai.rst <- rast(lai.dat)
ext(lai.rst) <- ext(min(lon.lai), max(lon.lai), min(lat.lai, na.rm=T), max(lat.lai, na.rm=T))
crs(lai.rst) <- "epsg:4326"  # set the coordinate reference system (typically WGS84)
nc_close(lai)
terra::plot(lai.rst)
lai.rst
lai.mn.rst <- app(lai.rst, fun = mean, na.rm = T)
terra::plot(lai.mn.rst)

lai.rst <- app(lai.mn.rst, function(x) {
  ifelse(x == 0, 0, log10(x))
})
lai.rst <- na.omit(lai.rst)
lai.rst

# extract
lai.Hgpts <- terra::extract(lai.mn.rst, Hgpts)
head(lai.Hgpts)

# add to data
datmrg$lai <- lai.Hgpts$mean


## radiometric data https://portal.ga.gov.au/persona/gadds
# KThU
KThU <- rast("radmap_v4_2019_filtered_ML_KThU_RGB_24bit.tif")
KThU.1 <- subset(KThU, 1) # r
KThU.2 <- subset(KThU, 2) # g
KThU.3 <- subset(KThU, 3) # b
KThUcmb <- (KThU.1 + KThU.2 + KThU.3)/3
terra::plot(KThUcmb)

KThU.rst <- app(KThUcmb, function(x) {
  ifelse(is.na(x) == T, NA, log10(x))
})
KThU.rst <- app(KThU.rst, function(x) {
  ifelse(is.infinite(x) == T, NA, x)
})
KThu.rst <- na.omit(KThU.rst)
KThu.rst
KThU.rsmp <- resample(KThU.rst, lai.rst)
plot(KThU.rsmp)
crsUse <- "+proj=longlat +datum=WGS84 +no_defs"
crs(KThU.rsmp) <- crsUse

# extract
KThU.Hgpts <- terra::extract(KThUcmb, Hgpts)
head(KThU.Hgpts)

# add to data
datmrg$KThU <- KThU.Hgpts$radmap_v4_2019_filtered_ML_KThU_RGB_24bit_1


# Th ppm
Th <- rast("radmap_v4_2019_filtered_ML_ppmTh_32bitfloat_grid.tif")
terra::plot(Th)

# extract
Th.Hgpts <- terra::extract(Th, Hgpts)
head(Th.Hgpts)

# add to data
datmrg$Th <- Th.Hgpts$radmap_v4_2019_filtered_ML_ppmTh_32bitfloat_grid


# U ppm
U <- rast("radmap_v4_2019_filtered_ML_ppmU_32bitfloat_grid.tif")
terra::plot(U)

# extract
U.Hgpts <- terra::extract(U, Hgpts)
head(U.Hgpts)
range(U.Hgpts$radmap_v4_2019_filtered_ML_ppmU_32bitfloat_grid)

# add to data
datmrg$U <- U.Hgpts$radmap_v4_2019_filtered_ML_ppmU_32bitfloat_grid


# pctk
pctk <- rast("radmap_v4_2019_filtered_ML_pctk_32bitfloat_grid.tif")
terra::plot(pctk)

# extract
pctk.Hgpts <- terra::extract(pctk, Hgpts)
head(pctk.Hgpts)
range(pctk.Hgpts$radmap_v4_2019_filtered_ML_pctk_32bitfloat_grid)

# add to data
datmrg$pctk <- pctk.Hgpts$radmap_v4_2019_filtered_ML_pctk_32bitfloat_grid


## soil nitrogen
# https://data.csiro.au/collection/csiro:61522?_st=browse&_str=2&_si=2&browseType=kw&browseValue=total%20soil%20nitrogen
soilN <- rast("NTO_000_005_EV_N_P_AU_NAT_C_20231101.tif")
terra::plot(soilN)
soilN.rst <- log10(soilN)
soilN.rsmp <- resample(soilN.rst, lai.rst)
plot(soilN.rsmp)
crs(soilN.rsmp) <- crsUse

# extract
soilN.Hgpts <- terra::extract(soilN, Hgpts)
head(soilN.Hgpts)

# add to data
datmrg$soilN <- soilN.Hgpts$focal_mean


## soil phosphorus
# https://data.csiro.au/collection/csiro:61526?_st=browse&_str=2&_si=1&browseType=kw&browseValue=total%20soil%20nitrogen
soilP <- rast("PTO_000_005_EV_N_P_AU_NAT_C_20231101.tif")
terra::plot(soilP)
soilP.rst <- log10(soilP)
soilP.rsmp <- resample(soilP.rst, lai.rst)
plot(soilP.rsmp)
crs(soilP.rsmp) <- crsUse

# extract
soilP.Hgpts <- terra::extract(soilP, Hgpts)
head(soilP.Hgpts)

# add to data
datmrg$soilP <- soilP.Hgpts$focal_mean


# soil pH (extracted)
soilpH <- rast("pHc_000_005_EV_N_P_AU_NAT_C_20140801.tif")
terra::plot(soilpH)
pH.rst <- soilpH
pH.rsmp <- resample(pH.rst, lai.rst)
plot(pH.rsmp)
crs(pH.rsmp) <- crsUse

# extract
soilpH.Hgpts <- terra::extract(soilpH, Hgpts)
head(soilpH.Hgpts)

# add to data
datmrg$soilpH <- soilpH.Hgpts$pHc_000_005_EV_N_P_AU_NAT_C_20140801


## rainfall
# https://ausenv.tern.org.au/aer/how-to-use-australias-environment-data-explorer/aer/australias-environment/index.html
# https://thredds.nci.org.au/thredds/catalog/ub8/au/OzWALD/annual/catalog.html
rain <- nc_open("OzWALD.annual.Pg.AnnualSums.nc")
print(rain)
rain.dat <- ncvar_get(rain, "AnnualSums")
lon.rain <- ncvar_get(rain, "longitude")
lat.rain <- ncvar_get(rain, "latitude")

rain.rst <- rast(rain.dat)
ext(rain.rst) <- ext(min(lon.rain), max(lon.rain), min(lat.rain, na.rm=T), max(lat.rain, na.rm=T))
crs(rain.rst) <- "epsg:4326"  # set the coordinate reference system (typically WGS84)
nc_close(rain)
terra::plot(rain.rst)
rain.rst
rain.mn.rst <- app(rain.rst, fun = mean, na.rm = T)
terra::plot(rain.mn.rst)

# extract
rain.Hgpts <- terra::extract(rain.mn.rst, Hgpts)
head(rain.Hgpts)

# add to data
datmrg$rain <- rain.Hgpts$mean


## Prescott index https://data.csiro.au/collection/csiro:9636v2
# a measure of water balance; designed to give an indication of the intensity of leaching by excess water
# and is calculated using long-term average precipitation P and potential evaporation E, both expressed as mean
# monthly values in mm (mean annual values divided by 12): PI = 0.445P / E^0.75
# Evaporation is estimated from temperature and net radiation; the net radiation is computed by the SRAD solar
# radiation model using the smoothed 1 arc-second resolution DEM-S (ANZCW0703014016) and includes both regional
# climatic influences and local topographic effects. Precipitation and temperature  obtained from national
# climate surfaces averaged over the same time period as the climatic information used in the radiation calculations
# (1981-2006). The Prescott Index has no units: larger values indicate wetter conditions.
prescott <- rast("PrescottIndex_01_3s_lzw.tif")
terra::plot(prescott)
Prescott.rst <- prescott
Prescott.rsmp <- resample(Prescott.rst, lai.rst)
plot(Prescott.rsmp)
crs(Prescott.rsmp) <- crsUse

# extract
prescott.Hgpts <- terra::extract(prescott, Hgpts)
head(prescott.Hgpts)

# add to data
datmrg$prescott <- prescott.Hgpts$PrescottIndex_01_3s_lzw


## soil water availability
# https://ausenv.tern.org.au/aer/how-to-use-australias-environment-data-explorer/aer/australias-environment/index.html
# https://thredds.nci.org.au/thredds/catalog/ub8/au/OzWALD/annual/catalog.html
soilH20 <- nc_open("OzWALD.Ssoil.AnnualMeans.nc")
print(soilH20)
soilH20.dat <- ncvar_get(soilH20, "AnnualMeans")
lon.soilH20 <- ncvar_get(soilH20, "longitude")
lat.soilH20 <- ncvar_get(soilH20, "latitude")

soilH20.rst <- rast(soilH20.dat)
ext(soilH20.rst) <- ext(min(lon.soilH20), max(lon.soilH20), min(lat.soilH20, na.rm=T), max(lat.soilH20, na.rm=T))
crs(soilH20.rst) <- "epsg:4326"  # set the coordinate reference system (typically WGS84)
nc_close(soilH20)
terra::plot(soilH20.rst)
soilH20.rst
soilH20.mn.rst <- app(soilH20.rst, fun = mean, na.rm = T)
terra::plot(soilH20.mn.rst)

soilH20.rst <- soilH20.mn.rst
soilH20.rsmp <- resample(soilH20.rst, lai.rst)
plot(soilH20.rsmp)

# extract
soilH20.Hgpts <- terra::extract(soilH20.mn.rst, Hgpts)
head(soilH20.Hgpts)

# add to data
datmrg$soilH20 <- soilH20.Hgpts$mean


# do lai.rsmp after defining soilH20.rsmp
lai.rsmp <- resample(lai.rst, soilH20.rsmp)
plot(lai.rsmp)
crs(lai.rsmp) <- crsUse


## vegetation carbon uptake (GPP)
# amount of carbon taken up by the vegetation through photosynthesis,
# as estimated by the OzWALD model-data fusion system
# https://ausenv.tern.org.au/aer/how-to-use-australias-environment-data-explorer/aer/australias-environment/index.html
# https://thredds.nci.org.au/thredds/catalog/ub8/au/OzWALD/annual/catalog.html
gpp <- nc_open("OzWALD.GPP.AnnualMeans.nc")
print(gpp)
gpp.dat <- ncvar_get(gpp, "AnnualMeans")
lon.gpp <- ncvar_get(gpp, "longitude")
lat.gpp <- ncvar_get(gpp, "latitude")

gpp.rst <- rast(gpp.dat)
ext(gpp.rst) <- ext(min(lon.gpp), max(lon.gpp), min(lat.gpp, na.rm=T), max(lat.gpp, na.rm=T))
crs(gpp.rst) <- "epsg:4326"  # set the coordinate reference system (typically WGS84)
nc_close(gpp)
terra::plot(gpp.rst)
gpp.rst
gpp.mn.rst <- app(gpp.rst, fun = mean, na.rm = T)
terra::plot(gpp.mn.rst)

# extract
gpp.Hgpts <- terra::extract(gpp.mn.rst, Hgpts)
head(gpp.Hgpts)

# add to data
datmrg$gpp <- gpp.Hgpts$mean


## soil clay content (%) extracted
# https://data.csiro.au/collection/csiro:55684
soilclay <- rast("CLY_000_005_EV_N_P_AU_TRN_N_20210902.tif")
terra::plot(soilclay)
clay.rst <- log(soilclay/100)/(1-(soilclay/100))
clay.rsmp <- resample(clay.rst, lai.rst)
plot(clay.rsmp)
crs(clay.rsmp) <- crsUse

# extract
soilclay.Hgpts <- terra::extract(soilclay, Hgpts)
head(soilclay.Hgpts)

# add to data
datmrg$soilclay <- soilclay.Hgpts$CLY_000_005_EV_N_P_AU_TRN_N_20210902


## soil silt content (%) extracted
# https://data.csiro.au/collection/csiro:10688?q=soil%20silt&_st=keyword&_str=12&_si=1
soilsilt <- rast("SLT_000_005_EV_N_P_AU_TRN_N_20210902.tif")
terra::plot(soilsilt)

# extract
soilsilt.Hgpts <- terra::extract(soilsilt, Hgpts)
head(soilsilt.Hgpts)

# add to data
datmrg$soilsilt <- soilsilt.Hgpts$SLT_000_005_EV_N_P_AU_TRN_N_20210902


## soil categories
## https://data.csiro.au/collection/csiro:40340
## https://doi.org/10.25919/5f1632a855c17
## citation: CSIRO; & National Resource Information Centre, BRS (1991): Atlas of Australian Soils (digital). v3. 
## CSIRO. Data Collection. https://doi.org/10.25919/5f1632a855c17
## Soil type classification based on the Australian Soil Classification (Isbell, 1996)
## https://www.soilscienceaustralia.org.au/asc/soilhome.htm
soil <- vect("soilAtlas2M.shp")
soil
names(soil)
attr(table(soil$MAP_UNIT), 'names')

# soil attributes
soil.key <- read.csv("asclut.csv")
head(soil.key)

soils <- merge(soil, soil.key, by="MAP_UNIT")
soils
attr(table(soils$TYPE), 'names')

# extract
soils.Hgpts <- terra::extract(soils, Hgpts)
head(soils.Hgpts)

# add to data
datmrg$soils <- soils.Hgpts$TYPE
table(datmrg$soils)
sum(as.data.frame(table(datmrg$soils))$Freq)


# enhanced barest earth (proxy for soil categories on continuous scale)
# https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/144231
# https://pid.geoscience.gov.au/dataset/ga/144231
# https://doi.org/10.26186/144231
# Wilford, J. and Roberts, D., 2020. Enhanced barest earth Landsat imagery for soil
# and lithological modelling. In: Czarnota, K., Roach, I., Abbott, S., Haynes, M.,
# Kositcin, N., Ray, A. and Slatter, E. (eds.) Exploring for the Future: Extended Abstracts,
# Geoscience Australia, Canberra, 1–4
# 1. BLUE; # 2. GREEN; # 3. RED; # 4. NIR; # 5. SWIR1; # 6. SWIR2
# Normalised ratios bands
# ((RED - BLUE) / (RED + BLUE)) == "-ND-RED-BLUE.tif"
# ((SWIR1 - NIR) / (SWIR1 + NIR)) == "-ND-SWIR1-NIR.tif"
# ((SWIR1 - SWIR2) / (SWIR1 + SWIR2)) == "-ND-SWIR1-SWIR2.tif"
# ((NIR - GREEN) / (NIR + GREEN)) == "-ND-NIR-GREEN.tif"
# ((SWIR1 - BLUE) / (SWIR1 + BLUE)) == "-ND-SWIR1-BLUE.tif"
# ((SWIR2 - NIR) / (SWIR2 + NIR)) == "-ND-SWIR2-NIR.tif"
# ((SWIR2 - RED) / (SWIR2 + RED)) =="-ND-SWIR2-RED.tif"
# ((RED - GREEN) / (RED + GREEN)) == "-ND-RED-GREEN.tif"
# ((SWIR2 - GREEN) / (SWIR2 + GREEN)) == "-ND-SWIR2-GREEN.tif"
# Ferric PC2 of BLUE and RED == "-FERRIC-PC2.tif"(based on Chavez and Kwarteng 1989)
# Ferric - PC4 of BLUE, RED, NIR, SWIR1 == "-FERRIC-PC4.tif" (based on Loughin, 1991)
# Hydroxyl (clay) PC 1 of SWIR1/SWIR2 and NIR/RED (Fraser and Green, 1987) == "-HYDROXYL-1-PC1.tif"
# (clay)
# Hydroxyl (clay) PC 2 of SWIR1/SWIR2 and NIR/RED (Fraser and Green, 1987) == "-HYDROXYL-1-PC2.tif"
# (mixed vegetation and clay)
# Hydroxyl - PC2 of SWIR1, SWIR2 == "-HYDROXYL-2-PC2.tif"(based on Chavez and Kwarteng 1989)
# Hydroxyl - PC3, PC4 of BLUE, NIR, SWIR1, SWIR2 == "-HYDROXYL-3-PC3.tif", "-HYDROXYL-3-PC4.tif"(based on Loughin, 1991)
# Band additions: # Carbonate/quartz (non-clays highly reflective) BLUE + SWIR2
ferricpc2 <- rast("FERRICPC2mrg.tif")

crs(ferricpc2)
crs(ferricpc2) <- "epsg:3577"
ext(ferricpc2)
terra::plot(ferricpc2)

aus.terr <- vect("aus.shp")
crs(aus.terr)
ext(aus.terr)

# project aus.terr to "epsg:3577"
aus.terr3577 <- project(aus.terr, "epsg:3577")
terra::plot(aus.terr3577)

# find projection boundaries for aus.terr9473
ext3577 <- ext(aus.terr3577)
ext3577

# clip to aus.terr
ferricpc2.clip <- crop(ferricpc2, aus.terr3577)

# project ferricpc4.clip to "epsg:9473"
ferricpc2clip9473 <- project(ferricpc2.clip, "epsg:9473")
terra::plot(ferricpc2clip9473)

# write raster
writeRaster(ferricpc2clip9473, "ferricpc2clip9473.tif", overwrite=T)
ferricpc2clip9473 <- rast("ferricpc2clip9473.tif")

# project back to WGS84
ferricpc2clip4326 <- project(ferricpc2clip9473, "epsg:4326")
ferric2.rst <- ferricpc2clip4326
ferric2.rsmp <- resample(ferric2.rst, lai.rst)
plot(ferric2.rsmp)
crs(ferric2.rsmp) <- crsUse

# extract
ferricpc2.Hgpts <- terra::extract(ferricpc2clip4326, Hgpts)
head(ferricpc2.Hgpts)

# add to data
datmrg$ferricpc2 <- ferricpc2.Hgpts$FERRICPC2mrg


## aluminium oxide
## https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/148587
## doi:10.26186/148587
Alox <- rast("Aluminium_oxide_prediction_median.tif")
terra::plot(Alox)
Alox9473 <- project(Alox, "epsg:9473")
terra::plot(Alox9473,col=rev(map.pal("sepia",200)))
Alox4326 <- project(Alox9473, crs(Hgpts))
Alox.rst <- Alox4326
Alox.rsmp <- resample(Alox.rst, lai.rst)
plot(Alox.rsmp)
crs(Alox.rsmp) <- crsUse

# extract
Alox.Hgpts <- terra::extract(Alox4326, Hgpts)
head(Alox.Hgpts)

# add to data
datmrg$Alox <- Alox.Hgpts$Aluminium_oxide_prediction_median


## iron oxide
## https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/148587
## doi:10.26186/148587
Feox9473 <- rast("Feox9473.tif")
terra::plot(Feox9473,col=(map.pal("gyr",200)))
Feox4326 <- project(Feox9473, crs(Hgpts))
Feox.rst <- Feox4326
Feox.rsmp <- resample(Feox.rst, lai.rst)
plot(Feox.rsmp)
crs(Feox.rsmp) <- crsUse

# extract
Feox.Hgpts <- terra::extract(Feox4326, Hgpts)
head(Feox.Hgpts)

# add to data
datmrg$Feox <- Feox.Hgpts$Feox9473


## phosphorus oxide
## https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/148587
## doi:10.26186/148587
Pox9473 <- rast("Pox9473.tif")
terra::plot(Pox9473,col=rev(map.pal("plasma",200)))
Pox4326 <- project(Pox9473, crs(Hgpts))
terra::saveRDS(Pox4326, "Pox4326.rds")
Pox.rst <- Pox4326
Pox.rsmp <- resample(Pox.rst, lai.rst)
plot(Pox.rsmp)
crs(Pox.rsmp) <- crsUse

# extract
Pox.Hgpts <- terra::extract(Pox4326, Hgpts)
head(Pox.Hgpts)

# add to data
datmrg$Pox <- Pox.Hgpts$Pox9473



############################################
############################################
## separate bottom from top of soil sample
############################################
# bottom
BOS.dat <- datmrg[datmrg$DEPTH == "BOS",]
head(BOS.dat)
dim(BOS.dat)

# top
TOS.dat <- datmrg[datmrg$DEPTH == "TOS",]
head(TOS.dat)
dim(TOS.dat)

## Hg columns
head(TOS.dat[,grep("Hg", colnames(TOS.dat))])

## bivariate plots
## time plots
# Hg vs. time
plot(datmrg$posixdate, log10(datmrg$HgCOMP), xlab="date", ylab="log10 [Hg]", pch=19, col="blue")

# Hg vs. mean grain size
plot(datmrg$gs.wmn, log10(datmrg$HgCOMP), xlab="weighted mean grain size", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(datmrg$HgCOMP) ~ datmrg$gs.wmn), col="red", lwd=2, lty=2)

# Hg vs. pH (field)
plot(TOS.dat$FIELDpH, log10(TOS.dat$HgCOMP), xlab="TOS field pH", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$FIELDpH), col="red", lwd=2, lty=2)
plot(BOS.dat$FIELDpH, log10(BOS.dat$HgCOMP), xlab="BOS field pH", ylab="log10 BOS [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$FIELDpH), col="red", lwd=2, lty=2)

# Hg vs. pH (1:5)
plot(TOS.dat$pH15, log10(TOS.dat$HgCOMP), xlab="TOS 1:5 pH", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$pH15), col="red", lwd=2, lty=2)
plot(BOS.dat$pH15, log10(BOS.dat$HgCOMP), xlab="BOS 1:5 pH", ylab="log10 BOS [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$pH15), col="red", lwd=2, lty=2)

# soil pH (extracted) vs. Hg
plot(TOS.dat$soilpH, log10(TOS.dat$HgCOMP), xlab="TOS soil pH", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$soilpH), col="red", lwd=2, lty=2)
plot(BOS.dat$soilpH, log10(BOS.dat$HgCOMP), xlab="BOS soil pH", ylab="log10 BOS [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$soilpH), col="red", lwd=2, lty=2)

# soil organic matter (LOI) vs. Hg
plot(TOS.dat$LOI, log10(TOS.dat$HgCOMP), xlab="TOS LOI", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$LOI), col="red", lwd=2, lty=2)
plot(BOS.dat$LOI, log10(BOS.dat$HgCOMP), xlab="BOS LOI", ylab="log10 BOS [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$LOI), col="red", lwd=2, lty=2)

# soil electroconductivity (EC15) vs. Hg
TOS.EC15.Hg <- na.omit(data.frame(EC15=TOS.dat$EC15, HgCOMP=log10(TOS.dat$HgCOMP)))
plot(log10(TOS.dat$EC15), log10(TOS.dat$HgCOMP), xlab="TOS soil electroconductivity", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ log10(TOS.dat$EC15)), col="red", lwd=2, lty=2)
plot(log10(BOS.dat$EC15), log10(BOS.dat$HgCOMP), xlab="BOS soil electroconductivity", ylab="log10 BOS [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ log10(BOS.dat$EC15)), col="red", lwd=2, lty=2)

# % clay vs. Hg
plot(TOS.dat$CLAY, log10(TOS.dat$HgCOMP), xlab="TOS % clay", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$CLAY), col="red", lwd=2, lty=2)
plot(BOS.dat$CLAY, log10(BOS.dat$HgCOMP), xlab="BOS % clay", ylab="log10 BOS [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$CLAY), col="red", lwd=2, lty=2)

# % soil clay vs. Hg (extracted)
plot(TOS.dat$soilclay, log10(TOS.dat$HgCOMP), xlab="TOS % soil clay", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$soilclay), col="red", lwd=2, lty=2)
plot(BOS.dat$soilclay, log10(BOS.dat$HgCOMP), xlab="BOS % soil clay", ylab="log10 BOS [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$soilclay), col="red", lwd=2, lty=2)

# % silt vs. Hg
plot(TOS.dat$SILT, log10(TOS.dat$HgCOMP), xlab="TOS % silt", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$SILT), col="red", lwd=2, lty=2)
plot(BOS.dat$SILT, log10(BOS.dat$HgCOMP), xlab="BOS % silt", ylab="log10 BOS [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$SILT), col="red", lwd=2, lty=2)

# % soil silt vs. Hg (extracted)
plot(TOS.dat$soilsilt, log10(TOS.dat$HgCOMP), xlab="TOS % soil silt", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$soilsilt), col="red", lwd=2, lty=2)
plot(BOS.dat$soilsilt, log10(BOS.dat$HgCOMP), xlab="BOS % soil silt", ylab="log10 BOS [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$soilsilt), col="red", lwd=2, lty=2)

# % sand vs. Hg
plot(TOS.dat$SAND, log10(TOS.dat$HgCOMP), xlab="TOS % sand", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$SAND), col="red", lwd=2, lty=2)
plot(BOS.dat$SAND, log10(BOS.dat$HgCOMP), xlab="BOS % sand", ylab="log10 BOS [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$SAND), col="red", lwd=2, lty=2)

# elements
# Al vs. Hg
plot(TOS.dat$Al, log10(TOS.dat$HgCOMP), xlab="Al", ylab="log10 [Hg]", pch=19, col="blue")
plot(BOS.dat$Al, log10(BOS.dat$HgCOMP), xlab="Al", ylab="log10 [Hg]", pch=19, col="blue")

# Fe vs. Hg
plot(log10(TOS.dat$Fe), log10(TOS.dat$HgCOMP), xlab="Fe", ylab="log10 [Hg]", pch=19, col="blue")
plot(log10(BOS.dat$Fe), log10(BOS.dat$HgCOMP), xlab="Fe", ylab="log10 [Hg]", pch=19, col="blue")

# Mn vs. Hg
plot(log10(TOS.dat$Mn), log10(TOS.dat$HgCOMP), xlab="Mn", ylab="log10 [Hg]", pch=19, col="blue")
plot(log10(BOS.dat$Mn), log10(BOS.dat$HgCOMP), xlab="Mn", ylab="log10 [Hg]", pch=19, col="blue")

# Cu vs. Hg
plot(log10(TOS.dat$CuICPMS), log10(TOS.dat$HgCOMP), xlab="Cu", ylab="log10 [Hg]", pch=19, col="blue")
plot(log10(BOS.dat$CuICPMS), log10(BOS.dat$HgCOMP), xlab="Cu", ylab="log10 [Hg]", pch=19, col="blue")

# Zn vs. Hg
plot(log10(as.numeric(TOS.dat$ZnICPMS)), log10(TOS.dat$HgCOMP), xlab="Zn", ylab="log10 [Hg]", pch=19, col="blue")
plot(log10(as.numeric(BOS.dat$ZnICPMS)), log10(BOS.dat$HgCOMP), xlab="Zn", ylab="log10 [Hg]", pch=19, col="blue")

# Pb vs. Hg
plot(log10(TOS.dat$PbICPMS), log10(TOS.dat$HgCOMP), xlab="Pb", ylab="log10 [Hg]", pch=19, col="blue")
plot(log10(BOS.dat$PbICPMS), log10(BOS.dat$HgCOMP), xlab="Pb", ylab="log10 [Hg]", pch=19, col="blue")

# Sb vs. Hg
plot(log10(TOS.dat$SbICPMS), log10(TOS.dat$HgCOMP), xlab="Sb", ylab="log10 [Hg]", pch=19, col="blue")
plot(log10(BOS.dat$SbICPMS), log10(BOS.dat$HgCOMP), xlab="Sb", ylab="log10 [Hg]", pch=19, col="blue")

# Ni vs. Hg
plot(log10(TOS.dat$NiCPMS), log10(TOS.dat$HgCOMP), xlab="Ni", ylab="log10 [Hg]", pch=19, col="blue")
plot(log10(BOS.dat$NiCPMS), log10(BOS.dat$HgCOMP), xlab="Ni", ylab="log10 [Hg]", pch=19, col="blue")

# Va vs. Hg
plot(log10(TOS.dat$V), log10(TOS.dat$HgCOMP), xlab="Va", ylab="log10 [Hg]", pch=19, col="blue")
plot(log10(BOS.dat$V), log10(BOS.dat$HgCOMP), xlab="Va", ylab="log10 [Hg]", pch=19, col="blue")

# oxides
# Alox
plot((TOS.dat$Alox), log10(TOS.dat$HgCOMP), xlab="log10 Al oxide", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ (TOS.dat$Alox)), col="red", lwd=2, lty=2)

# Feox
plot((TOS.dat$Feox), log10(TOS.dat$HgCOMP), xlab="log10 Fe oxide", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ (TOS.dat$Feox)), col="red", lwd=2, lty=2)

# Pox
plot((TOS.dat$Pox), log10(TOS.dat$HgCOMP), xlab="log10 P oxide", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ (TOS.dat$Pox)), col="red", lwd=2, lty=2)

# radiometric
# KThU
plot(log10(TOS.dat$KThU), log10(TOS.dat$HgCOMP), xlab="log10 KThU", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ log10(TOS.dat$KThU)), col="red", lwd=2, lty=2)
plot(log10(BOS.dat$KThU), log10(BOS.dat$HgCOMP), xlab="log10 KThU", ylab="log10 BOS [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ log10(BOS.dat$KThU)), col="red", lwd=2, lty=2)

# Th
plot(log10(TOS.dat$Th), log10(TOS.dat$HgCOMP), xlab="log10 Th", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ log10(TOS.dat$Th)), col="red", lwd=2, lty=2)

# U
plot(log10(TOS.dat$U), log10(TOS.dat$HgCOMP), xlab="log10 U", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ log10(TOS.dat$U)), col="red", lwd=2, lty=2)

# pctk
plot(log10(TOS.dat$pctk), log10(TOS.dat$HgCOMP), xlab="log10 pctk", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ log10(TOS.dat$pctk)), col="red", lwd=2, lty=2)

# soil N vs. Hg
plot(log10(TOS.dat$soilN), log10(TOS.dat$HgCOMP), xlab="soil N", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ log10(TOS.dat$soilN)), col="red", lwd=2, lty=2)
plot(log10(BOS.dat$soilN), log10(BOS.dat$HgCOMP), xlab="soil N", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ log10(BOS.dat$soilN)), col="red", lwd=2, lty=2)

# soil P vs. Hg
plot((TOS.dat$soilP), log10(TOS.dat$HgCOMP), xlab="soil P", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ (TOS.dat$soilP)), col="red", lwd=2, lty=2)
plot((BOS.dat$soilP), log10(BOS.dat$HgCOMP), xlab="soil P", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ (BOS.dat$soilP)), col="red", lwd=2, lty=2)

# ferric PC2 barest earth soil vs. Hg
plot((TOS.dat$ferricpc2), log10(TOS.dat$HgCOMP), xlab="ferric PC4", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$ferricpc2), col="red", lwd=2, lty=2)
hist(TOS.dat$ferricpc2)

# latitude vs. Hg
plot(datmrg$LAT, log10(datmrg$HgCOMP), xlab="latitude", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(datmrg$HgCOMP) ~ datmrg$LAT), col="red", lwd=2, lty=2)

# longitude vs. Hg
plot(datmrg$LON, log10(datmrg$HgCOMP), xlab="longitude", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(datmrg$HgCOMP) ~ datmrg$LON), col="red", lwd=2, lty=2)

# rainfall vs. Hg
plot(log10(TOS.dat$rain), log10(TOS.dat$HgCOMP), xlab="annual rainfall (mm)", ylab="TOS log10 [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ log10(TOS.dat$rain)), col="red", lwd=2, lty=2)
plot(log10(BOS.dat$rain), log10(BOS.dat$HgCOMP), xlab="annual rainfall (mm)", ylab="BOS log10 [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ log10(BOS.dat$rain)), col="red", lwd=2, lty=2)

# LAI vs. Hg
plot((TOS.dat$lai), log10(TOS.dat$HgCOMP), xlab="leaf area index", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$lai), col="red", lwd=2, lty=2)
plot((BOS.dat$lai), log10(BOS.dat$HgCOMP), xlab="leaf area index", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$lai), col="red", lwd=2, lty=2)

# GPP vs. Hg
plot((TOS.dat$gpp), log10(TOS.dat$HgCOMP), xlab="GPP", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$gpp), col="red", lwd=2, lty=2)
plot((BOS.dat$gpp), log10(BOS.dat$HgCOMP), xlab="GPP", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$gpp), col="red", lwd=2, lty=2)

# Prescott index vs. Hg
plot((TOS.dat$prescott), log10(TOS.dat$HgCOMP), xlab="Prescott index", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$prescott), col="red", lwd=2, lty=2)
plot((BOS.dat$prescott), log10(BOS.dat$HgCOMP), xlab="Prescott index", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$prescott), col="red", lwd=2, lty=2)

# Hg by state
TOS.HgXstate <- sort(xtabs(TOS.dat$HgCOMP ~ TOS.dat$STATE) / table(TOS.dat$STATE))
barplot(TOS.HgXstate, xlab="state", ylab="mean TOS [Hg]", col="blue")
BOS.HgXstate <- sort(xtabs(BOS.dat$HgCOMP ~ BOS.dat$STATE) / table(BOS.dat$STATE))
barplot(BOS.HgXstate, xlab="state", ylab="mean BOS [Hg]", col="blue")

TOSHgXstate.stats <- TOS.dat %>%
    group_by(STATE) %>%
    summarise(
      mean = mean(HgCOMP, na.rm = TRUE),
      median = median(HgCOMP, na.rm = TRUE), 
      var = var(HgCOMP, na.rm = TRUE),
      sd = sd(HgCOMP, na.rm = TRUE),
      se = sd/sqrt(n()),
      upper = quantile(HgCOMP, probs=0.975, na.rm = TRUE),
      lower = quantile(HgCOMP, probs=0.025, na.rm = TRUE),
      n = n()
  )

TOSHgXstate.sort <- TOSHgXstate.stats %>% arrange(mean)
print(TOSHgXstate.sort, digits = 5, width = Inf)

TOSHgXstate.sort$mean
TOSHgXstate.sort$se

ggplot(TOSHgXstate.sort, aes(x = reorder(STATE, mean), y = mean)) +
  geom_bar(stat='identity', fill = "steelblue") +
  geom_errorbar(aes(ymin = mean-(2*se), ymax = mean+(2*se)),
                width = 0.2) +
  labs(x = "", y = "mean [Hg] ± 2 se") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))

BOSHgXstate.stats <- BOS.dat %>%
  group_by(STATE) %>%
  summarise(
    mean = mean(HgCOMP, na.rm = TRUE),
    median = median(HgCOMP, na.rm = TRUE), 
    var = var(HgCOMP, na.rm = TRUE),
    sd = sd(HgCOMP, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(HgCOMP, probs=0.975, na.rm = TRUE),
    lower = quantile(HgCOMP, probs=0.025, na.rm = TRUE),
    n = n()
  )

BOSHgXstate.sort <- BOSHgXstate.stats %>% arrange(mean)
print(BOSHgXstate.sort, digits = 5, width = Inf)

BOSHgXstate.sort$mean
BOSHgXstate.sort$se

ggplot(BOSHgXstate.sort, aes(x = reorder(STATE, mean), y = mean)) +
  geom_bar(stat='identity', fill = "steelblue") +
  geom_errorbar(aes(ymin = mean-(2*se), ymax = mean+(2*se)),
                width = 0.2) +
  labs(x = "", y = "mean [Hg] ± 2 se") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))



# Hg by landuse
TOS.HgXlucat <- sort(xtabs(TOS.dat$HgCOMP ~ TOS.dat$lucat) / table(TOS.dat$lucat))
barplot(TOS.HgXlucat, xlab="landuse category", ylab="mean TOS [Hg]", col="blue")
BOS.HgXlucat <- sort(xtabs(BOS.dat$HgCOMP ~ BOS.dat$lucat) / table(BOS.dat$lucat))
barplot(BOS.HgXlucat, xlab="landuse category", ylab="mean BOS [Hg]", col="blue")

TOSHgXlucat.stats <- TOS.dat %>%
  group_by(lucat) %>%
  summarise(
    mean = mean(HgCOMP, na.rm = TRUE),
    median = median(HgCOMP, na.rm = TRUE), 
    var = var(HgCOMP, na.rm = TRUE),
    sd = sd(HgCOMP, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(HgCOMP, probs=0.975, na.rm = TRUE),
    lower = quantile(HgCOMP, probs=0.025, na.rm = TRUE),
    n = n()
  )

TOSHgXlucat.sort <- TOSHgXlucat.stats %>% arrange(mean)
TOSHgXlucat.sort
ggplot(TOSHgXlucat.sort, aes(x = reorder(lucat, mean), y = mean)) +
  geom_bar(stat='identity', fill = "steelblue") +
  geom_errorbar(aes(ymin = mean-(2*se), ymax = mean+(2*se)),
                width = 0.2) +
  labs(x = "land-use category", y = "mean [Hg] ± 2 se") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))

BOSHgXlucat.stats <- BOS.dat %>%
  group_by(lucat) %>%
  summarise(
    mean = mean(HgCOMP, na.rm = TRUE),
    median = median(HgCOMP, na.rm = TRUE), 
    var = var(HgCOMP, na.rm = TRUE),
    sd = sd(HgCOMP, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(HgCOMP, probs=0.975, na.rm = TRUE),
    lower = quantile(HgCOMP, probs=0.025, na.rm = TRUE),
    n = n()
  )

BOSHgXlucat.sort <- BOSHgXlucat.stats %>% arrange(mean)
BOSHgXlucat.sort
ggplot(BOSHgXlucat.sort, aes(x = reorder(lucat, mean), y = mean)) +
  geom_bar(stat='identity', fill = "steelblue") +
  geom_errorbar(aes(ymin = mean-(2*se), ymax = mean+(2*se)),
                width = 0.2) +
  labs(x = "land-use category", y = "mean [Hg] ± 2 se") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))



TOS.HgXlanduse <- sort(xtabs(TOS.dat$HgCOMP ~ TOS.dat$landuse) / table(TOS.dat$landuse))
barplot(TOS.HgXlanduse, xlab="primary landuse category", ylab="mean TOS [Hg]", col="blue")
BOS.HgXlanduse <- sort(xtabs(BOS.dat$HgCOMP ~ BOS.dat$landuse) / table(BOS.dat$landuse))
barplot(BOS.HgXlanduse, xlab="primary landuse category", ylab="mean BOS [Hg]", col="blue")

TOSHgXlanduse.stats <- TOS.dat %>%
  group_by(landuse) %>%
  summarise(
    mean = mean(HgCOMP, na.rm = TRUE),
    median = median(HgCOMP, na.rm = TRUE), 
    var = var(HgCOMP, na.rm = TRUE),
    sd = sd(HgCOMP, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(HgCOMP, probs=0.975, na.rm = TRUE),
    lower = quantile(HgCOMP, probs=0.025, na.rm = TRUE),
    n = n()
  )

TOSHgXlanduse.sort <- TOSHgXlanduse.stats %>% arrange(mean)
TOSHgXlanduse.sort

TOSHgXlanduse.sort$mean
TOSHgXlanduse.sort$se

ggplot(TOSHgXlanduse.sort, aes(x = reorder(landuse, mean), y = mean)) +
  geom_bar(stat='identity', fill = "steelblue") +
  geom_errorbar(aes(ymin = mean-(2*se), ymax = mean+(2*se)),
                width = 0.2) +
  labs(x = "land-use category", y = "mean [Hg] ± 2 se") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))

BOSHgXlanduse.stats <- BOS.dat %>%
  group_by(landuse) %>%
  summarise(
    mean = mean(HgCOMP, na.rm = TRUE),
    median = median(HgCOMP, na.rm = TRUE), 
    var = var(HgCOMP, na.rm = TRUE),
    sd = sd(HgCOMP, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(HgCOMP, probs=0.975, na.rm = TRUE),
    lower = quantile(HgCOMP, probs=0.025, na.rm = TRUE),
    n = n()
  )

BOSHgXlanduse.sort <- BOSHgXlanduse.stats %>% arrange(mean)
BOSHgXlanduse.sort

BOSHgXlanduse.sort$mean
BOSHgXlanduse.sort$se

ggplot(BOSHgXlanduse.sort, aes(x = reorder(landuse, mean), y = mean)) +
  geom_bar(stat='identity', fill = "steelblue") +
  geom_errorbar(aes(ymin = mean-(2*se), ymax = mean+(2*se)),
                width = 0.2) +
  labs(x = "land-use category", y = "mean [Hg] ± 2 se") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))


# Hg by landform
TOS.HgXlndfrm <- sort(xtabs(TOS.dat$HgCOMP ~ TOS.dat$LANDFORM.TYPE) / table(TOS.dat$LANDFORM.TYPE))
barplot(TOS.HgXlndfrm, xlab="landform", ylab="mean TOS [Hg]", col="blue")
BOS.HgXlndfrm <- sort(xtabs(BOS.dat$HgCOMP ~ BOS.dat$LANDFORM.TYPE) / table(BOS.dat$LANDFORM.TYPE))
barplot(BOS.HgXlndfrm, xlab="landform", ylab="mean BOS [Hg]", col="blue")

# Hg by biome
TOS.HgXbiome <- sort(xtabs(TOS.dat$HgCOMP ~ TOS.dat$biome) / table(TOS.dat$biome))
barplot(TOS.HgXbiome, xlab="biome", ylab="mean TOS [Hg]", col="blue")
BOS.HgXbiome <- sort(xtabs(BOS.dat$HgCOMP ~ BOS.dat$biome) / table(BOS.dat$biome))
barplot(BOS.HgXbiome, xlab="biome", ylab="mean BOS [Hg]", col="blue")

# Hg by biome type
TOS.HgXbtype <- sort(xtabs(TOS.dat$HgCOMP ~ TOS.dat$btype) / table(TOS.dat$btype))
barplot(TOS.HgXbtype, xlab="biome type", ylab="mean TOS [Hg]", col="blue")
BOS.HgXbtype <- sort(xtabs(BOS.dat$HgCOMP ~ BOS.dat$btype) / table(BOS.dat$btype))
barplot(BOS.HgXbtype, xlab="biome type", ylab="mean BOS [Hg]", col="blue")

TOSHgXbtype.stats <- TOS.dat %>%
  group_by(btype) %>%
  summarise(
    mean = mean(HgCOMP, na.rm = TRUE),
    median = median(HgCOMP, na.rm = TRUE), 
    var = var(HgCOMP, na.rm = TRUE),
    sd = sd(HgCOMP, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(HgCOMP, probs=0.975, na.rm = TRUE),
    lower = quantile(HgCOMP, probs=0.025, na.rm = TRUE),
    n = n()
  )

TOSHgXbtype.sort <- TOSHgXbtype.stats %>% arrange(mean)
TOSHgXbtype.sort

TOSHgXbtype.sort$mean
TOSHgXbtype.sort$se

ggplot(TOSHgXbtype.sort, aes(x = reorder(btype, mean), y = mean)) +
  geom_bar(stat='identity', fill = "steelblue") +
  geom_errorbar(aes(ymin = mean-(2*se), ymax = mean+(2*se)),
                width = 0.2) +
  labs(x = "biome type", y = "mean [Hg] ± 2 se") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))

BOSHgXbtype.stats <- BOS.dat %>%
  group_by(btype) %>%
  summarise(
    mean = mean(HgCOMP, na.rm = TRUE),
    median = median(HgCOMP, na.rm = TRUE), 
    var = var(HgCOMP, na.rm = TRUE),
    sd = sd(HgCOMP, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(HgCOMP, probs=0.975, na.rm = TRUE),
    lower = quantile(HgCOMP, probs=0.025, na.rm = TRUE),
    n = n()
  )

BOSHgXbtype.sort <- BOSHgXbtype.stats %>% arrange(mean)
BOSHgXbtype.sort

BOSHgXbtype.sort$mean
BOSHgXbtype.sort$se

ggplot(BOSHgXbtype.sort, aes(x = reorder(btype, mean), y = mean)) +
  geom_bar(stat='identity', fill = "steelblue") +
  geom_errorbar(aes(ymin = mean-(2*se), ymax = mean+(2*se)),
                width = 0.2) +
  labs(x = "biome type", y = "mean [Hg] ± 2 se") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))


# Hg by lithology
TOS.HgXlith <- sort(xtabs(TOS.dat$HgCOMP ~ TOS.dat$lithgrp) / table(TOS.dat$lithgrp))
barplot(TOS.HgXlith, xlab="lithology", ylab="mean TOS [Hg]", col="blue")
BOS.HgXlith <- sort(xtabs(BOS.dat$HgCOMP ~ BOS.dat$lithgrp) / table(BOS.dat$lithgrp))
barplot(BOS.HgXlith, xlab="lithology", ylab="mean BOS [Hg]", col="blue")

TOSHgXlith.stats <- TOS.dat %>%
  group_by(lithgrp) %>%
  summarise(
    mean = mean(HgCOMP, na.rm = TRUE),
    median = median(HgCOMP, na.rm = TRUE), 
    var = var(HgCOMP, na.rm = TRUE),
    sd = sd(HgCOMP, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(HgCOMP, probs=0.975, na.rm = TRUE),
    lower = quantile(HgCOMP, probs=0.025, na.rm = TRUE),
    n = n()
  )

TOSHgXlith.sort <- na.omit(TOSHgXlith.stats %>% arrange(mean))
TOSHgXlith.sort

TOSHgXlith.sort$mean
TOSHgXlith.sort$se

ggplot(TOSHgXlith.sort, aes(x = reorder(lithgrp, mean), y = mean)) +
  geom_bar(stat='identity', fill = "steelblue") +
  geom_errorbar(aes(ymin = mean-(2*se), ymax = mean+(2*se)),
                width = 0.2) +
  labs(x = "lithology class", y = "mean [Hg] ± 2 se") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))

BOSHgXlith.stats <- BOS.dat %>%
  group_by(lithgrp) %>%
  summarise(
    mean = mean(HgCOMP, na.rm = TRUE),
    median = median(HgCOMP, na.rm = TRUE), 
    var = var(HgCOMP, na.rm = TRUE),
    sd = sd(HgCOMP, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(HgCOMP, probs=0.975, na.rm = TRUE),
    lower = quantile(HgCOMP, probs=0.025, na.rm = TRUE),
    n = n()
  )

BOSHgXlith.sort <- na.omit(BOSHgXlith.stats %>% arrange(mean))
BOSHgXlith.sort

BOSHgXlith.sort$mean
BOSHgXlith.sort$se

ggplot(BOSHgXlith.sort, aes(x = reorder(lithgrp, mean), y = mean)) +
  geom_bar(stat='identity', fill = "steelblue") +
  geom_errorbar(aes(ymin = mean-(2*se), ymax = mean+(2*se)),
                width = 0.2) +
  labs(x = "lithology class", y = "mean [Hg] ± 2 se") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))


# Hg by soil type
TOS.HgXsoil <- sort(xtabs(TOS.dat$HgCOMP ~ TOS.dat$soils, na.rm=T) / table(TOS.dat$soils))
barplot(TOS.HgXsoil, xlab="soil type", ylab="mean TOS [Hg]", col="blue")
BOS.HgXsoil <- sort(xtabs(BOS.dat$HgCOMP ~ BOS.dat$soils, na.rm=T) / table(BOS.dat$soils))
barplot(BOS.HgXsoil, xlab="soil type", ylab="mean BOS [Hg]", col="blue")

table(TOS.dat$soils)
TOSHgXsoil.stats <- TOS.dat %>%
  group_by(soils) %>%
  summarise(
    mean = mean(HgCOMP, na.rm = TRUE),
    median = median(HgCOMP, na.rm = TRUE), 
    var = var(HgCOMP, na.rm = TRUE),
    sd = sd(HgCOMP, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(HgCOMP, probs=0.975, na.rm = TRUE),
    lower = quantile(HgCOMP, probs=0.025, na.rm = TRUE),
    n = n()
  )
TOSHgXsoil.sort <- na.omit(TOSHgXsoil.stats %>% arrange(mean))
TOSHgXsoil.sort

TOSHgXsoil.sort$mean
TOSHgXsoil.sort$se

soils.plot <- ggplot(TOSHgXsoil.sort, aes(x = reorder(soils, mean), y = mean)) +
  geom_bar(stat='identity', fill = "steelblue") +
  geom_errorbar(aes(ymin = mean-(2*se), ymax = mean+(2*se)),
                width = 0.2) +
  labs(x = "soil type", y = "mean [Hg] ± 2 se") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
soils.plot.flip <- soils.plot + coord_flip()
soils.plot.flip

BOSHgXsoil.stats <- BOS.dat %>%
  group_by(soils) %>%
  summarise(
    mean = mean(HgCOMP, na.rm = TRUE),
    median = median(HgCOMP, na.rm = TRUE), 
    var = var(HgCOMP, na.rm = TRUE),
    sd = sd(HgCOMP, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(HgCOMP, probs=0.975, na.rm = TRUE),
    lower = quantile(HgCOMP, probs=0.025, na.rm = TRUE),
    n = n()
  )
BOSHgXsoil.sort <- na.omit(BOSHgXsoil.stats %>% arrange(mean))
BOSHgXsoil.sort

BOSHgXsoil.sort$mean
BOSHgXsoil.sort$se

soils.plot <- ggplot(BOSHgXsoil.sort, aes(x = reorder(soils, mean), y = mean)) +
  geom_bar(stat='identity', fill = "steelblue") +
  geom_errorbar(aes(ymin = mean-(2*se), ymax = mean+(2*se)),
                width = 0.2) +
  labs(x = "soil type", y = "mean [Hg] ± 2 se") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
soils.plot.flip <- soils.plot + coord_flip()
soils.plot.flip



#######################################
## ratio of TOS:BOS Hg across biomes ##
#######################################
head(datmrg)
dim(datmrg)
SITEID.vec <- unique(datmrg$SITEID)
lenSITEID <- length(SITEID.vec)

Hgratio <- rep(NA, lenSITEID)
for (s in 1:lenSITEID) {
  site <- SITEID.vec[s]
  TOS.Hg <- datmrg$HgCOMP[datmrg$SITEID == site & datmrg$DEPTH == "TOS"]
  BOS.Hg <- datmrg$HgCOMP[datmrg$SITEID == site & datmrg$DEPTH == "BOS"]
  
  if (length(TOS.Hg) > 0 & length(BOS.Hg) > 0) {
    Hgratio[s] <- TOS.Hg / BOS.Hg
  } # end if
} # end s

## merge Hgratio with datmrg
Hgratio.dat <- data.frame(SITEID=SITEID.vec, Hgratio=Hgratio)
head(Hgratio.dat)
datmrg.ratio <- merge(datmrg, Hgratio.dat, by="SITEID", all.x=T)
head(datmrg.ratio)

# Hgratio stats by biome
HgratioXbiome.stats <- datmrg.ratio %>%
  group_by(biome) %>%
  summarise(
    mean = mean(Hgratio, na.rm = TRUE),
    median = median(Hgratio, na.rm = TRUE), 
    var = var(Hgratio, na.rm = TRUE),
    sd = sd(Hgratio, na.rm = TRUE),
    se = sd/sqrt(n()),
    sum = sum(Hgratio, na.rm = TRUE), 
    upper = quantile(Hgratio, probs=0.975, na.rm = TRUE),
    lower = quantile(Hgratio, probs=0.025, na.rm = TRUE),
    n = n()
  )
HgratioXbiome.sort <- na.omit(HgratioXbiome.stats %>% arrange(mean))
HgratioXbiome.sort
HgratioXbiome.sort$median

## make violin plot of Hgratio by biome
Hgratio.plot <- ggplot(datmrg.ratio, aes(x = biome, y = Hgratio)) +
  geom_violin(fill = "steelblue", alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  labs(x = "biome", y = "TOS:BOS Hg ratio") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    
    # Axis text (values)
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
Hgratio.plot


## create separate data.frames for each biome
Hgratio.DXS <- datmrg.ratio[datmrg.ratio$biome == "DXS", c("SITEID", "biome", "Hgratio")]
head(Hgratio.DXS)
Hgratio.MFW <- datmrg.ratio[datmrg.ratio$biome == "MFW", c("SITEID", "biome", "Hgratio")]
head(Hgratio.MFW)
Hgratio.TBMF <- datmrg.ratio[datmrg.ratio$biome == "TBMF", c("SITEID", "biome", "Hgratio")]
head(Hgratio.TBMF)
Hgratio.TGSS <- datmrg.ratio[datmrg.ratio$biome == "TGSS", c("SITEID", "biome", "Hgratio")]
head(Hgratio.TGSS)
Hgratio.TSGSS <- datmrg.ratio[datmrg.ratio$biome == "TSGSS", c("SITEID", "biome", "Hgratio")]
head(Hgratio.TSGSS)
Hgratio.TSMBF <- datmrg.ratio[datmrg.ratio$biome == "TSMBF", c("SITEID", "biome", "Hgratio")]
head(Hgratio.TSMBF)

# histograms of ratios by biome
# combine data
HgratioXbiome.dat <- rbind(na.omit(Hgratio.DXS[,c(2:3)]),
                            na.omit(Hgratio.MFW[,c(2:3)]),
                            na.omit(Hgratio.TBMF[,c(2:3)]),
                            na.omit(Hgratio.TGSS[,c(2:3)]),
                            na.omit(Hgratio.TSGSS[,c(2:3)]),
                            na.omit(Hgratio.TSMBF[,c(2:3)]))
head(HgratioXbiome.dat)
dim(HgratioXbiome.dat)

scaling_factors <- HgratioXbiome.dat %>%
  group_by(biome) %>%
  summarise(
    n = n(),
    bin_width = (max(Hgratio) - min(Hgratio)) / 20,  # assuming 20 bins
    scale_factor = n * bin_width,
    .groups = 'drop'
  )

# Join scaling factors back to original data
HgratioXbiome.dat.scaled <- HgratioXbiome.dat %>%
  left_join(scaling_factors, by = "biome")
head(HgratioXbiome.dat.scaled)

hist.plotXbiome <- ggplot(HgratioXbiome.dat.scaled, aes(x = Hgratio)) +
  geom_histogram(aes(y = after_stat(count)), bins = 20, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_density(aes(y = after_stat(count)), 
               color = "red", 
               size = 1.2, 
               alpha = 0.8) +
  facet_wrap(~ biome, scales = "free_y") +
  labs(x = "TOS [Hg]:BOS [Hg]",
    y = "frequency"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )
hist.plotXbiome


# export
write.csv(Hgratio.DXS, file="Hgratio.DXS.csv", row.names=F)
write.csv(Hgratio.MFW, file="Hgratio.MFW.csv", row.names=F)
write.csv(Hgratio.TBMF, file="Hgratio.TBMF.csv", row.names=F)
write.csv(Hgratio.TGSS, file="Hgratio.TGSS.csv", row.names=F)
write.csv(Hgratio.TSGSS, file="Hgratio.TSGSS.csv", row.names=F)
write.csv(Hgratio.TSMBF, file="Hgratio.TSMBF.csv", row.names=F)
biome.key


## resampling F-ratio test
# full dataset first
group_DXF <- na.omit(Hgratio.DXS$Hgratio)
group_MFW <- na.omit(Hgratio.MFW$Hgratio)
group_TBMF <- na.omit(Hgratio.TBMF$Hgratio)
group_TGSS <- na.omit(Hgratio.TGSS$Hgratio)
group_TSGSS <- na.omit(Hgratio.TSGSS$Hgratio)
group_TSMBF <- na.omit(Hgratio.TSMBF$Hgratio)

# combine data
all_data <- c(group_DXF, group_MFW, group_TBMF, group_TGSS, group_TSGSS, group_TSMBF)
groups <- factor(c(rep("DXF", length(group_DXF)),
                 rep("MFW", length(group_MFW)),
                 rep("TBMF", length(group_TBMF)),
                 rep("TGSS", length(group_TGSS)),
                 rep("TSGSS", length(group_TSGSS)),
                 rep("TSMBF", length(group_TSMBF))))

# means
grand_mean <- mean(all_data)
group_means <- tapply(all_data, groups, mean)

# SS between (SSB)
ssb <- sum(sapply(levels(groups), function(g) {
  n_g <- sum(groups == g)
  n_g * (group_means[g] - grand_mean)^2
}))

# SS within (SSW)
ssw <- sum(tapply(all_data, groups, function(x) sum((x - mean(x))^2)))

# degrees of freedom
df_between <- length(unique(groups)) - 1
df_within <- length(all_data) - length(unique(groups))

# mean squares
ms_between <- ssb / df_between
ms_within <- ssw / df_within

# F ratio
f_ratio <- ms_between / ms_within
f_ratio


## resample with replacement for each group
HgratioXbiome.sort$n
rsampn <- min(HgratioXbiome.sort$n) # number of resampled samples per group
Fiter <- 10000 # number of iterations
Fiter.div <- Fiter/10

f_ratio.obs.vec <- f_ratio.rnd.vec <- rndGTobs <- rep(NA, Fiter)
for (f in 1:Fiter) {

  # resample data
  DXF.rsmp <- sample(group_DXF, rsampn, replace = TRUE)
  MFW.rsmp <- sample(group_MFW, rsampn, replace = TRUE)
  TBMF.rsmp <- sample(group_TBMF, rsampn, replace = TRUE)
  TGSS.rsmp <- sample(group_TGSS, rsampn, replace = TRUE)
  TSGSS.rsmp <- sample(group_TSGSS, rsampn, replace = TRUE)
  TSMBF.rsmp <- sample(group_TSMBF, rsampn, replace = TRUE)

  # combine data
  all_data <- c(DXF.rsmp, MFW.rsmp, TBMF.rsmp, TGSS.rsmp, TSGSS.rsmp, TSMBF.rsmp)
  groups <- factor(c(rep("DXF", length(DXF.rsmp)),
                     rep("MFW", length(MFW.rsmp)),
                     rep("TBMF", length(TBMF.rsmp)),
                     rep("TGSS", length(TGSS.rsmp)),
                     rep("TSGSS", length(TSGSS.rsmp)),
                     rep("TSMBF", length(TSMBF.rsmp))))
  
  # means
  grand_mean <- mean(all_data)
  group_means <- tapply(all_data, groups, mean)
  
  # SS between (SSB)
  ssb <- sum(sapply(levels(groups), function(g) {
    n_g <- sum(groups == g)
    n_g * (group_means[g] - grand_mean)^2
  }))
  
  # SS within (SSW)
  ssw <- sum(tapply(all_data, groups, function(x) sum((x - mean(x))^2)))
  
  # degrees of freedom
  df_between <- length(unique(groups)) - 1
  df_within <- length(all_data) - length(unique(groups))
  
  # mean squares
  ms_between <- ssb / df_between
  ms_within <- ssw / df_within
  
  # observed F ratio
  f_ratio.obs.vec[f] <- ms_between / ms_within
  
  # randomise values among groups
  all_data.rnd <- sample(all_data, length(all_data), replace = FALSE)
  
  # means
  grand_mean.rnd <- mean(all_data.rnd)
  group_means.rnd <- tapply(all_data.rnd, groups, mean)
  
  # SS between (SSB)
  ssb.rnd <- sum(sapply(levels(groups), function(g) {
    n_g <- sum(groups == g)
    n_g * (group_means.rnd[g] - grand_mean.rnd)^2
  }))
  
  # SS within (SSW)
  ssw.rnd <- sum(tapply(all_data.rnd, groups, function(x) sum((x - mean(x))^2)))
  
  # degrees of freedom
  df_between.rnd <- length(unique(groups)) - 1
  df_within.rnd <- length(all_data.rnd) - length(unique(groups))
  
  # mean squares
  ms_between.rnd <- ssb.rnd / df_between.rnd
  ms_within.rnd <- ssw.rnd / df_within.rnd
  
  # observed F ratio
  f_ratio.rnd.vec[f] <- ms_between.rnd / ms_within.rnd
  
  # F ratio compare
  rndGTobs[f] <- ifelse(f_ratio.rnd.vec[f] >= f_ratio.obs.vec[f], 1, 0)
  
  # counter
  if (f %% Fiter.div == 0) {
    print(paste("iteration: ", f))
  }
} # end f loop

# histogram of F ratios
hist(f_ratio.vec, breaks=50, main="F ratio distribution", xlab="F ratio", col="lightblue", border="black")
median(f_ratio.vec, na.rm=T)
quantile(f_ratio.vec, probs=c(0.025, 0.975), na.rm=T)

# p random
p.rnd <- sum(rndGTobs, na.rm=T) / Fiter
print(paste("p random = ", p.rnd))


## repeat without tropical forest TSMBF (low sample size in this group)
HgratioXbiome.sort$n
rsampn <- min(HgratioXbiome.sort$n[-which(HgratioXbiome.sort$n == min(HgratioXbiome.sort$n))]) # number of resampled samples per group
Fiter <- 10000 # number of iterations
Fiter.div <- Fiter/10

f_ratio.obs.vec <- f_ratio.rnd.vec <- rndGTobs <- rep(NA, Fiter)
for (f in 1:Fiter) {
  
  # resample data
  DXF.rsmp <- sample(group_DXF, rsampn, replace = TRUE)
  MFW.rsmp <- sample(group_MFW, rsampn, replace = TRUE)
  TBMF.rsmp <- sample(group_TBMF, rsampn, replace = TRUE)
  TGSS.rsmp <- sample(group_TGSS, rsampn, replace = TRUE)
  TSGSS.rsmp <- sample(group_TSGSS, rsampn, replace = TRUE)

  # combine data
  all_data <- c(DXF.rsmp, MFW.rsmp, TBMF.rsmp, TGSS.rsmp, TSGSS.rsmp)
  groups <- factor(c(rep("DXF", length(DXF.rsmp)),
                     rep("MFW", length(MFW.rsmp)),
                     rep("TBMF", length(TBMF.rsmp)),
                     rep("TGSS", length(TGSS.rsmp)),
                     rep("TSGSS", length(TSGSS.rsmp))))
  
  # means
  grand_mean <- mean(all_data)
  group_means <- tapply(all_data, groups, mean)
  
  # SS between (SSB)
  ssb <- sum(sapply(levels(groups), function(g) {
    n_g <- sum(groups == g)
    n_g * (group_means[g] - grand_mean)^2
  }))
  
  # SS within (SSW)
  ssw <- sum(tapply(all_data, groups, function(x) sum((x - mean(x))^2)))
  
  # degrees of freedom
  df_between <- length(unique(groups)) - 1
  df_within <- length(all_data) - length(unique(groups))
  
  # mean squares
  ms_between <- ssb / df_between
  ms_within <- ssw / df_within
  
  # observed F ratio
  f_ratio.obs.vec[f] <- ms_between / ms_within
  
  # randomise values among groups
  all_data.rnd <- sample(all_data, length(all_data), replace = FALSE)
  
  # means
  grand_mean.rnd <- mean(all_data.rnd)
  group_means.rnd <- tapply(all_data.rnd, groups, mean)
  
  # SS between (SSB)
  ssb.rnd <- sum(sapply(levels(groups), function(g) {
    n_g <- sum(groups == g)
    n_g * (group_means.rnd[g] - grand_mean.rnd)^2
  }))
  
  # SS within (SSW)
  ssw.rnd <- sum(tapply(all_data.rnd, groups, function(x) sum((x - mean(x))^2)))
  
  # degrees of freedom
  df_between.rnd <- length(unique(groups)) - 1
  df_within.rnd <- length(all_data.rnd) - length(unique(groups))
  
  # mean squares
  ms_between.rnd <- ssb.rnd / df_between.rnd
  ms_within.rnd <- ssw.rnd / df_within.rnd
  
  # observed F ratio
  f_ratio.rnd.vec[f] <- ms_between.rnd / ms_within.rnd
  
  # F ratio compare
  rndGTobs[f] <- ifelse(f_ratio.rnd.vec[f] >= f_ratio.obs.vec[f], 1, 0)
  
  # counter
  if (f %% Fiter.div == 0) {
    print(paste("iteration: ", f))
  }
} # end f loop

# histogram of F ratios
hist(f_ratio.vec, breaks=50, main="F ratio distribution", xlab="F ratio", col="lightblue", border="black")
median(f_ratio.vec, na.rm=T)
quantile(f_ratio.vec, probs=c(0.025, 0.975), na.rm=T)

# p random
p.rnd <- sum(rndGTobs, na.rm=T) / Fiter
print(paste("p random = ", p.rnd))



## create dataset for testing relationships
# TOS
colnames(TOS.dat)
TOS.test1 <- data.frame(SITEID=TOS.dat$SITEID, LAT=TOS.dat$LAT, LON=TOS.dat$LON, lHg=log10(TOS.dat$HgCOMP), lgs=log10(TOS.dat$gs.wmn),
                       lgtclay=logit(TOS.dat$CLAY/100), lgtsoilclay=logit(TOS.dat$soilclay/100),
                       lgtsilt=logit(TOS.dat$SILT/100), lgtsoilsilt=logit(TOS.dat$soilsilt/100),  
                       pH=TOS.dat$pH, pH15=TOS.dat$pH15, soilpH=TOS.dat$soilpH, lelecond=log10(TOS.dat$EC15), lLOI=log10(TOS.dat$LOI),
                       lrain=log10(TOS.dat$rain), prescott=TOS.dat$prescott, soilH20=TOS.dat$soilH20,
                       llai=log10(TOS.dat$lai), lgpp=log10(TOS.dat$gpp),
                       soillN=log10(TOS.dat$soilN), soillP=log10(TOS.dat$soilP),
                       lAl=log10(TOS.dat$Al), lFe=log10(TOS.dat$Fe), lMn=log10(TOS.dat$Mn), lCu=log10(TOS.dat$CuICPMS),
                       lZn=log10(TOS.dat$ZnICPMS), lPb=log10(TOS.dat$PbICPMS), lSb=log10(TOS.dat$SbICPMS),
                       lNi=log10(TOS.dat$NiCPMS), lV=log10(TOS.dat$V), lKThU=log10(TOS.dat$KThU),
                       lucat=TOS.dat$lucat, biome=TOS.dat$biome, geol=TOS.dat$lithgrp, soils=TOS.dat$soils,
                       ferric2=TOS.dat$ferricpc2, Alox=TOS.dat$Alox, 
                       Feox=TOS.dat$Feox, Pox=TOS.dat$Pox)

head(TOS.test1)
dim(TOS.test1)
TOS.test <- distinct(TOS.test1)
dim(TOS.test)
head(TOS.test)

## any NAs?
length(which(is.na(TOS.test$lHg)==T)) # good
length(which(is.na(TOS.test$lgs)==T)) # good
length(which(is.na(TOS.test$lgtclay)==T)) # good
length(which(is.na(TOS.test$lgtsoilclay)==T)) # 7 missing
length(which(is.na(TOS.test$lgtsilt)==T)) # good
length(which(is.na(TOS.test$lgtsoilsilt)==T)) # 7 missing
length(which(is.na(TOS.test$pH15)==T)) # good
length(which(is.na(TOS.test$soilpH)==T)) # 5 missing
length(which(is.na(TOS.test$pH)==T)) # good
length(which(is.na(TOS.test$lelecond)==T)) # good
length(which(is.na(TOS.test$lLOI)==T)) # 3 missing
length(which(is.na(TOS.test$lrain)==T)) # 16 missing
length(which(is.na(TOS.test$prescott)==T)) # good
length(which(is.na(TOS.test$soilH20)==T)) # 1 missing
length(which(is.na(TOS.test$llai)==T)) # 1 missing
length(which(is.na(TOS.test$lgpp)==T)) # 1 missing
length(which(is.na(TOS.test$soillN)==T)) # 5 missing
length(which(is.na(TOS.test$soillP)==T)) # 5 missing
length(which(is.na(TOS.test$lAl)==T)) # 3 missing
length(which(is.na(TOS.test$lFe)==T)) # 3 missing
length(which(is.na(TOS.test$lMn)==T)) # 3 missing
length(which(is.na(TOS.test$lCu)==T)) # 3 missing
length(which(is.na(TOS.test$lZn)==T)) # 3 missing
length(which(is.na(TOS.test$lPb)==T)) # 3 missing
length(which(is.na(TOS.test$lSb)==T)) # 3 missing
length(which(is.na(TOS.test$lNi)==T)) # 3 missing
length(which(is.na(TOS.test$lV)==T)) # 3 missing
length(which(is.na(TOS.test$lKThU)==T)) # good
length(which(is.na(TOS.test$ferric2)==T)) # good
length(which(is.na(TOS.test$Alox)==T)) # good
length(which(is.na(TOS.test$Feox)==T)) # good
length(which(is.na(TOS.test$Pox)==T)) # good

## any infinites?
length(which(is.infinite(TOS.test$lHg)==T))
length(which(is.infinite(TOS.test$lgs)==T))
length(which(is.infinite(TOS.test$pH15)==T))
length(which(is.infinite(TOS.test$lgtclay)==T)) # 3 infinites
length(which(is.infinite(TOS.test$lgtsilt)==T)) # 1 infinite
length(which(is.infinite(TOS.test$lelecond)==T))
length(which(is.infinite(TOS.test$lLOI)==T))
length(which(is.infinite(TOS.test$lrain)==T))
length(which(is.infinite(TOS.test$prescott)==T))
length(which(is.infinite(TOS.test$soilH20)==T))
length(which(is.infinite(TOS.test$llai)==T)) # 1 infinite
length(which(is.infinite(TOS.test$lgpp)==T)) # 3 infinites
length(which(is.infinite(TOS.test$soillN)==T))
length(which(is.infinite(TOS.test$soillP)==T))
length(which(is.infinite(TOS.test$lAl)==T))
length(which(is.infinite(TOS.test$lFe)==T))
length(which(is.infinite(TOS.test$lMn)==T))
length(which(is.infinite(TOS.test$lCu)==T))
length(which(is.infinite(TOS.test$lZn)==T))
length(which(is.infinite(TOS.test$lPb)==T))
length(which(is.infinite(TOS.test$lSb)==T))
length(which(is.infinite(TOS.test$lNi)==T))
length(which(is.infinite(TOS.test$lV)==T))
length(which(is.infinite(TOS.test$lKThU)==T))
length(which(is.infinite(TOS.test$ferric2)==T))
length(which(is.infinite(TOS.test$Alox)==T))
length(which(is.infinite(TOS.test$Feox)==T))

testpts <- vect(cbind(TOS.test$LON, TOS.test$LAT), crs="+proj=longlat")
terra::plot(testpts)

# histograms
hist(TOS.test$lHg, xlab="log10 [Hg]", main="")
hist(TOS.test$lgs, xlab="log10 mean grain size", main="")
hist(TOS.test$pH, xlab="pH", main="")
hist(TOS.test$pH15, xlab="pH", main="")
hist(TOS.test$soilpH, xlab="soil pH", main="")
hist(TOS.test$lelecond, xlab="log10 soil electroconductivity", main="")
hist(TOS.test$lgtclay, xlab="logit proportion clay", main="")
hist(TOS.test$lgtsoilclay, xlab="logit proportion soil clay", main="")
hist(TOS.test$lgtsilt, xlab="logit proportion silt", main="")
hist(TOS.test$lgtsoilsilt, xlab="logit proportion soil silt", main="")
hist(TOS.test$lLOI, xlab="log10 LOI", main="")
hist(TOS.test$lrain, xlab="log10 annual rainfall", main="")
hist(TOS.test$prescott, xlab="Prescott index", main="")
hist(TOS.test$soilH20, xlab="soil H2O availability", main="")
hist(TOS.test$llai, xlab="log10 leaf area index", main="")
hist(TOS.test$lgpp, xlab="log10 vegetation C uptake (GPP)", main="")
hist(TOS.test$soillN, xlab="log10 soil N", main="")
hist(TOS.test$soillP, xlab="log10 soil P", main="")
hist(TOS.test$lAl, xlab="log10 [Al]", main="")
hist(TOS.test$lFe, xlab="log10 [Fe]", main="")
hist(TOS.test$lMn, xlab="log10 [Mn]", main="")
hist(TOS.test$lCu, xlab="log10 [Cu]", main="")
hist(TOS.test$lZn, xlab="log10 [Zn]", main="")
hist(TOS.test$lPb, xlab="log10 [Pb]", main="")
hist(TOS.test$lSb, xlab="log10 [Sb]", main="")
hist(TOS.test$lNi, xlab="log10 [Ni]", main="")
hist(TOS.test$lV, xlab="log10 [V]", main="")
hist(TOS.test$lKThU, xlab="log10 [KThU]", main="")
hist(TOS.test$ferric2, xlab="ferric PC2", main="")
hist(TOS.test$Alox, xlab="Al oxide", main="")
hist(TOS.test$Feox, xlab="Fe oxide", main="")
hist(TOS.test$Pox, xlab="P oxide", main="")

dim(TOS.test)

# replace infinites with NA
TOS.test$lgtclay[(which(is.infinite(TOS.test$lgtclay)==T))] <- NA
TOS.test$lgtsilt[(which(is.infinite(TOS.test$lgtsilt)==T))] <- NA
TOS.test$llai[(which(is.infinite(TOS.test$llai)==T))] <- NA
TOS.test$lgpp[(which(is.infinite(TOS.test$lgpp)==T))] <- NA

## full dataset correlation matrix
TOS.cor <- na.omit(TOS.test[,c("lHg","lgs","lgtclay","lgtsilt","pH15","lelecond","lLOI","prescott","soilH20","llai","lgpp","soillN","soillP",
                               "lAl","lFe","lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric2","Alox","Feox","Pox")])
class(TOS.cor)
dim(TOS.cor)
head(TOS.cor)
round(cor(TOS.cor, method="spearman"), 2)

## variance inflation
VIF.TOS <- usdm::vif(TOS.cor)
VIF.TOS
VIF.TOS[which(VIF.TOS$VIF > 4),]

# remove lFe, lgpp, lgtsilt, lAl

## limit variables after VIF inspection
TOS.test.brt <- TOS.cor[,c("lHg","lgs","lgtclay","pH15","lelecond","lLOI","prescott","soilH20","llai","soillN","soillP",
                           "lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric2","Alox","Feox","Pox")]
usdm::vif(TOS.test.brt)

dim(TOS.test.brt)
colnames(TOS.test.brt)

## boosted regression tree
TOS.brt <- gbm.step(TOS.test.brt, gbm.x = attr(TOS.test.brt, "names")[c(2:length(colnames(TOS.test.brt)))],
                    gbm.y = attr(TOS.test.brt, "names")[1], family="gaussian", max.trees=100000,
                    tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.75,
                    tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(TOS.brt)
barplot(summary(TOS.brt)$rel.inf, names.arg = summary(TOS.brt)$var, xlab="relative influence", ylab="", col="blue")
TOS.brt.summ <- summary(TOS.brt)

TOS.brt.plot <- ggplot(TOS.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
TOS.brt.plot.flip <- TOS.brt.plot + coord_flip()
TOS.brt.plot.flip


gbm.plot(TOS.brt)

gbm.plot(TOS.brt, smooth=T, rug=T, n.plots=12, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]")

gbm.plot(TOS.brt, variable.no=3)
colnames(TOS.test.brt)

gbm.plot(TOS.brt, variable.no=5, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]", x.label="log10 LOI", plot.layout=c(1,1))
gbm.plot(TOS.brt, variable.no=19, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]", x.label="ferric PC2", plot.layout=c(1,1))
gbm.plot(TOS.brt, variable.no=14, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]", x.label="log10 Pb", plot.layout=c(1,1))
gbm.plot(TOS.brt, variable.no=9, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]", x.label="log10 soil nitrogen", plot.layout=c(1,1))
gbm.plot(TOS.brt, variable.no=13, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]", x.label="log10 Zn", plot.layout=c(1,1))
gbm.plot(TOS.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]", x.label="pH1:5", plot.layout=c(1,1))

TOS.brt.partial.deps <- list()
TOS.pred.names <- names(TOS.test.brt)
TOS.pred.names <- TOS.pred.names[-1]
eq.sp.pts <- 100

for(i in seq_along(TOS.pred.names)) {
  # Use gbm.plot to get partial dependency values
  # Set plot.layout = c(1,1) to prevent automatic plotting
  pd_data <- plot.gbm(TOS.brt, i.var=i, continuous.resolution = eq.sp.pts, return.grid=T)
  TOS.brt.partial.deps[[TOS.pred.names[i]]] <- pd_data
  write.csv(TOS.brt.partial.deps[[i]], file = paste0("TOS_", TOS.pred.names[i], "_partial_deps.csv"))
}

TOS.brt.CV.cor <- 100 * TOS.brt$cv.statistics$correlation.mean
TOS.brt.CV.cor.se <- 100 * TOS.brt$cv.statistics$correlation.se
print(c(TOS.brt.CV.cor, TOS.brt.CV.cor.se))


#########################
## BOS
colnames(BOS.dat)
BOS.test1 <- data.frame(SITEID=BOS.dat$SITEID, LAT=BOS.dat$LAT, LON=BOS.dat$LON, lHg=log10(BOS.dat$HgCOMP), lgs=log10(BOS.dat$gs.wmn),
                        lgtclay=logit(BOS.dat$CLAY/100), lgtsoilclay=logit(BOS.dat$soilclay/100),
                        lgtsilt=logit(BOS.dat$SILT/100), lgtsoilsilt=logit(BOS.dat$soilsilt/100),  
                        pH=BOS.dat$pH, pH15=BOS.dat$pH15, soilpH=BOS.dat$soilpH, lelecond=log10(BOS.dat$EC15), lLOI=log10(BOS.dat$LOI),
                        lrain=log10(BOS.dat$rain), prescott=BOS.dat$prescott, soilH20=BOS.dat$soilH20,
                        llai=log10(BOS.dat$lai), lgpp=log10(BOS.dat$gpp),
                        soillN=log10(BOS.dat$soilN), soillP=log10(BOS.dat$soilP),
                        lAl=log10(BOS.dat$Al), lFe=log10(BOS.dat$Fe), lMn=log10(BOS.dat$Mn), lCu=log10(BOS.dat$CuICPMS),
                        lZn=log10(BOS.dat$ZnICPMS), lPb=log10(BOS.dat$PbICPMS), lSb=log10(BOS.dat$SbICPMS),
                        lNi=log10(BOS.dat$NiCPMS), lV=log10(BOS.dat$V), lKThU=log10(BOS.dat$KThU),
                        lucat=BOS.dat$lucat, biome=BOS.dat$biome, geol=BOS.dat$lithgrp, soils=BOS.dat$soils,
                        ferric2=BOS.dat$ferricpc2, Alox=BOS.dat$Alox, 
                        Feox=BOS.dat$Feox, Pox=BOS.dat$Pox)


head(BOS.test1)
dim(BOS.test1)
BOS.test <- distinct(BOS.test1)
dim(BOS.test)
head(BOS.test)

## any NAs?
length(which(is.na(BOS.test$lHg)==T)) # good
length(which(is.na(BOS.test$lgs)==T)) # good
length(which(is.na(BOS.test$lgtclay)==T)) # good
length(which(is.na(BOS.test$lgtsoilclay)==T)) # 7 missing
length(which(is.na(BOS.test$lgtsilt)==T)) # good
length(which(is.na(BOS.test$lgtsoilsilt)==T)) # 7 missing
length(which(is.na(BOS.test$pH15)==T)) # good
length(which(is.na(BOS.test$soilpH)==T)) # 5 missing
length(which(is.na(BOS.test$pH)==T)) # good
length(which(is.na(BOS.test$lelecond)==T)) # good
length(which(is.na(BOS.test$lLOI)==T)) # 3 missing
length(which(is.na(BOS.test$lrain)==T)) # 15 missing
length(which(is.na(BOS.test$prescott)==T)) # good
length(which(is.na(BOS.test$soilH20)==T)) # 1 missing
length(which(is.na(BOS.test$llai)==T)) # 1 missing
length(which(is.na(BOS.test$lgpp)==T)) # 1 missing
length(which(is.na(BOS.test$soillN)==T)) # 5 missing
length(which(is.na(BOS.test$soillP)==T)) # 5 missing
length(which(is.na(BOS.test$lAl)==T)) # 1 missing
length(which(is.na(BOS.test$lFe)==T)) # 1 missing
length(which(is.na(BOS.test$lMn)==T)) # 1 missing
length(which(is.na(BOS.test$lCu)==T)) # 1 missing
length(which(is.na(BOS.test$lZn)==T)) # 1 missing
length(which(is.na(BOS.test$lPb)==T)) # 1 missing
length(which(is.na(BOS.test$lSb)==T)) # 1 missing
length(which(is.na(BOS.test$lNi)==T)) # 1 missing
length(which(is.na(BOS.test$lV)==T)) # 1 missing
length(which(is.na(BOS.test$lKThU)==T)) # good
length(which(is.na(BOS.test$ferric2)==T)) # good
length(which(is.na(BOS.test$Alox)==T)) # good
length(which(is.na(BOS.test$Feox)==T)) # good
length(which(is.na(BOS.test$Pox)==T)) # good

## any infinites?
length(which(is.infinite(BOS.test$lHg)==T))
length(which(is.infinite(BOS.test$lgs)==T))
length(which(is.infinite(BOS.test$pH15)==T))
length(which(is.infinite(BOS.test$lgtclay)==T))
length(which(is.infinite(BOS.test$lgtsilt)==T))
length(which(is.infinite(BOS.test$lelecond)==T))
length(which(is.infinite(BOS.test$lLOI)==T))
length(which(is.infinite(BOS.test$lrain)==T))
length(which(is.infinite(BOS.test$prescott)==T))
length(which(is.infinite(BOS.test$soilH20)==T))
length(which(is.infinite(BOS.test$llai)==T)) # 1 infinite
length(which(is.infinite(BOS.test$lgpp)==T)) # 3 infinites
length(which(is.infinite(BOS.test$soillN)==T))
length(which(is.infinite(BOS.test$soillP)==T))
length(which(is.infinite(BOS.test$lAl)==T))
length(which(is.infinite(BOS.test$lFe)==T))
length(which(is.infinite(BOS.test$lMn)==T))
length(which(is.infinite(BOS.test$lCu)==T))
length(which(is.infinite(BOS.test$lZn)==T))
length(which(is.infinite(BOS.test$lPb)==T))
length(which(is.infinite(BOS.test$lSb)==T))
length(which(is.infinite(BOS.test$lNi)==T))
length(which(is.infinite(BOS.test$lV)==T))
length(which(is.infinite(BOS.test$lKThU)==T))
length(which(is.infinite(BOS.test$ferric2)==T))
length(which(is.infinite(BOS.test$Alox)==T))
length(which(is.infinite(BOS.test$Feox)==T))

testpts <- vect(cbind(BOS.test$LON, BOS.test$LAT), crs="+proj=longlat")
terra::plot(testpts)

# histograms
hist(BOS.test$lHg, xlab="log10 [Hg]", main="")
hist(BOS.test$lgs, xlab="log10 mean grain size", main="")
hist(BOS.test$pH, xlab="pH", main="")
hist(BOS.test$pH15, xlab="pH", main="")
hist(BOS.test$soilpH, xlab="soil pH", main="")
hist(BOS.test$lelecond, xlab="log10 soil electroconductivity", main="")
hist(BOS.test$lgtclay, xlab="logit proportion clay", main="")
hist(BOS.test$lgtsoilclay, xlab="logit proportion soil clay", main="")
hist(BOS.test$lgtsilt, xlab="logit proportion silt", main="")
hist(BOS.test$lgtsoilsilt, xlab="logit proportion soil silt", main="")
hist(BOS.test$lLOI, xlab="log10 LOI", main="")
hist(BOS.test$lrain, xlab="log10 annual rainfall", main="")
hist(BOS.test$prescott, xlab="Prescott index", main="")
hist(BOS.test$soilH20, xlab="soil H2O availability", main="")
hist(BOS.test$llai, xlab="log10 leaf area index", main="")
hist(BOS.test$lgpp, xlab="log10 vegetation C uptake (GPP)", main="")
hist(BOS.test$soillN, xlab="log10 soil N", main="")
hist(BOS.test$soillP, xlab="log10 soil P", main="")
hist(BOS.test$lAl, xlab="log10 [Al]", main="")
hist(BOS.test$lFe, xlab="log10 [Fe]", main="")
hist(BOS.test$lMn, xlab="log10 [Mn]", main="")
hist(BOS.test$lCu, xlab="log10 [Cu]", main="")
hist(BOS.test$lZn, xlab="log10 [Zn]", main="")
hist(BOS.test$lPb, xlab="log10 [Pb]", main="")
hist(BOS.test$lSb, xlab="log10 [Sb]", main="")
hist(BOS.test$lNi, xlab="log10 [Ni]", main="")
hist(BOS.test$lV, xlab="log10 [V]", main="")
hist(BOS.test$lKThU, xlab="log10 [KThU]", main="")
hist(BOS.test$ferric2, xlab="ferric PC2", main="")
hist(BOS.test$Alox, xlab="Al oxide", main="")
hist(BOS.test$Feox, xlab="Fe oxide", main="")
hist(BOS.test$Pox, xlab="P oxide", main="")

dim(BOS.test)

# replace infinites with NA
BOS.test$llai[(which(is.infinite(BOS.test$llai)==T))] <- NA
BOS.test$lgpp[(which(is.infinite(BOS.test$lgpp)==T))] <- NA

## correlation matrix
BOS.cor <- na.omit(BOS.test[,c("lHg","lgs","lgtclay","lgtsilt","pH15","lelecond","lLOI","prescott","soilH20","llai","lgpp","soillN","soillP",
                               "lAl","lFe","lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric2","Alox","Feox","Pox")])
class(BOS.cor)
dim(BOS.cor)
head(BOS.cor)
round(cor(BOS.cor, method="spearman"), 2)

## variance inflation
usdm::vif(BOS.cor)
usdm::vif(BOS.cor[,c("lHg","lgs","lgtclay","lgtsilt","pH15","lelecond","lLOI","prescott","soilH20","llai","lgpp","soillN","soillP",
                     "lAl","lFe","lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric2","Alox","Feox","Pox")])

# remove lgpp, lgtsilt, lFe
usdm::vif(BOS.cor[,c("lHg","lgs","lgtclay","pH15","lelecond","lLOI","prescott","soilH20","llai","soillN","soillP",
                     "lAl","lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric2","Alox","Feox","Pox")])

## limit variables after VIF inspection
BOS.test.brt <- BOS.cor[,c("lHg","lgs","lgtclay","pH15","lelecond","lLOI","prescott","soilH20","llai","soillN","soillP",
                           "lAl","lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric2","Alox","Feox","Pox")]
dim(BOS.test.brt)
colnames(BOS.test.brt)

## boosted regression tree
BOS.brt <- gbm.step(BOS.test.brt, gbm.x = attr(BOS.test.brt, "names")[c(2:length(colnames(BOS.test.brt)))],
                    gbm.y = attr(BOS.test.brt, "names")[1], family="gaussian", max.trees=100000,
                    tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.75,
                    tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(BOS.brt)
barplot(summary(BOS.brt)$rel.inf, names.arg = summary(BOS.brt)$var, xlab="relative influence", ylab="", col="blue")

gbm.plot(BOS.brt, smooth=T, rug=T, n.plots=12, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]")

BOS.brt.fitted <- predict.gbm(BOS.brt, BOS.test.brt, n.trees=BOS.brt$gbm.call$best.trees, type='response')
plot(BOS.test.brt$lHg, BOS.brt.fitted, xlab = "observed ", ylab = "fitted")
abline(0, 1, col = "red")
plot(BOS.brt, n.trees = BOS.brt$gbm.call$best.trees,
     write.title = FALSE)
BOS.brt.results <- data.frame(
  obs = BOS.test.brt$lHg,
  fit = BOS.brt.fitted
)

BOS.brt.partial.deps <- list()
BOS.pred.names <- names(BOS.test.brt)
BOS.pred.names <- BOS.pred.names[-1]

for(i in seq_along(BOS.pred.names)) {
  # Use gbm.plot to get partial dependency values
  # Set plot.layout = c(1,1) to prevent automatic plotting
  pd_data <- plot.gbm(BOS.brt, i.var=i, continuous.resolution = eq.sp.pts, return.grid=T)
  BOS.brt.partial.deps[[BOS.pred.names[i]]] <- pd_data
  write.csv(BOS.brt.partial.deps[[i]], file = paste0("BOS_", BOS.pred.names[i], "_partial_deps.csv"))
}

gbm.plot(BOS.brt)
gbm.plot.fits(BOS.brt)

BOS.brt.CV.cor <- 100 * BOS.brt$cv.statistics$correlation.mean
BOS.brt.CV.cor.se <- 100 * BOS.brt$cv.statistics$correlation.se
print(c(BOS.brt.CV.cor, BOS.brt.CV.cor.se))



###################################################### 
## multiscale spatial analysis of multivariate data ##
######################################################
## see https://cran.r-universe.dev/adespatial/doc/tutorial.html

# remove NAs
colnames(TOS.test)
TOS.noNA <- na.omit(TOS.test[,c("LON","LAT","lHg","lgs","lgtclay","pH15","lelecond","lLOI","prescott","soilH20","llai","soillN","soillP",
                                "lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric2","Alox","Feox","Pox")])
head(TOS.noNA)
dim(TOS.noNA)

TOS.xy <- data.frame(x=TOS.noNA$LON, y=TOS.noNA$LAT)
head(TOS.xy)
mxy <- as.matrix(TOS.xy)
rownames(mxy) <- NULL

# Aus map as a SpatialPolygons object
aus.sf <- sf::st_read("aus.shp")

s.label(mxy, ppoint.pch = 15, ppoint.col = "darkseagreen4")

# Gabriel spatial weighting matrix
nbgab <- graph2nb(gabrielneigh(mxy, nnmult=20), sym = TRUE)
nb2listw(nbgab)
distgab <- nbdists(nbgab, mxy)
nbgab[[1]]
fdist <- lapply(distgab, function(x) 1 - x/max(dist(mxy))) # spatial weights
listwgab <- nb2listw(nbgab, glist = fdist) # spatial weighting matrix
listwgab
names(listwgab)
listwgab$neighbours[[1]]
listwgab$weights[[1]]
print(listw2mat(listwgab)[1:10, 1:10], digits = 3)

# environmental data
TOS.env <- TOS.noNA[,3:dim(TOS.noNA)[2]]
head(TOS.env)



########################
## distribution model ##
########################

# store boundaries in a single extent object
geogr.extent <- ext(soilclay)

# TOS
# prepare data
TOS.dm.dat <- na.omit(TOS.test[,c("LON","LAT","lHg","lgs","lgtclay","pH15","lelecond","lLOI","prescott","soilH20","llai","soillN","soillP",
                                  "lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric2","Alox","Feox","Pox")])

TOSHgDat <- data.frame(
  x = TOS.dm.dat$LON,
  y = TOS.dm.dat$LAT,
  Hg = TOS.dm.dat$lHg,
  gs = TOS.dm.dat$lgs,
  clay = TOS.dm.dat$lgtclay,
  pH = TOS.dm.dat$pH15,
  Prescott = TOS.dm.dat$prescott,
  soilH20 = TOS.dm.dat$soilH20,
  lai = TOS.dm.dat$llai,
  soilN = TOS.dm.dat$soillN,
  soilP = TOS.dm.dat$soillP,
  KThU = TOS.dm.dat$lKThU,
  ferric2=TOS.dm.dat$ferric2,
  Alox=TOS.dm.dat$Alox,
  Feox=TOS.dm.dat$Feox,
  Pox=TOS.dm.dat$Pox
)

# convert to spatial points
TOSHgSpat <- TOSHgDat
coordinates(TOSHgSpat) <- ~x+y

# prepare training data
TOStraining_data <- data.frame(
  Hg = TOSHgSpat$Hg,
  # environmental predictors
  clay = TOSHgSpat$clay,
  pH = TOSHgSpat$pH,
  Prescott = TOSHgSpat$Prescott,
  soilH20 = TOSHgSpat$soilH20,
  lai = TOSHgSpat$lai,
  soilN = TOSHgSpat$soilN,
  soilP = TOSHgSpat$soilP,
  KThU = TOSHgSpat$KThU,
  ferric2=TOSHgSpat$ferric2,
  Alox=TOSHgSpat$Alox,
  Feox=TOSHgSpat$Feox,
  Pox=TOSHgSpat$Pox
)

crs(TOSHgSpat) <- crsUse

# choose size of spatial blocks for spatial autocorrelation
TOSblocksize <- cv_spatial_autocor(x=TOSHgSpat, column='Hg', plot=T, deg_to_metre=111139) 
TOSblocksize$range

# spatial blocks for cross-validation
TOSspatial_blocks <- cv_spatial(x=TOSHgSpat, k=5, hexagon=F, size=200000)
TOSspatial_blocks$folds_list

# custom indices for spatial cross-validation
TOSspatial_indices <- TOSspatial_blocks$folds_list

# define spatial cross-validation control
TOSsactrl <- trainControl(
  method = "cv",
  classProbs = F,
  number = 5,
  index = TOSspatial_indices,
  savePredictions = TRUE
)

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# train random forest model (no control for spatial autocorrelation)
TOSrf_model <- train(
  Hg ~ .,
  data = TOStraining_data,
  method = "rf",
  trControl = trainControl(method = "cv", number = 10)
)

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

TOSrf_model$results
TOSrf_model$finalModel$importance
str(TOSrf_model)


# create prediction raster stack of environmental variables
envVars <- c(clay.rsmp,pH.rsmp,Prescott.rsmp,soilH20.rsmp,lai.rsmp,soilN.rsmp,soilP.rsmp,KThU.rsmp,
             ferric2.rsmp,Alox.rsmp,Feox.rsmp,Pox.rsmp)
names(envVars) <- c("clay","pH","Prescott","soilH20","lai","soilN","soilP","KThU","ferric2","Alox","Feox","Pox")
names(envVars)
plot(envVars)
envVars

# spatial predictions
TOSHgPred <- terra::predict(envVars, TOSrf_model, na.rm=T)
plot(TOSHgPred)
points(TOSHgDat, pch=20, col="black", cex=0.5)
TOSHgPred.rf.bt <- TOSHgPred^10
plot(TOSHgPred.rf.bt)
writeRaster(TOSHgPred.rf.bt, "TOSHgPredRFbt.tif", overwrite=T)


# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# model cross-validation
TOScv <- train(
  Hg ~ .,
  data = TOStraining_data,
  method = "rf",
  trControl = trainControl(
    method = "cv",
    number = 10,
    savePredictions = TRUE
  )
)

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# performance metrics
TOSrmse <- sqrt(mean((TOScv$pred$obs - TOScv$pred$pred)^2))
TOSrmse
TOSr2 <- cor(TOScv$pred$obs, TOScv$pred$pred)^2
TOSr2

# variable importance
varImp(TOSrf_model)

# plot variable importance
plot(varImp(TOSrf_model))


# BOS
# prepare data
BOS.dm.dat <- na.omit(BOS.test[,c("LON","LAT","lHg","lgs","lgtclay","pH15","lelecond","lLOI","prescott","soilH20","llai","soillN","soillP",
                                  "lAl","lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric2","Alox","Feox","Pox")])

BOSHgDat <- data.frame(
  x = BOS.dm.dat$LON,
  y = BOS.dm.dat$LAT,
  Hg = BOS.dm.dat$lHg,
  gs = BOS.dm.dat$lgs,
  clay = BOS.dm.dat$lgtclay,
  pH = BOS.dm.dat$pH15,
  Prescott = BOS.dm.dat$prescott,
  soilH20 = BOS.dm.dat$soilH20,
  lai = BOS.dm.dat$llai,
  soilN = BOS.dm.dat$soillN,
  soilP = BOS.dm.dat$soillP,
  KThU = BOS.dm.dat$lKThU,
  ferric2=BOS.dm.dat$ferric2,
  Alox=BOS.dm.dat$Alox,
  Feox=BOS.dm.dat$Feox,
  Pox=BOS.dm.dat$Pox
)

# convert to spatial points
BOSHgSpat <- BOSHgDat
coordinates(BOSHgSpat) <- ~x+y

# prepare training data
BOStraining_data <- data.frame(
  Hg = BOSHgSpat$Hg,
  # environmental predictors
  clay = BOSHgSpat$clay,
  pH = BOSHgSpat$pH,
  Prescott = BOSHgSpat$Prescott,
  soilH20 = BOSHgSpat$soilH20,
  lai = BOSHgSpat$lai,
  soilN = BOSHgSpat$soilN,
  soilP = BOSHgSpat$soilP,
  KThU = BOSHgSpat$KThU,
  ferric2=BOSHgSpat$ferric2,
  Alox=BOSHgSpat$Alox,
  Feox=BOSHgSpat$Feox,
  Pox=BOSHgSpat$Pox
)

crs(BOSHgSpat) <- crsUse

# choose size of spatial blocks for spatial autocorrelation
BOSblocksize <- cv_spatial_autocor(x=BOSHgSpat, column='Hg', plot=T, deg_to_metre=111139) 
BOSblocksize$range

# spatial blocks for cross-validation
BOSspatial_blocks <- cv_spatial(x=BOSHgSpat, k=5, hexagon=F, size=200000)
BOSspatial_blocks$folds_list

# custom indices for spatial cross-validation
BOSspatial_indices <- BOSspatial_blocks$folds_list

# define spatial cross-validation control
BOSsactrl <- trainControl(
  method = "cv",
  classProbs = F,
  number = 5,
  index = BOSspatial_indices,
  savePredictions = TRUE
)

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# train random forest model (no control for spatial autocorrelation)
BOSrf_model <- train(
  Hg ~ .,
  data = BOStraining_data,
  method = "rf",
  trControl = trainControl(method = "cv", number = 10)
)

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

BOSrf_model$results
BOSrf_model$finalModel$importance
str(BOSrf_model)

# spatial predictions
BOSHgPred <- terra::predict(envVars, BOSrf_model, na.rm=T)
plot(BOSHgPred)
points(BOSHgDat, pch=20, col="black", cex=0.5)
BOSHgPred.rf.bt <- TOSHgPred^10
plot(BOSHgPred.rf.bt)
writeRaster(BOSHgPred.rf.bt, "BOSHgPredRFbt.tif", overwrite=T)

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# model cross-validation
BOScv <- train(
  Hg ~ .,
  data = BOStraining_data,
  method = "rf",
  trControl = trainControl(
    method = "cv",
    number = 10,
    savePredictions = TRUE
  )
)

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# performance metrics
BOSrmse <- sqrt(mean((BOScv$pred$obs - BOScv$pred$pred)^2))
BOSrmse
BOSr2 <- cor(BOScv$pred$obs, BOScv$pred$pred)^2
BOSr2

# variable importance
varImp(BOSrf_model)

# plot variable importance
plot(varImp(BOSrf_model))



###############
## spatialRF ##
###############

## TOS
# create distance matrix
# Haversine distance matrix function
# most accurate for geographic coordinates because it accounts for Earth's spherical shape
TOScoords <- data.frame(
  longitude = TOSHgDat$x,
  latitude = TOSHgDat$y
)

TOShaversine_matrix <- distm(
  x = TOScoords,
  fun = distHaversine
)

# format the distance matrix for spatialRF
#' spatialRF expects: square matrix; row and column names; symmetric distances
# add row and column names
rownames(TOShaversine_matrix) <- paste0("point_", 1:dim(TOSHgDat)[1])
colnames(TOShaversine_matrix) <- paste0("point_", 1:dim(TOSHgDat)[1])

# verify the distance matrix properties
#  checks ensure  matrix is properly formatted for spatialRF
stopifnot(
  nrow(TOShaversine_matrix) == ncol(TOShaversine_matrix),  # square matrix
  isSymmetric(TOShaversine_matrix),                          # symmetric
  all(diag(TOShaversine_matrix) == 0)                        # zero diagonal
)


#names of the response variable and the predictors
dependent.variable.name <- "Hg"
TOSpredictor.variable.names <- colnames(TOSHgDat)[5:dim(TOSHgDat)[2]]

#coordinates of the cases
TOSxy <- TOSHgDat[, c("x", "y")]

#distance matrix
TOSdistance.matrix <- TOShaversine_matrix
range(TOSdistance.matrix)

#distance thresholds (same units as distance_matrix)
distance.thresholds <- seq(from=0, to=.9e6, by=5e4)

world <- rnaturalearth::ne_countries(
  scale = "medium", 
  returnclass = "sf"
)

# plot data
ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf, 
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = TOSHgDat,
    ggplot2::aes(
      x = x,
      y = y,
      color = Hg
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(
    direction = -1, 
    option = "F"
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "[Hg]") +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::ggtitle("soil [Hg]") +
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

# scatterplots of the response variable (y axis) against each predictor (x axis)
spatialRF::plot_training_df(
  data = TOSHgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = TOSpredictor.variable.names,
  ncol = 3,
  point.color = rev(viridis::viridis(100, option = "F")),
  line.color = "gray30"
)

# Moran's I
TOSmoransIplot <- spatialRF::plot_training_df_moran(
  data = TOSHgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = TOSpredictor.variable.names,
  distance.matrix = TOSdistance.matrix,
  distance.thresholds = distance.thresholds,
  fill.color = viridis::viridis(
    100,
    option = "F",
    direction = -1
  ),
  point.color = "gray40"
)
TOSmoransIplot
str(TOSmoransIplot)

TOSmoransIplot$data$distance.threshold
TOSmoransIplot$data$p.value
write.csv(TOSmoransIplot$data, "TOSmoransIplot.csv", row.names=F)


## BOS
# create distance matrix
# Haversine distance matrix function
# most accurate for geographic coordinates because it accounts for Earth's spherical shape
BOScoords <- data.frame(
  longitude = BOSHgDat$x,
  latitude = BOSHgDat$y
)

BOShaversine_matrix <- distm(
  x = BOScoords,
  fun = distHaversine
)

# format the distance matrix for spatialRF
#' spatialRF expects: square matrix; row and column names; symmetric distances
# add row and column names
rownames(BOShaversine_matrix) <- paste0("point_", 1:dim(BOSHgDat)[1])
colnames(BOShaversine_matrix) <- paste0("point_", 1:dim(BOSHgDat)[1])

# verify the distance matrix properties
#  checks ensure  matrix is properly formatted for spatialRF
stopifnot(
  nrow(BOShaversine_matrix) == ncol(BOShaversine_matrix),  # square matrix
  isSymmetric(BOShaversine_matrix),                          # symmetric
  all(diag(BOShaversine_matrix) == 0)                        # zero diagonal
)


#names of the response variable and the predictors
dependent.variable.name <- "Hg"
BOSpredictor.variable.names <- colnames(BOSHgDat)[5:dim(BOSHgDat)[2]]

#coordinates of the cases
BOSxy <- BOSHgDat[, c("x", "y")]

#distance matrix
BOSdistance.matrix <- BOShaversine_matrix
range(BOSdistance.matrix)

# plot data
ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf, 
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = BOSHgDat,
    ggplot2::aes(
      x = x,
      y = y,
      color = Hg
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(
    direction = -1, 
    option = "F"
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "[Hg]") +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::ggtitle("soil [Hg]") +
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

# scatterplots of the response variable (y axis) against each predictor (x axis)
spatialRF::plot_training_df(
  data = BOSHgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = BOSpredictor.variable.names,
  ncol = 3,
  point.color = rev(viridis::viridis(100, option = "F")),
  line.color = "gray30"
)

# Moran's I
BOSmoransIplot <- spatialRF::plot_training_df_moran(
  data = BOSHgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = BOSpredictor.variable.names,
  distance.matrix = BOSdistance.matrix,
  distance.thresholds = distance.thresholds,
  fill.color = viridis::viridis(
    100,
    option = "F",
    direction = -1
  ),
  point.color = "gray40"
)
BOSmoransIplot
str(BOSmoransIplot)

BOSmoransIplot$data$distance.threshold
BOSmoransIplot$data$p.value
write.csv(BOSmoransIplot$data, "BOSmoransIplot.csv", row.names=F)



## TOS
# possible interactions
TOSinteractions <- spatialRF::the_feature_engineer(
  data = TOSHgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = TOSpredictor.variable.names,
  xy = TOSxy,
  importance.threshold = 0.50, # uses 50% best predictors
  cor.threshold = 0.60, # max corr between interactions and predictors
  repetitions = 100,
  verbose = TRUE
)

# **run only if promising interactions found**
kableExtra::kbl(
  head(TOSinteractions$screening, 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# non-spatial random forest
TOSmodel.non.spatial <- spatialRF::rf(
  data = TOSHgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = TOSpredictor.variable.names,
  distance.matrix = TOSdistance.matrix,
  distance.thresholds = distance.thresholds,
  xy = TOSxy, # not needed by rf, but other functions read it from the model
  verbose = FALSE,
  scaled.importance=T
)

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# residual diagnostics
spatialRF::plot_residuals_diagnostics(
  TOSmodel.non.spatial,
  verbose = FALSE
)

# output multi-scale Moran's I of non-spatial model residuals
TOSns.resids <- spatialRF::plot_residuals_diagnostics(
  TOSmodel.non.spatial,
  verbose = FALSE
)
str(TOSns.resids)
TOSns.MoranI.resid.dat <- TOSns.resids[[3]]$data
write.csv(TOSns.MoranI.resid.dat, "TOSnsMoransIresid.csv", row.names=F)

# variable importance
spatialRF::plot_importance(
  TOSmodel.non.spatial,
  verbose = FALSE
)

TOSnsm.imp <- spatialRF::plot_importance(
  TOSmodel.non.spatial,
  verbose = FALSE
)
TOSnsm.imp

TOSimportance.df <- randomForestExplainer::measure_importance(
  TOSmodel.non.spatial,
  measures = c("mean_min_depth", "no_of_nodes", "times_a_root", "p_value")
)

kableExtra::kbl(
  TOSimportance.df %>% 
    dplyr::arrange(mean_min_depth) %>% 
    dplyr::mutate(p_value = round(p_value, 4)),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# contribution of predictors to model transferability
# (predictive ability on independent spatial folds measured with rf_evaluate())
TOSmodel.non.spatial <- spatialRF::rf_importance(
  model = TOSmodel.non.spatial
)

TOSmodel.non.spatial$importance$per.variable %>% 
  ggplot2::ggplot() +
  ggplot2::aes(
    x = importance.oob,
    y = importance.cv
  ) + 
  ggplot2::geom_point(size = 3) + 
  ggplot2::theme_bw() +
  ggplot2::xlab("importance (out-of-bag)") + 
  ggplot2::ylab("contribution to transferability") + 
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "red4")

# local importance
TOSlocal.importance <- spatialRF::get_importance_local(TOSmodel.non.spatial)
kableExtra::kbl(
  round(TOSlocal.importance[1:10, 1:5], 3),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# map changing variable importance of space
# adding coordinates
TOSlocal.importance <- cbind(
  TOSxy,
  TOSlocal.importance
)
# colours
color.low <- viridis::viridis(
  3,
  option = "F"
)[2]
color.high <- viridis::viridis(
  3,
  option = "F"
)[1]

p1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = TOSlocal.importance,
    ggplot2::aes(
      x = x,
      y = y,
      color = ferric2
    )
  ) +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::scale_color_gradient2(
    low = color.low, 
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") + 
  ggplot2::ggtitle("ferric PC2") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) + 
  ggplot2::labs(color = "importance") + 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

p2 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = TOSlocal.importance,
    ggplot2::aes(
      x = x,
      y = y,
      color = soilN
    )
  ) +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
ggplot2::scale_color_gradient2(
    low = color.low, 
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") + 
  ggplot2::ggtitle("soil N") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) + 
  ggplot2::labs(color = "importance") + 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

p3 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = TOSlocal.importance,
    ggplot2::aes(
      x = x,
      y = y,
      color = lai
    )
  ) +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::scale_color_gradient2(
    low = color.low, 
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::ggtitle("LAI") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) + 
  ggplot2::labs(color = "importance") + 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

p4 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = TOSlocal.importance,
    ggplot2::aes(
      x = x,
      y = y,
      color = Pox
    )
  ) +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::scale_color_gradient2(
    low = color.low, 
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::ggtitle("phosphorus oxide") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) + 
  ggplot2::labs(color = "importance") + 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

p5 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = TOSlocal.importance,
    ggplot2::aes(
      x = x,
      y = y,
      color = Prescott
    )
  ) +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::scale_color_gradient2(
    low = color.low, 
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::ggtitle("Prescott index") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) + 
  ggplot2::labs(color = "importance") + 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

p6 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = TOSlocal.importance,
    ggplot2::aes(
      x = x,
      y = y,
      color = soilH20
    )
  ) +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::scale_color_gradient2(
    low = color.low, 
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::ggtitle("soil H2O availability") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) + 
  ggplot2::labs(color = "importance") + 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")


(p1 + p2 + p3) / (p4 + p5 + p6)


## BOS
# possible interactions
BOSinteractions <- spatialRF::the_feature_engineer(
  data = BOSHgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = BOSpredictor.variable.names,
  xy = BOSxy,
  importance.threshold = 0.50, # uses 50% best predictors
  cor.threshold = 0.60, # max corr between interactions and predictors
  repetitions = 100,
  verbose = TRUE
)

# **run only if promising interactions found**
kableExtra::kbl(
  head(BOSinteractions$screening, 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# non-spatial random forest
BOSmodel.non.spatial <- spatialRF::rf(
  data = BOSHgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = BOSpredictor.variable.names,
  distance.matrix = BOSdistance.matrix,
  distance.thresholds = distance.thresholds,
  xy = BOSxy, # not needed by rf, but other functions read it from the model
  verbose = FALSE,
  scaled.importance=T
)

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# residual diagnostics
spatialRF::plot_residuals_diagnostics(
  BOSmodel.non.spatial,
  verbose = FALSE
)

# output multi-scale Moran's I of non-spatial model residuals
BOSns.resids <- spatialRF::plot_residuals_diagnostics(
  BOSmodel.non.spatial,
  verbose = FALSE
)
str(BOSns.resids)
BOSns.MoranI.resid.dat <- BOSns.resids[[3]]$data
write.csv(BOSns.MoranI.resid.dat, "BOSnsMoransIresid.csv", row.names=F)

# variable importance
spatialRF::plot_importance(
  BOSmodel.non.spatial,
  verbose = FALSE
)

BOSnsm.imp <- spatialRF::plot_importance(
  BOSmodel.non.spatial,
  verbose = FALSE
)
BOSnsm.imp$data

BOSimportance.df <- randomForestExplainer::measure_importance(
  BOSmodel.non.spatial,
  measures = c("mean_min_depth", "no_of_nodes", "times_a_root", "p_value")
)

kableExtra::kbl(
  BOSimportance.df %>% 
    dplyr::arrange(mean_min_depth) %>% 
    dplyr::mutate(p_value = round(p_value, 4)),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# contribution of predictors to model transferability
# (predictive ability on independent spatial folds measured with rf_evaluate())
BOSmodel.non.spatial <- spatialRF::rf_importance(
  model = BOSmodel.non.spatial
)

BOSmodel.non.spatial$importance$per.variable %>% 
  ggplot2::ggplot() +
  ggplot2::aes(
    x = importance.oob,
    y = importance.cv
  ) + 
  ggplot2::geom_point(size = 3) + 
  ggplot2::theme_bw() +
  ggplot2::xlab("importance (out-of-bag)") + 
  ggplot2::ylab("contribution to transferability") + 
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "red4")

# local importance
BOSlocal.importance <- spatialRF::get_importance_local(BOSmodel.non.spatial)
kableExtra::kbl(
  round(BOSlocal.importance[1:10, 1:5], 3),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# map changing variable importance of space
# adding coordinates
BOSlocal.importance <- cbind(
  BOSxy,
  BOSlocal.importance
)
# colours
color.low <- viridis::viridis(
  3,
  option = "F"
)[2]
color.high <- viridis::viridis(
  3,
  option = "F"
)[1]

p1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = BOSlocal.importance,
    ggplot2::aes(
      x = x,
      y = y,
      color = ferric2
    )
  ) +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::scale_color_gradient2(
    low = color.low, 
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") + 
  ggplot2::ggtitle("ferric PC2") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) + 
  ggplot2::labs(color = "importance") + 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

p2 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = BOSlocal.importance,
    ggplot2::aes(
      x = x,
      y = y,
      color = soilN
    )
  ) +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::scale_color_gradient2(
    low = color.low, 
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") + 
  ggplot2::ggtitle("soil N") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) + 
  ggplot2::labs(color = "importance") + 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

p3 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = BOSlocal.importance,
    ggplot2::aes(
      x = x,
      y = y,
      color = lai
    )
  ) +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::scale_color_gradient2(
    low = color.low, 
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::ggtitle("LAI") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) + 
  ggplot2::labs(color = "importance") + 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

p4 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = BOSlocal.importance,
    ggplot2::aes(
      x = x,
      y = y,
      color = Pox
    )
  ) +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::scale_color_gradient2(
    low = color.low, 
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::ggtitle("phosphorus oxide") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) + 
  ggplot2::labs(color = "importance") + 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

p5 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = BOSlocal.importance,
    ggplot2::aes(
      x = x,
      y = y,
      color = Prescott
    )
  ) +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::scale_color_gradient2(
    low = color.low, 
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::ggtitle("Prescott index") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) + 
  ggplot2::labs(color = "importance") + 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

p6 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = aus.sf,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = BOSlocal.importance,
    ggplot2::aes(
      x = x,
      y = y,
      color = soilH20
    )
  ) +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::scale_color_gradient2(
    low = color.low, 
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::ggtitle("soil H2O availability") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) + 
  ggplot2::labs(color = "importance") + 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")


(p1 + p2 + p3) / (p4 + p5 + p6)



## TOS
# response curves and surfaces
spatialRF::plot_response_curves(
  TOSmodel.non.spatial,
  quantiles = c(0.1, 0.5, 0.9),
  line.color = viridis::viridis(
    3, # same number of colours as quantiles
    option = "F", 
    end = 0.9
  ),
  ncol = 4,
  show.data = TRUE
)
spatialRF::plot_response_curves(
  TOSmodel.non.spatial,
  quantiles = 0.5,
  ncol = 4
)

## BOS
# response curves and surfaces
spatialRF::plot_response_curves(
  BOSmodel.non.spatial,
  quantiles = c(0.1, 0.5, 0.9),
  line.color = viridis::viridis(
    3, # same number of colours as quantiles
    option = "F", 
    end = 0.9
  ),
  ncol = 4,
  show.data = TRUE
)
spatialRF::plot_response_curves(
  BOSmodel.non.spatial,
  quantiles = 0.5,
  ncol = 4
)


# model performance
# R squared (oob) and RMSE (oob) are the R squared of the model and its root mean squared error when predicting
# the out-of-bag data (fraction of data not used to train individual trees).
# From all the values available in the performance slot, probably these are the most honest ones, because
# it is the closest trying to get a performance estimate on independent data.
# However, out-of-bag data are not fully independent, and therefore will still be inflated,
# especially if the data is highly aggregated in space.
# R squared and pseudo R squared are computed from the observations and the predictions, and indicate to what
# extent model outcomes represent the input data. These values will usually be high if the data are
# highly aggregated in space.
# RMSE and its normalized version are computed via root_mean_squared_error(), and are linear with R squared
# and pseudo R squared.

spatialRF::print_performance(TOSmodel.non.spatial)
spatialRF::print_performance(BOSmodel.non.spatial)


## TOS
# spatial cross-validation
TOSmodel.non.spatial <- spatialRF::rf_evaluate(
  model = TOSmodel.non.spatial,
  xy = TOSxy,                  # data coordinates
  repetitions = 30,         # number of spatial folds
  training.fraction = 0.75, # training data fraction on each fold
  metrics = "r.squared",
  verbose = FALSE
)
spatialRF::plot_evaluation(TOSmodel.non.spatial)

# predicted values of [Hg]
TOSpredicted.non.spatial <- stats::predict(object = TOSmodel.non.spatial, data = TOSHgDat, type = "response")$predictions

# spatial random forest
# Moran's I
spatialRF::plot_moran(
  TOSmodel.non.spatial, 
  verbose = FALSE
)

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  #1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

TOSmodel.spatial <- spatialRF::rf_spatial(
   model = TOSmodel.non.spatial,
   method = "mem.moran.sequential", # default method
   verbose = FALSE
 )

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

spatialRF::print_performance(TOSmodel.spatial)

spatialRF::plot_moran(
  TOSmodel.spatial, 
  verbose = FALSE
)

TOSmodel.spatial$importance$per.variable$variable


## BOS
# spatial cross-validation
BOSmodel.non.spatial <- spatialRF::rf_evaluate(
  model = BOSmodel.non.spatial,
  xy = BOSxy,                  # data coordinates
  repetitions = 30,         # number of spatial folds
  training.fraction = 0.75, # training data fraction on each fold
  metrics = "r.squared",
  verbose = FALSE
)
spatialRF::plot_evaluation(BOSmodel.non.spatial)

# predicted values of [Hg]
BOSpredicted.non.spatial <- stats::predict(object = BOSmodel.non.spatial, data = BOSHgDat, type = "response")$predictions

# spatial random forest
# Moran's I
spatialRF::plot_moran(
  BOSmodel.non.spatial, 
  verbose = FALSE
)

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  #1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

BOSmodel.spatial <- spatialRF::rf_spatial(
  model = BOSmodel.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE
)

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

spatialRF::print_performance(BOSmodel.spatial)

spatialRF::plot_moran(
  BOSmodel.spatial, 
  verbose = FALSE
)

BOSmodel.spatial$importance$per.variable$variable


## redo environmental variables stack
names(envVars)
envVars2 <- c(soilN.rsmp,lai.rsmp,Prescott.rsmp,soilH20.rsmp,KThU.rsmp,soilP.rsmp,pH.rsmp,clay.rsmp,
              ferric2.rsmp,Alox.rsmp,Feox.rsmp,Pox.rsmp)
names(envVars2) <- c("soilN", "lai", "Prescott", "soilH20", "KThU", "soilP", "pH", "clay", "ferric2",
                     "Alox","Feox","Pox")
names(envVars2)
plot(envVars2)
envVars2


## TOS
TOSHgPredNoSpat <- terra::predict(envVars2, TOSmodel.non.spatial, na.rm=T)
plot(TOSHgPredNoSpat)
points(TOSHgDat, pch=20, col="black", cex=0.5)
TOSHgPredNoSpat.rf.bt <- TOSHgPredNoSpat^10
plot(TOSHgPredNoSpat.rf.bt)
writeRaster(TOSHgPredNoSpat, "TOSHgPredNoSpat.tif", overwrite=T)
writeRaster(TOSHgPredNoSpat.rf.bt, "TOSHgPredNoSpatRFbt.tif", overwrite=T)

TOSHgPredNoSpatImp <- rast('TOSHgPredNoSpat.tif')
writeRaster(TOSHgPredNoSpatImp, "TOSHgPredNoSpatImplog10.nc", filetype="netcdf", overwrite=TRUE) # write to NetCDF format
TOSHgPredNoSpatRFbtImp <- rast('TOSHgPredNoSpatRFbt.tif')
writeRaster(TOSHgPredNoSpatRFbtImp, "TOSHgPredNoSpatRFbtImp.nc", filetype="netcdf", overwrite=TRUE) # write to NetCDF format

p1 <- spatialRF::plot_importance(
  TOSmodel.non.spatial, 
  verbose = FALSE) + 
  ggplot2::ggtitle("TOS non-spatial model") 

p2 <- spatialRF::plot_importance(
  TOSmodel.spatial,
  verbose = FALSE) + 
  ggplot2::ggtitle("TOS spatial model")

p1 | p2 

kableExtra::kbl(
  head(TOSmodel.spatial$importance$per.variable, n = 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)


## BOS
BOSHgPredNoSpat <- terra::predict(envVars2, BOSmodel.non.spatial, na.rm=T)
plot(BOSHgPredNoSpat)
points(BOSHgDat, pch=20, col="black", cex=0.5)
BOSHgPredNoSpat.rf.bt <- BOSHgPredNoSpat^10
plot(BOSHgPredNoSpat.rf.bt)
writeRaster(BOSHgPredNoSpat, "BOSHgPredNoSpat.tif", overwrite=T)
writeRaster(BOSHgPredNoSpat.rf.bt, "BOSHgPredNoSpatRFbt.tif", overwrite=T)

BOSHgPredNoSpatImp <- rast('BOSHgPredNoSpat.tif')
writeRaster(BOSHgPredNoSpatImp, "BOSHgPredNoSpatImplog10.nc", filetype="netcdf", overwrite=TRUE) # write to NetCDF format
BOSHgPredNoSpatRFbtImp <- rast('BOSHgPredNoSpatRFbt.tif')
writeRaster(BOSHgPredNoSpatRFbtImp, "BOSHgPredNoSpatRFbtImp.nc", filetype="netcdf", overwrite=TRUE) # write to NetCDF format

p1 <- spatialRF::plot_importance(
  BOSmodel.non.spatial, 
  verbose = FALSE) + 
  ggplot2::ggtitle("BOS non-spatial model") 

p2 <- spatialRF::plot_importance(
  BOSmodel.spatial,
  verbose = FALSE) + 
  ggplot2::ggtitle("BOS spatial model")

p1 | p2 

kableExtra::kbl(
  head(BOSmodel.spatial$importance$per.variable, n = 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)


## TOS
# model tuning
TOSmodel.spatial <- spatialRF::rf_spatial(
  model = TOSmodel.spatial,
  method = "mem.moran.sequential", # default method
  ranger.arguments = list(
    mtry = 5,
    min.node.size = 20,
    num.trees = 500
  ),
  verbose = FALSE,
  seed = random.seed
)

TOSmodel.spatial <- rf_tuning(
  model = TOSmodel.spatial,
  xy = TOSxy,
  repetitions = 30,
  num.trees = c(500, 1000),
  mtry = seq(
    2,
    length(TOSmodel.spatial$ranger.arguments$predictor.variable.names), # number of predictors
    by = 9),
  min.node.size = c(5, 15),
  verbose = FALSE
)

## BOS
# model tuning
BOSmodel.spatial <- spatialRF::rf_spatial(
  model = BOSmodel.spatial,
  method = "mem.moran.sequential", # default method
  ranger.arguments = list(
    mtry = 5,
    min.node.size = 20,
    num.trees = 500
  ),
  verbose = FALSE,
  seed = random.seed
)

BOSmodel.spatial <- rf_tuning(
  model = BOSmodel.spatial,
  xy = BOSxy,
  repetitions = 30,
  num.trees = c(500, 1000),
  mtry = seq(
    2,
    length(BOSmodel.spatial$ranger.arguments$predictor.variable.names), # number of predictors
    by = 9),
  min.node.size = c(5, 15),
  verbose = FALSE
)


# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  #1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)


## TOS
# stochastic version
TOSmodel.spatial.repeat <- spatialRF::rf_repeat(
  model = TOSmodel.spatial, 
  repetitions = 250,
  verbose = FALSE
)

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# variable importance
TOSstochSpatialRFvarImp <- spatialRF::plot_importance(
  TOSmodel.spatial.repeat, 
  verbose = FALSE
)
TOSstochSpatialRFvarImp
str(TOSstochSpatialRFvarImp)
TOSstochSpatialRFvarImp.dat <- TOSstochSpatialRFvarImp$data

TOSviXvar.stats <- TOSstochSpatialRFvarImp.dat %>%
  group_by(variable) %>%
  summarise(
    mean = mean(importance, na.rm = TRUE),
    median = median(importance, na.rm = TRUE), 
    var = var(importance, na.rm = TRUE),
    sd = sd(importance, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(importance, probs=0.975, na.rm = TRUE),
    lower = quantile(importance, probs=0.025, na.rm = TRUE),
    n = n()
  )
TOSviXvar.stats
write.csv(TOSviXvar.stats, "TOSviXvarstats.csv", row.names=F)

# response curves
spatialRF::plot_response_curves(model=TOSmodel.spatial.repeat, quantiles = c(0.5, 0.975, 0.025), ncol = 3)


## BOS
# stochastic version
BOSmodel.spatial.repeat <- spatialRF::rf_repeat(
  model = BOSmodel.spatial, 
  repetitions = 250,
  verbose = FALSE
)

# variable importance
BOSstochSpatialRFvarImp <- spatialRF::plot_importance(
  BOSmodel.spatial.repeat, 
  verbose = FALSE
)
BOSstochSpatialRFvarImp
str(BOSstochSpatialRFvarImp)
BOSstochSpatialRFvarImp.dat <- BOSstochSpatialRFvarImp$data

BOSviXvar.stats <- BOSstochSpatialRFvarImp.dat %>%
  group_by(variable) %>%
  summarise(
    mean = mean(importance, na.rm = TRUE),
    median = median(importance, na.rm = TRUE), 
    var = var(importance, na.rm = TRUE),
    sd = sd(importance, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(importance, probs=0.975, na.rm = TRUE),
    lower = quantile(importance, probs=0.025, na.rm = TRUE),
    n = n()
  )
BOSviXvar.stats
write.csv(BOSviXvar.stats, "BOSviXvarstats.csv", row.names=F)

# response curves
spatialRF::plot_response_curves(model=BOSmodel.spatial.repeat, quantiles = c(0.5, 0.975, 0.025), ncol = 3)





#########################
## tuned spatial model ##
## stochastic          ##
#########################

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# redefine distance thresholds vector (shortened for reduced processing time)
distance.thresholds2 <- seq(from=1e5, to=5e5, by=5e4)

## TOS
  TOSmodel.tuned <- spatialRF::rf_spatial(
  data = TOSHgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = TOSpredictor.variable.names,
  distance.matrix = TOSdistance.matrix,
  distance.thresholds = distance.thresholds2,
  xy = TOSxy,
  cluster = local.cluster,
  scaled.importance=T,
  method = c("mem.moran.sequential"),
  max.spatial.predictors = 40, # main spatial predictors only
  ranger.arguments=list(
    num.trees=500,
    mtry=1,
    min.node.size=40)
) %>%
  rf_evaluate() %>%
  rf_repeat(repetitions = 250)

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# optimisation plot
TOSp <- spatialRF::plot_optimization(TOSmodel.tuned)

# importance
spatialRF::plot_importance(
  TOSmodel.tuned, 
  verbose = FALSE
)

# variable importance
TOSstochSpatTunedRFvarImp <- spatialRF::plot_importance(
  TOSmodel.tuned, 
  verbose = FALSE
)
TOSstochSpatTunedRFvarImp
str(TOSstochSpatTunedRFvarImp)
TOSstochSpatTunedRFvarImp.dat <- TOSstochSpatTunedRFvarImp$data

TOSviXvar.stats <- TOSstochSpatTunedRFvarImp.dat %>%
  group_by(variable) %>%
  summarise(
    mean = mean(importance, na.rm = TRUE),
    median = median(importance, na.rm = TRUE), 
    var = var(importance, na.rm = TRUE),
    sd = sd(importance, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(importance, probs=0.975, na.rm = TRUE),
    lower = quantile(importance, probs=0.025, na.rm = TRUE),
    n = n()
  )
TOSviXvar.stats
write.csv(TOSviXvar.stats, "TOSviXvarstatsTunedSpatial.csv", row.names=F)

# performance
spatialRF::print_performance(TOSmodel.tuned)
TOSmodel.tuned$importance$per.variable$variable

# spatial predictors
TOSspatial.predictors <- spatialRF::get_spatial_predictors(TOSmodel.tuned)
class(TOSspatial.predictors)
TOSspatpred.names <- colnames(TOSspatial.predictors)
TOSNspatpreds <- length(TOSspatpred.names)
TOSNspatpreds

TOSpr <- data.frame(TOSspatial.predictors, TOSHgDat[, c("x", "y")])
head(TOSpr)

kableExtra::kbl(
  head(TOSmodel.tuned$importance$per.variable, n = 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

TOSp1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = aus.sf, fill = "white") +
  ggplot2::geom_point(
    data = TOSpr,
    ggplot2::aes(
      x = x,
      y = y,
      color = spatial_predictor_100000_28	
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "eigenvalue") +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::ggtitle("spatial_predictor_100000_28") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

TOSp2 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = aus.sf, fill = "white") +
  ggplot2::geom_point(
    data = TOSpr,
    ggplot2::aes(
      x = x,
      y = y,
      color = spatial_predictor_100000_30
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "eigenvalue") +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::ggtitle("spatial_predictor_100000_30") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

TOSp3 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = aus.sf, fill = "white") +
  ggplot2::geom_point(
    data = TOSpr,
    ggplot2::aes(
      x = x,
      y = y,
      color = spatial_predictor_100000_24	
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "eigenvalue") +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::ggtitle("spatial_predictor_100000_24") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

TOSp4 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = aus.sf, fill = "white") +
  ggplot2::geom_point(
    data = TOSpr,
    ggplot2::aes(
      x = x,
      y = y,
      color = spatial_predictor_100000_37
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "eigenvalue") +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::ggtitle("spatial_predictor_100000_37") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

(TOSp1 + TOSp2) / (TOSp3 + TOSp4)


# response curves
TOSrespCurvRF <- spatialRF::plot_response_curves(model=TOSmodel.tuned, quantiles = c(0.5), ncol = 3)
TOSrespCurvRF
str(TOSrespCurvRF[[1]])

TOSrespCurvRF[[1]]$labels$x


# ferric2
TOSrespCurvRF[[1]]
TOSrespCurvRF[[1]]$labels$x
TOSferricrespCurvRF <- TOSrespCurvRF[[1]]
str(TOSferricrespCurvRF)

table(TOSferricrespCurvRF$layers[[1]]$data$id)

TOSferricpredRF.dat <- data.frame(iter=TOSferricrespCurvRF$layers[[1]]$data$id, 
                               ferric=TOSferricrespCurvRF$layers[[1]]$data$ferric,
                               Hg=TOSferricrespCurvRF$layers[[1]]$data$Hg,
                               quantile=TOSferricrespCurvRF$layers[[1]]$data$quantile)

head(TOSferricpredRF.dat)
TOSferricpredRFq.5.dat <- subset(TOSferricpredRF.dat, quantile==0.5)
dim(TOSferricpredRFq.5.dat)

TOSferricpredRF.5Xiter.stats <- TOSferricpredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(ferric, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
TOSferricpredRF.5Xiter.stats
range(TOSferricpredRF.5Xiter.stats$x, na.rm=T)
c(min(TOSferricpredRF.5Xiter.stats$lowery, na.rm=T), max(TOSferricpredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(TOSferricpredRF.5Xiter.stats, "TOSferricpredRF.5Xiter.stats.csv", row.names=F)


# soil N
TOSrespCurvRF[[2]]$labels$x
TOSsoilNrespCurvRF <- TOSrespCurvRF[[2]]
str(TOSsoilNrespCurvRF)

table(TOSsoilNrespCurvRF$layers[[1]]$data$id)

TOSsoilNpredRF.dat <- data.frame(iter=TOSsoilNrespCurvRF$layers[[1]]$data$id, 
           soilN=TOSsoilNrespCurvRF$layers[[1]]$data$soilN,
           Hg=TOSsoilNrespCurvRF$layers[[1]]$data$Hg,
           quantile=TOSsoilNrespCurvRF$layers[[1]]$data$quantile)
head(TOSsoilNpredRF.dat)
TOSsoilNpredRFq.5.dat <- subset(TOSsoilNpredRF.dat, quantile==0.5)
dim(TOSsoilNpredRFq.5.dat)

TOSsoilNpredRF.5Xiter.stats <- TOSsoilNpredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(soilN, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
TOSsoilNpredRF.5Xiter.stats
range(TOSsoilNpredRF.5Xiter.stats$x, na.rm=T)
c(min(TOSsoilNpredRF.5Xiter.stats$lowery, na.rm=T), max(TOSsoilNpredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(TOSsoilNpredRF.5Xiter.stats, "TOSsoilNpredRF.5Xiter.stats.csv", row.names=F)

# LAI
TOSrespCurvRF[[3]]$labels$x
TOSlairespCurvRF <- TOSrespCurvRF[[3]]
str(TOSlairespCurvRF)

table(TOSlairespCurvRF$layers[[1]]$data$id)

TOSlaipredRF.dat <- data.frame(iter=TOSlairespCurvRF$layers[[1]]$data$id, 
                            lai=TOSlairespCurvRF$layers[[1]]$data$lai,
                            Hg=TOSlairespCurvRF$layers[[1]]$data$Hg,
                            quantile=TOSlairespCurvRF$layers[[1]]$data$quantile)
head(TOSlaipredRF.dat)
TOSlaipredRFq.5.dat <- subset(TOSlaipredRF.dat, quantile==0.5)
dim(TOSlaipredRFq.5.dat)

TOSlaipredRF.5Xiter.stats <- TOSlaipredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(lai, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
TOSlaipredRF.5Xiter.stats
range(TOSlaipredRF.5Xiter.stats$x, na.rm=T)
c(min(TOSlaipredRF.5Xiter.stats$lowery, na.rm=T), max(TOSlaipredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(TOSlaipredRF.5Xiter.stats, "TOSlaipredRF.5Xiter.stats.csv", row.names=F)

# Pox
TOSrespCurvRF[[4]]$labels$x
TOSPoxrespCurvRF <- TOSrespCurvRF[[4]]
str(TOSPoxrespCurvRF)

table(TOSPoxrespCurvRF$layers[[1]]$data$id)

TOSPoxpredRF.dat <- data.frame(iter=TOSPoxrespCurvRF$layers[[1]]$data$id,
                            Pox=TOSPoxrespCurvRF$layers[[1]]$data$Pox,
                            Hg=TOSPoxrespCurvRF$layers[[1]]$data$Hg,
                            quantile=TOSPoxrespCurvRF$layers[[1]]$data$quantile)
head(TOSPoxpredRF.dat)
TOSPoxpredRFq.5.dat <- subset(TOSPoxpredRF.dat, quantile==0.5)
dim(TOSPoxpredRFq.5.dat)

TOSPoxpredRF.5Xiter.stats <- TOSPoxpredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(Pox, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
TOSPoxpredRF.5Xiter.stats
range(TOSPoxpredRF.5Xiter.stats$x, na.rm=T)
c(min(TOSPoxpredRF.5Xiter.stats$lowery, na.rm=T), max(TOSPoxpredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(TOSPoxpredRF.5Xiter.stats, "TOSPoxpredRF.5Xiter.stats.csv", row.names=F)

# Prescott
TOSrespCurvRF[[5]]$labels$x
TOSPrescottrespCurvRF <- TOSrespCurvRF[[5]]
str(TOSPrescottrespCurvRF)

table(TOSPrescottrespCurvRF$layers[[1]]$data$id)

TOSPrescottpredRF.dat <- data.frame(iter=TOSPrescottrespCurvRF$layers[[1]]$data$id, 
                              Prescott=TOSPrescottrespCurvRF$layers[[1]]$data$Prescott,
                              Hg=TOSPrescottrespCurvRF$layers[[1]]$data$Hg,
                              quantile=TOSPrescottrespCurvRF$layers[[1]]$data$quantile)
head(TOSPrescottpredRF.dat)
TOSPrescottpredRFq.5.dat <- subset(TOSPrescottpredRF.dat, quantile==0.5)
dim(TOSPrescottpredRFq.5.dat)

TOSPrescottpredRF.5Xiter.stats <- TOSPrescottpredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(Prescott, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
TOSPrescottpredRF.5Xiter.stats
range(TOSPrescottpredRF.5Xiter.stats$x, na.rm=T)
c(min(TOSPrescottpredRF.5Xiter.stats$lowery, na.rm=T), max(TOSPrescottpredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(TOSPrescottpredRF.5Xiter.stats, "TOSPrescottpredRF.5Xiter.stats.csv", row.names=F)


# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

## BOS
BOSmodel.tuned <- spatialRF::rf_spatial(
  data = BOSHgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = BOSpredictor.variable.names,
  distance.matrix = BOSdistance.matrix,
  distance.thresholds = distance.thresholds2,
  xy = BOSxy,
  cluster = local.cluster,
  scaled.importance=T,
  method = c("mem.moran.sequential"),
  max.spatial.predictors = 40, # main spatial predictors only
  ranger.arguments=list(
    num.trees=500,
    mtry=1,
    min.node.size=40)
) %>%
  rf_evaluate() %>%
  rf_repeat(repetitions = 250)


# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# optimisation plot
BOSp <- spatialRF::plot_optimization(BOSmodel.tuned)

# importance
spatialRF::plot_importance(
  BOSmodel.tuned, 
  verbose = FALSE
)

# variable importance
BOSstochSpatTunedRFvarImp <- spatialRF::plot_importance(
  BOSmodel.tuned, 
  verbose = FALSE
)
BOSstochSpatTunedRFvarImp
str(BOSstochSpatTunedRFvarImp)
BOSstochSpatTunedRFvarImp.dat <- BOSstochSpatTunedRFvarImp$data

BOSviXvar.stats <- BOSstochSpatTunedRFvarImp.dat %>%
  group_by(variable) %>%
  summarise(
    mean = mean(importance, na.rm = TRUE),
    median = median(importance, na.rm = TRUE), 
    var = var(importance, na.rm = TRUE),
    sd = sd(importance, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(importance, probs=0.975, na.rm = TRUE),
    lower = quantile(importance, probs=0.025, na.rm = TRUE),
    n = n()
  )
BOSviXvar.stats
write.csv(BOSviXvar.stats, "BOSviXvarstatsTunedSpatial.csv", row.names=F)

# performance
spatialRF::print_performance(BOSmodel.tuned)
BOSmodel.tuned$importance$per.variable$variable

# spatial predictors
BOSspatial.predictors <- spatialRF::get_spatial_predictors(BOSmodel.tuned)
class(BOSspatial.predictors)
BOSspatpred.names <- colnames(BOSspatial.predictors)
BOSNspatpreds <- length(BOSspatpred.names)
BOSNspatpreds

BOSpr <- data.frame(BOSspatial.predictors, BOSHgDat[, c("x", "y")])
head(BOSpr)

kableExtra::kbl(
  head(BOSmodel.tuned$importance$per.variable, n = 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

BOSp1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = aus.sf, fill = "white") +
  ggplot2::geom_point(
    data = BOSpr,
    ggplot2::aes(
      x = x,
      y = y,
      color = spatial_predictor_100000_16	
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "eigenvalue") +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::ggtitle("spatial_predictor_100000_16") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

BOSp2 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = aus.sf, fill = "white") +
  ggplot2::geom_point(
    data = BOSpr,
    ggplot2::aes(
      x = x,
      y = y,
      color = spatial_predictor_100000_18
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "eigenvalue") +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::ggtitle("spatial_predictor_100000_18") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

BOSp3 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = aus.sf, fill = "white") +
  ggplot2::geom_point(
    data = BOSpr,
    ggplot2::aes(
      x = x,
      y = y,
      color = spatial_predictor_100000_31	
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "eigenvalue") +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::ggtitle("spatial_predictor_100000_31") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

BOSp4 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = aus.sf, fill = "white") +
  ggplot2::geom_point(
    data = BOSpr,
    ggplot2::aes(
      x = x,
      y = y,
      color = spatial_predictor_100000_17
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "eigenvalue") +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::ggtitle("spatial_predictor_100000_17") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

(BOSp1 + BOSp2) / (BOSp3 + BOSp4)


# response curves
BOSrespCurvRF <- spatialRF::plot_response_curves(model=BOSmodel.tuned, quantiles = c(0.5), ncol = 3)
BOSrespCurvRF
str(BOSrespCurvRF[[1]])

BOSrespCurvRF[[1]]$labels$x


# ferric2
BOSrespCurvRF[[1]]
BOSrespCurvRF[[1]]$labels$x
BOSferricrespCurvRF <- BOSrespCurvRF[[1]]
str(BOSferricrespCurvRF)

table(BOSferricrespCurvRF$layers[[1]]$data$id)

BOSferricpredRF.dat <- data.frame(iter=BOSferricrespCurvRF$layers[[1]]$data$id, 
                                  ferric=BOSferricrespCurvRF$layers[[1]]$data$ferric,
                                  Hg=BOSferricrespCurvRF$layers[[1]]$data$Hg,
                                  quantile=BOSferricrespCurvRF$layers[[1]]$data$quantile)

head(BOSferricpredRF.dat)
BOSferricpredRFq.5.dat <- subset(BOSferricpredRF.dat, quantile==0.5)
dim(BOSferricpredRFq.5.dat)

BOSferricpredRF.5Xiter.stats <- BOSferricpredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(ferric, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
BOSferricpredRF.5Xiter.stats
range(BOSferricpredRF.5Xiter.stats$x, na.rm=T)
c(min(BOSferricpredRF.5Xiter.stats$lowery, na.rm=T), max(BOSferricpredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(BOSferricpredRF.5Xiter.stats, "BOSferricpredRF.5Xiter.stats.csv", row.names=F)


# soil N
BOSrespCurvRF[[2]]$labels$x
BOSsoilNrespCurvRF <- BOSrespCurvRF[[2]]
str(BOSsoilNrespCurvRF)

table(BOSsoilNrespCurvRF$layers[[1]]$data$id)

BOSsoilNpredRF.dat <- data.frame(iter=BOSsoilNrespCurvRF$layers[[1]]$data$id, 
                                 soilN=BOSsoilNrespCurvRF$layers[[1]]$data$soilN,
                                 Hg=BOSsoilNrespCurvRF$layers[[1]]$data$Hg,
                                 quantile=BOSsoilNrespCurvRF$layers[[1]]$data$quantile)
head(BOSsoilNpredRF.dat)
BOSsoilNpredRFq.5.dat <- subset(BOSsoilNpredRF.dat, quantile==0.5)
dim(BOSsoilNpredRFq.5.dat)

BOSsoilNpredRF.5Xiter.stats <- BOSsoilNpredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(soilN, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
BOSsoilNpredRF.5Xiter.stats
range(BOSsoilNpredRF.5Xiter.stats$x, na.rm=T)
c(min(BOSsoilNpredRF.5Xiter.stats$lowery, na.rm=T), max(BOSsoilNpredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(BOSsoilNpredRF.5Xiter.stats, "BOSsoilNpredRF.5Xiter.stats.csv", row.names=F)

# LAI
BOSrespCurvRF[[3]]$labels$x
BOSlairespCurvRF <- BOSrespCurvRF[[3]]
str(BOSlairespCurvRF)

table(BOSlairespCurvRF$layers[[1]]$data$id)

BOSlaipredRF.dat <- data.frame(iter=BOSlairespCurvRF$layers[[1]]$data$id, 
                               lai=BOSlairespCurvRF$layers[[1]]$data$lai,
                               Hg=BOSlairespCurvRF$layers[[1]]$data$Hg,
                               quantile=BOSlairespCurvRF$layers[[1]]$data$quantile)
head(BOSlaipredRF.dat)
BOSlaipredRFq.5.dat <- subset(BOSlaipredRF.dat, quantile==0.5)
dim(BOSlaipredRFq.5.dat)

BOSlaipredRF.5Xiter.stats <- BOSlaipredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(lai, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
BOSlaipredRF.5Xiter.stats
range(BOSlaipredRF.5Xiter.stats$x, na.rm=T)
c(min(BOSlaipredRF.5Xiter.stats$lowery, na.rm=T), max(BOSlaipredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(BOSlaipredRF.5Xiter.stats, "BOSlaipredRF.5Xiter.stats.csv", row.names=F)

# Prescott
BOSrespCurvRF[[4]]$labels$x
BOSPrescottrespCurvRF <- BOSrespCurvRF[[4]]
str(BOSPrescottrespCurvRF)

table(BOSPrescottrespCurvRF$layers[[1]]$data$id)

BOSPrescottpredRF.dat <- data.frame(iter=BOSPrescottrespCurvRF$layers[[1]]$data$id, 
                                    Prescott=BOSPrescottrespCurvRF$layers[[1]]$data$Prescott,
                                    Hg=BOSPrescottrespCurvRF$layers[[1]]$data$Hg,
                                    quantile=BOSPrescottrespCurvRF$layers[[1]]$data$quantile)
head(BOSPrescottpredRF.dat)
BOSPrescottpredRFq.5.dat <- subset(BOSPrescottpredRF.dat, quantile==0.5)
dim(BOSPrescottpredRFq.5.dat)

BOSPrescottpredRF.5Xiter.stats <- BOSPrescottpredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(Prescott, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
BOSPrescottpredRF.5Xiter.stats
range(BOSPrescottpredRF.5Xiter.stats$x, na.rm=T)
c(min(BOSPrescottpredRF.5Xiter.stats$lowery, na.rm=T), max(BOSPrescottpredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(BOSPrescottpredRF.5Xiter.stats, "BOSPrescottpredRF.5Xiter.stats.csv", row.names=F)

# Pox
BOSrespCurvRF[[5]]$labels$x
BOSPoxrespCurvRF <- BOSrespCurvRF[[5]]
str(BOSPoxrespCurvRF)

table(BOSPoxrespCurvRF$layers[[1]]$data$id)

BOSPoxpredRF.dat <- data.frame(iter=BOSPoxrespCurvRF$layers[[1]]$data$id,
                               Pox=BOSPoxrespCurvRF$layers[[1]]$data$Pox,
                               Hg=BOSPoxrespCurvRF$layers[[1]]$data$Hg,
                               quantile=BOSPoxrespCurvRF$layers[[1]]$data$quantile)
head(BOSPoxpredRF.dat)
BOSPoxpredRFq.5.dat <- subset(BOSPoxpredRF.dat, quantile==0.5)
dim(BOSPoxpredRFq.5.dat)

BOSPoxpredRF.5Xiter.stats <- BOSPoxpredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(Pox, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
BOSPoxpredRF.5Xiter.stats
range(BOSPoxpredRF.5Xiter.stats$x, na.rm=T)
c(min(BOSPoxpredRF.5Xiter.stats$lowery, na.rm=T), max(BOSPoxpredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(BOSPoxpredRF.5Xiter.stats, "BOSPoxpredRF.5Xiter.stats.csv", row.names=F)




# comparing spatial and non-spatial models stochastically
BOScomparison <- spatialRF::rf_compare(
  models = list(
    'non-spatial' = BOSmodel.non.spatial,
    'spatial' = BOSmodel.tuned
  ),
  xy = BOSxy,
  repetitions = 100,
  training.fraction = 0.8,
  metrics = "r.squared"
)

BOSx <- BOScomparison$comparison.df %>% 
  dplyr::group_by(model, metric) %>% 
  dplyr::summarise(value = round(median(value), 3)) %>% 
  dplyr::arrange(metric) %>% 
  as.data.frame()
colnames(BOSx) <- c("model", "metric", "median")
kableExtra::kbl(
  BOSx,
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)





#########################################
# randomised BRT using distance matrix ##
#########################################

## TOS
TOS.brt.dat.rsmp <- na.omit(TOS.test[,c("LON","LAT","lHg","lgs","lgtclay","pH15","lelecond","lLOI","prescott","soilH20","llai","soillN","soillP",
                               "lAl","lMn","lCu","lPb","lSb","lNi","lKThU","ferric2","Alox","Feox","Pox")])

# remove infinite values
TOS.brt.dat.rsmp <- TOS.brt.dat.rsmp[!is.infinite(rowSums(TOS.brt.dat.rsmp)), ]

dim(TOS.brt.dat.rsmp)
head(TOS.brt.dat.rsmp)
TOS.brt.dat.rsmp.coords <- as.data.frame(TOS.brt.dat.rsmp[,c("LON","LAT")])
dim(TOS.brt.dat.rsmp.coords)
head(TOS.brt.dat.rsmp.coords)

# recalculate Haversine distance matrix
coords.TOS.brt.dat.rsmp <- data.frame(
  longitude = TOS.brt.dat.rsmp$LON,
  latitude = TOS.brt.dat.rsmp$LAT
)

TOS.brt.dat.rsmp_haversine_matrix <- distm(
  x = coords.TOS.brt.dat.rsmp,
  fun = distHaversine
)

dim(TOS.brt.dat.rsmp_haversine_matrix)
range(TOS.brt.dat.rsmp_haversine_matrix)

# @param df data frame containing coordinate columns
# @param min_dist numeric value specifying the minimum distance required between points
# @param target_n optional integer specifying desired number of points (default: NULL)
# @param dist_matrix distance matrix between points
# @param x_col character string naming the x-coordinate column (default: "x")
# @param y_col character string naming the y-coordinate column (default: "y")
# @return data frame containing the selected points
select_distant_points <- function(df, min_dist, target_n = NULL, dist_matrix,
                                  x_col = "x", y_col = "y") {
  # input validation
  if (!is.data.frame(df)) {
    stop("input 'df' must be a data frame")
  }
  if (!all(c(x_col, y_col) %in% names(df))) {
    stop("specified coordinate columns not found in data frame")
  }
  if (!is.numeric(min_dist) || min_dist <= 0) {
    stop("min_dist must be a positive number")
  }
  if (!is.null(target_n)) {
    if (!is.numeric(target_n) || target_n <= 0 || target_n != round(target_n)) {
      stop("target_n must be a positive integer")
    }
    if (target_n > nrow(df)) {
      stop("target_n cannot be larger than the number of input points")
    }
  }
  
  # convert coordinates to matrix for distance calculation
  coords <- as.matrix(df[, c(x_col, y_col)])
  
  # Initialize variables
  n_points <- nrow(df)
  selected <- logical(n_points)
  
  # select first point randomly
  current <- sample(n_points, 1)
  selected[current] <- TRUE
  n_selected <- 1
  
  # maximum possible points given distance constraint
  available <- which(!selected)
  
  # selection loop
  while(length(available) > 0) {
    # Find points that satisfy distance constraint
    valid <- available[apply(dist_matrix[selected, available, drop = FALSE], 2,
                             function(x) all(x >= min_dist))]
    
    # break if no valid points remain or target reached
    if (length(valid) == 0 || (!is.null(target_n) && n_selected >= target_n)) {
      break
    }
    
    # randomly select next point from valid candidates
    next_point <- sample(valid, 1)
    selected[next_point] <- TRUE
    n_selected <- n_selected + 1
    
    # update available points
    available <- which(!selected)
  }
  
  # Return selected points
  return(df[selected, , drop = FALSE])
}

# select points with minimum distance of 350 km
min.dist <- 250000 # large enough to remove most spatial autocorrelation, but small enough to ensure adequate points for training BRT
npts2generate <- 100
TOScoords.ran <- select_distant_points(df = TOS.brt.dat.rsmp.coords, min_dist = min.dist, target_n = npts2generate,
                      dist_matrix = TOS.brt.dat.rsmp_haversine_matrix, x_col = "LON", y_col = "LAT")
dim(TOScoords.ran)
head(TOScoords.ran)
TOScoords.ran.pts <- vect(cbind(TOScoords.ran$LON, 
                                TOScoords.ran$LAT), crs="+proj=longlat")
terra::plot(TOScoords.ran.pts)

# subset random points for BRT training data
TOS.ran.subset <- na.omit(TOS.brt.dat.rsmp[as.numeric(row.names(TOScoords.ran)), ])
head(TOS.ran.subset)

## boosted regression tree
TOS.ran.brt <- gbm.step(TOS.ran.subset, gbm.x = attr(TOS.ran.subset, "names")[c(4:length(colnames(TOS.ran.subset)))],
                    gbm.y = attr(TOS.ran.subset, "names")[3], family="gaussian", max.trees=100000,
                    tolerance = 0.0002, learning.rate = 0.0005, bag.fraction=0.75,
                    tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(TOS.ran.brt)
gbm.plot(TOS.ran.brt, smooth=T, rug=T, n.plots=12, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]")

object_exists <- function(obj_name) {
  exists(obj_name, envir = .GlobalEnv)
}

# needed functions
modifyVecFunc <- function(obj_name, index, new_value) {
  tryCatch({
    if (!object_exists(obj_name)) {
      stop("object does not exist: ", obj_name)
    }
    
    obj <- get(obj_name)
    
    # check object type and handle accordingly
    if (is.vector(obj) && !is.list(obj)) {
      if (index > length(obj)) stop("index out of bounds")
      obj[index] <- new_value
    } else if (is.list(obj)) {
      if (index > length(obj)) stop("index out of bounds")
      obj[[index]] <- new_value
    } else if (is.data.frame(obj)) {
      if (!all(index <= dim(obj))) stop("index out of bounds")
      obj[index[1], index[2]] <- new_value
    } else {
      stop("unsupported object type")
    }
    
    # Save modified object back to its original name
    assign(obj_name, obj, envir = .GlobalEnv)
    
  }, error = function(e) {
    message("error: ", e$message)
    return(FALSE)
  })
}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# resampled BRT
biter <- 10000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()
eq.sp.pts <- 100
TOStraincols <- attr(TOS.ran.subset, "names")[c(4:length(colnames(TOS.ran.subset)))] # variable columns used to train data
TOSntraincols <- length(TOStraincols)

# create storage arrays
TOSval.arr <- TOSpred.arr <- array(data=NA, dim=c(eq.sp.pts, TOSntraincols, biter), dimnames=list(paste("x",1:eq.sp.pts,sep=""), TOStraincols, paste("b",1:biter,sep="")))

# create storage vectors
TOSri.vec.names <- paste(TOStraincols,".ri",sep="")
TOSCV.cor.vec <- TOSCV.cor.se.vec <- rep(NA,biter)
for (r in 1:TOSntraincols) {
  assign(TOSri.vec.names[r], rep(NA,biter))}

# b loop
for (b in 1:biter) {
  # generate random coords
  coords.ran <- select_distant_points(df = TOS.brt.dat.rsmp.coords, min_dist = min.dist, target_n = npts2generate,
                                      dist_matrix = TOS.brt.dat.rsmp_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for BRT training data
  TOS.ran.subset <- na.omit(TOS.brt.dat.rsmp[as.numeric(row.names(coords.ran)), ])
  
  ## boosted regression tree
  TOS.ran.brt <- gbm.step(TOS.ran.subset, gbm.x = attr(TOS.ran.subset, "names")[c(4:length(colnames(TOS.ran.subset)))],
                          gbm.y = attr(TOS.ran.subset, "names")[3], family="gaussian", max.trees=100000,
                          tolerance = 0.0002, learning.rate = 0.0005, bag.fraction=0.75,
                          tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  
  # error catch
  if (b == 1 & is.null(TOS.ran.brt)==F) {
    TOS.ran.brt.old <- TOS.ran.brt
  }
  if (is.null(TOS.ran.brt) == T) {
    TOS.ran.brt <- gbm.step(TOS.ran.subset, gbm.x = attr(TOS.ran.subset, "names")[c(4:length(colnames(TOS.ran.subset)))],
                            gbm.y = attr(TOS.ran.subset, "names")[3], family="gaussian", max.trees=100000,
                            tolerance = 0.0002, learning.rate = 0.0005, bag.fraction=0.75,
                            tree.complexity = 2, silent=T, step.size=20, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(TOS.ran.brt) == T) {
    TOS.ran.brt <- TOS.ran.brt.old
  }
  
  # summary
  TOSsumm.fit <- summary(TOS.ran.brt)
  
  if (is.null(TOS.ran.brt) == F) {
    TOS.ran.brt.old <- TOS.ran.brt
  }
  
  # variable relative importance
  for (ri in 1:TOSntraincols) {
    modifyVecFunc(TOSri.vec.names[ri], b, new_value=TOSsumm.fit$rel.inf[which(TOSsumm.fit$var == TOStraincols[ri])])
  }
  
  # goodness of fit
  TOSCV.cor.vec[b] <- 100*TOS.ran.brt$cv.statistics$correlation.mean
  TOSCV.cor.se.vec[b] <- 100*TOS.ran.brt$cv.statistics$correlation.se
  
  # response curves
  TOSRESP.val <- TOSRESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=TOSntraincols)
  for (p in 1:TOSntraincols) {
    TOSRESP.val[,p] <- plot.gbm(TOS.ran.brt, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    TOSRESP.pred[,p] <- plot.gbm(TOS.ran.brt, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  TOSRESP.val.dat <- as.data.frame(TOSRESP.val)
  colnames(TOSRESP.val.dat) <- TOS.ran.brt$var.names
  TOSRESP.pred.dat <- as.data.frame(TOSRESP.pred)
  colnames(TOSRESP.pred.dat) <- TOS.ran.brt$var.names
  
  # add to storage arrays
  TOSval.arr[, , b] <- as.matrix(TOSRESP.val.dat)
  TOSpred.arr[, , b] <- as.matrix(TOSRESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                       "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                            "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                            "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                      units = "mins"), 2), "minutes elapsed")))
} # end b

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- TOSntraincols
TOSpred.update <- TOSpred.arr[,,1:biter]

for (k in 1:kappa.n) {
  TOSboot.mean <- apply(TOSpred.update, MARGIN=c(1,2), mean, na.rm=T)
  TOSboot.sd <- apply(TOSpred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    TOSpred.update[,,z] <- ifelse((TOSpred.update[,,z] < (TOSboot.mean-kappa*TOSboot.sd) | TOSpred.update[,,z] >
                                  (TOSboot.mean+kappa*TOSboot.sd)), NA, TOSpred.update[,,z])
  } # end z
  print(k)
} # end k

TOSpred.med <- apply(TOSpred.update, MARGIN=c(1,2), median, na.rm=T)
TOSpred.lo <- apply(TOSpred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
TOSpred.up <- apply(TOSpred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
TOSval.med <- apply(TOSval.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
TOSCV.cor.update <- TOSCV.cor.vec[1:biter]
TOSCV.cor.se.update <- TOSCV.cor.se.vec[1:biter]

# update ri vectors
TOSri.vec.update.names <- paste(TOSri.vec.names,".update",sep="")
for (ri in 1:TOSntraincols) {
  assign(TOSri.vec.update.names[ri], get(TOSri.vec.names[ri])[1:biter])
}

TOSvec.mean.names <- paste(TOStraincols,".mean",sep="")
TOSvec.sd.names <- paste(TOStraincols,".sd",sep="")

for (k in 1:kappa.n) {
  TOSCV.cor.mean <- mean(TOSCV.cor.update, na.rm=T); TOSCV.cor.sd <- sd(TOSCV.cor.update, na.rm=T)
  TOSCV.cor.se.mean <- mean(TOSCV.cor.se.update, na.rm=T); TOSCV.cor.se.sd <- sd(TOSCV.cor.se.update, na.rm=T)
  
  for (v in 1:TOSntraincols) {
    assign(TOSvec.mean.names[v], mean(get(TOSri.vec.update.names[v]), na.rm=T))
    assign(TOSvec.sd.names[v], sd(get(TOSri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    TOSCV.cor.update[u] <- ifelse((TOSCV.cor.update[u] < (TOSCV.cor.mean-kappa*TOSCV.cor.sd) | TOSCV.cor.update[u] >
                                  (TOSCV.cor.mean+kappa*TOSCV.cor.sd)), NA, TOSCV.cor.update[u])
    TOSCV.cor.se.update[u] <- ifelse((TOSCV.cor.se.update[u] < (TOSCV.cor.se.mean-kappa*TOSCV.cor.se.sd) | TOSCV.cor.se.update[u] >
                                    (TOSCV.cor.se.mean+kappa*TOSCV.cor.se.sd)), NA, TOSCV.cor.se.update[u])
    for (ri in 1:TOSntraincols) {
      modifyVecFunc(TOSri.vec.update.names[ri], u, ifelse((get(TOSri.vec.update.names[ri])[u]) < 
                                                       (get(TOSvec.mean.names[ri]) - kappa*get(TOSvec.sd.names[ri])),
                                                        NA, get(TOSri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
TOSCV.cor.med <- median(TOSCV.cor.update, na.rm=T)
TOSCV.cor.lo <- quantile(TOSCV.cor.update, probs=0.025, na.rm=T)
TOSCV.cor.up <- quantile(TOSCV.cor.update, probs=0.975, na.rm=T)
print(c(TOSCV.cor.lo, TOSCV.cor.med, TOSCV.cor.up))

TOSri.vec.lo.names <- paste(TOStraincols,".ri.lo",sep="")
TOSri.vec.up.names <- paste(TOStraincols,".ri.up",sep="")
TOSri.vec.med.names <- paste(TOStraincols,".ri.med",sep="")

for (ri in 1:TOSntraincols) {
  assign(TOSri.vec.lo.names[ri], quantile(get(TOSri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(TOSri.vec.med.names[ri], median(get(TOSri.vec.update.names[ri]), na.rm=T))
  assign(TOSri.vec.up.names[ri], quantile(get(TOSri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

TOSri.lo <- as.numeric(mget(TOSri.vec.lo.names))
TOSri.med <- as.numeric(mget(TOSri.vec.med.names))
TOSri.up <- as.numeric(mget(TOSri.vec.up.names))
TOSri.out <- as.data.frame(cbind(TOSri.med, TOSri.up, TOSri.lo))
rownames(TOSri.out) <- TOStraincols
TOSri.sort <- TOSri.out[order(TOSri.out[,1], decreasing=T),]
TOSri.sort

# plot
TOSri.plt <- ggplot(TOSri.sort) +
  geom_bar(aes(x=reorder(row.names(TOSri.sort), TOSri.med), y=TOSri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(TOSri.sort), ymin=TOSri.lo, ymax=TOSri.up),
                 linewidth=0.4, colour="black", alpha=0.9)
TOSri.plt + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(TOSCV.cor.lo,TOSCV.cor.med,TOSCV.cor.up), 2))

## plot predicted relationships of top x variables
TOStopNvars <- 9 # x
head(TOSpred.med)
TOSri.sort
TOStop.ri.sort <- TOSri.sort[1:TOStopNvars,]
TOStopNvars.names <- rownames(TOStop.ri.sort)
TOSylims <- c(min(TOSpred.lo[,TOStopNvars.names], na.rm=T), max(TOSpred.up[,TOStopNvars.names], na.rm=T))

TOSplotNvec <- paste("TOSplt",1:TOStopNvars,sep="")
for (v in 1:TOStopNvars) {
  assign(TOSplotNvec[v], ggplot(data=as.data.frame(cbind(TOSval.med[,TOStopNvars.names[v]], TOSpred.med[,TOStopNvars.names[v]],
                                                         TOSpred.lo[,TOStopNvars.names[v]], TOSpred.up[,TOStopNvars.names[v]]))) +
    geom_line(aes(x=V1, y=V2), colour="blue") +
    geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
    lims(y=TOSylims) +
    xlab(TOStopNvars.names[v]) + ylab("log10 [Hg]"))
}
(TOSplt1 + TOSplt2 + TOSplt3) / (TOSplt4 + TOSplt5 + TOSplt6) / (TOSplt7 + TOSplt8 + TOSplt9)

# export results
for (v in 1:TOSntraincols) {
  TOSdata <- data.frame(x=TOSval.med[,TOStraincols[v]], mn=TOSpred.med[,TOStraincols[v]],
                         up=TOSpred.up[,TOStraincols[v]], lo=TOSpred.lo[,TOStraincols[v]])
  row.names(TOSdata) <- NULL
  write.csv(TOSdata, file=paste(TOStraincols[v], "TOSPred.csv", sep=""), row.names=F)
}



## BOS
BOS.brt.dat.rsmp <- na.omit(BOS.test[,c("LON","LAT","lHg","lgs","lgtclay","pH15","lelecond","lLOI","prescott","soilH20","llai","soillN","soillP",
                                        "lAl","lMn","lCu","lPb","lSb","lNi","lKThU","ferric2","Alox","Feox","Pox")])

# remove infinite values
BOS.brt.dat.rsmp <- BOS.brt.dat.rsmp[!is.infinite(rowSums(BOS.brt.dat.rsmp)), ]

dim(BOS.brt.dat.rsmp)
head(BOS.brt.dat.rsmp)
BOS.brt.dat.rsmp.coords <- as.data.frame(BOS.brt.dat.rsmp[,c("LON","LAT")])
dim(BOS.brt.dat.rsmp.coords)
head(BOS.brt.dat.rsmp.coords)

# recalculate Haversine distance matrix
coords.BOS.brt.dat.rsmp <- data.frame(
  longitude = BOS.brt.dat.rsmp$LON,
  latitude = BOS.brt.dat.rsmp$LAT
)

BOS.brt.dat.rsmp_haversine_matrix <- distm(
  x = coords.BOS.brt.dat.rsmp,
  fun = distHaversine
)

dim(BOS.brt.dat.rsmp_haversine_matrix)
range(BOS.brt.dat.rsmp_haversine_matrix)


# select points with minimum distance of 350 km
BOScoords.ran <- select_distant_points(df = BOS.brt.dat.rsmp.coords, min_dist = min.dist, target_n = npts2generate,
                                       dist_matrix = BOS.brt.dat.rsmp_haversine_matrix, x_col = "LON", y_col = "LAT")
dim(BOScoords.ran)
head(BOScoords.ran)
BOScoords.ran.pts <- vect(cbind(BOScoords.ran$LON, 
                                BOScoords.ran$LAT), crs="+proj=longlat")
terra::plot(BOScoords.ran.pts)

# subset random points for BRT training data
BOS.ran.subset <- na.omit(BOS.brt.dat.rsmp[as.numeric(row.names(BOScoords.ran)), ])
head(BOS.ran.subset)

## boosted regression tree
BOS.ran.brt <- gbm.step(BOS.ran.subset, gbm.x = attr(BOS.ran.subset, "names")[c(4:length(colnames(BOS.ran.subset)))],
                        gbm.y = attr(BOS.ran.subset, "names")[3], family="gaussian", max.trees=100000,
                        tolerance = 0.0002, learning.rate = 0.0005, bag.fraction=0.75,
                        tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(BOS.ran.brt)
gbm.plot(BOS.ran.brt, smooth=T, rug=T, n.plots=12, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]")

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# resampled BRT
biter <- 10000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()
eq.sp.pts <- 100
BOStraincols <- attr(BOS.ran.subset, "names")[c(4:length(colnames(BOS.ran.subset)))] # variable columns used to train data
BOSntraincols <- length(BOStraincols)

# create storage arrays
BOSval.arr <- BOSpred.arr <- array(data=NA, dim=c(eq.sp.pts, BOSntraincols, biter), dimnames=list(paste("x",1:eq.sp.pts,sep=""), BOStraincols, paste("b",1:biter,sep="")))

# create storage vectors
BOSri.vec.names <- paste(BOStraincols,".ri",sep="")
BOSCV.cor.vec <- BOSCV.cor.se.vec <- rep(NA,biter)
for (r in 1:BOSntraincols) {
  assign(BOSri.vec.names[r], rep(NA,biter))}

# b loop
for (b in 1:biter) {
  # generate random coords
  coords.ran <- select_distant_points(df = BOS.brt.dat.rsmp.coords, min_dist = min.dist, target_n = npts2generate,
                                      dist_matrix = BOS.brt.dat.rsmp_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for BRT training data
  BOS.ran.subset <- na.omit(BOS.brt.dat.rsmp[as.numeric(row.names(coords.ran)), ])
  
  ## boosted regression tree
  BOS.ran.brt <- gbm.step(BOS.ran.subset, gbm.x = attr(BOS.ran.subset, "names")[c(4:length(colnames(BOS.ran.subset)))],
                          gbm.y = attr(BOS.ran.subset, "names")[3], family="gaussian", max.trees=100000,
                          tolerance = 0.0002, learning.rate = 0.0005, bag.fraction=0.75,
                          tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  
  # error catch
  if (b == 1 & is.null(BOS.ran.brt)==F) {
    BOS.ran.brt.old <- BOS.ran.brt
  }
  if (is.null(BOS.ran.brt) == T) {
    BOS.ran.brt <- gbm.step(BOS.ran.subset, gbm.x = attr(BOS.ran.subset, "names")[c(4:length(colnames(BOS.ran.subset)))],
                            gbm.y = attr(BOS.ran.subset, "names")[3], family="gaussian", max.trees=100000,
                            tolerance = 0.0002, learning.rate = 0.0005, bag.fraction=0.75,
                            tree.complexity = 2, silent=T, step.size=20, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(BOS.ran.brt) == T) {
    BOS.ran.brt <- BOS.ran.brt.old
  }
  
  # summary
  BOSsumm.fit <- summary(BOS.ran.brt)
  
  if (is.null(BOS.ran.brt) == F) {
    BOS.ran.brt.old <- BOS.ran.brt
  }
  
  # variable relative importance
  for (ri in 1:BOSntraincols) {
    modifyVecFunc(BOSri.vec.names[ri], b, new_value=BOSsumm.fit$rel.inf[which(BOSsumm.fit$var == BOStraincols[ri])])
  }
  
  # goodness of fit
  BOSCV.cor.vec[b] <- 100*BOS.ran.brt$cv.statistics$correlation.mean
  BOSCV.cor.se.vec[b] <- 100*BOS.ran.brt$cv.statistics$correlation.se
  
  # response curves
  BOSRESP.val <- BOSRESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=BOSntraincols)
  for (p in 1:BOSntraincols) {
    BOSRESP.val[,p] <- plot.gbm(BOS.ran.brt, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    BOSRESP.pred[,p] <- plot.gbm(BOS.ran.brt, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  BOSRESP.val.dat <- as.data.frame(BOSRESP.val)
  colnames(BOSRESP.val.dat) <- BOS.ran.brt$var.names
  BOSRESP.pred.dat <- as.data.frame(BOSRESP.pred)
  colnames(BOSRESP.pred.dat) <- BOS.ran.brt$var.names
  
  # add to storage arrays
  BOSval.arr[, , b] <- as.matrix(BOSRESP.val.dat)
  BOSpred.arr[, , b] <- as.matrix(BOSRESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- BOSntraincols
BOSpred.update <- BOSpred.arr[,,1:biter]

for (k in 1:kappa.n) {
  BOSboot.mean <- apply(BOSpred.update, MARGIN=c(1,2), mean, na.rm=T)
  BOSboot.sd <- apply(BOSpred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    BOSpred.update[,,z] <- ifelse((BOSpred.update[,,z] < (BOSboot.mean-kappa*BOSboot.sd) | BOSpred.update[,,z] >
                                     (BOSboot.mean+kappa*BOSboot.sd)), NA, BOSpred.update[,,z])
  } # end z
  print(k)
} # end k

BOSpred.med <- apply(BOSpred.update, MARGIN=c(1,2), median, na.rm=T)
BOSpred.lo <- apply(BOSpred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
BOSpred.up <- apply(BOSpred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
BOSval.med <- apply(BOSval.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
BOSCV.cor.update <- BOSCV.cor.vec[1:biter]
BOSCV.cor.se.update <- BOSCV.cor.se.vec[1:biter]

# update ri vectors
BOSri.vec.update.names <- paste(BOSri.vec.names,".update",sep="")
for (ri in 1:BOSntraincols) {
  assign(BOSri.vec.update.names[ri], get(BOSri.vec.names[ri])[1:biter])
}

BOSvec.mean.names <- paste(BOStraincols,".mean",sep="")
BOSvec.sd.names <- paste(BOStraincols,".sd",sep="")

for (k in 1:kappa.n) {
  BOSCV.cor.mean <- mean(BOSCV.cor.update, na.rm=T); BOSCV.cor.sd <- sd(BOSCV.cor.update, na.rm=T)
  BOSCV.cor.se.mean <- mean(BOSCV.cor.se.update, na.rm=T); BOSCV.cor.se.sd <- sd(BOSCV.cor.se.update, na.rm=T)
  
  for (v in 1:BOSntraincols) {
    assign(BOSvec.mean.names[v], mean(get(BOSri.vec.update.names[v]), na.rm=T))
    assign(BOSvec.sd.names[v], sd(get(BOSri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    BOSCV.cor.update[u] <- ifelse((BOSCV.cor.update[u] < (BOSCV.cor.mean-kappa*BOSCV.cor.sd) | BOSCV.cor.update[u] >
                                     (BOSCV.cor.mean+kappa*BOSCV.cor.sd)), NA, BOSCV.cor.update[u])
    BOSCV.cor.se.update[u] <- ifelse((BOSCV.cor.se.update[u] < (BOSCV.cor.se.mean-kappa*BOSCV.cor.se.sd) | BOSCV.cor.se.update[u] >
                                        (BOSCV.cor.se.mean+kappa*BOSCV.cor.se.sd)), NA, BOSCV.cor.se.update[u])
    for (ri in 1:BOSntraincols) {
      modifyVecFunc(BOSri.vec.update.names[ri], u, ifelse((get(BOSri.vec.update.names[ri])[u]) < 
                                                            (get(BOSvec.mean.names[ri]) - kappa*get(BOSvec.sd.names[ri])),
                                                          NA, get(BOSri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
BOSCV.cor.med <- median(BOSCV.cor.update, na.rm=T)
BOSCV.cor.lo <- quantile(BOSCV.cor.update, probs=0.025, na.rm=T)
BOSCV.cor.up <- quantile(BOSCV.cor.update, probs=0.975, na.rm=T)
print(c(BOSCV.cor.lo, BOSCV.cor.med, BOSCV.cor.up))

BOSri.vec.lo.names <- paste(BOStraincols,".ri.lo",sep="")
BOSri.vec.up.names <- paste(BOStraincols,".ri.up",sep="")
BOSri.vec.med.names <- paste(BOStraincols,".ri.med",sep="")

for (ri in 1:BOSntraincols) {
  assign(BOSri.vec.lo.names[ri], quantile(get(BOSri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(BOSri.vec.med.names[ri], median(get(BOSri.vec.update.names[ri]), na.rm=T))
  assign(BOSri.vec.up.names[ri], quantile(get(BOSri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

BOSri.lo <- as.numeric(mget(BOSri.vec.lo.names))
BOSri.med <- as.numeric(mget(BOSri.vec.med.names))
BOSri.up <- as.numeric(mget(BOSri.vec.up.names))
BOSri.out <- as.data.frame(cbind(BOSri.med, BOSri.up, BOSri.lo))
rownames(BOSri.out) <- BOStraincols
BOSri.sort <- BOSri.out[order(BOSri.out[,1], decreasing=T),]
BOSri.sort

# plot
BOSri.plt <- ggplot(BOSri.sort) +
  geom_bar(aes(x=reorder(row.names(BOSri.sort), BOSri.med), y=BOSri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(BOSri.sort), ymin=BOSri.lo, ymax=BOSri.up),
                linewidth=0.4, colour="black", alpha=0.9)
BOSri.plt + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(BOSCV.cor.lo,BOSCV.cor.med,BOSCV.cor.up), 2))

## plot predicted relationships of top x variables
BOStopNvars <- 9 # x
head(BOSpred.med)
BOSri.sort
BOStop.ri.sort <- BOSri.sort[1:BOStopNvars,]
BOStopNvars.names <- rownames(BOStop.ri.sort)
BOSylims <- c(min(BOSpred.lo[,BOStopNvars.names], na.rm=T), max(BOSpred.up[,BOStopNvars.names], na.rm=T))

BOSplotNvec <- paste("BOSplt",1:BOStopNvars,sep="")
for (v in 1:BOStopNvars) {
  assign(BOSplotNvec[v], ggplot(data=as.data.frame(cbind(BOSval.med[,BOStopNvars.names[v]], BOSpred.med[,BOStopNvars.names[v]],
                                                         BOSpred.lo[,BOStopNvars.names[v]], BOSpred.up[,BOStopNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=BOSylims) +
           xlab(BOStopNvars.names[v]) + ylab("log10 [Hg]"))
}
(BOSplt1 + BOSplt2 + BOSplt3) / (BOSplt4 + BOSplt5 + BOSplt6) / (BOSplt7 + BOSplt8 + BOSplt9)

# export results
for (v in 1:BOSntraincols) {
  BOSdata <- data.frame(x=BOSval.med[,BOStraincols[v]], mn=BOSpred.med[,BOStraincols[v]],
                        up=BOSpred.up[,BOStraincols[v]], lo=BOSpred.lo[,BOStraincols[v]])
  row.names(BOSdata) <- NULL
  write.csv(BOSdata, file=paste(BOStraincols[v], "BOSPred.csv", sep=""), row.names=F)
}



##############################################################
## resampled for > min distance to factor-level differences ##
##############################################################

## functions
AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC

contains_zero_sign <- function(x_min, x_max) {
  return(sign(x_min) * sign(x_max) <= 0)
}

## TOS
colnames(TOS.dat)
TOS.cl.test <- data.frame(SITEID=TOS.dat$SITEID, LAT=TOS.dat$LAT, LON=TOS.dat$LON, state=TOS.dat$STATE,
                          lHg=log10(TOS.dat$HgCOMP), lucat=TOS.dat$lucat, landuse=TOS.dat$landuse,
                          biome=TOS.dat$biome, geol=TOS.dat$lithgrp, soil=TOS.dat$soils)
head(TOS.cl.test)

table(TOS.cl.test$state)
sum(table(TOS.cl.test$state))
table(TOS.cl.test$lucat)
sum(table(TOS.cl.test$lucat))
table(TOS.cl.test$landuse)
sum(table(TOS.cl.test$landuse))
table(TOS.cl.test$geol)
sum(table(TOS.cl.test$geol))
which(is.na(TOS.cl.test$geol)==T)
table(TOS.cl.test$soil)
sum(table(TOS.cl.test$soil))

dim(TOS.cl.test)
dim(na.omit(TOS.cl.test))

TOS.lm <- na.omit(TOS.cl.test)
dim(TOS.lm)
head(TOS.lm)
TOS.lm.coords <- as.data.frame(TOS.lm[,c("LON","LAT")])
dim(TOS.lm.coords)
head(TOS.lm.coords)
TOS.lm.coords.mat <- as.matrix(TOS.lm.coords)
colnames(TOS.lm.coords.mat) <- c("x","y")

# haversine matrix
TOS.lm_haversine_matrix <- distm(
  x = TOS.lm.coords,
  fun = distHaversine
)
dim(TOS.lm_haversine_matrix)
range(TOS.lm_haversine_matrix)

min.dist <- 250000

TOScoords.ran <- select_distant_points(df = TOS.lm.coords, min_dist = min.dist,
                                    dist_matrix = TOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
dim(TOScoords.ran)
head(TOScoords.ran)
TOScoords.ran.pts <- vect(cbind(TOScoords.ran$LON, 
                                TOScoords.ran$LAT), crs="+proj=longlat")
terra::plot(TOScoords.ran.pts)

# subset random points
TOS.ran.subset <- na.omit(TOS.lm[as.numeric(row.names(coords.ran)), ])
head(TOS.ran.subset)
dim(TOS.ran.subset)

# by state
table(TOS.ran.subset$state)
TOSHgXstate.stats.ran <- TOS.ran.subset %>%
  group_by(state) %>%
  summarise(
    mean = mean(lHg, na.rm = TRUE),
    median = median(lHg, na.rm = TRUE), 
    var = var(lHg, na.rm = TRUE),
    sd = sd(lHg, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(lHg, probs=0.975, na.rm = TRUE),
    lower = quantile(lHg, probs=0.025, na.rm = TRUE),
    n = n()
  )
TOSHgXstate.stats.ran
na.omit(TOSHgXstate.stats.ran)
hist(TOS.ran.subset$lHg)


TOSmod1 <- lm(lHg ~ state, data=TOS.ran.subset)
TOSmod.null <- lm(lHg ~ 1, data=TOS.ran.subset)
check_model(TOSmod1, detrend=T)
plot_model(TOSmod1)
TOSpmmod1 <- plot_model(TOSmod1)
TOSpmmod1.coef <- data.frame(state=as.character(TOSpmmod1[[1]]$term), lo=TOSpmmod1[[1]]$conf.low, up=TOSpmmod1[[1]]$conf.high)

TOSpmmod1.coef$nonzero <- ifelse(contains_zero_sign(TOSpmmod1.coef$lo, TOSpmmod1.coef$up)==F, 1, 0)
TOSpmmod1.coef
TOSnonzerosum <- sum(TOSpmmod1.coef$nonzero)
TOSpmmod1.coef[which(TOSpmmod1.coef$nonzero==1),]$state

TOSwAICc <- weight.IC(delta.IC(c(AICc(TOSmod1),AICc(TOSmod.null))))
TOSER <- TOSwAICc[1]/TOSwAICc[2]
TOSER


# resampling loop
iter <- 1000
itdiv <- iter/10

# storage matrix
table(TOS.lm$state)
TOSnlevels <- length(table(TOS.lm$state))
TOSstor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*TOSnlevels)
TOSstatedir.lab <- paste("state",attr(table(TOS.lm$state), "names"),"dir",sep="")
colnames(TOSstor.mat) <- c("ER", "nonzerosum",
                        paste("state",attr(table(TOS.lm$state), "names"),sep=""),TOSstatedir.lab)
head(TOSstor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = TOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = TOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  TOS.ran.subset <- na.omit(TOS.lm[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  TOSmod1 <- lm(lHg ~ state, data=TOS.ran.subset) # class level model
  TOSmod.null <- lm(lHg ~ 1, data=TOS.ran.subset) # null model
  
  # coefficient boundaries
  TOSpmmod1 <- plot_model(TOSmod1)
  TOSpmmod1.coef <- data.frame(state=as.character(TOSpmmod1[[1]]$term), lo=TOSpmmod1[[1]]$conf.low,
                            up=TOSpmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  TOSpmmod1.coef$nonzero <- ifelse(contains_zero_sign(TOSpmmod1.coef$lo, TOSpmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  TOSpmmod1.coef$dir <- ifelse(sign(TOSpmmod1.coef$lo) == -1 & sign(TOSpmmod1.coef$up) == -1, -1, 0)
  TOSpmmod1.coef$dir <- ifelse(sign(TOSpmmod1.coef$lo) == 1 & sign(TOSpmmod1.coef$up) == 1, 1, TOSpmmod1.coef$dir)
    
  # how many variables are non-zero?
  TOSstor.mat[i,"nonzerosum"] <- sum(TOSpmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  TOSnzstates <- TOSpmmod1.coef[which(TOSpmmod1.coef$nonzero==1),]$state
  TOSstor.mat[i, which(colnames(TOSstor.mat) %in% TOSnzstates)] <- 1
  
  # for non-zeros, what is the direction?
  TOSnzstatesdir <- paste(TOSnzstates,"dir",sep="")
  TOSstor.mat[i, which(colnames(TOSstor.mat) %in% TOSnzstatesdir)] <- TOSpmmod1.coef$dir[which(TOSpmmod1.coef$dir != 0)]
  
  # AICc model comparison
  TOSwAICc <- weight.IC(delta.IC(c(AICc(TOSmod1),AICc(TOSmod.null))))

  # store evidence ratio
  TOSstor.mat[i,"ER"] <- TOSwAICc[1]/TOSwAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

# ER
median(TOSstor.mat[,"ER"])
quantile(TOSstor.mat[,"ER"], probs=c(0.025,0.975))

head(TOSstor.mat)
state.labs <- c("NSW","NT","QLD","SA","TAS","VIC","WA")
TOSlennzNSW <- length(table(TOSstor.mat[,"stateNSW"]))
TOSlennzNT <- length(table(TOSstor.mat[,"stateNT"]))
TOSlennzQLD <- length(table(TOSstor.mat[,"stateQLD"]))
TOSlennzSA <- length(table(TOSstor.mat[,"stateSA"]))
TOSlennzTAS <- length(table(TOSstor.mat[,"stateTAS"]))
TOSlennzVIC <- length(table(TOSstor.mat[,"stateVIC"]))
TOSlennzWA <- length(table(TOSstor.mat[,"stateWA"]))

TOSnzrslts <- c(ifelse(TOSlennzNSW==1, as.numeric(table(TOSstor.mat[,"stateNSW"])), 0),
  ifelse(TOSlennzNT==1, as.numeric(table(TOSstor.mat[,"stateNT"])), 0),
  ifelse(TOSlennzQLD==1, as.numeric(table(TOSstor.mat[,"stateQLD"])), 0),
  ifelse(TOSlennzSA==1, as.numeric(table(TOSstor.mat[,"stateSA"])), 0),
  ifelse(TOSlennzTAS==1, as.numeric(table(TOSstor.mat[,"stateTAS"])), 0),
  ifelse(TOSlennzVIC==1, as.numeric(table(TOSstor.mat[,"stateVIC"])), 0),
  ifelse(TOSlennzWA==1, as.numeric(table(TOSstor.mat[,"stateWA"])), 0))

TOSdirNSWtab <- table(TOSstor.mat[,"stateNSWdir"])
TOSdirNTtab <- table(TOSstor.mat[,"stateNTdir"])
TOSdirQLDtab <- table(TOSstor.mat[,"stateQLDdir"])
TOSdirSAtab <- table(TOSstor.mat[,"stateSAdir"])
TOSdirTAStab <- table(TOSstor.mat[,"stateTASdir"])
TOSdirVICtab <- table(TOSstor.mat[,"stateVICdir"])
TOSdirWAtab <- table(TOSstor.mat[,"stateWAdir"])

TOSdirlenNSW <- length(TOSdirNSWtab)
TOSdirlenNT <- length(TOSdirNTtab)
TOSdirlenQLD <- length(TOSdirQLDtab)
TOSdirlenSA <- length(TOSdirSAtab)
TOSdirlenTAS <- length(TOSdirTAStab)
TOSdirlenVIC <- length(TOSdirVICtab)
TOSdirlenWA <- length(TOSdirWAtab)

if (TOSdirlenNSW == 0) {
  TOSNSWpos <- 0
  TOSNSWneg <- 0
} else if (TOSdirlenNSW == 1 & as.numeric(attr(TOSdirNSWtab,"names")[1]) == 1) {
  TOSNSWpos <- as.numeric(TOSdirNSWtab)[1]
  TOSNSWneg <- 0
} else if (TOSdirlenNSW == 1 & as.numeric(attr(TOSdirNSWtab,"names")[1]) == -1) {
  TOSNSWpos <- 0
  TOSNSWneg <- as.numeric(TOSdirNSWtab)[1]
} else if (TOSdirlenNSW == 2) {
  TOSNSWneg <- as.numeric(TOSdirNSWtab)[1]
  TOSNSWpos <- as.numeric(TOSdirNSWtab)[2]
} else {
  TOSNSWpos <- 0
  TOSNSWneg <- 0
}

if (TOSdirlenNT == 0) {
  TOSNTpos <- 0
  TOSNTneg <- 0
  } else if (TOSdirlenNT == 1 & as.numeric(attr(TOSdirNTtab,"names")[1]) == 1) {
    TOSNTpos <- as.numeric(TOSdirNTtab)[1]
    TOSNTneg <- 0
  } else if (TOSdirlenNT == 1 & as.numeric(attr(TOSdirNTtab,"names")[1]) == -1) {
    TOSNTpos <- 0
    TOSNTneg <- as.numeric(TOSdirNTtab)[1]
  } else if (TOSdirlenNT == 2) {
    TOSNTneg <- as.numeric(TOSdirNTtab)[1]
    TOSNTpos <- as.numeric(TOSdirNTtab)[2]
  } else {
    TOSNTpos <- 0
    TOSNTneg <- 0
  }

if (TOSdirlenQLD == 0) {
  TOSQLDpos <- 0
  TOSQLDneg <- 0
  } else if (TOSdirlenQLD == 1 & as.numeric(attr(TOSdirQLDtab,"names")[1]) == 1) {
    TOSQLDpos <- as.numeric(TOSdirQLDtab)[1]
    TOSQLDneg <- 0
  } else if (TOSdirlenQLD == 1 & as.numeric(attr(TOSdirQLDtab,"names")[1]) == -1) {
    TOSQLDpos <- 0
    TOSQLDneg <- as.numeric(TOSdirQLDtab)[1]
  } else if (TOSdirlenQLD == 2) {
    TOSQLDneg <- as.numeric(TOSdirQLDtab)[1]
    TOSQLDpos <- as.numeric(TOSdirQLDtab)[2]
  } else {
    TOSQLDpos <- 0
    TOSQLDneg <- 0
  }

if (TOSdirlenSA == 0) {
  TOSSApos <- 0
  TOSSAneg <- 0
  } else if (TOSdirlenSA == 1 & as.numeric(attr(TOSdirSAtab,"names")[1]) == 1) {
    TOSSApos <- as.numeric(TOSdirSAtab)[1]
    TOSSAneg <- 0
  } else if (TOSdirlenSA == 1 & as.numeric(attr(TOSdirSAtab,"names")[1]) == -1) {
    TOSSApos <- 0
    TOSSAneg <- as.numeric(TOSdirSAtab)[1]
  } else if (TOSdirlenSA == 2) {
    TOSSAneg <- as.numeric(TOSdirSAtab)[1]
    TOSSApos <- as.numeric(TOSdirSAtab)[2]
  } else {
    TOSSApos <- 0
    TOSSAneg <- 0
  }

if (TOSdirlenTAS == 0) {
  TOSTASpos <- 0
  TOSTASneg <- 0
  } else if (TOSdirlenTAS == 1 & as.numeric(attr(TOSdirTAStab,"names")[1]) == 1) {
    TOSTASpos <- as.numeric(TOSdirTAStab)[1]
    TOSTASneg <- 0
  } else if (TOSdirlenTAS == 1 & as.numeric(attr(TOSdirTAStab,"names")[1]) == -1) {
    TOSTASpos <- 0
    TOSTASneg <- as.numeric(TOSdirTAStab)[1]
  } else if (TOSdirlenTAS == 2) {
    TOSTASneg <- as.numeric(TOSdirTAStab)[1]
    TOSTASpos <- as.numeric(TOSdirTAStab)[2]
  } else {
    TOSTASpos <- 0
    TOSTASneg <- 0
  }

if (TOSdirlenVIC == 0) {
  TOSVICpos <- 0
  TOSVICneg <- 0
  } else if (TOSdirlenVIC == 1 & as.numeric(attr(TOSdirVICtab,"names")[1]) == 1) {
    TOSVICpos <- as.numeric(TOSdirVICtab)[1]
    TOSVICneg <- 0
  } else if (TOSdirlenVIC == 1 & as.numeric(attr(TOSdirVICtab,"names")[1]) == -1) {
    TOSVICpos <- 0
    TOSVICneg <- as.numeric(TOSdirVICtab)[1]
  } else if (TOSdirlenVIC == 2) {
    TOSVICneg <- as.numeric(TOSdirVICtab)[1]
    TOSVICpos <- as.numeric(TOSdirVICtab)[2]
  } else {
    TOSVICpos <- 0
    TOSVICneg <- 0
  }

if (TOSdirlenWA == 0) {
  TOSWApos <- 0
  TOSWAneg <- 0
  } else if (TOSdirlenWA == 1 & as.numeric(attr(TOSdirWAtab,"names")[1]) == 1) {
    TOSWApos <- as.numeric(TOSdirWAtab)[1]
    TOSWAneg <- 0
  } else if (TOSdirlenWA == 1 & as.numeric(attr(TOSdirWAtab,"names")[1]) == -1) {
    TOSWApos <- 0
    TOSWAneg <- as.numeric(TOSdirWAtab)[1]
  } else if (TOSdirlenWA == 2) {
    TOSWAneg <- as.numeric(TOSdirWAtab)[1]
    TOSWApos <- as.numeric(TOSdirWAtab)[2]
  } else {
    TOSWApos <- 0
    TOSWAneg <- 0
  }


TOSnznegdirrslts <- c(TOSNSWneg, TOSNTneg, TOSQLDneg, TOSSAneg, TOSTASneg, TOSVICneg, TOSWAneg)
TOSnzposdirrslts <- c(TOSNSWpos, TOSNTpos, TOSQLDpos, TOSSApos, TOSTASpos, TOSVICpos, TOSWApos)

TOSresults.out <- data.frame(state=state.labs, nonzero=TOSnzrslts/iter, negdir=TOSnznegdirrslts/iter,
                          posdir=TOSnzposdirrslts/iter)
TOSresults.out

# ER
10^median(log10(TOSstor.mat[,"ER"]))
quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(TOSstor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(TOSstor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(TOSstor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(TOSstor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)


## BOS
colnames(BOS.dat)
BOS.cl.test <- data.frame(SITEID=BOS.dat$SITEID, LAT=BOS.dat$LAT, LON=BOS.dat$LON, state=BOS.dat$STATE,
                          lHg=log10(BOS.dat$HgCOMP), lucat=BOS.dat$lucat, landuse=BOS.dat$landuse,
                          biome=BOS.dat$biome, geol=BOS.dat$lithgrp, soil=BOS.dat$soils)
head(BOS.cl.test)

table(BOS.cl.test$state)
sum(table(BOS.cl.test$state))
table(BOS.cl.test$lucat)
sum(table(BOS.cl.test$lucat))
table(BOS.cl.test$landuse)
sum(table(BOS.cl.test$landuse))
table(BOS.cl.test$geol)
sum(table(BOS.cl.test$geol))
which(is.na(BOS.cl.test$geol)==T)
table(BOS.cl.test$soil)
sum(table(BOS.cl.test$soil))

dim(BOS.cl.test)
dim(na.omit(BOS.cl.test))

BOS.lm <- na.omit(BOS.cl.test)
dim(BOS.lm)
head(BOS.lm)
BOS.lm.coords <- as.data.frame(BOS.lm[,c("LON","LAT")])
dim(BOS.lm.coords)
head(BOS.lm.coords)
BOS.lm.coords.mat <- as.matrix(BOS.lm.coords)
colnames(BOS.lm.coords.mat) <- c("x","y")

# haversine matrix
BOS.lm_haversine_matrix <- distm(
  x = BOS.lm.coords,
  fun = distHaversine
)
dim(BOS.lm_haversine_matrix)
range(BOS.lm_haversine_matrix)

BOScoords.ran <- select_distant_points(df = BOS.lm.coords, min_dist = min.dist,
                                       dist_matrix = BOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
dim(BOScoords.ran)
head(BOScoords.ran)
BOScoords.ran.pts <- vect(cbind(BOScoords.ran$LON, 
                                BOScoords.ran$LAT), crs="+proj=longlat")
terra::plot(BOScoords.ran.pts)

# subset random points
BOS.ran.subset <- na.omit(BOS.lm[as.numeric(row.names(coords.ran)), ])
head(BOS.ran.subset)
dim(BOS.ran.subset)

# by state
table(BOS.ran.subset$state)
BOSHgXstate.stats.ran <- BOS.ran.subset %>%
  group_by(state) %>%
  summarise(
    mean = mean(lHg, na.rm = TRUE),
    median = median(lHg, na.rm = TRUE), 
    var = var(lHg, na.rm = TRUE),
    sd = sd(lHg, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(lHg, probs=0.975, na.rm = TRUE),
    lower = quantile(lHg, probs=0.025, na.rm = TRUE),
    n = n()
  )
BOSHgXstate.stats.ran
na.omit(BOSHgXstate.stats.ran)
hist(BOS.ran.subset$lHg)


BOSmod1 <- lm(lHg ~ state, data=BOS.ran.subset)
BOSmod.null <- lm(lHg ~ 1, data=BOS.ran.subset)
check_model(BOSmod1, detrend=T)
plot_model(BOSmod1)
BOSpmmod1 <- plot_model(BOSmod1)
BOSpmmod1.coef <- data.frame(state=as.character(BOSpmmod1[[1]]$term), lo=BOSpmmod1[[1]]$conf.low, up=BOSpmmod1[[1]]$conf.high)

BOSpmmod1.coef$nonzero <- ifelse(contains_zero_sign(BOSpmmod1.coef$lo, BOSpmmod1.coef$up)==F, 1, 0)
BOSpmmod1.coef
BOSnonzerosum <- sum(BOSpmmod1.coef$nonzero)
BOSpmmod1.coef[which(BOSpmmod1.coef$nonzero==1),]$state

BOSwAICc <- weight.IC(delta.IC(c(AICc(BOSmod1),AICc(BOSmod.null))))
BOSER <- BOSwAICc[1]/BOSwAICc[2]
BOSER


# resampling loop
iter <- 1000
itdiv <- iter/10

# storage matrix
table(BOS.lm$state)
BOSnlevels <- length(table(BOS.lm$state))
BOSstor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*BOSnlevels)
BOSstatedir.lab <- paste("state",attr(table(BOS.lm$state), "names"),"dir",sep="")
colnames(BOSstor.mat) <- c("ER", "nonzerosum",
                           paste("state",attr(table(BOS.lm$state), "names"),sep=""),BOSstatedir.lab)
head(BOSstor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = BOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = BOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  BOS.ran.subset <- na.omit(BOS.lm[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  BOSmod1 <- lm(lHg ~ state, data=BOS.ran.subset) # class level model
  BOSmod.null <- lm(lHg ~ 1, data=BOS.ran.subset) # null model
  
  # coefficient boundaries
  BOSpmmod1 <- plot_model(BOSmod1)
  BOSpmmod1.coef <- data.frame(state=as.character(BOSpmmod1[[1]]$term), lo=BOSpmmod1[[1]]$conf.low,
                               up=BOSpmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  BOSpmmod1.coef$nonzero <- ifelse(contains_zero_sign(BOSpmmod1.coef$lo, BOSpmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  BOSpmmod1.coef$dir <- ifelse(sign(BOSpmmod1.coef$lo) == -1 & sign(BOSpmmod1.coef$up) == -1, -1, 0)
  BOSpmmod1.coef$dir <- ifelse(sign(BOSpmmod1.coef$lo) == 1 & sign(BOSpmmod1.coef$up) == 1, 1, BOSpmmod1.coef$dir)
  
  # how many variables are non-zero?
  BOSstor.mat[i,"nonzerosum"] <- sum(BOSpmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  BOSnzstates <- BOSpmmod1.coef[which(BOSpmmod1.coef$nonzero==1),]$state
  BOSstor.mat[i, which(colnames(BOSstor.mat) %in% BOSnzstates)] <- 1
  
  # for non-zeros, what is the direction?
  BOSnzstatesdir <- paste(BOSnzstates,"dir",sep="")
  BOSstor.mat[i, which(colnames(BOSstor.mat) %in% BOSnzstatesdir)] <- BOSpmmod1.coef$dir[which(BOSpmmod1.coef$dir != 0)]
  
  # AICc model comparison
  BOSwAICc <- weight.IC(delta.IC(c(AICc(BOSmod1),AICc(BOSmod.null))))
  
  # store evidence ratio
  BOSstor.mat[i,"ER"] <- BOSwAICc[1]/BOSwAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

# ER
median(BOSstor.mat[,"ER"])
quantile(BOSstor.mat[,"ER"], probs=c(0.025,0.975))

head(BOSstor.mat)
state.labs <- c("NSW","NT","QLD","SA","TAS","VIC","WA")
BOSlennzNSW <- length(table(BOSstor.mat[,"stateNSW"]))
BOSlennzNT <- length(table(BOSstor.mat[,"stateNT"]))
BOSlennzQLD <- length(table(BOSstor.mat[,"stateQLD"]))
BOSlennzSA <- length(table(BOSstor.mat[,"stateSA"]))
BOSlennzTAS <- length(table(BOSstor.mat[,"stateTAS"]))
BOSlennzVIC <- length(table(BOSstor.mat[,"stateVIC"]))
BOSlennzWA <- length(table(BOSstor.mat[,"stateWA"]))

BOSnzrslts <- c(ifelse(BOSlennzNSW==1, as.numeric(table(BOSstor.mat[,"stateNSW"])), 0),
                ifelse(BOSlennzNT==1, as.numeric(table(BOSstor.mat[,"stateNT"])), 0),
                ifelse(BOSlennzQLD==1, as.numeric(table(BOSstor.mat[,"stateQLD"])), 0),
                ifelse(BOSlennzSA==1, as.numeric(table(BOSstor.mat[,"stateSA"])), 0),
                ifelse(BOSlennzTAS==1, as.numeric(table(BOSstor.mat[,"stateTAS"])), 0),
                ifelse(BOSlennzVIC==1, as.numeric(table(BOSstor.mat[,"stateVIC"])), 0),
                ifelse(BOSlennzWA==1, as.numeric(table(BOSstor.mat[,"stateWA"])), 0))

BOSdirNSWtab <- table(BOSstor.mat[,"stateNSWdir"])
BOSdirNTtab <- table(BOSstor.mat[,"stateNTdir"])
BOSdirQLDtab <- table(BOSstor.mat[,"stateQLDdir"])
BOSdirSAtab <- table(BOSstor.mat[,"stateSAdir"])
BOSdirTAStab <- table(BOSstor.mat[,"stateTASdir"])
BOSdirVICtab <- table(BOSstor.mat[,"stateVICdir"])
BOSdirWAtab <- table(BOSstor.mat[,"stateWAdir"])

BOSdirlenNSW <- length(BOSdirNSWtab)
BOSdirlenNT <- length(BOSdirNTtab)
BOSdirlenQLD <- length(BOSdirQLDtab)
BOSdirlenSA <- length(BOSdirSAtab)
BOSdirlenTAS <- length(BOSdirTAStab)
BOSdirlenVIC <- length(BOSdirVICtab)
BOSdirlenWA <- length(BOSdirWAtab)

if (BOSdirlenNSW == 0) {
  BOSNSWpos <- 0
  BOSNSWneg <- 0
} else if (BOSdirlenNSW == 1 & as.numeric(attr(BOSdirNSWtab,"names")[1]) == 1) {
  BOSNSWpos <- as.numeric(BOSdirNSWtab)[1]
  BOSNSWneg <- 0
} else if (BOSdirlenNSW == 1 & as.numeric(attr(BOSdirNSWtab,"names")[1]) == -1) {
  BOSNSWpos <- 0
  BOSNSWneg <- as.numeric(BOSdirNSWtab)[1]
} else if (BOSdirlenNSW == 2) {
  BOSNSWneg <- as.numeric(BOSdirNSWtab)[1]
  BOSNSWpos <- as.numeric(BOSdirNSWtab)[2]
} else {
  BOSNSWpos <- 0
  BOSNSWneg <- 0
}

if (BOSdirlenNT == 0) {
  BOSNTpos <- 0
  BOSNTneg <- 0
} else if (BOSdirlenNT == 1 & as.numeric(attr(BOSdirNTtab,"names")[1]) == 1) {
  BOSNTpos <- as.numeric(BOSdirNTtab)[1]
  BOSNTneg <- 0
} else if (BOSdirlenNT == 1 & as.numeric(attr(BOSdirNTtab,"names")[1]) == -1) {
  BOSNTpos <- 0
  BOSNTneg <- as.numeric(BOSdirNTtab)[1]
} else if (BOSdirlenNT == 2) {
  BOSNTneg <- as.numeric(BOSdirNTtab)[1]
  BOSNTpos <- as.numeric(BOSdirNTtab)[2]
} else {
  BOSNTpos <- 0
  BOSNTneg <- 0
}

if (BOSdirlenQLD == 0) {
  BOSQLDpos <- 0
  BOSQLDneg <- 0
} else if (BOSdirlenQLD == 1 & as.numeric(attr(BOSdirQLDtab,"names")[1]) == 1) {
  BOSQLDpos <- as.numeric(BOSdirQLDtab)[1]
  BOSQLDneg <- 0
} else if (BOSdirlenQLD == 1 & as.numeric(attr(BOSdirQLDtab,"names")[1]) == -1) {
  BOSQLDpos <- 0
  BOSQLDneg <- as.numeric(BOSdirQLDtab)[1]
} else if (BOSdirlenQLD == 2) {
  BOSQLDneg <- as.numeric(BOSdirQLDtab)[1]
  BOSQLDpos <- as.numeric(BOSdirQLDtab)[2]
} else {
  BOSQLDpos <- 0
  BOSQLDneg <- 0
}

if (BOSdirlenSA == 0) {
  BOSSApos <- 0
  BOSSAneg <- 0
} else if (BOSdirlenSA == 1 & as.numeric(attr(BOSdirSAtab,"names")[1]) == 1) {
  BOSSApos <- as.numeric(BOSdirSAtab)[1]
  BOSSAneg <- 0
} else if (BOSdirlenSA == 1 & as.numeric(attr(BOSdirSAtab,"names")[1]) == -1) {
  BOSSApos <- 0
  BOSSAneg <- as.numeric(BOSdirSAtab)[1]
} else if (BOSdirlenSA == 2) {
  BOSSAneg <- as.numeric(BOSdirSAtab)[1]
  BOSSApos <- as.numeric(BOSdirSAtab)[2]
} else {
  BOSSApos <- 0
  BOSSAneg <- 0
}

if (BOSdirlenTAS == 0) {
  BOSTASpos <- 0
  BOSTASneg <- 0
} else if (BOSdirlenTAS == 1 & as.numeric(attr(BOSdirTAStab,"names")[1]) == 1) {
  BOSTASpos <- as.numeric(BOSdirTAStab)[1]
  BOSTASneg <- 0
} else if (BOSdirlenTAS == 1 & as.numeric(attr(BOSdirTAStab,"names")[1]) == -1) {
  BOSTASpos <- 0
  BOSTASneg <- as.numeric(BOSdirTAStab)[1]
} else if (BOSdirlenTAS == 2) {
  BOSTASneg <- as.numeric(BOSdirTAStab)[1]
  BOSTASpos <- as.numeric(BOSdirTAStab)[2]
} else {
  BOSTASpos <- 0
  BOSTASneg <- 0
}

if (BOSdirlenVIC == 0) {
  BOSVICpos <- 0
  BOSVICneg <- 0
} else if (BOSdirlenVIC == 1 & as.numeric(attr(BOSdirVICtab,"names")[1]) == 1) {
  BOSVICpos <- as.numeric(BOSdirVICtab)[1]
  BOSVICneg <- 0
} else if (BOSdirlenVIC == 1 & as.numeric(attr(BOSdirVICtab,"names")[1]) == -1) {
  BOSVICpos <- 0
  BOSVICneg <- as.numeric(BOSdirVICtab)[1]
} else if (BOSdirlenVIC == 2) {
  BOSVICneg <- as.numeric(BOSdirVICtab)[1]
  BOSVICpos <- as.numeric(BOSdirVICtab)[2]
} else {
  BOSVICpos <- 0
  BOSVICneg <- 0
}

if (BOSdirlenWA == 0) {
  BOSWApos <- 0
  BOSWAneg <- 0
} else if (BOSdirlenWA == 1 & as.numeric(attr(BOSdirWAtab,"names")[1]) == 1) {
  BOSWApos <- as.numeric(BOSdirWAtab)[1]
  BOSWAneg <- 0
} else if (BOSdirlenWA == 1 & as.numeric(attr(BOSdirWAtab,"names")[1]) == -1) {
  BOSWApos <- 0
  BOSWAneg <- as.numeric(BOSdirWAtab)[1]
} else if (BOSdirlenWA == 2) {
  BOSWAneg <- as.numeric(BOSdirWAtab)[1]
  BOSWApos <- as.numeric(BOSdirWAtab)[2]
} else {
  BOSWApos <- 0
  BOSWAneg <- 0
}


BOSnznegdirrslts <- c(BOSNSWneg, BOSNTneg, BOSQLDneg, BOSSAneg, BOSTASneg, BOSVICneg, BOSWAneg)
BOSnzposdirrslts <- c(BOSNSWpos, BOSNTpos, BOSQLDpos, BOSSApos, BOSTASpos, BOSVICpos, BOSWApos)

BOSresults.out <- data.frame(state=state.labs, nonzero=BOSnzrslts/iter, negdir=BOSnznegdirrslts/iter,
                             posdir=BOSnzposdirrslts/iter)
BOSresults.out

# ER
10^median(log10(BOSstor.mat[,"ER"]))
quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(BOSstor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(BOSstor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(BOSstor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(BOSstor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)


############
# landuse ##
############

## TOS
table(TOS.ran.subset$landuse)

table(TOS.ran.subset$landuse)
TOSHgXlanduse.stats.ran <- TOS.ran.subset %>%
  group_by(landuse) %>%
  summarise(
    mean = mean(lHg, na.rm = TRUE),
    median = median(lHg, na.rm = TRUE), 
    var = var(lHg, na.rm = TRUE),
    sd = sd(lHg, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(lHg, probs=0.975, na.rm = TRUE),
    lower = quantile(lHg, probs=0.025, na.rm = TRUE),
    n = n()
  )
TOSHgXlanduse.stats.ran
na.omit(TOSHgXlanduse.stats.ran)

# resampling loop
# storage matrix
table(TOS.lm$landuse)
TOSnlevels <- length(table(TOS.lm$landuse))
TOSstor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*TOSnlevels)
TOSlandusedir.lab <- paste("landuse",attr(table(TOS.lm$landuse), "names"),"dir",sep="")
colnames(TOSstor.mat) <- c("ER", "nonzerosum",
                           paste("landuse",attr(table(TOS.lm$landuse), "names"),sep=""),TOSlandusedir.lab)
head(TOSstor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = TOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = TOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  TOS.ran.subset <- na.omit(TOS.lm[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  TOSmod1 <- lm(lHg ~ landuse, data=TOS.ran.subset) # class level model
  TOSmod.null <- lm(lHg ~ 1, data=TOS.ran.subset) # null model
  
  # coefficient boundaries
  TOSpmmod1 <- plot_model(TOSmod1)
  TOSpmmod1.coef <- data.frame(landuse=as.character(TOSpmmod1[[1]]$term), lo=TOSpmmod1[[1]]$conf.low,
                               up=TOSpmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  TOSpmmod1.coef$nonzero <- ifelse(contains_zero_sign(TOSpmmod1.coef$lo, TOSpmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  TOSpmmod1.coef$dir <- ifelse(sign(TOSpmmod1.coef$lo) == -1 & sign(TOSpmmod1.coef$up) == -1, -1, 0)
  TOSpmmod1.coef$dir <- ifelse(sign(TOSpmmod1.coef$lo) == 1 & sign(TOSpmmod1.coef$up) == 1, 1, TOSpmmod1.coef$dir)
  
  # how many variables are non-zero?
  TOSstor.mat[i,"nonzerosum"] <- sum(TOSpmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  TOSnzlanduses <- TOSpmmod1.coef[which(TOSpmmod1.coef$nonzero==1),]$landuse
  TOSstor.mat[i, which(colnames(TOSstor.mat) %in% TOSnzlanduses)] <- 1
  
  # for non-zeros, what is the direction?
  TOSnzlandusessdir <- paste(TOSnzlanduses,"dir",sep="")
  TOSstor.mat[i, which(colnames(TOSstor.mat) %in% TOSnzlandusessdir)] <- TOSpmmod1.coef$dir[which(TOSpmmod1.coef$dir != 0)]
  
  # AICc model comparison
  TOSwAICc <- weight.IC(delta.IC(c(AICc(TOSmod1),AICc(TOSmod.null))))
  
  # store evidence ratio
  TOSstor.mat[i,"ER"] <- TOSwAICc[1]/TOSwAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

head(TOSstor.mat)
TOSlanduse.labs <- c("CN","INT","PDA","PIA","PRN","WAT","WA")
TOSlennzCN <- length(table(TOSstor.mat[,"landuseconservation/natural"]))
TOSlennzINT <- length(table(TOSstor.mat[,"landuseintensive"]))
TOSlennzPDA <- length(table(TOSstor.mat[,"landuseproduction-dryland agr"]))
TOSlennzPIA <- length(table(TOSstor.mat[,"landuseproduction-irrigated agr"]))
TOSlennzPRN <- length(table(TOSstor.mat[,"landuseproduction-relatively natural"]))
TOSlennzWAT <- length(table(TOSstor.mat[,"landusewater"]))

TOSnzrslts <- c(ifelse(TOSlennzCN==1, as.numeric(table(TOSstor.mat[,"landuseconservation/natural"])), 0),
                ifelse(TOSlennzINT==1, as.numeric(table(TOSstor.mat[,"landuseintensive"])), 0),
                ifelse(TOSlennzPDA==1, as.numeric(table(TOSstor.mat[,"landuseproduction-dryland agr"])), 0),
                ifelse(TOSlennzPIA==1, as.numeric(table(TOSstor.mat[,"landuseproduction-irrigated agr"])), 0),
                ifelse(TOSlennzPRN==1, as.numeric(table(TOSstor.mat[,"landuseproduction-relatively natural"])), 0),
                ifelse(TOSlennzWAT==1, as.numeric(table(TOSstor.mat[,"landusewater"])), 0))

TOSdirCNtab <- table(TOSstor.mat[,"landuseconservation/naturaldir"])
TOSdirINTtab <- table(TOSstor.mat[,"landuseintensivedir"])
TOSdirPDAtab <- table(TOSstor.mat[,"landuseproduction-dryland agrdir"])
TOSdirPIAtab <- table(TOSstor.mat[,"landuseproduction-irrigated agrdir"])
TOSdirPRNtab <- table(TOSstor.mat[,"landuseproduction-relatively naturaldir"])
TOSdirWATtab <- table(TOSstor.mat[,"landusewater"])

TOSdirlenCN <- length(TOSdirCNtab)
TOSdirlenINT <- length(TOSdirINTtab)
TOSdirlenPDA <- length(TOSdirPDAtab)
TOSdirlenPIA <- length(TOSdirPIAtab)
TOSdirlenPRN <- length(TOSdirPRNtab)
TOSdirlenWAT <- length(TOSdirWATtab)


if (TOSdirlenCN == 0) {
  TOSCNpos <- 0
  TOSCNneg <- 0
} else if (TOSdirlenCN == 1 & as.numeric(attr(TOSdirCNtab,"names")[1]) == 1) {
  TOSCNpos <- as.numeric(TOSdirCNtab)[1]
  TOSCNneg <- 0
} else if (TOSdirlenCN == 1 & as.numeric(attr(TOSdirCNtab,"names")[1]) == -1) {
  TOSCNpos <- 0
  TOSCNneg <- as.numeric(TOSTOSdirCNtab)[1]
} else if (TOSdirlenCN == 2) {
  TOSCNneg <- as.numeric(TOSdirCNtab)[1]
  TOSCNpos <- as.numeric(TOSdirCNtab)[2]
} else {
  TOSCNpos <- 0
  TOSCNneg <- 0
}

if (TOSdirlenINT == 0) {
  TOSINTpos <- 0
  TOSINTneg <- 0
} else if (TOSdirlenINT == 1 & as.numeric(attr(TOSdirINTtab,"names")[1]) == 1) {
  TOSINTpos <- as.numeric(TOSdirINTtab)[1]
  TOSINTneg <- 0
} else if (TOSdirlenINT == 1 & as.numeric(attr(TOSdirINTtab,"names")[1]) == -1) {
  TOSINTpos <- 0
  TOSINTneg <- as.numeric(TOSdirINTtab)[1]
} else if (TOSdirlenINT == 2) {
  TOSINTneg <- as.numeric(TOSdirINTtab)[1]
  TOSINTpos <- as.numeric(TOSdirINTtab)[2]
} else {
  TOSINTpos <- 0
  TOSINTneg <- 0
}

if (TOSdirlenPDA == 0) {
  TOSPDApos <- 0
  TOSPDAneg <- 0
} else if (TOSdirlenPDA == 1 & as.numeric(attr(TOSdirPDAtab,"names")[1]) == 1) {
  TOSPDApos <- as.numeric(TOSdirPDAtab)[1]
  TOSPDAneg <- 0
} else if (TOSdirlenPDA == 1 & as.numeric(attr(TOSdirPDAtab,"names")[1]) == -1) {
  TOSPDApos <- 0
  TOSPDAneg <- as.numeric(TOSdirPDAtab)[1]
} else if (TOSdirlenPDA == 2) {
  TOSPDAneg <- as.numeric(TOSdirPDAtab)[1]
  TOSPDApos <- as.numeric(TOSdirPDAtab)[2]
} else {
  TOSPDApos <- 0
  TOSPDAneg <- 0
}

if (TOSdirlenPIA == 0) {
  TOSPIApos <- 0
  TOSPIAneg <- 0
} else if (TOSdirlenPIA == 1 & as.numeric(attr(TOSdirPIAtab,"names")[1]) == 1) {
  TOSPIApos <- as.numeric(TOSdirPIAtab)[1]
  TOSPIAneg <- 0
} else if (TOSdirlenPIA == 1 & as.numeric(attr(TOSdirPIAtab,"names")[1]) == -1) {
  TOSPIApos <- 0
  TOSPIAneg <- as.numeric(TOSdirPIAtab)[1]
} else if (TOSdirlenPIA == 2) {
  TOSPIAneg <- as.numeric(TOSdirPIAtab)[1]
  TOSPIApos <- as.numeric(TOSdirPIAtab)[2]
} else {
  TOSPIApos <- 0
  TOSPIAneg <- 0
}

if (TOSdirlenPRN == 0) {
  TOSPRNpos <- 0
  TOSPRNneg <- 0
} else if (TOSdirlenPRN == 1 & as.numeric(attr(TOSdirPRNtab,"names")[1]) == 1) {
  TOSPRNpos <- as.numeric(TOSdirPRNtab)[1]
  TOSPRNneg <- 0
} else if (TOSdirlenPRN == 1 & as.numeric(attr(TOSdirPRNtab,"names")[1]) == -1) {
  TOSPRNpos <- 0
  TOSPRNneg <- as.numeric(TOSdirPRNtab)[1]
} else if (TOSdirlenPRN == 2) {
  TOSPRNneg <- as.numeric(TOSdirPRNtab)[1]
  TOSPRNpos <- as.numeric(TOSdirPRNtab)[2]
} else {
  TOSPRNpos <- 0
  TOSPRNneg <- 0
}

if (TOSdirlenWAT == 0) {
  TOSWATpos <- 0
  TOSWATneg <- 0
} else if (TOSdirlenWAT == 1 & as.numeric(attr(TOSdirWATtab,"names")[1]) == 1) {
  TOSWATpos <- as.numeric(TOSdirWATtab)[1]
  TOSWATneg <- 0
} else if (TOSdirlenWAT == 1 & as.numeric(attr(TOSdirWATtab,"names")[1]) == -1) {
  TOSWATpos <- 0
  TOSWATneg <- as.numeric(TOSdirWATtab)[1]
} else if (TOSdirlenWAT == 2) {
  TOSWATneg <- as.numeric(TOSdirWATtab)[1]
  TOSWATpos <- as.numeric(TOSdirWATtab)[2]
} else {
  TOSWATpos <- 0
  TOSWATneg <- 0
}

TOSnznegdirrslts <- c(TOSCNneg, TOSINTneg, TOSPDAneg, TOSPIAneg, TOSPRNneg, TOSWATneg)
TOSnzposdirrslts <- c(TOSCNpos, TOSINTpos, TOSPDApos, TOSPIApos, TOSPRNpos, TOSWATpos)

TOSresults.out <- data.frame(landuse=TOSlandusedir.lab, nonzero=TOSnzrslts/iter, negdir=TOSnznegdirrslts/iter,
                             posdir=TOSnzposdirrslts/iter)
TOSresults.out

# ER
10^median(log10(TOSstor.mat[,"ER"]))
quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(TOSstor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(TOSstor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(TOSstor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(TOSstor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)


## BOS
table(BOS.ran.subset$landuse)

table(BOS.ran.subset$landuse)
BOSHgXlanduse.stats.ran <- BOS.ran.subset %>%
  group_by(landuse) %>%
  summarise(
    mean = mean(lHg, na.rm = TRUE),
    median = median(lHg, na.rm = TRUE), 
    var = var(lHg, na.rm = TRUE),
    sd = sd(lHg, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(lHg, probs=0.975, na.rm = TRUE),
    lower = quantile(lHg, probs=0.025, na.rm = TRUE),
    n = n()
  )
BOSHgXlanduse.stats.ran
na.omit(BOSHgXlanduse.stats.ran)

# resampling loop
# storage matrix
table(BOS.lm$landuse)
BOSnlevels <- length(table(BOS.lm$landuse))
BOSstor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*BOSnlevels)
BOSlandusedir.lab <- paste("landuse",attr(table(BOS.lm$landuse), "names"),"dir",sep="")
colnames(BOSstor.mat) <- c("ER", "nonzerosum",
                           paste("landuse",attr(table(BOS.lm$landuse), "names"),sep=""),BOSlandusedir.lab)
head(BOSstor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = BOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = BOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  BOS.ran.subset <- na.omit(BOS.lm[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  BOSmod1 <- lm(lHg ~ landuse, data=BOS.ran.subset) # class level model
  BOSmod.null <- lm(lHg ~ 1, data=BOS.ran.subset) # null model
  
  # coefficient boundaries
  BOSpmmod1 <- plot_model(BOSmod1)
  BOSpmmod1.coef <- data.frame(landuse=as.character(BOSpmmod1[[1]]$term), lo=BOSpmmod1[[1]]$conf.low,
                               up=BOSpmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  BOSpmmod1.coef$nonzero <- ifelse(contains_zero_sign(BOSpmmod1.coef$lo, BOSpmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  BOSpmmod1.coef$dir <- ifelse(sign(BOSpmmod1.coef$lo) == -1 & sign(BOSpmmod1.coef$up) == -1, -1, 0)
  BOSpmmod1.coef$dir <- ifelse(sign(BOSpmmod1.coef$lo) == 1 & sign(BOSpmmod1.coef$up) == 1, 1, BOSpmmod1.coef$dir)
  
  # how many variables are non-zero?
  BOSstor.mat[i,"nonzerosum"] <- sum(BOSpmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  BOSnzlanduses <- BOSpmmod1.coef[which(BOSpmmod1.coef$nonzero==1),]$landuse
  BOSstor.mat[i, which(colnames(BOSstor.mat) %in% BOSnzlanduses)] <- 1
  
  # for non-zeros, what is the direction?
  BOSnzlandusessdir <- paste(BOSnzlanduses,"dir",sep="")
  BOSstor.mat[i, which(colnames(BOSstor.mat) %in% BOSnzlandusessdir)] <- BOSpmmod1.coef$dir[which(BOSpmmod1.coef$dir != 0)]
  
  # AICc model comparison
  BOSwAICc <- weight.IC(delta.IC(c(AICc(BOSmod1),AICc(BOSmod.null))))
  
  # store evidence ratio
  BOSstor.mat[i,"ER"] <- BOSwAICc[1]/BOSwAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

# ER
median(BOSstor.mat[,"ER"])
quantile(BOSstor.mat[,"ER"], probs=c(0.025,0.975))

head(BOSstor.mat)
BOSlanduse.labs <- c("CN","INT","PDA","PIA","PRN","WAT","WA")
BOSlennzCN <- length(table(BOSstor.mat[,"landuseconservation/natural"]))
BOSlennzINT <- length(table(BOSstor.mat[,"landuseintensive"]))
BOSlennzPDA <- length(table(BOSstor.mat[,"landuseproduction-dryland agr"]))
BOSlennzPIA <- length(table(BOSstor.mat[,"landuseproduction-irrigated agr"]))
BOSlennzPRN <- length(table(BOSstor.mat[,"landuseproduction-relatively natural"]))
BOSlennzWAT <- length(table(BOSstor.mat[,"landusewater"]))

BOSnzrslts <- c(ifelse(BOSlennzCN==1, as.numeric(table(BOSstor.mat[,"landuseconservation/natural"])), 0),
                ifelse(BOSlennzINT==1, as.numeric(table(BOSstor.mat[,"landuseintensive"])), 0),
                ifelse(BOSlennzPDA==1, as.numeric(table(BOSstor.mat[,"landuseproduction-dryland agr"])), 0),
                ifelse(BOSlennzPIA==1, as.numeric(table(BOSstor.mat[,"landuseproduction-irrigated agr"])), 0),
                ifelse(BOSlennzPRN==1, as.numeric(table(BOSstor.mat[,"landuseproduction-relatively natural"])), 0),
                ifelse(BOSlennzWAT==1, as.numeric(table(BOSstor.mat[,"landusewater"])), 0))

BOSdirCNtab <- table(BOSstor.mat[,"landuseconservation/naturaldir"])
BOSdirINTtab <- table(BOSstor.mat[,"landuseintensivedir"])
BOSdirPDAtab <- table(BOSstor.mat[,"landuseproduction-dryland agrdir"])
BOSdirPIAtab <- table(BOSstor.mat[,"landuseproduction-irrigated agrdir"])
BOSdirPRNtab <- table(BOSstor.mat[,"landuseproduction-relatively naturaldir"])
BOSdirWATtab <- table(BOSstor.mat[,"landusewaterdir"])

BOSdirlenCN <- length(BOSdirCNtab)
BOSdirlenINT <- length(BOSdirINTtab)
BOSdirlenPDA <- length(BOSdirPDAtab)
BOSdirlenPIA <- length(BOSdirPIAtab)
BOSdirlenPRN <- length(BOSdirPRNtab)
BOSdirlenWAT <- length(BOSdirWATtab)


if (BOSdirlenCN == 0) {
  BOSCNpos <- 0
  BOSCNneg <- 0
} else if (BOSdirlenCN == 1 & as.numeric(attr(BOSdirCNtab,"names")[1]) == 1) {
  BOSCNpos <- as.numeric(BOSdirCNtab)[1]
  BOSCNneg <- 0
} else if (BOSdirlenCN == 1 & as.numeric(attr(BOSdirCNtab,"names")[1]) == -1) {
  BOSCNpos <- 0
  BOSCNneg <- as.numeric(BOSBOSdirCNtab)[1]
} else if (BOSdirlenCN == 2) {
  BOSCNneg <- as.numeric(BOSdirCNtab)[1]
  BOSCNpos <- as.numeric(BOSdirCNtab)[2]
} else {
  BOSCNpos <- 0
  BOSCNneg <- 0
}

if (BOSdirlenINT == 0) {
  BOSINTpos <- 0
  BOSINTneg <- 0
} else if (BOSdirlenINT == 1 & as.numeric(attr(BOSdirINTtab,"names")[1]) == 1) {
  BOSINTpos <- as.numeric(BOSdirINTtab)[1]
  BOSINTneg <- 0
} else if (BOSdirlenINT == 1 & as.numeric(attr(BOSdirINTtab,"names")[1]) == -1) {
  BOSINTpos <- 0
  BOSINTneg <- as.numeric(BOSdirINTtab)[1]
} else if (BOSdirlenINT == 2) {
  BOSINTneg <- as.numeric(BOSdirINTtab)[1]
  BOSINTpos <- as.numeric(BOSdirINTtab)[2]
} else {
  BOSINTpos <- 0
  BOSINTneg <- 0
}

if (BOSdirlenPDA == 0) {
  BOSPDApos <- 0
  BOSPDAneg <- 0
} else if (BOSdirlenPDA == 1 & as.numeric(attr(BOSdirPDAtab,"names")[1]) == 1) {
  BOSPDApos <- as.numeric(BOSdirPDAtab)[1]
  BOSPDAneg <- 0
} else if (BOSdirlenPDA == 1 & as.numeric(attr(BOSdirPDAtab,"names")[1]) == -1) {
  BOSPDApos <- 0
  BOSPDAneg <- as.numeric(BOSdirPDAtab)[1]
} else if (BOSdirlenPDA == 2) {
  BOSPDAneg <- as.numeric(BOSdirPDAtab)[1]
  BOSPDApos <- as.numeric(BOSdirPDAtab)[2]
} else {
  BOSPDApos <- 0
  BOSPDAneg <- 0
}

if (BOSdirlenPIA == 0) {
  BOSPIApos <- 0
  BOSPIAneg <- 0
} else if (BOSdirlenPIA == 1 & as.numeric(attr(BOSdirPIAtab,"names")[1]) == 1) {
  BOSPIApos <- as.numeric(BOSdirPIAtab)[1]
  BOSPIAneg <- 0
} else if (BOSdirlenPIA == 1 & as.numeric(attr(BOSdirPIAtab,"names")[1]) == -1) {
  BOSPIApos <- 0
  BOSPIAneg <- as.numeric(BOSdirPIAtab)[1]
} else if (BOSdirlenPIA == 2) {
  BOSPIAneg <- as.numeric(BOSdirPIAtab)[1]
  BOSPIApos <- as.numeric(BOSdirPIAtab)[2]
} else {
  BOSPIApos <- 0
  BOSPIAneg <- 0
}

if (BOSdirlenPRN == 0) {
  BOSPRNpos <- 0
  BOSPRNneg <- 0
} else if (BOSdirlenPRN == 1 & as.numeric(attr(BOSdirPRNtab,"names")[1]) == 1) {
  BOSPRNpos <- as.numeric(BOSdirPRNtab)[1]
  BOSPRNneg <- 0
} else if (BOSdirlenPRN == 1 & as.numeric(attr(BOSdirPRNtab,"names")[1]) == -1) {
  BOSPRNpos <- 0
  BOSPRNneg <- as.numeric(BOSdirPRNtab)[1]
} else if (BOSdirlenPRN == 2) {
  BOSPRNneg <- as.numeric(BOSdirPRNtab)[1]
  BOSPRNpos <- as.numeric(BOSdirPRNtab)[2]
} else {
  BOSPRNpos <- 0
  BOSPRNneg <- 0
}

if (BOSdirlenWAT == 0) {
  BOSWATpos <- 0
  BOSWATneg <- 0
} else if (BOSdirlenWAT == 1 & as.numeric(attr(BOSdirWATtab,"names")[1]) == 1) {
  BOSWATpos <- as.numeric(BOSdirWATtab)[1]
  BOSWATneg <- 0
} else if (BOSdirlenWAT == 1 & as.numeric(attr(BOSdirWATtab,"names")[1]) == -1) {
  BOSWATpos <- 0
  BOSWATneg <- as.numeric(BOSdirWATtab)[1]
} else if (BOSdirlenWAT == 2) {
  BOSWATneg <- as.numeric(BOSdirWATtab)[1]
  BOSWATpos <- as.numeric(BOSdirWATtab)[2]
} else {
  BOSWATpos <- 0
  BOSWATneg <- 0
}

BOSnznegdirrslts <- c(BOSCNneg, BOSINTneg, BOSPDAneg, BOSPIAneg, BOSPRNneg, BOSWATneg)
BOSnzposdirrslts <- c(BOSCNpos, BOSINTpos, BOSPDApos, BOSPIApos, BOSPRNpos, BOSWATpos)

BOSresults.out <- data.frame(landuse=BOSlandusedir.lab, nonzero=BOSnzrslts/iter, negdir=BOSnznegdirrslts/iter,
                             posdir=BOSnzposdirrslts/iter)
BOSresults.out

# ER
10^median(log10(BOSstor.mat[,"ER"]))
quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(BOSstor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(BOSstor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(BOSstor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(BOSstor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)



#############
# by biome ##
#############
## TOS
table(TOS.ran.subset$biome)

table(TOS.ran.subset$biome)
TOSHgXbiome.stats.ran <- TOS.ran.subset %>%
  group_by(biome) %>%
  summarise(
    mean = mean(lHg, na.rm = TRUE),
    median = median(lHg, na.rm = TRUE), 
    var = var(lHg, na.rm = TRUE),
    sd = sd(lHg, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(lHg, probs=0.975, na.rm = TRUE),
    lower = quantile(lHg, probs=0.025, na.rm = TRUE),
    n = n()
  )
TOSHgXbiome.stats.ran
na.omit(TOSHgXbiome.stats.ran)

# resampling loop
# storage matrix
table(TOS.lm$biome)
TOSnlevels <- length(table(TOS.lm$biome))
TOSstor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*TOSnlevels)
TOSbiomedir.lab <- paste("biome",attr(table(TOS.lm$biome), "names"),"dir",sep="")
colnames(TOSstor.mat) <- c("ER", "nonzerosum",
                        paste("biome",attr(table(TOS.lm$biome), "names"),sep=""),TOSbiomedir.lab)
head(TOSstor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = TOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = TOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  TOS.ran.subset <- na.omit(TOS.lm[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  TOSmod1 <- lm(lHg ~ biome, data=TOS.ran.subset) # class level model
  TOSmod.null <- lm(lHg ~ 1, data=TOS.ran.subset) # null model
  
  # coefficient boundaries
  TOSpmmod1 <- plot_model(TOSmod1)
  TOSpmmod1.coef <- data.frame(biome=as.character(TOSpmmod1[[1]]$term), lo=TOSpmmod1[[1]]$conf.low,
                            up=TOSpmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  TOSpmmod1.coef$nonzero <- ifelse(contains_zero_sign(TOSpmmod1.coef$lo, TOSpmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  TOSpmmod1.coef$dir <- ifelse(sign(TOSpmmod1.coef$lo) == -1 & sign(TOSpmmod1.coef$up) == -1, -1, 0)
  TOSpmmod1.coef$dir <- ifelse(sign(TOSpmmod1.coef$lo) == 1 & sign(TOSpmmod1.coef$up) == 1, 1, TOSpmmod1.coef$dir)
  
  # how many variables are non-zero?
  TOSstor.mat[i,"nonzerosum"] <- sum(TOSpmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  TOSnzbiomes <- TOSpmmod1.coef[which(TOSpmmod1.coef$nonzero==1),]$biome
  TOSstor.mat[i, which(colnames(TOSstor.mat) %in% TOSnzbiomes)] <- 1
  
  # for non-zeros, what is the direction?
  TOSnzbiomessdir <- paste(TOSnzbiomes,"dir",sep="")
  TOSstor.mat[i, which(colnames(TOSstor.mat) %in% TOSnzbiomessdir)] <- TOSpmmod1.coef$dir[which(TOSpmmod1.coef$dir != 0)]
  
  # AICc model comparison
  TOSwAICc <- weight.IC(delta.IC(c(AICc(TOSmod1),AICc(TOSmod.null))))
  
  # store evidence ratio
  TOSstor.mat[i,"ER"] <- TOSwAICc[1]/TOSwAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

head(TOSstor.mat)
biome.labs <- c("DXS","MFW","TBMF","TGSS","TSGSS","TSMBF")
TOSlennzDXS <- length(table(TOSstor.mat[,"biomeDXS"]))
TOSlennzMFW <- length(table(TOSstor.mat[,"biomeMFW"]))
TOSlennzTBMF <- length(table(TOSstor.mat[,"biomeTBMF"]))
TOSlennzTGSS <- length(table(TOSstor.mat[,"biomeTGSS"]))
TOSlennzTSGSS <- length(table(TOSstor.mat[,"biomeTSGSS"]))
TOSlennzTSMBF <- length(table(TOSstor.mat[,"biomeTSMBF"]))

TOSnzrslts <- c(ifelse(TOSlennzDXS==1, as.numeric(table(TOSstor.mat[,"biomeDXS"])), 0),
             ifelse(TOSlennzMFW==1, as.numeric(table(TOSstor.mat[,"biomeMFW"])), 0),
             ifelse(TOSlennzTBMF==1, as.numeric(table(TOSstor.mat[,"biomeTBMF"])), 0),
             ifelse(TOSlennzTGSS==1, as.numeric(table(TOSstor.mat[,"biomeTGSS"])), 0),
             ifelse(TOSlennzTSGSS==1, as.numeric(table(TOSstor.mat[,"biomeTSGSS"])), 0),
             ifelse(TOSlennzTSMBF==1, as.numeric(table(TOSstor.mat[,"biomeTSMBF"])), 0))

TOSdirDXStab <- table(TOSstor.mat[,"biomeDXSdir"])
TOSdirMFWtab <- table(TOSstor.mat[,"biomeMFWdir"])
TOSdirTBMFtab <- table(TOSstor.mat[,"biomeTBMFdir"])
TOSdirTGSStab <- table(TOSstor.mat[,"biomeTGSSdir"])
TOSdirTSGSStab <- table(TOSstor.mat[,"biomeTSGSSdir"])
TOSdirTSMBFtab <- table(TOSstor.mat[,"biomeTSMBFdir"])

TOSdirlenDXS <- length(TOSdirDXStab)
TOSdirlenMFW <- length(TOSdirMFWtab)
TOSdirlenTBMF <- length(TOSdirTBMFtab)
TOSdirlenTGSS <- length(TOSdirTGSStab)
TOSdirlenTSGSS <- length(TOSdirTSGSStab)
TOSdirlenTSMBF <- length(TOSdirTSMBFtab)


if (TOSdirlenDXS == 0) {
  TOSDXSpos <- 0
  TOSDXSneg <- 0
} else if (TOSdirlenDXS == 1 & as.numeric(attr(TOSdirDXStab,"names")[1]) == 1) {
  TOSDXSpos <- as.numeric(TOSdirDXStab)[1]
  TOSDXSneg <- 0
} else if (TOSdirlenDXS == 1 & as.numeric(attr(TOSdirDXStab,"names")[1]) == -1) {
  TOSDXSpos <- 0
  TOSDXSneg <- as.numeric(TOSdirDXStab)[1]
} else if (TOSdirlenDXS == 2) {
  TOSDXSneg <- as.numeric(TOSdirDXStab)[1]
  TOSDXSpos <- as.numeric(TOSdirDXStab)[2]
} else {
  TOSDXSpos <- 0
  TOSDXSneg <- 0
}

if (TOSdirlenMFW == 0) {
  TOSMFWpos <- 0
  TOSMFWneg <- 0
} else if (TOSdirlenMFW == 1 & as.numeric(attr(TOSdirMFWtab,"names")[1]) == 1) {
  TOSMFWpos <- as.numeric(TOSdirMFWtab)[1]
  TOSMFWneg <- 0
} else if (TOSdirlenMFW == 1 & as.numeric(attr(TOSdirMFWtab,"names")[1]) == -1) {
  TOSMFWpos <- 0
  TOSMFWneg <- as.numeric(TOSdirMFWtab)[1]
} else if (TOSdirlenMFW == 2) {
  TOSMFWneg <- as.numeric(TOSdirMFWtab)[1]
  TOSMFWpos <- as.numeric(TOSdirMFWtab)[2]
} else {
  TOSMFWpos <- 0
  TOSMFWneg <- 0
}

if (TOSdirlenTBMF == 0) {
  TOSTBMFpos <- 0
  TOSTBMFneg <- 0
} else if (TOSdirlenTBMF == 1 & as.numeric(attr(TOSdirTBMFtab,"names")[1]) == 1) {
  TOSTBMFpos <- as.numeric(TOSdirTBMFtab)[1]
  TOSTBMFneg <- 0
} else if (TOSdirlenTBMF == 1 & as.numeric(attr(TOSdirTBMFtab,"names")[1]) == -1) {
  TOSTBMFpos <- 0
  TOSTBMFneg <- as.numeric(TOSdirTBMFtab)[1]
} else if (TOSdirlenTBMF == 2) {
  TOSTBMFneg <- as.numeric(TOSdirTBMFtab)[1]
  TOSTBMFpos <- as.numeric(TOSdirTBMFtab)[2]
} else {
  TOSTBMFpos <- 0
  TOSTBMFneg <- 0
}

if (TOSdirlenTGSS == 0) {
  TOSTGSSpos <- 0
  TOSTGSSneg <- 0
} else if (TOSdirlenTGSS == 1 & as.numeric(attr(TOSdirTGSStab,"names")[1]) == 1) {
  TOSTGSSpos <- as.numeric(TOSdirTGSStab)[1]
  TOSTGSSneg <- 0
} else if (TOSdirlenTGSS == 1 & as.numeric(attr(TOSdirTGSStab,"names")[1]) == -1) {
  TOSTGSSpos <- 0
  TOSTGSSneg <- as.numeric(TOSdirTGSStab)[1]
} else if (TOSdirlenTGSS == 2) {
  TOSTGSSneg <- as.numeric(TOSdirTGSStab)[1]
  TOSTGSSpos <- as.numeric(TOSdirTGSStab)[2]
} else {
  TOSTGSSpos <- 0
  TOSTGSSneg <- 0
}

if (TOSdirlenTSGSS == 0) {
  TOSTSGSSpos <- 0
  TOSTSGSSneg <- 0
} else if (TOSdirlenTSGSS == 1 & as.numeric(attr(TOSdirTSGSStab,"names")[1]) == 1) {
  TOSTSGSSpos <- as.numeric(TOSdirTSGSStab)[1]
  TOSTSGSSneg <- 0
} else if (TOSdirlenTSGSS == 1 & as.numeric(attr(TOSdirTSGSStab,"names")[1]) == -1) {
  TOSTSGSSpos <- 0
  TOSTSGSSneg <- as.numeric(TOSdirTSGSStab)[1]
} else if (TOSdirlenTSGSS == 2) {
  TOSTSGSSneg <- as.numeric(TOSdirTSGSStab)[1]
  TOSTSGSSpos <- as.numeric(TOSdirTSGSStab)[2]
} else {
  TOSTSGSSpos <- 0
  TOSTSGSSneg <- 0
}

if (TOSdirlenTSMBF == 0) {
  TOSTSMBFpos <- 0
  TOSTSMBFneg <- 0
} else if (TOSdirlenTSMBF == 1 & as.numeric(attr(TOSdirTSMBFtab,"names")[1]) == 1) {
  TOSTSMBFpos <- as.numeric(TOSdirTSMBFtab)[1]
  TOSTSMBFneg <- 0
} else if (TOSdirlenTSMBF == 1 & as.numeric(attr(TOSdirTSMBFtab,"names")[1]) == -1) {
  TOSTSMBFpos <- 0
  TOSTSMBFneg <- as.numeric(TOSdirTSMBFtab)[1]
} else if (TOSdirlenTSMBF == 2) {
  TOSTSMBFneg <- as.numeric(TOSdirTSMBFtab)[1]
  TOSTSMBFpos <- as.numeric(TOSdirTSMBFtab)[2]
} else {
  TOSTSMBFpos <- 0
  TOSTSMBFneg <- 0
}

TOSnznegdirrslts <- c(TOSDXSneg, TOSMFWneg, TOSTBMFneg, TOSTGSSneg, TOSTSGSSneg, TOSTSMBFneg)
TOSnzposdirrslts <- c(TOSDXSpos, TOSMFWpos, TOSTBMFpos, TOSTGSSpos, TOSTSGSSpos, TOSTSMBFpos)

TOSresults.out <- data.frame(biome=TOSbiomedir.lab, nonzero=TOSnzrslts/iter, negdir=TOSnznegdirrslts/iter,
                          posdir=TOSnzposdirrslts/iter)
TOSresults.out

# ER
10^median(log10(TOSstor.mat[,"ER"]))
quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(TOSstor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(TOSstor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(TOSstor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(TOSstor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)


## BOS
table(BOS.ran.subset$biome)

table(BOS.ran.subset$biome)
BOSHgXbiome.stats.ran <- BOS.ran.subset %>%
  group_by(biome) %>%
  summarise(
    mean = mean(lHg, na.rm = TRUE),
    median = median(lHg, na.rm = TRUE), 
    var = var(lHg, na.rm = TRUE),
    sd = sd(lHg, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(lHg, probs=0.975, na.rm = TRUE),
    lower = quantile(lHg, probs=0.025, na.rm = TRUE),
    n = n()
  )
BOSHgXbiome.stats.ran
na.omit(BOSHgXbiome.stats.ran)

# resampling loop
# storage matrix
table(BOS.lm$biome)
BOSnlevels <- length(table(BOS.lm$biome))
BOSstor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*BOSnlevels)
BOSbiomedir.lab <- paste("biome",attr(table(BOS.lm$biome), "names"),"dir",sep="")
colnames(BOSstor.mat) <- c("ER", "nonzerosum",
                           paste("biome",attr(table(BOS.lm$biome), "names"),sep=""),BOSbiomedir.lab)
head(BOSstor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = BOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = BOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  BOS.ran.subset <- na.omit(BOS.lm[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  BOSmod1 <- lm(lHg ~ biome, data=BOS.ran.subset) # class level model
  BOSmod.null <- lm(lHg ~ 1, data=BOS.ran.subset) # null model
  
  # coefficient boundaries
  BOSpmmod1 <- plot_model(BOSmod1)
  BOSpmmod1.coef <- data.frame(biome=as.character(BOSpmmod1[[1]]$term), lo=BOSpmmod1[[1]]$conf.low,
                               up=BOSpmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  BOSpmmod1.coef$nonzero <- ifelse(contains_zero_sign(BOSpmmod1.coef$lo, BOSpmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  BOSpmmod1.coef$dir <- ifelse(sign(BOSpmmod1.coef$lo) == -1 & sign(BOSpmmod1.coef$up) == -1, -1, 0)
  BOSpmmod1.coef$dir <- ifelse(sign(BOSpmmod1.coef$lo) == 1 & sign(BOSpmmod1.coef$up) == 1, 1, BOSpmmod1.coef$dir)
  
  # how many variables are non-zero?
  BOSstor.mat[i,"nonzerosum"] <- sum(BOSpmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  BOSnzbiomes <- BOSpmmod1.coef[which(BOSpmmod1.coef$nonzero==1),]$biome
  BOSstor.mat[i, which(colnames(BOSstor.mat) %in% BOSnzbiomes)] <- 1
  
  # for non-zeros, what is the direction?
  BOSnzbiomessdir <- paste(BOSnzbiomes,"dir",sep="")
  BOSstor.mat[i, which(colnames(BOSstor.mat) %in% BOSnzbiomessdir)] <- BOSpmmod1.coef$dir[which(BOSpmmod1.coef$dir != 0)]
  
  # AICc model comparison
  BOSwAICc <- weight.IC(delta.IC(c(AICc(BOSmod1),AICc(BOSmod.null))))
  
  # store evidence ratio
  BOSstor.mat[i,"ER"] <- BOSwAICc[1]/BOSwAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

head(BOSstor.mat)
biome.labs <- c("DXS","MFW","TBMF","TGSS","TSGSS","TSMBF")
BOSlennzDXS <- length(table(BOSstor.mat[,"biomeDXS"]))
BOSlennzMFW <- length(table(BOSstor.mat[,"biomeMFW"]))
BOSlennzTBMF <- length(table(BOSstor.mat[,"biomeTBMF"]))
BOSlennzTGSS <- length(table(BOSstor.mat[,"biomeTGSS"]))
BOSlennzTSGSS <- length(table(BOSstor.mat[,"biomeTSGSS"]))
BOSlennzTSMBF <- length(table(BOSstor.mat[,"biomeTSMBF"]))

BOSnzrslts <- c(ifelse(BOSlennzDXS==1, as.numeric(table(BOSstor.mat[,"biomeDXS"])), 0),
                ifelse(BOSlennzMFW==1, as.numeric(table(BOSstor.mat[,"biomeMFW"])), 0),
                ifelse(BOSlennzTBMF==1, as.numeric(table(BOSstor.mat[,"biomeTBMF"])), 0),
                ifelse(BOSlennzTGSS==1, as.numeric(table(BOSstor.mat[,"biomeTGSS"])), 0),
                ifelse(BOSlennzTSGSS==1, as.numeric(table(BOSstor.mat[,"biomeTSGSS"])), 0),
                ifelse(BOSlennzTSMBF==1, as.numeric(table(BOSstor.mat[,"biomeTSMBF"])), 0))

BOSdirDXStab <- table(BOSstor.mat[,"biomeDXSdir"])
BOSdirMFWtab <- table(BOSstor.mat[,"biomeMFWdir"])
BOSdirTBMFtab <- table(BOSstor.mat[,"biomeTBMFdir"])
BOSdirTGSStab <- table(BOSstor.mat[,"biomeTGSSdir"])
BOSdirTSGSStab <- table(BOSstor.mat[,"biomeTSGSSdir"])
BOSdirTSMBFtab <- table(BOSstor.mat[,"biomeTSMBFdir"])

BOSdirlenDXS <- length(BOSdirDXStab)
BOSdirlenMFW <- length(BOSdirMFWtab)
BOSdirlenTBMF <- length(BOSdirTBMFtab)
BOSdirlenTGSS <- length(BOSdirTGSStab)
BOSdirlenTSGSS <- length(BOSdirTSGSStab)
BOSdirlenTSMBF <- length(BOSdirTSMBFtab)


if (BOSdirlenDXS == 0) {
  BOSDXSpos <- 0
  BOSDXSneg <- 0
} else if (BOSdirlenDXS == 1 & as.numeric(attr(BOSdirDXStab,"names")[1]) == 1) {
  BOSDXSpos <- as.numeric(BOSdirDXStab)[1]
  BOSDXSneg <- 0
} else if (BOSdirlenDXS == 1 & as.numeric(attr(BOSdirDXStab,"names")[1]) == -1) {
  BOSDXSpos <- 0
  BOSDXSneg <- as.numeric(BOSdirDXStab)[1]
} else if (BOSdirlenDXS == 2) {
  BOSDXSneg <- as.numeric(BOSdirDXStab)[1]
  BOSDXSpos <- as.numeric(BOSdirDXStab)[2]
} else {
  BOSDXSpos <- 0
  BOSDXSneg <- 0
}

if (BOSdirlenMFW == 0) {
  BOSMFWpos <- 0
  BOSMFWneg <- 0
} else if (BOSdirlenMFW == 1 & as.numeric(attr(BOSdirMFWtab,"names")[1]) == 1) {
  BOSMFWpos <- as.numeric(BOSdirMFWtab)[1]
  BOSMFWneg <- 0
} else if (BOSdirlenMFW == 1 & as.numeric(attr(BOSdirMFWtab,"names")[1]) == -1) {
  BOSMFWpos <- 0
  BOSMFWneg <- as.numeric(BOSdirMFWtab)[1]
} else if (BOSdirlenMFW == 2) {
  BOSMFWneg <- as.numeric(BOSdirMFWtab)[1]
  BOSMFWpos <- as.numeric(BOSdirMFWtab)[2]
} else {
  BOSMFWpos <- 0
  BOSMFWneg <- 0
}

if (BOSdirlenTBMF == 0) {
  BOSTBMFpos <- 0
  BOSTBMFneg <- 0
} else if (BOSdirlenTBMF == 1 & as.numeric(attr(BOSdirTBMFtab,"names")[1]) == 1) {
  BOSTBMFpos <- as.numeric(BOSdirTBMFtab)[1]
  BOSTBMFneg <- 0
} else if (BOSdirlenTBMF == 1 & as.numeric(attr(BOSdirTBMFtab,"names")[1]) == -1) {
  BOSTBMFpos <- 0
  BOSTBMFneg <- as.numeric(BOSdirTBMFtab)[1]
} else if (BOSdirlenTBMF == 2) {
  BOSTBMFneg <- as.numeric(BOSdirTBMFtab)[1]
  BOSTBMFpos <- as.numeric(BOSdirTBMFtab)[2]
} else {
  BOSTBMFpos <- 0
  BOSTBMFneg <- 0
}

if (BOSdirlenTGSS == 0) {
  BOSTGSSpos <- 0
  BOSTGSSneg <- 0
} else if (BOSdirlenTGSS == 1 & as.numeric(attr(BOSdirTGSStab,"names")[1]) == 1) {
  BOSTGSSpos <- as.numeric(BOSdirTGSStab)[1]
  BOSTGSSneg <- 0
} else if (BOSdirlenTGSS == 1 & as.numeric(attr(BOSdirTGSStab,"names")[1]) == -1) {
  BOSTGSSpos <- 0
  BOSTGSSneg <- as.numeric(BOSdirTGSStab)[1]
} else if (BOSdirlenTGSS == 2) {
  BOSTGSSneg <- as.numeric(BOSdirTGSStab)[1]
  BOSTGSSpos <- as.numeric(BOSdirTGSStab)[2]
} else {
  BOSTGSSpos <- 0
  BOSTGSSneg <- 0
}

if (BOSdirlenTSGSS == 0) {
  BOSTSGSSpos <- 0
  BOSTSGSSneg <- 0
} else if (BOSdirlenTSGSS == 1 & as.numeric(attr(BOSdirTSGSStab,"names")[1]) == 1) {
  BOSTSGSSpos <- as.numeric(BOSdirTSGSStab)[1]
  BOSTSGSSneg <- 0
} else if (BOSdirlenTSGSS == 1 & as.numeric(attr(BOSdirTSGSStab,"names")[1]) == -1) {
  BOSTSGSSpos <- 0
  BOSTSGSSneg <- as.numeric(BOSdirTSGSStab)[1]
} else if (BOSdirlenTSGSS == 2) {
  BOSTSGSSneg <- as.numeric(BOSdirTSGSStab)[1]
  BOSTSGSSpos <- as.numeric(BOSdirTSGSStab)[2]
} else {
  BOSTSGSSpos <- 0
  BOSTSGSSneg <- 0
}

if (BOSdirlenTSMBF == 0) {
  BOSTSMBFpos <- 0
  BOSTSMBFneg <- 0
} else if (BOSdirlenTSMBF == 1 & as.numeric(attr(BOSdirTSMBFtab,"names")[1]) == 1) {
  BOSTSMBFpos <- as.numeric(BOSdirTSMBFtab)[1]
  BOSTSMBFneg <- 0
} else if (BOSdirlenTSMBF == 1 & as.numeric(attr(BOSdirTSMBFtab,"names")[1]) == -1) {
  BOSTSMBFpos <- 0
  BOSTSMBFneg <- as.numeric(BOSdirTSMBFtab)[1]
} else if (BOSdirlenTSMBF == 2) {
  BOSTSMBFneg <- as.numeric(BOSdirTSMBFtab)[1]
  BOSTSMBFpos <- as.numeric(BOSdirTSMBFtab)[2]
} else {
  BOSTSMBFpos <- 0
  BOSTSMBFneg <- 0
}

BOSnznegdirrslts <- c(BOSDXSneg, BOSMFWneg, BOSTBMFneg, BOSTGSSneg, BOSTSGSSneg, BOSTSMBFneg)
BOSnzposdirrslts <- c(BOSDXSpos, BOSMFWpos, BOSTBMFpos, BOSTGSSpos, BOSTSGSSpos, BOSTSMBFpos)

BOSresults.out <- data.frame(biome=BOSbiomedir.lab, nonzero=BOSnzrslts/iter, negdir=BOSnznegdirrslts/iter,
                             posdir=BOSnzposdirrslts/iter)
BOSresults.out

# ER
10^median(log10(BOSstor.mat[,"ER"]))
quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(BOSstor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(BOSstor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(BOSstor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(BOSstor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)


##############
# lithology ##
##############

## TOS
table(TOS.ran.subset$geol)
colnames(TOS.ran.subset)

table(TOS.ran.subset$geol)
TOSHgXgeol.stats.ran <- TOS.ran.subset %>%
  group_by(geol) %>%
  summarise(
    mean = mean(lHg, na.rm = TRUE),
    median = median(lHg, na.rm = TRUE), 
    var = var(lHg, na.rm = TRUE),
    sd = sd(lHg, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(lHg, probs=0.975, na.rm = TRUE),
    lower = quantile(lHg, probs=0.025, na.rm = TRUE),
    n = n()
  )
TOSHgXgeol.stats.ran
na.omit(TOSHgXgeol.stats.ran)

# resampling loop
# storage matrix
table(TOS.lm$geol)
TOSnlevels <- length(table(TOS.lm$geol))
TOSstor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*TOSnlevels)
TOSgeoldir.lab <- paste("geol",attr(table(TOS.lm$geol), "names"),"dir",sep="")
colnames(TOSstor.mat) <- c("ER", "nonzerosum",
                        paste("geol",attr(table(TOS.lm$geol), "names"),sep=""),TOSgeoldir.lab)
head(TOSstor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = TOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = TOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  TOS.ran.subset <- na.omit(TOS.lm[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  TOSmod1 <- lm(lHg ~ geol, data=TOS.ran.subset) # class level model
  TOSmod.null <- lm(lHg ~ 1, data=TOS.ran.subset) # null model
  
  # coefficient boundaries
  TOSpmmod1 <- plot_model(TOSmod1)
  TOSpmmod1.coef <- data.frame(geol=as.character(TOSpmmod1[[1]]$term), lo=TOSpmmod1[[1]]$conf.low,
                            up=TOSpmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  TOSpmmod1.coef$nonzero <- ifelse(contains_zero_sign(TOSpmmod1.coef$lo, TOSpmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  TOSpmmod1.coef$dir <- ifelse(sign(TOSpmmod1.coef$lo) == -1 & sign(TOSpmmod1.coef$up) == -1, -1, 0)
  TOSpmmod1.coef$dir <- ifelse(sign(TOSpmmod1.coef$lo) == 1 & sign(TOSpmmod1.coef$up) == 1, 1, TOSpmmod1.coef$dir)
  
  # how many variables are non-zero?
  TOSstor.mat[i,"nonzerosum"] <- sum(TOSpmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  TOSnzgeol <- TOSpmmod1.coef[which(TOSpmmod1.coef$nonzero==1),]$geol
  TOSstor.mat[i, which(colnames(TOSstor.mat) %in% TOSnzgeol)] <- 1
  
  # for non-zeros, what is the direction?
  TOSnzgeoldir <- paste(TOSnzgeol,"dir",sep="")
  TOSstor.mat[i, which(colnames(TOSstor.mat) %in% TOSnzgeoldir)] <- TOSpmmod1.coef$dir[which(TOSpmmod1.coef$dir != 0)]
  
  # AICc model comparison
  TOSwAICc <- weight.IC(delta.IC(c(AICc(TOSmod1),AICc(TOSmod.null))))
  
  # store evidence ratio
  TOSstor.mat[i,"ER"] <- TOSwAICc[1]/TOSwAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

head(TOSstor.mat)
geol.labs <- c("CAR","FEL","INT","MAF","MET","OTH","SIL")
TOSlennzCAR <- length(table(TOSstor.mat[,"geolCAR"]))
TOSlennzFEL <- length(table(TOSstor.mat[,"geolFEL"]))
TOSlennzINT <- length(table(TOSstor.mat[,"geolINT"]))
TOSlennzMAF <- length(table(TOSstor.mat[,"geolMAF"]))
TOSlennzMET <- length(table(TOSstor.mat[,"geolMET"]))
TOSlennzOTH <- length(table(TOSstor.mat[,"geolOTH"]))
TOSlennzSIL <- length(table(TOSstor.mat[,"geolSIL"]))

TOSnzrslts <- c(ifelse(TOSlennzCAR==1, as.numeric(table(TOSstor.mat[,"geolCAR"])), 0),
             ifelse(TOSlennzFEL==1, as.numeric(table(TOSstor.mat[,"geolFEL"])), 0),
             ifelse(TOSlennzINT==1, as.numeric(table(TOSstor.mat[,"geolINT"])), 0),
             ifelse(TOSlennzMAF==1, as.numeric(table(TOSstor.mat[,"geolMAF"])), 0),
             ifelse(TOSlennzMET==1, as.numeric(table(TOSstor.mat[,"geolMET"])), 0),
             ifelse(TOSlennzOTH==1, as.numeric(table(TOSstor.mat[,"geolOTH"])), 0),
             ifelse(TOSlennzSIL==1, as.numeric(table(TOSstor.mat[,"geolSIL"])), 0))

TOSdirCARtab <- table(TOSstor.mat[,"geolCARdir"])
TOSdirFELtab <- table(TOSstor.mat[,"geolFELdir"])
TOSdirINTtab <- table(TOSstor.mat[,"geolINTdir"])
TOSdirMAFtab <- table(TOSstor.mat[,"geolMAFdir"])
TOSdirMETtab <- table(TOSstor.mat[,"geolMETdir"])
TOSdirOTHtab <- table(TOSstor.mat[,"geolOTHdir"])
TOSdirSILtab <- table(TOSstor.mat[,"geolSILdir"])

TOSdirlenCAR <- length(TOSdirCARtab)
TOSdirlenFEL <- length(TOSdirFELtab)
TOSdirlenINT <- length(TOSdirINTtab)
TOSdirlenMAF <- length(TOSdirMAFtab)
TOSdirlenMET <- length(TOSdirMETtab)
TOSdirlenOTH <- length(TOSdirOTHtab)
TOSdirlenSIL <- length(TOSdirSILtab)


if (TOSdirlenCAR == 0) {
  TOSCARpos <- 0
  TOSCARneg <- 0
} else if (TOSdirlenCAR == 1 & as.numeric(attr(TOSdirCARtab,"names")[1]) == 1) {
  TOSCARpos <- as.numeric(TOSdirCARtab)[1]
  TOSCARneg <- 0
} else if (TOSdirlenCAR == 1 & as.numeric(attr(TOSdirCARtab,"names")[1]) == -1) {
  TOSCARpos <- 0
  TOSCARneg <- as.numeric(TOSdirCARtab)[1]
} else if (dirlenCAR == 2) {
  TOSCARneg <- as.numeric(TOSdirCARtab)[1]
  TOSCARpos <- as.numeric(TOSdirCARtab)[2]
} else {
  TOSCARpos <- 0
  TOSCARneg <- 0
}

if (TOSdirlenFEL == 0) {
  TOSFELpos <- 0
  TOSFELneg <- 0
} else if (TOSdirlenFEL == 1 & as.numeric(attr(TOSdirFELtab,"names")[1]) == 1) {
  TOSFELpos <- as.numeric(TOSdirFELtab)[1]
  TOSFELneg <- 0
} else if (TOSdirlenFEL == 1 & as.numeric(attr(TOSdirFELtab,"names")[1]) == -1) {
  TOSFELpos <- 0
  TOSFELneg <- as.numeric(TOSdirFELtab)[1]
} else if (TOSdirlenFEL == 2) {
  TOSFELneg <- as.numeric(TOSdirFELtab)[1]
  TOSFELpos <- as.numeric(TOSdirFELtab)[2]
} else {
  TOSFELpos <- 0
  TOSFELneg <- 0
}

if (TOSdirlenINT == 0) {
  TOSINTpos <- 0
  TOSINTneg <- 0
} else if (TOSdirlenINT == 1 & as.numeric(attr(TOSdirINTtab,"names")[1]) == 1) {
  TOSINTpos <- as.numeric(TOSdirINTtab)[1]
  TOSINTneg <- 0
} else if (TOSdirlenINT == 1 & as.numeric(attr(TOSdirINTtab,"names")[1]) == -1) {
  TOSINTpos <- 0
  TOSINTneg <- as.numeric(TOSdirINTtab)[1]
} else if (dirlenINT == 2) {
  TOSINTneg <- as.numeric(TOSdirINTtab)[1]
  TOSINTpos <- as.numeric(TOSdirINTtab)[2]
} else {
  TOSINTpos <- 0
  TOSINTneg <- 0
}

if (TOSdirlenMAF == 0) {
  TOSMAFpos <- 0
  TOSMAFneg <- 0
} else if (TOSdirlenMAF == 1 & as.numeric(attr(TOSdirMAFtab,"names")[1]) == 1) {
  TOSMAFpos <- as.numeric(TOSdirMAFtab)[1]
  TOSMAFneg <- 0
} else if (TOSdirlenMAF == 1 & as.numeric(attr(TOSdirMAFtab,"names")[1]) == -1) {
  TOSMAFpos <- 0
  TOSMAFneg <- as.numeric(TOSdirMAFtab)[1]
} else if (TOSdirlenMAF == 2) {
  TOSMAFneg <- as.numeric(TOSdirMAFtab)[1]
  TOSMAFpos <- as.numeric(TOSdirMAFtab)[2]
} else {
  TOSMAFpos <- 0
  TOSMAFneg <- 0
}

if (TOSdirlenMET == 0) {
  TOSMETpos <- 0
  TOSMETneg <- 0
} else if (TOSdirlenMET == 1 & as.numeric(attr(TOSdirMETtab,"names")[1]) == 1) {
  TOSMETpos <- as.numeric(TOSdirMETtab)[1]
  TOSMETneg <- 0
} else if (TOSdirlenMET == 1 & as.numeric(attr(TOSdirMETtab,"names")[1]) == -1) {
  TOSMETpos <- 0
  TOSMETneg <- as.numeric(TOSdirMETtab)[1]
} else if (TOSdirlenMET == 2) {
  TOSMETneg <- as.numeric(TOSdirMETtab)[1]
  TOSMETpos <- as.numeric(TOSdirMETtab)[2]
} else {
  TOSMETpos <- 0
  TOSMETneg <- 0
}

if (TOSdirlenOTH == 0) {
  TOSOTHpos <- 0
  TOSOTHneg <- 0
} else if (TOSdirlenOTH == 1 & as.numeric(attr(TOSdirOTHtab,"names")[1]) == 1) {
  TOSOTHpos <- as.numeric(TOSdirOTHtab)[1]
  TOSOTHneg <- 0
} else if (TOSdirlenOTH == 1 & as.numeric(attr(TOSdirOTHtab,"names")[1]) == -1) {
  TOSOTHpos <- 0
  TOSOTHneg <- as.numeric(TOSdirOTHtab)[1]
} else if (TOSdirlenOTH == 2) {
  TOSOTHneg <- as.numeric(TOSdirOTHtab)[1]
  TOSOTHpos <- as.numeric(TOSdirOTHtab)[2]
} else {
  TOSOTHpos <- 0
  TOSOTHneg <- 0
}

if (TOSdirlenSIL == 0) {
  TOSSILpos <- 0
  TOSSILneg <- 0
} else if (TOSdirlenSIL == 1 & as.numeric(attr(TOSdirSILtab,"names")[1]) == 1) {
  TOSSILpos <- as.numeric(TOSdirSILtab)[1]
  TOSSILneg <- 0
} else if (TOSdirlenSIL == 1 & as.numeric(attr(TOSdirSILtab,"names")[1]) == -1) {
  TOSSILpos <- 0
  TOSSILneg <- as.numeric(TOSdirSILtab)[1]
} else if (TOSdirlenSIL == 2) {
  TOSSILneg <- as.numeric(TOSdirSILtab)[1]
  TOSSILpos <- as.numeric(TOSdirSILtab)[2]
} else {
  TOSSILpos <- 0
  TOSSILneg <- 0
}

TOSnznegdirrslts <- c(TOSCARneg, TOSFELneg, TOSINTneg, TOSMAFneg, TOSMETneg, TOSOTHneg, TOSSILneg)
TOSnzposdirrslts <- c(TOSCARpos, TOSFELpos, TOSINTpos, TOSMAFpos, TOSMETpos, TOSOTHpos, TOSSILpos)

TOSresults.out <- data.frame(geol=TOSgeoldir.lab, nonzero=TOSnzrslts/iter, negdir=TOSnznegdirrslts/iter,
                          posdir=TOSnzposdirrslts/iter)
TOSresults.out

# ER
10^median(log10(TOSstor.mat[,"ER"]))
quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(TOSstor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(TOSstor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(TOSstor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(TOSstor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)


## BOS
table(BOS.ran.subset$geol)
colnames(BOS.ran.subset)

table(BOS.ran.subset$geol)
BOSHgXgeol.stats.ran <- BOS.ran.subset %>%
  group_by(geol) %>%
  summarise(
    mean = mean(lHg, na.rm = TRUE),
    median = median(lHg, na.rm = TRUE), 
    var = var(lHg, na.rm = TRUE),
    sd = sd(lHg, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(lHg, probs=0.975, na.rm = TRUE),
    lower = quantile(lHg, probs=0.025, na.rm = TRUE),
    n = n()
  )
BOSHgXgeol.stats.ran
na.omit(BOSHgXgeol.stats.ran)

# resampling loop
# storage matrix
table(BOS.lm$geol)
BOSnlevels <- length(table(BOS.lm$geol))
BOSstor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*BOSnlevels)
BOSgeoldir.lab <- paste("geol",attr(table(BOS.lm$geol), "names"),"dir",sep="")
colnames(BOSstor.mat) <- c("ER", "nonzerosum",
                           paste("geol",attr(table(BOS.lm$geol), "names"),sep=""),BOSgeoldir.lab)
head(BOSstor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = BOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = BOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  BOS.ran.subset <- na.omit(BOS.lm[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  BOSmod1 <- lm(lHg ~ geol, data=BOS.ran.subset) # class level model
  BOSmod.null <- lm(lHg ~ 1, data=BOS.ran.subset) # null model
  
  # coefficient boundaries
  BOSpmmod1 <- plot_model(BOSmod1)
  BOSpmmod1.coef <- data.frame(geol=as.character(BOSpmmod1[[1]]$term), lo=BOSpmmod1[[1]]$conf.low,
                               up=BOSpmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  BOSpmmod1.coef$nonzero <- ifelse(contains_zero_sign(BOSpmmod1.coef$lo, BOSpmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  BOSpmmod1.coef$dir <- ifelse(sign(BOSpmmod1.coef$lo) == -1 & sign(BOSpmmod1.coef$up) == -1, -1, 0)
  BOSpmmod1.coef$dir <- ifelse(sign(BOSpmmod1.coef$lo) == 1 & sign(BOSpmmod1.coef$up) == 1, 1, BOSpmmod1.coef$dir)
  
  # how many variables are non-zero?
  BOSstor.mat[i,"nonzerosum"] <- sum(BOSpmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  BOSnzgeol <- BOSpmmod1.coef[which(BOSpmmod1.coef$nonzero==1),]$geol
  BOSstor.mat[i, which(colnames(BOSstor.mat) %in% BOSnzgeol)] <- 1
  
  # for non-zeros, what is the direction?
  BOSnzgeoldir <- paste(BOSnzgeol,"dir",sep="")
  BOSstor.mat[i, which(colnames(BOSstor.mat) %in% BOSnzgeoldir)] <- BOSpmmod1.coef$dir[which(BOSpmmod1.coef$dir != 0)]
  
  # AICc model comparison
  BOSwAICc <- weight.IC(delta.IC(c(AICc(BOSmod1),AICc(BOSmod.null))))
  
  # store evidence ratio
  BOSstor.mat[i,"ER"] <- BOSwAICc[1]/BOSwAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

head(BOSstor.mat)
geol.labs <- c("CAR","FEL","INT","MAF","MET","OTH","SIL")
BOSlennzCAR <- length(table(BOSstor.mat[,"geolCAR"]))
BOSlennzFEL <- length(table(BOSstor.mat[,"geolFEL"]))
BOSlennzINT <- length(table(BOSstor.mat[,"geolINT"]))
BOSlennzMAF <- length(table(BOSstor.mat[,"geolMAF"]))
BOSlennzMET <- length(table(BOSstor.mat[,"geolMET"]))
BOSlennzOTH <- length(table(BOSstor.mat[,"geolOTH"]))
BOSlennzSIL <- length(table(BOSstor.mat[,"geolSIL"]))

BOSnzrslts <- c(ifelse(BOSlennzCAR==1, as.numeric(table(BOSstor.mat[,"geolCAR"])), 0),
                ifelse(BOSlennzFEL==1, as.numeric(table(BOSstor.mat[,"geolFEL"])), 0),
                ifelse(BOSlennzINT==1, as.numeric(table(BOSstor.mat[,"geolINT"])), 0),
                ifelse(BOSlennzMAF==1, as.numeric(table(BOSstor.mat[,"geolMAF"])), 0),
                ifelse(BOSlennzMET==1, as.numeric(table(BOSstor.mat[,"geolMET"])), 0),
                ifelse(BOSlennzOTH==1, as.numeric(table(BOSstor.mat[,"geolOTH"])), 0),
                ifelse(BOSlennzSIL==1, as.numeric(table(BOSstor.mat[,"geolSIL"])), 0))

BOSdirCARtab <- table(BOSstor.mat[,"geolCARdir"])
BOSdirFELtab <- table(BOSstor.mat[,"geolFELdir"])
BOSdirINTtab <- table(BOSstor.mat[,"geolINTdir"])
BOSdirMAFtab <- table(BOSstor.mat[,"geolMAFdir"])
BOSdirMETtab <- table(BOSstor.mat[,"geolMETdir"])
BOSdirOTHtab <- table(BOSstor.mat[,"geolOTHdir"])
BOSdirSILtab <- table(BOSstor.mat[,"geolSILdir"])

BOSdirlenCAR <- length(BOSdirCARtab)
BOSdirlenFEL <- length(BOSdirFELtab)
BOSdirlenINT <- length(BOSdirINTtab)
BOSdirlenMAF <- length(BOSdirMAFtab)
BOSdirlenMET <- length(BOSdirMETtab)
BOSdirlenOTH <- length(BOSdirOTHtab)
BOSdirlenSIL <- length(BOSdirSILtab)


if (BOSdirlenCAR == 0) {
  BOSCARpos <- 0
  BOSCARneg <- 0
} else if (BOSdirlenCAR == 1 & as.numeric(attr(BOSdirCARtab,"names")[1]) == 1) {
  BOSCARpos <- as.numeric(BOSdirCARtab)[1]
  BOSCARneg <- 0
} else if (BOSdirlenCAR == 1 & as.numeric(attr(BOSdirCARtab,"names")[1]) == -1) {
  BOSCARpos <- 0
  BOSCARneg <- as.numeric(BOSdirCARtab)[1]
} else if (dirlenCAR == 2) {
  BOSCARneg <- as.numeric(BOSdirCARtab)[1]
  BOSCARpos <- as.numeric(BOSdirCARtab)[2]
} else {
  BOSCARpos <- 0
  BOSCARneg <- 0
}

if (BOSdirlenFEL == 0) {
  BOSFELpos <- 0
  BOSFELneg <- 0
} else if (BOSdirlenFEL == 1 & as.numeric(attr(BOSdirFELtab,"names")[1]) == 1) {
  BOSFELpos <- as.numeric(BOSdirFELtab)[1]
  BOSFELneg <- 0
} else if (BOSdirlenFEL == 1 & as.numeric(attr(BOSdirFELtab,"names")[1]) == -1) {
  BOSFELpos <- 0
  BOSFELneg <- as.numeric(BOSdirFELtab)[1]
} else if (BOSdirlenFEL == 2) {
  BOSFELneg <- as.numeric(BOSdirFELtab)[1]
  BOSFELpos <- as.numeric(BOSdirFELtab)[2]
} else {
  BOSFELpos <- 0
  BOSFELneg <- 0
}

if (BOSdirlenINT == 0) {
  BOSINTpos <- 0
  BOSINTneg <- 0
} else if (BOSdirlenINT == 1 & as.numeric(attr(BOSdirINTtab,"names")[1]) == 1) {
  BOSINTpos <- as.numeric(BOSdirINTtab)[1]
  BOSINTneg <- 0
} else if (BOSdirlenINT == 1 & as.numeric(attr(BOSdirINTtab,"names")[1]) == -1) {
  BOSINTpos <- 0
  BOSINTneg <- as.numeric(BOSdirINTtab)[1]
} else if (dirlenINT == 2) {
  BOSINTneg <- as.numeric(BOSdirINTtab)[1]
  BOSINTpos <- as.numeric(BOSdirINTtab)[2]
} else {
  BOSINTpos <- 0
  BOSINTneg <- 0
}

if (BOSdirlenMAF == 0) {
  BOSMAFpos <- 0
  BOSMAFneg <- 0
} else if (BOSdirlenMAF == 1 & as.numeric(attr(BOSdirMAFtab,"names")[1]) == 1) {
  BOSMAFpos <- as.numeric(BOSdirMAFtab)[1]
  BOSMAFneg <- 0
} else if (BOSdirlenMAF == 1 & as.numeric(attr(BOSdirMAFtab,"names")[1]) == -1) {
  BOSMAFpos <- 0
  BOSMAFneg <- as.numeric(BOSdirMAFtab)[1]
} else if (BOSdirlenMAF == 2) {
  BOSMAFneg <- as.numeric(BOSdirMAFtab)[1]
  BOSMAFpos <- as.numeric(BOSdirMAFtab)[2]
} else {
  BOSMAFpos <- 0
  BOSMAFneg <- 0
}

if (BOSdirlenMET == 0) {
  BOSMETpos <- 0
  BOSMETneg <- 0
} else if (BOSdirlenMET == 1 & as.numeric(attr(BOSdirMETtab,"names")[1]) == 1) {
  BOSMETpos <- as.numeric(BOSdirMETtab)[1]
  BOSMETneg <- 0
} else if (BOSdirlenMET == 1 & as.numeric(attr(BOSdirMETtab,"names")[1]) == -1) {
  BOSMETpos <- 0
  BOSMETneg <- as.numeric(BOSdirMETtab)[1]
} else if (BOSdirlenMET == 2) {
  BOSMETneg <- as.numeric(BOSdirMETtab)[1]
  BOSMETpos <- as.numeric(BOSdirMETtab)[2]
} else {
  BOSMETpos <- 0
  BOSMETneg <- 0
}

if (BOSdirlenOTH == 0) {
  BOSOTHpos <- 0
  BOSOTHneg <- 0
} else if (BOSdirlenOTH == 1 & as.numeric(attr(BOSdirOTHtab,"names")[1]) == 1) {
  BOSOTHpos <- as.numeric(BOSdirOTHtab)[1]
  BOSOTHneg <- 0
} else if (BOSdirlenOTH == 1 & as.numeric(attr(BOSdirOTHtab,"names")[1]) == -1) {
  BOSOTHpos <- 0
  BOSOTHneg <- as.numeric(BOSdirOTHtab)[1]
} else if (BOSdirlenOTH == 2) {
  BOSOTHneg <- as.numeric(BOSdirOTHtab)[1]
  BOSOTHpos <- as.numeric(BOSdirOTHtab)[2]
} else {
  BOSOTHpos <- 0
  BOSOTHneg <- 0
}

if (BOSdirlenSIL == 0) {
  BOSSILpos <- 0
  BOSSILneg <- 0
} else if (BOSdirlenSIL == 1 & as.numeric(attr(BOSdirSILtab,"names")[1]) == 1) {
  BOSSILpos <- as.numeric(BOSdirSILtab)[1]
  BOSSILneg <- 0
} else if (BOSdirlenSIL == 1 & as.numeric(attr(BOSdirSILtab,"names")[1]) == -1) {
  BOSSILpos <- 0
  BOSSILneg <- as.numeric(BOSdirSILtab)[1]
} else if (BOSdirlenSIL == 2) {
  BOSSILneg <- as.numeric(BOSdirSILtab)[1]
  BOSSILpos <- as.numeric(BOSdirSILtab)[2]
} else {
  BOSSILpos <- 0
  BOSSILneg <- 0
}

BOSnznegdirrslts <- c(BOSCARneg, BOSFELneg, BOSINTneg, BOSMAFneg, BOSMETneg, BOSOTHneg, BOSSILneg)
BOSnzposdirrslts <- c(BOSCARpos, BOSFELpos, BOSINTpos, BOSMAFpos, BOSMETpos, BOSOTHpos, BOSSILpos)

BOSresults.out <- data.frame(geol=BOSgeoldir.lab, nonzero=BOSnzrslts/iter, negdir=BOSnznegdirrslts/iter,
                             posdir=BOSnzposdirrslts/iter)
BOSresults.out

# ER
10^median(log10(BOSstor.mat[,"ER"]))
quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(BOSstor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(BOSstor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(BOSstor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(BOSstor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)


#########
# soil ##
#########

## TOS
table(TOS.ran.subset$soil)

table(TOS.ran.subset$soil)
TOSHgXsoil.stats.ran <- TOS.ran.subset %>%
  group_by(soil) %>%
  summarise(
    mean = mean(lHg, na.rm = TRUE),
    median = median(lHg, na.rm = TRUE), 
    var = var(lHg, na.rm = TRUE),
    sd = sd(lHg, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(lHg, probs=0.975, na.rm = TRUE),
    lower = quantile(lHg, probs=0.025, na.rm = TRUE),
    n = n()
  )
TOSHgXsoil.stats.ran
na.omit(TOSHgXsoil.stats.ran)

# resampling loop
# storage matrix
table(TOS.lm$soil)

TOS.lm.soils.nolake <- TOS.lm[which(TOS.lm$soil != "lake"),]
table(TOS.lm.soils.nolake$soil)

# remove ferrosols (n=1) and organosols (n=2)
# TOS.lm.soils.nolakemod <- TOS.lm[which(TOS.lm.soils.nolake$soil != "ferrosol" & TOS.lm.soils.nolake$soil != "organosol"),]
# dim(TOS.lm.soils.nolakemod)
# table(TOS.lm.soils.nolakemod$soil)
TOS.lm.soils.nolakemod <- TOS.lm.soils.nolake

TOSHgXsoil.stats.nolake <- TOS.lm.soils.nolakemod %>%
  group_by(soil) %>%
  summarise(
    mean = mean(10^lHg, na.rm = TRUE),
    median = median(10^lHg, na.rm = TRUE), 
    var = var(10^lHg, na.rm = TRUE),
    sd = sd(10^lHg, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(10^lHg, probs=0.975, na.rm = TRUE),
    lower = quantile(10^lHg, probs=0.025, na.rm = TRUE),
    n = n()
  )
TOSHgXsoil.stats.nolake
data.frame(soil=TOSHgXsoil.stats.nolake$soil, mean=TOSHgXsoil.stats.nolake$mean, se=TOSHgXsoil.stats.nolake$se,
           upper=TOSHgXsoil.stats.nolake$upper, lower=TOSHgXsoil.stats.nolake$lower)

TOSnlevels <- length(table(TOS.lm.soils.nolakemod$soil))
TOSstor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*TOSnlevels)
TOSsoildir.lab <- paste("soil",attr(table(TOS.lm.soils.nolakemod$soil), "names"),"dir",sep="")
colnames(TOSstor.mat) <- c("ER", "nonzerosum",
                        paste("soil",attr(table(TOS.lm.soils.nolakemod$soil), "names"),sep=""),TOSsoildir.lab)
head(TOSstor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = TOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = TOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  TOS.ran.subset <- na.omit(TOS.lm.soils.nolakemod[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  TOSmod1 <- lm(lHg ~ soil, data=TOS.ran.subset) # class level model
  TOSmod.null <- lm(lHg ~ 1, data=TOS.ran.subset) # null model
  
  # coefficient boundaries
  TOSpmmod1 <- plot_model(TOSmod1)
  TOSpmmod1.coef <- data.frame(soil=as.character(TOSpmmod1[[1]]$term), lo=TOSpmmod1[[1]]$conf.low,
                            up=TOSpmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  TOSpmmod1.coef$nonzero <- ifelse(contains_zero_sign(TOSpmmod1.coef$lo, TOSpmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  TOSpmmod1.coef$dir <- ifelse(sign(TOSpmmod1.coef$lo) == -1 & sign(TOSpmmod1.coef$up) == -1, -1, 0)
  TOSpmmod1.coef$dir <- ifelse(sign(TOSpmmod1.coef$lo) == 1 & sign(TOSpmmod1.coef$up) == 1, 1, TOSpmmod1.coef$dir)
  
  # how many variables are non-zero?
  TOSstor.mat[i,"nonzerosum"] <- sum(TOSpmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  TOSnzsoils <- TOSpmmod1.coef[which(TOSpmmod1.coef$nonzero==1),]$soil
  TOSstor.mat[i, which(colnames(TOSstor.mat) %in% TOSnzsoils)] <- 1
  
  # for non-zeros, what is the direction?
  TOSnzsoilssdir <- paste(TOSnzsoils,"dir",sep="")
  TOSstor.mat[i, which(colnames(TOSstor.mat) %in% TOSnzsoilssdir)] <- TOSpmmod1.coef$dir[which(TOSpmmod1.coef$dir != 0)]
  
  # AICc model comparison
  TOSwAICc <- weight.IC(delta.IC(c(AICc(TOSmod1),AICc(TOSmod.null))))
  
  # store evidence ratio
  TOSstor.mat[i,"ER"] <- TOSwAICc[1]/TOSwAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

head(TOSstor.mat)
colnames(TOSstor.mat)
TOSsoil.labs <- c("CALC","CHROMO","DERMO","FERRO","HYDRO","KANDO","KURO","ORGANO","PODO","RUDO","SODO","TENO","VERTO")
#TOSsoil.labs <- c("CALC","CHROMO","DERMO","HYDRO","KANDO","KURO","PODO","RUDO","SODO","TENO","VERTO")

TOSlennzCALC <- length(table(TOSstor.mat[,"soilcalcarosol"]))
TOSlennzCHROMO <- length(table(TOSstor.mat[,"soilchromosol"]))
TOSlennzDERMO <- length(table(TOSstor.mat[,"soildermosol"]))
TOSlennzFERRO <- length(table(TOSstor.mat[,"soilferrosol"]))
TOSlennzHYDRO <- length(table(TOSstor.mat[,"soilhydrosol"]))
TOSlennzKANDO <- length(table(TOSstor.mat[,"soilkandosol"]))
TOSlennzKURO <- length(table(TOSstor.mat[,"soilkurosol"]))
TOSlennzORGANO <- length(table(TOSstor.mat[,"soilorganosol"]))
TOSlennzPODO <- length(table(TOSstor.mat[,"soilpodosol"]))
TOSlennzRUDO <- length(table(TOSstor.mat[,"soilrudosol"]))
TOSlennzSODO <- length(table(TOSstor.mat[,"soilsodosol"]))
TOSlennzTENO <- length(table(TOSstor.mat[,"soiltenosol"]))
TOSlennzVERTO<- length(table(TOSstor.mat[,"soilvertosol"]))

TOSnzrslts <- c(ifelse(TOSlennzCALC==1, as.numeric(table(TOSstor.mat[,"soilcalcarosol"])), 0),
             ifelse(TOSlennzCHROMO==1, as.numeric(table(TOSstor.mat[,"soilchromosol"])), 0),
             ifelse(TOSlennzDERMO==1, as.numeric(table(TOSstor.mat[,"soildermosol"])), 0),
             ifelse(TOSlennzFERRO==1, as.numeric(table(TOSstor.mat[,"soilferrosol"])), 0),
             ifelse(TOSlennzHYDRO==1, as.numeric(table(TOSstor.mat[,"soilhydrosol"])), 0),
             ifelse(TOSlennzKANDO==1, as.numeric(table(TOSstor.mat[,"soilkandosol"])), 0),
             ifelse(TOSlennzKURO==1, as.numeric(table(TOSstor.mat[,"soilkurosol"])), 0),
             ifelse(TOSlennzORGANO==1, as.numeric(table(TOSstor.mat[,"soilorganosol"])), 0),
             ifelse(TOSlennzPODO==1, as.numeric(table(TOSstor.mat[,"soilpodosol"])), 0),
             ifelse(TOSlennzRUDO==1, as.numeric(table(TOSstor.mat[,"soilrudosol"])), 0),
             ifelse(TOSlennzSODO==1, as.numeric(table(TOSstor.mat[,"soilsodosol"])), 0),
             ifelse(TOSlennzTENO==1, as.numeric(table(TOSstor.mat[,"soiltenosol"])), 0),
             ifelse(TOSlennzVERTO==1, as.numeric(table(TOSstor.mat[,"soilvertosol"])), 0))

TOSdirCALCtab <- table(TOSstor.mat[,"soilcalcarosoldir"])
TOSdirCHROMOtab <- table(TOSstor.mat[,"soilchromosoldir"])
TOSdirDERMOtab <- table(TOSstor.mat[,"soildermosoldir"])
TOSdirFERROtab <- table(TOSstor.mat[,"soilferrosoldir"])
TOSdirHYDROtab <- table(TOSstor.mat[,"soilhydrosoldir"])
TOSdirKANDOtab <- table(TOSstor.mat[,"soilkandosoldir"])
TOSdirKUROtab <- table(TOSstor.mat[,"soilkurosoldir"])
TOSdirORGANOTab <- table(TOSstor.mat[,"soilorganosoldir"])
TOSdirPODOtab <- table(TOSstor.mat[,"soilpodosoldir"])
TOSdirRUDOtab <- table(TOSstor.mat[,"soilrudosoldir"])
TOSdirSODOtab <- table(TOSstor.mat[,"soilsodosoldir"])
TOSdirTENOTab <- table(TOSstor.mat[,"soiltenosoldir"])
TOSdirVERTOTab <- table(TOSstor.mat[,"soilvertosoldir"])

TOSdirlenCALC <- length(TOSdirCALCtab)
TOSdirlenCHROMO <- length(TOSdirCHROMOtab)
TOSdirlenDERMO <- length(TOSdirDERMOtab)
TOSdirlenFERRO <- length(TOSdirFERROtab)
TOSdirlenHYDRO <- length(TOSdirHYDROtab)
TOSdirlenKANDO <- length(TOSdirKANDOtab)
TOSdirlenKURO <- length(TOSdirKUROtab)
TOSdirlenORGANO <- length(TOSdirORGANOTab)
TOSdirlenPODO <- length(TOSdirPODOtab)
TOSdirlenRUDO <- length(TOSdirRUDOtab)
TOSdirlenSODO <- length(TOSdirSODOtab)
TOSdirlenTENO <- length(TOSdirTENOTab)
TOSdirlenVERTO <- length(TOSdirVERTOTab)

if (TOSdirlenCALC == 0) {
  TOSCALCpos <- 0
  TOSCALCneg <- 0
} else if (TOSdirlenCALC == 1 & as.numeric(attr(TOSdirCALCtab,"names")[1]) == 1) {
  TOSCALCpos <- as.numeric(TOSdirCALCtab)[1]
  TOSCALCneg <- 0
} else if (TOSdirlenCALC == 1 & as.numeric(attr(TOSdirCALCtab,"names")[1]) == -1) {
  TOSCALCpos <- 0
  TOSCALCneg <- as.numeric(TOSdirCALCtab)[1]
} else if (dirlenCALC == 2) {
  TOSCALCneg <- as.numeric(TOSdirCALCtab)[1]
  TOSCALCpos <- as.numeric(TOSdirCALCtab)[2]
} else {
  TOSCALCpos <- 0
  TOSCALCneg <- 0
}

if (TOSdirlenCHROMO == 0) {
  TOSCHROMOpos <- 0
  TOSCHROMOneg <- 0
} else if (TOSdirlenCHROMO == 1 & as.numeric(attr(TOSdirCHROMOtab,"names")[1]) == 1) {
  TOSCHROMOpos <- as.numeric(TOSdirCHROMOtab)[1]
  TOSCHROMOneg <- 0
} else if (TOSdirlenCHROMO == 1 & as.numeric(attr(TOSdirCHROMOtab,"names")[1]) == -1) {
  TOSCHROMOpos <- 0
  TOSCHROMOneg <- as.numeric(TOSdirCHROMOtab)[1]
} else if (TOSdirlenCHROMO == 2) {
  TOSCHROMOneg <- as.numeric(TOSdirCHROMOtab)[1]
  TOSCHROMOpos <- as.numeric(TOSdirCHROMOtab)[2]
} else {
  TOSCHROMOpos <- 0
  TOSCHROMOneg <- 0
}

if (TOSdirlenDERMO == 0) {
  TOSDERMOpos <- 0
  TOSDERMOneg <- 0
} else if (TOSdirlenDERMO == 1 & as.numeric(attr(TOSdirDERMOtab,"names")[1]) == 1) {
  TOSDERMOpos <- as.numeric(TOSdirDERMOtab)[1]
  TOSDERMOneg <- 0
} else if (TOSdirlenDERMO == 1 & as.numeric(attr(TOSdirDERMOtab,"names")[1]) == -1) {
  TOSDERMOpos <- 0
  TOSDERMOneg <- as.numeric(TOSdirDERMOtab)[1]
} else if (TOSdirlenDERMO == 2) {
  TOSDERMOneg <- as.numeric(TOSdirDERMOtab)[1]
  TOSDERMOpos <- as.numeric(TOSdirDERMOtab)[2]
} else {
  TOSDERMOpos <- 0
  TOSDERMOneg <- 0
}

if (TOSdirlenFERRO == 0) {
  TOSFERROpos <- 0
  TOSFERROneg <- 0
} else if (TOSdirlenFERRO == 1 & as.numeric(attr(TOSdirFERROtab,"names")[1]) == 1) {
  TOSFERROpos <- as.numeric(TOSdirFERROtab)[1]
  TOSFERROneg <- 0
} else if (TOSdirlenFERRO == 1 & as.numeric(attr(TOSdirFERROtab,"names")[1]) == -1) {
  TOSFERROpos <- 0
  TOSFERROneg <- as.numeric(TOSdirFERROtab)[1]
} else if (TOSdirlenFERRO == 2) {
  TOSFERROneg <- as.numeric(TOSdirFERROtab)[1]
  TOSFERROpos <- as.numeric(TOSdirFERROtab)[2]
} else {
  TOSFERROpos <- 0
  TOSFERROneg <- 0
}

if (TOSdirlenHYDRO == 0) {
  TOSHYDROpos <- 0
  TOSHYDROneg <- 0
} else if (TOSdirlenHYDRO == 1 & as.numeric(attr(TOSdirHYDROtab,"names")[1]) == 1) {
  TOSHYDROpos <- as.numeric(TOSdirHYDROtab)[1]
  TOSHYDROneg <- 0
} else if (TOSdirlenHYDRO == 1 & as.numeric(attr(TOSdirHYDROtab,"names")[1]) == -1) {
  TOSHYDROpos <- 0
  TOSHYDROneg <- as.numeric(TOSdirHYDROtab)[1]
} else if (TOSdirlenHYDRO == 2) {
  TOSHYDROneg <- as.numeric(TOSdirHYDROtab)[1]
  TOSHYDROpos <- as.numeric(TOSdirHYDROtab)[2]
} else {
  TOSHYDROpos <- 0
  TOSHYDROneg <- 0
}

if (TOSdirlenKANDO == 0) {
  TOSKANDOpos <- 0
  TOSKANDOneg <- 0
} else if (TOSdirlenKANDO == 1 & as.numeric(attr(TOSdirKANDOtab,"names")[1]) == 1) {
  TOSKANDOpos <- as.numeric(TOSdirKANDOtab)[1]
  TOSKANDOneg <- 0
} else if (TOSdirlenKANDO == 1 & as.numeric(attr(TOSdirKANDOtab,"names")[1]) == -1) {
  TOSKANDOpos <- 0
  TOSKANDOneg <- as.numeric(TOSdirKANDOtab)[1]
} else if (TOSdirlenKANDO == 2) {
  TOSKANDOneg <- as.numeric(TOSdirKANDOtab)[1]
  TOSKANDOpos <- as.numeric(TOSdirKANDOtab)[2]
} else {
  TOSKANDOpos <- 0
  TOSKANDOneg <- 0
}

if (TOSdirlenKURO == 0) {
  TOSKUROpos <- 0
  TOSKUROneg <- 0
} else if (TOSdirlenKURO == 1 & as.numeric(attr(TOSdirKUROtab,"names")[1]) == 1) {
  TOSKUROpos <- as.numeric(TOSdirKUROtab)[1]
  TOSKUROneg <- 0
} else if (TOSdirlenKURO == 1 & as.numeric(attr(TOSdirKUROtab,"names")[1]) == -1) {
  TOSKUROpos <- 0
  TOSKUROneg <- as.numeric(TOSdirKUROtab)[1]
} else if (TOSdirlenKURO == 2) {
  TOSKUROneg <- as.numeric(TOSdirKUROtab)[1]
  TOSKUROpos <- as.numeric(TOSdirKUROtab)[2]
} else {
  TOSKUROpos <- 0
  TOSKUROneg <- 0
}

if (TOSdirlenORGANO == 0) {
  TOSORGANOpos <- 0
  TOSORGANOneg <- 0
} else if (TOSdirlenORGANO == 1 & as.numeric(attr(TOSdirORGANOTab,"names")[1]) == 1) {
  TOSORGANOpos <- as.numeric(TOSdirORGANOTab)[1]
  TOSORGANOneg <- 0
} else if (TOSdirlenORGANO == 1 & as.numeric(attr(TOSdirORGANOTab,"names")[1]) == -1) {
  TOSORGANOpos <- 0
  TOSORGANOneg <- as.numeric(TOSdirORGANOTab)[1]
} else if (TOSdirlenORGANO == 2) {
  TOSORGANOneg <- as.numeric(TOSdirORGANOTab)[1]
  TOSORGANOpos <- as.numeric(TOSdirORGANOTab)[2]
} else {
  TOSORGANOpos <- 0
  TOSORGANOneg <- 0
}

if (TOSdirlenPODO == 0) {
  TOSPODOpos <- 0
  TOSPODOneg <- 0
} else if (TOSdirlenPODO == 1 & as.numeric(attr(TOSdirPODOtab,"names")[1]) == 1) {
  TOSPODOpos <- as.numeric(TOSdirPODOtab)[1]
  TOSPODOneg <- 0
} else if (TOSdirlenPODO == 1 & as.numeric(attr(TOSdirPODOtab,"names")[1]) == -1) {
  TOSPODOpos <- 0
  TOSPODOneg <- as.numeric(TOSdirPODOtab)[1]
} else if (TOSdirlenPODO == 2) {
  TOSPODOneg <- as.numeric(TOSdirPODOtab)[1]
  TOSPODOpos <- as.numeric(TOSdirPODOtab)[2]
} else {
  TOSPODOpos <- 0
  TOSPODOneg <- 0
}

if (TOSdirlenRUDO == 0) {
  TOSRUDOpos <- 0
  TOSRUDOneg <- 0
} else if (TOSdirlenRUDO == 1 & as.numeric(attr(TOSdirRUDOtab,"names")[1]) == 1) {
  TOSRUDOpos <- as.numeric(TOSdirRUDOtab)[1]
  TOSRUDOneg <- 0
} else if (TOSdirlenRUDO == 1 & as.numeric(attr(TOSdirRUDOtab,"names")[1]) == -1) {
  TOSRUDOpos <- 0
  TOSRUDOneg <- as.numeric(TOSdirRUDOtab)[1]
} else if (TOSdirlenRUDO == 2) {
  TOSRUDOneg <- as.numeric(TOSdirRUDOtab)[1]
  TOSRUDOpos <- as.numeric(TOSdirRUDOtab)[2]
} else {
  TOSRUDOpos <- 0
  TOSRUDOneg <- 0
}

if (TOSdirlenSODO == 0) {
  TOSSODOpos <- 0
  TOSSODOneg <- 0
} else if (TOSdirlenSODO == 1 & as.numeric(attr(TOSdirSODOtab,"names")[1]) == 1) {
  TOSSODOpos <- as.numeric(TOSdirSODOtab)[1]
  TOSSODOneg <- 0
} else if (TOSdirlenSODO == 1 & as.numeric(attr(TOSdirSODOtab,"names")[1]) == -1) {
  TOSSODOpos <- 0
  TOSSODOneg <- as.numeric(TOSdirSODOtab)[1]
} else if (TOSdirlenSODO == 2) {
  TOSSODOneg <- as.numeric(TOSdirSODOtab)[1]
  TOSSODOpos <- as.numeric(TOSdirSODOtab)[2]
} else {
  TOSSODOpos <- 0
  TOSSODOneg <- 0
}

if (TOSdirlenTENO == 0) {
  TOSTENOpos <- 0
  TOSTENOneg <- 0
} else if (TOSdirlenTENO == 1 & as.numeric(attr(TOSdirTENOTab,"names")[1]) == 1) {
  TOSTENOpos <- as.numeric(TOSdirTENOTab)[1]
  TOSTENOneg <- 0
} else if (TOSdirlenTENO == 1 & as.numeric(attr(TOSdirTENOTab,"names")[1]) == -1) {
  TOSTENOpos <- 0
  TOSTENOneg <- as.numeric(TOSdirTENOTab)[1]
} else if (TOSdirlenTENO == 2) {
  TOSTENOneg <- as.numeric(TOSdirTENOTab)[1]
  TOSTENOpos <- as.numeric(TOSdirTENOTab)[2]
} else {
  TOSTENOpos <- 0
  TOSTENOneg <- 0
}

if (TOSdirlenVERTO == 0) {
  TOSVERTOpos <- 0
  TOSVERTOneg <- 0
} else if (TOSdirlenVERTO == 1 & as.numeric(attr(TOSdirVERTOTab,"names")[1]) == 1) {
  TOSVERTOpos <- as.numeric(TOSdirVERTOTab)[1]
  TOSVERTOneg <- 0
} else if (TOSdirlenVERTO == 1 & as.numeric(attr(TOSdirVERTOTab,"names")[1]) == -1) {
  TOSVERTOpos <- 0
  TOSVERTOneg <- as.numeric(TOSdirVERTOTab)[1]
} else if (TOSdirlenVERTO == 2) {
  TOSVERTOneg <- as.numeric(TOSdirVERTOTab)[1]
  TOSVERTOpos <- as.numeric(TOSdirVERTOTab)[2]
} else {
  TOSVERTOpos <- 0
  TOSVERTOneg <- 0
}

TOSnznegdirrslts <- c(TOSCALCneg, TOSCHROMOneg, TOSDERMOneg, TOSFERROneg, TOSHYDROneg, TOSKANDOneg, TOSKUROneg, TOSORGANOneg, TOSPODOneg, TOSRUDOneg, TOSSODOneg, TOSTENOneg, TOSVERTOneg)
TOSnzposdirrslts <- c(TOSCALCpos, TOSCHROMOpos, TOSDERMOpos, TOSFERROpos, TOSHYDROpos, TOSKANDOpos, TOSKUROpos, TOSORGANOpos, TOSPODOpos, TOSRUDOpos, TOSSODOpos, TOSTENOpos, TOSVERTOpos)
#TOSnznegdirrslts <- c(TOSCALCneg, TOSCHROMOneg, TOSDERMOneg, TOSHYDROneg, TOSKANDOneg, TOSKUROneg, TOSPODOneg, TOSRUDOneg, TOSSODOneg, TOSTENOneg, TOSVERTOneg)
#TOSnzposdirrslts <- c(TOSCALCpos, TOSCHROMOpos, TOSDERMOpos, TOSHYDROpos, TOSKANDOpos, TOSKUROpos, TOSPODOpos, TOSRUDOpos, TOSSODOpos, TOSTENOpos, TOSVERTOpos)

TOSresults.out <- data.frame(soil=TOSsoildir.lab, nonzero=TOSnzrslts/iter, negdir=TOSnznegdirrslts/iter,
                          posdir=TOSnzposdirrslts/iter)
TOSresults.out

# ER
10^median(log10(TOSstor.mat[,"ER"]))
quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(TOSstor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(TOSstor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(TOSstor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(TOSstor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(TOSstor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)


## BOS
table(BOS.ran.subset$soil)

table(BOS.ran.subset$soil)
BOSHgXsoil.stats.ran <- BOS.ran.subset %>%
  group_by(soil) %>%
  summarise(
    mean = mean(lHg, na.rm = TRUE),
    median = median(lHg, na.rm = TRUE), 
    var = var(lHg, na.rm = TRUE),
    sd = sd(lHg, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(lHg, probs=0.975, na.rm = TRUE),
    lower = quantile(lHg, probs=0.025, na.rm = TRUE),
    n = n()
  )
BOSHgXsoil.stats.ran
na.omit(BOSHgXsoil.stats.ran)

# resampling loop
# storage matrix
table(BOS.lm$soil)
BOS.lm.soils.nolake <- BOS.lm[which(BOS.lm$soil != "lake"),]
table(BOS.lm.soils.nolake$soil)

BOSHgXsoil.stats.nolake <- BOS.lm.soils.nolake %>%
  group_by(soil) %>%
  summarise(
    mean = mean(10^lHg, na.rm = TRUE),
    median = median(10^lHg, na.rm = TRUE), 
    var = var(10^lHg, na.rm = TRUE),
    sd = sd(10^lHg, na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(10^lHg, probs=0.975, na.rm = TRUE),
    lower = quantile(10^lHg, probs=0.025, na.rm = TRUE),
    n = n()
  )
BOSHgXsoil.stats.nolake
data.frame(soil=BOSHgXsoil.stats.nolake$soil, mean=BOSHgXsoil.stats.nolake$mean, se=BOSHgXsoil.stats.nolake$se,
           upper=BOSHgXsoil.stats.nolake$upper, lower=BOSHgXsoil.stats.nolake$lower)

BOSnlevels <- length(table(BOS.lm.soils.nolake$soil))
BOSstor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*BOSnlevels)
BOSsoildir.lab <- paste("soil",attr(table(BOS.lm.soils.nolake$soil), "names"),"dir",sep="")
colnames(BOSstor.mat) <- c("ER", "nonzerosum",
                           paste("soil",attr(table(BOS.lm.soils.nolake$soil), "names"),sep=""),BOSsoildir.lab)
head(BOSstor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = BOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = BOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  BOS.ran.subset <- na.omit(BOS.lm.soils.nolake[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  BOSmod1 <- lm(lHg ~ soil, data=BOS.ran.subset) # class level model
  BOSmod.null <- lm(lHg ~ 1, data=BOS.ran.subset) # null model
  
  # coefficient boundaries
  BOSpmmod1 <- plot_model(BOSmod1)
  BOSpmmod1.coef <- data.frame(soil=as.character(BOSpmmod1[[1]]$term), lo=BOSpmmod1[[1]]$conf.low,
                               up=BOSpmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  BOSpmmod1.coef$nonzero <- ifelse(contains_zero_sign(BOSpmmod1.coef$lo, BOSpmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  BOSpmmod1.coef$dir <- ifelse(sign(BOSpmmod1.coef$lo) == -1 & sign(BOSpmmod1.coef$up) == -1, -1, 0)
  BOSpmmod1.coef$dir <- ifelse(sign(BOSpmmod1.coef$lo) == 1 & sign(BOSpmmod1.coef$up) == 1, 1, BOSpmmod1.coef$dir)
  
  # how many variables are non-zero?
  BOSstor.mat[i,"nonzerosum"] <- sum(BOSpmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  BOSnzsoils <- BOSpmmod1.coef[which(BOSpmmod1.coef$nonzero==1),]$soil
  BOSstor.mat[i, which(colnames(BOSstor.mat) %in% BOSnzsoils)] <- 1
  
  # for non-zeros, what is the direction?
  BOSnzsoilssdir <- paste(BOSnzsoils,"dir",sep="")
  BOSstor.mat[i, which(colnames(BOSstor.mat) %in% BOSnzsoilssdir)] <- BOSpmmod1.coef$dir[which(BOSpmmod1.coef$dir != 0)]
  
  # AICc model comparison
  BOSwAICc <- weight.IC(delta.IC(c(AICc(BOSmod1),AICc(BOSmod.null))))
  
  # store evidence ratio
  BOSstor.mat[i,"ER"] <- BOSwAICc[1]/BOSwAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

head(BOSstor.mat)
colnames(BOSstor.mat)
BOSsoil.labs <- c("CALC","CHROMO","DERMO","FERRO","HYDRO","KANDO","KURO","ORGANO","PODO","RUDO","SODO","TENO","VERTO")
BOSlennzCALC <- length(table(BOSstor.mat[,"soilcalcarosol"]))
BOSlennzCHROMO <- length(table(BOSstor.mat[,"soilchromosol"]))
BOSlennzDERMO <- length(table(BOSstor.mat[,"soildermosol"]))
BOSlennzFERRO <- length(table(BOSstor.mat[,"soilferrosol"]))
BOSlennzHYDRO <- length(table(BOSstor.mat[,"soilhydrosol"]))
BOSlennzKANDO <- length(table(BOSstor.mat[,"soilkandosol"]))
BOSlennzKURO <- length(table(BOSstor.mat[,"soilkurosol"]))
BOSlennzORGANO <- length(table(BOSstor.mat[,"soilorganosol"]))
BOSlennzPODO <- length(table(BOSstor.mat[,"soilpodosol"]))
BOSlennzRUDO <- length(table(BOSstor.mat[,"soilrudosol"]))
BOSlennzSODO <- length(table(BOSstor.mat[,"soilsodosol"]))
BOSlennzTENO <- length(table(BOSstor.mat[,"soiltenosol"]))
BOSlennzVERTO<- length(table(BOSstor.mat[,"soilvertosol"]))

BOSnzrslts <- c(ifelse(BOSlennzCALC==1, as.numeric(table(BOSstor.mat[,"soilcalcarosol"])), 0),
                ifelse(BOSlennzCHROMO==1, as.numeric(table(BOSstor.mat[,"soilchromosol"])), 0),
                ifelse(BOSlennzDERMO==1, as.numeric(table(BOSstor.mat[,"soildermosol"])), 0),
                ifelse(BOSlennzFERRO==1, as.numeric(table(BOSstor.mat[,"soilferrosol"])), 0),
                ifelse(BOSlennzHYDRO==1, as.numeric(table(BOSstor.mat[,"soilhydrosol"])), 0),
                ifelse(BOSlennzKANDO==1, as.numeric(table(BOSstor.mat[,"soilkandosol"])), 0),
                ifelse(BOSlennzKURO==1, as.numeric(table(BOSstor.mat[,"soilkurosol"])), 0),
                ifelse(BOSlennzORGANO==1, as.numeric(table(BOSstor.mat[,"soilorganosol"])), 0),
                ifelse(BOSlennzPODO==1, as.numeric(table(BOSstor.mat[,"soilpodosol"])), 0),
                ifelse(BOSlennzRUDO==1, as.numeric(table(BOSstor.mat[,"soilrudosol"])), 0),
                ifelse(BOSlennzSODO==1, as.numeric(table(BOSstor.mat[,"soilsodosol"])), 0),
                ifelse(BOSlennzTENO==1, as.numeric(table(BOSstor.mat[,"soiltenosol"])), 0),
                ifelse(BOSlennzVERTO==1, as.numeric(table(BOSstor.mat[,"soilvertosol"])), 0))

BOSdirCALCtab <- table(BOSstor.mat[,"soilcalcarosoldir"])
BOSdirCHROMOtab <- table(BOSstor.mat[,"soilchromosoldir"])
BOSdirDERMOtab <- table(BOSstor.mat[,"soildermosoldir"])
BOSdirFERROtab <- table(BOSstor.mat[,"soilferrosoldir"])
BOSdirHYDROtab <- table(BOSstor.mat[,"soilhydrosoldir"])
BOSdirKANDOtab <- table(BOSstor.mat[,"soilkandosoldir"])
BOSdirKUROtab <- table(BOSstor.mat[,"soilkurosoldir"])
BOSdirORGANOTab <- table(BOSstor.mat[,"soilorganosoldir"])
BOSdirPODOtab <- table(BOSstor.mat[,"soilpodosoldir"])
BOSdirRUDOtab <- table(BOSstor.mat[,"soilrudosoldir"])
BOSdirSODOtab <- table(BOSstor.mat[,"soilsodosoldir"])
BOSdirTENOTab <- table(BOSstor.mat[,"soiltenosoldir"])
BOSdirVERTOTab <- table(BOSstor.mat[,"soilvertosoldir"])

BOSdirlenCALC <- length(BOSdirCALCtab)
BOSdirlenCHROMO <- length(BOSdirCHROMOtab)
BOSdirlenDERMO <- length(BOSdirDERMOtab)
BOSdirlenFERRO <- length(BOSdirFERROtab)
BOSdirlenHYDRO <- length(BOSdirHYDROtab)
BOSdirlenKANDO <- length(BOSdirKANDOtab)
BOSdirlenKURO <- length(BOSdirKUROtab)
BOSdirlenORGANO <- length(BOSdirORGANOTab)
BOSdirlenPODO <- length(BOSdirPODOtab)
BOSdirlenRUDO <- length(BOSdirRUDOtab)
BOSdirlenSODO <- length(BOSdirSODOtab)
BOSdirlenTENO <- length(BOSdirTENOTab)
BOSdirlenVERTO <- length(BOSdirVERTOTab)

if (BOSdirlenCALC == 0) {
  BOSCALCpos <- 0
  BOSCALCneg <- 0
} else if (BOSdirlenCALC == 1 & as.numeric(attr(BOSdirCALCtab,"names")[1]) == 1) {
  BOSCALCpos <- as.numeric(BOSdirCALCtab)[1]
  BOSCALCneg <- 0
} else if (BOSdirlenCALC == 1 & as.numeric(attr(BOSdirCALCtab,"names")[1]) == -1) {
  BOSCALCpos <- 0
  BOSCALCneg <- as.numeric(BOSdirCALCtab)[1]
} else if (dirlenCALC == 2) {
  BOSCALCneg <- as.numeric(BOSdirCALCtab)[1]
  BOSCALCpos <- as.numeric(BOSdirCALCtab)[2]
} else {
  BOSCALCpos <- 0
  BOSCALCneg <- 0
}

if (BOSdirlenCHROMO == 0) {
  BOSCHROMOpos <- 0
  BOSCHROMOneg <- 0
} else if (BOSdirlenCHROMO == 1 & as.numeric(attr(BOSdirCHROMOtab,"names")[1]) == 1) {
  BOSCHROMOpos <- as.numeric(BOSdirCHROMOtab)[1]
  BOSCHROMOneg <- 0
} else if (BOSdirlenCHROMO == 1 & as.numeric(attr(BOSdirCHROMOtab,"names")[1]) == -1) {
  BOSCHROMOpos <- 0
  BOSCHROMOneg <- as.numeric(BOSdirCHROMOtab)[1]
} else if (BOSdirlenCHROMO == 2) {
  BOSCHROMOneg <- as.numeric(BOSdirCHROMOtab)[1]
  BOSCHROMOpos <- as.numeric(BOSdirCHROMOtab)[2]
} else {
  BOSCHROMOpos <- 0
  BOSCHROMOneg <- 0
}

if (BOSdirlenDERMO == 0) {
  BOSDERMOpos <- 0
  BOSDERMOneg <- 0
} else if (BOSdirlenDERMO == 1 & as.numeric(attr(BOSdirDERMOtab,"names")[1]) == 1) {
  BOSDERMOpos <- as.numeric(BOSdirDERMOtab)[1]
  BOSDERMOneg <- 0
} else if (BOSdirlenDERMO == 1 & as.numeric(attr(BOSdirDERMOtab,"names")[1]) == -1) {
  BOSDERMOpos <- 0
  BOSDERMOneg <- as.numeric(BOSdirDERMOtab)[1]
} else if (BOSdirlenDERMO == 2) {
  BOSDERMOneg <- as.numeric(BOSdirDERMOtab)[1]
  BOSDERMOpos <- as.numeric(BOSdirDERMOtab)[2]
} else {
  BOSDERMOpos <- 0
  BOSDERMOneg <- 0
}

if (BOSdirlenFERRO == 0) {
  BOSFERROpos <- 0
  BOSFERROneg <- 0
} else if (BOSdirlenFERRO == 1 & as.numeric(attr(BOSdirFERROtab,"names")[1]) == 1) {
  BOSFERROpos <- as.numeric(BOSdirFERROtab)[1]
  BOSFERROneg <- 0
} else if (BOSdirlenFERRO == 1 & as.numeric(attr(BOSdirFERROtab,"names")[1]) == -1) {
  BOSFERROpos <- 0
  BOSFERROneg <- as.numeric(BOSdirFERROtab)[1]
} else if (BOSdirlenFERRO == 2) {
  BOSFERROneg <- as.numeric(BOSdirFERROtab)[1]
  BOSFERROpos <- as.numeric(BOSdirFERROtab)[2]
} else {
  BOSFERROpos <- 0
  BOSFERROneg <- 0
}

if (BOSdirlenHYDRO == 0) {
  BOSHYDROpos <- 0
  BOSHYDROneg <- 0
} else if (BOSdirlenHYDRO == 1 & as.numeric(attr(BOSdirHYDROtab,"names")[1]) == 1) {
  BOSHYDROpos <- as.numeric(BOSdirHYDROtab)[1]
  BOSHYDROneg <- 0
} else if (BOSdirlenHYDRO == 1 & as.numeric(attr(BOSdirHYDROtab,"names")[1]) == -1) {
  BOSHYDROpos <- 0
  BOSHYDROneg <- as.numeric(BOSdirHYDROtab)[1]
} else if (BOSdirlenHYDRO == 2) {
  BOSHYDROneg <- as.numeric(BOSdirHYDROtab)[1]
  BOSHYDROpos <- as.numeric(BOSdirHYDROtab)[2]
} else {
  BOSHYDROpos <- 0
  BOSHYDROneg <- 0
}

if (BOSdirlenKANDO == 0) {
  BOSKANDOpos <- 0
  BOSKANDOneg <- 0
} else if (BOSdirlenKANDO == 1 & as.numeric(attr(BOSdirKANDOtab,"names")[1]) == 1) {
  BOSKANDOpos <- as.numeric(BOSdirKANDOtab)[1]
  BOSKANDOneg <- 0
} else if (BOSdirlenKANDO == 1 & as.numeric(attr(BOSdirKANDOtab,"names")[1]) == -1) {
  BOSKANDOpos <- 0
  BOSKANDOneg <- as.numeric(BOSdirKANDOtab)[1]
} else if (BOSdirlenKANDO == 2) {
  BOSKANDOneg <- as.numeric(BOSdirKANDOtab)[1]
  BOSKANDOpos <- as.numeric(BOSdirKANDOtab)[2]
} else {
  BOSKANDOpos <- 0
  BOSKANDOneg <- 0
}

if (BOSdirlenKURO == 0) {
  BOSKUROpos <- 0
  BOSKUROneg <- 0
} else if (BOSdirlenKURO == 1 & as.numeric(attr(BOSdirKUROtab,"names")[1]) == 1) {
  BOSKUROpos <- as.numeric(BOSdirKUROtab)[1]
  BOSKUROneg <- 0
} else if (BOSdirlenKURO == 1 & as.numeric(attr(BOSdirKUROtab,"names")[1]) == -1) {
  BOSKUROpos <- 0
  BOSKUROneg <- as.numeric(BOSdirKUROtab)[1]
} else if (BOSdirlenKURO == 2) {
  BOSKUROneg <- as.numeric(BOSdirKUROtab)[1]
  BOSKUROpos <- as.numeric(BOSdirKUROtab)[2]
} else {
  BOSKUROpos <- 0
  BOSKUROneg <- 0
}

if (BOSdirlenORGANO == 0) {
  BOSORGANOpos <- 0
  BOSORGANOneg <- 0
} else if (BOSdirlenORGANO == 1 & as.numeric(attr(BOSdirORGANOTab,"names")[1]) == 1) {
  BOSORGANOpos <- as.numeric(BOSdirORGANOTab)[1]
  BOSORGANOneg <- 0
} else if (BOSdirlenORGANO == 1 & as.numeric(attr(BOSdirORGANOTab,"names")[1]) == -1) {
  BOSORGANOpos <- 0
  BOSORGANOneg <- as.numeric(BOSdirORGANOTab)[1]
} else if (BOSdirlenORGANO == 2) {
  BOSORGANOneg <- as.numeric(BOSdirORGANOTab)[1]
  BOSORGANOpos <- as.numeric(BOSdirORGANOTab)[2]
} else {
  BOSORGANOpos <- 0
  BOSORGANOneg <- 0
}

if (BOSdirlenPODO == 0) {
  BOSPODOpos <- 0
  BOSPODOneg <- 0
} else if (BOSdirlenPODO == 1 & as.numeric(attr(BOSdirPODOtab,"names")[1]) == 1) {
  BOSPODOpos <- as.numeric(BOSdirPODOtab)[1]
  BOSPODOneg <- 0
} else if (BOSdirlenPODO == 1 & as.numeric(attr(BOSdirPODOtab,"names")[1]) == -1) {
  BOSPODOpos <- 0
  BOSPODOneg <- as.numeric(BOSdirPODOtab)[1]
} else if (BOSdirlenPODO == 2) {
  BOSPODOneg <- as.numeric(BOSdirPODOtab)[1]
  BOSPODOpos <- as.numeric(BOSdirPODOtab)[2]
} else {
  BOSPODOpos <- 0
  BOSPODOneg <- 0
}

if (BOSdirlenRUDO == 0) {
  BOSRUDOpos <- 0
  BOSRUDOneg <- 0
} else if (BOSdirlenRUDO == 1 & as.numeric(attr(BOSdirRUDOtab,"names")[1]) == 1) {
  BOSRUDOpos <- as.numeric(BOSdirRUDOtab)[1]
  BOSRUDOneg <- 0
} else if (BOSdirlenRUDO == 1 & as.numeric(attr(BOSdirRUDOtab,"names")[1]) == -1) {
  BOSRUDOpos <- 0
  BOSRUDOneg <- as.numeric(BOSdirRUDOtab)[1]
} else if (BOSdirlenRUDO == 2) {
  BOSRUDOneg <- as.numeric(BOSdirRUDOtab)[1]
  BOSRUDOpos <- as.numeric(BOSdirRUDOtab)[2]
} else {
  BOSRUDOpos <- 0
  BOSRUDOneg <- 0
}

if (BOSdirlenSODO == 0) {
  BOSSODOpos <- 0
  BOSSODOneg <- 0
} else if (BOSdirlenSODO == 1 & as.numeric(attr(BOSdirSODOtab,"names")[1]) == 1) {
  BOSSODOpos <- as.numeric(BOSdirSODOtab)[1]
  BOSSODOneg <- 0
} else if (BOSdirlenSODO == 1 & as.numeric(attr(BOSdirSODOtab,"names")[1]) == -1) {
  BOSSODOpos <- 0
  BOSSODOneg <- as.numeric(BOSdirSODOtab)[1]
} else if (BOSdirlenSODO == 2) {
  BOSSODOneg <- as.numeric(BOSdirSODOtab)[1]
  BOSSODOpos <- as.numeric(BOSdirSODOtab)[2]
} else {
  BOSSODOpos <- 0
  BOSSODOneg <- 0
}

if (BOSdirlenTENO == 0) {
  BOSTENOpos <- 0
  BOSTENOneg <- 0
} else if (BOSdirlenTENO == 1 & as.numeric(attr(BOSdirTENOTab,"names")[1]) == 1) {
  BOSTENOpos <- as.numeric(BOSdirTENOTab)[1]
  BOSTENOneg <- 0
} else if (BOSdirlenTENO == 1 & as.numeric(attr(BOSdirTENOTab,"names")[1]) == -1) {
  BOSTENOpos <- 0
  BOSTENOneg <- as.numeric(BOSdirTENOTab)[1]
} else if (BOSdirlenTENO == 2) {
  BOSTENOneg <- as.numeric(BOSdirTENOTab)[1]
  BOSTENOpos <- as.numeric(BOSdirTENOTab)[2]
} else {
  BOSTENOpos <- 0
  BOSTENOneg <- 0
}

if (BOSdirlenVERTO == 0) {
  BOSVERTOpos <- 0
  BOSVERTOneg <- 0
} else if (BOSdirlenVERTO == 1 & as.numeric(attr(BOSdirVERTOTab,"names")[1]) == 1) {
  BOSVERTOpos <- as.numeric(BOSdirVERTOTab)[1]
  BOSVERTOneg <- 0
} else if (BOSdirlenVERTO == 1 & as.numeric(attr(BOSdirVERTOTab,"names")[1]) == -1) {
  BOSVERTOpos <- 0
  BOSVERTOneg <- as.numeric(BOSdirVERTOTab)[1]
} else if (BOSdirlenVERTO == 2) {
  BOSVERTOneg <- as.numeric(BOSdirVERTOTab)[1]
  BOSVERTOpos <- as.numeric(BOSdirVERTOTab)[2]
} else {
  BOSVERTOpos <- 0
  BOSVERTOneg <- 0
}

BOSnznegdirrslts <- c(BOSCALCneg, BOSCHROMOneg, BOSDERMOneg, BOSFERROneg, BOSHYDROneg, BOSKANDOneg, BOSKUROneg, BOSORGANOneg, BOSPODOneg, BOSRUDOneg, BOSSODOneg, BOSTENOneg, BOSVERTOneg)
BOSnzposdirrslts <- c(BOSCALCpos, BOSCHROMOpos, BOSDERMOpos, BOSFERROpos, BOSHYDROpos, BOSKANDOpos, BOSKUROpos, BOSORGANOpos, BOSPODOpos, BOSRUDOpos, BOSSODOpos, BOSTENOpos, BOSVERTOpos)

BOSresults.out <- data.frame(soil=BOSsoildir.lab, nonzero=BOSnzrslts/iter, negdir=BOSnznegdirrslts/iter,
                             posdir=BOSnzposdirrslts/iter)
BOSresults.out

# ER
10^median(log10(BOSstor.mat[,"ER"]))
quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(BOSstor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(BOSstor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(BOSstor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(BOSstor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(BOSstor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)
