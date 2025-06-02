# soil Hg analysis
# determinants of Hg in soil for Australia
# Corey Bradshaw
# June 2025

# rm(list = ls())

# increase memory to max
mem.maxVSize(v = Inf)

# libraries
library(ade4) # analyse ecological/environmental data
library(adegraphics) # visualisation of ade4 objects
library(adespatial) 
library(blockCV) # blocking of spatially structured data
library(boot) # bootstrapping
library(caret) # for model training and prediction
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
library(rayshader) # for 3D visualization
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
geochem.dat$CLAY

# find duplicates
geochem.dupl.vec <- which(is.na(geochem.dat$DUPL) == F)
length(geochem.dupl.vec)

# list duplicate pairs
geodupl.pairs <- data.frame(ID1=NA, ID2=NA)
for (d in 1:length(geochem.dupl.vec)) {
  geodupl.pairs[d, ] <- (geochem.dat[geochem.dupl.vec[d], c("SITEID", "DUPL")])
}
head(geodupl.pairs)
dim(geodupl.pairs)

# average duplicates
for (p in 1:dim(geodupl.pairs)[1]) {
  l1 <- length(which(geochem.dat$SITEID == geodupl.pairs[p,1]))
  l2 <- length(which(geochem.dat$SITEID == geodupl.pairs[p,2]))

  if (l1 == 2 & l2 == 2) {
    sub1 <- which(geochem.dat$SITEID == geodupl.pairs[p,1])[1]
    sub2 <- which(geochem.dat$SITEID == geodupl.pairs[p,1])[2]
    sub3 <- which(geochem.dat$SITEID == geodupl.pairs[p,2])[1]
    sub4 <- which(geochem.dat$SITEID == geodupl.pairs[p,2])[2]
  } # end if
  
  if (l1 == 2 & l2 == 2) {
    geochem.dat[sub1, ]$Al <- mean(c(geochem.dat[sub1, ]$Al, geochem.dat[sub2, ]$Al,
                                     geochem.dat[sub3, ]$Al, geochem.dat[sub4, ]$Al), na.rm=T)
    geochem.dat[sub2, ]$Al <- geochem.dat[sub1, ]$Al; geochem.dat[sub3, ]$Al <- geochem.dat[sub1, ]$Al
    geochem.dat[sub4, ]$Al <- geochem.dat[sub1, ]$Al 
    
    geochem.dat[sub1, ]$As <- mean(c(geochem.dat[sub1, ]$As, geochem.dat[sub2, ]$As,
                                     geochem.dat[sub3, ]$As, geochem.dat[sub4, ]$As), na.rm=T)
    geochem.dat[sub2, ]$As <- geochem.dat[sub1, ]$As; geochem.dat[sub3, ]$As <- geochem.dat[sub1, ]$As
    geochem.dat[sub4, ]$As <- geochem.dat[sub1, ]$As 
    geochem.dat[sub1, ]$Ba <- mean(c(geochem.dat[sub1, ]$Ba, geochem.dat[sub2, ]$Ba,
                                     geochem.dat[sub3, ]$Ba, geochem.dat[sub4, ]$Ba), na.rm=T)
    geochem.dat[sub2, ]$Ba <- geochem.dat[sub1, ]$Ba; geochem.dat[sub3, ]$Ba <- geochem.dat[sub1, ]$Ba
    geochem.dat[sub4, ]$Ba <- geochem.dat[sub1, ]$Ba 
    geochem.dat[sub1, ]$Bi <- mean(c(geochem.dat[sub1, ]$Bi, geochem.dat[sub2, ]$Bi,
                                     geochem.dat[sub3, ]$Bi, geochem.dat[sub4, ]$Bi), na.rm=T)
    geochem.dat[sub2, ]$Bi <- geochem.dat[sub1, ]$Bi; geochem.dat[sub3, ]$Bi <- geochem.dat[sub1, ]$Bi
    geochem.dat[sub4, ]$Bi <- geochem.dat[sub1, ]$Bi 
    geochem.dat[sub1, ]$Ca <- mean(c(geochem.dat[sub1, ]$Ca, geochem.dat[sub2, ]$Ca,
                                     geochem.dat[sub3, ]$Ca, geochem.dat[sub4, ]$Ca), na.rm=T)
    geochem.dat[sub2, ]$Ca <- geochem.dat[sub1, ]$Ca; geochem.dat[sub3, ]$Ca <- geochem.dat[sub1, ]$Ca
    geochem.dat[sub4, ]$Ca <- geochem.dat[sub1, ]$Ca 
    geochem.dat[sub1, ]$CdICPMS <- mean(c(geochem.dat[sub1, ]$CdICPMS, geochem.dat[sub2, ]$CdICPMS,
                                     geochem.dat[sub3, ]$CdICPMS, geochem.dat[sub4, ]$CdICPMS), na.rm=T)
    geochem.dat[sub2, ]$CdICPMS <- geochem.dat[sub1, ]$CdICPMS; geochem.dat[sub3, ]$CdICPMS <- geochem.dat[sub1, ]$CdICPMS
    geochem.dat[sub4, ]$CdICPMS <- geochem.dat[sub1, ]$CdICPMS 
    geochem.dat[sub1, ]$CdAR <- mean(c(geochem.dat[sub1, ]$CdAR, geochem.dat[sub2, ]$CdAR,
                                     geochem.dat[sub3, ]$CdAR, geochem.dat[sub4, ]$CdAR), na.rm=T)
    geochem.dat[sub2, ]$CdAR <- geochem.dat[sub1, ]$CdAR; geochem.dat[sub3, ]$CdAR <- geochem.dat[sub1, ]$CdAR
    geochem.dat[sub4, ]$CdAR <- geochem.dat[sub1, ]$CdAR 
    geochem.dat[sub1, ]$Ce <- mean(c(geochem.dat[sub1, ]$Ce, geochem.dat[sub2, ]$Ce,
                                     geochem.dat[sub3, ]$Ce, geochem.dat[sub4, ]$Ce), na.rm=T)
    geochem.dat[sub2, ]$Ce <- geochem.dat[sub1, ]$Ce; geochem.dat[sub3, ]$Ce <- geochem.dat[sub1, ]$Ce
    geochem.dat[sub4, ]$Ce <- geochem.dat[sub1, ]$Ce 
    geochem.dat[sub1, ]$Co <- mean(c(geochem.dat[sub1, ]$Co, geochem.dat[sub2, ]$Co,
                                     geochem.dat[sub3, ]$Co, geochem.dat[sub4, ]$Co), na.rm=T)
    geochem.dat[sub2, ]$Co <- geochem.dat[sub1, ]$Co; geochem.dat[sub3, ]$Co <- geochem.dat[sub1, ]$Co
    geochem.dat[sub4, ]$Co <- geochem.dat[sub1, ]$Co 
    geochem.dat[sub1, ]$Cr <- mean(c(geochem.dat[sub1, ]$Cr, geochem.dat[sub2, ]$Cr,
                                     geochem.dat[sub3, ]$Cr, geochem.dat[sub4, ]$Cr), na.rm=T)
    geochem.dat[sub2, ]$Cr <- geochem.dat[sub1, ]$Cr; geochem.dat[sub3, ]$Cr <- geochem.dat[sub1, ]$Cr
    geochem.dat[sub4, ]$Cr <- geochem.dat[sub1, ]$Cr 
    geochem.dat[sub1, ]$Cs <- mean(c(geochem.dat[sub1, ]$Cs, geochem.dat[sub2, ]$Cs,
                                     geochem.dat[sub3, ]$Cs, geochem.dat[sub4, ]$Cs), na.rm=T)
    geochem.dat[sub2, ]$Cs <- geochem.dat[sub1, ]$Cs; geochem.dat[sub3, ]$Cs <- geochem.dat[sub1, ]$Cs
    geochem.dat[sub4, ]$Cs <- geochem.dat[sub1, ]$Cs 
    geochem.dat[sub1, ]$CuICPMS <- mean(c(geochem.dat[sub1, ]$CuICPMS, geochem.dat[sub2, ]$CuICPMS,
                                     geochem.dat[sub3, ]$CuICPMS, geochem.dat[sub4, ]$CuICPMS), na.rm=T)
    geochem.dat[sub2, ]$CuICPMS <- geochem.dat[sub1, ]$CuICPMS; geochem.dat[sub3, ]$CuICPMS <- geochem.dat[sub1, ]$CuICPMS
    geochem.dat[sub4, ]$CuICPMS <- geochem.dat[sub1, ]$CuICPMS 
    geochem.dat[sub1, ]$CuAR <- mean(c(geochem.dat[sub1, ]$CuAR, geochem.dat[sub2, ]$CuAR,
                                     geochem.dat[sub3, ]$CuAR, geochem.dat[sub4, ]$CuAR), na.rm=T)
    geochem.dat[sub2, ]$CuAR <- geochem.dat[sub1, ]$CuAR; geochem.dat[sub3, ]$CuAR <- geochem.dat[sub1, ]$CuAR
    geochem.dat[sub4, ]$CuAR <- geochem.dat[sub1, ]$CuAR 
    geochem.dat[sub1, ]$F <- mean(c(geochem.dat[sub1, ]$F, geochem.dat[sub2, ]$F,
                                     geochem.dat[sub3, ]$F, geochem.dat[sub4, ]$F), na.rm=T)
    geochem.dat[sub2, ]$F <- geochem.dat[sub1, ]$F; geochem.dat[sub3, ]$F <- geochem.dat[sub1, ]$F
    geochem.dat[sub4, ]$F <- geochem.dat[sub1, ]$F 
    geochem.dat[sub1, ]$Fe <- mean(c(geochem.dat[sub1, ]$Fe, geochem.dat[sub2, ]$Fe,
                                     geochem.dat[sub3, ]$Fe, geochem.dat[sub4, ]$Fe), na.rm=T)
    geochem.dat[sub2, ]$Fe <- geochem.dat[sub1, ]$Fe; geochem.dat[sub3, ]$Fe <- geochem.dat[sub1, ]$Fe
    geochem.dat[sub4, ]$Fe <- geochem.dat[sub1, ]$Fe 
    geochem.dat[sub1, ]$Ga <- mean(c(geochem.dat[sub1, ]$Ga, geochem.dat[sub2, ]$Ga,
                                     geochem.dat[sub3, ]$Ga, geochem.dat[sub4, ]$Ga), na.rm=T)
    geochem.dat[sub2, ]$Ga <- geochem.dat[sub1, ]$Ga; geochem.dat[sub3, ]$Ga <- geochem.dat[sub1, ]$Ga
    geochem.dat[sub4, ]$Ga <- geochem.dat[sub1, ]$Ga 
    geochem.dat[sub1, ]$Ge <- mean(c(geochem.dat[sub1, ]$Ge, geochem.dat[sub2, ]$Ge,
                                     geochem.dat[sub3, ]$Ge, geochem.dat[sub4, ]$Ge), na.rm=T)
    geochem.dat[sub2, ]$Ge <- geochem.dat[sub1, ]$Ge; geochem.dat[sub3, ]$Ge <- geochem.dat[sub1, ]$Ge
    geochem.dat[sub4, ]$Ge <- geochem.dat[sub1, ]$Ge 
    geochem.dat[sub1, ]$Hg <- mean(c(geochem.dat[sub1, ]$Hg, geochem.dat[sub2, ]$Hg,
                                     geochem.dat[sub3, ]$Hg, geochem.dat[sub4, ]$Hg), na.rm=T)
    geochem.dat[sub2, ]$Hg <- geochem.dat[sub1, ]$Hg; geochem.dat[sub3, ]$Hg <- geochem.dat[sub1, ]$Hg
    geochem.dat[sub4, ]$Hg <- geochem.dat[sub1, ]$Hg 
    geochem.dat[sub1, ]$K <- mean(c(geochem.dat[sub1, ]$La, geochem.dat[sub2, ]$La,
                                     geochem.dat[sub3, ]$La, geochem.dat[sub4, ]$La), na.rm=T)
    geochem.dat[sub2, ]$La <- geochem.dat[sub1, ]$La; geochem.dat[sub3, ]$La <- geochem.dat[sub1, ]$La
    geochem.dat[sub4, ]$La <- geochem.dat[sub1, ]$La 
    geochem.dat[sub1, ]$Li <- mean(c(geochem.dat[sub1, ]$Li, geochem.dat[sub2, ]$Li,
                                     geochem.dat[sub3, ]$Li, geochem.dat[sub4, ]$Li), na.rm=T)
    geochem.dat[sub2, ]$Li <- geochem.dat[sub1, ]$Li; geochem.dat[sub3, ]$Li <- geochem.dat[sub1, ]$Li
    geochem.dat[sub4, ]$Li <- geochem.dat[sub1, ]$Li 
    geochem.dat[sub1, ]$Mg <- mean(c(geochem.dat[sub1, ]$Mg, geochem.dat[sub2, ]$Mg,
                                     geochem.dat[sub3, ]$Mg, geochem.dat[sub4, ]$Mg), na.rm=T)
    geochem.dat[sub2, ]$Mg <- geochem.dat[sub1, ]$Mg; geochem.dat[sub3, ]$Mg <- geochem.dat[sub1, ]$Mg
    geochem.dat[sub4, ]$Mg <- geochem.dat[sub1, ]$Mg 
    geochem.dat[sub1, ]$Mn <- mean(c(geochem.dat[sub1, ]$Mn, geochem.dat[sub2, ]$Mn,
                                     geochem.dat[sub3, ]$Mn, geochem.dat[sub4, ]$Mn), na.rm=T)
    geochem.dat[sub2, ]$Mn <- geochem.dat[sub1, ]$Mn; geochem.dat[sub3, ]$Mn <- geochem.dat[sub1, ]$Mn
    geochem.dat[sub4, ]$Mn <- geochem.dat[sub1, ]$Mn 
    geochem.dat[sub1, ]$MoICPMS <- mean(c(geochem.dat[sub1, ]$MoICPMS, geochem.dat[sub2, ]$MoICPMS,
                                     geochem.dat[sub3, ]$MoICPMS, geochem.dat[sub4, ]$MoICPMS), na.rm=T)
    geochem.dat[sub2, ]$MoICPMS <- geochem.dat[sub1, ]$MoICPMS; geochem.dat[sub3, ]$MoICPMS <- geochem.dat[sub1, ]$MoICPMS
    geochem.dat[sub4, ]$MoICPMS <- geochem.dat[sub1, ]$MoICPMS 
    geochem.dat[sub1, ]$MoAR <- mean(c(geochem.dat[sub1, ]$MoAR, geochem.dat[sub2, ]$MoAR,
                                     geochem.dat[sub3, ]$MoAR, geochem.dat[sub4, ]$MoAR), na.rm=T)
    geochem.dat[sub2, ]$MoAR <- geochem.dat[sub1, ]$MoAR; geochem.dat[sub3, ]$MoAR <- geochem.dat[sub1, ]$MoAR
    geochem.dat[sub4, ]$MoAR <- geochem.dat[sub1, ]$MoAR 
    geochem.dat[sub1, ]$Na <- mean(c(geochem.dat[sub1, ]$Na, geochem.dat[sub2, ]$Na,
                                     geochem.dat[sub3, ]$Na, geochem.dat[sub4, ]$Na), na.rm=T)
    geochem.dat[sub2, ]$Na <- geochem.dat[sub1, ]$Na; geochem.dat[sub3, ]$Na <- geochem.dat[sub1, ]$Na
    geochem.dat[sub4, ]$Na <- geochem.dat[sub1, ]$Na 
    geochem.dat[sub1, ]$Nb <- mean(c(geochem.dat[sub1, ]$Nb, geochem.dat[sub2, ]$Nb,
                                     geochem.dat[sub3, ]$Nb, geochem.dat[sub4, ]$Nb), na.rm=T)
    geochem.dat[sub2, ]$Nb <- geochem.dat[sub1, ]$Nb; geochem.dat[sub3, ]$Nb <- geochem.dat[sub1, ]$Nb
    geochem.dat[sub4, ]$Nb <- geochem.dat[sub1, ]$Nb 
    geochem.dat[sub1, ]$NiCPMS <- mean(c(geochem.dat[sub1, ]$NiCPMS, geochem.dat[sub2, ]$NiCPMS,
                                     geochem.dat[sub3, ]$NiCPMS, geochem.dat[sub4, ]$NiCPMS), na.rm=T)
    geochem.dat[sub2, ]$NiCPMS <- geochem.dat[sub1, ]$NiCPMS; geochem.dat[sub3, ]$NiCPMS <- geochem.dat[sub1, ]$NiCPMS
    geochem.dat[sub4, ]$NiCPMS <- geochem.dat[sub1, ]$NiCPMS 
    geochem.dat[sub1, ]$NiAR <- mean(c(geochem.dat[sub1, ]$NiAR, geochem.dat[sub2, ]$NiAR,
                                     geochem.dat[sub3, ]$NiAR, geochem.dat[sub4, ]$NiAR), na.rm=T)
    geochem.dat[sub2, ]$NiAR <- geochem.dat[sub1, ]$NiAR; geochem.dat[sub3, ]$NiAR <- geochem.dat[sub1, ]$NiAR
    geochem.dat[sub4, ]$NiAR <- geochem.dat[sub1, ]$NiAR 
    geochem.dat[sub1, ]$P <- mean(c(geochem.dat[sub1, ]$P, geochem.dat[sub2, ]$P,
                                     geochem.dat[sub3, ]$P, geochem.dat[sub4, ]$P), na.rm=T)
    geochem.dat[sub2, ]$P <- geochem.dat[sub1, ]$P; geochem.dat[sub3, ]$P <- geochem.dat[sub1, ]$P
    geochem.dat[sub4, ]$P <- geochem.dat[sub1, ]$P 
    geochem.dat[sub1, ]$PbICPMS <- mean(c(geochem.dat[sub1, ]$PbICPMS, geochem.dat[sub2, ]$PbICPMS,
                                     geochem.dat[sub3, ]$PbICPMS, geochem.dat[sub4, ]$PbICPMS), na.rm=T)
    geochem.dat[sub2, ]$PbICPMS <- geochem.dat[sub1, ]$PbICPMS; geochem.dat[sub3, ]$PbICPMS <- geochem.dat[sub1, ]$PbICPMS
    geochem.dat[sub4, ]$PbICPMS <- geochem.dat[sub1, ]$PbICPMS 
    geochem.dat[sub1, ]$PbAR <- mean(c(geochem.dat[sub1, ]$PbAR, geochem.dat[sub2, ]$PbAR,
                                     geochem.dat[sub3, ]$PbAR, geochem.dat[sub4, ]$PbAR), na.rm=T)
    geochem.dat[sub2, ]$PbAR <- geochem.dat[sub1, ]$PbAR; geochem.dat[sub3, ]$PbAR <- geochem.dat[sub1, ]$PbAR
    geochem.dat[sub4, ]$PbAR <- geochem.dat[sub1, ]$PbAR 
    geochem.dat[sub1, ]$Pr <- mean(c(geochem.dat[sub1, ]$Pr, geochem.dat[sub2, ]$Pr,
                                     geochem.dat[sub3, ]$Pr, geochem.dat[sub4, ]$Pr), na.rm=T)
    geochem.dat[sub2, ]$Pr <- geochem.dat[sub1, ]$Pr; geochem.dat[sub3, ]$Pr <- geochem.dat[sub1, ]$Pr
    geochem.dat[sub4, ]$Pr <- geochem.dat[sub1, ]$Pr 
    geochem.dat[sub1, ]$Rb <- mean(c(geochem.dat[sub1, ]$Rb, geochem.dat[sub2, ]$Rb,
                                     geochem.dat[sub3, ]$Rb, geochem.dat[sub4, ]$Rb), na.rm=T)
    geochem.dat[sub2, ]$Rb <- geochem.dat[sub1, ]$Rb; geochem.dat[sub3, ]$Rb <- geochem.dat[sub1, ]$Rb
    geochem.dat[sub4, ]$Rb <- geochem.dat[sub1, ]$Rb 
    geochem.dat[sub1, ]$SbICPMS <- mean(c(geochem.dat[sub1, ]$SbICPMS, geochem.dat[sub2, ]$SbICPMS,
                                     geochem.dat[sub3, ]$SbICPMS, geochem.dat[sub4, ]$SbICPMS), na.rm=T)
    geochem.dat[sub2, ]$SbICPMS <- geochem.dat[sub1, ]$SbICPMS; geochem.dat[sub3, ]$SbICPMS <- geochem.dat[sub1, ]$SbICPMS
    geochem.dat[sub4, ]$SbICPMS <- geochem.dat[sub1, ]$SbICPMS 
    geochem.dat[sub1, ]$SbAR <- mean(c(geochem.dat[sub1, ]$SbAR, geochem.dat[sub2, ]$SbAR,
                                     geochem.dat[sub3, ]$SbAR, geochem.dat[sub4, ]$SbAR), na.rm=T)
    geochem.dat[sub2, ]$SbAR <- geochem.dat[sub1, ]$SbAR; geochem.dat[sub3, ]$SbAR <- geochem.dat[sub1, ]$SbAR
    geochem.dat[sub4, ]$SbAR <- geochem.dat[sub1, ]$SbAR 
    geochem.dat[sub1, ]$Sc <- mean(c(geochem.dat[sub1, ]$Sc, geochem.dat[sub2, ]$Sc,
                                     geochem.dat[sub3, ]$Sc, geochem.dat[sub4, ]$Sc), na.rm=T)
    geochem.dat[sub2, ]$Sc <- geochem.dat[sub1, ]$Sc; geochem.dat[sub3, ]$Sc <- geochem.dat[sub1, ]$Sc
    geochem.dat[sub4, ]$Sc <- geochem.dat[sub1, ]$Sc 
    geochem.dat[sub1, ]$Se <- mean(c(geochem.dat[sub1, ]$Se, geochem.dat[sub2, ]$Se,
                                     geochem.dat[sub3, ]$Se, geochem.dat[sub4, ]$Se), na.rm=T)
    geochem.dat[sub2, ]$Se <- geochem.dat[sub1, ]$Se; geochem.dat[sub3, ]$Se <- geochem.dat[sub1, ]$Se
    geochem.dat[sub4, ]$Se <- geochem.dat[sub1, ]$Se 
    geochem.dat[sub1, ]$Si <- mean(c(geochem.dat[sub1, ]$Si, geochem.dat[sub2, ]$Si,
                                     geochem.dat[sub3, ]$Si, geochem.dat[sub4, ]$Si), na.rm=T)
    geochem.dat[sub2, ]$Si <- geochem.dat[sub1, ]$Si; geochem.dat[sub3, ]$Si <- geochem.dat[sub1, ]$Si
    geochem.dat[sub4, ]$Si <- geochem.dat[sub1, ]$Si 
    geochem.dat[sub1, ]$Sm <- mean(c(geochem.dat[sub1, ]$Sm, geochem.dat[sub2, ]$Sm,
                                     geochem.dat[sub3, ]$Sm, geochem.dat[sub4, ]$Sm), na.rm=T)
    geochem.dat[sub2, ]$Sm <- geochem.dat[sub1, ]$Sm; geochem.dat[sub3, ]$Sm <- geochem.dat[sub1, ]$Sm
    geochem.dat[sub4, ]$Sm <- geochem.dat[sub1, ]$Sm 
    geochem.dat[sub1, ]$SnICPMS <- mean(c(geochem.dat[sub1, ]$SnICPMS, geochem.dat[sub2, ]$SnICPMS,
                                     geochem.dat[sub3, ]$SnICPMS, geochem.dat[sub4, ]$SnICPMS), na.rm=T)
    geochem.dat[sub2, ]$SnICPMS <- geochem.dat[sub1, ]$SnICPMS; geochem.dat[sub3, ]$SnICPMS <- geochem.dat[sub1, ]$SnICPMS
    geochem.dat[sub4, ]$SnICPMS <- geochem.dat[sub1, ]$SnICPMS 
    geochem.dat[sub1, ]$SnAR <- mean(c(geochem.dat[sub1, ]$SnAR, geochem.dat[sub2, ]$SnAR,
                                     geochem.dat[sub3, ]$SnAR, geochem.dat[sub4, ]$SnAR), na.rm=T)
    geochem.dat[sub2, ]$SnAR <- geochem.dat[sub1, ]$SnAR; geochem.dat[sub3, ]$SnAR <- geochem.dat[sub1, ]$SnAR
    geochem.dat[sub4, ]$SnAR <- geochem.dat[sub1, ]$SnAR 
    geochem.dat[sub1, ]$Sr <- mean(c(geochem.dat[sub1, ]$Sr, geochem.dat[sub2, ]$Sr,
                                     geochem.dat[sub3, ]$Sr, geochem.dat[sub4, ]$Sr), na.rm=T)
    geochem.dat[sub2, ]$Sr <- geochem.dat[sub1, ]$Sr; geochem.dat[sub3, ]$Sr <- geochem.dat[sub1, ]$Sr
    geochem.dat[sub4, ]$Sr <- geochem.dat[sub1, ]$Sr 
    geochem.dat[sub1, ]$Ta <- mean(c(geochem.dat[sub1, ]$Ta, geochem.dat[sub2, ]$Ta,
                                     geochem.dat[sub3, ]$Ta, geochem.dat[sub4, ]$Ta), na.rm=T)
    geochem.dat[sub2, ]$Ta <- geochem.dat[sub1, ]$Ta; geochem.dat[sub3, ]$Ta <- geochem.dat[sub1, ]$Ta
    geochem.dat[sub4, ]$Ta <- geochem.dat[sub1, ]$Ta 
    geochem.dat[sub1, ]$Tb <- mean(c(geochem.dat[sub1, ]$Tb, geochem.dat[sub2, ]$Tb,
                                     geochem.dat[sub3, ]$Tb, geochem.dat[sub4, ]$Tb), na.rm=T)
    geochem.dat[sub2, ]$Tb <- geochem.dat[sub1, ]$Tb; geochem.dat[sub3, ]$Tb <- geochem.dat[sub1, ]$Tb
    geochem.dat[sub4, ]$Tb <- geochem.dat[sub1, ]$Tb 
    geochem.dat[sub1, ]$Te <- mean(c(geochem.dat[sub1, ]$Te, geochem.dat[sub2, ]$Te,
                                     geochem.dat[sub3, ]$Te, geochem.dat[sub4, ]$Te), na.rm=T)
    geochem.dat[sub2, ]$Te <- geochem.dat[sub1, ]$Te; geochem.dat[sub3, ]$Te <- geochem.dat[sub1, ]$Te
    geochem.dat[sub4, ]$Te <- geochem.dat[sub1, ]$Te 
    geochem.dat[sub1, ]$Th <- mean(c(geochem.dat[sub1, ]$Th, geochem.dat[sub2, ]$Th,
                                     geochem.dat[sub3, ]$Th, geochem.dat[sub4, ]$Th), na.rm=T)
    geochem.dat[sub2, ]$Th <- geochem.dat[sub1, ]$Th; geochem.dat[sub3, ]$Th <- geochem.dat[sub1, ]$Th
    geochem.dat[sub4, ]$Th <- geochem.dat[sub1, ]$Th 
    geochem.dat[sub1, ]$Ti <- mean(c(geochem.dat[sub1, ]$Ti, geochem.dat[sub2, ]$Ti,
                                     geochem.dat[sub3, ]$Ti, geochem.dat[sub4, ]$Ti), na.rm=T)
    geochem.dat[sub2, ]$Ti <- geochem.dat[sub1, ]$Ti; geochem.dat[sub3, ]$Ti <- geochem.dat[sub1, ]$Ti
    geochem.dat[sub4, ]$Ti <- geochem.dat[sub1, ]$Ti 
    geochem.dat[sub1, ]$Tl <- mean(c(geochem.dat[sub1, ]$Tl, geochem.dat[sub2, ]$Tl,
                                     geochem.dat[sub3, ]$Tl, geochem.dat[sub4, ]$Tl), na.rm=T)
    geochem.dat[sub2, ]$Tl <- geochem.dat[sub1, ]$Tl; geochem.dat[sub3, ]$Tl <- geochem.dat[sub1, ]$Tl
    geochem.dat[sub4, ]$Tl <- geochem.dat[sub1, ]$Tl 
    geochem.dat[sub1, ]$Tm <- mean(c(geochem.dat[sub1, ]$Tm, geochem.dat[sub2, ]$Tm,
                                     geochem.dat[sub3, ]$Tm, geochem.dat[sub4, ]$Tm), na.rm=T)
    geochem.dat[sub2, ]$Tm <- geochem.dat[sub1, ]$Tm; geochem.dat[sub3, ]$Tm <- geochem.dat[sub1, ]$Tm
    geochem.dat[sub4, ]$Tm <- geochem.dat[sub1, ]$Tm 
    geochem.dat[sub1, ]$U <- mean(c(geochem.dat[sub1, ]$U, geochem.dat[sub2, ]$U,
                                     geochem.dat[sub3, ]$U, geochem.dat[sub4, ]$U), na.rm=T)
    geochem.dat[sub2, ]$U <- geochem.dat[sub1, ]$U; geochem.dat[sub3, ]$U <- geochem.dat[sub1, ]$U
    geochem.dat[sub4, ]$U <- geochem.dat[sub1, ]$U 
    geochem.dat[sub1, ]$V <- mean(c(geochem.dat[sub1, ]$V, geochem.dat[sub2, ]$V,
                                     geochem.dat[sub3, ]$V, geochem.dat[sub4, ]$V), na.rm=T)
    geochem.dat[sub2, ]$V <- geochem.dat[sub1, ]$V; geochem.dat[sub3, ]$V <- geochem.dat[sub1, ]$V
    geochem.dat[sub4, ]$V <- geochem.dat[sub1, ]$V 
    geochem.dat[sub1, ]$W <- mean(c(geochem.dat[sub1, ]$W, geochem.dat[sub2, ]$W,
                                     geochem.dat[sub3, ]$W, geochem.dat[sub4, ]$W), na.rm=T)
    geochem.dat[sub2, ]$W <- geochem.dat[sub1, ]$W; geochem.dat[sub3, ]$W <- geochem.dat[sub1, ]$W
    geochem.dat[sub4, ]$W <- geochem.dat[sub1, ]$W 
    geochem.dat[sub1, ]$Y <- mean(c(geochem.dat[sub1, ]$Y, geochem.dat[sub2, ]$Y,
                                     geochem.dat[sub3, ]$Y, geochem.dat[sub4, ]$Y), na.rm=T)
    geochem.dat[sub2, ]$Y <- geochem.dat[sub1, ]$Y; geochem.dat[sub3, ]$Y <- geochem.dat[sub1, ]$Y
    geochem.dat[sub4, ]$Y <- geochem.dat[sub1, ]$Y 
    geochem.dat[sub1, ]$Yb <- mean(c(geochem.dat[sub1, ]$Yb, geochem.dat[sub2, ]$Yb,
                                     geochem.dat[sub3, ]$Yb, geochem.dat[sub4, ]$Yb), na.rm=T)
    geochem.dat[sub2, ]$Yb <- geochem.dat[sub1, ]$Yb; geochem.dat[sub3, ]$Yb <- geochem.dat[sub1, ]$Yb
    geochem.dat[sub4, ]$Yb <- geochem.dat[sub1, ]$Yb 
    geochem.dat[sub1, ]$ZnICPMS <- mean(c(geochem.dat[sub1, ]$ZnICPMS, geochem.dat[sub2, ]$ZnICPMS,
                                     geochem.dat[sub3, ]$ZnICPMS, geochem.dat[sub4, ]$ZnICPMS), na.rm=T)
    geochem.dat[sub2, ]$ZnICPMS <- geochem.dat[sub1, ]$ZnICPMS; geochem.dat[sub3, ]$ZnICPMS <- geochem.dat[sub1, ]$ZnICPMS
    geochem.dat[sub4, ]$ZnICPMS <- geochem.dat[sub1, ]$ZnICPMS 
    geochem.dat[sub1, ]$ZnAR <- mean(c(geochem.dat[sub1, ]$ZnAR, geochem.dat[sub2, ]$ZnAR,
                                     geochem.dat[sub3, ]$ZnAR, geochem.dat[sub4, ]$ZnAR), na.rm=T)
    geochem.dat[sub2, ]$ZnAR <- geochem.dat[sub1, ]$ZnAR; geochem.dat[sub3, ]$ZnAR <- geochem.dat[sub1, ]$ZnAR
    geochem.dat[sub4, ]$ZnAR <- geochem.dat[sub1, ]$ZnAR 
    geochem.dat[sub1, ]$FIELDpH <- mean(c(geochem.dat[sub1, ]$FIELDpH, geochem.dat[sub2, ]$FIELDpH,
                                     geochem.dat[sub3, ]$FIELDpH, geochem.dat[sub4, ]$FIELDpH), na.rm=T)
    geochem.dat[sub2, ]$FIELDpH <- geochem.dat[sub1, ]$FIELDpH; geochem.dat[sub3, ]$FIELDpH <- geochem.dat[sub1, ]$FIELDpH
    geochem.dat[sub4, ]$FIELDpH <- geochem.dat[sub1, ]$FIELDpH 
    geochem.dat[sub1, ]$pH15 <- mean(c(geochem.dat[sub1, ]$pH15, geochem.dat[sub2, ]$pH15,
                                     geochem.dat[sub3, ]$pH15, geochem.dat[sub4, ]$pH15), na.rm=T)
    geochem.dat[sub2, ]$pH15 <- geochem.dat[sub1, ]$pH15; geochem.dat[sub3, ]$pH15 <- geochem.dat[sub1, ]$pH15
    geochem.dat[sub4, ]$pH15 <- geochem.dat[sub1, ]$pH15 
    geochem.dat[sub1, ]$EC15 <- mean(c(geochem.dat[sub1, ]$EC15, geochem.dat[sub2, ]$EC15,
                                     geochem.dat[sub3, ]$EC15, geochem.dat[sub4, ]$EC15), na.rm=T)
    geochem.dat[sub2, ]$EC15 <- geochem.dat[sub1, ]$EC15; geochem.dat[sub3, ]$EC15 <- geochem.dat[sub1, ]$EC15
    geochem.dat[sub4, ]$EC15 <- geochem.dat[sub1, ]$EC15 
    geochem.dat[sub1, ]$SAND <- mean(c(geochem.dat[sub1, ]$SAND, geochem.dat[sub2, ]$SAND,
                                     geochem.dat[sub3, ]$SAND, geochem.dat[sub4, ]$SAND), na.rm=T)
    geochem.dat[sub2, ]$SAND <- geochem.dat[sub1, ]$SAND; geochem.dat[sub3, ]$SAND <- geochem.dat[sub1, ]$SAND
    geochem.dat[sub4, ]$SAND <- geochem.dat[sub1, ]$SAND 
    geochem.dat[sub1, ]$SILT <- mean(c(geochem.dat[sub1, ]$SILT, geochem.dat[sub2, ]$SILT,
                                     geochem.dat[sub3, ]$SILT, geochem.dat[sub4, ]$SILT), na.rm=T)
    geochem.dat[sub2, ]$SILT <- geochem.dat[sub1, ]$SILT; geochem.dat[sub3, ]$SILT <- geochem.dat[sub1, ]$SILT
    geochem.dat[sub4, ]$SILT <- geochem.dat[sub1, ]$SILT 
    geochem.dat[sub1, ]$CLAY <- mean(c(geochem.dat[sub1, ]$CLAY, geochem.dat[sub2, ]$CLAY,
                                     geochem.dat[sub3, ]$CLAY, geochem.dat[sub4, ]$CLAY), na.rm=T)
    geochem.dat[sub2, ]$CLAY <- geochem.dat[sub1, ]$CLAY; geochem.dat[sub3, ]$CLAY <- geochem.dat[sub1, ]$CLAY
    geochem.dat[sub4, ]$CLAY <- geochem.dat[sub1, ]$CLAY 
    geochem.dat[sub1, ]$LOI <- mean(c(geochem.dat[sub1, ]$LOI, geochem.dat[sub2, ]$LOI,
                                     geochem.dat[sub3, ]$LOI, geochem.dat[sub4, ]$LOI), na.rm=T)
    geochem.dat[sub2, ]$LOI <- geochem.dat[sub1, ]$LOI; geochem.dat[sub3, ]$LOI <- geochem.dat[sub1, ]$LOI
    geochem.dat[sub4, ]$LOI <- geochem.dat[sub1, ]$LOI 
    } # end if
  
} # end p loop

# remove duplicates
result <- geodupl.pairs %>%
  mutate(min_id = pmin(ID1, ID2),
         max_id = pmax(ID1, ID2)) %>%
  distinct(min_id, max_id)

removes1 <- which((geodupl.pairs$ID1 %in% result$min_id) == F)
ID.rem1 <- geodupl.pairs$ID1[removes1]
length(ID.rem1)

removes2 <- which((geodupl.pairs$ID2 %in% result$max_id) == F)
ID.rem2 <- geodupl.pairs$ID2[removes2]
length(ID.rem2)

geochem.dat2 <- geochem.dat[which(geochem.dat$SITEID %in% ID.rem1 == F), ]
dim(geochem.dat2)
dim(geochem.dat)

geochem.dat3 <- geochem.dat2[which(geochem.dat2$SITEID %in% ID.rem2 == F), ]
head(geochem.dat3)
dim(geochem.dat3)
dim(geochem.dat)


# check no more duplicates
geochem.dupl2.vec <- which(is.na(geochem.dat3$DUPL) == F)
length(geochem.dupl2.vec)

## field dataset
field.dat <- read.csv("field.csv", header = T)
head(field.dat)

# find duplicates
field.dupl.vec <- which(is.na(field.dat$DUPL) == F)
length(field.dupl.vec)

# list duplicate pairs
fielddupl.pairs <- data.frame(ID1=NA, ID2=NA)
for (d in 1:length(field.dupl.vec)) {
  fielddupl.pairs[d, ] <- (field.dat[field.dupl.vec[d], c("SITEID", "DUPL")])
}
head(fielddupl.pairs)
dim(fielddupl.pairs)

# average duplicates
for (p in 1:dim(fielddupl.pairs)[1]) {
  sub1 <- which(field.dat$SITEID == fielddupl.pairs[p,1])
  sub2 <- which(field.dat$SITEID == fielddupl.pairs[p,2])
  
  field.dat[sub1, ]$TOSTOPD <- mean(c(field.dat[sub1, ]$TOSTOPD, field.dat[sub2, ]$TOSTOPD), na.rm=T)
  field.dat[sub2, ]$TOSTOPD <- mean(c(field.dat[sub1, ]$TOSTOPD, field.dat[sub2, ]$TOSTOPD), na.rm=T)
  field.dat[sub1, ]$TOSBD <- mean(c(field.dat[sub1, ]$TOSBD, field.dat[sub2, ]$TOSBD), na.rm=T)
  field.dat[sub2, ]$TOSBD <- mean(c(field.dat[sub1, ]$TOSBD, field.dat[sub2, ]$TOSBD), na.rm=T)
  field.dat[sub1, ]$TOSpH <- mean(c(field.dat[sub1, ]$TOSpH, field.dat[sub2, ]$TOSpH), na.rm=T)
  field.dat[sub2, ]$TOSpH <- mean(c(field.dat[sub1, ]$TOSpH, field.dat[sub2, ]$TOSpH), na.rm=T)
  field.dat[sub1, ]$BOSTOPD <- mean(c(field.dat[sub1, ]$BOSTOPD, field.dat[sub2, ]$BOSTOPD), na.rm=T)
  field.dat[sub2, ]$BOSTOPD <- mean(c(field.dat[sub1, ]$BOSTOPD, field.dat[sub2, ]$BOSTOPD), na.rm=T)
  field.dat[sub1, ]$BOSBD <- mean(c(field.dat[sub1, ]$BOSBD, field.dat[sub2, ]$BOSBD), na.rm=T)
  field.dat[sub2, ]$BOSBD <- mean(c(field.dat[sub1, ]$BOSBD, field.dat[sub2, ]$BOSBD), na.rm=T)
  field.dat[sub1, ]$BOSpH <- mean(c(field.dat[sub1, ]$BOSpH, field.dat[sub2, ]$BOSpH), na.rm=T)
  field.dat[sub2, ]$BOSpH <- mean(c(field.dat[sub1, ]$BOSpH, field.dat[sub2, ]$BOSpH), na.rm=T)
}

# remove duplicates
result <- fielddupl.pairs %>%
  mutate(min_id = pmin(ID1, ID2),
         max_id = pmax(ID1, ID2)) %>%
         distinct(min_id, max_id)

removes <- which((fielddupl.pairs$ID1 %in% result$min_id) == F)
ID.rem <- fielddupl.pairs$ID1[removes]
length(ID.rem)

field.dat2 <- field.dat[which(field.dat$SITEID %in% ID.rem == F), ]
dim(field.dat2)
dim(field.dat)

removes2 <- which((fielddupl.pairs$ID2 %in% result$max_id) == F)
ID.rem2 <- fielddupl.pairs$ID2[removes2]
length(ID.rem2)

field.dat3 <- field.dat2[which(field.dat2$SITEID %in% ID.rem2 == F), ]
dim(field.dat3)
dim(field.dat)

# check no more duplicates
field.dupl2.vec <- which(is.na(field.dat3$DUPL) == F)
length(field.dupl2.vec)


## Hg dataset
hg.dat <- read.csv("hgTSID.csv", header = T)
head(hg.dat)


## grain size dataset
gs.dat <- read.csv("gs.csv", header = T)
head(gs.dat)

# grain size categories
gs.cat <- read.csv("gscats.csv", header = T)
head(gs.cat)

# create weighted mean grain size column
dim(gs.dat)
head(gs.dat[10:dim(gs.dat)[2]])
apply(gs.dat[,10:dim(gs.dat)[2]], 1, sum, na.rm=T)
gs.dat[1, 10:dim(gs.dat)[2]]
sum(gs.dat[1, 10:dim(gs.dat)[2]])
sum(gs.cat$mn * gs.dat[1, 10:dim(gs.dat)[2]]/100)

gs.dat$gs.wmn <- NA
for (i in 1:dim(gs.dat)[1]) {
  gs.dat$gs.wmn[i] <- sum(gs.cat$mn * as.numeric(gs.dat[i, 10:(dim(gs.dat)[2]-1)])/100, na.rm=T)
}
gs.dat$gs.wmn
head(gs.dat[1:3,])

# remove rows with now mean grain size elements
gs.dat2 <- subset(gs.dat, GS=='bulk')
dim(gs.dat2)
head(gs.dat2)
gs.dat2$gs.wmn

# merge datasets
datmrg1 <- merge(geochem.dat3, field.dat3, by="SITEID", all.y=F, all.x=F)
head(datmrg1)

datmrg2 <- merge(hg.dat, datmrg1, by="SITEID", all.y=F, all.x=F)
head(datmrg2)

datmrg3 <- merge(gs.dat2, datmrg2, by="SITEID", all.y=F, all.x=F)
head(datmrg3)
dim(datmrg3)

# remove duplicate rows
datmrg <- distinct(datmrg3)
head(datmrg)
dim(datmrg)

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

## Hg sample location coordinates
Hgpts <- vect(cbind(datmrg$LON.x, datmrg$LAT.x), crs="+proj=longlat")
terra::plot(Hgpts)

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
datmrg$posixdate <- as.POSIXct(datmrg$DATE.x, format="%d.%m.%Y")

## overlay data onto WWF ecoregions
# WWF ecoregion polygon # download shapefile from:
# https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
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

# plot (warning: can take long time to plot)
terra::plot(WWFecoregions, 'BIOME')

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


## radiometric data https://portal.ga.gov.au/persona/gadds
# KThU
KThU <- rast("radmap_v4_2019_filtered_ML_KThU_RGB_24bit.tif")
KThU.1 <- subset(KThU, 1) # r
KThU.2 <- subset(KThU, 2) # g
KThU.3 <- subset(KThU, 3) # b
KThUcmb <- (KThU.1 + KThU.2 + KThU.3)/3
terra::plot(KThUcmb)

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

# extract
soilN.Hgpts <- terra::extract(soilN, Hgpts)
head(soilN.Hgpts)

# add to data
datmrg$soilN <- soilN.Hgpts$focal_mean

## soil phosphorus
# https://data.csiro.au/collection/csiro:61526?_st=browse&_str=2&_si=1&browseType=kw&browseValue=total%20soil%20nitrogen
soilP <- rast("PTO_000_005_EV_N_P_AU_NAT_C_20231101.tif")
terra::plot(soilP)

# extract
soilP.Hgpts <- terra::extract(soilP, Hgpts)
head(soilP.Hgpts)

# add to data
datmrg$soilP <- soilP.Hgpts$focal_mean

# soil pH (extracted)
soilpH <- rast("pHc_000_005_EV_N_P_AU_NAT_C_20140801.tif")
terra::plot(soilpH)

# extract
soilpH.Hgpts <- terra::extract(soilpH, Hgpts)
head(soilpH.Hgpts)

# add to data
datmrg$soilpH <- soilpH.Hgpts$pHc_000_005_EV_N_P_AU_NAT_C_20140801

## rainfall
# https://ausenv.tern.org.au/aer/how-to-use-australias-environment-data-explorer/aer/australias-environment/index.html
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
#writeRaster(soilH20.mn.rst, "soilH20Mn.tif")

# extract
soilH20.Hgpts <- terra::extract(soilH20.mn.rst, Hgpts)
head(soilH20.Hgpts)

# add to data
datmrg$soilH20 <- soilH20.Hgpts$mean

## leaf area index
# https://ausenv.tern.org.au/aer/how-to-use-australias-environment-data-explorer/aer/australias-environment/index.html
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

# extract
lai.Hgpts <- terra::extract(lai.mn.rst, Hgpts)
head(lai.Hgpts)

# add to data
datmrg$lai <- lai.Hgpts$mean

## vegetation carbon uptake (GPP)
# amount of carbon taken up by the vegetation through photosynthesis,
# as estimated by the OzWALD model-data fusion system
# https://ausenv.tern.org.au/aer/how-to-use-australias-environment-data-explorer/aer/australias-environment/index.html
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
saveRDS(soil, "soil.rds")

# soil attributes
soil.key <- read.csv("asclut.csv")
head(soil.key)

soils <- merge(soil, soil.key, by="MAP_UNIT")
soils
attr(table(soils$TYPE), 'names')
saveRDS(soils, "soils.rds")

# export
#writeVector(soils, filename = "soilsMain.shp", overwrite = TRUE)

# extract
soils.Hgpts <- terra::extract(soils, Hgpts)
head(soils.Hgpts)

# add to data
datmrg$soils <- soils.Hgpts$TYPE

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
ferricpc4 <- rast("FERRICPC4mrg.tif")
crs(ferricpc4)
crs(ferricpc4) <- "epsg:3577"
ext(ferricpc4)
terra::plot(ferricpc4)
saveRDS(ferricpc4, "ferricpc4.rds")

ferricpc2 <- rast("FERRICPC2mrg.tif")
crs(ferricpc2)
crs(ferricpc2) <- "epsg:3577"
ext(ferricpc2)
terra::plot(ferricpc2)
saveRDS(ferricpc2, "ferricpc2.rds")

aus.terr <- vect("aus.shp")
crs(aus.terr)
ext(aus.terr)

# project aus.terr to "epsg:3577"
aus.terr3577 <- project(aus.terr, "epsg:3577")
terra::plot(aus.terr3577)
saveRDS(aus.terr3577, "austerr3577.rds")

# find projection boundaries for aus.terr9473
ext3577 <- ext(aus.terr3577)

# clip to aus.terr
ferricpc4.clip <- crop(ferricpc4, aus.terr3577)
ferricpc2.clip <- crop(ferricpc2, aus.terr3577)

# project ferricpc4.clip to "epsg:9473"
ferricpc4clip9473 <- project(ferricpc4.clip, "epsg:9473")
terra::plot(ferricpc4clip9473)
saveRDS(ferricpc4clip9473, "ferricpc4clip9473.rds")
ferricpc2clip9473 <- project(ferricpc2.clip, "epsg:9473")
terra::plot(ferricpc2clip9473)
saveRDS(ferricpc2clip9473, "ferricpc2clip9473.rds")

# write raster
writeRaster(ferricpc4clip9473, "ferricpc4clip9473.tif")
writeRaster(ferricpc2clip9473, "ferricpc2clip9473.tif")
ferricpc4clip9473 <- rast("ferricpc4clip9473.tif")
ferricpc2clip9473 <- rast("ferricpc2clip9473.tif")

# project back to WGS84
ferricpc4clip4326 <- project(ferricpc4clip9473, "epsg:4326")
ferricpc2clip4326 <- project(ferricpc2clip9473, "epsg:4326")
ferric2.rst <- ferricpc2clip4326
ferric4.rst <- ferricpc4clip4326
ferric2.rsmp <- resample(ferric2.rst, lai.rst)
ferric4.rsmp <- resample(ferric4.rst, lai.rst)
plot(ferric2.rsmp)
plot(ferric4.rsmp)
crs(ferric2.rsmp) <- crsUse
crs(ferric4.rsmp) <- crsUse

# extract
ferricpc4.Hgpts <- terra::extract(ferricpc4clip4326, Hgpts)
head(ferricpc4.Hgpts)
ferricpc2.Hgpts <- terra::extract(ferricpc2clip4326, Hgpts)
head(ferricpc2.Hgpts)

# add to data
datmrg$ferricpc4 <- ferricpc4.Hgpts$FERRICPC4mrg
datmrg$ferricpc2 <- ferricpc2.Hgpts$FERRICPC2mrg

## aluminium oxide
## https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/148587
## doi:10.26186/148587
Alox <- rast("Aluminium_oxide_prediction_median.tif")  # need to project original GeoTIFF into 9473 projection first
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
Feox9473 <- rast("Feox9473.tif")  # need to project original GeoTIFF into 9473 projection first
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
Pox9473 <- rast("Pox9473.tif") # need to project original GeoTIFF into 9473 projection first
terra::plot(Pox9473,col=rev(map.pal("plasma",200)))
Pox4326 <- project(Pox9473, crs(Hgpts))
Pox.rst <- Pox4326
Pox.rsmp <- resample(Pox.rst, lai.rst)
plot(Pox.rsmp)
crs(Pox.rsmp) <- crsUse

# extract
Pox.Hgpts <- terra::extract(Pox4326, Hgpts)
head(Pox.Hgpts)

# add to data
datmrg$Pox <- Pox.Hgpts$Pox9473


## merged file elements
head(datmrg)

## separate bottom from top soil
# bottom
BOS.dat <- datmrg[datmrg$DEPTH.x == "BOS",]
head(BOS.dat)
dim(BOS.dat)

# remove TOS fields from BOS.dat
str(BOS.dat)
BOS.collabs <- colnames(BOS.dat)
BOS.dat <- BOS.dat[, -(grep("TOS", BOS.collabs))]
head(BOS.dat)

# top
TOS.dat <- datmrg[datmrg$DEPTH.x == "TOS",]
head(TOS.dat)
dim(TOS.dat)

# remove BOS fields from TOS.dat
str(TOS.dat)
TOS.collabs <- colnames(TOS.dat)
TOS.dat <- TOS.dat[, -(grep("BOS", TOS.collabs))]
head(TOS.dat)

## Hg columns
head(TOS.dat[,grep("Hg", colnames(TOS.dat))])

## time plots
# Hg vs. time
plot(datmrg$posixdate, log10(datmrg$HgCOMP), xlab="date", ylab="log10 [Hg]", pch=19, col="blue")

## bivariate plots
# Hg vs. mean grain size
plot(datmrg$gs.wmn, log10(datmrg$HgCOMP), xlab="mean grain size", ylab="log10 [Hg]", pch=19, col="blue")
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

# Hg vs. pH (1:5)
plot(TOS.dat$TOSpH, log10(TOS.dat$HgCOMP), xlab="TOS pH", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$TOSpH), col="red", lwd=2, lty=2)
plot(BOS.dat$BOSpH, log10(BOS.dat$HgCOMP), xlab="BOS pH", ylab="log10 BOS [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$BOSpH), col="red", lwd=2, lty=2)

# soil pH (extracted) vs. Hg
plot(TOS.dat$soilpH, log10(TOS.dat$HgCOMP), xlab="TOS soil pH", ylab="log10 TOS [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$soilpH), col="red", lwd=2, lty=2)
plot(BOS.dat$soilpH, log10(BOS.dat$HgCOMP), xlab="BOS soil pH", ylab="log10 BOS [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ BOS.dat$soilpH), col="red", lwd=2, lty=2)

# compare pH measures
plot(TOS.dat$soilpH, TOS.dat$TOSpH, xlab="extracted pH", ylab="field pH", pch=19, col="blue")
abline(lm(TOS.dat$TOSpH ~ TOS.dat$soilpH), col="red", lwd=2, lty=2)
plot(BOS.dat$soilpH, BOS.dat$BOSpH, xlab="extracted pH", ylab="field pH", pch=19, col="blue")
abline(lm(BOS.dat$BOSpH ~ BOS.dat$soilpH), col="red", lwd=2, lty=2)

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
plot(log10(TOS.dat$ZnICPMS), log10(TOS.dat$HgCOMP), xlab="Zn", ylab="log10 [Hg]", pch=19, col="blue")
plot(log10(BOS.dat$ZnICPMS), log10(BOS.dat$HgCOMP), xlab="Zn", ylab="log10 [Hg]", pch=19, col="blue")

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

# ferric PC4 barest earth soil vs. Hg
plot((TOS.dat$ferricpc4), log10(TOS.dat$HgCOMP), xlab="ferric PC4", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$ferricpc4), col="red", lwd=2, lty=2)
hist(TOS.dat$ferricpc4)

# ferric PC2 barest earth soil vs. Hg
plot((TOS.dat$ferricpc2), log10(TOS.dat$HgCOMP), xlab="ferric PC4", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ TOS.dat$ferricpc2), col="red", lwd=2, lty=2)
hist(TOS.dat$ferricpc2)

# latitude vs. Hg
plot(datmrg$LAT.x, log10(datmrg$HgCOMP), xlab="latitude", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(datmrg$HgCOMP) ~ datmrg$LAT.x), col="red", lwd=2, lty=2)

# longitude vs. Hg
plot(datmrg$LON.x, log10(datmrg$HgCOMP), xlab="longitude", ylab="log10 [Hg]", pch=19, col="blue")
abline(lm(log10(datmrg$HgCOMP) ~ datmrg$LON.x), col="red", lwd=2, lty=2)

# rainfall vs. Hg
plot(log10(TOS.dat$rain), log10(TOS.dat$HgCOMP), xlab="annual rainfall (mm)", ylab="TOS log10 [Hg]", pch=19, col="blue")
abline(lm(log10(TOS.dat$HgCOMP) ~ log10(TOS.dat$rain), na.rm=T), col="red", lwd=2, lty=2)
plot(log10(BOS.dat$rain), log10(BOS.dat$HgCOMP), xlab="annual rainfall (mm)", ylab="BOS log10 [Hg]", pch=19, col="blue")
abline(lm(log10(BOS.dat$HgCOMP) ~ log10(BOS.dat$rain), na.rm=T), col="red", lwd=2, lty=2)

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
TOS.HgXstate <- sort(xtabs(TOS.dat$HgCOMP ~ TOS.dat$STATE.y) / table(TOS.dat$STATE.y))
barplot(TOS.HgXstate, xlab="state", ylab="mean TOS [Hg]", col="blue")
BOS.HgXstate <- sort(xtabs(BOS.dat$HgCOMP ~ BOS.dat$STATE.y) / table(BOS.dat$STATE.y))
barplot(BOS.HgXstate, xlab="state", ylab="mean BOS [Hg]", col="blue")

HgXstate.stats <- TOS.dat %>%
    group_by(STATE.x) %>%
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

HgXstate.sort <- HgXstate.stats %>% arrange(mean)
HgXstate.sort
ggplot(HgXstate.sort, aes(x = reorder(STATE.x, mean), y = mean)) +
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

HgXlucat.stats <- TOS.dat %>%
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

HgXlucat.sort <- HgXlucat.stats %>% arrange(mean)
HgXlucat.sort
ggplot(HgXlucat.sort, aes(x = reorder(lucat, mean), y = mean)) +
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

HgXlanduse.stats <- TOS.dat %>%
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

HgXlanduse.sort <- HgXlanduse.stats %>% arrange(mean)
HgXlanduse.sort
ggplot(HgXlanduse.sort, aes(x = reorder(landuse, mean), y = mean)) +
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

HgXbtype.stats <- TOS.dat %>%
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

HgXbtype.sort <- HgXbtype.stats %>% arrange(mean)
HgXbtype.sort
ggplot(HgXbtype.sort, aes(x = reorder(btype, mean), y = mean)) +
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


# Hg by geology
# lithology
TOS.HgXlith <- sort(xtabs(TOS.dat$HgCOMP ~ TOS.dat$lithgrp) / table(TOS.dat$lithgrp))
barplot(TOS.HgXlith, xlab="lithology", ylab="mean TOS [Hg]", col="blue")
BOS.HgXlith <- sort(xtabs(BOS.dat$HgCOMP ~ BOS.dat$lithgrp) / table(BOS.dat$lithgrp))
barplot(BOS.HgXlith, xlab="lithology", ylab="mean BOS [Hg]", col="blue")

HgXlith.stats <- TOS.dat %>%
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

HgXlith.sort <- na.omit(HgXlith.stats %>% arrange(mean))
HgXlith.sort
ggplot(HgXlith.sort, aes(x = reorder(lithgrp, mean), y = mean)) +
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

HgXsoil.stats <- TOS.dat %>%
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
HgXsoil.sort <- na.omit(HgXsoil.stats %>% arrange(mean))
HgXsoil.sort

soils.plot <- ggplot(HgXsoil.sort, aes(x = reorder(soils, mean), y = mean)) +
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



## create dataset for testing relationships
colnames(TOS.dat)
TOS.test1 <- data.frame(SITEID=TOS.dat$SITEID, LAT=TOS.dat$LAT.x, LON=TOS.dat$LON.x, lHg=log10(TOS.dat$HgCOMP), lgs=log10(TOS.dat$gs.wmn),
                       lgtclay=logit(TOS.dat$CLAY/100), lgtsoilclay=logit(TOS.dat$soilclay/100),
                       lgtsilt=logit(TOS.dat$SILT/100), lgtsoilsilt=logit(TOS.dat$soilsilt/100),  
                       pH=TOS.dat$TOSpH, pH15=TOS.dat$pH15, soilpH=TOS.dat$soilpH, lelecond=log10(TOS.dat$EC15), lLOI=log10(TOS.dat$LOI),
                       lrain=log10(TOS.dat$rain), prescott=TOS.dat$prescott, soilH20=TOS.dat$soilH20,
                       llai=log10(TOS.dat$lai), lgpp=log10(TOS.dat$gpp),
                       soillN=log10(TOS.dat$soilN), soillP=log10(TOS.dat$soilP),
                       lAl=log10(TOS.dat$Al), lFe=log10(TOS.dat$Fe), lMn=log10(TOS.dat$Mn), lCu=log10(TOS.dat$CuICPMS),
                       lZn=log10(TOS.dat$ZnICPMS), lPb=log10(TOS.dat$PbICPMS), lSb=log10(TOS.dat$SbICPMS),
                       lNi=log10(TOS.dat$NiCPMS), lV=log10(TOS.dat$V), lKThU=log10(TOS.dat$KThU),
                       lucat=TOS.dat$lucat, biome=TOS.dat$biome, geol=TOS.dat$lithgrp, soils=TOS.dat$soils,
                       ferric4=TOS.dat$ferricpc4, ferric2=TOS.dat$ferricpc2, Alox=TOS.dat$Alox, 
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
length(which(is.na(TOS.test$lgtsoilclay)==T)) # 28 missing
length(which(is.na(TOS.test$lgtsilt)==T)) # good
length(which(is.na(TOS.test$lgtsoilsilt)==T)) # 28 missing
length(which(is.na(TOS.test$pH15)==T)) # good
length(which(is.na(TOS.test$soilpH)==T)) # 12 missing
length(which(is.na(TOS.test$pH)==T)) # good
length(which(is.na(TOS.test$lelecond)==T)) # good
length(which(is.na(TOS.test$lLOI)==T)) # 9 missing
length(which(is.na(TOS.test$lrain)==T)) # 54 missing
length(which(is.na(TOS.test$prescott)==T)) # good
length(which(is.na(TOS.test$soilH20)==T)) # 4 missing
length(which(is.na(TOS.test$llai)==T)) # 4 missing
length(which(is.na(TOS.test$lgpp)==T)) # 4 missing
length(which(is.na(TOS.test$soillN)==T)) # 20 missing
length(which(is.na(TOS.test$soillP)==T)) # 20 missing
length(which(is.na(TOS.test$lAl)==T)) # 9 missing
length(which(is.na(TOS.test$lFe)==T)) # 9 missing
length(which(is.na(TOS.test$lMn)==T)) # 9 missing
length(which(is.na(TOS.test$lCu)==T)) # 9 missing
length(which(is.na(TOS.test$lZn)==T)) # 9 missing
length(which(is.na(TOS.test$lPb)==T)) # 9 missing
length(which(is.na(TOS.test$lSb)==T)) # 9 missing
length(which(is.na(TOS.test$lNi)==T)) # 9 missing
length(which(is.na(TOS.test$lV)==T)) # 9 missing
length(which(is.na(TOS.test$lKThU)==T)) # good
length(which(is.na(TOS.test$ferric4)==T)) # good
length(which(is.na(TOS.test$ferric2)==T)) # good
length(which(is.na(TOS.test$Alox)==T)) # good
length(which(is.na(TOS.test$Feox)==T)) # good
length(which(is.na(TOS.test$Pox)==T)) # good

## any infinites?
length(which(is.infinite(TOS.test$lHg)==T))
length(which(is.infinite(TOS.test$lgs)==T))
length(which(is.infinite(TOS.test$pH15)==T))
length(which(is.infinite(TOS.test$lgtclay)==T))
length(which(is.infinite(TOS.test$lgtsilt)==T)) # 3 missing
length(which(is.infinite(TOS.test$lelecond)==T))
length(which(is.infinite(TOS.test$lLOI)==T)) # 1 missing
length(which(is.infinite(TOS.test$lrain)==T))
length(which(is.infinite(TOS.test$prescott)==T))
length(which(is.infinite(TOS.test$soilH20)==T))
length(which(is.infinite(TOS.test$llai)==T)) # 4 missing
length(which(is.infinite(TOS.test$lgpp)==T)) # 4 missing
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
length(which(is.infinite(TOS.test$ferric4)==T))
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
hist(TOS.test$ferric4, xlab="ferric PC4", main="")
hist(TOS.test$ferric2, xlab="ferric PC2", main="")
hist(TOS.test$Alox, xlab="Al oxide", main="")
hist(TOS.test$Feox, xlab="Fe oxide", main="")
hist(TOS.test$Pox, xlab="P oxide", main="")

dim(TOS.test)

# replace infinites with NA
TOS.test$lgtclay[(which(is.infinite(TOS.test$lgtclay)==T))] <- NA
TOS.test$lgtsilt[(which(is.infinite(TOS.test$lgtsilt)==T))] <- NA
TOS.test$lLOI[(which(is.infinite(TOS.test$lLOI)==T))] <- NA
TOS.test$llai[(which(is.infinite(TOS.test$llai)==T))] <- NA
TOS.test$lgpp[(which(is.infinite(TOS.test$lgpp)==T))] <- NA

## split by biome & re-examine vif
table(TOS.test$biome)
TOS.test.tempfor <- subset(TOS.test, biome=="TBMF")
TOS.test.tempgrass <- subset(TOS.test, biome=="TGSS")
TOS.test.Medfor <- subset(TOS.test, biome=="MFW")
TOS.test.tropfor <- subset(TOS.test, biome=="TSMBF")
TOS.test.tropgrass <- subset(TOS.test, biome=="TSGSS")
TOS.test.dryxeric <- subset(TOS.test, biome=="DXS")

## correlation matrix
# temperate forest
TOS.cor.tempfor <- na.omit(TOS.test.tempfor[,c("lHg","lgs","lgtclay","lgtsilt","pH15","lelecond","lLOI","prescott","soilH20","llai","lgpp","soillN","soillP",
                               "lAl","lFe","lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric4","ferric2","Alox","Feox","Pox")])
VIF.tempfor <- usdm::vif(TOS.cor.tempfor)
VIF.tempfor[which(VIF.tempfor$VIF > 4),]
usdm::vifstep(TOS.cor.tempfor, th=5)

# remove multicollinear variables
TOS.cor.tempfor.brt <- TOS.cor.tempfor %>% dplyr::select(-c('lFe', 'lV', 'lgtsilt', 'ferric4'))
usdm::vif(TOS.cor.tempfor.brt)

# Mediterranean forest
TOS.cor.Medfor <- na.omit(TOS.test.Medfor[,c("lHg","lgs","lgtclay","lgtsilt","pH15","lelecond","lLOI","prescott","soilH20","llai","lgpp","soillN","soillP",
                               "lAl","lFe","lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric4","ferric2","Alox","Feox","Pox")])
VIF.Medfor <- usdm::vif(TOS.cor.Medfor)
VIF.Medfor[which(VIF.Medfor$VIF > 4),]
usdm::vifstep(TOS.cor.tempfor, th=5)

# remove multicollinear variables
TOS.cor.Medfor.brt <- TOS.cor.Medfor %>% dplyr::select(-c('lFe', 'lV', 'lgpp','lAl','prescott','ferric4'))
usdm::vif(TOS.cor.Medfor.brt)

# temperate grassland/savanna
TOS.cor.tempgrass <- na.omit(TOS.test.tempgrass[,c("lHg","lgs","lgtclay","lgtsilt","pH15","lelecond","lLOI","prescott","soilH20","llai","lgpp","soillN","soillP",
                               "lAl","lFe","lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric4","ferric2","Alox","Feox","Pox")])
VIF.tempgrass <- usdm::vif(TOS.cor.tempgrass)
VIF.tempgrass[which(VIF.tempgrass$VIF > 4),]

# remove multicollinear variables
TOS.cor.tempgrass.brt <- TOS.cor.tempgrass %>% dplyr::select(-c('lFe', 'lV', 'lgpp','llai', "ferric4"))
usdm::vif(TOS.cor.tempgrass.brt)

# tropical forest (n = 28 only; insufficient for VIF)
TOS.cor.tropfor <- na.omit(TOS.test.tropfor[,c("lHg","lgs","lgtclay","lgtsilt","pH15","lelecond","lLOI","prescott","soilH20","llai","lgpp","soillN","soillP",
                               "lAl","lFe","lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric4","ferric2","Alox","Feox","Pox")])
VIF.tropfor <- usdm::vif(TOS.cor.tropfor)
VIF.tropfor[which(VIF.tropfor$VIF > 4),]
dim(TOS.cor.tropfor)

# tropical grassland/savanna
TOS.cor.tropgrass <- na.omit(TOS.test.tropgrass[,c("lHg","lgs","lgtclay","lgtsilt","pH15","lelecond","lLOI","prescott","soilH20","llai","lgpp","soillN","soillP",
                               "lAl","lFe","lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric4","ferric2","Alox","Feox","Pox")])
# remove infinite values
TOS.cor.tropgrass$lgtclay[(which(is.infinite(TOS.cor.tropgrass$lgtclay)==T))] <- NA

VIF.tropgrass <- usdm::vif(TOS.cor.tropgrass)
VIF.tropgrass[which(VIF.tropgrass$VIF > 4),]

# remove multicollinear variables
TOS.cor.tropgrass.brt <- TOS.cor.tropgrass %>% dplyr::select(-c('lFe', 'lV', 'lgpp','lgtsilt', "ferric4"))
usdm::vif(TOS.cor.tropgrass.brt)

# dryland/xeric
TOS.cor.dryxeric <- na.omit(TOS.test.dryxeric[,c("lHg","lgs","lgtclay","lgtsilt","pH15","lelecond","lLOI","prescott","soilH20","llai","lgpp","soillN","soillP",
                               "lAl","lFe","lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric4","ferric2","Alox","Feox","Pox")])
VIF.dryxeric <- usdm::vif(TOS.cor.dryxeric)
VIF.dryxeric[which(VIF.dryxeric$VIF > 4),]

# remove multicollinear variables
TOS.cor.dryxeric.brt <- TOS.cor.dryxeric %>% dplyr::select(-c('lFe', 'lCu','lZn', 'lAl','lgtsilt','lgpp', "ferric4"))
usdm::vif(TOS.cor.dryxeric.brt)


## full dataset correlation matrix
TOS.cor <- na.omit(TOS.test[,c("lHg","lgs","lgtclay","lgtsilt","pH15","lelecond","lLOI","prescott","soilH20","llai","lgpp","soillN","soillP",
                "lAl","lFe","lMn","lCu","lZn","lPb","lSb","lNi","lV","lKThU","ferric4","ferric2","Alox","Feox","Pox")])
class(TOS.cor)
dim(TOS.cor)
head(TOS.cor)
round(cor(TOS.cor, method="spearman"), 2)

## variance inflation
VIF.TOS <- usdm::vif(TOS.cor)
VIF.TOS
VIF.TOS[which(VIF.TOS$VIF > 4),]

## combine VIF tables
VIF.biome <- cbind(all=VIF.TOS$VIF, tempfor=VIF.tempfor$VIF, Medfor=VIF.Medfor$VIF, tempgrass=VIF.tempgrass$VIF, tropfor=VIF.tropfor$VIF,
                   tropgrass=VIF.tropgrass$VIF, dryxeric=VIF.dryxeric$VIF)
row.names(VIF.biome) <- VIF.tempfor$Variable
VIF.biome

# plot heatmap of VIF
heatmap(VIF.biome[,-5], col=rev(heat.colors(256)), scale="none", na.rm=T, margins=c(4,4))


## limit variables after VIF inspection
TOS.test.brt <- TOS.cor[,c("lHg","lgs","lgtclay","pH15","lelecond","lLOI","prescott","soilH20","llai","soillN","soillP",
                            "lAl","lMn","lCu","lPb","lSb","lNi","lKThU","ferric2","Alox","Feox","Pox")]
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

gbm.plot(TOS.brt, variable.no=18, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]", x.label="ferric PC4", plot.layout=c(1,1))
gbm.plot(TOS.brt, variable.no=5, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]", x.label="log10 LOI", plot.layout=c(1,1))
gbm.plot(TOS.brt, variable.no=9, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]", x.label="log10 soil nitrogen", plot.layout=c(1,1))
gbm.plot(TOS.brt, variable.no=16, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]", x.label="log10 nickel", plot.layout=c(1,1))
gbm.plot(TOS.brt, variable.no=4, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]", x.label="log10 electro-conductivity", plot.layout=c(1,1))
gbm.plot(TOS.brt, variable.no=14, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 [Hg]", x.label="log10 lead", plot.layout=c(1,1))

gbm.plot.fits(TOS.brt)

TOS.brt.partial.deps <- list()
TOS.pred.names <- names(TOS.test.brt)[2:19]
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


## biome-specific BRTs
# temperate forests
TOS.cor.tempfor.brt
dim(TOS.cor.tempfor.brt)
colnames(TOS.cor.tempfor.brt)

## boosted regression tree
tempfor.brt <- gbm.step(TOS.cor.tempfor.brt, gbm.x = attr(TOS.cor.tempfor.brt, "names")[c(2:length(colnames(TOS.cor.tempfor.brt)))],
                    gbm.y = attr(TOS.cor.tempfor.brt, "names")[1], family="gaussian", max.trees=100000,
                    tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.75,
                    tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(tempfor.brt)
tempfor.brt.summ <- summary(tempfor.brt)

# plots
tempfor.ri.plt <- ggplot(tempfor.brt.summ) +
  geom_bar(aes(x=reorder(row.names(tempfor.brt.summ), rel.inf), y=rel.inf), stat="identity", fill="blue", alpha=0.7)
tempfor.ri.plt + coord_flip() +
  ylab("relative influence (temperate forest)") + xlab("")
tempfor.ri.plt.flip <- tempfor.ri.plt + coord_flip() +
  ylab("relative influence (temperate forest)") + xlab("")

## plot predicted relationships of top x variables
topNvars <- 9 # x
topNvar.names <- tempfor.brt.summ[1:topNvars,]$var
plotNvec <- paste("tempfor.plt",1:topNvars,sep="")

for (v in 1:topNvars) {
  var.dat <- as.data.frame(plot.gbm(tempfor.brt, i.var=v, continuous.resolution = 100, return.grid=T))
  colnames(var.dat) <- c("V1", "V2")
  assign(plotNvec[v], ggplot(data=var.dat) +
           geom_line(aes(x = V1, y = V2), col="blue") +
           xlab(topNvar.names[v]) +
           ylab("log10 [Hg]"))
}
grid.arrange(tempfor.plt1, tempfor.plt2, tempfor.plt3, tempfor.plt4, tempfor.plt5, tempfor.plt6,
             tempfor.plt7, tempfor.plt8, tempfor.plt9, ncol=3)

tempfor.CV <- 100 * tempfor.brt$cv.statistics$correlation.mean
tempfor.CV.se <- 100 * tempfor.brt$cv.statistics$correlation.se
print(c(tempfor.CV, tempfor.CV.se))


# Mediterranean forests
TOS.cor.Medfor.brt
dim(TOS.cor.Medfor.brt)
colnames(TOS.cor.Medfor.brt)

## boosted regression tree
Medfor.brt <- gbm.step(TOS.cor.Medfor.brt, gbm.x = attr(TOS.cor.Medfor.brt, "names")[c(2:length(colnames(TOS.cor.Medfor.brt)))],
                    gbm.y = attr(TOS.cor.Medfor.brt, "names")[1], family="gaussian", max.trees=100000,
                    tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.75,
                    tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(Medfor.brt)
Medfor.brt.summ <- summary(Medfor.brt)

# plots
Medfor.ri.plt <- ggplot(Medfor.brt.summ) +
  geom_bar(aes(x=reorder(row.names(Medfor.brt.summ), rel.inf), y=rel.inf), stat="identity", fill="blue", alpha=0.7)
Medfor.ri.plt + coord_flip() +
  ylab("relative influence (Mediterranean forest)") + xlab("")
Medfor.ri.plt.flip <- Medfor.ri.plt + coord_flip() +
  ylab("relative influence (Mediterranean forest)") + xlab("")

## plot predicted relationships of top x variables
topNvars <- 9 # x
topNvar.names <- Medfor.brt.summ[1:topNvars,]$var
plotNvec <- paste("Medfor.plt",1:topNvars,sep="")
for (v in 1:topNvars) {
  var.dat <- as.data.frame(plot.gbm(Medfor.brt, i.var=v, continuous.resolution = 100, return.grid=T))
  colnames(var.dat) <- c("V1", "V2")
  assign(plotNvec[v], ggplot(data=var.dat) +
           geom_line(aes(x = V1, y = V2), col="blue") +
           xlab(topNvar.names[v]) +
           ylab("log10 [Hg]"))
}
grid.arrange(Medfor.plt1, Medfor.plt2, Medfor.plt3, Medfor.plt4, Medfor.plt5, Medfor.plt6,
             Medfor.plt7, Medfor.plt8, Medfor.plt9, ncol=3)

Medfor.CV <- 100 * Medfor.brt$cv.statistics$correlation.mean
Medfor.CV.se <- 100 * Medfor.brt$cv.statistics$correlation.se
print(c(Medfor.CV, Medfor.CV.se))


# temperate grassland/savanna
TOS.cor.tempgrass.brt
dim(TOS.cor.tempgrass.brt)
colnames(TOS.cor.tempgrass.brt)

## boosted regression tree
tempgrass.brt <- gbm.step(TOS.cor.tempgrass.brt, gbm.x = attr(TOS.cor.tempgrass.brt, "names")[c(2:length(colnames(TOS.cor.tempgrass.brt)))],
                    gbm.y = attr(TOS.cor.tempgrass.brt, "names")[1], family="gaussian", max.trees=100000,
                    tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.75,
                    tree.complexity = 2, silent=F, tolerance.method = "auto")

summary(tempgrass.brt)
tempgrass.brt.summ <- summary(tempgrass.brt)

# plots
tempgrass.ri.plt <- ggplot(tempgrass.brt.summ) +
  geom_bar(aes(x=reorder(row.names(tempgrass.brt.summ), rel.inf), y=rel.inf), stat="identity", fill="blue", alpha=0.7)
tempgrass.ri.plt + coord_flip() +
  ylab("relative influence (temperate grassland/savanna)") + xlab("")
tempgrass.ri.plt.flip <- tempgrass.ri.plt + coord_flip() +
  ylab("relative influence (temperate grassland/savanna)") + xlab("")

## plot predicted relationships of top x variables
topNvars <- 9 # x
topNvar.names <- tempgrass.brt.summ[1:topNvars,]$var
plotNvec <- paste("tempgrass.plt",1:topNvars,sep="")
for (v in 1:topNvars) {
  var.dat <- as.data.frame(plot.gbm(tempgrass.brt, i.var=v, continuous.resolution = 100, return.grid=T))
  colnames(var.dat) <- c("V1", "V2")
  assign(plotNvec[v], ggplot(data=var.dat) +
           geom_line(aes(x = V1, y = V2), col="blue") +
           xlab(topNvar.names[v]) +
           ylab("log10 [Hg]"))
}

grid.arrange(tempgrass.plt1, tempgrass.plt2, tempgrass.plt3, tempgrass.plt4, tempgrass.plt5, tempgrass.plt6,
             tempgrass.plt7, tempgrass.plt8, tempgrass.plt9, ncol=3)


tempgrass.CV <- 100 * tempgrass.brt$cv.statistics$correlation.mean
tempgrass.CV.se <- 100 * tempgrass.brt$cv.statistics$correlation.se
print(c(tempgrass.CV, tempgrass.CV.se))


# tropical grasslands/savannas
TOS.cor.tropgrass.brt
dim(TOS.cor.tropgrass.brt)
colnames(TOS.cor.tropgrass.brt)

## boosted regression tree
tropgrass.brt <- gbm.step(TOS.cor.tropgrass.brt, gbm.x = attr(TOS.cor.tropgrass.brt, "names")[c(2:length(colnames(TOS.cor.tropgrass.brt)))],
                    gbm.y = attr(TOS.cor.tropgrass.brt, "names")[1], family="gaussian", max.trees=100000,
                    tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.75,
                    tree.complexity = 2, silent=F, tolerance.method = "auto")

summary(tropgrass.brt)
tropgrass.brt.summ <- summary(tropgrass.brt)

# plots
tropgrass.ri.plt <- ggplot(tropgrass.brt.summ) +
  geom_bar(aes(x=reorder(row.names(tropgrass.brt.summ), rel.inf), y=rel.inf), stat="identity", fill="blue", alpha=0.7)
tropgrass.ri.plt + coord_flip() +
  ylab("relative influence (tropical grassland/savanna)") + xlab("")
tropgrass.ri.plt.flip <- tropgrass.ri.plt + coord_flip() +
  ylab("relative influence (tropical grassland/savanna)") + xlab("")

## plot predicted relationships of top x variables
topNvars <- 9 # x
topNvar.names <- tropgrass.brt.summ[1:topNvars,]$var
plotNvec <- paste("tropgrass.plt",1:topNvars,sep="")
for (v in 1:topNvars) {
  var.dat <- as.data.frame(plot.gbm(tropgrass.brt, i.var=v, continuous.resolution = 100, return.grid=T))
  colnames(var.dat) <- c("V1", "V2")
  assign(plotNvec[v], ggplot(data=var.dat) +
           geom_line(aes(x = V1, y = V2), col="blue") +
           xlab(topNvar.names[v]) +
           ylab("log10 [Hg]"))
}

grid.arrange(tropgrass.plt1, tropgrass.plt2, tropgrass.plt3, tropgrass.plt4, tropgrass.plt5, tropgrass.plt6,
             tropgrass.plt7, tropgrass.plt8, tropgrass.plt9, ncol=3)

tropgrass.CV <- 100 * tropgrass.brt$cv.statistics$correlation.mean
tropgrass.CV.se <- 100 * tropgrass.brt$cv.statistics$correlation.se
print(c(tropgrass.CV, tropgrass.CV.se))


# dryland/xeric
TOS.cor.dryxeric.brt
dim(TOS.cor.dryxeric.brt)
colnames(TOS.cor.dryxeric.brt)

## boosted regression tree
dryxeric.brt <- gbm.step(TOS.cor.dryxeric.brt, gbm.x = attr(TOS.cor.dryxeric.brt, "names")[c(2:length(colnames(TOS.cor.dryxeric.brt)))],
                    gbm.y = attr(TOS.cor.dryxeric.brt, "names")[1], family="gaussian", max.trees=100000,
                    tolerance = 0.002, learning.rate = 0.001, bag.fraction=0.75,
                    tree.complexity = 2, silent=F, tolerance.method = "auto")

summary(dryxeric.brt)
dryxeric.brt.summ <- summary(dryxeric.brt)

# plots
dryxeric.ri.plt <- ggplot(dryxeric.brt.summ) +
  geom_bar(aes(x=reorder(row.names(dryxeric.brt.summ), rel.inf), y=rel.inf), stat="identity", fill="blue", alpha=0.7)
dryxeric.ri.plt + coord_flip() +
  ylab("relative influence (dryland/xeric)") + xlab("")
dryxeric.ri.plt.flip <- dryxeric.ri.plt + coord_flip() +
  ylab("relative influence (dryland/xeric)") + xlab("")

## plot predicted relationships of top x variables
topNvars <- 9 # x
topNvar.names <- dryxeric.brt.summ[1:topNvars,]$var
plotNvec <- paste("dryxeric.plt",1:topNvars,sep="")
for (v in 1:topNvars) {
  var.dat <- as.data.frame(plot.gbm(dryxeric.brt, i.var=v, continuous.resolution = 100, return.grid=T))
  colnames(var.dat) <- c("V1", "V2")
  assign(plotNvec[v], ggplot(data=var.dat) +
           geom_line(aes(x = V1, y = V2), col="blue") +
           xlab(topNvar.names[v]) +
           ylab("log10 [Hg]"))
}

grid.arrange(dryxeric.plt1, dryxeric.plt2, dryxeric.plt3, dryxeric.plt4, dryxeric.plt5, dryxeric.plt6,
             dryxeric.plt7, dryxeric.plt8, dryxeric.plt9, ncol=3)

dryxeric.CV <- 100 * dryxeric.brt$cv.statistics$correlation.mean
dryxeric.CV.se <- 100 * dryxeric.brt$cv.statistics$correlation.se
print(c(dryxeric.CV, dryxeric.CV.se))


# plot rel infl plots for each biome together
grid.arrange(tempfor.ri.plt.flip, Medfor.ri.plt.flip, tempgrass.ri.plt.flip, tropgrass.ri.plt.flip,
             dryxeric.ri.plt.flip, ncol=3)
biome.brt.CVs <- data.frame(CV=c(tempfor.CV, Medfor.CV, tempgrass.CV, tropgrass.CV, dryxeric.CV),
                            CV.se=c(tempfor.CV.se, Medfor.CV.se, tempgrass.CV.se, tropgrass.CV.se, dryxeric.CV.se))
rownames(biome.brt.CVs) <- c("tempfor", "Medfor", "tempgrass", "tropgrass", "dryxeric")
biome.brt.CVs

# merge brt outputs per biome to single table
biome.brt1 <- merge(tempfor.brt.summ, Medfor.brt.summ, by="var", all.x=T, all.y=T)
biome.brt2 <- merge(biome.brt1, tempgrass.brt.summ, by="var", all.x=T, all.y=T)
biome.brt3 <- merge(biome.brt2, tropgrass.brt.summ, by="var", all.x=T, all.y=T)
biome.brt <- merge(biome.brt3, dryxeric.brt.summ, by="var", all.x=T, all.y=T)
colnames(biome.brt)[2:6] <- rownames(biome.brt.CVs)
biome.brt$relInfl.mn <- rowMeans(biome.brt[,2:6], na.rm=T)
biome.brt
biome.sort <- biome.brt[order(biome.brt$relInfl.mn, decreasing=T),]
biome.sort

# plot relative influence in multi-column plot
library(data.table)
biome.sort2 <- biome.sort[,c('var','tempfor','Medfor','tempgrass','tropgrass','dryxeric')]
biome.mlt <- reshape2::melt(biome.sort2, id.vars = 1, verbose=T)
biome.mlt
head(biome.mlt)
write.csv(biome.sort, file="biomes_brt_relInfl.csv", row.names=F)

ggplot(biome.mlt, aes(x = var, y = value)) + 
  geom_bar(aes(fill = variable), stat = "identity", position = "dodge") +
  coord_flip() + xlab("relative influence") + ylab("") + theme(axis.text.y=element_text(size=8)) +
  labs(fill = 'biome')

ggplot(biome.mlt, aes(x = var, y = value)) + 
  geom_bar(aes(fill = variable), stat = "identity", position = "stack") +
  coord_flip() + ylab("relative influence") + xlab("") + theme(axis.text.y=element_text(size=8)) +
  labs(fill = 'biome')


#########################
## BOS
colnames(BOS.dat)
BOS.test1 <- data.frame(SITEID=BOS.dat$SITEID, LAT=BOS.dat$LAT.x, LON=BOS.dat$LON.x, lHg=log10(BOS.dat$HgCOMP), lgs=log10(BOS.dat$gs.wmn),
                        lgtclay=logit(BOS.dat$CLAY/100), lgtsoilclay=logit(BOS.dat$soilclay/100),
                        lgtsilt=logit(BOS.dat$SILT/100), lgtsoilsilt=logit(BOS.dat$soilsilt/100),  
                        pH=BOS.dat$BOSpH, pH15=BOS.dat$pH15, soilpH=BOS.dat$soilpH, lelecond=log10(BOS.dat$EC15), lLOI=log10(BOS.dat$LOI),
                        lrain=log10(BOS.dat$rain), prescott=BOS.dat$prescott, soilH20=BOS.dat$soilH20,
                        llai=log10(BOS.dat$lai), lgpp=log10(BOS.dat$gpp),
                        soillN=log10(BOS.dat$soilN), soillP=log10(BOS.dat$soilP),
                        lAl=log10(BOS.dat$Al), lFe=log10(BOS.dat$Fe), lMn=log10(BOS.dat$Mn), lCu=log10(BOS.dat$CuICPMS),
                        lZn=log10(BOS.dat$ZnICPMS), lPb=log10(BOS.dat$PbICPMS), lSb=log10(BOS.dat$SbICPMS),
                        lNi=log10(BOS.dat$NiCPMS), lV=log10(BOS.dat$V), lKThU=log10(BOS.dat$KThU),
                        lucat=BOS.dat$lucat, biome=BOS.dat$biome, geol=BOS.dat$lithgrp, soils=BOS.dat$soils,
                        ferric4=BOS.dat$ferricpc4, ferric2=BOS.dat$ferricpc2, Alox=BOS.dat$Alox, 
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
length(which(is.na(BOS.test$lgtsoilclay)==T)) # 28 missing
length(which(is.na(BOS.test$lgtsilt)==T)) # good
length(which(is.na(BOS.test$lgtsoilsilt)==T)) # 28 missing
length(which(is.na(BOS.test$pH15)==T)) # good
length(which(is.na(BOS.test$soilpH)==T)) # 12 missing
length(which(is.na(BOS.test$pH)==T)) # good
length(which(is.na(BOS.test$lelecond)==T)) # good
length(which(is.na(BOS.test$lLOI)==T)) # 9 missing
length(which(is.na(BOS.test$lrain)==T)) # 54 missing
length(which(is.na(BOS.test$prescott)==T)) # good
length(which(is.na(BOS.test$soilH20)==T)) # 4 missing
length(which(is.na(BOS.test$llai)==T)) # 4 missing
length(which(is.na(BOS.test$lgpp)==T)) # 4 missing
length(which(is.na(BOS.test$soillN)==T)) # 20 missing
length(which(is.na(BOS.test$soillP)==T)) # 20 missing
length(which(is.na(BOS.test$lAl)==T)) # 9 missing
length(which(is.na(BOS.test$lFe)==T)) # 9 missing
length(which(is.na(BOS.test$lMn)==T)) # 9 missing
length(which(is.na(BOS.test$lCu)==T)) # 9 missing
length(which(is.na(BOS.test$lZn)==T)) # 9 missing
length(which(is.na(BOS.test$lPb)==T)) # 9 missing
length(which(is.na(BOS.test$lSb)==T)) # 9 missing
length(which(is.na(BOS.test$lNi)==T)) # 9 missing
length(which(is.na(BOS.test$lV)==T)) # 9 missing
length(which(is.na(BOS.test$lKThU)==T)) # good

## any infinites?
length(which(is.infinite(BOS.test$lHg)==T))
length(which(is.infinite(BOS.test$lgs)==T))
length(which(is.infinite(BOS.test$pH15)==T))
length(which(is.infinite(BOS.test$lgtclay)==T)) # 7 missing
length(which(is.infinite(BOS.test$lgtsilt)==T)) # 3 missing
length(which(is.infinite(BOS.test$lelecond)==T))
length(which(is.infinite(BOS.test$lLOI)==T)) # 1 missing
length(which(is.infinite(BOS.test$lrain)==T))
length(which(is.infinite(BOS.test$prescott)==T))
length(which(is.infinite(BOS.test$soilH20)==T))
length(which(is.infinite(BOS.test$llai)==T)) # 4 missing
length(which(is.infinite(BOS.test$lgpp)==T)) # 4 missing
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

dim(BOS.test)

# replace infinties with NA
BOS.test$lgtclay[(which(is.infinite(BOS.test$lgtclay)==T))] <- NA
BOS.test$lgtsilt[(which(is.infinite(BOS.test$lgtsilt)==T))] <- NA
BOS.test$lLOI[(which(is.infinite(BOS.test$lLOI)==T))] <- NA
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
usdm::vif(BOS.cor[,c("lHg","lgs","lgtclay","pH15","lelecond","lLOI","prescott","soilH20","llai","soillN","soillP",
                     "lAl","lMn","lCu","lPb","lSb","lNi","lKThU","ferric2","Alox","Feox","Pox")])

## limit variables after VIF inspection
BOS.test.brt <- BOS.cor[,c("lHg","lgs","lgtclay","pH15","lelecond","lLOI","prescott","soilH20","llai","soillN","soillP",
                           "lAl","lMn","lCu","lPb","lSb","lNi","lKThU","ferric2","Alox","Feox","Pox")]
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
BOS.pred.names <- names(BOS.test.brt)[2:18]

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


#########################################
# randomised BRT using distance matrix ##
#########################################
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
coords.ran <- select_distant_points(df = TOS.brt.dat.rsmp.coords, min_dist = min.dist, target_n = npts2generate,
                      dist_matrix = TOS.brt.dat.rsmp_haversine_matrix, x_col = "LON", y_col = "LAT")
dim(coords.ran)
head(coords.ran)
coords.ran.pts <- vect(cbind(coords.ran$LON, 
                                           coords.ran$LAT), crs="+proj=longlat")
terra::plot(coords.ran.pts)

# subset random points for BRT training data
TOS.ran.subset <- na.omit(TOS.brt.dat.rsmp[as.numeric(row.names(coords.ran)), ])
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

# resampled BRT
biter <- 10000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()
eq.sp.pts <- 100
traincols <- attr(TOS.ran.subset, "names")[c(4:length(colnames(TOS.ran.subset)))] # variable columns used to train data
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter), dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

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
  summ.fit <- summary(TOS.ran.brt)
  
  if (is.null(TOS.ran.brt) == F) {
    TOS.ran.brt.old <- TOS.ran.brt
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*TOS.ran.brt$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*TOS.ran.brt$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(TOS.ran.brt, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(TOS.ran.brt, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- TOS.ran.brt$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- TOS.ran.brt$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
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


# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                    (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                       (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                        NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med <- median(CV.cor.update, na.rm=T)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo, CV.cor.med, CV.cor.up))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort

# plot
ri.plt <- ggplot(ri.sort) +
  geom_bar(aes(x=reorder(row.names(ri.sort), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar( aes(x=row.names(ri.sort), ymin=ri.lo, ymax=ri.up),
                 linewidth=0.4, colour="black", alpha=0.9)
ri.plt + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo,CV.cor.med,CV.cor.up), 2))

## plot predicted relationships of top x variables
topNvars <- 9 # x
head(pred.med)
ri.sort
top.ri.sort <- ri.sort[1:topNvars,]
topNvars.names <- rownames(top.ri.sort)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                           pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
    geom_line(aes(x=V1, y=V2), colour="blue") +
    geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
    lims(y=ylims) +
    xlab(topNvars.names[v]) + ylab("log10 [Hg]"))
}
(plt1 + plt2 + plt3) / (plt4 + plt5 + plt6) / (plt7 + plt8 + plt9)

# export results
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                         up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "Pred.csv", sep=""), row.names=F)
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

colnames(TOS.dat)
TOS.cl.test <- data.frame(SITEID=TOS.dat$SITEID, LAT=TOS.dat$LAT.x, LON=TOS.dat$LON.x, state=TOS.dat$STATE,
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

min.dist <- 175000

coords.ran <- select_distant_points(df = TOS.lm.coords, min_dist = min.dist,
                                    dist_matrix = TOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
dim(coords.ran)
head(coords.ran)
coords.ran.pts <- vect(cbind(coords.ran$LON, 
                             coords.ran$LAT), crs="+proj=longlat")
terra::plot(coords.ran.pts)

# subset random points
TOS.ran.subset <- na.omit(TOS.lm[as.numeric(row.names(coords.ran)), ])
head(TOS.ran.subset)
dim(TOS.ran.subset)

# by state
table(TOS.ran.subset$state)
HgXstate.stats.ran <- TOS.ran.subset %>%
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
HgXstate.stats.ran
na.omit(HgXstate.stats.ran)
hist(TOS.ran.subset$lHg)


mod1 <- lm(lHg ~ state, data=TOS.ran.subset)
mod.null <- lm(lHg ~ 1, data=TOS.ran.subset)
check_model(mod1, detrend=T)
plot_model(mod1)
pmmod1 <- plot_model(mod1)
pmmod1.coef <- data.frame(state=as.character(pmmod1[[1]]$term), lo=pmmod1[[1]]$conf.low, up=pmmod1[[1]]$conf.high)

pmmod1.coef$nonzero <- ifelse(contains_zero_sign(pmmod1.coef$lo, pmmod1.coef$up)==F, 1, 0)
pmmod1.coef
nonzerosum <- sum(pmmod1.coef$nonzero)
pmmod1.coef[which(pmmod1.coef$nonzero==1),]$state

wAICc <- weight.IC(delta.IC(c(AICc(mod1),AICc(mod.null))))
ER <- wAICc[1]/wAICc[2]
ER


# resampling loop
iter <- 1000
itdiv <- iter/10

# storage matrix
table(TOS.lm$state)
nlevels <- length(table(TOS.lm$state))
stor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*nlevels)
statedir.lab <- paste("state",attr(table(TOS.lm$state), "names"),"dir",sep="")
colnames(stor.mat) <- c("ER", "nonzerosum",
                        paste("state",attr(table(TOS.lm$state), "names"),sep=""),statedir.lab)
head(stor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = TOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = TOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  TOS.ran.subset <- na.omit(TOS.lm[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  mod1 <- lm(lHg ~ state, data=TOS.ran.subset) # class level model
  mod.null <- lm(lHg ~ 1, data=TOS.ran.subset) # null model
  
  # coefficient boundaries
  pmmod1 <- plot_model(mod1)
  pmmod1.coef <- data.frame(state=as.character(pmmod1[[1]]$term), lo=pmmod1[[1]]$conf.low,
                            up=pmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  pmmod1.coef$nonzero <- ifelse(contains_zero_sign(pmmod1.coef$lo, pmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  pmmod1.coef$dir <- ifelse(sign(pmmod1.coef$lo) == -1 & sign(pmmod1.coef$up) == -1, -1, 0)
  pmmod1.coef$dir <- ifelse(sign(pmmod1.coef$lo) == 1 & sign(pmmod1.coef$up) == 1, 1, pmmod1.coef$dir)
    
  # how many variables are non-zero?
  stor.mat[i,"nonzerosum"] <- sum(pmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  nzstates <- pmmod1.coef[which(pmmod1.coef$nonzero==1),]$state
  stor.mat[i, which(colnames(stor.mat) %in% nzstates)] <- 1
  
  # for non-zeros, what is the direction?
  nzstatesdir <- paste(nzstates,"dir",sep="")
  stor.mat[i, which(colnames(stor.mat) %in% nzstatesdir)] <- pmmod1.coef$dir[which(pmmod1.coef$dir != 0)]
  
  # AICc model comparison
  wAICc <- weight.IC(delta.IC(c(AICc(mod1),AICc(mod.null))))

  # store evidence ratio
  stor.mat[i,"ER"] <- wAICc[1]/wAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

# ER
median(stor.mat[,"ER"])
quantile(stor.mat[,"ER"], probs=c(0.025,0.975))

head(stor.mat)
state.labs <- c("NSW","NT","QLD","SA","TAS","VIC","WA")
lennzNSW <- length(table(stor.mat[,"stateNSW"]))
lennzNT <- length(table(stor.mat[,"stateNT"]))
lennzQLD <- length(table(stor.mat[,"stateQLD"]))
lennzSA <- length(table(stor.mat[,"stateSA"]))
lennzTAS <- length(table(stor.mat[,"stateTAS"]))
lennzVIC <- length(table(stor.mat[,"stateVIC"]))
lennzWA <- length(table(stor.mat[,"stateWA"]))


nzrslts <- c(ifelse(lennzNSW==1, as.numeric(table(stor.mat[,"stateNSW"])), 0),
  ifelse(lennzNT==1, as.numeric(table(stor.mat[,"stateNT"])), 0),
  ifelse(lennzQLD==1, as.numeric(table(stor.mat[,"stateQLD"])), 0),
  ifelse(lennzSA==1, as.numeric(table(stor.mat[,"stateSA"])), 0),
  ifelse(lennzTAS==1, as.numeric(table(stor.mat[,"stateTAS"])), 0),
  ifelse(lennzVIC==1, as.numeric(table(stor.mat[,"stateVIC"])), 0),
  ifelse(lennzWA==1, as.numeric(table(stor.mat[,"stateWA"])), 0))


dirNSWtab <- table(stor.mat[,"stateNSWdir"])
dirNTtab <- table(stor.mat[,"stateNTdir"])
dirQLDtab <- table(stor.mat[,"stateQLDdir"])
dirSAtab <- table(stor.mat[,"stateSAdir"])
dirTAStab <- table(stor.mat[,"stateTASdir"])
dirVICtab <- table(stor.mat[,"stateVICdir"])
dirWAtab <- table(stor.mat[,"stateWAdir"])

dirlenNSW <- length(dirNSWtab)
dirlenNT <- length(dirNTtab)
dirlenQLD <- length(dirQLDtab)
dirlenSA <- length(dirSAtab)
dirlenTAS <- length(dirTAStab)
dirlenVIC <- length(dirVICtab)
dirlenWA <- length(dirWAtab)

if (dirlenNSW == 0) {
  NSWpos <- 0
  NSWneg <- 0
} else if (dirlenNSW == 1 & as.numeric(attr(dirNSWtab,"names")[1]) == 1) {
  NSWpos <- as.numeric(dirNSWtab)[1]
  NSWneg <- 0
} else if (dirlenNSW == 1 & as.numeric(attr(dirNSWtab,"names")[1]) == -1) {
  NSWpos <- 0
  NSWneg <- as.numeric(dirNSWtab)[1]
} else if (dirlenNSW == 2) {
  NSWneg <- as.numeric(dirNSWtab)[1]
  NSWpos <- as.numeric(dirNSWtab)[2]
} else {
  NSWpos <- 0
  NSWneg <- 0
}

if (dirlenNT == 0) {
  NTpos <- 0
  NTneg <- 0
  } else if (dirlenNT == 1 & as.numeric(attr(dirNTtab,"names")[1]) == 1) {
    NTpos <- as.numeric(dirNTtab)[1]
    NTneg <- 0
  } else if (dirlenNT == 1 & as.numeric(attr(dirNTtab,"names")[1]) == -1) {
    NTpos <- 0
    NTneg <- as.numeric(dirNTtab)[1]
  } else if (dirlenNT == 2) {
    NTneg <- as.numeric(dirNTtab)[1]
    NTpos <- as.numeric(dirNTtab)[2]
  } else {
    NTpos <- 0
    NTneg <- 0
  }

if (dirlenQLD == 0) {
  QLDpos <- 0
  QLDneg <- 0
  } else if (dirlenQLD == 1 & as.numeric(attr(dirQLDtab,"names")[1]) == 1) {
    QLDpos <- as.numeric(dirQLDtab)[1]
    QLDneg <- 0
  } else if (dirlenQLD == 1 & as.numeric(attr(dirQLDtab,"names")[1]) == -1) {
    QLDpos <- 0
    QLDneg <- as.numeric(dirQLDtab)[1]
  } else if (dirlenQLD == 2) {
    QLDneg <- as.numeric(dirQLDtab)[1]
    QLDpos <- as.numeric(dirQLDtab)[2]
  } else {
    QLDpos <- 0
    QLDneg <- 0
  }

if (dirlenSA == 0) {
  SApos <- 0
  SAneg <- 0
  } else if (dirlenSA == 1 & as.numeric(attr(dirSAtab,"names")[1]) == 1) {
    SApos <- as.numeric(dirSAtab)[1]
    SAneg <- 0
  } else if (dirlenSA == 1 & as.numeric(attr(dirSAtab,"names")[1]) == -1) {
    SApos <- 0
    SAneg <- as.numeric(dirSAtab)[1]
  } else if (dirlenSA == 2) {
    SAneg <- as.numeric(dirSAtab)[1]
    SApos <- as.numeric(dirSAtab)[2]
  } else {
    SApos <- 0
    SAneg <- 0
  }

if (dirlenTAS == 0) {
  TASpos <- 0
  TASneg <- 0
  } else if (dirlenTAS == 1 & as.numeric(attr(dirTAStab,"names")[1]) == 1) {
    TASpos <- as.numeric(dirTAStab)[1]
    TASneg <- 0
  } else if (dirlenTAS == 1 & as.numeric(attr(dirTAStab,"names")[1]) == -1) {
    TASpos <- 0
    TASneg <- as.numeric(dirTAStab)[1]
  } else if (dirlenTAS == 2) {
    TASneg <- as.numeric(dirTAStab)[1]
    TASpos <- as.numeric(dirTAStab)[2]
  } else {
    TASpos <- 0
    TASneg <- 0
  }

if (dirlenVIC == 0) {
  VICpos <- 0
  VICneg <- 0
  } else if (dirlenVIC == 1 & as.numeric(attr(dirVICtab,"names")[1]) == 1) {
    VICpos <- as.numeric(dirVICtab)[1]
    VICneg <- 0
  } else if (dirlenVIC == 1 & as.numeric(attr(dirVICtab,"names")[1]) == -1) {
    VICpos <- 0
    VICneg <- as.numeric(dirVICtab)[1]
  } else if (dirlenVIC == 2) {
    VICneg <- as.numeric(dirVICtab)[1]
    VICpos <- as.numeric(dirVICtab)[2]
  } else {
    VICpos <- 0
    VICneg <- 0
  }

if (dirlenWA == 0) {
  WApos <- 0
  WAneg <- 0
  } else if (dirlenWA == 1 & as.numeric(attr(dirWAtab,"names")[1]) == 1) {
    WApos <- as.numeric(dirWAtab)[1]
    WAneg <- 0
  } else if (dirlenWA == 1 & as.numeric(attr(dirWAtab,"names")[1]) == -1) {
    WApos <- 0
    WAneg <- as.numeric(dirWAtab)[1]
  } else if (dirlenWA == 2) {
    WAneg <- as.numeric(dirWAtab)[1]
    WApos <- as.numeric(dirWAtab)[2]
  } else {
    WApos <- 0
    WAneg <- 0
  }


nznegdirrslts <- c(NSWneg, NTneg, QLDneg, SAneg, TASneg, VICneg, WAneg)
nzposdirrslts <- c(NSWpos, NTpos, QLDpos, SApos, TASpos, VICpos, WApos)

results.out <- data.frame(state=state.labs, nonzero=nzrslts/iter, negdir=nznegdirrslts/iter,
                          posdir=nzposdirrslts/iter)
results.out

# ER
10^median(log10(stor.mat[,"ER"]))
quantile(log10(stor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(stor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(stor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(stor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(stor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(stor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(stor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(stor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)



# by biome
table(TOS.ran.subset$biome)

table(TOS.ran.subset$biome)
HgXbiome.stats.ran <- TOS.ran.subset %>%
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
HgXbiome.stats.ran
na.omit(HgXbiome.stats.ran)

# resampling loop
# storage matrix
table(TOS.lm$biome)
nlevels <- length(table(TOS.lm$biome))
stor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*nlevels)
biomedir.lab <- paste("biome",attr(table(TOS.lm$biome), "names"),"dir",sep="")
colnames(stor.mat) <- c("ER", "nonzerosum",
                        paste("biome",attr(table(TOS.lm$biome), "names"),sep=""),biomedir.lab)
head(stor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = TOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = TOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  TOS.ran.subset <- na.omit(TOS.lm[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  mod1 <- lm(lHg ~ biome, data=TOS.ran.subset) # class level model
  mod.null <- lm(lHg ~ 1, data=TOS.ran.subset) # null model
  
  # coefficient boundaries
  pmmod1 <- plot_model(mod1)
  pmmod1.coef <- data.frame(biome=as.character(pmmod1[[1]]$term), lo=pmmod1[[1]]$conf.low,
                            up=pmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  pmmod1.coef$nonzero <- ifelse(contains_zero_sign(pmmod1.coef$lo, pmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  pmmod1.coef$dir <- ifelse(sign(pmmod1.coef$lo) == -1 & sign(pmmod1.coef$up) == -1, -1, 0)
  pmmod1.coef$dir <- ifelse(sign(pmmod1.coef$lo) == 1 & sign(pmmod1.coef$up) == 1, 1, pmmod1.coef$dir)
  
  # how many variables are non-zero?
  stor.mat[i,"nonzerosum"] <- sum(pmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  nzbiomes <- pmmod1.coef[which(pmmod1.coef$nonzero==1),]$biome
  stor.mat[i, which(colnames(stor.mat) %in% nzbiomes)] <- 1
  
  # for non-zeros, what is the direction?
  nzbiomessdir <- paste(nzbiomes,"dir",sep="")
  stor.mat[i, which(colnames(stor.mat) %in% nzbiomessdir)] <- pmmod1.coef$dir[which(pmmod1.coef$dir != 0)]
  
  # AICc model comparison
  wAICc <- weight.IC(delta.IC(c(AICc(mod1),AICc(mod.null))))
  
  # store evidence ratio
  stor.mat[i,"ER"] <- wAICc[1]/wAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

head(stor.mat)
biome.labs <- c("DXS","MFW","TBMF","TGSS","TSGSS","TSMBF")
lennzDXS <- length(table(stor.mat[,"biomeDXS"]))
lennzMFW <- length(table(stor.mat[,"biomeMFW"]))
lennzTBMF <- length(table(stor.mat[,"biomeTBMF"]))
lennzTGSS <- length(table(stor.mat[,"biomeTGSS"]))
lennzTSGSS <- length(table(stor.mat[,"biomeTSGSS"]))
lennzTSMBF <- length(table(stor.mat[,"biomeTSMBF"]))

nzrslts <- c(ifelse(lennzDXS==1, as.numeric(table(stor.mat[,"biomeDXS"])), 0),
             ifelse(lennzMFW==1, as.numeric(table(stor.mat[,"biomeMFW"])), 0),
             ifelse(lennzTBMF==1, as.numeric(table(stor.mat[,"biomeTBMF"])), 0),
             ifelse(lennzTGSS==1, as.numeric(table(stor.mat[,"biomeTGSS"])), 0),
             ifelse(lennzTSGSS==1, as.numeric(table(stor.mat[,"biomeTSGSS"])), 0),
             ifelse(lennzTSMBF==1, as.numeric(table(stor.mat[,"biomeTSMBF"])), 0))


dirDXStab <- table(stor.mat[,"biomeDXS"])
dirMFWtab <- table(stor.mat[,"biomeMFW"])
dirTBMFtab <- table(stor.mat[,"biomeTBMF"])
dirTGSStab <- table(stor.mat[,"biomeTGSS"])
dirTSGSStab <- table(stor.mat[,"biomeTSGSS"])
dirTSMBFtab <- table(stor.mat[,"biomeTSMBF"])

dirlenDXS <- length(dirDXStab)
dirlenMFW <- length(dirMFWtab)
dirlenTBMF <- length(dirTBMFtab)
dirlenTGSS <- length(dirTGSStab)
dirlenTSGSS <- length(dirTSGSStab)
dirlenTSMBF <- length(dirTSMBFtab)


if (dirlenDXS == 0) {
  DXSpos <- 0
  DXSneg <- 0
} else if (dirlenDXS == 1 & as.numeric(attr(dirDXStab,"names")[1]) == 1) {
  DXSpos <- as.numeric(dirDXStab)[1]
  DXSneg <- 0
} else if (dirlenDXS == 1 & as.numeric(attr(dirDXStab,"names")[1]) == -1) {
  DXSpos <- 0
  DXSneg <- as.numeric(dirDXStab)[1]
} else if (dirlenDXS == 2) {
  DXSneg <- as.numeric(dirDXStab)[1]
  DXSpos <- as.numeric(dirDXStab)[2]
} else {
  DXSpos <- 0
  DXSneg <- 0
}

if (dirlenMFW == 0) {
  MFWpos <- 0
  MFWneg <- 0
} else if (dirlenMFW == 1 & as.numeric(attr(dirMFWtab,"names")[1]) == 1) {
  MFWpos <- as.numeric(dirMFWtab)[1]
  MFWneg <- 0
} else if (dirlenMFW == 1 & as.numeric(attr(dirMFWtab,"names")[1]) == -1) {
  MFWpos <- 0
  MFWneg <- as.numeric(dirMFWtab)[1]
} else if (dirlenMFW == 2) {
  MFWneg <- as.numeric(dirMFWtab)[1]
  MFWpos <- as.numeric(dirMFWtab)[2]
} else {
  MFWpos <- 0
  MFWneg <- 0
}

if (dirlenTBMF == 0) {
  TBMFpos <- 0
  TBMFneg <- 0
} else if (dirlenTBMF == 1 & as.numeric(attr(dirTBMFtab,"names")[1]) == 1) {
  TBMFpos <- as.numeric(dirTBMFtab)[1]
  TBMFneg <- 0
} else if (dirlenTBMF == 1 & as.numeric(attr(dirTBMFtab,"names")[1]) == -1) {
  TBMFpos <- 0
  TBMFneg <- as.numeric(dirTBMFtab)[1]
} else if (dirlenTBMF == 2) {
  TBMFneg <- as.numeric(dirTBMFtab)[1]
  TBMFpos <- as.numeric(dirTBMFtab)[2]
} else {
  TBMFpos <- 0
  TBMFneg <- 0
}

if (dirlenTGSS == 0) {
  TGSSpos <- 0
  TGSSneg <- 0
} else if (dirlenTGSS == 1 & as.numeric(attr(dirTGSStab,"names")[1]) == 1) {
  TGSSpos <- as.numeric(dirTGSStab)[1]
  TGSSneg <- 0
} else if (dirlenTGSS == 1 & as.numeric(attr(dirTGSStab,"names")[1]) == -1) {
  TGSSpos <- 0
  TGSSneg <- as.numeric(dirTGSStab)[1]
} else if (dirlenTGSS == 2) {
  TGSSneg <- as.numeric(dirTGSStab)[1]
  TGSSpos <- as.numeric(dirTGSStab)[2]
} else {
  TGSSpos <- 0
  TGSSneg <- 0
}

if (dirlenTSGSS == 0) {
  TSGSSpos <- 0
  TSGSSneg <- 0
} else if (dirlenTSGSS == 1 & as.numeric(attr(dirTSGSStab,"names")[1]) == 1) {
  TSGSSpos <- as.numeric(dirTSGSStab)[1]
  TSGSSneg <- 0
} else if (dirlenTSGSS == 1 & as.numeric(attr(dirTSGSStab,"names")[1]) == -1) {
  TSGSSpos <- 0
  TSGSSneg <- as.numeric(dirTSGSStab)[1]
} else if (dirlenTSGSS == 2) {
  TSGSSneg <- as.numeric(dirTSGSStab)[1]
  TSGSSpos <- as.numeric(dirTSGSStab)[2]
} else {
  TSGSSpos <- 0
  TSGSSneg <- 0
}

if (dirlenTSMBF == 0) {
  TSMBFpos <- 0
  TSMBFneg <- 0
} else if (dirlenTSMBF == 1 & as.numeric(attr(dirTSMBFtab,"names")[1]) == 1) {
  TSMBFpos <- as.numeric(dirTSMBFtab)[1]
  TSMBFneg <- 0
} else if (dirlenTSMBF == 1 & as.numeric(attr(dirTSMBFtab,"names")[1]) == -1) {
  TSMBFpos <- 0
  TSMBFneg <- as.numeric(dirTSMBFtab)[1]
} else if (dirlenTSMBF == 2) {
  TSMBFneg <- as.numeric(dirTSMBFtab)[1]
  TSMBFpos <- as.numeric(dirTSMBFtab)[2]
} else {
  TSMBFpos <- 0
  TSMBFneg <- 0
}

nznegdirrslts <- c(DXSneg, MFWneg, TBMFneg, TGSSneg, TSGSSneg, TSMBFneg)
nzposdirrslts <- c(DXSpos, MFWpos, TBMFpos, TGSSpos, TSGSSpos, TSMBFpos)

results.out <- data.frame(biome=biomedir.lab, nonzero=nzrslts/iter, negdir=nznegdirrslts/iter,
                          posdir=nzposdirrslts/iter)
results.out

# ER
10^median(log10(stor.mat[,"ER"]))
quantile(log10(stor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(stor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(stor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(stor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(stor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(stor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(stor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(stor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)



# lithology
table(TOS.ran.subset$geol)

table(TOS.ran.subset$geol)
HgXgeol.stats.ran <- TOS.ran.subset %>%
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
HgXgeol.stats.ran
na.omit(HgXgeol.stats.ran)

# resampling loop
# storage matrix
table(TOS.lm$geol)
nlevels <- length(table(TOS.lm$geol))
stor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*nlevels)
geoldir.lab <- paste("geol",attr(table(TOS.lm$geol), "names"),"dir",sep="")
colnames(stor.mat) <- c("ER", "nonzerosum",
                        paste("geol",attr(table(TOS.lm$geol), "names"),sep=""),geoldir.lab)
head(stor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = TOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = TOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  TOS.ran.subset <- na.omit(TOS.lm[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  mod1 <- lm(lHg ~ geol, data=TOS.ran.subset) # class level model
  mod.null <- lm(lHg ~ 1, data=TOS.ran.subset) # null model
  
  # coefficient boundaries
  pmmod1 <- plot_model(mod1)
  pmmod1.coef <- data.frame(geol=as.character(pmmod1[[1]]$term), lo=pmmod1[[1]]$conf.low,
                            up=pmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  pmmod1.coef$nonzero <- ifelse(contains_zero_sign(pmmod1.coef$lo, pmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  pmmod1.coef$dir <- ifelse(sign(pmmod1.coef$lo) == -1 & sign(pmmod1.coef$up) == -1, -1, 0)
  pmmod1.coef$dir <- ifelse(sign(pmmod1.coef$lo) == 1 & sign(pmmod1.coef$up) == 1, 1, pmmod1.coef$dir)
  
  # how many variables are non-zero?
  stor.mat[i,"nonzerosum"] <- sum(pmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  nzgeol <- pmmod1.coef[which(pmmod1.coef$nonzero==1),]$geol
  stor.mat[i, which(colnames(stor.mat) %in% nzgeol)] <- 1
  
  # for non-zeros, what is the direction?
  nzgeoldir <- paste(nzgeol,"dir",sep="")
  stor.mat[i, which(colnames(stor.mat) %in% nzgeoldir)] <- pmmod1.coef$dir[which(pmmod1.coef$dir != 0)]
  
  # AICc model comparison
  wAICc <- weight.IC(delta.IC(c(AICc(mod1),AICc(mod.null))))
  
  # store evidence ratio
  stor.mat[i,"ER"] <- wAICc[1]/wAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

head(stor.mat)
geol.labs <- c("CAR","FEL","INT","MAF","MET","OTH","SIL")
lennzCAR <- length(table(stor.mat[,"geolCAR"]))
lennzFEL <- length(table(stor.mat[,"geolFEL"]))
lennzINT <- length(table(stor.mat[,"geolINT"]))
lennzMAF <- length(table(stor.mat[,"geolMAF"]))
lennzMET <- length(table(stor.mat[,"geolMET"]))
lennzOTH <- length(table(stor.mat[,"geolOTH"]))
lennzSIL <- length(table(stor.mat[,"geolSIL"]))

nzrslts <- c(ifelse(lennzCAR==1, as.numeric(table(stor.mat[,"geolCAR"])), 0),
             ifelse(lennzFEL==1, as.numeric(table(stor.mat[,"geolFEL"])), 0),
             ifelse(lennzINT==1, as.numeric(table(stor.mat[,"geolINT"])), 0),
             ifelse(lennzMAF==1, as.numeric(table(stor.mat[,"geolMAF"])), 0),
             ifelse(lennzMET==1, as.numeric(table(stor.mat[,"geolMET"])), 0),
             ifelse(lennzOTH==1, as.numeric(table(stor.mat[,"geolOTH"])), 0),
             ifelse(lennzSIL==1, as.numeric(table(stor.mat[,"geolSIL"])), 0))

dirCARtab <- table(stor.mat[,"geolCAR"])
dirFELtab <- table(stor.mat[,"geolFEL"])
dirINTtab <- table(stor.mat[,"geolINT"])
dirMAFtab <- table(stor.mat[,"geolMAF"])
dirMETtab <- table(stor.mat[,"geolMET"])
dirOTHtab <- table(stor.mat[,"geolOTH"])
dirSILtab <- table(stor.mat[,"geolSIL"])

dirlenCAR <- length(dirCARtab)
dirlenFEL <- length(dirFELtab)
dirlenINT <- length(dirINTtab)
dirlenMAF <- length(dirMAFtab)
dirlenMET <- length(dirMETtab)
dirlenOTH <- length(dirOTHtab)
dirlenSIL <- length(dirSILtab)


if (dirlenCAR == 0) {
  CARpos <- 0
  CARneg <- 0
} else if (dirlenCAR == 1 & as.numeric(attr(dirCARtab,"names")[1]) == 1) {
  CARpos <- as.numeric(dirCARtab)[1]
  CARneg <- 0
} else if (dirlenCAR == 1 & as.numeric(attr(dirCARtab,"names")[1]) == -1) {
  CARpos <- 0
  CARneg <- as.numeric(dirCARtab)[1]
} else if (dirlenCAR == 2) {
  CARneg <- as.numeric(dirCARtab)[1]
  CARpos <- as.numeric(dirCARtab)[2]
} else {
  CARpos <- 0
  CARneg <- 0
}

if (dirlenFEL == 0) {
  FELpos <- 0
  FELneg <- 0
} else if (dirlenFEL == 1 & as.numeric(attr(dirFELtab,"names")[1]) == 1) {
  FELpos <- as.numeric(dirFELtab)[1]
  FELneg <- 0
} else if (dirlenFEL == 1 & as.numeric(attr(dirFELtab,"names")[1]) == -1) {
  FELpos <- 0
  FELneg <- as.numeric(dirFELtab)[1]
} else if (dirlenFEL == 2) {
  FELneg <- as.numeric(dirFELtab)[1]
  FELpos <- as.numeric(dirFELtab)[2]
} else {
  FELpos <- 0
  FELneg <- 0
}

if (dirlenINT == 0) {
  INTpos <- 0
  INTneg <- 0
} else if (dirlenINT == 1 & as.numeric(attr(dirINTtab,"names")[1]) == 1) {
  INTpos <- as.numeric(dirINTtab)[1]
  INTneg <- 0
} else if (dirlenINT == 1 & as.numeric(attr(dirINTtab,"names")[1]) == -1) {
  INTpos <- 0
  INTneg <- as.numeric(dirINTtab)[1]
} else if (dirlenINT == 2) {
  INTneg <- as.numeric(dirINTtab)[1]
  INTpos <- as.numeric(dirINTtab)[2]
} else {
  INTpos <- 0
  INTneg <- 0
}

if (dirlenMAF == 0) {
  MAFpos <- 0
  MAFneg <- 0
} else if (dirlenMAF == 1 & as.numeric(attr(dirMAFtab,"names")[1]) == 1) {
  MAFpos <- as.numeric(dirMAFtab)[1]
  MAFneg <- 0
} else if (dirlenMAF == 1 & as.numeric(attr(dirMAFtab,"names")[1]) == -1) {
  MAFpos <- 0
  MAFneg <- as.numeric(dirMAFtab)[1]
} else if (dirlenMAF == 2) {
  MAFneg <- as.numeric(dirMAFtab)[1]
  MAFpos <- as.numeric(dirMAFtab)[2]
} else {
  MAFpos <- 0
  MAFneg <- 0
}

if (dirlenMET == 0) {
  METpos <- 0
  METneg <- 0
} else if (dirlenMET == 1 & as.numeric(attr(dirMETtab,"names")[1]) == 1) {
  METpos <- as.numeric(dirMETtab)[1]
  METneg <- 0
} else if (dirlenMET == 1 & as.numeric(attr(dirMETtab,"names")[1]) == -1) {
  METpos <- 0
  METneg <- as.numeric(dirMETtab)[1]
} else if (dirlenMET == 2) {
  METneg <- as.numeric(dirMETtab)[1]
  METpos <- as.numeric(dirMETtab)[2]
} else {
  METpos <- 0
  METneg <- 0
}

if (dirlenOTH == 0) {
  OTHpos <- 0
  OTHneg <- 0
} else if (dirlenOTH == 1 & as.numeric(attr(dirOTHtab,"names")[1]) == 1) {
  OTHpos <- as.numeric(dirOTHtab)[1]
  OTHneg <- 0
} else if (dirlenOTH == 1 & as.numeric(attr(dirOTHtab,"names")[1]) == -1) {
  OTHpos <- 0
  OTHneg <- as.numeric(dirOTHtab)[1]
} else if (dirlenOTH == 2) {
  OTHneg <- as.numeric(dirOTHtab)[1]
  OTHpos <- as.numeric(dirOTHtab)[2]
} else {
  OTHpos <- 0
  OTHneg <- 0
}

if (dirlenSIL == 0) {
  SILpos <- 0
  SILneg <- 0
} else if (dirlenSIL == 1 & as.numeric(attr(dirSILtab,"names")[1]) == 1) {
  SILpos <- as.numeric(dirSILtab)[1]
  SILneg <- 0
} else if (dirlenSIL == 1 & as.numeric(attr(dirSILtab,"names")[1]) == -1) {
  SILpos <- 0
  SILneg <- as.numeric(dirSILtab)[1]
} else if (dirlenSIL == 2) {
  SILneg <- as.numeric(dirSILtab)[1]
  SILpos <- as.numeric(dirSILtab)[2]
} else {
  SILpos <- 0
  SILneg <- 0
}

nznegdirrslts <- c(CARneg, FELneg, INTneg, MAFneg, METneg, OTHneg, SILneg)
nzposdirrslts <- c(CARpos, FELpos, INTpos, MAFpos, METpos, OTHpos, SILpos)

results.out <- data.frame(geol=geoldir.lab, nonzero=nzrslts/iter, negdir=nznegdirrslts/iter,
                          posdir=nzposdirrslts/iter)
results.out

# ER
10^median(log10(stor.mat[,"ER"]))
quantile(log10(stor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(stor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(stor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(stor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(stor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(stor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(stor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(stor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)


# soil
table(TOS.ran.subset$soil)

table(TOS.ran.subset$soil)
HgXsoil.stats.ran <- TOS.ran.subset %>%
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
HgXsoil.stats.ran
na.omit(HgXsoil.stats.ran)

# resampling loop
# storage matrix
table(TOS.lm$soil)
TOS.lm.soils.nolake <- TOS.lm[which(TOS.lm$soil != "lake"),]
table(TOS.lm.soils.nolake$soil)
nlevels <- length(table(TOS.lm.soils.nolake$soil))
stor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*nlevels)
soildir.lab <- paste("soil",attr(table(TOS.lm.soils.nolake$soil), "names"),"dir",sep="")
colnames(stor.mat) <- c("ER", "nonzerosum",
                        paste("soil",attr(table(TOS.lm.soils.nolake$soil), "names"),sep=""),soildir.lab)
head(stor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = TOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = TOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  TOS.ran.subset <- na.omit(TOS.lm.soils.nolake[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  mod1 <- lm(lHg ~ soil, data=TOS.ran.subset) # class level model
  mod.null <- lm(lHg ~ 1, data=TOS.ran.subset) # null model
  
  # coefficient boundaries
  pmmod1 <- plot_model(mod1)
  pmmod1.coef <- data.frame(soil=as.character(pmmod1[[1]]$term), lo=pmmod1[[1]]$conf.low,
                            up=pmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  pmmod1.coef$nonzero <- ifelse(contains_zero_sign(pmmod1.coef$lo, pmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  pmmod1.coef$dir <- ifelse(sign(pmmod1.coef$lo) == -1 & sign(pmmod1.coef$up) == -1, -1, 0)
  pmmod1.coef$dir <- ifelse(sign(pmmod1.coef$lo) == 1 & sign(pmmod1.coef$up) == 1, 1, pmmod1.coef$dir)
  
  # how many variables are non-zero?
  stor.mat[i,"nonzerosum"] <- sum(pmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  nzsoils <- pmmod1.coef[which(pmmod1.coef$nonzero==1),]$soil
  stor.mat[i, which(colnames(stor.mat) %in% nzsoils)] <- 1
  
  # for non-zeros, what is the direction?
  nzsoilssdir <- paste(nzsoils,"dir",sep="")
  stor.mat[i, which(colnames(stor.mat) %in% nzsoilssdir)] <- pmmod1.coef$dir[which(pmmod1.coef$dir != 0)]
  
  # AICc model comparison
  wAICc <- weight.IC(delta.IC(c(AICc(mod1),AICc(mod.null))))
  
  # store evidence ratio
  stor.mat[i,"ER"] <- wAICc[1]/wAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

# ER
median(stor.mat[,"ER"], na.rm=T)
quantile(stor.mat[,"ER"], probs=c(0.025,0.975), na.rm=T)

head(stor.mat)
soil.labs <- c("CALC","CHROMO","DERMO","FERRO","HYDRO","KANDO","KURO","ORGANO","PODO","RUDO","SODO","TENO","VERTO")
lennzCALC <- length(table(stor.mat[,"soilcalcarosol"]))
lennzCHROMO <- length(table(stor.mat[,"soilchromosol"]))
lennzDERMO <- length(table(stor.mat[,"soildermosol"]))
lennzFERRO <- length(table(stor.mat[,"soilferrosol"]))
lennzHYDRO <- length(table(stor.mat[,"soilhydrosol"]))
lennzKANDO <- length(table(stor.mat[,"soilkandosol"]))
lennzKURO <- length(table(stor.mat[,"soilkurosol"]))
lennzORGANO <- length(table(stor.mat[,"soilorganosol"]))
lennzPODO <- length(table(stor.mat[,"soilpodosol"]))
lennzRUDO <- length(table(stor.mat[,"soilrudosol"]))
lennzSODO <- length(table(stor.mat[,"soilsodosol"]))
lennzTENO <- length(table(stor.mat[,"soiltenosol"]))
lennzVERTO<- length(table(stor.mat[,"soilvertosol"]))

nzrslts <- c(ifelse(lennzCALC==1, as.numeric(table(stor.mat[,"soilcalcarosol"])), 0),
             ifelse(lennzCHROMO==1, as.numeric(table(stor.mat[,"soilchromosol"])), 0),
             ifelse(lennzDERMO==1, as.numeric(table(stor.mat[,"soildermosol"])), 0),
             ifelse(lennzFERRO==1, as.numeric(table(stor.mat[,"soilferrosol"])), 0),
             ifelse(lennzHYDRO==1, as.numeric(table(stor.mat[,"soilhydrosol"])), 0),
             ifelse(lennzKANDO==1, as.numeric(table(stor.mat[,"soilkandosol"])), 0),
             ifelse(lennzKURO==1, as.numeric(table(stor.mat[,"soilkurosol"])), 0),
             ifelse(lennzORGANO==1, as.numeric(table(stor.mat[,"soilorganosol"])), 0),
             ifelse(lennzPODO==1, as.numeric(table(stor.mat[,"soilpodosol"])), 0),
             ifelse(lennzRUDO==1, as.numeric(table(stor.mat[,"soilrudosol"])), 0),
             ifelse(lennzSODO==1, as.numeric(table(stor.mat[,"soilsodosol"])), 0),
             ifelse(lennzTENO==1, as.numeric(table(stor.mat[,"soiltenosol"])), 0),
             ifelse(lennzVERTO==1, as.numeric(table(stor.mat[,"soilvertosol"])), 0))

dirCALCtab <- table(stor.mat[,"soilcalcarosol"])
dirCHROMOtab <- table(stor.mat[,"soilchromosol"])
dirDERMOtab <- table(stor.mat[,"soildermosol"])
dirFERROtab <- table(stor.mat[,"soilferrosol"])
dirHYDROtab <- table(stor.mat[,"soilhydrosol"])
dirKANDOtab <- table(stor.mat[,"soilkandosol"])
dirKUROtab <- table(stor.mat[,"soilkurosol"])
dirORGANOTab <- table(stor.mat[,"soilorganosol"])
dirPODOtab <- table(stor.mat[,"soilpodosol"])
dirRUDOtab <- table(stor.mat[,"soilrudosol"])
dirSODOtab <- table(stor.mat[,"soilsodosol"])
dirTENOTab <- table(stor.mat[,"soiltenosol"])
dirVERTOTab <- table(stor.mat[,"soilvertosol"])

dirlenCALC <- length(dirCALCtab)
dirlenCHROMO <- length(dirCHROMOtab)
dirlenDERMO <- length(dirDERMOtab)
dirlenFERRO <- length(dirFERROtab)
dirlenHYDRO <- length(dirHYDROtab)
dirlenKANDO <- length(dirKANDOtab)
dirlenKURO <- length(dirKUROtab)
dirlenORGANO <- length(dirORGANOTab)
dirlenPODO <- length(dirPODOtab)
dirlenRUDO <- length(dirRUDOtab)
dirlenSODO <- length(dirSODOtab)
dirlenTENO <- length(dirTENOTab)
dirlenVERTO <- length(dirVERTOTab)

if (dirlenCALC == 0) {
  CALCpos <- 0
  CALCneg <- 0
} else if (dirlenCALC == 1 & as.numeric(attr(dirCALCtab,"names")[1]) == 1) {
  CALCpos <- as.numeric(dirCALCtab)[1]
  CALCneg <- 0
} else if (dirlenCALC == 1 & as.numeric(attr(dirCALCtab,"names")[1]) == -1) {
  CALCpos <- 0
  CALCneg <- as.numeric(dirCALCtab)[1]
} else if (dirlenCALC == 2) {
  CALCneg <- as.numeric(dirCALCtab)[1]
  CALCpos <- as.numeric(dirCALCtab)[2]
} else {
  CALCpos <- 0
  CALCneg <- 0
}

if (dirlenCHROMO == 0) {
  CHROMOpos <- 0
  CHROMOneg <- 0
} else if (dirlenCHROMO == 1 & as.numeric(attr(dirCHROMOtab,"names")[1]) == 1) {
  CHROMOpos <- as.numeric(dirCHROMOtab)[1]
  CHROMOneg <- 0
} else if (dirlenCHROMO == 1 & as.numeric(attr(dirCHROMOtab,"names")[1]) == -1) {
  CHROMOpos <- 0
  CHROMOneg <- as.numeric(dirCHROMOtab)[1]
} else if (dirlenCHROMO == 2) {
  CHROMOneg <- as.numeric(dirCHROMOtab)[1]
  CHROMOpos <- as.numeric(dirCHROMOtab)[2]
} else {
  CHROMOpos <- 0
  CHROMOneg <- 0
}

if (dirlenDERMO == 0) {
  DERMOpos <- 0
  DERMOneg <- 0
} else if (dirlenDERMO == 1 & as.numeric(attr(dirDERMOtab,"names")[1]) == 1) {
  DERMOpos <- as.numeric(dirDERMOtab)[1]
  DERMOneg <- 0
} else if (dirlenDERMO == 1 & as.numeric(attr(dirDERMOtab,"names")[1]) == -1) {
  DERMOpos <- 0
  DERMOneg <- as.numeric(dirDERMOtab)[1]
} else if (dirlenDERMO == 2) {
  DERMOneg <- as.numeric(dirDERMOtab)[1]
  DERMOpos <- as.numeric(dirDERMOtab)[2]
} else {
  DERMOpos <- 0
  DERMOneg <- 0
}

if (dirlenFERRO == 0) {
  FERROpos <- 0
  FERROneg <- 0
} else if (dirlenFERRO == 1 & as.numeric(attr(dirFERROtab,"names")[1]) == 1) {
  FERROpos <- as.numeric(dirFERROtab)[1]
  FERROneg <- 0
} else if (dirlenFERRO == 1 & as.numeric(attr(dirFERROtab,"names")[1]) == -1) {
  FERROpos <- 0
  FERROneg <- as.numeric(dirFERROtab)[1]
} else if (dirlenFERRO == 2) {
  FERROneg <- as.numeric(dirFERROtab)[1]
  FERROpos <- as.numeric(dirFERROtab)[2]
} else {
  FERROpos <- 0
  FERROneg <- 0
}

if (dirlenHYDRO == 0) {
  HYDROpos <- 0
  HYDROneg <- 0
} else if (dirlenHYDRO == 1 & as.numeric(attr(dirHYDROtab,"names")[1]) == 1) {
  HYDROpos <- as.numeric(dirHYDROtab)[1]
  HYDROneg <- 0
} else if (dirlenHYDRO == 1 & as.numeric(attr(dirHYDROtab,"names")[1]) == -1) {
  HYDROpos <- 0
  HYDROneg <- as.numeric(dirHYDROtab)[1]
} else if (dirlenHYDRO == 2) {
  HYDROneg <- as.numeric(dirHYDROtab)[1]
  HYDROpos <- as.numeric(dirHYDROtab)[2]
} else {
  HYDROpos <- 0
  HYDROneg <- 0
}

if (dirlenKANDO == 0) {
  KANDOpos <- 0
  KANDOneg <- 0
} else if (dirlenKANDO == 1 & as.numeric(attr(dirKANDOtab,"names")[1]) == 1) {
  KANDOpos <- as.numeric(dirKANDOtab)[1]
  KANDOneg <- 0
} else if (dirlenKANDO == 1 & as.numeric(attr(dirKANDOtab,"names")[1]) == -1) {
  KANDOpos <- 0
  KANDOneg <- as.numeric(dirKANDOtab)[1]
} else if (dirlenKANDO == 2) {
  KANDOneg <- as.numeric(dirKANDOtab)[1]
  KANDOpos <- as.numeric(dirKANDOtab)[2]
} else {
  KANDOpos <- 0
  KANDOneg <- 0
}

if (dirlenKURO == 0) {
  KUROpos <- 0
  KUROneg <- 0
} else if (dirlenKURO == 1 & as.numeric(attr(dirKUROtab,"names")[1]) == 1) {
  KUROpos <- as.numeric(dirKUROtab)[1]
  KUROneg <- 0
} else if (dirlenKURO == 1 & as.numeric(attr(dirKUROtab,"names")[1]) == -1) {
  KUROpos <- 0
  KUROneg <- as.numeric(dirKUROtab)[1]
} else if (dirlenKURO == 2) {
  KUROneg <- as.numeric(dirKUROtab)[1]
  KUROpos <- as.numeric(dirKUROtab)[2]
} else {
  KUROpos <- 0
  KUROneg <- 0
}

if (dirlenORGANO == 0) {
  ORGANOpos <- 0
  ORGANOneg <- 0
} else if (dirlenORGANO == 1 & as.numeric(attr(dirORGANOTab,"names")[1]) == 1) {
  ORGANOpos <- as.numeric(dirORGANOTab)[1]
  ORGANOneg <- 0
} else if (dirlenORGANO == 1 & as.numeric(attr(dirORGANOTab,"names")[1]) == -1) {
  ORGANOpos <- 0
  ORGANOneg <- as.numeric(dirORGANOTab)[1]
} else if (dirlenORGANO == 2) {
  ORGANOneg <- as.numeric(dirORGANOTab)[1]
  ORGANOpos <- as.numeric(dirORGANOTab)[2]
} else {
  ORGANOpos <- 0
  ORGANOneg <- 0
}

if (dirlenPODO == 0) {
  PODOpos <- 0
  PODOneg <- 0
} else if (dirlenPODO == 1 & as.numeric(attr(dirPODOtab,"names")[1]) == 1) {
  PODOpos <- as.numeric(dirPODOtab)[1]
  PODOneg <- 0
} else if (dirlenPODO == 1 & as.numeric(attr(dirPODOtab,"names")[1]) == -1) {
  PODOpos <- 0
  PODOneg <- as.numeric(dirPODOtab)[1]
} else if (dirlenPODO == 2) {
  PODOneg <- as.numeric(dirPODOtab)[1]
  PODOpos <- as.numeric(dirPODOtab)[2]
} else {
  PODOpos <- 0
  PODOneg <- 0
}

if (dirlenRUDO == 0) {
  RUDOpos <- 0
  RUDOneg <- 0
} else if (dirlenRUDO == 1 & as.numeric(attr(dirRUDOtab,"names")[1]) == 1) {
  RUDOpos <- as.numeric(dirRUDOtab)[1]
  RUDOneg <- 0
} else if (dirlenRUDO == 1 & as.numeric(attr(dirRUDOtab,"names")[1]) == -1) {
  RUDOpos <- 0
  RUDOneg <- as.numeric(dirRUDOtab)[1]
} else if (dirlenRUDO == 2) {
  RUDOneg <- as.numeric(dirRUDOtab)[1]
  RUDOpos <- as.numeric(dirRUDOtab)[2]
} else {
  RUDOpos <- 0
  RUDOneg <- 0
}

if (dirlenSODO == 0) {
  SODOpos <- 0
  SODOneg <- 0
} else if (dirlenSODO == 1 & as.numeric(attr(dirSODOtab,"names")[1]) == 1) {
  SODOpos <- as.numeric(dirSODOtab)[1]
  SODOneg <- 0
} else if (dirlenSODO == 1 & as.numeric(attr(dirSODOtab,"names")[1]) == -1) {
  SODOpos <- 0
  SODOneg <- as.numeric(dirSODOtab)[1]
} else if (dirlenSODO == 2) {
  SODOneg <- as.numeric(dirSODOtab)[1]
  SODOpos <- as.numeric(dirSODOtab)[2]
} else {
  SODOpos <- 0
  SODOneg <- 0
}

if (dirlenTENO == 0) {
  TENOpos <- 0
  TENOneg <- 0
} else if (dirlenTENO == 1 & as.numeric(attr(dirTENOTab,"names")[1]) == 1) {
  TENOpos <- as.numeric(dirTENOTab)[1]
  TENOneg <- 0
} else if (dirlenTENO == 1 & as.numeric(attr(dirTENOTab,"names")[1]) == -1) {
  TENOpos <- 0
  TENOneg <- as.numeric(dirTENOTab)[1]
} else if (dirlenTENO == 2) {
  TENOneg <- as.numeric(dirTENOTab)[1]
  TENOpos <- as.numeric(dirTENOTab)[2]
} else {
  TENOpos <- 0
  TENOneg <- 0
}

if (dirlenVERTO == 0) {
  VERTOpos <- 0
  VERTOneg <- 0
} else if (dirlenVERTO == 1 & as.numeric(attr(dirVERTOTab,"names")[1]) == 1) {
  VERTOpos <- as.numeric(dirVERTOTab)[1]
  VERTOneg <- 0
} else if (dirlenVERTO == 1 & as.numeric(attr(dirVERTOTab,"names")[1]) == -1) {
  VERTOpos <- 0
  VERTOneg <- as.numeric(dirVERTOTab)[1]
} else if (dirlenVERTO == 2) {
  VERTOneg <- as.numeric(dirVERTOTab)[1]
  VERTOpos <- as.numeric(dirVERTOTab)[2]
} else {
  VERTOpos <- 0
  VERTOneg <- 0
}

nznegdirrslts <- c(CALCneg, CHROMOneg, DERMOneg, FERROneg, HYDROneg, KANDOneg, KUROneg, ORGANOneg, PODOneg, RUDOneg, SODOneg, TENOneg, VERTOneg)
nzposdirrslts <- c(CALCpos, CHROMOpos, DERMOpos, FERROpos, HYDROpos, KANDOpos, KUROpos, ORGANOpos, PODOpos, RUDOpos, SODOpos, TENOpos, VERTOpos)

results.out <- data.frame(soil=soildir.lab, nonzero=nzrslts/iter, negdir=nznegdirrslts/iter,
                          posdir=nzposdirrslts/iter)
results.out

# ER
10^median(log10(stor.mat[,"ER"]))
quantile(log10(stor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(stor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(stor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(stor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(stor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(stor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(stor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(stor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)





# landuse
table(TOS.ran.subset$landuse)

table(TOS.ran.subset$landuse)
HgXlanduse.stats.ran <- TOS.ran.subset %>%
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
HgXlanduse.stats.ran
na.omit(HgXlanduse.stats.ran)


# resampling loop
# storage matrix
table(TOS.lm$landuse)
nlevels <- length(table(TOS.lm$landuse))
stor.mat <- matrix(data=NA, nrow=iter, ncol=2+2*nlevels)
landusedir.lab <- paste("landuse",attr(table(TOS.lm$landuse), "names"),"dir",sep="")
colnames(stor.mat) <- c("ER", "nonzerosum",
                        paste("landuse",attr(table(TOS.lm$landuse), "names"),sep=""),landusedir.lab)
head(stor.mat)

for (i in 1:iter) {
  # generate random coords
  coords.ran <- select_distant_points(df = TOS.lm.coords, min_dist = min.dist,
                                      dist_matrix = TOS.lm_haversine_matrix, x_col = "LON", y_col = "LAT")
  
  # subset random points for data summaries
  TOS.ran.subset <- na.omit(TOS.lm[as.numeric(row.names(coords.ran)), ])
  
  # fit linear models
  mod1 <- lm(lHg ~ landuse, data=TOS.ran.subset) # class level model
  mod.null <- lm(lHg ~ 1, data=TOS.ran.subset) # null model
  
  # coefficient boundaries
  pmmod1 <- plot_model(mod1)
  pmmod1.coef <- data.frame(landuse=as.character(pmmod1[[1]]$term), lo=pmmod1[[1]]$conf.low,
                            up=pmmod1[[1]]$conf.high)
  
  # determine if zero is in the confidence interval
  pmmod1.coef$nonzero <- ifelse(contains_zero_sign(pmmod1.coef$lo, pmmod1.coef$up)==F, 1, 0)
  
  # determine if negative or positive relationship
  pmmod1.coef$dir <- ifelse(sign(pmmod1.coef$lo) == -1 & sign(pmmod1.coef$up) == -1, -1, 0)
  pmmod1.coef$dir <- ifelse(sign(pmmod1.coef$lo) == 1 & sign(pmmod1.coef$up) == 1, 1, pmmod1.coef$dir)
  
  # how many variables are non-zero?
  stor.mat[i,"nonzerosum"] <- sum(pmmod1.coef$nonzero)
  
  # which category levels are non-zero?
  nzlanduses <- pmmod1.coef[which(pmmod1.coef$nonzero==1),]$landuse
  stor.mat[i, which(colnames(stor.mat) %in% nzlanduses)] <- 1
  
  # for non-zeros, what is the direction?
  nzlandusessdir <- paste(nzlanduses,"dir",sep="")
  stor.mat[i, which(colnames(stor.mat) %in% nzlandusessdir)] <- pmmod1.coef$dir[which(pmmod1.coef$dir != 0)]
  
  # AICc model comparison
  wAICc <- weight.IC(delta.IC(c(AICc(mod1),AICc(mod.null))))
  
  # store evidence ratio
  stor.mat[i,"ER"] <- wAICc[1]/wAICc[2]
  
  if (i %% itdiv == 0) print(paste("iter = ", i, sep=""))
}

# ER
median(stor.mat[,"ER"])
quantile(stor.mat[,"ER"], probs=c(0.025,0.975))

head(stor.mat)
landuse.labs <- c("CN","INT","PDA","PIA","PRN","WAT","WA")
lennzCN <- length(table(stor.mat[,"landuseconservation/natural"]))
lennzINT <- length(table(stor.mat[,"landuseintensive"]))
lennzPDA <- length(table(stor.mat[,"landuseproduction-dryland agr"]))
lennzPIA <- length(table(stor.mat[,"landuseproduction-irrigated agr"]))
lennzPRN <- length(table(stor.mat[,"landuseproduction-relatively natural"]))
lennzWAT <- length(table(stor.mat[,"landusewater"]))

nzrslts <- c(ifelse(lennzCN==1, as.numeric(table(stor.mat[,"landuseconservation/natural"])), 0),
             ifelse(lennzINT==1, as.numeric(table(stor.mat[,"landuseintensive"])), 0),
             ifelse(lennzPDA==1, as.numeric(table(stor.mat[,"landuseproduction-dryland agr"])), 0),
             ifelse(lennzPIA==1, as.numeric(table(stor.mat[,"landuseproduction-irrigated agr"])), 0),
             ifelse(lennzPRN==1, as.numeric(table(stor.mat[,"landuseproduction-relatively natural"])), 0),
             ifelse(lennzWAT==1, as.numeric(table(stor.mat[,"landusewater"])), 0))


dirCNtab <- table(stor.mat[,"landuseconservation/natural"])
dirINTtab <- table(stor.mat[,"landuseintensive"])
dirPDAtab <- table(stor.mat[,"landuseproduction-dryland agr"])
dirPIAtab <- table(stor.mat[,"landuseproduction-irrigated agr"])
dirPRNtab <- table(stor.mat[,"landuseproduction-relatively natural"])
dirWATtab <- table(stor.mat[,"landusewater"])

dirlenCN <- length(dirCNtab)
dirlenINT <- length(dirINTtab)
dirlenPDA <- length(dirPDAtab)
dirlenPIA <- length(dirPIAtab)
dirlenPRN <- length(dirPRNtab)
dirlenWAT <- length(dirWATtab)


if (dirlenCN == 0) {
  CNpos <- 0
  CNneg <- 0
} else if (dirlenCN == 1 & as.numeric(attr(dirCNtab,"names")[1]) == 1) {
  CNpos <- as.numeric(dirCNtab)[1]
  CNneg <- 0
} else if (dirlenCN == 1 & as.numeric(attr(dirCNtab,"names")[1]) == -1) {
  CNpos <- 0
  CNneg <- as.numeric(dirCNtab)[1]
} else if (dirlenCN == 2) {
  CNneg <- as.numeric(dirCNtab)[1]
  CNpos <- as.numeric(dirCNtab)[2]
} else {
  CNpos <- 0
  CNneg <- 0
}

if (dirlenINT == 0) {
  INTpos <- 0
  INTneg <- 0
  } else if (dirlenINT == 1 & as.numeric(attr(dirINTtab,"names")[1]) == 1) {
    INTpos <- as.numeric(dirINTtab)[1]
    INTneg <- 0
  } else if (dirlenINT == 1 & as.numeric(attr(dirINTtab,"names")[1]) == -1) {
    INTpos <- 0
    INTneg <- as.numeric(dirINTtab)[1]
  } else if (dirlenINT == 2) {
    INTneg <- as.numeric(dirINTtab)[1]
    INTpos <- as.numeric(dirINTtab)[2]
  } else {
    INTpos <- 0
    INTneg <- 0
}

if (dirlenPDA == 0) {
  PDApos <- 0
  PDAneg <- 0
  } else if (dirlenPDA == 1 & as.numeric(attr(dirPDAtab,"names")[1]) == 1) {
    PDApos <- as.numeric(dirPDAtab)[1]
    PDAneg <- 0
  } else if (dirlenPDA == 1 & as.numeric(attr(dirPDAtab,"names")[1]) == -1) {
    PDApos <- 0
    PDAneg <- as.numeric(dirPDAtab)[1]
  } else if (dirlenPDA == 2) {
    PDAneg <- as.numeric(dirPDAtab)[1]
    PDApos <- as.numeric(dirPDAtab)[2]
  } else {
    PDApos <- 0
    PDAneg <- 0
}

if (dirlenPIA == 0) {
  PIApos <- 0
  PIAneg <- 0
  } else if (dirlenPIA == 1 & as.numeric(attr(dirPIAtab,"names")[1]) == 1) {
    PIApos <- as.numeric(dirPIAtab)[1]
    PIAneg <- 0
  } else if (dirlenPIA == 1 & as.numeric(attr(dirPIAtab,"names")[1]) == -1) {
    PIApos <- 0
    PIAneg <- as.numeric(dirPIAtab)[1]
  } else if (dirlenPIA == 2) {
    PIAneg <- as.numeric(dirPIAtab)[1]
    PIApos <- as.numeric(dirPIAtab)[2]
  } else {
    PIApos <- 0
    PIAneg <- 0
}

if (dirlenPRN == 0) {
  PRNpos <- 0
  PRNneg <- 0
  } else if (dirlenPRN == 1 & as.numeric(attr(dirPRNtab,"names")[1]) == 1) {
    PRNpos <- as.numeric(dirPRNtab)[1]
    PRNneg <- 0
  } else if (dirlenPRN == 1 & as.numeric(attr(dirPRNtab,"names")[1]) == -1) {
    PRNpos <- 0
    PRNneg <- as.numeric(dirPRNtab)[1]
  } else if (dirlenPRN == 2) {
    PRNneg <- as.numeric(dirPRNtab)[1]
    PRNpos <- as.numeric(dirPRNtab)[2]
  } else {
    PRNpos <- 0
    PRNneg <- 0
}

if (dirlenWAT == 0) {
  WATpos <- 0
  WATneg <- 0
  } else if (dirlenWAT == 1 & as.numeric(attr(dirWATtab,"names")[1]) == 1) {
    WATpos <- as.numeric(dirWATtab)[1]
    WATneg <- 0
  } else if (dirlenWAT == 1 & as.numeric(attr(dirWATtab,"names")[1]) == -1) {
    WATpos <- 0
    WATneg <- as.numeric(dirWATtab)[1]
  } else if (dirlenWAT == 2) {
    WATneg <- as.numeric(dirWATtab)[1]
    WATpos <- as.numeric(dirWATtab)[2]
  } else {
    WATpos <- 0
    WATneg <- 0
}

nznegdirrslts <- c(CNneg, INTneg, PDAneg, PIAneg, PRNneg, WATneg)
nzposdirrslts <- c(CNpos, INTpos, PDApos, PIApos, PRNpos, WATpos)

results.out <- data.frame(landuse=landusedir.lab, nonzero=nzrslts/iter, negdir=nznegdirrslts/iter,
                          posdir=nzposdirrslts/iter)
results.out

# ER
10^median(log10(stor.mat[,"ER"]))
quantile(log10(stor.mat[,"ER"]), probs=c(0.025,0.975))
10^quantile(log10(stor.mat[,"ER"]), probs=c(0.025,0.975))
quantile(log10(stor.mat[,"ER"]), probs=c(0.1,0.9))
10^quantile(log10(stor.mat[,"ER"]), probs=c(0.1,0.9))

hist(log10(stor.mat[,"ER"]), main="", xlab="log10 ER")
abline(v=median(log10(stor.mat[,"ER"])), col="red", lwd=2, lty=2)
abline(v=quantile(log10(stor.mat[,"ER"]), probs=0.1), col="blue", lwd=2, lty=2)
abline(v=quantile(log10(stor.mat[,"ER"]), probs=0.9), col="blue", lwd=2, lty=2)


# by biome
WWFecoregions <- vect("wwf_terr_ecos.shp") 
WWFecoregions.crop <- crop(WWFecoregions, aus.ext)
WWFecoregions.sf <- st_as_sf(WWFecoregions.crop)

# join to biome key
biome.key
biome <- inner_join(WWFecoregions.sf, biome.key, by="BIOME")
biome

lHgXbiome.stats <- TOS.dat %>%
  group_by(biome) %>%
  summarise(
    mean = mean(log10(HgCOMP), na.rm = TRUE),
    median = median(log10(HgCOMP), na.rm = TRUE), 
    var = var(log10(HgCOMP), na.rm = TRUE),
    sd = sd(log10(HgCOMP), na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(log10(HgCOMP), probs=0.975, na.rm = TRUE),
    lower = quantile(log10(HgCOMP), probs=0.025, na.rm = TRUE),
    n = n()
  )
lHgXbiome.stats

# add mean Hg to biome shapefile
biome$lHg <- 0

table(biome$abbr)
biome$lHg <- ifelse(biome$abbr == "DXS",
                    lHgXbiome.stats[lHgXbiome.stats$biome=="DXS",]$mean, 0)
biome$lHg <- ifelse(biome$abbr == "MFW",
                    lHgXbiome.stats[lHgXbiome.stats$biome=="MFW",]$mean, biome$lHg)
biome$lHg <- ifelse(biome$abbr == "TBMF",
                    lHgXbiome.stats[lHgXbiome.stats$biome=="TBMF",]$mean, biome$lHg)
biome$lHg <- ifelse(biome$abbr == "TGSS",
                    lHgXbiome.stats[lHgXbiome.stats$biome=="TGSS",]$mean, biome$lHg)
biome$lHg <- ifelse(biome$abbr == "TSGSS",
                    lHgXbiome.stats[lHgXbiome.stats$biome=="TSGSS",]$mean, biome$lHg)
biome$lHg <- ifelse(biome$abbr == "TSMBF",
                    lHgXbiome.stats[lHgXbiome.stats$biome=="TSMBF",]$mean, biome$lHg)

names(biome)
summary(biome[,"lHg"])
class(biome)

# write
st_write(biome, "biomeHg.shp", append=F)

## reimport as spatVector to dissolve
biome.sv <- vect("biomeHg.shp")
biome.diss <- aggregate(biome.sv, "abbr")
biome.sf <- st_as_sf(biome.diss)
names(biome.sf)

# write
st_write(biome.sf, "biomeHg.shp", append=F)

# 3D plot with ggplot2 and plot_gg
biome3Dplot = ggplot(data=biome.sf, (aes(fill = mean_lHg))) +
  scale_fill_viridis_c(direction=-1, alpha = .7, limits = c(1, 1.7)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.ticks = element_blank(), axis.text = element_blank(),
                     legend.title = element_blank(), legend.position="none") + 
  geom_sf()
biome3Dplot
plot_gg(biome3Dplot, multicore = TRUE, width=5, height=5,
        scale=200, windowsize=c(1400,866), zoom = 0.4, phi = 30, shadow_intensity = 0.3,
        offset_edges = T, raytrace = T, triangulate = T, anglebreaks=seq(50,60,0.1))
render_snapshot()


# by lithology
lith.orig <- st_read("GeologicUnitPolygons1M.shp")
names(lith.orig)
attr(table(lith.orig$LITHOLOGY), 'names')

lith.key
lith <- inner_join(lith.orig, lith.key, by="LITHOLOGY")
head(lith)
class(lith)

lHgXlith.stats <- TOS.dat %>%
  group_by(lithgrp) %>%
  summarise(
    mean = mean(log10(HgCOMP), na.rm = TRUE),
    median = median(log10(HgCOMP), na.rm = TRUE), 
    var = var(log10(HgCOMP), na.rm = TRUE),
    sd = sd(log10(HgCOMP), na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(log10(HgCOMP), probs=0.975, na.rm = TRUE),
    lower = quantile(log10(HgCOMP), probs=0.025, na.rm = TRUE),
    n = n()
  )
lHgXlith.stats

# add mean Hg to biome shapefile
lith$lHg <- 0

table(lith$LITHGRP)
lith$lHg <- ifelse(lith$LITHGRP == "CAR",
                    lHgXlith.stats[lHgXlith.stats$lithgrp=="CAR",]$mean, 0)
lith$lHg <- ifelse(lith$LITHGRP == "FEL",
                    lHgXlith.stats[lHgXlith.stats$lithgrp=="FEL",]$mean, lith$lHg)
lith$lHg <- ifelse(lith$LITHGRP == "INT",
                    lHgXlith.stats[lHgXlith.stats$lithgrp=="INT",]$mean, lith$lHg)
lith$lHg <- ifelse(lith$LITHGRP == "MAF",
                    lHgXlith.stats[lHgXlith.stats$lithgrp=="MAF",]$mean, lith$lHg)
lith$lHg <- ifelse(lith$LITHGRP == "MET",
                    lHgXlith.stats[lHgXlith.stats$lithgrp=="MET",]$mean, lith$lHg)
lith$lHg <- ifelse(lith$LITHGRP == "OTH",
                    lHgXlith.stats[lHgXlith.stats$lithgrp=="OTH",]$mean, lith$lHg)
lith$lHg <- ifelse(lith$LITHGRP == "SIL",
                    lHgXlith.stats[lHgXlith.stats$lithgrp=="SIL",]$mean, lith$lHg)

# save to shapefile
st_write(lith, "lithHg.shp", append=F)

# re-import as spatVector
lithHg <- vect("lithHg.shp")
names(lithHg)

# aggregate spatVector (dissolve)
# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

lithHg.diss <- aggregate(lithHg, fun='mean', by="LITHGRP", cores=cl, na.rm=T)
parallel::stopCluster(cl = local.cluster)

lithHg.diss
names(lithHg.diss)
lithHg.diss$mean_lHg

# convert back to sf object
lithHg.sf <- st_as_sf(lithHg.diss)
lithHg.sf$mean_lHg
str(lithHg.sf)
str(lithHg.sf$mean_lHg)

st_write(lithHg.sf, "lithHgDiss.shp", append=F)

# change lHg to factor to avoid colour-scale errors
lithHg.sf$lHgfact <- as.factor(lithHg.sf$mean_lHg)

# 3D plot with ggplot2 and plot_gg
lith3Dplot = ggplot(data=lithHg.sf, (aes(fill = lHgfact))) +
  scale_fill_viridis_d(direction=-1, alpha = .7) +
  scale_y_continuous(limits = c(1, 1.7)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.ticks = element_blank(), axis.text = element_blank(),
                     legend.title = element_blank(), legend.position="none") + 
  geom_sf()
lith3Dplot
plot_gg(lith3Dplot, multicore = TRUE, width=5, height=5,
        scale=200, windowsize=c(1400,866), zoom = 0.4, phi = 30, shadow_intensity = 0.3,
        offset_edges = T, raytrace = T, triangulate = T, anglebreaks=seq(50,60,0.1))
render_snapshot()



## by land use
lnduse.orig <- rast('NLUM_v7_250_ALUMV8_2020_21_alb.tif')
names(lnduse.orig)

# convert to vector
lnduse.sv <- as.polygons(lnduse.orig, aggregate=T, na.rm=T)
names(lnduse.sv)
lnduse.orig$PRIMV8N

# convert to sf object
lnduse.sf <- st_as_sf(lnduse.sv)
names(lnduse.sf)
table(lnduse.sf$TERTV8)

# identify only first-order land-use categories
lnduse.sf$plu <- as.numeric(substr(lnduse.sf$TERTV8, 1, 1))
table(lnduse.sf$plu)

# join to lu.key
lu.key
lnduse <- inner_join(lnduse.sf, lu.key, by="plu")
head(lnduse)
class(lnduse)

lHgXlu.stats <- TOS.dat %>%
  group_by(landuse) %>%
  summarise(
    mean = mean(log10(HgCOMP), na.rm = TRUE),
    median = median(log10(HgCOMP), na.rm = TRUE), 
    var = var(log10(HgCOMP), na.rm = TRUE),
    sd = sd(log10(HgCOMP), na.rm = TRUE),
    se = sd/sqrt(n()),
    upper = quantile(log10(HgCOMP), probs=0.975, na.rm = TRUE),
    lower = quantile(log10(HgCOMP), probs=0.025, na.rm = TRUE),
    n = n()
  )
lHgXlu.stats

# add mean Hg to landuse
lnduse$lHg <- 0

table(lnduse$plu)
lnduse$lHg <- ifelse(lnduse$plu == 1,
                   lHgXlu.stats[lHgXlu.stats$landuse=="conservation/natural",]$mean, 0)
lnduse$lHg <- ifelse(lnduse$plu == 2,
                     lHgXlu.stats[lHgXlu.stats$landuse=="intensive",]$mean, lnduse$lHg)
lnduse$lHg <- ifelse(lnduse$plu == 3,
                     lHgXlu.stats[lHgXlu.stats$landuse=="production-dryland agr",]$mean, lnduse$lHg)
lnduse$lHg <- ifelse(lnduse$plu == 4,
                     lHgXlu.stats[lHgXlu.stats$landuse=="production-irrigated agr",]$mean, lnduse$lHg)
lnduse$lHg <- ifelse(lnduse$plu == 5,
                     lHgXlu.stats[lHgXlu.stats$landuse=="production-relatively natural",]$mean, lnduse$lHg)
lnduse$lHg <- ifelse(lnduse$plu == 6,
                     lHgXlu.stats[lHgXlu.stats$landuse=="water",]$mean, lnduse$lHg)

# save to shapefile
st_write(lnduse, "lnduseHg.shp", append=F)

# re-import as spatVector
lnduseHg <- vect("lnduseHg.shp")
names(lnduseHg)

# aggregate spatVector (dissolve)
# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

lithHg.diss <- aggregate(lithHg, fun='mean', by="LITHGRP", cores=cl, na.rm=T)
parallel::stopCluster(cl = local.cluster)

lithHg.diss
names(lithHg.diss)
lithHg.diss$mean_lHg

# convert back to sf object
lithHg.sf <- st_as_sf(lithHg.diss)
lithHg.sf$mean_lHg
str(lithHg.sf)
str(lithHg.sf$mean_lHg)

st_write(lithHg.sf, "lithHgDiss.shp", append=F)


## sugarcane https://www.abs.gov.au/statistics/industry/agriculture/sugarcane-experimental-regional-estimates-using-new-data-sources-and-methods/latest-release#data-downloads

## Hg measures from geochem
scXgeochem <- read.csv("geochemXsugarcane.csv", header=T)
head(scXgeochem)
scXgeochemProdgtZero <- scXgeochem[which(scXgeochem$sugarcane_scp2021 >= 0 & is.na(scXgeochem$sugarcane_scp2021)==F),]
names(scXgeochemProdgtZero)

# average by SA2
scXgeochemProdgtZeroSA2 <- scXgeochemProdgtZero %>%
  group_by(SA2_CODE21) %>%
  summarise(
    meanHg = mean(Hg, na.rm=T),
    seHg = sd(Hg, na.rm=T)/sqrt(n()),
    loHg = quantile(Hg, probs=0.025, na.rm=T),
    upHg = quantile(Hg, probs=0.975, na.rm=T),
    meanscp2021 = mean(sugarcane_scp2021, na.rm=T),
    sescp2021 = sd(sugarcane_scp2021, na.rm=T)/sqrt(n()),
    loscp2021 = quantile(sugarcane_scp2021, probs=0.025, na.rm=T),
    upscp2021 = quantile(sugarcane_scp2021, probs=0.975, na.rm=T)
)
scXgeochemProdgtZeroSA2

plot(log10(scXgeochemProdgtZeroSA2$meanscp2021), log10(scXgeochemProdgtZeroSA2$meanHg), xlab="log10 sugarcane production 2020/2021 (t)",
     ylab="log10 Hg (ng/g)", pch=19)
abline(lm(log10(scXgeochemProdgtZeroSA2$meanHg) ~ log10(scXgeochemProdgtZeroSA2$meanscp2021)), col="red", lty=2)

ggplot(data=scXgeochemProdgtZeroSA2, aes(x=log10(meanscp2021), y=log10(meanHg))) +
  geom_point() +
  
  # add y error bars
  geom_errorbar(aes(ymin=log10(meanHg-seHg), ymax=log10(meanHg+seHg)), width=.1) +
  
  geom_smooth(method="lm", se=T) +
  xlab("log10 sugarcane production 2020/2021 (t)") +
  ylab("log10 mean measured Hg (ng/g)")

## Hg predicted from RF model
scpoly <- vect("sugarcane.shp")
scpoly
names(scpoly)
scpoly$sugarcan_4 # 2020-2021 production in t

lHg.rst <- rast("HgPredSpatRFlog10.tif")
lHg.rst

# overlay scpoly on lHg
scHg <- terra::extract(lHg.rst, scpoly, fun='mean', na.rm=T)
scHg

scpoly$lHgPred <- scHg$prediction
scpoly

plot(log10(scpoly$sugarcan_4), scpoly$lHgPred, pch=19,
     xlab="log10 sugarcane production 2020/2021 (t)", ylab="log10 mean predicted Hg (ng/g)")
abline(lm(scpoly$lHgPred ~ log10(scpoly$sugarcan_4)), col="red", lty=2)

plotpreddat <- data.frame(scprod=scpoly$sugarcan_4, lHg=scpoly$lHgPred)
ggplot(data=plotpreddat, aes(x=log10(scprod), y=lHg)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  xlab("log10 sugarcane production 2020/2021 (t)") +
  ylab("log10 mean predicted Hg (ng/g)")

                             


###################################################### 
## multiscale spatial analysis of multivariate data ##
######################################################

# remove NAs
colnames(TOS.test)
TOS.noNA <- na.omit(TOS.test[,c("LON","LAT","lgs","lgtclay","pH15","lelecond","lLOI","prescott","soilH20",
                            "llai","soillN","soillP","lAl","lMn","lCu","lPb","lSb","lNi","lKThU","ferric2",
                            "Alox", "Feox", "Pox")])
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

# Moran's eigenvector map
mem.gab <- mem(listwgab)
mem.gab
class(mem.gab)
names(attributes(mem.gab))
plot(mem.gab[,c(1, 5, 10, 20, 30, 40, 50, 60, 70)], SpORcoords = mxy)
s.value(mxy, mem.gab[,c(1, 5, 10, 20, 30, 40, 50, 60, 70)], symbol = "circle", ppoint.cex = 0.6)
barplot(attr(mem.gab, "values"), 
        main = "eigenvalues of the spatial weighting matrix", cex.main = 0.7)
moranI <- moran.randtest(mem.gab, listwgab, 99)
moranI
head(attr(mem.gab, "values") / moranI$obs)

# environmental data
TOS.env <- TOS.noNA[,3:dim(TOS.noNA)[2]]
head(TOS.env)

# Moran’s coefficient of spatial autocorrelation
#sp.env <- SpatialPolygonsDataFrame(Sr = aus.sp, data = TOS.env, match.ID = FALSE)
MC.env <- moran.randtest(TOS.env, listwgab, nrepet = 500)
MC.env
mc.bounds <- moran.bounds(listwgab)
mc.bounds
env.maps <- s1d.barchart(MC.env$obs, labels = MC.env$names, plot = FALSE, xlim = 1.1 * mc.bounds, paxes.draw = TRUE, 
                         pgrid.draw = FALSE)
addline(env.maps, v = mc.bounds, plot = TRUE, pline.col = 'red', pline.lty = 3)

colnames(TOS.env)

NP.lgs <- moranNP.randtest(TOS.env[,1], listwgab, nrepet = 999, alter = "two-sided") 
NP.lgs
plot(NP.lgs)
sum(NP.lgs$obs)
MC.env$obs[1]

NP.clay <- moranNP.randtest(TOS.env[,2], listwgab, nrepet = 999, alter = "two-sided") 
NP.clay # positive autocorrelation only
plot(NP.clay)
sum(NP.clay$obs)
MC.env$obs[2]

NP.pH15 <- moranNP.randtest(TOS.env[,3], listwgab, nrepet = 999, alter = "two-sided")
NP.pH15
plot(NP.pH15)
sum(NP.pH15$obs)
MC.env$obs[3]

NP.elecond <- moranNP.randtest(TOS.env[,4], listwgab, nrepet = 999, alter = "two-sided")
NP.elecond
plot(NP.elecond)
sum(NP.elecond$obs)
MC.env$obs[4]

NP.lLOI <- moranNP.randtest(TOS.env[,5], listwgab, nrepet = 999, alter = "two-sided")
NP.lLOI
plot(NP.lLOI)
sum(NP.lLOI$obs)
MC.env$obs[5]

NP.prescott <- moranNP.randtest(TOS.env[,6], listwgab, nrepet = 999, alter = "two-sided")
NP.prescott
plot(NP.prescott)
sum(NP.prescott$obs)
MC.env$obs[6]

NP.soilH20 <- moranNP.randtest(TOS.env[,7], listwgab, nrepet = 999, alter = "two-sided")
NP.soilH20
plot(NP.soilH20)
sum(NP.soilH20$obs)
MC.env$obs[7]

NP.llai <- moranNP.randtest(TOS.env[,8], listwgab, nrepet = 999, alter = "two-sided")
NP.llai
plot(NP.llai)
sum(NP.llai$obs)
MC.env$obs[8]

NP.soillN <- moranNP.randtest(TOS.env[,9], listwgab, nrepet = 999, alter = "two-sided")
NP.soillN
plot(NP.soillN)
sum(NP.soillN$obs)
MC.env$obs[9]

NP.soillP <- moranNP.randtest(TOS.env[,10], listwgab, nrepet = 999, alter = "two-sided")
NP.soillP
plot(NP.soillP)
sum(NP.soillP$obs)
MC.env$obs[10]

NP.lAl <- moranNP.randtest(TOS.env[,11], listwgab, nrepet = 999, alter = "two-sided")
NP.lAl
plot(NP.lAl)
sum(NP.lAl$obs)
MC.env$obs[11]

NP.lMn <- moranNP.randtest(TOS.env[,12], listwgab, nrepet = 999, alter = "two-sided")
NP.lMn
plot(NP.lMn)
sum(NP.lMn$obs)
MC.env$obs[12]

NP.lCu <- moranNP.randtest(TOS.env[,13], listwgab, nrepet = 999, alter = "two-sided")
NP.lCu
plot(NP.lCu)
sum(NP.lCu$obs)
MC.env$obs[13]

NP.lPb <- moranNP.randtest(TOS.env[,14], listwgab, nrepet = 999, alter = "two-sided")
NP.lPb
plot(NP.lPb)
sum(NP.lPb$obs)
MC.env$obs[14]

NP.lSb <- moranNP.randtest(TOS.env[,15], listwgab, nrepet = 999, alter = "two-sided")
NP.lSb
plot(NP.lSb)
sum(NP.lSb$obs)
MC.env$obs[15]

NP.lNi <- moranNP.randtest(TOS.env[,16], listwgab, nrepet = 999, alter = "two-sided")
NP.lNi
plot(NP.lNi)
sum(NP.lNi$obs)
MC.env$obs[16]

NP.lKThU <- moranNP.randtest(TOS.env[,17], listwgab, nrepet = 999, alter = "two-sided")
NP.lKThU
plot(NP.lKThU)
sum(NP.lKThU$obs)
MC.env$obs[17]





########################
## distribution model ##
## based on train rf  ##
## (superceded)       ##
########################

# store boundaries in a single extent object
geogr.extent <- ext(soilclay)

# prepare data
TOS.dm.dat <- na.omit(TOS.test[,c("LON","LAT","lHg","lgs","lgtclay","pH15","lelecond","lLOI","prescott","soilH20",
                                "llai","soillN","soillP","lAl","lMn","lCu","lPb","lSb","lNi","lKThU","ferric2",
                                "Alox","Feox","Pox")])

HgDat <- data.frame(
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
saveRDS(HgDat, file = "HgDat.rds")

# convert to spatial points
HgSpat <- HgDat
coordinates(HgSpat) <- ~x+y
saveRDS(HgSpat, file = "HgSpat.rds")


# prepare training data
training_data <- data.frame(
  Hg = HgSpat$Hg,
  # environmental predictors
  clay = HgSpat$clay,
  pH = HgSpat$pH,
  Prescott = HgSpat$Prescott,
  soilH20 = HgSpat$soilH20,
  lai = HgSpat$lai,
  soilN = HgSpat$soilN,
  soilP = HgSpat$soilP,
  KThU = HgSpat$KThU,
  ferric2=HgSpat$ferric2,
  Alox=HgSpat$Alox,
  Feox=HgSpat$Feox,
  Pox=HgSpat$Pox
)

crs(HgSpat) <- crsUse

# choose size of spatial blocks for spatial autocorrelation
blocksize <- cv_spatial_autocor(x=HgSpat, column='Hg', plot=T, deg_to_metre=111139) 
blocksize$range

cv_block_size(x = HgSpat, column = 'Hg') # interactive block-size visalisation

# spatial blocks for cross-validation
spatial_blocks <- cv_spatial(x=HgSpat, k=5, hexagon=F, size=200000)
spatial_blocks$folds_list

# custom indices for spatial cross-validation
spatial_indices <- spatial_blocks$folds_list

# define spatial cross-validation control
sactrl <- trainControl(
  method = "cv",
  classProbs = F,
  number = 5,
  index = spatial_indices,
  savePredictions = TRUE
)

# train random forest model (no control for spatial autocorrelation)
rf_model <- train(
  Hg ~ .,
  data = training_data,
  method = "rf",
  trControl = trainControl(method = "cv", number = 10)
)
rf_model$results
rf_model$finalModel$importance
str(rf_model)


# create prediction raster stack of environmental variables
envVars <- c(clay.rsmp,pH.rsmp,Prescott.rsmp,soilH20.rsmp,lai.rsmp,soilN.rsmp,soilP.rsmp,KThU.rsmp,
             ferric2.rsmp,Alox.rsmp,Feox.rsmp,Pox.rsmp)
names(envVars) <- c("clay","pH","Prescott","soilH20","lai","soilN","soilP","KThU","ferric2","Alox","Feox","Pox")
names(envVars)
plot(envVars)
envVars
saveRDS(envVars, file = "envVars.rds")

# spatial predictions
HgPred <- terra::predict(envVars, rf_model, na.rm=T)
plot(HgPred)
points(HgDat, pch=20, col="black", cex=0.5)
HgPred.rf.bt <- HgPred^10
plot(HgPred.rf.bt)
writeRaster(HgPred.rf.bt, "HgPredRFbt.tif", overwrite=T)

# model cross-validation
cv <- train(
  Hg ~ .,
  data = training_data,
  method = "rf",
  trControl = trainControl(
    method = "cv",
    number = 10,
    savePredictions = TRUE
  )
)

# performance metrics
rmse <- sqrt(mean((cv$pred$obs - cv$pred$pred)^2))
rmse
r2 <- cor(cv$pred$obs, cv$pred$pred)^2
r2

# variable importance
varImp(rf_model)

# plot variable importance
plot(varImp(rf_model))





###############
## spatialRF ##
###############

# create distance matrix
# Haversine distance matrix function
# most accurate for geographic coordinates because it accounts for Earth's spherical shape
coords <- data.frame(
  longitude = HgDat$x,
  latitude = HgDat$y
)

haversine_matrix <- distm(
  x = coords,
  fun = distHaversine
)

# format the distance matrix for spatialRF
#' spatialRF expects: square matrix; row and column names; symmetric distances
# add row and column names
rownames(haversine_matrix) <- paste0("point_", 1:dim(HgDat)[1])
colnames(haversine_matrix) <- paste0("point_", 1:dim(HgDat)[1])

# verify the distance matrix properties
#  checks ensure  matrix is properly formatted for spatialRF
stopifnot(
  nrow(haversine_matrix) == ncol(haversine_matrix),  # square matrix
  isSymmetric(haversine_matrix),                          # symmetric
  all(diag(haversine_matrix) == 0)                        # zero diagonal
)


#names of the response variable and the predictors
dependent.variable.name <- "Hg"
predictor.variable.names <- colnames(HgDat)[5:dim(HgDat)[2]]

#coordinates of the cases
xy <- HgDat[, c("x", "y")]

#distance matrix
distance.matrix <- haversine_matrix
range(distance.matrix)

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
    data = HgDat,
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
  data = HgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  ncol = 3,
  point.color = rev(viridis::viridis(100, option = "F")),
  line.color = "gray30"
)

# Moran's I
moransIplot <- spatialRF::plot_training_df_moran(
  data = HgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  fill.color = viridis::viridis(
    100,
    option = "F",
    direction = -1
  ),
  point.color = "gray40"
)
moransIplot
str(moransIplot)

moransIplot$data$distance.threshold
moransIplot$data$p.value
write.csv(moransIplot$data, "moransIplot.csv", row.names=F)


# possible interactions
interactions <- spatialRF::the_feature_engineer(
  data = HgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  xy = xy,
  importance.threshold = 0.50, # uses 50% best predictors
  cor.threshold = 0.60, # max corr between interactions and predictors
  repetitions = 100,
  verbose = TRUE
)

# **run only if promising interactions found**
kableExtra::kbl(
  head(interactions$screening, 10),
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
model.non.spatial <- spatialRF::rf(
  data = HgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, # not needed by rf, but other functions read it from the model
  verbose = FALSE,
  scaled.importance=T
)

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# residual diagnostics
spatialRF::plot_residuals_diagnostics(
  model.non.spatial,
  verbose = FALSE
)

# output multi-scale Moran's I of non-spatial model residuals
ns.resids <- spatialRF::plot_residuals_diagnostics(
  model.non.spatial,
  verbose = FALSE
)
str(ns.resids)
ns.MoranI.resid.dat <- ns.resids[[3]]$data
write.csv(ns.MoranI.resid.dat, "nsMoransIresid.csv", row.names=F)


# variable importance
spatialRF::plot_importance(
  model.non.spatial,
  verbose = FALSE
)

nsm.imp <- spatialRF::plot_importance(
  model.non.spatial,
  verbose = FALSE
)
nsm.imp

importance.df <- randomForestExplainer::measure_importance(
  model.non.spatial,
  measures = c("mean_min_depth", "no_of_nodes", "times_a_root", "p_value")
)

kableExtra::kbl(
  importance.df %>% 
    dplyr::arrange(mean_min_depth) %>% 
    dplyr::mutate(p_value = round(p_value, 4)),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# contribution of predictors to model transferability
# (predictive ability on independent spatial folds measured with rf_evaluate())
model.non.spatial <- spatialRF::rf_importance(
  model = model.non.spatial
)

model.non.spatial$importance$per.variable %>% 
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
local.importance <- spatialRF::get_importance_local(model.non.spatial)
kableExtra::kbl(
  round(local.importance[1:10, 1:5], 3),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# map changing variable importance of space
# adding coordinates
local.importance <- cbind(
  xy,
  local.importance
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
    data = local.importance,
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
    data = local.importance,
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
    data = local.importance,
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
    data = local.importance,
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
    data = local.importance,
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
    data = local.importance,
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

# no real evidence for changing variable importance over space

# response curves and surfaces
spatialRF::plot_response_curves(
  model.non.spatial,
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
  model.non.spatial,
  quantiles = 0.5,
  ncol = 4
)

# ferric2
pdp::partial(
  model.non.spatial, 
  train = HgDat, 
  pred.var = "ferric2", 
  plot = TRUE, 
  grid.resolution = 200
)

# soil N
pdp::partial(
  model.non.spatial, 
  train = HgDat, 
  pred.var = "soilN", 
  plot = TRUE, 
  grid.resolution = 200
)

# Prescott index
pdp::partial(
  model.non.spatial, 
  train = HgDat, 
  pred.var = "Prescott", 
  plot = TRUE, 
  grid.resolution = 200
)

# Pox
pdp::partial(
  model.non.spatial, 
  train = HgDat, 
  pred.var = "Pox", 
  plot = TRUE, 
  grid.resolution = 200
)

# lai
pdp::partial(
  model.non.spatial, 
  train = HgDat, 
  pred.var = "lai", 
  plot = TRUE, 
  grid.resolution = 200
)

# soil H2O
pdp::partial(
  model.non.spatial, 
  train = HgDat, 
  pred.var = "soilH20", 
  plot = TRUE, 
  grid.resolution = 200
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

spatialRF::print_performance(model.non.spatial)

# spatial cross-validation
model.non.spatial <- spatialRF::rf_evaluate(
  model = model.non.spatial,
  xy = xy,                  # data coordinates
  repetitions = 30,         # number of spatial folds
  training.fraction = 0.75, # training data fraction on each fold
  metrics = "r.squared",
  verbose = FALSE
)
spatialRF::plot_evaluation(model.non.spatial)

# predicted values of [Hg]
predicted.non.spatial <- stats::predict(object = model.non.spatial, data = HgDat, type = "response")$predictions

# spatial random forest
# Moran's I
spatialRF::plot_moran(
  model.non.spatial, 
  verbose = FALSE
)

 model.spatial <- spatialRF::rf_spatial(
   model = model.non.spatial,
   method = "mem.moran.sequential", #default method
   verbose = FALSE
 )

spatialRF::print_performance(model.spatial)

spatialRF::plot_moran(
  model.spatial, 
  verbose = FALSE
)

model.spatial$importance$per.variable$variable
names(envVars)
#envVars2 <- c(soilN.rsmp,lai.rsmp,Prescott.rsmp,soilH20.rsmp,KThU.rsmp,soilP.rsmp,pH.rsmp,clay.rsmp,
#              ferric2.rsmp,Alox.rsmp,Feox.rsmp,Pox.rsmp)
#names(envVars2) <- c("soilN", "lai", "Prescott", "soilH20", "KThU", "soilP", "pH", "clay", "ferric2",
#                     "Alox","Feox","Pox")
names(envVars2)
plot(envVars2)
envVars2
saveRDS(envVars2, file = "envVars2.rds")


HgPredNoSpat <- terra::predict(envVars2, model.non.spatial, na.rm=T)
plot(HgPredNoSpat)
points(HgDat, pch=20, col="black", cex=0.5)
HgPredNoSpat.rf.bt <- HgPredNoSpat^10
plot(HgPredNoSpat.rf.bt)
writeRaster(HgPredNoSpat.rf.bt, "HgPredNoSpatRFbt.tif", overwrite=T)

HgPredSpat <- terra::predict(envVars2, model.spatial, na.rm=T)
plot(HgPredSpat)
writeRaster(HgPredSpat, "HgPredSpatRFlog10.tif", overwrite=T)

HgPredSpatImp <- rast('HgPredSpatRFlog10.tif')
writeRaster(HgPredSpatImp, "HgPredSpatRFlog10.nc", filetype="netcdf", overwrite=TRUE) # write to NetCDF format

points(HgDat, pch=20, col="black", cex=0.5)
HgPredSpat.rf.bt <- HgPredSpat^10
plot(HgPredSpat.rf.bt)
writeRaster(HgPredSpat.rf.bt, "HgPredSpatRFbt.tif", overwrite=T)

HgPredSpatBTImp <- rast('HgPredSpatRFbt.tif')
writeRaster(HgPredSpatBTImp, "HgPredSpatRF.nc", filetype="netcdf", overwrite=TRUE) # write to NetCDF format

p1 <- spatialRF::plot_importance(
  model.non.spatial, 
  verbose = FALSE) + 
  ggplot2::ggtitle("non-spatial model") 

p2 <- spatialRF::plot_importance(
  model.spatial,
  verbose = FALSE) + 
  ggplot2::ggtitle("spatial model")

p1 | p2 

kableExtra::kbl(
  head(model.spatial$importance$per.variable, n = 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# model tuning
model.spatial <- spatialRF::rf_spatial(
  model = model.spatial,
  method = "mem.moran.sequential", #default method
  ranger.arguments = list(
    mtry = 5,
    min.node.size = 20,
    num.trees = 500
  ),
  verbose = FALSE,
  seed = random.seed
)

model.spatial <- rf_tuning(
  model = model.spatial,
  xy = xy,
  repetitions = 30,
  num.trees = c(500, 1000),
  mtry = seq(
    2,
    length(model.spatial$ranger.arguments$predictor.variable.names), # number of predictors
    by = 9),
  min.node.size = c(5, 15),
  verbose = FALSE
)

# increase memory to max
mem.maxVSize(v = Inf)

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  #1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# stochastic version
model.spatial.repeat <- spatialRF::rf_repeat(
  model = model.spatial, 
  repetitions = 250,
  verbose = FALSE
)

# stopping the cluster
parallel::stopCluster(cl = local.cluster)


# variable importance
stochSpatialRFvarImp <- spatialRF::plot_importance(
  model.spatial.repeat, 
  verbose = FALSE
)
stochSpatialRFvarImp
str(stochSpatialRFvarImp)
stochSpatialRFvarImp.dat <- stochSpatialRFvarImp$data

viXvar.stats <- stochSpatialRFvarImp.dat %>%
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
viXvar.stats
write.csv(viXvar.stats, "viXvarstats.csv", row.names=F)

# response curves
spatialRF::plot_response_curves(model=model.spatial.repeat, quantiles = c(0.5, 0.975, 0.025), ncol = 3)


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

  model.tuned <- spatialRF::rf_spatial(
  data = HgDat,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds2,
  xy = xy,
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
p <- spatialRF::plot_optimization(model.tuned)

# importance
spatialRF::plot_importance(
  model.tuned, 
  verbose = FALSE
)

# variable importance
stochSpatTunedRFvarImp <- spatialRF::plot_importance(
  model.tuned, 
  verbose = FALSE
)
stochSpatTunedRFvarImp
str(stochSpatTunedRFvarImp)
stochSpatTunedRFvarImp.dat <- stochSpatTunedRFvarImp$data

viXvar.stats <- stochSpatTunedRFvarImp.dat %>%
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
viXvar.stats
write.csv(viXvar.stats, "viXvarstatsTunedSpatial.csv", row.names=F)

# performance
spatialRF::print_performance(model.tuned)
model.tuned$importance$per.variable$variable

# spatial predictors
spatial.predictors <- spatialRF::get_spatial_predictors(model.tuned)
class(spatial.predictors)
spatpred.names <- colnames(spatial.predictors)
Nspatpreds <- length(spatpred.names)
Nspatpreds

pr <- data.frame(spatial.predictors, HgDat[, c("x", "y")])
head(pr)

kableExtra::kbl(
  head(model.tuned$importance$per.variable, n = 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

p1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = aus.sf, fill = "white") +
  ggplot2::geom_point(
    data = pr,
    ggplot2::aes(
      x = x,
      y = y,
      color = spatial_predictor_100000_27	
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "eigenvalue") +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::ggtitle("spatial_predictor_100000_27	") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

p2 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = aus.sf, fill = "white") +
  ggplot2::geom_point(
    data = pr,
    ggplot2::aes(
      x = x,
      y = y,
      color = spatial_predictor_100000_15
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "eigenvalue") +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::ggtitle("spatial_predictor_100000_15") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

p3 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = aus.sf, fill = "white") +
  ggplot2::geom_point(
    data = pr,
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
  ggplot2::ggtitle("spatial_predictor_100000_18	") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

p4 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = aus.sf, fill = "white") +
  ggplot2::geom_point(
    data = pr,
    ggplot2::aes(
      x = x,
      y = y,
      color = spatial_predictor_100000_29
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(option = "F") +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "eigenvalue") +
  ggplot2::scale_x_continuous(limits = c(112, 155)) +
  ggplot2::scale_y_continuous(limits = c(-45, -10))  +
  ggplot2::ggtitle("spatial_predictor_100000_29") + 
  ggplot2::theme(legend.position = "bottom")+ 
  ggplot2::xlab("longitude") + 
  ggplot2::ylab("latitude")

(p1 + p2) / (p3 + p4)


# response curves
respCurvRF <- spatialRF::plot_response_curves(model=model.tuned, quantiles = c(0.5), ncol = 3)
respCurvRF
str(respCurvRF[[1]])

respCurvRF[[1]]$labels$x


# ferric2
respCurvRF[[1]]
respCurvRF[[1]]$labels$x
ferricrespCurvRF <- respCurvRF[[1]]
str(ferricrespCurvRF)

table(ferricrespCurvRF$layers[[1]]$data$id)

ferricpredRF.dat <- data.frame(iter=ferricrespCurvRF$layers[[1]]$data$id, 
                               ferric=ferricrespCurvRF$layers[[1]]$data$ferric,
                               Hg=ferricrespCurvRF$layers[[1]]$data$Hg,
                               quantile=ferricrespCurvRF$layers[[1]]$data$quantile)

head(ferricpredRF.dat)
ferricpredRFq.5.dat <- subset(ferricpredRF.dat, quantile==0.5)
dim(ferricpredRFq.5.dat)

ferricpredRF.5Xiter.stats <- ferricpredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(ferric, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
ferricpredRF.5Xiter.stats
range(ferricpredRF.5Xiter.stats$x, na.rm=T)
c(min(ferricpredRF.5Xiter.stats$lowery, na.rm=T), max(ferricpredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(ferricpredRF.5Xiter.stats, "ferricpredRF.5Xiter.stats.csv", row.names=F)


# soil N
respCurvRF[[2]]$labels$x
soilNrespCurvRF <- respCurvRF[[2]]
str(soilNrespCurvRF)

table(soilNrespCurvRF$layers[[1]]$data$id)

soilNpredRF.dat <- data.frame(iter=soilNrespCurvRF$layers[[1]]$data$id, 
           soilN=soilNrespCurvRF$layers[[1]]$data$soilN,
           Hg=soilNrespCurvRF$layers[[1]]$data$Hg,
           quantile=soilNrespCurvRF$layers[[1]]$data$quantile)
head(soilNpredRF.dat)
soilNpredRFq.5.dat <- subset(soilNpredRF.dat, quantile==0.5)
dim(soilNpredRFq.5.dat)

soilNpredRF.5Xiter.stats <- soilNpredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(soilN, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
soilNpredRF.5Xiter.stats
range(soilNpredRF.5Xiter.stats$x, na.rm=T)
c(min(soilNpredRF.5Xiter.stats$lowery, na.rm=T), max(soilNpredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(soilNpredRF.5Xiter.stats, "soilNpredRF.5Xiter.stats.csv", row.names=F)

# Prescott
respCurvRF[[3]]$labels$x
PrescottrespCurvRF <- respCurvRF[[3]]
str(PrescottrespCurvRF)

table(PrescottrespCurvRF$layers[[1]]$data$id)

PrescottpredRF.dat <- data.frame(iter=PrescottrespCurvRF$layers[[1]]$data$id, 
                              Prescott=PrescottrespCurvRF$layers[[1]]$data$Prescott,
                              Hg=PrescottrespCurvRF$layers[[1]]$data$Hg,
                              quantile=PrescottrespCurvRF$layers[[1]]$data$quantile)
head(PrescottpredRF.dat)
PrescottpredRFq.5.dat <- subset(PrescottpredRF.dat, quantile==0.5)
dim(PrescottpredRFq.5.dat)

PrescottpredRF.5Xiter.stats <- PrescottpredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(Prescott, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
PrescottpredRF.5Xiter.stats
range(PrescottpredRF.5Xiter.stats$x, na.rm=T)
c(min(PrescottpredRF.5Xiter.stats$lowery, na.rm=T), max(PrescottpredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(PrescottpredRF.5Xiter.stats, "PrescottpredRF.5Xiter.stats.csv", row.names=F)

# LAI
respCurvRF[[4]]$labels$x
lairespCurvRF <- respCurvRF[[4]]
str(lairespCurvRF)

table(lairespCurvRF$layers[[1]]$data$id)

laipredRF.dat <- data.frame(iter=lairespCurvRF$layers[[1]]$data$id, 
                                 lai=lairespCurvRF$layers[[1]]$data$lai,
                                 Hg=lairespCurvRF$layers[[1]]$data$Hg,
                                 quantile=lairespCurvRF$layers[[1]]$data$quantile)
head(laipredRF.dat)
laipredRFq.5.dat <- subset(laipredRF.dat, quantile==0.5)
dim(laipredRFq.5.dat)

laipredRF.5Xiter.stats <- laipredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(lai, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
laipredRF.5Xiter.stats
range(laipredRF.5Xiter.stats$x, na.rm=T)
c(min(laipredRF.5Xiter.stats$lowery, na.rm=T), max(laipredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(laipredRF.5Xiter.stats, "laipredRF.5Xiter.stats.csv", row.names=F)

# Pox
respCurvRF[[5]]$labels$x
PoxrespCurvRF <- respCurvRF[[5]]
str(PoxrespCurvRF)

table(PoxrespCurvRF$layers[[1]]$data$id)

PoxpredRF.dat <- data.frame(iter=PoxrespCurvRF$layers[[1]]$data$id,
                            Pox=PoxrespCurvRF$layers[[1]]$data$Pox,
                            Hg=PoxrespCurvRF$layers[[1]]$data$Hg,
                            quantile=PoxrespCurvRF$layers[[1]]$data$quantile)
head(PoxpredRF.dat)
PoxpredRFq.5.dat <- subset(PoxpredRF.dat, quantile==0.5)
dim(PoxpredRFq.5.dat)

PoxpredRF.5Xiter.stats <- PoxpredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(Pox, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
PoxpredRF.5Xiter.stats
range(PoxpredRF.5Xiter.stats$x, na.rm=T)
c(min(PoxpredRF.5Xiter.stats$lowery, na.rm=T), max(PoxpredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(PoxpredRF.5Xiter.stats, "PoxpredRF.5Xiter.stats.csv", row.names=F)


# soil P
soilPrespCurvRF <- spatialRF::plot_response_curves(model=model.tuned, variables="soilP", quantiles = c(0.5), ncol = 3)
soilPrespCurvRF$labels$x
str(soilPrespCurvRF)

table(soilPrespCurvRF$layers[[1]]$data$id)

soilPpredRF.dat <- data.frame(iter=soilPrespCurvRF$layers[[1]]$data$id,
                            soilP=soilPrespCurvRF$layers[[1]]$data$soilP,
                            Hg=soilPrespCurvRF$layers[[1]]$data$Hg,
                            quantile=soilPrespCurvRF$layers[[1]]$data$quantile)
head(soilPpredRF.dat)
soilPpredRFq.5.dat <- subset(soilPpredRF.dat, quantile==0.5)
dim(soilPpredRFq.5.dat)

soilPpredRF.5Xiter.stats <- soilPpredRFq.5.dat %>%
  group_by(iter) %>%
  summarise(
    x = mean(soilP, na.rm = TRUE),
    meany = mean(Hg, na.rm = TRUE),
    uppery = quantile(Hg, probs=0.975, na.rm = TRUE),
    lowery = quantile(Hg, probs=0.025, na.rm = TRUE),
  )
soilPpredRF.5Xiter.stats
range(soilPpredRF.5Xiter.stats$x, na.rm=T)
c(min(soilPpredRF.5Xiter.stats$lowery, na.rm=T), max(soilPpredRF.5Xiter.stats$uppery, na.rm=T))
write.csv(soilPpredRF.5Xiter.stats, "soilPpredRF.5Xiter.stats.csv", row.names=F)



# comparing spatial and non-spatial models stochastically
comparison <- spatialRF::rf_compare(
  models = list(
    'non-spatial' = model.non.spatial,
    'spatial' = model.tuned
  ),
  xy = xy,
  repetitions = 100,
  training.fraction = 0.8,
  metrics = "r.squared"
)

x <- comparison$comparison.df %>% 
  dplyr::group_by(model, metric) %>% 
  dplyr::summarise(value = round(median(value), 3)) %>% 
  dplyr::arrange(metric) %>% 
  as.data.frame()
colnames(x) <- c("model", "metric", "median")
kableExtra::kbl(
  x,
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

