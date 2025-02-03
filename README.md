# Predicting continental distribution of soil mercury concentration in Australia
<img align="right" src="www/HgPredSpatRFbt.jpg" alt="predicted [Hg]" width="300" style="margin-top: 20px">

Professor <a href="https://globalecologyflinders.com/people">Corey J. A. Bradshaw</a><br>
Global Ecology, Flinders University<br>
<a href="mailto:corey.bradshaw@flinders.edu.au">e-mail</a><br>
<br>
<strong>Team</strong>:<br>
- Associate Professor <a href="https://researchportalplus.anu.edu.au/en/persons/larissa-schneider">Larissa Schneider</a>, Australian National University
- Adjunct Professor <a href="https://scholar.google.com.au/citations?user=O3mHBygAAAAJ&hl=en">Patrice de Caritat</a>, Curtin University

## Aims
- identify external and indirect determinants of mercury (Hg)
- understand environmental conditions that influence mercury retention and mobility
- predict continental distribution of soil mercury

## Scripts

## Data
### Sample point
- <em>geochem.csv</em>: geochemical data
- <em>field.csv</em>: sample point characteristics
- <em>hg.csv</em>: re-analysed [Hg] estimates
- <em>gs.csv</em>: grain-size category percentages
 
### Spatial
- <em><a href="data/spatial/aus.zip">aus.shp</a></em>: Australia boundary shapefile (zipped)
- <em><a href="https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip">wwf_terr_ecos.shp</a></em>: WWF ecoregions shapefile (zipped)
- <em>GeologicUnitPolygons1M.shp</em>: 1:5e6 geological unit polygon shapefile
- <em><a href="data/spatial/lithreclass.csv">lithreclass.csv</a></em>: reclassified lithology groups text file
- <em>radmap_v4_2019_filtered_ML_KThU_RGB_24bit.tif</em>: potassium:thorium:uranium geotif raster
- <em>radmap_v4_2019_filtered_ML_ppmTh_32bitfloat_grid.tif</em>: thorium ppm geotif raster
- <em>radmap_v4_2019_filtered_ML_ppmU_32bitfloat_grid.tif</em>: uranium ppm geotif raster
- <em>radmap_v4_2019_filtered_ML_pctk_32bitfloat_grid.tif</em>: % potassium geotif raster
- <em>NTO_000_005_EV_N_P_AU_NAT_C_20231101.tif</em>: soil nitrogen (0-5 cm) geotif raster
- <em>PTO_000_005_EV_N_P_AU_NAT_C_20231101.tif</em>: soil phosphorus (0-5 cm) geotif raster
- <em>pHc_000_005_EV_N_P_AU_NAT_C_20140801.tif</em>: soil pH (0-5 cm) geotif raster
- <em>OzWALD.annual.Pg.AnnualSums.nc</em>: annual precipitation NetCDF
- <em>PrescottIndex_01_3s_lzw.tif</em>: Prescott index geotif raster
- <em>OzWALD.Ssoil.AnnualMeans.nc</em>: soil water availability NetCDF
- <em>OzWALD.LAI.AnnualMeans.nc</em>: leaf area index NetCDF
- <em>OzWALD.GPP.AnnualMeans.nc</em>: vegetation carbon uptake (gross primary production) NetCDF
- <em>CLY_000_005_EV_N_P_AU_TRN_N_20210902.tif</em>: % soil clay content geotif raster
- <em>SLT_000_005_EV_N_P_AU_TRN_N_20210902.tif</em>: % soil silt content geotif raster

## required R libraries
- <code>ade4</code>,<code>adegraphics</code>,<code>adespatial</code>,<code>boot</code>,<code>dismo</code>,<code>dplyr</code>,<code>fields</code>,<code>gbm</code>,<code>geodata</code>,<code>geosphere</code>,<code>ggplot2</code>,<code>kableExtra</code>,<code>ncdf4</code>,<code>pdp</code>,<code>randomForestExplainer</code>,<code>rnaturalearthdata</code>,<code>sf</code>,<code>sp</code>,<code>spatialRF</code>,<code>spdep</code>,<code>terra</code>,<code>tidyverse</code>,<code>usdm</code>

<p><a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Stacked_RGB_Master.jpg" alt="Flinders University logo" width="80" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="130" style="margin-top: 20px"></a>  &nbsp; &nbsp;
 <a href="https://ciehf.au"><img align="bottom-left" src="www/CIEHF_Logo_Email_Version Transparent.png" alt="CIEHF logo" width="200" style="margin-top: 20px"></a>  &nbsp; &nbsp; <a href="https://www.anu.edu.au"><img align="bottom-left" src="www/ANUlogo.png" alt="ANU logo" width="80" style="margin-top: 20px"></a>  &nbsp; &nbsp; <a href="https://www.curtin.edu.au"><img align="bottom-left" src="www/CUlogo.png" alt="ANU logo" width="50" style="margin-top: 20px"></a></p>
