# Predicting continental distribution of soil mercury concentration in Australia
<img align="right" src="www/HgPredSpatRFbt.jpg" alt="predicted [Hg]" width="300" style="margin-top: 20px">

Professor <a href="https://globalecologyflinders.com/people">Corey J. A. Bradshaw</a><br>
<a href="https://globalecologyflinders.com/">Global Ecology | <em>Partuyarta Ngadluku Wardli Kuu</em></a>, <a href="https://flinders.edu.au">Flinders University</a><br>
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
- <code><a href="scripts/HgGH.R">HgGH.R</a></code>: all required R code combined

## <a href="data">Data</a>
### <a href="data/samplept">Sample point</a>
- <em><a href="data/samplept/geochem.csv">geochem.csv</a></em>: geochemical data
- <em><a href="data/samplept/field.csv">field.csv</a></em>: sample point characteristics
- <em><a href="data/samplept/Hg.csv">Hg.csv</a></em>: re-analysed [Hg] estimates
- <em><a href="data/samplept/gs.csv">gs.csv</a></em>: grain-size category percentages
 
### <a href="data/spatial">Spatial</a>
Most of these files are too large to store in this repository directly, so in most cases the links refer to the original repository URLs where you can download the datasets.<br>
- <em><a href="data/spatial/aus.zip">aus.shp</a></em>: Australia boundary shapefile (zipped)
- <em><a href="https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip">wwf_terr_ecos.shp</a></em>: WWF ecoregions shapefile (download zipped from original site)
- <em><a href="https://d28rz98at9flks.cloudfront.net/74619/74619_1M_shapefiles.zip">GeologicUnitPolygons1M.shp</a></em>: 1:1,000,000 geological unit polygon shapefile (download zipped from original site)
- <em><a href="data/spatial/lithreclass.csv">lithreclass.csv</a></em>: reclassified lithology groups text file
- <em><a href="https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/144413">radmap_v4_2019_filtered_ML_KThU_RGB_24bit.tif</a></em>: potassium:thorium:uranium geotif raster (download the data package from original site)
- <em><a href="https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/144413">radmap_v4_2019_filtered_ML_ppmTh_32bitfloat_grid.tif</a></em>: thorium ppm geotif raster (download the data package from original site)
- <em><a href="https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/144413">radmap_v4_2019_filtered_ML_ppmU_32bitfloat_grid.tif</a></em>: uranium ppm geotif raster (download the data package from original site)
- <em><a href="https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/144413">radmap_v4_2019_filtered_ML_pctk_32bitfloat_grid.tif</a></em>: % potassium geotif raster (download the data package from original site)
- <em><a href="https://data.csiro.au/collection/csiro:61522?_st=browse&_str=2&_si=2&browseType=kw&browseValue=total%20soil%20nitrogen">NTO_000_005_EV_N_P_AU_NAT_C_20231101.tif</a></em>: soil nitrogen (0-5 cm) geotif raster (download from original site)
- <em><a href="https://data.csiro.au/collection/csiro:61526?_st=browse&_str=2&_si=1&browseType=kw&browseValue=total%20soil%20nitrogen">PTO_000_005_EV_N_P_AU_NAT_C_20231101.tif</a></em>: soil phosphorus (0-5 cm) geotif raster (download from original site)
- <em><a href="https://data.csiro.au/collection/csiro%3A11030v4">pHc_000_005_EV_N_P_AU_NAT_C_20140801.tif</a></em>: soil pH (0-5 cm) geotif raster (download from original site)
- <em><a href="http://dap.nci.org.au/thredds/remoteCatalogService'command=subset&catalog=http://dapds00.nci.org.au/thredds/catalog/ub8/au/OzWALD/annual//catalog.xml&dataset=ub8-au/OzWALD/annual/OzWALD.annual.Pg.AnnualSums.nc">OzWALD.annual.Pg.AnnualSums.nc</a></em>: annual precipitation NetCDF (download from original site)
- <em><a href="https://data.csiro.au/collection/csiro:9636v2">PrescottIndex_01_3s_lzw.tif</a></em>: Prescott index geotif raster (download from original site)
- <em><a href="http://dap.nci.org.au/thredds/remoteCatalogService'command=subset&catalog=http://dapds00.nci.org.au/thredds/catalog/ub8/au/OzWALD/annual//catalog.xml&dataset=ub8-au/OzWALD/annual/OzWALD.Ssoil.AnnualMeans.nc">OzWALD.Ssoil.AnnualMeans.nc</a></em>: soil water availability NetCDF (download from original site)
- <em><a href="http://dap.nci.org.au/thredds/remoteCatalogService'command=subset&catalog=http://dapds00.nci.org.au/thredds/catalog/ub8/au/OzWALD/annual//catalog.xml&dataset=ub8-au/OzWALD/annual/OzWALD.LAI.AnnualMeans.nc">OzWALD.LAI.AnnualMeans.nc</a></em>: leaf area index NetCDF (download from original site)
- <em><a href="http://dap.nci.org.au/thredds/remoteCatalogService'command=subset&catalog=http://dapds00.nci.org.au/thredds/catalog/ub8/au/OzWALD/annual//catalog.xml&dataset=ub8-au/OzWALD/annual/OzWALD.GPP.AnnualSums.nc">OzWALD.GPP.AnnualMeans.nc</a></em>: vegetation carbon uptake (gross primary production) NetCDF (download from original site)
- <em><a href="https://data.csiro.au/collection/csiro:55684">CLY_000_005_EV_N_P_AU_TRN_N_20210902.tif</a></em>: % soil clay content geotif raster (download from original site)
- <em><a href="https://data.csiro.au/collection/csiro:10688?q=soil%20silt&_st=keyword&_str=12&_si=1">SLT_000_005_EV_N_P_AU_TRN_N_20210902.tif</a></em>: % soil silt content geotif raster (download from original site)

## required R libraries
- <code>ade4</code>,<code>adegraphics</code>,<code>adespatial</code>,<code>boot</code>,<code>dismo</code>,<code>dplyr</code>,<code>fields</code>,<code>gbm</code>,<code>geodata</code>,<code>geosphere</code>,<code>ggplot2</code>,<code>kableExtra</code>,<code>ncdf4</code>,<code>patchwork</code>,<code>pdp</code>,<code>randomForestExplainer</code>,<code>rnaturalearthdata</code>,<code>sf</code>,<code>sp</code>,<code>spatialRF</code>,<code>spdep</code>,<code>terra</code>,<code>tidyverse</code>,<code>usdm</code>

<p><a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Stacked_RGB_Master.jpg" alt="Flinders University logo" width="80" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="130" style="margin-top: 20px"></a>  &nbsp; &nbsp;
 <a href="https://ciehf.au"><img align="bottom-left" src="www/CIEHF_Logo_Email_Version Transparent.png" alt="CIEHF logo" width="200" style="margin-top: 20px"></a>  &nbsp; &nbsp; <a href="https://www.anu.edu.au"><img align="bottom-left" src="www/ANUlogo.png" alt="ANU logo" width="80" style="margin-top: 20px"></a>  &nbsp; &nbsp; <a href="https://www.curtin.edu.au"><img align="bottom-left" src="www/CUlogo.png" alt="ANU logo" width="50" style="margin-top: 20px"></a></p>
