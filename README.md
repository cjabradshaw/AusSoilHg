# Predicting continental distribution of soil mercury concentration in Australia
<a href="https://doi.org/10.5281/zenodo.14824785"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14824785.svg" alt="DOI"></a>
<img align="right" src="www/HgPredSpatRFbt.jpg" alt="predicted [Hg]" width="300" style="margin-top: 20px">

Professor <a href="https://globalecologyflinders.com/people">Corey J. A. Bradshaw</a><br>
<a href="https://globalecologyflinders.com/">Global Ecology | <em>Partuyarta Ngadluku Wardli Kuu</em></a>, <a href="https://flinders.edu.au">Flinders University</a><br>
<a href="mailto:corey.bradshaw@flinders.edu.au">e-mail</a><br>
<br>
<strong>Team</strong>:<br>
- Associate Professor <a href="https://researchportalplus.anu.edu.au/en/persons/larissa-schneider">Larissa Schneider</a>, Australian National University
- Adjunct Professor <a href="https://scholar.google.com.au/citations?user=O3mHBygAAAAJ&hl=en">Patrice de Caritat</a>, Curtin University
- Professor <a href="https://researchportalplus.anu.edu.au/en/persons/simon-haberle">Simon Haberle</a>, Australian National University

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
- <em><a href="data/samplept/Hg.csv">Hg.csv</a></em>: re-analysed [Hg] estimates (ng/g)
- <em><a href="data/samplept/gs.csv">gs.csv</a></em>: grain-size category percentages
 
### <a href="data/spatial">Spatial</a>
Most of these files are too large to store in this repository directly, so in most cases the links refer to the original repository URLs where you can download the datasets.<br>
- <em><a href="data/spatial/aus.zip">aus.shp</a></em>: Australia boundary shapefile (zipped)
- <em><a href="https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip">wwf_terr_ecos.shp</a></em>: WWF ecoregions shapefile (download zipped file from original site)
- <em><a href="https://d28rz98at9flks.cloudfront.net/74619/74619_1M_shapefiles.zip">GeologicUnitPolygons1M.shp</a></em>: 1:1,000,000 geological unit polygon shapefile (download zipped file from original site)
- <em><a href="data/spatial/lithreclass.csv">lithreclass.csv</a></em>: reclassified lithology groups text file
- <em><a href="https://www.agriculture.gov.au/abares/aclump/land-use/data-download">NLUM_v7_250_ALUMV8_2020_21_alb.tif</a></em>: land use of Australia 2010–11 to 2020–21 (download geotif raster from original site)
- <em><a href="https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/144413">radmap_v4_2019_filtered_ML_KThU_RGB_24bit.tif</a></em>: potassium:thorium:uranium geotif raster (download data package from original site)
- <em><a href="https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/144413">radmap_v4_2019_filtered_ML_ppmTh_32bitfloat_grid.tif</a></em>: thorium ppm geotif raster (download data package from original site)
- <em><a href="https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/144413">radmap_v4_2019_filtered_ML_ppmU_32bitfloat_grid.tif</a></em>: uranium ppm geotif raster (download data package from original site)
- <em><a href="https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/144413">radmap_v4_2019_filtered_ML_pctk_32bitfloat_grid.tif</a></em>: % potassium geotif raster (download data package from original site)
- <em><a href="https://data.csiro.au/collection/csiro:61522?_st=browse&_str=2&_si=2&browseType=kw&browseValue=total%20soil%20nitrogen">NTO_000_005_EV_N_P_AU_NAT_C_20231101.tif</a></em>: soil nitrogen (0-5 cm) geotif raster (download from original site)
- <em><a href="https://data.csiro.au/collection/csiro:61526?_st=browse&_str=2&_si=1&browseType=kw&browseValue=total%20soil%20nitrogen">PTO_000_005_EV_N_P_AU_NAT_C_20231101.tif</a></em>: soil phosphorus (0-5 cm) geotif raster (download from original site)
- <em><a href="https://data.csiro.au/collection/csiro%3A11030v4">pHc_000_005_EV_N_P_AU_NAT_C_20140801.tif</a></em>: soil pH (0-5 cm) geotif raster (download from original site)
- <em><a href="http://dap.nci.org.au/thredds/remoteCatalogService'command=subset&catalog=http://dapds00.nci.org.au/thredds/catalog/ub8/au/OzWALD/annual//catalog.xml&dataset=ub8-au/OzWALD/annual/OzWALD.annual.Pg.AnnualSums.nc">OzWALD.annual.Pg.AnnualSums.nc</a></em>: annual precipitation <a href="https://www.unidata.ucar.edu/software/netcdf/">NetCDF</a> (download from original site)
- <em><a href="https://data.csiro.au/collection/csiro:9636v2">PrescottIndex_01_3s_lzw.tif</a></em>: Prescott index geotif raster (download from original site)
- <em><a href="http://dap.nci.org.au/thredds/remoteCatalogService'command=subset&catalog=http://dapds00.nci.org.au/thredds/catalog/ub8/au/OzWALD/annual//catalog.xml&dataset=ub8-au/OzWALD/annual/OzWALD.Ssoil.AnnualMeans.nc">OzWALD.Ssoil.AnnualMeans.nc</a></em>: soil water availability NetCDF (download from original site)
- <em><a href="http://dap.nci.org.au/thredds/remoteCatalogService'command=subset&catalog=http://dapds00.nci.org.au/thredds/catalog/ub8/au/OzWALD/annual//catalog.xml&dataset=ub8-au/OzWALD/annual/OzWALD.LAI.AnnualMeans.nc">OzWALD.LAI.AnnualMeans.nc</a></em>: leaf area index NetCDF (download from original site)
- <em><a href="http://dap.nci.org.au/thredds/remoteCatalogService'command=subset&catalog=http://dapds00.nci.org.au/thredds/catalog/ub8/au/OzWALD/annual//catalog.xml&dataset=ub8-au/OzWALD/annual/OzWALD.GPP.AnnualSums.nc">OzWALD.GPP.AnnualMeans.nc</a></em>: vegetation carbon uptake (gross primary production) NetCDF (download from original site)
- <em><a href="https://data.csiro.au/collection/csiro:55684">CLY_000_005_EV_N_P_AU_TRN_N_20210902.tif</a></em>: % soil clay content geotif raster (download from original site)
- <em><a href="https://data.csiro.au/collection/csiro:10688?q=soil%20silt&_st=keyword&_str=12&_si=1">SLT_000_005_EV_N_P_AU_TRN_N_20210902.tif</a></em>: % soil silt content geotif raster (download from original site)

### <a href="data/predicted">Predicted</a>
- <em>HgPredSpatRFlog10.nc</em>: This is the <a href="https://www.unidata.ucar.edu/software/netcdf/">NetCDF</a> file for the output map of Australia-wide Hg concentration predicted from the random forest model. Due to Github file-size constraints, we have broken the file into similar-sized chunks, and then compressed them using the fast, lossless compression algorithm <code><a href="https://man.archlinux.org/man/zstd.1.en">zstd</a></code>. First, decompress each .zst chunk using the following command in Terminal: <code>zstd -d 'HgPredSpatRFlog10.nc_chunk_a*.zst'</code>, and then combine chunks <em>a</em> to <em>i</em> using the following <a href="https://support.apple.com/en-au/guide/terminal/welcome/mac">Terminal</a> (or equivalent) command: <code>cat HgPredSpatRFlog10.nc_chunk_* > HgPredSpatRFlog10.nc</code>. You can import the NetCDF file in R using the <code>ncdf4</code> package and its function <a href="https://search.r-project.org/CRAN/refmans/ncdf4/html/nc_open.html"><em>nc_open</em></a>. This produces a <code><a href="https://cran.r-project.org/web/packages/ncdf4/index.html">ncdf4</a></code> object that can be converted to a <code><a href="https://rdrr.io/cran/terra/man/SpatRaster-class.html">SpatRaster</a></code> object using the <em><a href="https://rdrr.io/cran/terra/man/rast.html">rast</a></em> function in package <code><a href="https://rspatial.org/pkg/index.html">terra</a></code>.

## required R libraries
- <code><a href="https://adeverse.github.io/ade4/">ade4</a></code>,<code><a href="https://cran.r-project.org/web/packages/adegraphics/index.html">adegraphics</a></code>,<code><a href="https://cran.r-project.org/web/packages/adespatial/index.html">adespatial</a></code>,<code><a href="https://cran.r-project.org/web/packages/boot/index.html">boot</a></code>,<code><a href="https://cran.r-project.org/web/packages/dismo/index.html">dismo</a></code>,<code><a href="https://dplyr.tidyverse.org">dplyr</a></code>,<code><a href="https://cran.r-project.org/web/packages/fields/index.html">fields</a></code>,<code><a href="https://github.com/gbm-developers/gbm">gbm</a></code>,<code><a href="https://cran.r-project.org/web/packages/geodata/index.html">geodata</a></code>,<code><a href="https://cran.r-project.org/web/packages/geosphere/index.html">geosphere</a></code>,<code><a href="https://ggplot2.tidyverse.org">ggplot2</a></code>,<code><a href="https://github.com/haozhu233/kableExtra">kableExtra</a></code>,<code><a href="https://cran.r-project.org/web/packages/ncdf4/index.html">ncdf4</a></code>,<code><a href="https://patchwork.data-imaginist.com">patchwork</a></code>,<code><a href="https://github.com/bgreenwell/pdp">pdp</a></code>,<code><a href="https://cran.rstudio.com/web/packages/randomForestExplainer/vignettes/randomForestExplainer.html">randomForestExplainer</a></code>,<code><a href="https://docs.ropensci.org/rnaturalearthdata/">rnaturalearthdata</a></code>,<code><a href="https://r-spatial.github.io/sf/">sf</a></code>,<code>sp</code>,<code><a href="https://github.com/BlasBenito/spatialRF">spatialRF</a></code>,<code><a href="https://spatstat.org">spatstat</a></code>,<code><a href="https://r-spatial.github.io/spdep/">spdep</a></code>,<code><a href="https://rspatial.org/pkg/index.html">terra</a></code>,<code><a href="https://www.tidyverse.org">tidyverse</a></code>,<code><a href="https://cran.r-project.org/web/packages/usdm/index.html">usdm</a></code>

<p><a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Stacked_RGB_Master.jpg" alt="Flinders University logo" width="80" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="130" style="margin-top: 20px"></a>  &nbsp; &nbsp;
 <a href="https://ciehf.au"><img align="bottom-left" src="www/CIEHF_Logo_Email_Version Transparent.png" alt="CIEHF logo" width="200" style="margin-top: 20px"></a>  &nbsp; &nbsp; <a href="https://www.anu.edu.au"><img align="bottom-left" src="www/ANUlogo.png" alt="ANU logo" width="80" style="margin-top: 20px"></a>  &nbsp; &nbsp; <a href="https://www.curtin.edu.au"><img align="bottom-left" src="www/CUlogo.png" alt="ANU logo" width="50" style="margin-top: 20px"></a></p>
