# Extract HadISST data for sediment-trap comparison
# Lukas Jonkers, last edit 23 Oct 2018
library(raster)

# download HadSST data from https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html

# prepare data ####
# load HadISST SST data
HadSST <- brick('HadISST_sst.nc')

# replace -1000 by NA
na.fun <- function(x) {replace(x, x == -1000, NA)}
SST <- calc(HadSST, na.fun)
rm(HadSST)

# get dates in HadSST
Had.date <- substr(names(SST),2,11)
Had.date <- as.Date(gsub('\\.', '\\-', Had.date))
Had.month <- format(Had.date, '%Y-%m')
Had.year <- format(Had.date, '%Y')

# mean SST over specified years
start.year <- 1870
end.year <- 1899
sub.indx <- which(Had.year >= start.year & Had.year <= end.year)
MSST <- mean(subset(SST, sub.indx), na.rm = TRUE)

# function to extract SST values using polygons (circle) around point
# polygons made using make_polygons.R
extractSST.fun <- function(dat, circle){
  tem <- extract(dat, circle, small = TRUE, weights= TRUE, df = TRUE)
  tem <- subset(tem, select = -ID)
  weights <- tem$weight
  vals <- subset(tem, select = -weight)
  if(nrow(tem)>1){
    apply(vals, 2, function(x) sum(x[!is.na(x)]*weights[!is.na(x)])/sum(weights[!is.na(x)]))}
  else{
    as.numeric(vals)
  }
} 

# historical SST forcens ####
ForCenS_polys <- readRDS('ForCens_polys.RDS')
ForCenS.SST <- sapply(ForCenS_polys, function(x) extractSST.fun(MSST, x))
saveRDS(ForCenS.SST, 'forcens_HadSST_1870-1899.RDS')

# historical SST sediment traps ####
# load polygons
trap_polys <- readRDS('trap_polys.RDS')
trap.SST <- sapply(trap_polys, function(x) extractSST.fun(MSST, x))
names(trap.SST) <- names(trap_polys)
saveRDS(trap.SST, 'traps_HadSST_1870-1899.RDS')

# get mean temperature over collecting period of sediment trap ####
# needs annual assemblages (TS.assem) for length of the time series
get.mean.temp <- function(trap){
  begin <- format(as.Date(dat.sel[[which(names(dat.sel) == trap)]]$dat$open[1]), '%Y-%m')
  n.years <- nrow(TS.assem[[which(names(TS.assem) == trap)]])
  begin.indx <- which(Had.month == begin)
  end.indx <- (begin.indx + n.years*12)-1
  trap_poly <- trap_polys[[which(names(trap_polys) == trap)]]
  temperature <- extractSST.fun(SST[[begin.indx:end.indx]], trap_poly)
  mean.TS.SST <- mean(temperature, na.rm = TRUE)
  an.mean.TS.SST <- colMeans(matrix(temperature, nrow = 12), na.rm = TRUE)
  list(mean.TS.SST = mean.TS.SST, an.mean.TS.SST = an.mean.TS.SST)
}

dat.sel <- readRDS('dat_sel.RDS')
TS.assem <- readRDS('TS.assem.rds')

TS.real.SST <- lapply(names(TS.assem), get.mean.temp)
names(TS.real.SST) <- names(TS.assem)
saveRDS(TS.real.SST, 'traps_HadSST_period.RDS')


### get linear temperature trend over time frame specified in subset for Figure 1 ####
SST_sub <- subset(SST, which(Had.year >= 1870 & Had.year <= 2015))
# function to extract linear trend
lm_fun = function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[2] }}
# get time, each step is 1 month
time <- 1:nlayers(SST_sub)
# get decadal trend
SST.slope <- calc(SST_sub, lm_fun) * 120
# convert to tibble using ggplot_data function
gplot_data <- function(x, maxpixels = 100000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x)))
  names(dat) <- c('value', 'variable')
  dat <- dplyr::as.tbl(data.frame(coords, dat))
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]],
                            by = c("value" = "ID"))
  }
  dat
}
SST.slope.tibble <- gplot_data(SST.slope)
# save
saveRDS(SST.slope.tibble, file = 'hadisst_trend_1870-2015.RDS')

### get mean temperature  over time frame specified in subset to use as background for Figure 3 ####
# time frame here spans all years in which sediment traps were active
SST_sub <- subset(SST, which(Had.year >= 1978 & Had.year <= 2013))
SST.mean <- mean(SST_sub, na.rm = TRUE)
# convert to tibble using ggplot_data function
gplot_data <- function(x, maxpixels = 100000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x)))
  names(dat) <- c('value', 'variable')
  
  dat <- dplyr::as.tbl(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]],
                            by = c("value" = "ID"))
  }
  dat
}
SST.mean.tibble <- gplot_data(SST.mean)
# save
saveRDS(SST.mean.tibble, file = 'hadisst_mean_1978-2013.RDS')
