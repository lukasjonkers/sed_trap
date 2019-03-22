# extract ERSST data for sediment-trap comparison
# using circular polygons created using get_polygons_function

library(raster)

# download ERSST v5 data from https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.ersst.v5.html
# To cite use of dataset: Boyin Huang, Peter W. Thorne, Viva F. Banzon, Tim Boyer, Gennady Chepurin, Jay H. Lawrimore, Matthew J. Menne, Thomas M. Smith, Russell S. Vose, and Huai-Min Zhang (2017): NOAA Extended Reconstructed Sea Surface Temperature (ERSST), Version 5. [indicate subset used]. NOAA National Centers for Environmental Information. doi:10.7289/V5T72FNM [access date].
# Please note: If you acquire NOAA_ERSST_V5 data products from PSD, we ask that you acknowledge us in your use of the data. This may be done by including text such as NOAA_ERSST_V5 data provided by the NOAA/OAR/ESRL PSD, Boulder, Colorado, USA, from their Web site at https://www.esrl.noaa.gov/psd/ in any documents or publications using these data. We would also appreciate receiving a copy of the relevant publications. This will help PSD to justify keeping the NOAA_ERSST_V5 data set freely available online in the future. Thank you!

ERSST <- rotate(brick('~/Seafile/SST_data/ERSSTv5/sst.mnmean.nc'))

# get dates in ERSST
ERSST.date <- getZ(ERSST)
ERSST.month <- format(ERSST.date, '%Y-%m')
ERSST.year <- format(ERSST.date, '%Y')

# mean annual SST over specified years
start.year <- 1854
end.year <- 1883
sub.indx <- which(ERSST.year >= start.year & ERSST.year <= end.year)
MSST <- mean(subset(ERSST, sub.indx), na.rm = TRUE)


# function to extract SST values using polygons
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
saveRDS(ForCenS.SST, 'forcens_ERSST_1854-1883.RDS')


# historical SST sediment traps ####
# load polygons
trap_polys <- readRDS('trap_polys.RDS')
trap.SST <- sapply(trap_polys, function(x) extractSST.fun(MSST, x))
names(trap.SST) <- names(trap_polys)
saveRDS(trap.SST, 'traps_ERSST_1854-1883.RDS')


# get mean temperature over collecting period of sediment trap ####
# needs annual assemblages (TS.assem) for length of the time series
get.mean.temp <- function(trap){
  begin <- format(as.Date(dat.sel[[which(names(dat.sel) == trap)]]$dat$open[1]), '%Y-%m')
  n.years <- nrow(TS.assem[[which(names(TS.assem) == trap)]])
  begin.indx <- which(ERSST.month == begin)
  end.indx <- (begin.indx + n.years*12)-1
  trap_poly <- trap_polys[[which(names(trap_polys) == trap)]]
  temperature <- extractSST.fun(ERSST[[begin.indx:end.indx]], trap_poly)
  mean.TS.SST <- mean(temperature, na.rm = TRUE)
  an.mean.TS.SST <- colMeans(matrix(temperature, nrow = 12), na.rm = TRUE)
  list(mean.TS.SST = mean.TS.SST, an.mean.TS.SST = an.mean.TS.SST)
}

dat.sel <- readRDS('dat_sel.RDS')
TS.assem <- readRDS('TS.assem.rds')

TS.real.SST <- lapply(names(TS.assem), get.mean.temp)
names(TS.real.SST) <- names(TS.assem)
saveRDS(TS.real.SST, 'traps_ERSST_period.RDS')
