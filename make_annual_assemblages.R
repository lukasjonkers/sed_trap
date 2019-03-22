# create annual assemblages from sediment trap time series
# Lukas Jonkers, last edit 20 Mar 2019

# load all shell flux time series
dat.sel <- readRDS('dat_sel.RDS')

# load function to make annual time series
source('make_annual_fluxes_function.R')

# imputation threshold, fraction of year
it <- 1/4
# get annual assemblages
TS.ann <- lapply(dat.sel, function(i) process.TS(i, it))
saveRDS(TS.ann, file = 'TS.ann.RDS')

# load ForCenS data split by region
species_domains <- readRDS('species_domains_compare.RDS')

# calculate assemblages from trap time series for each year
# keep only species that are in ForCenS, i.e. T. parkerae removed from Arabian Sea, G. ungulata from GOM
TS.assem <- lapply(TS.ann, function(x){
  domain <- x$meta$ocean
  spp <- names(species_domains[[which(names(species_domains) %in% domain)]]$species)
  flux <- lapply(x$flux, function(i) i[names(i) %in% spp])
  as.data.frame(t(sapply(flux, function(j) colSums(j)/sum(j))))
})
saveRDS(TS.assem, file = 'TS.assem.RDS')

# SSTs during collection interval can now be extracted using extract_HadSST