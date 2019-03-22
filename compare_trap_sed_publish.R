# compare assemblages from sedimnet trap shell flux time series with ForCenS assemblages
# script by Lukas Jonkers, last edit: 22 Mar 2019

rm(list = ls())
library(rioja)
library(reshape2)
library(geosphere)
library(Hmisc)

# load all shell flux time series. Change to dat_small for sensitivity test with different shell sizes
dat.sel <- readRDS('dat_sel.RDS')
names(dat.sel) <- substr(names(dat.sel), 1, 3)

# load core top data from ForCenS, get lat and lon
ForCens <- readRDS('forcens_trimmed_compare.RDS')
ForCens.xy <- cbind.data.frame(x = ForCens$meta$Longitude, y = ForCens$meta$Latitude)

# get temperature data. Change here between HadISST and ERSST
# ForCenS SST
ForCens.SST <- readRDS('forcens_HadSST_1870-1899.RDS')
# traps historical SST
TS.SST <- readRDS('traps_HadSST_1870-1899.RDS')
TS.SST <- TS.SST[names(TS.SST) %in% names(dat.sel)]
# traps SST during collection
TS.real.SST <- readRDS('traps_HadSST_period.RDS')
TS.real.SST <- TS.real.SST[names(TS.real.SST) %in% names(dat.sel)]

# load ForCenS data split by region
species_domains <- readRDS('species_domains_compare.RDS')
domain_indeces <- readRDS('domain_indeces_compare.RDS')

# load annual fluxes and assemblages 
# made using make_annual_assemblages.R
TS.ann <- readRDS('TS.ann.RDS')
# assemblage for each year
TS.assem <- readRDS('TS.assem.RDS')

# get meta for each flux time series
TS.xy <- as.data.frame(t(sapply(dat.sel, function(i) c(x = i$meta$lon_dec[1], y= i$meta$lat_dec[1]))))
TS.region <- sapply(dat.sel, function(x) x$meta$ocean[1])
TS.waterdepth <- sapply(dat.sel, function(x) x$meta$water_depth_m[1])
TS.upwelling <- sapply(dat.sel, function(x) x$meta$upwelling[1])
TS.lat <- sapply(dat.sel, function(x) x$meta$lat_dec[1])
TS.trapdepth <- sapply(dat.sel, function(x) x$meta$trap_depth_m[1])
TS.trapdepth <- sapply(lapply(regmatches(TS.trapdepth, gregexpr("[[:digit:]]+", TS.trapdepth)), as.numeric), mean) # get average
TS.size <- ifelse(sapply(dat.sel, function(x) x$meta$min_size_micron[1]) <150, 'small', 'large')


# assemblage of complete time series, long term integrated assemblage
LT.assem <- lapply(TS.ann, function(x){
  domain <- x$meta$ocean
  spp <- names(species_domains[[which(names(species_domains) %in% domain)]]$species)
  flux <- lapply(x$flux, function(i) i[names(i) %in% spp])
  total.flux <- do.call(rbind, flux)
  as.data.frame(t(colSums(total.flux)/sum(total.flux)))
})

# calculate physical distances from traps to each point in ForCenS
TS.dist <- apply(TS.xy, 1, function(x) distGeo(p1 = ForCens.xy, p2 = x)/1000)
# what is minimum distance
TS.dist.min <- apply(TS.dist, 2, min)
# index of nearest sample
TS.dist.min.indx <- apply(TS.dist, 2, which.min)

# square chord distance to nearest sample
# Merge TS and sed assemblages, assume 0 for species not in traps (standard MAT approach)
# use Merge from rioja, not from HMisc
sed.assem.nearest <- lapply(TS.dist.min.indx, function(x) ForCens$species[x,])
f <- mapply(function(sed, trap) rioja::Merge(sed, trap, split = TRUE), sed = sed.assem.nearest, trap = TS.assem, SIMPLIFY = FALSE)
sqcd.nearest <- lapply(f, function(x) as.vector(paldist2(x$sed, x$trap, dist.method = 'sq.chord')))

# square chord distance from LT assem to nearest sample
# Merge LT trap and sed assemblages, assume 0 for species not in traps (standard MAT approach)
# can use Merge from rioja
f.lt <- mapply(function(sed, trap) rioja::Merge(sed, trap, split = TRUE), sed = sed.assem.nearest, trap = LT.assem, SIMPLIFY = FALSE)
sqcd.nearest.lt <- sapply(f.lt, function(x) as.vector(paldist2(x$sed, x$trap, dist.method = 'sq.chord')))

# minumum square chord distance
get.sqcd <- function(spp.trap, domain){
  spp.forcens <- species_domains[[which(names(species_domains) == domain)]]$species
  f <- rioja::Merge(spp.trap, spp.forcens, split = TRUE)
  sqcd <- paldist2(f$spp.trap, f$spp.forcens, dist.method = 'sq.chord')
  result <- as.data.frame(t(apply(sqcd, 1, function(y) c(min(y), which.min(y)))))
  names(result) <- c('min.sqcd', 'indx.sqcd')
  result
}
sqcd.min <- mapply(function(x, y) get.sqcd(x, y), x = TS.assem, y = TS.region, SIMPLIFY = FALSE)
sqcd.min.lt <- mapply(function(x, y) get.sqcd(x, y), x = LT.assem, y = TS.region)

# extract sample species and metadata for sqcd.min ####
# for long term only
samples.sqcd.min <- mapply(function(indx, domain){
  list(
  spp = species_domains[[which(names(species_domains) == domain)]]$species[indx,],
  meta = species_domains[[which(names(species_domains) == domain)]]$meta[indx,])
}, indx = sqcd.min.lt[2,], domain = TS.region, SIMPLIFY = FALSE)

# position of most similar sediment samples
sqcd.min.xy <- t(sapply(samples.sqcd.min, function(x) c(x$meta$Longitude, x$meta$Latitude)))

# distance between trap and most similar sediment sample
spatial.distance.sqcd.min <- distGeo(TS.xy, sqcd.min.xy)/1000
latitudinal.distance.sqcd.min <- distGeo(TS.xy, cbind(TS.xy[,1], sqcd.min.xy[,2]))/1000

# relate dissimilarity to distance ####
# dissimilarity  and distance between nearest sediment sample and all forcens samples
sqcd.forcens.nearest <- sapply(sed.assem.nearest, function(x) paldist2(x, ForCens$species, dist.method = 'sq.chord'))
nearest.xy <- ForCens.xy[TS.dist.min.indx, ]
nearest.dist <- apply(nearest.xy, 1, function(x) distGeo(p1 = ForCens.xy, p2 = x)/1000)

# latitudinal distance between nearest sediment sample and forcens
nearest.y <- cbind.data.frame(x = 0, y = ForCens.xy$y[TS.dist.min.indx])
ForCens.y <- cbind.data.frame(x = 0, y = ForCens$meta$Latitude)
nearest.dist.lat <- apply(nearest.y, 1, function(x) distGeo(p1 = ForCens.y, p2 = x)/1000)

# for model fits
new.x <- seq(0, 2, by = 0.05)

# get data that need modelling
# select within radius all latitudinal distances and sqcd
radius <- 5000
lat.extent <- 2500
mod.dat <- list()
for(i in 1:ncol(nearest.dist.lat)){
  site.indx <- intersect(which(nearest.dist[,i] <= radius), which(nearest.dist.lat[,i] <= lat.extent))
  mod.dat[[i]] <- cbind.data.frame(x = sqcd.forcens.nearest[site.indx,i],
                                   y = nearest.dist.lat[site.indx,i])
}

# linear model, forced through origin
lin.mod <- lapply(mod.dat, function(z) lm(y~x -1, data = z))
lin.dist <- sapply(lin.mod, function(x) predict(x, data.frame(x = new.x)))

lat.displacement <- mapply(function(x,y) predict(x, data.frame(x = y)), x = lin.mod, y = sqcd.nearest.lt)

round(median(lat.displacement), 0)
round(range(lat.displacement), 0)

# SST differences ####
SST.nearest <- ForCens.SST[TS.dist.min.indx]
dSST.closest <- SST.nearest - TS.SST

get.SST.most.sim <- function(trap, domain){
  SST.indx <- trap$indx.sqcd
  ForCens.SST[domain_indeces[[which(names(domain_indeces) == domain)]]][SST.indx]
}
SST.most.sim <- mapply(function(x, y) get.SST.most.sim(trap = x, domain = y), x = sqcd.min, y = TS.region, SIMPLIFY = FALSE)
dSST.most.sim <- mapply(function(x, y) x - y, x = SST.most.sim, y = TS.SST)

# for the long term integrated flux
get.SST.most.sim.lt <- function(trap, domain){
  SST.indx <- trap
  ForCens.SST[domain_indeces[[which(names(domain_indeces) == domain)]]][SST.indx]
}

SST.most.sim.lt <- mapply(function(x, y) get.SST.most.sim.lt(trap = x, domain = y), x = sqcd.min.lt[2,], y = TS.region)
dSST.most.sim.lt <- SST.most.sim.lt - TS.SST

# get temperature difference
real.dT <- sapply(TS.real.SST, function(x) x$mean.TS.SST) - TS.SST

# get proportion of annual dT consistent with long term mean ####
indx.long <- which(sapply(dSST.most.sim, function(x) length(x))>1)
cons.an.lt <- mapply(function(i, j) i == j, i =lapply(dSST.most.sim[indx.long], function(x) x>0), j = dSST.most.sim.lt[indx.long]>0)
prop.cons.an.lt <- sapply(cons.an.lt, function(x) sum(x, na.rm = TRUE)/length(x[!is.na(x)]))*100
length.long <- sapply(cons.an.lt, length)
melt.prop.cons <- cbind.data.frame(consistent = prop.cons.an.lt, trap = names(prop.cons.an.lt), length = length.long)
melt.prop.cons <- melt.prop.cons[order(melt.prop.cons$length),]

# weighted mean of proportion consistent with long term
#wtm.prop.cons <- wtd.mean(prop.cons.an.lt, length.long)
#wtsd.prop.cons <-sqrt(wtd.var(prop.cons.an.lt, length.long))

# proportion consistent as a function of minimum time series length
interAnnual <- rbind.data.frame(
  cbind.data.frame(minLength = 2:5,
                   Average = sapply(2:5, function(x)
                     wtd.mean(subset(melt.prop.cons, length >= x)$consistent, subset(melt.prop.cons, length >= x)$length)),
                   SD = sapply(2:5, function(x)
                     sqrt(wtd.var(subset(melt.prop.cons, length >= x)$consistent, subset(melt.prop.cons, length >= x)$length))),
                   N = sapply(2:5, function(x)
                     length(subset(melt.prop.cons, length >= x)$consistent)),
                   group = 'weighted'
  ),
  cbind.data.frame(minLength = 2:5,
                   Average = sapply(2:5, function(x) 
                     mean(subset(melt.prop.cons, length >= x)$consistent)),
                   SD = sapply(2:5, function(x)
                     sd(subset(melt.prop.cons, length >= x)$consistent)),
                   N = sapply(2:5, function(x)
                     length(subset(melt.prop.cons, length >= x)$consistent)),
                   group = 'unweighted')
)

# assess inter annual variability excluding years with >25 % imputation
# proportion of annual flux imputed <= 0.25 year
indx.no.impute.all <- lapply(TS.ann[indx.long], function(x){
  sapply(x$flux, function(y) sum(y$imputed)/length(y$imputed) <= 0.25)
})
# which are still > 1 year
indx.no.impute <- indx.no.impute.all[sapply(indx.no.impute.all, sum)>1]

# proportion consistent with long term mean, excluding the years where > 25% was imputed
cons.an.lt.no.impute <- mapply(function(x, y) x[y], x = cons.an.lt[sapply(indx.no.impute.all, sum)>1], y = indx.no.impute)
length.long.no.impute <- sapply(cons.an.lt.no.impute, length)
prop.cons.an.lt.no.impute <- sapply(cons.an.lt.no.impute, function(x) sum(x, na.rm = TRUE)/length(x[!is.na(x)]))*100
melt.prop.cons.no.impute <- cbind.data.frame(consistent = prop.cons.an.lt.no.impute, trap = names(prop.cons.an.lt.no.impute), length = length.long.no.impute)
melt.prop.cons.no.impute <- melt.prop.cons.no.impute[order(melt.prop.cons.no.impute$length),]

wtm.prop.cons.no.impute <- wtd.mean(prop.cons.an.lt.no.impute, length.long.no.impute)
wtsd.prop.cons.no.impute <- sqrt(wtd.var(prop.cons.an.lt.no.impute, length.long.no.impute))

T1 <- subset(melt.prop.cons.no.impute, length >= 2)
T1$N <- 2
T2 <- subset(melt.prop.cons.no.impute, length >= 3)
T2$N <- 3
T3 <- subset(melt.prop.cons.no.impute, length >= 4)
T3$N <- 4
T4 <- subset(melt.prop.cons.no.impute, length >= 5)
T4$N <- 5

# these are the data needed for figure 4A
interAnnual.no.impute <- rbind.data.frame(T1, T2, T3, T4)

# get number of consistent time series per year ####
an.real.dT <- mapply(function(x, y) x$an.mean.TS.SST- y, x = TS.real.SST, y = TS.SST)
an.cons <- mapply(function(x, y) x=y, x= unlist(an.real.dT)>0, y = unlist(dSST.most.sim)>0)
n.cons.year <- table(cbind.data.frame(year = factor(unlist(sapply(TS.ann, function(x) x$year)), levels = seq(1978, 2013)),
                       consistent = an.cons))

# set up plotting data
plot.LT <- cbind.data.frame(trap = names(dSST.most.sim),
                            x = TS.xy$x,
                            y = TS.xy$y,
                            lt.dSST = dSST.most.sim.lt,
                            real.dT = real.dT,
                            deep = factor(sapply(TS.waterdepth, function(x) if(x > 2000){'deep'}else{'shallow'})),
                            far = factor(sapply(TS.dist.min, function(x) if(x > 250){'far'}else{'near'})),
                            trap.trend = factor(sapply(dSST.most.sim.lt, function(x) if(x > 0){'warming'}else{'cooling'})),
                            real.trend = factor(sapply(real.dT, function(x) if(x > 0){'warming'}else{'cooling'})),
                            sqcd.nearest = sqcd.nearest.lt,
                            sqcd.min = unlist(sqcd.min.lt[1,]),
                            duration = sapply(TS.assem, nrow),
                            lat.displacement = lat.displacement,
                            shell_size = TS.size
)
plot.LT$group <- 'WW'
plot.LT$group[plot.LT$real.dT > 0 & plot.LT$lt.dSST <= 0] <- 'WC'
plot.LT$group[plot.LT$real.dT <= 0 & plot.LT$lt.dSST > 0] <- 'CW'
plot.LT$group[plot.LT$real.dT <= 0 & plot.LT$lt.dSST <= 0] <- 'CC'
plot.LT$consistent <- plot.LT$group == 'WW' | plot.LT$group == 'CC'

saveRDS(plot.LT, 'plot_LT_HadSST_1870-1899.RDS')
