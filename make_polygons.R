# make circles with 100 km radius centered at trap and ForCenS sites
# used to extract SST from raster data
# Lukas Jonkers, last edit 22 Mar 2019

source('make_polygon_function.R')

dat.sel <- readRDS('dat_sel.RDS')
trap.xy <- as.data.frame(t(sapply(dat.sel, function(i) c(x = i$meta$lon_dec[1], y= i$meta$lat_dec[1]))))
trap.polys <- apply(trap.xy, 1, function(x) get_polygon_function(x, 1e5))

plot(NA, ylim = c(-90, 90), xlim = c(-180, 180))
lapply(trap.polys, function(x) plot(x, add= TRUE))

saveRDS(trap.polys, 'trap_polys.RDS')

ForCenS <- readRDS('forcens_trimmed_compare.RDS')
ForCenS.xy <- cbind.data.frame(x = ForCenS$meta$Longitude, y = ForCenS$meta$Latitude)
rm(ForCenS)

ForCenS.polys <- apply(ForCenS.xy, 1, function(i) get_polygon_function(i, 1e5))

plot(NA, ylim = c(-90, 90), xlim = c(-180, 180))
lapply(ForCenS.polys, function(x) plot(x, add= TRUE))

saveRDS(ForCenS.polys, 'ForCenS_polys.RDS')