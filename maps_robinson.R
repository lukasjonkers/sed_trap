# make map with ggplot and Robinson projection
# based on http://rpsychologist.com/working-with-shapefiles-projections-and-world-maps-in-ggplot
# Lukas Jonkers, last edit 21 Mar 2019

rm(list = ls())
library(rgdal)
library(ggplot2)
library(sp)
library(reshape2)
library(viridis)

# get necessary shapefiles and reproject ####
# download shapefiles from http://www.naturalearthdata.com/downloads/110m-physical-vectors/

# read the shapefile for the simple worldmap
wmap <- readOGR(dsn = path.expand('~/Documents/worldshapefiles/ne_110m_land'), layer = 'ne_110m_land')
# convert to dataframe
wmap_df <- fortify(wmap)

# same for graticules
grat <- readOGR(path.expand('~/Documents/worldshapefiles/ne_110m_graticules_all'), layer="ne_110m_graticules_15")
grat_df <- fortify(grat)

# and bounding box
bbox <- readOGR(path.expand('~/Documents/worldshapefiles/ne_110m_graticules_all'), layer="ne_110m_wgs84_bounding_box")
bbox_df<- fortify(bbox)

# reproject all
wmap_robin <- spTransform(wmap, CRS("+proj=robin")) # reproject wmap
wmap_df_robin <- fortify(wmap_robin)
grat_robin <- spTransform(grat, CRS("+proj=robin"))  # reproject graticule
grat_df_robin <- fortify(grat_robin)
bbox_robin <- spTransform(bbox, CRS("+proj=robin"))  # reproject bounding box
bbox_robin_df <- fortify(bbox_robin)

# get xy for traps and ForCens and reproject ####
ForCens <- readRDS('forcens_trimmed_compare.RDS')
Fxy <- cbind.data.frame(x = ForCens$meta$Longitude, y = ForCens$meta$Latitude)
coordinates(Fxy) <- c("x", "y") # convert to spatialpointsdataframe
proj4string(Fxy) <- CRS("+proj=longlat + datum=WGS84")
F_robin <- spTransform(Fxy, CRS('+proj=robin')) # reproject points
F_robin_df <- as.data.frame(F_robin) # and convert coordinates to df

Txy <- readRDS('plot_LT_HadSST_1870-1899.RDS')
coordinates(Txy) <- c("x", "y") # convert to spatialpointsdataframe
proj4string(Txy) <- CRS("+proj=longlat + datum=WGS84")
T_robin <- spTransform(Txy, CRS('+proj=robin'))
T_robin_df <- as.data.frame(T_robin)


# Figure 1 ####
library(maptools)
gpclibPermit()
gpclibPermitStatus()

# create a blank ggplot theme
theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         panel.border = element_blank(),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         plot.title = element_text(size=12),
                         legend.key = element_blank(),
                         rect = element_rect(fill = 'transparent')
)
)

# plot decadal linear SST trend 1870-2015
# based on http://mazamascience.com/WorkingWithData/?p=1494
# get SST trend from HadISST
SST.slope <- as.data.frame(readRDS('hadisst_trend_1870-2015.RDS')[,-4])
# get coordinates
coordinates(SST.slope) <- c("x", "y")
proj4string(SST.slope) <- CRS("+proj=longlat + datum=WGS84")
# make polygons
polys = as(SpatialPixelsDataFrame(SST.slope, SST.slope@data, tolerance = 0.149842),"SpatialPolygonsDataFrame")
# reproject map
polys <- spTransform(polys, CRS("+proj=robin"))
# add to data a new column termed "id" composed of the rownames of data
polys@data$id <- rownames(polys@data)
# create a data frame from the spatial object
polysPoints <- fortify(polys, region = "id")
# merge the "fortified" data with the data from our spatial object
slopeDF <- merge(polysPoints, polys@data, by = "id")

# replace NA with zero and add continents for plotting
temp <- slopeDF
temp$value[is.na(temp$value)] <- 0
wmap_df_robin$value <- NA

tem <- rbind.data.frame(temp, wmap_df_robin)
#saveRDS(tem, 'polys4map.RDS')

Fig1 <- ggplot(bbox_robin_df, aes(long,lat)) +
          geom_polygon(data = tem, aes(x = long, y = lat, group = group, fill = value, colour = value)) +
          scale_fill_gradient2(low = 'dodgerblue4', high = 'firebrick3', na.value = 'grey80', name = expression(~degree~C~decade^-1)) +
          scale_colour_gradient2(low = 'dodgerblue4', high = 'firebrick3', na.value = 'grey80', guide = 'none') +
          geom_point(data = F_robin_df, aes(x, y), colour = 'grey40', alpha = 0.5, size = 0.05, shape = 16) +
          geom_point(data = T_robin_df, aes(x, y), colour = 'grey20', size = 1, alpha = 0.8) +
          geom_polygon(linetype = 'solid', fill = NA, colour = 'grey80', size = 0.3) +
          coord_equal() +
          theme(legend.position = 'right') +
          theme_opts +
          theme(text = element_text(size = 7),
                legend.key.width = unit(0.2, 'cm'),
                legend.key.height = unit(0.5, 'cm'),
                legend.margin = margin(0,0,0,0),
                legend.box.margin = margin(-10, 0,0, -15))

ggsave('~/Dropbox/projects_ongoing/traps_sed/Figs_paper/Fig1.pdf',
       plot = Fig1,
       units = 'cm',
       width = 9,
       height = 6)



# Figure 3A ####
# plot trap sites on map, colour by consistency in trend
library(maptools)
gpclibPermit()
gpclibPermitStatus()

# plot mean SST over sediment trap period (1978-2013)
# based on http://mazamascience.com/WorkingWithData/?p=1494
# get SST from HadISST
SST.mean <- as.data.frame(readRDS('hadisst_mean_1978-2013.RDS')[,-4])
# get coordinates
coordinates(SST.mean) <- c("x", "y")
proj4string(SST.mean) <- CRS("+proj=longlat + datum=WGS84")
# make polygons
polysMean = as(SpatialPixelsDataFrame(SST.mean, SST.mean@data, tolerance = 0.149842),"SpatialPolygonsDataFrame")
# reproject map
polysMean <- spTransform(polysMean, CRS("+proj=robin"))
# add to data a new column termed "id" composed of the rownames of data
polysMean@data$id <- rownames(polysMean@data)
# create a data.frame from the spatial object
polysPoints <- fortify(polysMean, region = "id")
# merge the "fortified" data with the data from our spatial object
meanDF <- merge(polysPoints, polysMean@data, by = "id")

# add continents for plotting
wmap_df_robin$value <- NA
tem <- rbind.data.frame(meanDF, wmap_df_robin)
#saveRDS(tem, 'polys4map.RDS')

# tweak T_robin to allow colour and fill
T_robin_df$fill <- ifelse(T_robin_df$trap.trend == 'warming', 'firebrick3', 'dodgerblue4')
T_robin_df$colour <- ifelse(T_robin_df$real.trend == 'warming', 'firebrick3', 'dodgerblue4')

# set colourscheme
cols <- viridis(5, option = 'viridis')

Fig3A <- ggplot(bbox_robin_df, aes(long,lat)) +
  geom_polygon(data = tem, aes(x = long, y = lat, group = group, fill = value, colour = value)) +
  scale_fill_gradientn(colours = cols, na.value = 'grey80', name = 'SST [degC]') +
  scale_colour_gradientn(colours = cols, na.value = 'grey80', guide = 'none') +
  geom_point(data = T_robin_df, aes(x, y), colour = 'white', fill = NA, size = 1.5, alpha = 0.8) +
  geom_point(data = T_robin_df, aes(x, y), colour = T_robin_df$colour, fill = alpha(T_robin_df$fill, 0.8), shape = 21, size = 1, stroke = 0.5) +
  geom_polygon(linetype = 'solid', fill = NA, colour = 'black', size = 0.3) +
  coord_equal() +
  theme(legend.position = 'bottom') +
  theme_opts +
  theme(text = element_text(size = 7),
        #legend.box.spacing = unit(1, 'cm'),
        legend.key.width = unit(0.7, 'cm'),
        legend.key.height = unit(0.1, 'cm'),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10, 0,0, 0))



ggsave('~/Dropbox/projects_ongoing/traps_sed/Figs_paper/Fig3A_viridis.pdf',
       plot = Fig3A,
       units = 'cm',
       width = 9,
       height = 6)