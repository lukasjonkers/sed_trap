# function to make polygons from longitude latitude pairs to use in raster::extract
# Lukas Jonkers, 22 Mar 2019

# in order to make circles coordinates need to be transformed to UTM, which has units = metres
# UTM needs zone and hemisphere
find_UTM_zone <- function(longitude, latitude) {
  
  # Special zones for Svalbard and Norway
  if (latitude >= 72.0 && latitude < 84.0 ) 
    if (longitude >= 0.0  && longitude <  9.0) 
      return(31);
  if (longitude >= 9.0  && longitude < 21.0)
    return(33)
  if (longitude >= 21.0 && longitude < 33.0)
    return(35)
  if (longitude >= 33.0 && longitude < 42.0) 
    return(37)
  
  (floor((longitude + 180) / 6) %% 60) + 1
}

get_polygon_function <- function(x, radius){
  require(sf)
  p <- as.numeric(as.matrix(x, nrow = 1))
  point <- st_point(p)
  point_WS84 <- st_sfc(point, crs = 4236)
  zone <- find_UTM_zone(x[1], x[2])
  EPSG_code <- as.integer(if(x[2]>0) {32600+zone} else{32700+zone})
  point_UTM <- st_transform(point_WS84, EPSG_code)
  poly_UTM <- st_buffer(point_UTM, radius)
  poly_WS84_tem <- st_transform(poly_UTM, 4326)
  poly_WS84 <- st_wrap_dateline(poly_WS84_tem, options = c("WRAPDATELINE=YES"))
  as(poly_WS84, "Spatial")
}





  
