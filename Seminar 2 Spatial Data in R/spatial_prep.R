#Layer prep for Spatial Session#

rivers <- readOGR('NHDNamedRiversStreams.shp')
plot(rivers)
proj4string(rivers)
elk_utm@bbox

rivers <- spTransform(rivers, CRS('+init=epsg:32612'))
x_min <- 604780.8 - 20000
x_max <- 631893.1 + 20000
y_min <- 4813177.7 - 20000
y_max <- 4850037.7 + 20000
coords = matrix(c(x_min, y_min, x_max, y_min, x_max, y_max, x_min, y_max, x_min, y_min), 
                ncol = 2, byrow = TRUE)


P1 = Polygon(coords)
Ps1 = SpatialPolygons(list(Polygons(list(P1), ID = "a")), proj4string=CRS('+init=epsg:32612'))
plot(Ps1, axes = TRUE); plot(rivers, col = 'blue', add=T); plot(elk_utm, col = 'red', add=T)

sub_rivers <- crop(rivers, Ps1)
plot(sub_rivers, col = 'green', add=T)
writeOGR(sub_rivers, "wiggins_rivers.shp", driver = "ESRI Shapefile", layer = "rivers")

#polygon
pub_lands <- readOGR('PublicLands_MTWY.shp')
pub_lands <- spTransform(pub_lands, 
                         CRS("+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

proj4string(pub_lands)
sub_lands <- crop(pub_lands, Ps1)
plot(Ps1, axes = TRUE); plot(pub_lands, col = 'blue', add=T); plot(elk_utm, col = 'red', add=T)
plot(sub_lands, col = 'green', add=T)
writeOGR(sub_lands, "wiggins_lands.shp", driver = "ESRI Shapefile", layer = "lands")

#elk_conus
elk_conus <- spTransform(elk_data, CRS('+init=epsg:5070'))
elk_conus <- as.data.frame(elk_conus)
elk_conus <- elk_conus[,-1]
colnames(elk_conus)[4] <- 'X'
colnames(elk_conus)[5] <- 'Y'
write.csv(elk_conus, 'elk_other_proj.csv')
