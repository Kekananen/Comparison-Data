// Load all packages first
library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales)   #for transparency
install.packages("raster")
library(raster)

// Load shape file 
crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
states=readShapePoly("/Users/User/Documents/Hufford/Hufnagel/mexstates/mexstates.shp",proj4string=crswgs84,verbose=TRUE)
plot(states)

// Load tiles of elevation graphics
block1 <- getData('SRTM', lon=-105, lat=25)
block2 <- getData('SRTM', lon=-100, lat=20)
block3 <- getData('SRTM', lon=-95, lat=20)
block4 <- getData('SRTM', lon=-90, lat=20)
block5 <- getData('SRTM', lon=-110, lat=25)
block6 <- getData('SRTM', lon=-100, lat=25)
block7 <- getData('SRTM', lon=-95, lat=25)
block8 <- getData('SRTM', lon=-90, lat=25)
block9 <- getData('SRTM', lon=-120, lat=30)
block10 <- getData('SRTM', lon=-115, lat=30)
block11 <- getData('SRTM', lon=-105, lat=20)
block12 <- getData('SRTM', lon=-110, lat=30)
block13 <- getData('SRTM', lon=-100, lat=30)
block14 <- getData('SRTM', lon=-115, lat=35)
block15 <- getData('SRTM', lon=-110, lat=35)
block16 <- getData('SRTM', lon=-115, lat=25)
block17 <- getData('SRTM', lon=-105, lat=30)
block18 <- getData('SRTM', lon=-120, lat=35)
block19 <- getData('SRTM', lon=-105, lat=35)

// Mash tiles into one graphic 
blockMosaic <- mosaic(block1, block2, block3, block4, block5, block6, block7, block8, block9, block10, block11, block12, block13, block14, block15, block16, block17, block18, block19, fun=mean)

// Fit tiles to shape file
image <- crop(blockMosaic, states)
image <- mask(blockMosaic, states)

// Plot shape file and graphics
plot(image)
plot(states, add = TRUE)

//Plots the data points for our data
allData = read.table("nonConfPArvPops.info.pops") 
df = data.frame(accession=allData$V1, lat=allData$V2, long=allData$V3)
accession = df[,1]
lat = df[,2]
long = df[,3]

points(allData$V7, allData$V6, pch=as.numeric(allData$V2), col="slateblue4", cex=.1) #plot the sample sites
points(df$long, df$lat, pch=1, col="slateblue4", cex=.5)

allData = read.table("parvCoords.txt")
df = data.frame(accession=allData$V1, lat=allData$V2, long=allData$V3)
accession = df[,1]
lat = df[,2]
long = df[,3]

points(allData$V7, allData$V6, pch=as.numeric(allData$V2), col="skyblue", cex=.1) #plot the sample sites
points(df$long, df$lat, pch=1, col="skyblue", cex=.5)

allData = read.table("hybridCoords.txt")
df = data.frame(accession=allData$V1, lat=allData$V2, long=allData$V3)
accession = df[,1]
lat = df[,2]
long = df[,3]

points(allData$V7, allData$V6, pch=as.numeric(allData$V2), col="violetred1", cex=.1) #plot the sample sites
points(df$long, df$lat, pch=1, col="violetred1", cex=.5)

allData = read.table("nonConfMexPops.info.pops")
df = data.frame(accession=allData$V1, lat=allData$V2, long=allData$V3)
accession = df[,1]
lat = df[,2]
long = df[,3]

points(allData$V7, allData$V6, pch=as.numeric(allData$V2), col="green", cex=.1) #plot the sample sites
points(df$long, df$lat, pch=1, col="green4", cex=.5)

allData = read.table("mexCoords.txt")
df = data.frame(accession=allData$V1, lat=allData$V2, long=allData$V3)
accession = df[,1]
lat = df[,2]
long = df[,3]

points(allData$V7, allData$V6, pch=as.numeric(allData$V2), col="pink", cex=.1) #plot the sample sites
points(df$long, df$lat, pch=1, col="pink", cex=.5)

// Set legend
legend("bottomleft", c("Conf Parv","Ambig Parv","Conf Mex", "Ambig Mex", "hybrids"), col = c("skyblue","slateblue4","pink1", "green4", "violetred1"), pch = c(16,16,16))