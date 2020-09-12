// Load all packages first
library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales)   #for transparency

// Load shape file 
crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
states=readShapePoly("/Users/User/Documents/Hufford/Hufnagel/mexstates/mexstates.shp",proj4string=crswgs84,verbose=TRUE)
plot(states)

// For non-confidence parv pops
setwd("~/Documents/Hufford/Hufnagel/Chapter2/parvMap")

allData = read.table("nonConfPArvPops.info.pops") 
df = data.frame(accession=allData$V1, lat=allData$V2, long=allData$V3)
accession = df[,1]
lat = df[,2]
long = df[,3]

points(allData$V7, allData$V6, pch=as.numeric(allData$V2), col="slateblue4", cex=.1) #plot the sample sites
points(df$long, df$lat, pch=1, col="slateblue4", cex=.5)

// For high confidence parv. pops
allData = read.table("parvCoords.txt") #import data
df = data.frame(accession=allData$V1, lat=allData$V2, long=allData$V3)
accession = df[,1]
lat = df[,2]
long = df[,3]

points(allData$V7, allData$V6, pch=as.numeric(allData$V2), col="skyblue", cex=.1) #plot the sample sites
points(df$long, df$lat, pch=1, col="skyblue", cex=.5)

// For non-confidence mex populations
allData = read.table("nonConfMexPops.info.pops")
df = data.frame(accession=allData$V1, lat=allData$V2, long=allData$V3)
accession = df[,1]
lat = df[,2]
long = df[,3]

points(allData$V7, allData$V6, pch=as.numeric(allData$V2), col="yellow", cex=.1) #plot the sample sites
points(df$long, df$lat, pch=1, col="yellow", cex=.5)

// For mex parent populations
allData = read.table("mexCoords.txt")
df = data.frame(accession=allData$V1, lat=allData$V2, long=allData$V3)
accession = df[,1]
lat = df[,2]
long = df[,3]

points(allData$V7, allData$V6, pch=as.numeric(allData$V2), col="yellow4", cex=.1) #plot the sample sites
points(df$long, df$lat, pch=1, col="yellow4", cex=.5)

// For hybrids
allData = read.table("hybridCoords.txt")
df = data.frame(accession=allData$V1, lat=allData$V2, long=allData$V3)
accession = df[,1]
lat = df[,2]
long = df[,3]

points(allData$V7, allData$V6, pch=as.numeric(allData$V2), col="violetred1", cex=.1) #plot the sample sites
points(df$long, df$lat, pch=1, col="violetred1", cex=.5)


// Set legend
legend("bottomleft", c("Conf Parv","Ambig Parv","Conf Mex", "Ambig Mex", "hybrids"), col = c("slateblue4","skyblue","yellow4", "yellow", "violetred1"), pch = c(16,16,16))