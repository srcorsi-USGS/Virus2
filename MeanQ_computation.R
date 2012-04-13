#################################################################################
## Compute mean Q for specified period of time for Menomonee sites
## Retrieve data from EnDDaT and compute mean
#################################################################################


bdate <- "2009-03-01"
edate <- "2011-03-31"
site <- c("04087030","04087050","04087088","04087119","04087120","04087142")
abb <- c("MF","LD","UW","HW","MW","MC")
parameter <- "00060"
stat <- "00003" # output statistic: 3 for daily mean, 0 for unit values

Qmeans <- numeric()

for (i in 1:length(site)){
  myurl <- paste("http://cida.usgs.gov/enddat/service/execute?style=tab&DateFormat=Excel&fill=&beginPosition=",bdate,"&endPosition=",edate,"&Lake=michigan&TZ=0_GMT&BeachName=UW&BeachLat=43.054785&BeachLon=-88.046417&shapefile=&shapefileFeature=null&precipOrigin=&maxTime0=&maxTime1=&maxTime2=&filterId=&timeInt=6&NWIS=",site[i],"%3A",parameter,"%3A",stat,sep="")
  df <- read.delim(myurl)
  Qmeans[i] <- mean(df[,2],na.rm=T)
}



Qmeans
Qmeans/DA
