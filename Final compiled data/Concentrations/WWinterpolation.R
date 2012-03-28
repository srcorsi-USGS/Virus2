##################################################################
# Interpolate wastewater data to match samples
##################################################################

###################
#Set R directory
SNsys <- system("wmic bios get serialnumber", intern = TRUE, show.output.on.console = FALSE)
SN <- gsub("\\s","",SNsys)[2]
if(SN == "R97R0NW") {Wlocal <- "D:/srcldata" #Laptop
                     Project <- "D://srcldata"
}else {Wlocal <- "//igsarmewwssrc/SRCdata"
       Project <- "M:/QW Monitoring Team"} #Network
Rlocal <- paste(Wlocal,"/R",sep="")

#Read in files
WW.dir <- paste(Project,"/MMSD/Viruses/Marshfield results/Jan 2012",sep="")
conc.dir <- paste(Rlocal,"/MMSD_virus/Virus2/Final compiled data/Concentrations",sep="")
dfWW <- read.delim(paste(WW.dir,"/DataSummarySewageAdVseparated.txt",sep=""))
dfPathFIB <- read.delim(paste(conc.dir,"/Path_FIB_WQ_Hy_All.txt",sep=""))


# Take mean of all samples collected on the same date: one from JI and one from SS.
dfWW$Collection.Start.Date <- as.character(dfWW$Collection.Start.Date)
WWdates <- unique(dfWW$Collection.Start.Date)
for (i in 1:length(WWdates)){
  subdf <- subset(dfWW,dfWW$Collection.Start.Date==WWdates[i])
  if(i==1) WWmean <- as.data.frame(t(mean(subdf[,5:26])))
  else WWmean <- rbind(WWmean,as.data.frame(t(mean(subdf[,5:26]))))
}
WWmean <- cbind(WWdates,WWmean)
names(WWmean)[1] <- "date"

#Generate POSIXct dates and add 3.5 days for mean of sample period
WWmean$pdate <- as.POSIXct(strptime(WWmean$date,format= "%d-%b-%y",tz="CST6CDT")) +3.5*24*3600
#convert to GMT
WWmean$pdate <- as.POSIXct(strptime(format(WWmean$pdate,tz="GMT",usetz=TRUE),format="%Y-%m-%d %H:%M",tz="GMT"))

dfPathFIB$Ebpdate <- as.POSIXct(strptime(dfPathFIB$Ebpdate,format=  "%m/%d/%Y %H:%M",tz="GMT"))
nsamples <- length(dfPathFIB$Ebpdate)
WWcols <- ncol(WWmean)
for (i in 1:nsamples){
  dtime <- difftime(dfPathFIB$Ebpdate[i],WWmean$pdate)
  minWWdate <- which(WWmean$pdate==min(WWmean$pdate))
  if(dfPathFIB$Ebpdate[i]<min(WWmean[minWWdate,"pdate"]) WWvalues <- WWmean[minWWdate,2:WWcols]
  else {
    WWsub <- WWmean[dtime>0,]
    WWsub[which(WWsub$pdate==max(WWsub$pdate)),"pdate"

  



