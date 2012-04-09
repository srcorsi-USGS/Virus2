##################################################################
# Interpolate wastewater FIB data to match samples
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
WWFIB.dir <- paste(Rlocal,"/MMSD_virus/Virus2/FIB",sep="")
conc.dir <- paste(Rlocal,"/MMSD_virus/Virus2/Final compiled data/Concentrations",sep="")
dfWWFIB <- read.delim(paste(WWFIB.dir,"/FIBWW.txt",sep=""))
dfPathFIB <- read.delim(paste(conc.dir,"/PathWWFIB.txt",sep=""))


# Take mean of all samples collected on the same date: one from JI and one from SS.
WWdates <- unique(as.character(dfWWFIB$Sample.Date))
FIBcols <- names(dfWWFIB)[c(4,5,7,8)]
for (i in 1:length(WWdates)){
  subdf <- subset(dfWWFIB,dfWWFIB$Sample.Date==WWdates[i])
  if(i==1) WWmean <- as.data.frame(t(mean(subdf[,FIBcols])))
  else WWmean <- rbind(WWmean,as.data.frame(t(mean(subdf[,FIBcols]))))
}
WWmean <- cbind(WWdates,WWmean)

#Generate POSIXct dates and subtract 3.5 days for mean of sample period. PATHOGEN FILE DATES
#WERE BEGINNING OF SAMPLING AND FIB FILE DATES WERE END. SO SUBTRACT HERE AND ADD FOR PATHOGENS
WWmean$pdate <- as.POSIXct(strptime(WWmean$WWdates,format= "%d-%b-%y",tz="CST6CDT")) -3.5*24*3600
dfPathFIB$Ebpdate <- as.POSIXct(strptime(dfPathFIB$Ebpdate,format=  "%Y-%m-%d %H:%M",tz="GMT"))
WWmean$pdate <- strptime(format(as.POSIXct(WWmean$pdate),tz="GMT",usetz=TRUE),format="%Y-%m-%d %H:%M",tz="GMT")
nsamples <- length(dfPathFIB$Ebpdate)

for (i in 1:nsamples){

  WWcols <- ncol(WWmean)
  if(dfPathFIB$Ebpdate[i]<min(WWmean$pdate)) {
    WWinterp1 <- WWmean[which(WWmean$pdate==min(WWmean$pdate)),2:(WWcols-1)]
  } else{
    #Determine closest date of WW sample previous to stream sample
    subWW <- subset(WWmean,pdate<dfPathFIB$Ebpdate[i])
    WWmean1 <- subWW[which(subWW$pdate==max(subWW$pdate)),2:WWcols]

    #Determine closest date of WW sample after to stream sample
    subWW <- subset(WWmean,pdate>=dfPathFIB$Ebpdate[i])
    WWmean2 <- subWW[which(subWW$pdate==min(subWW$pdate)),2:WWcols]
    Wcols2 <- ncol(WWmean2)
    
    #Interpolate concentrations of WW to sample date
    dtimeWW <- as.numeric(difftime(WWmean2$pdate,WWmean1$pdate,units="days"))
    dtimeS1 <- as.numeric(difftime(dfPathFIB$Ebpdate[i],WWmean1$pdate,units="days"))
    timefraction <- dtimeS1/dtimeWW
    WWinterp1 <- WWmean1[1:(Wcols2-1)] + (WWmean2[1:(Wcols2-1)]-WWmean1[1:(Wcols2-1)])*timefraction
  }
  
  #Add interpolated values to data fram
  if(i == 1) WWinterp <- WWinterp1
  else WWinterp <- rbind(WWinterp,WWinterp1)
}

names(WWinterp) <- paste("WWFIB.",names(WWinterp),sep="")
dfPathFIBWW <- cbind(dfPathFIB,WWinterp)

write.table(dfPathFIBWW,paste(conc.dir,"/PathWWFIBWW.txt",sep=""),sep="\t",row.names=F)

df<-dfPathFIBWW
plot(log10(df$BacHumancen),log10(df$WWFIB.BacHuman..CN.100ml.))
abline(lm(log10(df$WWFIB.BacHuman..CN.100ml.)~log10(df$BacHumancen)))

datlowess <- lowess(log10(df$BacHumancen),log10(df$WWFIB.BacHuman..CN.100ml.),f=0.5)
lines(datlowess,col=colors()[26], lwd=2) #lowess smooth (red) with 10% window
