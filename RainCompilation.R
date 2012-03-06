# Menomonee virus Rainfall compilation
#


###################
#Set R directory
SNsys <- system("wmic bios get serialnumber", intern = TRUE, show.output.on.console = FALSE)
SN <- gsub("\\s","",SNsys)[2]
if(SN == "R97R0NW") {Wlocal <- "D:/srcldata" #Laptop
                     Project <- "D://srcldata"
}else {Wlocal <- "//igsarmewwssrc/SRCdata"
       Project <- "M:/QW Monitoring Team"} #Network
Rlocal <- paste(Wlocal,"/R",sep="")
excel <- "%m/%d/%Y %H:%M"

source(paste(Rlocal,"/Rainmaker/RRainmaker.R",sep=""))
source(paste(Rlocal,"/Rainmaker/fxn_RMeventsSamples.R",sep=""))
source(paste(Rlocal,"/Rainmaker/fxn_RMevents.plotQ.R",sep=""))

setwd(paste(Rlocal,"/MMSD_virus/Virus2/Rain",sep=""))
#setwd(paste(Project,sep=""))

rain.site <- c("MenoFalls","LMDonges","Underwood","Honey","70th","16th")
site <- c("Meno Falls","Donges Bay","Underwood","Honey","70th St.","16th St.")

dir.prec <- paste(Project,"/MMSD/Viruses/Menomonee Viruses/Data compilation/EnDDat precip/",sep="")

for (i in 1:length(rain.site)){
  dir.Q <- paste(Project,"/MMSD/Viruses/Menomonee Viruses/Data compilation/UnitValues/mmsd/",sep="")
  
  Sitefile <- paste(rain.site[i],"_Final.txt",sep="")
  Qfile <- paste(rain.site[i],".rdb",sep="")

  dfRain <- read.delim(paste(dir.prec,Sitefile,sep=""))
  
  dfQ <- read.delim(paste(dir.Q,Qfile,sep=""))
  dfQ <- subset(dfQ,substr(dfQ$NAME,1,1)=="Q" & Q != -1.23456e+25)


  hr <- trunc(dfQ$MINUTE/60)
  min <- dfQ$MINUTE-hr*60
  
  GMToffset <- 5 # hours to central standard time
  dfQ$pdate <- strptime(paste(dfQ$YEAR,"-",dfQ$MONTH,"-",dfQ$DAY," ",hr,":",min,sep=""),format="%Y-%m-%d %H:%M")
  dfRain$pdate <- strptime(dfRain$GMT.Time,format="%m/%d/%Y %H:%M",tz="") - GMToffset*60*60
#  dfRain$pdate2 <- as.POSIXlt(dfRain$pdate,isdst=1)
  
  dfRain <- dfRain[which(!is.na(dfRain$pdate)),]

  dfsamples.all <- read.delim(paste(Rlocal,"/MMSD_virus/Virus2/virus2Prelim5.txt",sep=""))
  dfsamples <- subset(dfsamples.all,Site.Name==site[i])
  dfhydro <- read.delim(paste(dir.Q,rain.site[i],"/",rain.site[i],"-storm-volume-2009-11_FINAL.txt",sep=""))
  dfsamples2 <- merge(dfsamples,dfhydro,by.x="Study.Sample.ID",by.y="USGS.ID",all=T)
  
  dfsamples2$start_date <- as.character(dfsamples2$start_date)
  dfsamples2$end_date <- as.character(dfsamples2$end_date)
  
  for (i in 1:nrow(dfsamples2)){
    if(dfsamples2$Sample.Type.x[i]=="Baseflow"){
      dfsamples2$start_date[i] <- as.character(dfsamples2[i,"SSdate"])
      dfsamples2$end_date[i] <- as.character(dfsamples2[i,"SEdate"])
    }
  }
  
  dfsamples2$pdate <- strptime(dfsamples2$SSdate, format="%m/%d/%Y %H:%M")
  dfsamples2$ddate <- as.Date(dfsamples2$pdate)
  dfsamples2$Hbpdate <- strptime(dfsamples2$start_date, format="%m/%d/%Y %H:%M")
  dfsamples2$Hepdate <- strptime(dfsamples2$end_date, format="%m/%d/%Y %H:%M")

  dfsamples3 <- Hydrovol(dfQ=dfQ,Q="VALUE",time="pdate",df.dates=dfsamples2,bdate="Hbpdate",edate="Hepdate")
  
  #events <- RMevents(dfRain,ieHr=6,rainthresh=5.1,rain="rain",time="pdate")[1]

}


# Menomonee test data set
events.MF<- read.delim("D:/srcldata/R/MMSD_virus/Virus2/virus2Prelim5.txt")
events.MF$pdate <- strptime(events.MF$SSdate, format="%m/%d/%Y %H:%M")
events.MF$ddate <- as.Date(events.MF$pdate)
write.table(events.MF,"D:/srcldata/R/MMSD_virus/Virus2/virus2Prelim6.txt",sep="\t",row.names=F)

RMevents.plotQ(df=dfRain,dfQ=dfQ,date="pdate",Qdate="pdate",rain="rain",Q="VALUE",df.events=events.MF[1:2,],sdate="Braindate",edate="Eraindate",erain="depth",plot.buffer=6,site.name="MFtest3",)

str(events)

rm(events.MF)

dfsamples <- subset(dfVirus,Site.Name=="Meno Falls")
dfsamples$bpdate <- strptime(dfsamples[,"SSdate"], format="%m/%d/%Y %H:%M")
dfsamples$epdate <- strptime(dfsamples[,"SEdate"], format="%m/%d/%Y %H:%M")

SampleVols <- read.delim(paste(Project,"/MMSD/Viruses/Menomonee Viruses/Data compilation/UnitValues/mmsd/MenoFalls/menofalls-storm-volume-2009-11.txt",sep=""))
SampleVols$Ebpdate <- strptime(SampleVols$start_date, format="%m/%d/%Y %H:%M")
SampleVols$Eepdate <- strptime(SampleVols$end_date, format="%m/%d/%Y %H:%M")
#dfsamples <- dfsamples[,1:(ncol(dfsamples)-1)]
dfsamples2 <- merge(dfsamples,SampleVols,by.x="Study.Sample.ID",by.y="USGS.ID",all=T)

for (i in 1:nrow(dfsamples2)){
 if(dfsamples2$Sample.Type.x[i]=="Baseflow"){
   dfsamples2$Ebpdate[i] <- dfsamples2[i,bdate]
   dfsamples2$Eepdate[i] <- dfsamples2[i,edate]
 }
}


#   df <- dfRain
#   ieHr <- 6
#   rain <- "rain"
#   time <- "pdate"
#   dfsamples <- dfsamples
#   bdate <- "Ebpdate"
#   edate <- "Eepdate"

dfQ<-dfQ
Q<-"VALUE"
time<-"pdate"
df.dates<-dfsamples2
bdate<-"Hbpdate"
edate<-"Hepdate"

events.MF <- RMeventsSamples(df=dfRain,ieHr=12,rain="rain",time="pdate",dfsamples=dfsamples2,bdate="Ebpdate",edate="Eepdate")
write.table(events.MF,"D:/srcldata/R/MMSD_virus/Virus2/virus2Prelim7.txt",sep="\t",row.names=F)






for(i in 1:length(rws))
  as.character(dfsamples2[rws[i],"Event.Type"])==as.character(dfsamples2[rws[i],"Sample.Type.y"])
data.frame(x=dfsamples2[,"Event.Type"],y=dfsamples2[,"Sample.Type.y"])