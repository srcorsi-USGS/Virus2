# Menomonee virus Rainfall compilation
#


###################
#Set R directory
SNsys <- system("wmic bios get serialnumber", intern = TRUE, show.output.on.console = FALSE)
SN <- gsub("\\s","",SNsys)[2]
if(SN == "R97R0NW") {Wlocal <- "D:/srcldata" #Laptop
}else {Wlocal <- "//igsarmewwssrc/SRCdata"} #Network
Rlocal <- paste(Wlocal,"/R",sep="")

source(paste(Rlocal,"/Rainmaker/RRainmaker.R",sep=""))

setwd(paste(Rlocal,"/MMSD_virus/Virus2/Rain",sep=""))

rain.site <- c("MenoFalls","LMDonges","Underwood","Honey","70th","16th")

dir.prec <- paste(Wlocal,"/MMSD/Viruses/Menomonee Viruses/Data compilation/EnDDat precip/",sep="")

for (i in 1:length(rain.site))[
  dir.Q <- paste(Wlocal,"/MMSD/Viruses/Menomonee Viruses/Data compilation/UnitValues/mmsd/",sep="")
  
  Sitefile <- paste(rain.site[i],"_Final.txt",sep="")
  Qfile <- paste(rain.site[i],".rdb",sep="")

dfRain <- read.delim(paste(dir.prec,Sitefile,sep=""))
dfQ <- read.delim(paste(dir.Q,Qfile,sep=""))
dfQ <- subset(dfQ,substr(dfQ$NAME,1,1)=="Q")

  hr <- trunc(dfQ$MINUTE/60)
  min <- dfQ$MINUTE-hr*60
  
  GMToffset <- 5 # hours to central standard time
  dfQ$pdate <- strptime(paste(dfQ$YEAR,"-",dfQ$MONTH,"-",dfQ$DAY," ",hr,":",min,sep=""),format="%Y-%m-%d %H:%M")
  dfRain$pdate <- strptime(dfRain$GMT.Time,format="%m/%d/%Y %H:%M",tz="") - GMToffset*60*60
  dfRain$pdate2 <- as.POSIXlt(dfRain$pdate,isdst=1)
  
  dfRain <- dfRain[which(!is.na(dfRain$pdate)),]


events <- RMevents(dfRain,ieHr=6,rainthresh=5.1,rain="rain",time="pdate")[1]

}


RMevents.plotQ(df=dfRain,df.events=event.test,sdate="StartDate",edate="EndDate",plot.buffer=6,site.name="MFtest",dfQ=dfQ,)

str(events)

