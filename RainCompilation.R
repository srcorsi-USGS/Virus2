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
source(paste(Rlocal,"/Hydrovol/Fxn_WQcompos.R",sep=""))

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
  
  for (j in 1:nrow(dfsamples2)){
    if(dfsamples2$Sample.Type.x[j]=="Baseflow"){
      dfsamples2$start_date[j] <- as.character(dfsamples2[j,"SSdate"])
      dfsamples2$end_date[j] <- as.character(dfsamples2[j,"SEdate"])
    }
  }
  
  dfsamples2$pdate <- strptime(dfsamples2$SSdate, format="%m/%d/%Y %H:%M")
  dfsamples2$ddate <- as.Date(dfsamples2$pdate)
  dfsamples2$Hbpdate <- strptime(dfsamples2$start_date, format="%m/%d/%Y %H:%M")
  dfsamples2$Hepdate <- strptime(dfsamples2$end_date, format="%m/%d/%Y %H:%M")

  #Run scripts: 1. Hydrovol, 2. RMeventsSamples, 3. WQcompos
  dfsamples3 <- Hydrovol(dfQ=dfQ,Q="VALUE",time="pdate",df.dates=dfsamples2,bdate="Hbpdate",edate="Hepdate")
  dfsamples4 <- RMeventsSamples(df=dfRain,ieHr=12,rain="rain",time="pdate",dfsamples=dfsamples3,bdate="Hbpdate",edate="Hepdate")
  dfsamples5 <- WQcompos(df.samples=dfsamples4,sampleID="USGS_ID",parms="Human_virus",volume="event.vol")
  cf2liters <- 28.316846592
  dfsamples5$Load_Human_virus <- dfsamples5$Human_virus*dfsamples5$event.vol*cf2liters
  dfsamples5$Daily_Human_virus <- dfsamples5$Load_Human_virus/dfsamples5$duration*24
    
  # write file with date in name
  filenm <- paste(Rlocal,"/MMSD_virus/Virus2/",rain.site[i],"/",rain.site[i],"_",format(Sys.time(), "%b%d%Y_%H%M"),".txt",sep="")
  write.table(dfsamples5,filenm,sep="\t",row.names=F)
  
}

## End of data compilation
################################################################################################################################

# Graph events vs baseflow for individual sites

df <- subset(dfsamples5,!is.na(Human_virus))

bwplot((df$Human_virus+0.01)~df$Site.Name|df$Sample.Type.x,data=df,
       xlab="",ylab="Pathogen concentrations (gc/L)",
       main="Event type: Sum of human viruses",
       scales = list(y = list(log=10)))

bwplot((df$Load_Human_virus+0.01)~df$Site.Name|df$Sample.Type.x,data=df,
       xlab="",ylab="Pathogen load (gc)",
       main="Event type: Sum of human virus load",
       scales = list(y = list(log=10)))

bwplot((df$Load_Human_virus+0.01)~df$Site.Name|df$Event.Type,data=df,
       xlab="",ylab="Pathogen load (gc)",
       main="Event type: Sum of human virus load",
       scales = list(y = list(log=10)))

bwplot((df$Daily_Human_virus+0.01)~df$Site.Name|df$Event.Type,data=df,
       xlab="",ylab="Pathogen load (gc)",
       main="Event type: Daily mean human virus load",
       scales = list(y = list(log=10)))


######## Graph loadings for conditional graphs based on sites and event types  #########

rain.site <- c("MenoFalls","LMDonges","Underwood","Honey","70th","16th")
site <- c("Meno Falls","Donges Bay","Underwood","Honey","70th St.","16th St.")

# Combine data from all sites
for (i in 1:length(site)){

  filenm <- paste(Rlocal,"/MMSD_virus/Virus2/",rain.site[i],"/",rain.site[i],"_Mar05.txt",sep="")
  dfsite <- read.delim(filenm)
  if(i==1) df <- dfsite
  if(i>1) df <- rbind(df,dfsite)
  
}

bwplot((df$Daily_Human_virus+0.01)~df$Site.Name|df$Event.Type,data=df,
       xlab="",ylab="Pathogen load (gc)",
       main="Event type: Daily mean human virus load",
       scales = list(y = list(log=10),x=list(rot=90)))

bwplot((df$Daily_Human_virus+0.01)~df$Event.Type|df$Site.Name,data=df,
       xlab="",ylab="Pathogen load (gc)",
       main="Event type: Daily mean human virus load",
       scales = list(y = list(log=10),x=list(rot=90)))

bwplot((df$Daily_Human_virus+0.01)~df$Sample.Type.x|df$Site.Name,data=df,
       xlab="",ylab="Pathogen load (gc)",
       main="Event type: Daily mean human virus load",
       scales = list(y = list(log=10),x=list(rot=90)))

bwplot((df$Load_Human_virus+0.01)~df$Sample.Type.x|df$Site.Name,data=df,
       xlab="",ylab="Pathogen load (gc)",
       main="Event type: Daily mean human virus load",
       scales = list(y = list(log=10),x=list(rot=90)))

plot(df$depth,(df$Human_virus+0.01),log="y")

m1 <- lm((Human_virus+0.01)~depth+Qmax,data=df)
summary(m1)

subdf <- subset(df,Event.Type != "Baseflow")
m2 <- lm((Human_virus+0.01)~depth+Qmax,data=subdf)
summary(m2)
plot(subdf$depth,(subdf$Human_virus+0.01),log="y")

subdf <- subset(subdf,Human_virus > 0)
m3 <- lm((Human_virus+0.01)~depth+Qmax,data=subdf)
summary(m3)
plot(subdf$depth,(subdf$Human_virus+0.01),log="y")

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