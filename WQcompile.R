## Script to compile real-time water quality and Q data into virus/FIB data from 
## Menomonee R Study

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

#Initialize needed functions
source(paste(Rlocal,"/Rainmaker/RRainmaker.R",sep=""))
source(paste(Rlocal,"/Rainmaker/fxn_RMeventsSamples.R",sep=""))
source(paste(Rlocal,"/Rainmaker/fxn_RMevents.plotQ.R",sep=""))
source(paste(Rlocal,"/Hydrovol/Fxn_WQcompos.R",sep=""))
source(paste(Rlocal,"/Scripts/fxn_TSstormstats.R",sep=""))

UV.site <- c("MenoFalls","LMDonges","Underwood","Honey","70th","16th")
site <- c("Meno Falls","Donges Bay","Underwood","Honey","70th St.","16th St.")

dir.UV <- "MMSD/Viruses/Menomonee Viruses/Data compilation/UnitValues/mmsd/"


for(site.num in 1:length(site)){

UVdf <- read.delim(paste(Project,"/",dir.UV,UV.site[site.num],".RDB",sep=""))
UVdf <- subset(UVdf,VALUE>-100000)
UVdf <- RMprep(df=UVdf,prep.type=1,date.type=4)

Pathdf <- read.delim(paste(Rlocal,"/MMSD_virus/virus2/Final compiled data/PathFIBMar192012reconciled.txt",sep=""))

Pathdf$Ebpdate <- strptime(Pathdf$SSdate,format="%m/%d/%Y %H:%M")
Pathdf$Eepdate <- strptime(Pathdf$SEdate,format="%m/%d/%Y %H:%M")

#subUV <- subset(UVdf,NAME==levels(UVdf$NAME)[2])
#subdates <- Pathdf[1:3,]

WQparms <- character()
for (i in 1:length(levels(UVdf$NAME))){
  Vname <- levels(UVdf$NAME)[i]
  letter <- ""
  l <- 0
  while(letter !=" "){
    l <- l+1
    letter <- substr(Vname,l,l)
    }
  WQparms[i] <- substr(Vname,1,l-1)
}
levels(UVdf$NAME) <- WQparms

Pathdf.stats <- Pathdf
for (parm.num in 1:length(WQparms)){
Pathdf.stats <- TSstormstats(df=UVdf,date="pdate",varname="VALUE",
                     dates=Pathdf.stats,starttime="Ebpdate",endtime="Eepdate",
                     stats.return=c("mean","median","max","min"),
                     subdfvar="NAME", subdfvalue=WQparms[parm.num],
                     out.varname=WQparms[parm.num],
                     subdatesvar="Site.Name", subdatesvalue=site[site.num])
}

dir.data <- paste(Rlocal,"/MMSD_virus/Virus2/Final compiled data/",sep="")

write.table(Pathdf.stats,paste(dir.data,"PathFIBWQ_",UV.site[site.num],".txt",sep=""),sep="\t",row.names=F)

} # End of loop for all sites

plot(Pathdf.stats$TURB_median,Pathdf.stats$lachnocen,log="xy")
plot(as.numeric(as.character(Pathdf.stats$Entero.CN.100ml)),Pathdf.stats$lachnocen)

mBH <- with(Pathdf,lm(log10(BacHumancen)~log10(TURB_median)+SC_min+SC_max+WTEMP_min+log10(Q_mean)+AntBase+sin(nmonth/11*2*pi)+cos(nmonth/11*2*pi)))
summary(mBH)

mL2 <- with(Pathdf,lm(log10(lachnocen)~log10(TURB_median)+SC_min+SC_max+WTEMP_min+log10(Q_mean)+AntBase+sin(nmonth/11*2*pi)+cos(nmonth/11*2*pi)))
summary(mL2)
plot(mL2)


mEn <- with(Pathdf,lm(as.numeric(as.character(Pathdf.stats$Entero.CN.100ml))~TURB_median+SC_min+SC_max+WTEMP_min+Q_max+Q_mean))
summary(mEn)

PathSpring <- subset(Pathdf.stats,nmonth>5 | nmonth < 9)

mHV <- with(PathSpring,lm(Human_virus~TURB_median+SC_min+SC_max+WTEMP_min+Q_max+Q_mean+sin(nmonth/11*2*pi)+cos(nmonth/11*2*pi)))
summary(mHV)
plot(residuals(mHV))

parm.num <- 2

df<-UVdf
date<-"pdate"
varname<-"VALUE"
dates<-Pathdf.stats
starttime<-"Ebpdate"
endtime<-"Eepdate"
stats.return<-c("mean","max")
subdfvar<-"NAME"
subdfvalue<-WQparms[parm.num]
out.varname<-WQparms[parm.num]
subdatesvar<-"Site.Name"
subdatesvalue<-Pathdf.stats$Site.Name[1]

# QA on results
#UV.site <- c("MenoFalls","LMDonges","Underwood","Honey","70th","16th")
#site <- c("Meno Falls","Donges Bay","Underwood","Honey","70th St.","16th St.")
#WQparms: "DO"    "Q"     "SC"    "TURB"  "WTEMP"

site.num <- 6
parm.num <- 4
  
  UVdf <- read.delim(paste(Project,"/",dir.UV,UV.site[site.num],".RDB",sep=""))
  UVdf <- subset(UVdf,VALUE>-100000)
  UVdf <- RMprep(df=UVdf,prep.type=1,date.type=4)

WQparms <- character()
for (i in 1:length(levels(UVdf$NAME))){
  Vname <- levels(UVdf$NAME)[i]
  letter <- ""
  l <- 0
  while(letter !=" "){
    l <- l+1
    letter <- substr(Vname,l,l)
  }
  WQparms[i] <- substr(Vname,1,l-1)
}
levels(UVdf$NAME) <- WQparms

  subUV <- subset(UVdf,NAME==WQparms[parm.num])

plot(subUV$pdate, subUV[,"VALUE"],pch=20,cex=0.5)

range(subUV$pdate,na.rm=T)
