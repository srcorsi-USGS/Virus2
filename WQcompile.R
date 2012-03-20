

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
site.num <-5
parm.num <- 4 #turbidity
UVdf <- read.delim(paste(Project,"/",dir.UV,UV.site[site.num],".RDB",sep=""))
UVdf <- subset(UVdf,VALUE>-100000)
UVdf <- RMprep(df=UVdf,prep.type=1,date.type=4)

Pathdf <- read.delim(paste(Rlocal,"/MMSD_virus/virus2/virus2Prelim5.txt",sep=""))

Pathdf$Ebpdate <- strptime(Pathdf$SSdate,format="%m/%d/%Y %H:%M")
Pathdf$Eepdate <- strptime(Pathdf$SEdate,format="%m/%d/%Y %H:%M")

#subUV <- subset(UVdf,NAME==levels(UVdf$NAME)[2])
#subdates <- Pathdf[1:3,]
Vname <- levels(UVdf$NAME)[parm.num]

letter <- ""
l <- 0
while(letter !=" "){
  l <- l+1
  letter <- substr(Vname,l,l)
}
Vname <- substr(Vname,1,l-1)


Pathdf.stats <- TSstormstats(df=UVdf,date="pdate",varname="VALUE",
                     dates=Pathdf,starttime="Ebpdate",endtime="Eepdate",
                     stats.return=c("mean","median","max","min"),
                     subdfvar="NAME", subdfvalue=levels(UVdf$NAME)[parm.num],
                     out.varname=Vname,
                     subdatesvar="Site.Name", subdatesvalue=site[site.num])

plot(Pathdf.stats$TURB_median,Pathdf.stats$Human_virus+0.01,log="xy")

df<-subUV
date<-"pdate"
varname<-"VALUE"
dates<-subdates
starttime<-"Ebpdate"
endtime<-"Eepdate"
stats.return<-"mean"
subdfvar<-"NAME"
subdfvalue<-levels(UVdf$NAME)[2]
subdatesvar<-"Site.Name"
subdatesvalue<-"Meno Falls"
