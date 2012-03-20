

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

UV.site <- c("MenoFalls","LMDonges","Underwood","Honey","70th","16th")
site <- c("Meno Falls","Donges Bay","Underwood","Honey","70th St.","16th St.")

dir.UV <- "MMSD/Viruses/Menomonee Viruses/Data compilation/UnitValues/mmsd/"
i <-5
UVdf <- read.delim(paste(Project,"/",dir.UV,UV.site[i],".RDB",sep=""))
UVdf <- RMprep(df=UVdf,prep.type=1,date.type=4)

df <- read.delim(paste(Rlocal,"/MMSD_virus/virus2/virus2Prelim8.txt",sep=""))

df$Ebpdate <- strptime(df$Ebpdate,format="%Y-%m-%d %H:%M")
df$Eepdate <- strptime(df$Eepdate,format="%Y-%m-%d %H:%M")


test <- TSstormstats(UVdf,date="pdate",varnames="NAME",dates=df,
             starttime=Ebpdate,endtime=Eepdate,
             stats.return="mean",
             subdfvalue=levels(UVdf)[4],
             subdatesvar="Site.Name",
             subsdatesvalue="Meno Falls")
