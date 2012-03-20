## Script to add antecedent baseflow into the virus/FIB files

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

dir.Hysep <- paste(Rlocal,"/MMSD_virus/Virus2/Hysep",sep="")

dir.data <- paste(Rlocal,"/MMSD_virus/Virus2/Final compiled data/",sep="")


for (site.num in 1:length(site)){
  Pathdf <- read.delim(paste(dir.data,"PathFIBWQ_",UV.site[site.num],".txt",sep=""))
  Pathdf$Ebpdate <- strptime(Pathdf$SSdate,format="%m/%d/%Y %H:%M")
  Pathdf$Eepdate <- strptime(Pathdf$SEdate,format="%m/%d/%Y %H:%M")
  Pathdf$Ebddate <- as.Date(Pathdf$Ebpdate)
  
  Hdf <- read.delim(paste(dir.Hysep,"/",UV.site[site.num],".hysep.txt",sep=""))
  Hdf$pdate <- strptime(Hdf$Date,format="%m/%d/%Y")
  Hdf <- subset(Hdf, pdate > as.POSIXlt("2009-01-01"))
  Hdf$ddate <- as.Date(Hdf$pdate)
  
  AntBase <- numeric()
  for(i in 1:nrow(Pathdf)){
    if (sum(Hdf$ddate == (Pathdf$Ebddate[i]-1))>0){
      AntBase[i] <- Hdf[which(Hdf$ddate == (Pathdf$Ebddate[i]-1)),"Baseflow..cfs."]
    }
  }
  Pathdf$AntBase <- AntBase
  
  write.table(Pathdf,paste(dir.data,"Path_FIB_WQ_Hy_",UV.site[site.num],".txt",sep=""),sep="\t",row.names=F)
}
  
    
  
