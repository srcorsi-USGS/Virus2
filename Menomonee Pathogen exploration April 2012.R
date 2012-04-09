##################################################################
# Explore Menomonee virus results further
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



sum(dfPathFIBWW$Adenovirus.C.D.F>0,na.rm=T)/sum(!is.na(dfPathFIBWW$Adenovirus.C.D.F))
sum(dfPathFIBWW$Adenovirus.A >0,na.rm=T)/sum(!is.na(dfPathFIBWW$Adenovirus.A))
sum(dfPathFIBWW$Enterovirus >0,na.rm=T)/sum(!is.na(dfPathFIBWW$Enterovirus))
