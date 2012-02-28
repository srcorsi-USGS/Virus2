#####################################################################################
# Script to compile ancillary data into Memonomee watershed 2009-2011 pathogen data #
#####################################################################################
  

###################
#Set R directory
SNsys <- system("wmic bios get serialnumber", intern = TRUE, show.output.on.console = FALSE)
SN <- gsub("\\s","",SNsys)[2]
if(SN == "R97R0NW") {Rlocal <- "D:/srcldata/R" #Laptop
}else {Rlocal <- "//igsarmewwssrc/SRCdata/R"} #Network

setwd(paste(Rlocal,"/MMSD_virus/Virus2",sep=""))

###################
#Read data file and functions
source(paste(Rlocal,"/Scripts/fxn_categorize.R",sep=""))
source(paste(Rlocal,"/Scripts/fxns_rmNA_rmZeros.R",sep=""))
source(paste(Rlocal,"/Scripts/fxn_drop.unused.factors.R",sep=""))

#df <- read.delim("D:/srcldata/MMSD/Marshfield results/Jan 2012/menopathTEST.txt")
df <- read.delim("virus2Prelim5.txt")
df.orig <- df # archive original copy

FIB <- read.delim(paste(Rlocal,"/MMSD_virus/Virus2/FIB/FIBFeb222012.txt",sep=""),
                  as.is = TRUE)
FIB.orig <- FIB
#FIB <- FIB[,-which(names(FIB)=="Notes")]

#Remove QA/QC samples
FIB <- FIB[-which(FIB$Event=="QA/QC"|is.na(FIB$FT.)),]

FIB$lachno <- ifelse(FIB$RERUN.Lachno.2.CN.100ml=="",  #Combine original and rerun columns
                     FIB$Lachno.2.CN.100ml,
                     FIB$RERUN.Lachno.2.CN.100ml)
FIB$lachno <- ifelse(FIB$lachno=="0","BLD",FIB$lachno) #change 0's to BLD

#change all BLDs to detection level
dl <- 60
FIB$lachnocen<- ifelse(FIB[,"lachno"]=="BLD",dl,FIB[,"lachno"])

FIB$BacHuman <- ifelse(FIB$RERUN.Lachno.2.CN.100ml=="",  #Combine original and rerun columns
                     FIB$BacHuman.CN.100ml,
                     FIB$RERUN.BacHuman.CN.100ml)

FIB$BacHuman <- ifelse(FIB$BacHuman=="0","BLD",FIB$BacHuman) #change 0's to BLD

dl <- 60
FIB$BacHumancen <- ifelse(FIB[,"BacHuman"]=="BLD",dl,FIB[,"BacHuman"])

min(as.numeric(FIB$lachno),na.rm=T)
min(as.numeric(FIB$BacHuman),na.rm=T)

range(as.numeric(FIB$lachno),na.rm=T)
range(as.numeric(FIB$BacHuman),na.rm=T)

summary(as.numeric(FIB$lachnocen),na.rm=T)
summary(as.numeric(FIB$BacHumancen),na.rm=T)

df$pdate <- strptime(df$Collection.Start.Date,format="%d-%b-%y")
df$Ddate <- as.Date(df$pdate)
FIB$pdate <- strptime(FIB$Sample.Date,format="%d-%b-%y")
FIB$Ddate <- as.Date(FIB$pdate)

#Add factor variables for site abbreviation and site name to match those in pathogen data frame
Sites <- c("16th St.","70th St.","Honey","Underwood","Donges Bay","Meno Falls" )
abbrev <- c("MC","MW","HW","UW","LD","MF")
FIB$Abb2 <- factor(FIB$Abb.,levels=abbrev)
FIB$Site.Name <- FIB$Abb2
levels(FIB$Site.Name) <- Sites
df$Abb2 <- df$Site.Name
df$Abb2 <- factor(df$Abb2,levels=Sites)
levels(df$Abb2) <- abbrev
levels(df$Abb2)


test <- merge(df,FIB,by=c("Ddate","Abb2"),all.y=F)

23*.13

