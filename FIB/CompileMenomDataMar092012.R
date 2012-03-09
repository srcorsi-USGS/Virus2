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

FIB <- read.delim(paste(Rlocal,"/MMSD_virus/Virus2/FIB/FIBMar092012.txt",sep=""),
                  as.is = TRUE)
FIB.orig <- FIB
FIB <- FIB[,-which(names(FIB)=="Notes")]

#Remove QA/QC samples
FIB <- FIB[-which(FIB$Event=="QA/QC"|is.na(FIB$FT.)),]

FIB$lachno <- ifelse(FIB$RERUN.Lachno.2.CN.100ml==""
                     | FIB$RERUN.Lachno.2.CN.100ml=="*",  #Combine original and rerun columns
                     FIB$Lachno.2.CN.100ml,
                     FIB$RERUN.Lachno.2.CN.100ml)
FIB$lachno <- ifelse(FIB$lachno=="0","BLD",FIB$lachno) #change 0's to BLD

#change all BLDs to detection level
dl <- 60
FIB$lachnocen<- ifelse(FIB[,"lachno"]=="BLD",dl,FIB[,"lachno"])

FIB$BacHuman <- ifelse(FIB$RERUN.BacHuman.CN.100ml==""
                     | FIB$RERUN.BacHuman.CN.100ml=="*",  #Combine original and rerun columns
                     FIB$BacHuman.CN.100ml,
                     FIB$RERUN.BacHuman.CN.100ml)

FIB$BacHuman <- ifelse(FIB$BacHuman=="0","BLD",FIB$BacHuman) #change 0's to BLD

dl <- 60
FIB$BacHumancen <- ifelse(FIB[,"BacHuman"]=="BLD",dl,FIB[,"BacHuman"])

min(as.numeric(FIB$lachno),na.rm=T)
min(as.numeric(FIB$BacHuman),na.rm=T)

min(as.numeric(FIB$lachnocen),na.rm=T)
min(as.numeric(FIB$BacHumancen),na.rm=T)

range(as.numeric(FIB$lachno),na.rm=T)
range(as.numeric(FIB$BacHuman),na.rm=T)

summary(as.numeric(FIB$lachnocen),na.rm=T)
summary(as.numeric(FIB$BacHumancen),na.rm=T)


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

IDabbrev <- substr(FIB$Site.Code,5,6)
IDchar <- nchar(FIB$Site.Code)
IDnum <- substr(FIB$Site.Code,IDchar-2,IDchar)
FIB$Study.Sample.ID <- paste(IDabbrev,IDnum)
length(FIB$Study.Sample.ID)

df.orig2 <- df
FIB.orig2 <- FIB
df <- df.orig2
FIB <- FIB.orig2


df$pdate <- strptime(df$Collection.Start.Date,format="%d-%b-%y")
df$SEpdate <- strptime(df$SEdate, format="%m/%d/%Y %H:%M")
df$Ddate <- as.Date(df$pdate)
FIB$pdate <- strptime(FIB$Sample.Date,format="%d-%b-%y")
FIB$Ddate <- as.Date(FIB$pdate)


FIBpdate <- FIB$pdate + 24*3600-1
for (i in 1:nrow(FIB)) {
  timediff <- difftime(FIBpdate[i],df$SEpdate,units="days") # difference in times from FIB and virus samples

  nearby <- which(timediff > 0) # which samples were collected before delivery and not more than two days
  if (min(timediff[nearby]) > 4){
    final.sample[i] <- NA
  } else{
  site.match <- nearby[which(df[nearby,"Abb2"]==FIB$Abb.[i])] # which samples match site from FIB sample
  final.sample[i] <- site.match[which(timediff[site.match] == min(timediff[site.match]))]
  }
}

# Mannually check and fill in NAs
which(is.na(final.sample)) #5  46  81 121 159
FIB[which(is.na(final.sample)),"pdate"] #"2009-06-15" "2009-06-15" "2009-06-15" "2009-06-15" "2009-06-15"

final.sample[5] <- 5 # HW 006
final.sample[46] <- 47 # LD 005
final.sample[81] <- 99 # MC 007
final.sample[121] <- 133 # MF 003
final.sample[159] <- 175 # MW 006


df2 <- df
FIBnames <- names(FIB)
FIBnames.orig <- FIBnames
for (i in 1:length(FIBnames)) {
  
  if(sum(FIBnames[i]==names(df))>0) FIBnames[i] <- paste(FIBnames[i],"FIB",sep="")
}
     
for (i in 1:length(FIBnames)) df2[,FIBnames[i]] <- NA

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!   NEED TO FIX THIS -- NEED TO SET NEW VARIABLES IN DF2 TO SAME CLASS AS FIB VARIABLES
for (i in 1:length(FIBnames)) class(df2[,FIBnames[i]]) <- class(FIB[,FIBnames.orig[i]])

for (i in 1:nrow(FIB)) df2[final.sample[i],FIBnames] <- FIB[i,FIBnames.orig]

i <- 0

i <- i + 1
df2[final.sample[i],FIBnames] <- FIB[i,FIBnames.orig]




value <- numeric(nrow(df))
value <- NA
data.frame(value)

df2 <- cbind(df,data.frame(
test <- merge(df,FIB,by=c("Ddate","Abb2"),all.y=F)

