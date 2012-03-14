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
df <- read.delim(paste(Rlocal,"/MMSD_virus/Virus2/FIB/PathFIBMar122012.txt",sep=""))
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


# Determine which sample data from FIB is closest, but the same day or after the
# ending sample date.
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
final.sample[81] <- 88 # MC 007
final.sample[121] <- 133 # MF 003
final.sample[159] <- 175 # MW 006



# Resolve duplicates in final.sample vector
which(duplicated (final.sample))
final.sample[16:17] <- c(19,20) # from 20
final.sample[21:22] <- c(24,25) # from 25
final.sample[26:28] <- c(29,NA,NA) # NO GOOD MATCH FOR 27 AND 28 AT HONEY, BUT SAMPLES AT OTHERS) # from 29
final.sample[60:61] <- c(65,66) # from 66
final.sample[94:95] <- c(105,106) #from 106
final.sample[98:99] <- c(109,NA) #did not send in 99 MC037 for viruses
final.sample[106:107] <-c(116,117) # from 116
final.sample[125:126] <- c(137,138) # MF010 for 138, but it does not match well. end date virus is after received FIB
final.sample[137:138] <- c(151,152) # from 152
final.sample[171:172] <- c(190,191) # from 191
final.sample[176:177] <- c(195,196) # from 196
final.sample[188:189] <- c(207,208) # from 208

FIB.no.Virus <- c(27,28,99)



Vtimes <- df[final.sample,"SEpdate"]
dftimes <- data.frame(Vtimes,FIBpdate)
t2 <- round(with(dftimes,(difftime(FIBpdate,Vtimes,units="days"))),1)


df2 <- df
# FIBnames <- names(FIB)
# FIBnames.orig <- FIBnames
# for (i in 1:length(FIBnames)) {
#   
#   if(sum(FIBnames[i]==names(df))>0) FIBnames[i] <- paste(FIBnames[i],"FIB",sep="")
# }
#      
# for (i in 1:length(FIBnames)) df2[,FIBnames[i]] <- NA

FIBdf1 <- cbind(df[final.sample,],FIB)
FIBdf1 <- data.frame(FIBdf1)
tdiff <- difftime(FIBdf1$pdate.1,FIBdf1$SEpdate,units="days")

min(tdiff,na.rm=T)
max(tdiff,na.rm=T)




FIB.rows <- numeric()
FIB.NA <- numeric()
j <- 1
for (i in 1:nrow(df2)){
  if(length(which(final.sample==i) > 0) > 0){ # if there is an i that matches something in final.sample
    FIB.rows[i] <- which(final.sample==i)
    j <- j + 1
  } else{
    FIB.rows[i] <- 1
    FIB.NA <- c(FIB.NA,i)
  }
}



FIBdf <- FIB[FIB.rows,]
  
  
  
FIBdf <- FIB[FIB.rows,]
FIBdf[FIB.NA,] <- NA

df3 <- cbind(df,FIBdf)

#resolve duplicate column header names making them unique by adding ".1" on the 
#end of the second of the duplicates
df3 <- data.frame(df3) # could also use names(df3) <- make.names(names(df3),unique=TRUE)

NA1 <- (nrow(df3)+1)
NA2 <- (nrow(df3)+length(FIB.no.Virus))
df3[NA1:NA2,] <- NA
### START HERE--NEED TO DEAL WITH NA COLUMNS THAT DO NOT HAVE VIRUSES, BUT DO HAVE FIB
#ADD IN COLUMNS FROM fib NEED TO DETERMINE WHICH COLUMNS FIRST
bFIB <- ncol(df2)+1
eFIB <- ncol(df3)
df3[NA1:NA2,bFIB:eFIB]  <- FIB[FIB.no.Virus,]


tdiff <- difftime(df3$pdate.1,df3$SEpdate,units="days")

max(tdiff,na.rm=T)
min(tdiff,na.rm=T)

test <- cbind(df[final.sample,],FIB)
test <- data.frame(test)
difftime(test$pdate.1,test$SEpdate,units="days")

write.table(df3,paste(Rlocal,"/MMSD_virus/Virus2/FIB/PathFIB.txt",sep=""),sep="\t",row.names=F)

# Check on IDs
IDs <- data.frame(VirusID=df3$Study.Sample.ID,FIBID=df3$Site.Code,FIBID2=df3$Study.Sample.ID.1)
diffIDs <- IDs[which(as.character(IDs[,1])!=as.character(IDs[,3])),]
row.names(diffIDs)


df.premerge <- read.delim(paste(Rlocal,"/MMSD_virus/Virus2/virus2Prelim5.txt",sep=""))
IDmergedf <- merge(x=df.premerge,y=FIB,by.x="Study.Sample.ID",by.y="Study.Sample.ID",all=T)
IDs2 <- data.frame(VirusID=IDmergedf$Study.Sample.ID,FIBID=IDmergedf$Site.Code,FIBID2=IDmergedf$Study.Sample.ID)
diffIDs2 <- IDs2[which(as.character(IDs2[,1])!=as.character(IDs2[,3])),]
row.names(diffIDs2)

IDmergedf[which(duplicated(IDmergedf$Study.Sample.ID)),]
str(IDmergedf)
IDmergedf$SEpdate <- strptime(IDmergedf$SEdate, format="%m/%d/%Y %H:%M")

tdiff2 <- difftime((IDmergedf$pdate.y+24*60*60),IDmergedf$SEpdate,units="days")

large.tdiff2 <- IDmergedf[which(tdiff2>4),c("Study.Sample.ID","Site.Code")]
