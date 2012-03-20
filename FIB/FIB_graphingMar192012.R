
###################
#Set R directory
SNsys <- system("wmic bios get serialnumber", intern = TRUE, show.output.on.console = FALSE)
SN <- gsub("\\s","",SNsys)[2]
if(SN == "R97R0NW") {Rlocal <- "D:/srcldata/R" #Laptop
       }else {Rlocal <- "//igsarmewwssrc/SRCdata/R"} #Network
       
library(latticeExtra)
#library(quantreg)


FIB <- read.delim(paste(Rlocal,"/MMSD_virus/virus2/FIB/PathFIBMar192012reconciled.txt",sep=""))
names.to.change <- match(c("Enterococci.CFU.100ml","E..coli.CFU.100ml","BacHumancen","lachnocen",
                     "BacSp.CN.100ml", "Entero.CN.100ml"),names(FIB))
names(FIB)[names.to.change] <- c("Enterococci.CFU","E.coli.CFU","Human.Bacteroides","Lachno2",
                                 "Bacteroides.Sp","Enterococci.CN")

#NEED TO MODIFY THIS FOR FINAL DATA
# FIBww <- read.delim(paste(Rlocal,"/MMSD_virus/FIB_influent.txt",sep=""))
# names(FIBww) <- c("FT","Sample.Date","Site.Code","Entero.Culture","EC.Culture",
# "BacHuman","BacSp","Ratio.Human.Sp","Entero","Lachno","Notes")
# FIBww$Abb <- factor("South Shore",c("Jones Island","South Shore"))
# FIBww$Abb[grep("J",FIBww$Site.Code)] <- "Jones Island"
  
# ## Convert FIB data to character strings from factors ##
# for (i in 6:13) FIB[,i] <- as.character(FIB[,i])
# (min(as.numeric(FIB[which(FIB[,i]!="BLD"),i]))/10)

## NOT USED Set samples below Detection levels to 10% of minimum actual value ##
#for (i in 6:12) FIB[,i] <- ifelse(FIB[,i]=="BLD",
#                          (min(as.numeric(FIB[which(FIB[,i]!="BLD"),i]))/10),FIB[,i])

# ## Set left censored values to 1 and convert from alpha to numeric##
# for (i in 6:12) FIB[,i] <- as.numeric(ifelse(FIB[,i]=="BLD",
#                           1,FIB[,i]))

# ## Set left censored values to 100 for BacHuman and Lachno 2 ##
# i=8; FIB[,i] <- ifelse(FIB[,i]==1,100,FIB[,i])
# i=12; FIB[,i] <- ifelse(FIB[,i]==1,100,FIB[,i])

# FIB$Ratio.Human.Sp <- FIB$BacHuman/FIB$BacSp
# FIB$Ratio.Human.Sp <- ifelse(FIB$BacHuman==1 | FIB$BacSp==1,NA,FIB$Ratio.Human.Sp)
levels(FIB$USGS.Event.Type.2)
levels(FIB$Abb2)
FIB$Abb2 <- factor(FIB$Abb2,c("LD","MF","UW","HW","MW","MC"))
FIB$USGS.Event.Type.2 <- factor(FIB$USGS.Event.Type.2,c("Baseflow","Rainfall","CSO","Mix","Snowmelt"))

## Boxplot of Human Bacteroides  ##
options(scipen=10)
bx <- boxplot(Human.Bacteroides~USGS.Event.Type.2+Abb2,data=FIB,xaxt="n", main = "",
        ylab="Human Bacteroides (CN/100ml)",log="y",ylim=c(10,1000000))
#mtext(text=paste(levels(FIB$Abb),collapse="               "),side=1,line=2)
axis(side=1,labels=rep(c("Base","Rainfall","CSO","Mix","Snowmelt"),times=6),at=c(1:30),cex.axis=0.8,las=3)
#mtext(text=rep(c("Base","Event","CSO"),times=6),at=c(1:18),side=1,line=1)
mtext(text=levels(FIB$Abb2),at=c(3,8,13,18,23,28),side=3,line=1.)
mtext(text="Human Bacteroides at Menomonee R. Sites",side=3,line=3.,cex=1.5)
abline(v=((1:5)*5+0.5),lty=3)
abline(h=60,lty=2)
mtext("Detection level",line=-2.7,side=1)
  mtext(cex=0.7,paste(bx$n, sep = ""),
  at = seq_along(bx$n), line = 0.2, side = 3)

## Boxplot of Human Bacteroides--all sites  ##
options(scipen=10)
bx <- boxplot(Human.Bacteroides~USGS.Event.Type.2,data=FIB,xaxt="n", main = "",
              ylab="Human Bacteroides (CN/100ml)",log="y",ylim=c(10,1000000))
#mtext(text=paste(levels(FIB$Abb),collapse="               "),side=1,line=2)
axis(side=1,labels=levels(FIB$USGS.Event.Type.2),at=c(1:5),cex.axis=0.8,las=3)
#mtext(text=rep(c("Base","Event","CSO"),times=6),at=c(1:18),side=1,line=1)
mtext(text="Human Bacteroides at Menomonee R. Sites",side=3,line=3.,cex=1.5)
abline(h=60,lty=2)
mtext("Detection level",line=-2.7,side=1)
mtext(cex=0.7,paste(bx$n, sep = ""),
      at = seq_along(bx$n), line = 0.2, side = 3)

## Boxplot of Lachno 2  ##
options(scipen=10)
bx <- boxplot(Lachno~Event+Abb,data=FIB,xaxt="n", main = "",
        ylab="Lachnospiraceae (CN/100ml)",log="y",ylim=c(10,1000000))
#mtext(text=paste(levels(FIB$Abb),collapse="               "),side=1,line=2)
axis(side=1,labels=rep(c("Base","Event","CSO"),times=6),at=c(1:18),cex.axis=0.8,las=3)
#mtext(text=rep(c("Base","Event","CSO"),times=6),at=c(1:18),side=1,line=1)
mtext(text=levels(FIB$Abb),at=c(2,5,8,11,14,17),side=3,line=1.)
mtext(text="Lachnospiraceae at Menomonee R. Sites",side=3,line=3.,cex=1.5)
abline(v=((1:5)*3+0.5),lty=3)
abline(h=100,lty=2)
mtext("Detection level",line=-2.7,side=1)
  mtext(cex=0.7,paste(bx$n, sep = ""),
  at = seq_along(bx$n), line = 0.2, side = 3)


## Boxplot of Human Bacteroides including all sites  ##
options(scipen=10)
par(mar=c(5,5,4,2)+0.01)
bx <- boxplot(BacHuman~Event,data=FIB,xaxt="n", yaxt="n",main = "",
        log="y",ylim=c(10,1000000),las=2)
#mtext(text=paste(levels(FIB$Abb),collapse="               "),side=1,line=2)
axis(side=1,labels=c("Base","Event","CSO"),at=c(1:3),cex.axis=0.8,las=2)
axis(side=2,cex.axis=0.8,las=2)
#mtext(text=rep(c("Base","Event","CSO"),times=6),at=c(1:18),side=1,line=1)
mtext(text="Human Bacteroides (CN/100ml)",side=2,line=3.5,cex=1.0)
mtext(text="Human Bacteroides",side=3,line=2.7,cex=1.5)
#abline(v=((1:5)*3+0.5),lty=3)
abline(h=100,lty=2)
mtext("Detection level",line=-2.7,side=1)
  mtext(cex=0.7,paste(bx$n, sep = ""),
  at = seq_along(bx$n), line = 0.2, side = 3)


## Boxplot of Lachno 2 including all sites  ##
options(scipen=10)
par(mar=c(5,5,4,2)+0.01)
bx <- boxplot(Lachno~Event,data=FIB,xaxt="n", yaxt="n",main = "",
        log="y",ylim=c(10,1000000),las=2)
#mtext(text=paste(levels(FIB$Abb),collapse="               "),side=1,line=2)
axis(side=1,labels=c("Base","Event","CSO"),at=c(1:3),cex.axis=0.8,las=2)
axis(side=2,cex.axis=0.8,las=2)
#mtext(text=rep(c("Base","Event","CSO"),times=6),at=c(1:18),side=1,line=1)
mtext(text="Lachnospiraceae (CN/100ml)",side=2,line=3.5,cex=1.0)
mtext(text="Lachnospiraceae",side=3,line=2.7,cex=1.5)
#abline(v=((1:5)*3+0.5),lty=3)
abline(h=100,lty=2)
mtext("Detection level",line=-2.7,side=1)
  mtext(cex=0.7,paste(bx$n, sep = ""),
  at = seq_along(bx$n), line = 0.2, side = 3)


## Boxplot of Entero PCR including all sites  ##
options(scipen=10)
par(mar=c(5,5,4,2)+0.01)
bx <- boxplot(Entero~Event,data=FIB,xaxt="n", yaxt="n",main = "",
        log="y",ylim=c(10,1000000),las=2)
#mtext(text=paste(levels(FIB$Abb),collapse="               "),side=1,line=2)
axis(side=1,labels=c("Base","Event","CSO"),at=c(1:3),cex.axis=0.8,las=2)
axis(side=2,cex.axis=0.8,las=2)
#mtext(text=rep(c("Base","Event","CSO"),times=6),at=c(1:18),side=1,line=1)
mtext(text="Enterococci PCR (CN/100ml)",side=2,line=3.5,cex=1.0)
mtext(text="Enterococci PCR",side=3,line=2.7,cex=1.5)
#abline(v=((1:5)*3+0.5),lty=3)
abline(h=100,lty=2)
# mtext("Detection level",line=-2.7,side=1)
   mtext(cex=0.7,paste(bx$n, sep = ""),
   at = seq_along(bx$n), line = 0.2, side = 3)




## Boxplot of Entero Culture including all sites  ##
options(scipen=10)
par(mar=c(5,5,4,2)+0.01)
bx <- boxplot(Entero.Culture~Event,data=FIB,xaxt="n", yaxt="n",main = "",
        log="y",ylim=c(10,1000000),las=2)
#mtext(text=paste(levels(FIB$Abb),collapse="               "),side=1,line=2)
axis(side=1,labels=c("Base","Event","CSO"),at=c(1:3),cex.axis=0.8,las=2)
axis(side=2,cex.axis=0.8,las=2)
#mtext(text=rep(c("Base","Event","CSO"),times=6),at=c(1:18),side=1,line=1)
mtext(text="Enterococci Culture (CN/100ml)",side=2,line=3.5,cex=1.0)
mtext(text="Enterococci Culture",side=3,line=2.7,cex=1.5)
#abline(v=((1:5)*3+0.5),lty=3)
abline(h=100,lty=2)
# mtext("Detection level",line=-2.7,side=1)
   mtext(cex=0.7,paste(bx$n, sep = ""),
   at = seq_along(bx$n), line = 0.2, side = 3)

############## WASTEWATER INFLUENT FIB RESULTS ###################
###### Set left censored data from wastewater to 100 ###########
## Convert FIB data to character strings from factors ##
for (i in 4:8) FIB[,i] <- as.character(FIB[,i])
for (i in 4:8) FIB[,i] <- as.numeric(ifelse(FIB[,i]=="BLD",
                          1,FIB[,i]))
## Boxplot of Lachno 2 for wastewater sites  ##
options(scipen=10)
par(mar=c(5,6,4,2)+0.01)
bx <- boxplot(Lachno~Abb,data=FIBww,xaxt="n", yaxt="n",main = "",
        log="y",las=2)
#mtext(text=paste(levels(FIB$Abb),collapse="               "),side=1,line=2)
axis(side=1,labels=levels(FIBww$Abb),at=c(1:2),cex.axis=0.8,las=2)
axis(side=2,cex.axis=0.8,las=2)
#mtext(text=rep(c("Base","Event","CSO"),times=6),at=c(1:18),side=1,line=1)
mtext(text="Lachnospiraceae (CN/100ml)",side=2,line=4.5,cex=1.0)
mtext(text="Lachnospiraceae",side=3,line=2.7,cex=1.5)
#abline(v=((1:5)*3+0.5),lty=3)
#abline(h=100,lty=2)
#mtext("Detection level",line=-2.7,side=1)
  mtext(cex=0.7,paste(bx$n, sep = ""),
  at = seq_along(bx$n), line = 0.2, side = 3)

## Boxplot of Human Bacteroides for wastewater sites  ##
options(scipen=10)
par(mar=c(5,6,4,2)+0.01)
bx <- boxplot(BacHuman~Abb,data=FIBww,xaxt="n", yaxt="n",main = "",
        log="y",las=2)
#mtext(text=paste(levels(FIB$Abb),collapse="               "),side=1,line=2)
axis(side=1,labels=levels(FIBww$Abb),at=c(1:2),cex.axis=0.8,las=2)
axis(side=2,cex.axis=0.8,las=2)
#mtext(text=rep(c("Base","Event","CSO"),times=6),at=c(1:18),side=1,line=1)
mtext(text="Human Bacteroides (CN/100ml)",side=2,line=4.5,cex=1.0)
mtext(text="Human Bacteroides",side=3,line=2.7,cex=1.5)
#abline(v=((1:5)*3+0.5),lty=3)
#abline(h=100,lty=2)
#mtext("Detection level",line=-2.7,side=1)
  mtext(cex=0.7,paste(bx$n, sep = ""),
  at = seq_along(bx$n), line = 0.2, side = 3)


## Boxplot of PCR Entero for wastewater sites  ##
options(scipen=10)
par(mar=c(5,6,4,2)+0.01)
bx <- boxplot(Entero~Abb,data=FIBww,xaxt="n", yaxt="n",main = "",
        log="y",las=2)
#mtext(text=paste(levels(FIB$Abb),collapse="               "),side=1,line=2)
axis(side=1,labels=levels(FIBww$Abb),at=c(1:2),cex.axis=0.8,las=2)
axis(side=2,cex.axis=0.8,las=2)
#mtext(text=rep(c("Base","Event","CSO"),times=6),at=c(1:18),side=1,line=1)
mtext(text="Enterococci PCR (CN/100ml)",side=2,line=4.5,cex=1.0)
mtext(text="Enterococci PCR",side=3,line=2.7,cex=1.5)
#abline(v=((1:5)*3+0.5),lty=3)
#abline(h=100,lty=2)
#mtext("Detection level",line=-2.7,side=1)
  mtext(cex=0.7,paste(bx$n, sep = ""),
  at = seq_along(bx$n), line = 0.2, side = 3)


## Boxplot of BacSp for wastewater sites  ##
options(scipen=10)
par(mar=c(5,6,4,2)+0.01)
bx <- boxplot(BacSp~Abb,data=FIBww,xaxt="n", yaxt="n",main = "",
        log="y",las=2)
#mtext(text=paste(levels(FIB$Abb),collapse="               "),side=1,line=2)
axis(side=1,labels=levels(FIBww$Abb),at=c(1:2),cex.axis=0.8,las=2)
axis(side=2,cex.axis=0.8,las=2)
#mtext(text=rep(c("Base","Event","CSO"),times=6),at=c(1:18),side=1,line=1)
mtext(text="Bacteroides Sp (CN/100ml)",side=2,line=4.5,cex=1.0)
mtext(text="Bacteroides Sp",side=3,line=2.7,cex=1.5)
#abline(v=((1:5)*3+0.5),lty=3)
#abline(h=100,lty=2)
#mtext("Detection level",line=-2.7,side=1)
  mtext(cex=0.7,paste(bx$n, sep = ""),
  at = seq_along(bx$n), line = 0.2, side = 3)

