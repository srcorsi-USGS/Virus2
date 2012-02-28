
# Process menomonee virus file
# 
# Include bar chart of occurrence for events vs baseflow

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
df <- read.delim("virus2Prelim4.txt")
df.orig <- df # archive original copy

#write.table(names(meno), file="menoNames.txt",sep="\t",row.names=F)
dfNames <- read.delim("menoNames.txt")

for(i in 1:ncol(dfNames)) dfNames[,i] <- as.character(dfNames[,i])

## define columns in names file ##
variables <- "Variables"  # define colum in names file with variable names
groups <- "Groups"  # define column in names file with group definitions
VarType <- "VarType"

## define some critical columns in df ##
Site.Name <- "Site.Name"
Sample.Type  <- "Sample.Type"
Event.Type  <- "Event.Type"



# Define groups of pathogens (response variables)
pathType <- unique(dfNames[dfNames[,VarType]=="Response",groups])

#levels(dfNames$Groups) <- unique(dfNames$Groups)
#sum groups of viruses


for (i in 1:length(pathType)){
  pathRows <- which(dfNames[,groups]==pathType[i])
  pathNames <- dfNames[pathRows,variables]
  pathRowsdf <- match(pathNames,names(df))
  if(length(pathRows)!=1){ df[,pathType[i]] <- apply(df[,pathRowsdf],1,sum)
  }else df[,pathType[i]] <- df[,pathNames]
}

#sum all human pathogens
df$Human <- df$Human_virus + df$Bacteria

#sum of all pathogens
df$All_pathogens <- df$Human + df$Bovine_virus

pathType <- c(pathType,"Human","All_pathogens")


##########################################################
## Create barchart of occurrence by pathogen class
##########################################################
#

#Compute % occurrence for each of the pathogen groups
pathCols <- numeric()
for (i in 1:length(pathType)) pathCols[i] <- which(names(df)==pathType[i])
      
bcol <- min(pathCols)
ecol <- max(pathCols)
occurrence <- function(x) {sum(x != 0)/sum(!is.na(x))}
occurrence <- function(x) {sum(x != 0,na.rm=T)/sum(!is.na(x),na.rm=T)}
path_occur <- apply(df[,bcol:ecol],2,occurrence)
samplecount <- function(x) {sum(!is.na(x),na.rm=T)}
nvals <- apply(df[,bcol:ecol],2,samplecount)


# Create bar chart with occurrence for the three pathogen classes
#
# Turn off scientific notation on axes:
options(scipen=10)

#########
# Create bar chart with percent occurrence of pathogens by pathogen class
par(mar=c(8,4,4,2))

mp <- barplot(path_occur, beside = TRUE,
              col = c("lightblue"),      
              main = "Menomonee River Watershed Pathogen Occurrence", 
              font.main = 4,
              sub = "", col.sub = "black",
              cex.names = 1., ylim=c(0,1),
              ylab = "Pathogen Occurrence",
              xlab = "",
              las=2)


mtext(paste(round(path_occur,2), 
            sep = ""), 
      at = mp, line = -1, side = 3,
      cex=0.8)
mtext(paste(" (n=",nvals,")",sep = ""), 
      at = mp, line = -2, side = 3,
      cex=0.8)
mtext("Pathogen Class",side=1,line=6.5,cex=1.2)
box()



##########################################################
## Create barchart of occurrence by site for all pathogens
##########################################################
#

#Compute % occurrence for each of the pathogen groups for subsets by site


# Reorder factors in Site.Name column in downstream order
df[,Site.Name] <- factor(df[,Site.Name],levels=c("Meno Falls","Donges Bay","Underwood","Honey","70th St.","16th St."))

#pdf("pathOccurBySite.pdf")

sites <- as.character(levels(df[,Site.Name]))
numsites <- length(sites)
par(mfrow=c(3,2))
for (j in 1:numsites){
subdf <- subset(df,Site.Name==sites[j])
pathCols <- numeric()
for (i in 1:length(pathType)) pathCols[i] <- which(names(subdf)==pathType[i])

bcol <- min(pathCols)
ecol <- max(pathCols)
occurrence <- function(x) {sum(x != 0,na.rm=T)/sum(!is.na(x),na.rm=T)}
path_occur <- apply(subdf[,bcol:ecol],2,occurrence)
samplecount <- function(x) {sum(!is.na(x),na.rm=T)}
nvals <- apply(subdf[,bcol:ecol],2,samplecount)

# Create bar chart with occurrence for the three pathogen classes
#
# Turn off scientific notation on axes:
options(scipen=10)

#########
# Create bar chart with percent occurrence of pathogens by pathogen class
par(mar=c(5,4,4,2))

mp <- barplot(path_occur, beside = TRUE,
              col = c("lightblue"),      
              main = paste(sites[j],"Pathogen Occurrence"),
              font.main = 4,
              sub = "", col.sub = "black",
              cex.names = 1., ylim=c(0,1),
              ylab = "% Occurrence",
              xlab = "",
              las=2)


mtext(paste(round(path_occur,2), 
            sep = ""), 
      at = mp, line = -1, side = 3,
      cex=0.8)
mtext(paste(" (n=",nvals,")",sep = ""), 
      at = mp, line = -2, side = 3,
      cex=0.8)

box()
}
mtext("Pathogen Class",side=1,line=6.5,cex=1.2)
#dev.off()


##########################################################
## Create barchart of occurrence by site for all pathogens
##########################################################
#

#Compute % occurrence for each of the pathogen groups for subsets by site
levels(df[,Sample.Type])


# Reorder factors in Site.Name column in downstream order
df[,Site.Name] <- factor(df[,Site.Name],levels=c("Meno Falls","Donges Bay","Underwood","Honey","70th St.","16th St."))

#pdf("pathOccurBySite.pdf")

sites <- as.character(levels(df[,Sample.Type]))


numsites <- length(sites)
par(mfrow=c(2,1))
for (j in 1:numsites){
  subdf <- subset(df,Sample.Type==sites[j])
  pathCols <- numeric()
  for (i in 1:length(pathType)) pathCols[i] <- which(names(subdf)==pathType[i])
  
  bcol <- min(pathCols)
  ecol <- max(pathCols)
  occurrence <- function(x) {sum(x != 0,na.rm=T)/sum(!is.na(x),na.rm=T)}
  path_occur <- apply(subdf[,bcol:ecol],2,occurrence)
  samplecount <- function(x) {sum(!is.na(x),na.rm=T)}
  nvals <- apply(subdf[,bcol:ecol],2,samplecount)
  
  # Create bar chart with occurrence for the three pathogen classes
  #
  # Turn off scientific notation on axes:
  options(scipen=10)
  
  #########
  # Create bar chart with percent occurrence of pathogens by pathogen class
  par(mar=c(5,4,4,2))
  
  mp <- barplot(path_occur, beside = TRUE,
                col = c("lightblue"),      
                main = paste(sites[j],"Pathogen Occurrence"),
                font.main = 4,
                sub = "", col.sub = "black",
                cex.names = 1., ylim=c(0,1),
                ylab = "% Occurrence",
                xlab = "",
                las=2)
  
  
  mtext(paste(round(path_occur,2), 
              sep = ""), 
        at = mp, line = -1, side = 3,
        cex=0.8)
  mtext(paste(" (n=",nvals,")",sep = ""), 
        at = mp, line = -2, side = 3,
        cex=0.8)
  
  box()
}



##########################################################
## Create barchart of occurrence by site for all pathogens
##########################################################
#

#Compute % occurrence for each of the pathogen groups for subsets by site
levels(df[,Sample.Type])


# Reorder factors in Site.Name column in downstream order
df[,Site.Name] <- factor(df[,Site.Name],levels=c("Meno Falls","Donges Bay","Underwood","Honey","70th St.","16th St."))

#pdf("pathOccurBySite.pdf")

sites <- as.character(levels(df[,Sample.Type]))


numsites <- length(sites)
par(mfrow=c(2,1))
for (j in 1:numsites){
  subdf <- subset(df,Sample.Type==sites[j])
  pathCols <- numeric()
  for (i in 1:length(pathType)) pathCols[i] <- which(names(subdf)==pathType[i])
  
  bcol <- min(pathCols)
  ecol <- max(pathCols)
  occurrence <- function(x) {sum(x != 0,na.rm=T)/sum(!is.na(x),na.rm=T)}
  path_occur <- apply(subdf[,bcol:ecol],2,occurrence)
  samplecount <- function(x) {sum(!is.na(x),na.rm=T)}
  nvals <- apply(subdf[,bcol:ecol],2,samplecount)
  
  # Create bar chart with occurrence for the three pathogen classes
  #
  # Turn off scientific notation on axes:
  options(scipen=10)
  
  #########
  # Create bar chart with percent occurrence of pathogens by pathogen class
  par(mar=c(5,4,4,2))
  y <- subdf$Bovine_virus
  y <- y[!is.na(y)]
  mp <- stripchart(y,
                log="y",ylim=c(0.01,10000))
  ,
                col = c("lightblue"),      
                main = paste(sites[j],"Pathogen Occurrence"),
                font.main = 4,
                sub = "", col.sub = "black",
                cex.names = 1.,
                ylab = "% Occurrence",
                xlab = "",
                las=2)
  
  
  mtext(paste(round(path_occur,2), 
              sep = ""), 
        at = mp, line = -1, side = 3,
        cex=0.8)
  mtext(paste(" (n=",nvals,")",sep = ""), 
        at = mp, line = -2, side = 3,
        cex=0.8)
  
  box()
}

dev.off()


###############################################################################
## 


#######################################################################
################   Seasonal for individual viruses ####################
#######################################################################

## determine occurrence for individual human viruses  and bacteria ##

# define POSIX date and month for each sample
df$pdate <- strptime(df$Collection.Start.Date,format="%d-%b-%y")
df$nmonth <- unclass(df$pdate)$mon

# Define individual pathogens to graph
pathNames <- dfNames[which(dfNames$Groups=="Human_virus" | dfNames$Groups=="Bacteria"),Variables]
pathColsdf <- match(pathNames,names(df))
sumdf <- apply(df[,pathColsdf],2,sum,na.rm=T)

# Only choose the pathogens that actually have occurrence
Present <- names(which(sumdf>0))
Presentcols <- match(Present,names(df))

#Compute % occurrence for each of the individual human and bacterial pathogens 
pathCols <- numeric()
for (i in 1:length(pathType)) pathCols[i] <- which(names(df)==pathType[i])

occurrence <- function(x) {sum(x != 0,na.rm=T)/sum(!is.na(x),na.rm=T)}
path_occur <- apply(df[,Presentcols],2,occurrence)
samplecount <- function(x) {sum(!is.na(x),na.rm=T)}
nvals <- apply(df[,Presentcols],2,samplecount)


for (i in 1:length(Presentcols)){
  varname <- paste(names(df)[Presentcols[i]],"occur",sep="_")
  df[,varname] <- ifelse(df[,Presentcols[i]]>0,1,0)
}

Cmonth <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

columns <- ceiling(length(Present)/3)
rows <- ceiling(length(Present)/columns)
par(mfrow=c(rows,columns))     #layout for 2 rows and 2 columns

for (i in 1:length(Presentcols)){
  varname <- paste(Present[i],"occur",sep="_")
  if(sum(df[,varname],na.rm=T)>=1) {
    V_month <- tapply(df[,varname],df$nmonth,mean,na.rm=T) # % occurrence
    V_month.df <- data.frame(V_month)
    V_month.df$month <- row.names(V_month.df)
    V_month <- rbind(rbind(rbind(rbind(V_month.df[1:6,],c(NA,6)),V_month.df[7:9,]),c(NA,10)),V_month.df[10,])
    V_month.final <- as.data.frame(t(V_month[,1]))
    names(V_month.final) <- Cmonth
    barplot(as.matrix(V_month.final), main="",
            xlab="", col=c("darkblue"),ylim=c(0,1),
            beside=TRUE,xaxt="n",yaxt="n")
    axis(side=2,at=c(0,0.5,1),labels=c("0","50%","100%"),cex.axis=1.5,las=2)
    axis(side=1,at=(seq(1:12)*2-0.5),labels=Cmonth,cex.axis=1.75,las=3)
    abline(h=0,v=0)
    mtext(names(df)[Presentcols[i]],side=3,line=0.5,cex=1.2)
  }
  
}
par(mfrow=c(1,1))     #change the graph window back to one figure
mtext("Human Viruses and Bacteria",side=3,line=2.9,cex=1.5)


#######################################################################
################   Seasonal for individual viruses ####################
#######################################################################

## determine occurrence for individual bovine viruses  ##

# define POSIX date and month for each sample
df$pdate <- strptime(df$Collection.Start.Date,format="%d-%b-%y")
df$nmonth <- unclass(df$pdate)$mon

# Define individual pathogens to graph
pathNames <- dfNames[which(dfNames$Groups=="Bovine_virus"),Variables]
pathColsdf <- match(pathNames,names(df))
sumdf <- apply(df[,pathColsdf],2,sum,na.rm=T)

# Only choose the pathogens that actually have occurrence
Present <- names(which(sumdf>0))
Presentcols <- match(Present,names(df))

#Compute % occurrence for each of the individual human and bacterial pathogens 
pathCols <- numeric()
for (i in 1:length(pathType)) pathCols[i] <- which(names(df)==pathType[i])

occurrence <- function(x) {sum(x != 0,na.rm=T)/sum(!is.na(x),na.rm=T)}
path_occur <- apply(df[,Presentcols],2,occurrence)
samplecount <- function(x) {sum(!is.na(x),na.rm=T)}
nvals <- apply(df[,Presentcols],2,samplecount)


for (i in 1:length(Presentcols)){
  varname <- paste(names(df)[Presentcols[i]],"occur",sep="_")
  df[,varname] <- ifelse(df[,Presentcols[i]]>0,1,0)
}

Cmonth <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

columns <- ceiling(length(Present)/3)
rows <- ceiling(length(Present)/columns)
par(mfrow=c(rows,columns))     #layout for 2 rows and 2 columns

for (i in 1:length(Presentcols)){
  varname <- paste(Present[i],"occur",sep="_")
  if(sum(df[,varname],na.rm=T)>=1) {
    V_month <- tapply(df[,varname],df$nmonth,mean,na.rm=T) # % occurrence
    V_month.df <- data.frame(V_month)
    V_month.df$month <- row.names(V_month.df)
    V_month <- rbind(rbind(rbind(rbind(V_month.df[1:6,],c(NA,6)),V_month.df[7:9,]),c(NA,10)),V_month.df[10,])
    V_month.final <- as.data.frame(t(V_month[,1]))
    names(V_month.final) <- Cmonth
    barplot(as.matrix(V_month.final), main="",
            xlab="", col=c("darkblue"),ylim=c(0,1),
            beside=TRUE,xaxt="n",yaxt="n")
    axis(side=2,at=c(0,0.5,1),labels=c("0","50%","100%"),cex.axis=1.5,las=2)
    axis(side=1,at=(seq(1:12)*2-0.5),labels=Cmonth,cex.axis=1.75,las=3)
    abline(h=0,v=0)
    mtext(names(df)[Presentcols[i]],side=3,line=0.5,cex=1.2)
  }
  
}
par(mfrow=c(1,1))     #change the graph window back to one figure
mtext("Bovine Viruses",side=3,line=2.9,cex=1.5)


# Lattice testing:

df[,Site.Name] <- factor(df[,Site.Name],levels=c("16th St.","70th St.","Honey","Underwood","Donges Bay","Meno Falls"))
df[,Site.Name] <- factor(df[,Site.Name],levels=c("Honey","70th St.","16th St.","Meno Falls","Donges Bay","Underwood"))

xyplot((df$Human_virus+.01)~df$nmonth|df$Site.Name,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Seasonal pathogens by site: Sum of Human viruses",
       scales = list(y = list(log=10)))

xyplot((df$Bacteria+.01)~df$nmonth|df$Site.Name,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Seasonal pathogens by site: Sum of Pathogenic Bacteria",
       scales = list(y = list(log=10)))

xyplot((df$All_pathogens+0.01)~df$nmonth|df$Site.Name,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Seasonal pathogens by site: Sum of All pathogens",
       scales = list(y = list(log=10)))


df[,"Event.Type"] <- factor(df[,"Event.Type"],levels=levels(df[,"Event.Type"])[c(1,4,2,3)])

bwplot((df$Human_virus+0.01)~df$Event.Type|df$Site.Name,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Event type by site: Sum of human viruses",
       scales = list(y = list(log=10)))

bwplot((df$All_pathogens+0.01)~df$Event.Type|df$Site.Name,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Event type by site: Sum of All pathogens",
       scales = list(y = list(log=10)))

bwplot((df$Human_virus+0.01)~df$Event.Type,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Event type: Sum of All pathogens",
       scales = list(y = list(log=10)))

bwplot((df$All_pathogens+0.01)~df$Event.Type,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Event type by site: Sum of All pathogens",
       scales = list(y = list(log=10)))

df[,Site.Name] <- factor(df[,Site.Name],levels=c("16th St.","70th St.","Honey","Underwood","Donges Bay","Meno Falls"))


bwplot((df$Human_virus+0.01)~df$Site.Name|df$Sample.Type,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Event type by site: Sum of human viruses",
       scales = list(y = list(log=10)))

bwplot((df$Bacteria+0.01)~df$Site.Name|df$Sample.Type,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Event type by site: Sum of human viruses",
       scales = list(y = list(log=10)))

bwplot((df$Human_virus+0.01)~df$Sample.Type,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Event type by site: Sum of human viruses",
       scales = list(y = list(log=10)))

bwplot((df$Bovine_virus+0.01)~df$Site.Name|df$Sample.Type,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Event type by site: Sum of bovine viruses",
       scales = list(y = list(log=10)))


bwplot((df$Human_virus+0.01)~df$Site.Name|df$Event.Type,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Event type by site: Sum of human viruses",
       scales = list(y = list(log=10)))

df[,Site.Name] <- factor(df[,Site.Name],levels=c("16th St.","70th St.","Honey","Underwood","Donges Bay","Meno Falls"))

bwplot((df$Human_virus+0.01)~df$Sample.Type|df$Site.Name,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Event type by site: Sum of human viruses",
       layout = c(2, 3),
       scales = list(y = list(log=10)))

bwplot((df$Human_virus+0.01)~df$Event.Type|df$Site.Name,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Event type by site: Sum of human viruses",
       layout = c(2, 3),
       scales = list(y = list(log=10)))


# generate new variable based on location of analysis: Marshfield clinic (MCRF) or USDA
df$Analysis_loc <- factor(ifelse(is.na(df$Cryptosporidium),"MCRF","USDA"),levels=c("MCRF","USDA"))
df$pdate <- strptime(df$Collection.Start.Date,"%d-%b-%y")

#one NA in 2011--need to refine
firstUSDA <- min(df[df$Analysis_loc=="USDA","pdate"])
df$Analysis_loc <- factor(ifelse(df$pdate<firstUSDA,"MCRF","USDA"),levels=c("MCRF","USDA"))

plot(df$pdate,as.numeric(df$Analysis_loc))

bwplot((df$Human_virus+0.01)~df$Analysis_loc,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Event type by site: Sum of human viruses",
       scales = list(y = list(log=10)))


############################################################################
## Plot occurrence for specific human viruses 
############################################################################

df[,Site.Name] <- factor(df[,Site.Name],levels=c("16th St.","70th St.","Honey","Underwood","Donges Bay","Meno Falls"))
Sites <- levels(df[,Site.Name])
nvals <- vector()
occur.df <- data.frame()
for (j in 1:length(Sites)) {
  
  
subdf <- subset(df,Site.Name==Sites[j])
pathogens <- dfNames[which(dfNames$Groups==pathType[1]),Variables]
for (i in 1:length(pathogens)) pathCols[i] <- which(names(subdf)==pathogens[i])

bcol <- min(pathCols)
ecol <- max(pathCols)
occurrence <- function(x) {sum(x != 0,na.rm=T)/sum(!is.na(x),na.rm=T)}
path_occur <- apply(subdf[,bcol:ecol],2,occurrence)
samplecount <- function(x) {sum(!is.na(x),na.rm=T)}
nvals <- apply(subdf[,bcol:ecol],2,samplecount)

if(j==1) occur.df <- data.frame(path_occur)
if(j>1) occur.df <- cbind(occur.df,path_occur)

if(j==1) nvals.df <- data.frame(nvals)
if(j>1) nvals.df <- cbind(nvals.df,nvals)
}


## Plot time series of human viruses by site

xyplot((df$Human_virus+.01)~df$ddate|df$Site.Name,data=df,
       xlab="date",ylab="Human Viruses (gc/L)",
       main="Seasonal pathogens by site: Sum of Human viruses",
       scales = list(y = list(log=10),x=list(rot=90)))



barchart((df$Human_virus+.01)~df$ddate|df$Site.Name,data=df,
       xlab="date",ylab="Human Viruses (gc/L)",
       main="Seasonal pathogens by site: Sum of Human viruses",
       scales = list(y = list(log=10),x=list(rot=90)))



##!!!!!!!!!!!!!!!!!NEED TO WORK ON THIS BAR GRAPH NOW--2/20/2012!!!!!!!!!!!!!!!!!!!!!!!!!!
names(occur.df) <- Sites
barplot(t(as.matrix(occur.df)), main="",
        xlab="", col=c("darkblue"),ylim=c(0,1),
        beside=TRUE)   #,xaxt="n",yaxt="n")
axis(side=2,at=c(0,0.5,1),labels=c("0","50%","100%"),cex.axis=1.5,las=2)
axis(side=1,at=(seq(1:length(Sites))*2+1)/2,labels=Sites,cex.axis=.5,las=3)
abline(h=0,v=0)
mtext(names(df)[Presentcols[i]],side=3,line=0.5,cex=1.2)


############################################################################
## Consider Human virus data before and after the lab move
############################################################################
qrange <- function(x) {
  x <- x[!is.na(x)]
  quantile(x,seq(0.6,1,0.01))
}

LocQArray <- tapply(df$Human_virus,INDEX=df$Analysis_loc,FUN=qrange)

LocQuantiles <- data.frame(MCRF=LocQArray$MCRF,USDA=LocQArray$USDA)

ylim <-with(LocQuantiles, range(c(MCRF,USDA))+0.01)
with(LocQuantiles,plot(LocQuantiles$MCRF,ylim=ylim,log="y"))
with(LocQuantiles,points(USDA,col=3))

,   scales=list(rot=c(45,45),abbreviate=T)






ylim <- range(log10(df$Human_virus+0.01))
xyplot((Human_virus+0.01)~as.Date(pdate)|df$Analysis_loc,data=df,
       xlab="Month",ylab="Pathogens (gc/L)",
       main="Event type by site: Sum of human viruses",
       scales = list(y = list(log=10),
                     rot=c(45,45),
                     ylim=ylim))

MCRFd <- sum(!df$Human_virus[df$Analysis_loc =="MCRF"]==0,na.rm=T)
MCRFnd <- sum(df$Human_virus[df$Analysis_loc =="MCRF"]==0,na.rm=T)

USDAd <- sum(!df$Human_virus[df$Analysis_loc =="USDA"]==0,na.rm=T)
USDAnd <- sum(df$Human_virus[df$Analysis_loc =="USDA"]==0,na.rm=T)

MCRFoccur <- MCRFd/(MCRFd+MCRFnd)
USDAoccur <- USDAd/(USDAd+USDAnd)
MCRFoccur;USDAoccur



x <- 1:10
y <- 21:30
d <- data.frame(dx=x,dy=y)
str(d)

plot(d$dx,d$dy)
plot(dx,dy,data=d)
