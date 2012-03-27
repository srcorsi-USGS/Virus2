# Menomonee virus Rainfall compilation
# and preliminary loading computations
#


###################
#Set R directory
SNsys <- system("wmic bios get serialnumber", intern = TRUE, show.output.on.console = FALSE)
SN <- gsub("\\s","",SNsys)[2]
if(SN == "R97R0NW") {Wlocal <- "D:/srcldata" #Laptop
                     Project <- "D://srcldata"
}else {Wlocal <- "//igsarmewwssrc/SRCdata"
       Project <- "M:/QW Monitoring Team"} #Network
Rlocal <- paste(Wlocal,"/R",sep="")

source(paste(Rlocal,"/Scripts/fxn_compute_seasonal.R",sep=""))

dfPathFIB <- read.delim(paste(Rlocal,"/MMSD_virus/Virus2/Final compiled data/Concentrations/Path_FIB_WQ_Hy_All.txt",sep=""))
PathFIBnames <- read.delim(paste(Rlocal,"/MMSD_virus/Virus2/Final compiled data/Concentrations/PathFIBnames.txt",sep=""))

dfPathFIB$Ebpdate <- strptime(dfPathFIB$Ebpdate,format="%m/%d/%Y %H:%M", tz="GMT")

dfPathFIB <- compute_seasonal(df=dfPathFIB,date="Ebpdate",return.var="jday")

logvars <- as.character(PathFIBnames[which(PathFIBnames$groups=="log"),"variables"])
logvars2 <- paste("log_",logvars,sep="")
dfPathFIB[,logvars2] <- log10(dfPathFIB[,logvars])

predictors <- PathFIBnames[which(PathFIBnames$groups=="Predictor" | PathFIBnames$groups=="log_no16"),"variables"]

predictors <- c(as.character(predictors),logvars2)
predictors <- c(predictors,"sin_jday","cos_jday")




##############################################################################
# Explore model for Human Virus
##############################################################################
#log transform Human virus data after adding 0.01 to remove zeros

min(dfPathFIB$Human_virus[dfPathFIB$Human_virus!=0],na.rm=T)
dfPathFIB$logHV <- log10(dfPathFIB$Human_virus + 0.01)

form <- formula(paste("logHV~",paste(predictors,collapse="+",sep=""),sep=""))

mHV <- lm(form,data=dfPathFIB)
summary(mHV)

mHVstep <- step(mHV,direction="both",trace=2,k=2) # Using BIC for selection k=log(nrow(bw3))
summary(mHVstep)

mHVall <- regsubsets(form,data=dfPathFIB, nvmax = 5,
                     method=c("exhaustive", "backward", "forward", "seqrep"))
summary(mHVall)
mHVall

predicted <- predict(mHVstep,newdata=dfPathFIB)
length(predicted)
logHV <- dfPathFIB$logHV
plot(predicted,logHV)

#########################################
# Explore model for Human Virus at 16th #
#########################################

subdf <- subset(dfPathFIB,Abb2=="MC")
form <- formula(paste("logHV~",paste(predictors,collapse="+",sep=""),sep=""))

mHV16 <- lm(form,data=subdf)
summary(mHV16)

mHVstep16 <- step(mHV,direction="both",trace=2,k=2) # Using BIC for selection k=log(nrow(bw3))
summary(mHVstep16)

mHVall16 <- regsubsets(form,data=subdf, nvmax = 5,
                     method=c("exhaustive", "backward", "forward", "seqrep"))
summary(mHVall16)


predicted <- predict(mHVstep16,newdata=subdf)
length(predicted)
logHV <- subdf$logHV
plot(predicted,logHV)
abline(0,1)


## RPART
mHVRP16 <- rpart(formula = form, data = subdf, 
      control = rpart.control(minsplit = 10, maxcompete = 10))
summary(mHVRP16)
plot(as.party(mHVRP16), tp_args=list(id=FALSE))
mHVRP16p <- prune(mHVRP16,cp=mHVRP16$cptable[which.min(mHVRP16$cptable[,"xerror"])])
plot(as.party(mHVRP16p), tp_args=list(id=FALSE),main="16th Street Human Virus")



##############################################
# Explore model for Human Virus at all sites #
##############################################

# Subset to one site
i <- 5
site.abb <- c("MF","LD","UW","HW","MW","MC")
subdf <- subset(dfPathFIB,Abb2==site.abb[i])

# Log transform Q and turbidity values
logvars <- as.character(PathFIBnames[which(PathFIBnames$groups=="log" | PathFIBnames$groups=="log_no16"),"variables"])
logvars2 <- paste("log_",logvars,sep="")
subdf[,logvars2] <- log10(subdf[,logvars])

#Define predictors
predictors <- PathFIBnames[which(PathFIBnames$groups2=="Predictor" | PathFIBnames$groups=="log_no16"),"variables"]
predictors <- c(as.character(predictors),logvars2)
predictors <- c(predictors,"sin_jday","cos_jday")

#log transform human virus
min(subdf$Human_virus[subdf$Human_virus!=0],na.rm=T)
subdf$logHV <- log10(subdf$Human_virus + 0.01)

#Define formula for constructing models
form <- formula(paste("logHV~",paste(predictors,collapse="+",sep=""),sep=""))

mHV <- lm(form,data=subdf)
summary(mHV)

mHVstep <- step(mHV,direction="both",trace=2,k=2) # Using BIC for selection k=log(nrow(bw3))
summary(mHVstep)

mHVall <- regsubsets(form,data=subdf, nvmax = 5,
                       method=c("exhaustive", "backward", "forward", "seqrep"))
summary(mHVall)


predicted <- predict(mHVstep,newdata=subdf)
length(predicted)
logHV <- subdf$logHV
plot(predicted,logHV)
abline(0,1)


## RPART
mHVRP <- rpart(formula = form, data = subdf, 
                 control = rpart.control(minsplit = 10, maxcompete = 10))
summary(mHVRP)
plot(as.party(mHVRP, tp_args=list(id=FALSE)))
mHVRPp <- prune(mHVRP,cp=mHVRP$cptable[which.min(mHVRP$cptable[,"xerror"])])
plot(as.party(mHVRPp), tp_args=list(id=FALSE),main="16th Street Human Virus")




##############################################
# Explore model for Enterovirus at all sites #
##############################################

# Subset to one site
i <- 5
site.abb <- c("MF","LD","UW","HW","MW","MC")
subdf <- subset(dfPathFIB,Abb2==site.abb[i])

# Log transform Q and turbidity values
logvars <- as.character(PathFIBnames[which(PathFIBnames$groups=="log" | PathFIBnames$groups=="log_no16"),"variables"])
logvars2 <- paste("log_",logvars,sep="")
subdf[,logvars2] <- log10(subdf[,logvars])

#Define predictors
predictors <- PathFIBnames[which(PathFIBnames$groups2=="Predictor" | PathFIBnames$groups=="log_no16"),"variables"]
predictors <- c(as.character(predictors),logvars2)
predictors <- c(predictors,"sin_jday","cos_jday")

#log transform human virus
min(subdf$Enterovirus[subdf$Enterovirus!=0],na.rm=T)
subdf$logEV <- log10(subdf$Enterovirus + 0.01)

#Define formula for constructing models
form <- formula(paste("logEV~",paste(predictors,collapse="+",sep=""),sep=""))

mEV <- lm(form,data=subdf)
summary(mEV)

mEVstep <- step(mEV,direction="both",trace=2,k=2) # Using BIC for selection k=log(nrow(bw3))
summary(mEVstep)

mEVall <- regsubsets(form,data=subdf, nvmax = 5,
                     method=c("exhaustive"), "backward", "forward", "seqrep"))
summary(mEVall)


predicted <- predict(mEVstep,newdata=subdf)
length(predicted)
logEV <- subdf$logEV
plot(predicted,logEV)
abline(0,1)


## RPART
mEVRP <- rpart(formula = form, data = subdf, 
               control = rpart.control(minsplit = 10, maxcompete = 10))
summary(mEVRP)
plot(as.party(mEVRP, tp_args=list(id=FALSE)))
mHVRPp <- prune(mEVRP,cp=mEVRP$cptable[which.min(mEVRP$cptable[,"xerror"])])
plot(as.party(mEVRPp), tp_args=list(id=FALSE),main="16th Street Human Virus")


##############################################################################
# Explore model for Lachno
##############################################################################
#log transform Human virus data after adding 0.01 to remove zeros

# Subset to one site
i <- 5
site.abb <- c("MF","LD","UW","HW","MW","MC")
subdf <- subset(dfPathFIB,Abb2==site.abb[i])

# Log transform Q and turbidity values
logvars <- as.character(PathFIBnames[which(PathFIBnames$groups=="log" | PathFIBnames$groups=="log_no16"),"variables"])
logvars2 <- paste("log_",logvars,sep="")
subdf[,logvars2] <- log10(subdf[,logvars])

#Define predictors
predictors <- PathFIBnames[which(PathFIBnames$groups2=="Predictor" | PathFIBnames$groups=="log_no16"),"variables"]
predictors <- c(as.character(predictors),logvars2)
predictors <- c(predictors,"sin_jday","cos_jday")

#log transform human virus
min(subdf$lachnocen[subdf$lachnocen!=0],na.rm=T)
subdf$loglachno <- log10(subdf$lachnocen)

#Define formula for constructing models
form <- formula(paste("loglachno~",paste(predictors,collapse="+",sep=""),sep=""))

mlachno <- lm(form,data=subdf)
summary(mlachno)

mlachnostep <- step(mEV,direction="both",trace=2,k=2) # Using BIC for selection k=log(nrow(bw3))
summary(mlachnostep)

mlachnoall <- regsubsets(form,data=subdf, nvmax = 5,
                     method=c("exhaustive", "backward", "forward", "seqrep"))
summary(mlachnoall)


predicted <- predict(mlachnostep,newdata=subdf)
length(predicted)
loglachno <- subdf$loglachno
plot(predicted,loglachno)
abline(0,1)


## RPART
mlachnoRP <- rpart(formula = form, data = subdf, 
               control = rpart.control(minsplit = 10, maxcompete = 10))
summary(mlachnoRP)
plot(as.party(mlachnoRP, tp_args=list(id=FALSE)))
mlachnoRPp <- prune(mlachnoRP,cp=mlachnoRP$cptable[which.min(mlachnoRP$cptable[,"xerror"])])
plot(as.party(mlachnoRPp), tp_args=list(id=FALSE),main="16th Street Human Virus")

