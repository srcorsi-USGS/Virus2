## Script to merge the six menomonee files into one combined file

###################
#Set R directory
SNsys <- system("wmic bios get serialnumber", intern = TRUE, show.output.on.console = FALSE)
SN <- gsub("\\s","",SNsys)[2]
if(SN == "R97R0NW") {Wlocal <- "D:/srcldata" #Laptop
                     Project <- "D://srcldata"
}else {Wlocal <- "//igsarmewwssrc/SRCdata"
       Project <- "M:/QW Monitoring Team"} #Network
Rlocal <- paste(Wlocal,"/R",sep="")

#write.table(Pathdf,paste(dir.data,"Path_FIB_WQ_Hy_",UV.site[site.num],".txt",sep=""),sep="\t",row.names=F)

UV.site <- c("MenoFalls","LMDonges","Underwood","Honey","70th","16th")
site <- c("Meno Falls","Donges Bay","Underwood","Honey","70th St.","16th St.")

dir.data <- paste(Rlocal,"/MMSD_virus/Virus2/Final compiled data/Concentrations/",sep="")

name.check <- numeric()
rows <- 0
for (site.num in 1:length(site)){
  Pathdf.site <- read.delim(paste(dir.data,"Path_FIB_WQ_Hy_",UV.site[site.num],".txt",sep=""))
  if(site.num==1) Pathdf <- Pathdf.site
  else Pathdf <- merge(Pathdf,Pathdf.site,all=T)
  name.check[site.num] <- sum(names(Pathdf.site)!=names(Pathdf))
  rows <- rows + nrow(Pathdf.site)
}
  
write.table(Pathdf,paste(dir.data,"Path_FIB_WQ_Hy_All.txt",sep=""),sep="\t",row.names=F)
