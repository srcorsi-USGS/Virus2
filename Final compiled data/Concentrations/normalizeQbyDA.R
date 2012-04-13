
df <- read.delim("D:/SRCData/R/MMSD_virus/Virus2/Final compiled data/Concentrations/PathWWFIBWW_stats_master.txt")

abb <- c("MF","LD","UW","HW","MW","MC")
DA <- as.numeric(c(34.7,8.0,18.2,10.3,123.,146))


df$Q_mean.DA <- numeric(nrow(df))
df$Q_max.DA <- numeric(nrow(df))
df$AntBase.DA <- numeric(nrow(df))

df$Q_mean <- as.numeric(as.character(df$Q_mean))
df$Q_max <- as.numeric(as.character(df$Q_max))

for (i in 1:length(abb)){
  site.rows <- which(df$Abbrev==abb[i])
  for (j in 1:length(site.rows)){
    df$Q_mean.DA[site.rows[j]] <- df$Q_mean[site.rows[j]]/DA[i]
    df$Q_max.DA[site.rows[j]] <- df$Q_max[site.rows[j]]/DA[i]
    df$AntBase.DA[site.rows[j]] <- df$AntBase[site.rows[j]]/DA[i]
  }
}

boxplot(df$Q_mean.DA~df$Abbrev,log="y")


boxplot(df$Q_max.DA~df$Abbrev,log="y")

boxplot(df$AntBase.DA~df$Abbrev,log="y")

write.table(df,"D:/SRCData/R/MMSD_virus/Virus2/Final compiled data/Concentrations/PathWWFIBWW_stats_master2.txt",sep="\t",row.names=F)