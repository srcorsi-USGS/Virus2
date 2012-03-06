str(df)
str(mendates)

mendates$USGSID <- paste(substr(mensum$SampleID,1,2),substr(mensum$SampleID,3,6))
df.merged <- merge(df,mendates,by.x="Study.Sample.ID",by.y="USGSID")

write.table(df.merged,"virus2.with.datestimes.txt",sep="\t",row.names=FALSE)

df.merged$ESdate <- NA
df.merged$EEdate <- NA

FIB$USGSID <- paste(substr(FIB$Site.Code,5,6),substr(FIB$Site.Code,7,9))

abbrev <- substr(FIB$Site.Code,5,6)

which(!abbrev==FIB$Abb.)