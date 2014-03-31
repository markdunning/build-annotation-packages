
library(lumiHumanIDMapping)
library(lumiMouseIDMapping)
library(lumiRatIDMapping)

library(AnnotationForge)

conn =lumiHumanIDMapping_dbconn()





source("annotFunctions.R")

####################################Human Packages###########################################

outDir ="."

###V4
manTemplate = "NewMappings.Rd"

prefix = "illuminaHumanv4"

tmp = read.table("Oct2011/human/Annotation_Illumina_Human_HT12_V4_hg19_Sept2011.txt",sep="\t", fill=TRUE, header=TRUE,as.is=TRUE,quote="", comment.char="~")


lumiTable = dbGetQuery(conn, "SELECT * FROM HumanHT12_V4_0_R2_15002873_B")

#array_address = lumiTable$Array_Address[match(tmp[,1], lumiTable$Probe_Id)]
nu_id = lumiTable$nuID[match(tmp[,1], lumiTable$Probe_Id)]


repInfo <- tmp$Repeat_Mask
repInfo[which(repInfo != "")] = paste(tmp$Repeat_Mask[which(repInfo != "")], tmp$Nr_nuc_RMasked[which(repInfo != "")], sep=":")


extraInfo = data.frame(IlluminaID = tmp[,1], ArrayAddress = tmp$Array_Address_Id,  NuID = nu_id,ProbeQuality = tmp$"Quality_score", CodingZone = tmp$Coding_zone, ProbeSequence = tmp$Probe_sequence, SecondMatches = tmp$Genomic_2nd_matches, OtherGenomicMatches = tmp$Other_genomic_matches, RepeatMask = repInfo, OverlappingSNP = tmp$SNPs, EntrezReannotated = tmp$Entrez, GenomicLocation = tmp$Genomic_location, SymbolReannotated = tmp$Gene_symbol, ReporterGroupName = tmp$Reporter_Group_Name, ReporterGroupID = tmp$Reporter_Group_ID, EnsemblReannotated = tmp$Ensembl_gene)

for(i in 1:ncol(extraInfo)){

extraInfo[which(extraInfo[,i] == ""),i] = NA

}
extraInfo$ArrayAddress <- sub("^[0]+", "", extraInfo$ArrayAddress)


makeBioconductorAnnotation("Humanv4", "HumanHT12v4", refseq  = tmp$RefSeq_transcripts, IlluminaID = tmp[,1], extraInfo = extraInfo, outDir =".", version = "1.22.0", manTemplate = manTemplate, "human")


#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD INSTALL ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)


#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD build ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)

#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD check ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)



##illuminaHumanWGDASLv4

prefix <- "illuminaHumanWGDASLv4"

tmp = read.table("Oct2011/human/Annotation_Illumina_Human_HT12_V4_WGDASL_hg19_Sept2011.txt",sep="\t", fill=TRUE, header=TRUE,as.is=TRUE,quote="", comment.char="~")


lumiTable = dbGetQuery(conn, "SELECT * FROM HumanHT12_V4_0_R2_15002873_B_WGDASL")

#array_address = lumiTable$Array_Address[match(tmp[,1], lumiTable$Probe_Id)]
nu_id = lumiTable$nuID[match(tmp[,1], lumiTable$Probe_Id)]


repInfo <- tmp$Repeat_Mask
repInfo[which(repInfo != "")] = paste(tmp$Repeat_Mask[which(repInfo != "")], tmp$Nr_nuc_RMasked[which(repInfo != "")], sep=":")



extraInfo = data.frame(IlluminaID = tmp[,1], ArrayAddress = tmp$Array_Address_Id,  NuID = nu_id,ProbeQuality = tmp$"Quality_score", CodingZone = tmp$Coding_zone, ProbeSequence = tmp$Probe_sequence, SecondMatches = tmp$Genomic_2nd_matches, OtherGenomicMatches = tmp$Other_genomic_matches, RepeatMask = repInfo, OverlappingSNP = tmp$SNPs, EntrezReannotated = tmp$Entrez, GenomicLocation = tmp$Genomic_location, SymbolReannotated = tmp$Gene_symbol, ReporterGroupName = tmp$Reporter_Group_Name, ReporterGroupID = tmp$Reporter_Group_ID, EnsemblReannotated = tmp$Ensembl_gene)


for(i in 1:ncol(extraInfo)){

extraInfo[which(extraInfo[,i] == ""),i] = NA

}
extraInfo$ArrayAddress <- sub("^[0]+", "", extraInfo$ArrayAddress)


makeBioconductorAnnotation("HumanWGDASLv4", "HumanWGDASLv4", refseq  = tmp$RefSeq_transcripts, IlluminaID = tmp[,1], extraInfo = extraInfo, outDir =".", version = "1.22.0", manTemplate = manTemplate, "human")


#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD INSTALL ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)


#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD build ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)

#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD check ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)








###V3


tmp = read.table("Oct2011/human/Annotation_Illumina_Human_HT12_V3_hg19_Sept2011.txt",sep="\t", fill=TRUE, header=TRUE,as.is=TRUE,quote="", comment.char="~")

prefix <- "illuminaHumanv3"


lumiTable = dbGetQuery(conn, "SELECT * FROM HumanHT12_V3_0_R3_11283641_A")

nu_id = lumiTable$nuID[match(tmp[,1], lumiTable$Probe_Id)]

repInfo <- tmp$Repeat_Mask
repInfo[which(repInfo != "")] = paste(tmp$Repeat_Mask[which(repInfo != "")], tmp$Nr_nuc_RMasked[which(repInfo != "")], sep=":")


extraInfo = data.frame(IlluminaID = tmp[,1], ArrayAddress = tmp$Array_Address_Id,  NuID = nu_id,ProbeQuality = tmp$"Quality_score", CodingZone = tmp$Coding_zone, ProbeSequence = tmp$Probe_sequence, SecondMatches = tmp$Genomic_2nd_matches, OtherGenomicMatches = tmp$Other_genomic_matches, RepeatMask = repInfo, OverlappingSNP = tmp$SNPs, EntrezReannotated = tmp$Entrez, GenomicLocation = tmp$Genomic_location, SymbolReannotated = tmp$Gene_symbol, ReporterGroupName = tmp$Reporter_Group_Name, ReporterGroupID = tmp$Reporter_Group_ID, EnsemblReannotated = tmp$Ensembl_gene)
for(i in 1:ncol(extraInfo)){

extraInfo[which(extraInfo[,i] == ""),i] = NA

}
extraInfo$ArrayAddress <- sub("^[0]+", "", extraInfo$ArrayAddress)

makeBioconductorAnnotation("Humanv3", "HumanHT12v3", refseq  = tmp$RefSeq_transcripts, IlluminaID = tmp[,1], extraInfo = extraInfo, outDir =".", version = "1.22.0", manTemplate = manTemplate, "human")


##Build and check in the background

#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD INSTALL ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)


#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD build ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)

#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD check ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)



##illuminaHumanWGDASLv3

tmp = read.table("Oct2011/human/Annotation_Illumina_Human_Ref-8_V3_WGDASL_hg19_Sept2011.txt",sep="\t", fill=TRUE, header=TRUE,as.is=TRUE,quote="", comment.char="~")

prefix <- "illuminaHumanWGDASLv3"

tmp$Genomic_location = gsub("t", "", tmp$Genomic_location)

lumiTable = dbGetQuery(conn, "SELECT * FROM HUMANREF8_V3_0_R1_11282963_A_WGDASL")

#array_address = lumiTable$Array_Address[match(tmp[,1], lumiTable$Probe_Id)]
nu_id = lumiTable$nuID[match(tmp[,1], lumiTable$Probe_Id)]


repInfo <- tmp$Repeat_Mask
repInfo[which(repInfo != "")] = paste(tmp$Repeat_Mask[which(repInfo != "")], tmp$Nr_nuc_RMasked[which(repInfo != "")], sep=":")



extraInfo = data.frame(IlluminaID = tmp[,1], ArrayAddress = tmp$Array_Address_Id,  NuID = nu_id,ProbeQuality = tmp$"Quality_score", CodingZone = tmp$Coding_zone, ProbeSequence = tmp$Probe_sequence, SecondMatches = tmp$Genomic_2nd_matches, OtherGenomicMatches = tmp$Other_genomic_matches, RepeatMask = repInfo, OverlappingSNP = tmp$SNPs, EntrezReannotated = tmp$Entrez, GenomicLocation = tmp$Genomic_location, SymbolReannotated = tmp$Gene_symbol, ReporterGroupName = tmp$Reporter_Group_Name, ReporterGroupID = tmp$Reporter_Group_ID, EnsemblReannotated = tmp$Ensembl_gene)

for(i in 1:ncol(extraInfo)){

extraInfo[which(extraInfo[,i] == ""),i] = NA

}
extraInfo$ArrayAddress <- sub("^[0]+", "", extraInfo$ArrayAddress)

makeBioconductorAnnotation("HumanWGDASLv3", "HumanHT12WGDASLv3", refseq  = tmp$RefSeq_transcripts, IlluminaID = tmp[,1], extraInfo = extraInfo, outDir =".", version = "1.22.0", manTemplate = manTemplate, "human")


#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD INSTALL ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)


#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD build ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)

#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD check ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)






###V2



tmp = read.table("Oct2011/human/Annotation_Illumina_Human_WG-6_V2_hg19_Sept2011.txt",sep="\t", fill=TRUE, header=TRUE,as.is=TRUE,quote="", comment.char="~")

prefix <- "illuminaHumanv2"

tmp$Genomic_location = gsub("t", "", tmp$Genomic_location)

lumiTable = dbGetQuery(conn, "SELECT * FROM HumanWG6_V2_0_R4_11223189_A")

#array_address = lumiTable$Array_Address[match(tmp[,1], lumiTable$Probe_Id)]
nu_id = lumiTable$nuID[match(tmp[,1], lumiTable$Probe_Id)]

repInfo <- tmp$Repeat_Mask
repInfo[which(repInfo != "")] = paste(tmp$Repeat_Mask[which(repInfo != "")], tmp$Nr_nuc_RMasked[which(repInfo != "")], sep=":")


extraInfo = data.frame(IlluminaID = tmp[,1], ArrayAddress = tmp$Array_Address_Id,  NuID = nu_id,ProbeQuality = tmp$"Quality_score", CodingZone = tmp$Coding_zone, ProbeSequence = tmp$Probe_sequence, SecondMatches = tmp$Genomic_2nd_matches, OtherGenomicMatches = tmp$Other_genomic_matches, RepeatMask = repInfo, OverlappingSNP = tmp$SNPs, EntrezReannotated = tmp$Entrez, GenomicLocation = tmp$Genomic_location, SymbolReannotated = tmp$Gene_symbol, ReporterGroupName = tmp$Reporter_Group_Name, ReporterGroupID = tmp$Reporter_Group_ID, EnsemblReannotated = tmp$Ensembl_gene)
for(i in 1:ncol(extraInfo)){

extraInfo[which(extraInfo[,i] == ""),i] = NA

}
extraInfo$ArrayAddress <- sub("^[0]+", "", extraInfo$ArrayAddress)


makeBioconductorAnnotation("Humanv2", "HumanWG6v2", refseq  = tmp$RefSeq_transcripts, IlluminaID = tmp[,1], extraInfo = extraInfo, outDir =".", version = "1.22.0", manTemplate = manTemplate, "human")

##Build and check in the background

#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD INSTALL ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)


#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD build ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)

#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD check ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)


###V1



tmp = read.table("Oct2011/human/Annotation_Illumina_Human_WG-6_V1_hg19_Sept2011.txt",sep="\t", fill=TRUE, header=TRUE,as.is=TRUE,quote="", comment.char="~")

prefix <- "illuminaHumanv1"

tmp$Genomic_location = gsub("t", "", tmp$Genomic_location)

lumiTable = dbGetQuery(conn, "SELECT * FROM HumanWG6_V1")

array_address = lumiTable$Array_Address[match(tmp[,1], lumiTable$Probe_Id)]
nu_id = lumiTable$nuID[match(tmp[,1], lumiTable$Probe_Id)]

repInfo <- tmp$Repeat_Mask
repInfo[which(repInfo != "")] = paste(tmp$Repeat_Mask[which(repInfo != "")], tmp$Nr_nuc_RMasked[which(repInfo != "")], sep=":")


extraInfo = data.frame(IlluminaID = tmp[,1], ArrayAddress = tmp$Array_Address_Id,  NuID = nu_id,ProbeQuality = tmp$"Quality_score", CodingZone = tmp$Coding_zone, ProbeSequence = tmp$Probe_sequence, SecondMatches = tmp$Genomic_2nd_matches, OtherGenomicMatches = tmp$Other_genomic_matches, RepeatMask = repInfo, OverlappingSNP = tmp$SNPs, EntrezReannotated = tmp$Entrez, GenomicLocation = tmp$Genomic_location, SymbolReannotated = tmp$Gene_symbol, ReporterGroupName = tmp$Reporter_Group_Name, ReporterGroupID = tmp$Reporter_Group_ID, EnsemblReannotated = tmp$Ensembl_gene)



###Load the control info

controlz <- read.csv("Humanv1Controls.csv")
controlz <- data.frame(IlluminaID = paste("CONTROL",controlz[,1],sep="_"), ArrayAddress = controlz[,1], ReporterGroupName = controlz[,2], ReporterGroupID = controlz[,3], ProbeSequence = controlz[,4])

extraInfo <- merge(extraInfo, controlz,all=T)

for(i in 1:ncol(extraInfo)){

extraInfo[which(extraInfo[,i] == ""),i] = NA

}
extraInfo$ArrayAddress <- sub("^[0]+", "", extraInfo$ArrayAddress)

makeBioconductorAnnotation("Humanv1", "HumanWG6v1", refseq  = tmp$RefSeq_transcripts, IlluminaID = tmp[,1], extraInfo = extraInfo, outDir =".", version = "1.22.0", manTemplate = manTemplate, "human")


##Build and check in the background

#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD INSTALL ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)


#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD build ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)
#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD check ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)



####################################Mouse Packages###########################################


conn =lumiMouseIDMapping_dbconn()


tmp = read.table("Oct2011/mouse/Annotation_Illumina_Mouse_WG-6_V2_mm9_Sept2011.txt",sep="\t", fill=TRUE, header=TRUE,as.is=TRUE,quote="", comment.char="~")

prefix <- "illuminaMousev2"


lumiTable = dbGetQuery(conn, "SELECT * FROM MouseWG6_V2_0_R3_11278593_A")

#array_address = lumiTable$ProbeId[match(tmp[,1], lumiTable$ProbeId)]
nu_id = lumiTable$nuID[match(tmp[,1], lumiTable$Probe_Id)]

repInfo <- tmp$Repeat_Mask
repInfo[which(repInfo != "")] = paste(tmp$Repeat_Mask[which(repInfo != "")], tmp$Nr_nuc_RMasked[which(repInfo != "")], sep=":")


extraInfo = data.frame(IlluminaID = tmp[,1], ArrayAddress = tmp$Array_Address_Id,  NuID = nu_id,ProbeQuality = tmp$"Quality_score", CodingZone = tmp$Coding_zone, ProbeSequence = tmp$Probe_sequence, SecondMatches = tmp$Genomic_2nd_matches, OtherGenomicMatches = tmp$Other_genomic_matches, RepeatMask = repInfo, OverlappingSNP = tmp$SNPs, EntrezReannotated = tmp$Entrez, GenomicLocation = tmp$Genomic_location, SymbolReannotated = tmp$Gene_symbol, ReporterGroupName = tmp$Reporter_Group_Name, ReporterGroupID = tmp$Reporter_Group_ID, EnsemblReannotated = tmp$Ensembl_gene)
for(i in 1:ncol(extraInfo)){

extraInfo[which(extraInfo[,i] == ""),i] = NA

}
extraInfo$ArrayAddress <- sub("^[0]+", "", extraInfo$ArrayAddress)

makeBioconductorAnnotation("Mousev2", "MouseWG6v2", refseq  = tmp$RefSeq_transcripts, IlluminaID = tmp[,1], extraInfo = extraInfo, outDir =".", version = "1.22.0", manTemplate = manTemplate, "mouse")




#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD INSTALL ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)


#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD build ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)
#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD check ", outDir, "/", prefix, ".db",sep=""))







tmp = read.table("Oct2011/mouse/Annotation_Illumina_Mouse_WG-6_V1_1_mm9_Sept2011.txt",sep="\t", fill=TRUE, header=TRUE,as.is=TRUE,quote="", comment.char="~")

prefix <- "illuminaMousev1p1"


lumiTable = dbGetQuery(conn, "SELECT * FROM MouseWG6_V1_B")

array_address = lumiTable$ProbeId[match(tmp[,1], lumiTable$ProbeId)]
nu_id = lumiTable$nuID[match(tmp[,1], lumiTable$ProbeId)]

repInfo <- tmp$Repeat_Mask
repInfo[which(repInfo != "")] = paste(tmp$Repeat_Mask[which(repInfo != "")], tmp$Nr_nuc_RMasked[which(repInfo != "")], sep=":")


extraInfo = data.frame(IlluminaID = tmp[,1], ArrayAddress = tmp$Array_Address_Id,  NuID = nu_id,ProbeQuality = tmp$"Quality_score", CodingZone = tmp$Coding_zone, ProbeSequence = tmp$Probe_sequence, SecondMatches = tmp$Genomic_2nd_matches, OtherGenomicMatches = tmp$Other_genomic_matches, RepeatMask = repInfo, OverlappingSNP = tmp$SNPs, EntrezReannotated = tmp$Entrez, GenomicLocation = tmp$Genomic_location, SymbolReannotated = tmp$Gene_symbol, ReporterGroupName = tmp$Reporter_Group_Name, ReporterGroupID = tmp$Reporter_Group_ID, EnsemblReannotated = tmp$Ensembl_gene)
for(i in 1:ncol(extraInfo)){

extraInfo[which(extraInfo[,i] == ""),i] = NA

}
extraInfo$ArrayAddress <- sub("^[0]+", "", extraInfo$ArrayAddress)

makeBioconductorAnnotation("Mousev1p1", "MouseWG6v1p1", refseq  = tmp$RefSeq_transcripts, IlluminaID = tmp[,1], extraInfo = extraInfo, outDir =".", version = "1.22.0", manTemplate = manTemplate, "mouse")




#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD INSTALL ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)


#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD build ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)
#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD check ", outDir, "/", prefix, ".db",sep=""))







tmp = read.table("Oct2011/mouse/Annotation_Illumina_Mouse_WG-6_V1_mm9_Sept2011.txt",sep="\t", fill=TRUE, header=TRUE,as.is=TRUE,quote="", comment.char="~")

prefix <- "illuminaMousev1"


lumiTable = dbGetQuery(conn, "SELECT * FROM MouseWG6_V1_B")

array_address = lumiTable$ProbeId[match(tmp[,1], lumiTable$ProbeId)]
nu_id = lumiTable$nuID[match(tmp[,1], lumiTable$ProbeId)]

repInfo <- tmp$Repeat_Mask
repInfo[which(repInfo != "")] = paste(tmp$Repeat_Mask[which(repInfo != "")], tmp$Nr_nuc_RMasked[which(repInfo != "")], sep=":")


extraInfo = data.frame(IlluminaID = tmp[,1], ArrayAddress = tmp$Array_Address_Id,  NuID = nu_id,ProbeQuality = tmp$"Quality_score", CodingZone = tmp$Coding_zone, ProbeSequence = tmp$Probe_sequence, SecondMatches = tmp$Genomic_2nd_matches, OtherGenomicMatches = tmp$Other_genomic_matches, RepeatMask = repInfo, OverlappingSNP = tmp$SNPs, EntrezReannotated = tmp$Entrez, GenomicLocation = tmp$Genomic_location, SymbolReannotated = tmp$Gene_symbol, ReporterGroupName = tmp$Reporter_Group_Name, ReporterGroupID = tmp$Reporter_Group_ID, EnsemblReannotated = tmp$Ensembl_gene)


controlz <- read.csv("Mousev1Controls.csv")
controlz <- data.frame(IlluminaID = paste("CONTROL",controlz[,1],sep="_"), ArrayAddress = controlz[,1], ReporterGroupName = controlz[,2], ReporterGroupID = controlz[,3], ProbeSequence = controlz[,4])

extraInfo <- merge(extraInfo, controlz,all=T)


for(i in 1:ncol(extraInfo)){

extraInfo[which(extraInfo[,i] == ""),i] = NA

}
extraInfo$ArrayAddress <- sub("^[0]+", "", extraInfo$ArrayAddress)

makeBioconductorAnnotation("Mousev1", "MouseWG6v1", refseq  = tmp$RefSeq_transcripts, IlluminaID = tmp[,1], extraInfo = extraInfo, outDir =".", version = "1.22.0", manTemplate = manTemplate, "mouse")




#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD INSTALL ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)


#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD build ", outDir, "/", prefix, ".db",sep=""),wait=FALSE)
#system(paste("/home/dunnin01/software/R-3.0.2/bin/R CMD check ", outDir, "/", prefix, ".db",sep=""))






####################################Rat Packages###########################################

conn =lumiRatIDMapping_dbconn()


tmp = read.table("Oct2011/rat/Annotation_Illumina_Rat_Ref-12_V1_rn4_Sept2011.txt",sep="\t", fill=TRUE, header=TRUE,as.is=TRUE,quote="", comment.char="~")

prefix <- "illuminaRatv1"


lumiTable = dbGetQuery(conn, "SELECT * FROM RatRef12_V1_0_R5_11222119_A")

array_address = lumiTable$ProbeId[match(tmp[,1], lumiTable$ProbeId)]
nu_id = lumiTable$nuID[match(tmp[,1], lumiTable$ProbeId)]

repInfo <- tmp$Repeat_Mask
repInfo[which(repInfo != "")] = paste(tmp$Repeat_Mask[which(repInfo != "")], tmp$Nr_nuc_RMasked[which(repInfo != "")], sep=":")


extraInfo = data.frame(IlluminaID = tmp[,1], ArrayAddress = tmp$Array_Address_Id,  NuID = nu_id,ProbeQuality = tmp$"Quality_score", CodingZone = tmp$Coding_zone, ProbeSequence = tmp$Probe_sequence, SecondMatches = tmp$Genomic_2nd_matches, OtherGenomicMatches = tmp$Other_genomic_matches, RepeatMask = repInfo, OverlappingSNP = tmp$SNPs, EntrezReannotated = tmp$Entrez, GenomicLocation = tmp$Genomic_location, SymbolReannotated = tmp$Gene_symbol, ReporterGroupName = tmp$Reporter_Group_Name, ReporterGroupID = tmp$Reporter_Group_ID, EnsemblReannotated = tmp$Ensembl_gene)
for(i in 1:ncol(extraInfo)){

extraInfo[which(extraInfo[,i] == ""),i] = NA

}
extraInfo$ArrayAddress <- sub("^[0]+", "", extraInfo$ArrayAddress)


makeBioconductorAnnotation("Ratv1", "Ratv1", refseq  = tmp$RefSeq_transcripts, IlluminaID = tmp[,1], extraInfo = extraInfo, outDir =".", version = "1.22.0", manTemplate = manTemplate, "rat")



######Manual steps
##1. Remove NAMESPACE.original from each package
##2. Add imports AnnotationForge to DESCRIPTION file
##3. Add 'import(AnnotationForge)' to NAMESPACE
##4. Remove zzz.R.original from R/
##5. Add display function to R/zzz.R
##e.g.
##illuminaHumanv4 <- function() {
##cat("####Mappings based on RefSeqID####\n")
##showQCData("illuminaHumanv4", datacache)
##cat("####Custom Mappings based on probe sequence####\n")
##illuminaHumanv4listNewMappings()

##}



