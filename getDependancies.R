source("http://www.bioconductor.org/biocLite.R")

##Make sure that these are the latest versions
biocLite(c("AnnotationDbi","AnnotationForge", "GO.db","human.db0", "org.Hs.eg.db", "mouse.db0", "org.Mm.eg.db", "rat.db0", "org.Rn.eg.db"))

biocLite(c("lumiHumanIDMapping", "lumiMouseIDMapping", "lumiRatIDMapping"))

