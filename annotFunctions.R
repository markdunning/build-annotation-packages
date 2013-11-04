AddNewBimap = function(file, name, Index, prefix){
	

	cmd = paste(prefix, toupper(name), " <- createSimpleBimap(\"ExtraInfo\", \"", Index, "\",","\"", name, "\"",", datacache, \"", toupper(name), "\",", "\"",prefix, ".db\")\n",sep="")

	cat(cmd, file= file, append=TRUE)
	cat("\n", file= file, append=TRUE)

}



insertExtraInfo = function(dbcon, extraInfo){	

	names = colnames(extraInfo)


	cmd = "INSERT INTO ExtraInfo VALUES ("

	for(i in 1:(length(names)-1)){
	
		cmd = paste(cmd, "$", names[i], ", ", sep="")
	}

	cmd = paste(cmd, "$", names[length(names)], ")",sep="")

bval <- dbBeginTransaction(dbcon)

gval <- dbGetPreparedQuery(dbcon, cmd, bind.data = extraInfo)

cval <- dbCommit(dbcon)


}

makeSqlTable= function(dbcon, names){

	
	cmd = "CREATE Table ExtraInfo ("

	for(i in 1:(length(names)-1)){

		cmd = paste(cmd, names[i], " TEXT, ",sep="")
	}



	cmd = paste(cmd, names[length(names)], " TEXT)",sep="")

	dbGetQuery(dbcon, cmd)


}



makeBioconductorAnnotation = function(baseName, chipName, refseq, IlluminaID, extraInfo, outDir, version, manTemplate, organism = "human"){

	##transform refseq into form required by AnnotationDbi	
	refseq = lapply(as.character(refseq), function(x) gsub('[[:space:]]', ';', x)) 

	rs = paste(baseName, "_refseq.txt",sep="")
	
	###write it out to a file

	write.table(cbind(IlluminaID, refseq), file=rs, 
         sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

	prefix=paste("illumina", baseName,sep="")

	if(organism == "human"){

          makeDBPackage("HUMANCHIP_DB",affy=FALSE,
               prefix=prefix,
               fileName=rs,
               baseMapType="refseq",
               outputDir = outDir,
               version=version,
               manufacturer = "Illumina",
               chipName = chipName,
               manufacturerUrl = "http://www.illumina.com",
		author ="Mark Dunning, Andy Lynch, Matthew Eldridge",
		maintainer="Mark Dunning <mark.dunning@cruk.cam.ac.uk>"
	)

	}

	else if (organism == "mouse"){

          makeDBPackage("MOUSECHIP_DB",affy=FALSE,
               prefix=prefix,
               fileName=rs,
               baseMapType="refseq",
               outputDir = outDir,
               version=version,
               manufacturer = "Illumina",
               chipName = chipName,
               manufacturerUrl = "http://www.illumina.com",
		author ="Mark Dunning, Andy Lynch, Matthew Eldridge",
		maintainer="Mark Dunning <mark.dunning@cruk.cam.ac.uk>"
	)

	}

	else if (organism == "rat"){

	makeDBPackage("RATCHIP_DB",affy=FALSE,
               prefix=prefix,
               fileName=rs,
               baseMapType="refseq",
               outputDir = outDir,
               version=version,
               manufacturer = "Illumina",
               chipName = chipName,
               manufacturerUrl = "http://www.illumina.com",
		author ="Mark Dunning, Andy Lynch, Matthew Eldridge",
		maintainer="Mark Dunning <mark.dunning@cruk.cam.ac.uk>"
	)
	
	


	}

	else stop("Invalid organism definiton\n")


	newPkgPath = paste(outDir, "/",prefix, ".db",sep="")

	newSQL = paste(newPkgPath,"/inst/extdata/", prefix,".sqlite", sep="")	

	###Make the new SQL file writable

	system(paste("chmod 755", newSQL))

	drv = dbDriver("SQLite")
	dbcon = dbConnect(drv, dbname=newSQL)


	makeSqlTable(dbcon, colnames(extraInfo))


	insertExtraInfo(dbcon, extraInfo) 


	cat("Checking that insert worked\n")

	dbGetQuery(dbcon, "SELECT * FROM ExtraInfo LIMIT 10")

	

	#sqlCreate = "CREATE Table ExtraInfo (IlluminaID TEXT, ArrayAddress TEXT, ProbeQuality TEXT, CodingZone TEXT, ProbeSequence TEXT, OtherMatches TEXT)"

	#dbGetQuery(dbcon, sqlCreate)

	#sqlInsert <- "INSERT INTO ExtraInfo VALUES ($IlluminaID, $ArrayAddress, $ProbeQuality, $CodingZone, $ProbeSequence, $OtherMatches)"

	#bval <- dbBeginTransaction(dbcon)

	#gval <- dbGetPreparedQuery(dbcon, sqlInsert, bind.data = extraInfo)

	#cval <- dbCommit(dbcon)

	zzzFile = paste(newPkgPath, "/R/zzz.R",sep="")

	##make a copy of zzz.r 

	system(paste("cp ", zzzFile, " ", zzzFile,".original",sep="")) 


	###Need to add these kinds of definitions to the zzz file

	##illuminaHumanv3PROBEQUALITY <- createSimpleBimap("probeinfo","ProbeID","ProbeQuality",datacache,"PROBEQUALITY","illuminaHumanv3.db")
	##illuminaHumanv3CODINGZONE <- createSimpleBimap("probeinfo","ProbeID","CodingZone",datacache,"CODINGZONE","illuminaHumanv3.db")

	cat("##Custom Bimaps for the package\n\n", file=zzzFile,append=TRUE)

	###We assume that column 2 of ExtraInfo is ArrayAddressID


	for(i in 2:ncol(extraInfo)){
	
		AddNewBimap(zzzFile, colnames(extraInfo)[i], colnames(extraInfo)[1], prefix)

	}		


	##Export them in the namespace
	
	nspace = paste(newPkgPath, "/NAMESPACE",sep="")

	system(paste("cp ", nspace, " ", nspace,".original",sep="")) 


	cat("##Custom Bimaps exported\n\n", file=nspace,append=TRUE)
	
	newFns = NULL

	for(i in 2:ncol(extraInfo)){

		newFns[i-1] = paste(prefix, toupper(colnames(extraInfo)[i]),sep="")
	}

	

	
	cat(paste("export(", paste(newFns, collapse=" , "), ", ", prefix,"listNewMappings ,", prefix, "fullReannotation)\n",sep=""),file=nspace, append=TRUE)
	

	cat(zzzFile, "##Define a utility function that lists the funtions we've just made\n",append=TRUE)

	cmd = paste(prefix,"listNewMappings = function(){\n",sep="")

	for(i in 1:length(newFns)){
		cmd = paste(cmd, "cat(\"", newFns[i], "()",sep="")
		
		cmd = paste(cmd, "\\n\")",sep="")

		cmd = paste(cmd,"\n") 
	}

	cmd = paste(cmd, "}\n",sep="")
	
	cat(cmd, file= zzzFile, append=TRUE)


	###utility function to retrieve the full table

	cmd = paste(prefix,"fullReannotation = function(){\n",sep="")

	cmd = paste(cmd, "dbGetQuery(", prefix, "_dbconn(), \"SELECT * FROM ExtraInfo\")\n}\n",sep="")


	cat(cmd, file= zzzFile, append=TRUE)


	###Now copy the template doc

	newManPage = paste(outDir, "/",prefix,".db/man/", prefix, "NewMappings.Rd",sep="")
		
	system(paste("cp ", manTemplate, " ", newManPage,sep=""))
	tmp = system(paste("sed \'s/PKGNAME/",prefix, "/g' ", newManPage,sep=""),intern=TRUE)
	write(tmp, newManPage)

}

