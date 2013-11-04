\name{PKGNAMElistNewMappings}
\alias{PKGNAMEARRAYADDRESS}
\alias{PKGNAMENUID}
\alias{PKGNAMEPROBESEQUENCE}
\alias{PKGNAMEPROBEQUALITY}
\alias{PKGNAMECODINGZONE}
\alias{PKGNAMEGENOMICLOCATION}
\alias{PKGNAMEGENOMICMATCHSIMILARITY}
\alias{PKGNAMESECONDMATCHES}
\alias{PKGNAMESECONDMATCHSIMILARITY}
\alias{PKGNAMETRANSCRIPTOMICMATCHSIMILARITY}
\alias{PKGNAMEOTHERGENOMICMATCHES}
\alias{PKGNAMEREPEATMASK}
\alias{PKGNAMEOVERLAPPINGSNP}
\alias{PKGNAMEENTREZREANNOTATED}
\alias{PKGNAMEENSEMBLREANNOTATED}
\alias{PKGNAMESYMBOLREANNOTATED}
\alias{PKGNAMElistNewMappings}
\alias{PKGNAMEfullReannotation}
\alias{PKGNAMEREPORTERGROUPNAME}
\alias{PKGNAMEREPORTERGROUPID}

\title{Custom mappings added to the package}
\description{
We have used an extensive re-annotation of the PKGNAME probe sequences to provide additional information that is not captured in the standard Bioconductor packages. Whereas Bioconductor annotations are based on the RefSeq ID that each probe maps to, our additional mappings provide data specific to each probe on the platform. See below for details. We recommend using the probe quality as a form of filtering, and retaining only perfect or good probes for an analysis. 

}

\section{Details of custom mappings}{

\describe{

\item{PKGNAMElistNewMappings}{List all the custom re-annotation mappings provided by the package}

\item{PKGNAMEfullReannotation}{Return all the re-annotation information as a matrix}

\item{PKGNAMEARRAYADDRESS}{Array Address code used to identify the probe at the bead-level}

\item{PKGNAMENUID}{Lumi's nuID (universal naming scheme for oligonucleotides) Reference: Du et al. (2007), Biol Direct 2:16}

\item{PKGNAMEPROBESEQUENCE}{The 50 base sequence for the probe}

\item{PKGNAMEPROBEQUALITY}{Quality grade assigned to the probe: \dQuote{Perfect} if it perfectly and uniquely matches the target transcript; \dQuote{Good} if the probe, although imperfectly matching the target transcript, is still likely to provide considerably sensitive signal (up to two mismatches are allowed, based on empirical evidence that the signal intensity for 50-mer probes with less than 95\% identity to the respective targets is less than 50\% of the signal associated with perfect matches *); \dQuote{Bad} if the probe matches repeat sequences, intergenic or intronic regions, or is unlikely to provide specific signal for any transcript; \dQuote{No match} if it does not match any genomic region or transcript.}

\item{PKGNAMECODINGZONE}{Coding status of target sequence: intergenic / intronic / Transcriptomic (\dQuote{Transcriptomic} when the target transcript is non-coding or there is no information on the coding sequence)}

\item{PKGNAMEGENOMICLOCATION}{Probe's genomic coordinates (hg19 for human, mm9 for mouse or rn4 for rat)}

\item{PKGNAMEGENOMICMATCHSIMILARITY}{Percentage of similarity between the probe and its best genomic match in the alignable region, taking the probe as reference}

\item{PKGNAMESECONDMATCHES}{Genomic coordinates of second best matches between the probe and the genome}

\item{PKGNAMESECONDMATCHSIMILARITY}{Percentage of similarity between the probe and its second best genomic match in the alignable region, taking the probe as reference}

\item{PKGNAMETRANSCRIPTOMICMATCHSIMILARITY}{Percentage of similarity between the probe and its target transcript in the alignable region, taking the probe as reference}

\item{PKGNAMEOTHERGENOMICMATCHES}{Genomic coordinates of sequences as alignable with the probe (in terms of number of matching nucleotides) as its main target}

\item{PKGNAMEREPEATMASK}{Overlapping RepeatMasked sequences, with number of bases overlapped by the repeat}

\item{PKGNAMEOVERLAPPINGSNP}{Overlapping annotated SNPs}

\item{PKGNAMEENTREZREANNOTATED}{Entrez IDs}

\item{PKGNAMEENSEMBLREANNOTATED}{Ensembl IDs}

\item{PKGNAMESYMBOLREANNOTATED}{Gene symbol derived by re-annotation}

\item{PKGNAMEREPORTERGROUPID}{For probes marked as controls in Illuminas annotation file, these gives the type of control}

\item{PKGNAMEREPORTERGROUPNAME}{Usually a more informative name for the control type}



}

}

\references{
  \url{http://remoat.sysbiol.cam.ac.uk}


 Barbosa-Morais et al. (2010) A re-annotation pipeline for Illumina BeadArrays: improving the interpretation of gene expression data. Nucleic Acids Research

}
\examples{

	##See what new mappings are available

	PKGNAMElistNewMappings()
	

        x <- PKGNAMEPROBEQUALITY

        mapped_probes <- mappedkeys(x)
        # Convert to a list
        xx <- as.list(x[mapped_probes])
        if(length(xx) > 0) {
          # Get the PROBEQUALITY for the first five probes
          xx[1:5]
          # Get the first one
          xx[[1]]
        }


	##Overall table of qualities
	table(unlist(xx))

	

        x <- PKGNAMEARRAYADDRESS

        mapped_probes <- mappedkeys(x)
        # Convert to a list
        xx <- as.list(x[mapped_probes])
        if(length(xx) > 0) {
          # Get the ARRAYADDRESS for the first five probes
          xx[1:5]
          # Get the first one
          xx[[1]]
        }

	##Can do the mapping from array address to illumina ID using a revmap
	
	y<- revmap(PKGNAMEARRAYADDRESS)
	
        mapped_probes <- mappedkeys(y)
        # Convert to a list
        yy <- as.list(y[mapped_probes])
        if(length(yy) > 0) {
          # Get the ARRAYADDRESS for the first five probes
          yy[1:5]
          # Get the first one
          yy[[1]]
        }
	


        x <- PKGNAMECODINGZONE

        mapped_probes <- mappedkeys(x)
        # Convert to a list
        xx <- as.list(x[mapped_probes])
        if(length(xx) > 0) {
          # Get the CODINGZONE for the first five probes
          xx[1:5]
          # Get the first one
          xx[[1]]
        }

        x <- PKGNAMEPROBESEQUENCE

        mapped_probes <- mappedkeys(x)
        # Convert to a list
        xx <- as.list(x[mapped_probes])
        if(length(xx) > 0) {
          # Get the PROBESEQUENCE for the first five probes
          xx[1:5]
          # Get the first one
          xx[[1]]
        }


}
\keyword{datasets}
