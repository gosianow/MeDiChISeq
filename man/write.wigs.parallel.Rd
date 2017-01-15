\name{write.wigs.parallel}
\alias{write.wigs.parallel}
\title{
WIG generator
}
\description{
Produces WIG files from files containing aligned reads.
}
\usage{
write.wigs.parallel(file, output.dir, chromosomes, sample.type, 
wig.res, reads.elong, split = FALSE, genome, format, zeros = F, 
nr.cores = 1, overwrite.wigs = FALSE, verbose = T)
}
\arguments{
  \item{file}{
Path to the file containing aligned reads.
}
  \item{output.dir}{
Directory to the folder where all the WIGS will be saved.
}
  \item{chromosomes}{
List of chromosomes for which the WIGS will be produced.
}
  \item{sample.type}{
Name prefix of the new WIG files.
}
  \item{wig.res}{
Resolution of the wigs.
}
  \item{reads.elong}{
Numeric. Reads are elongated by this value in order to create intensity wig files.
}
  \item{split}{
Logical. If \code{TRUE}, separet wigs for forward and reverse strand will be produced. 
}
  \item{genome}{
Chatacter specifying the genome \code{"hg18"}, \code{"hg19"}, \code{"mm8"}, \code{"mm9"}, \code{"dm2"}, \code{"dm3"}, \code{"ce4"}, \code{"ce6"}, \code{"rn3"}, \code{"rn4"}, \code{"danRer4"}, \code{"danRer6"}.
}
  \item{format}{
Format of the ChIP-Seq file. Accepted formats: \code{"sam"}, \code{"bam"}, \code{"bowtie"}, \code{"soap"}.
}
  \item{zeros}{
Logical. Whether to produce the entries with intensity zero.
}
  \item{nr.cores}{
Numeric. The number of cores to be used for parallel processing.
}
  \item{overwrite.wigs}{
Logical. Whether to overwrite wig files with the existing ones. When \code{FALSE} one can save time on wig construction.  
}
  \item{verbose}{
Logical. Whether to display processing status.
}
}
\value{
Produces the WIG files in pointed directory. Output of \code{write.wigs.parallel} is a list of paths to the generated WIGS.
}
\author{
Malgorzata Nowicka <gosia.nowicka@uzh.ch> and Wouter Van Gool under the direction of Marco Antonio Mendoza Parra <marco@igbmc.fr> and Hinrich Gronemeyer.
}
\examples{

dir <- system.file("extdata/CTCF", package="MeDiChISeq")
file <- file.path(dir, "GSM646314_GM12878_CTCF_rep1_rcd.bed")

write.wigs.parallel(file, output.dir=getwd(), chromosomes="chr19", sample.type="IP", 
wig.res=10, reads.elong=150, split = FALSE, genome="hg18", format="bed", zeros = FALSE, 
nr.cores = 1, overwrite.wigs = FALSE, verbose = TRUE)
}
