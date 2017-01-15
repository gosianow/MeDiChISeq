\name{remove.clons}
\alias{remove.clons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove clonal reads
}
\description{
Removes reads that starts in exactly the same position.
}
\usage{
remove.clons(file, reads.to.keep, output.dir, format)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
Path to the file containing alligned reads.
}
  \item{reads.to.keep}{
Numeric. Defines number of reads that are allowed to have exactly the same position. 
}
  \item{output.dir}{
Directory to the folder where new BED file will be saved.
}
  \item{format}{
Format of the ChIP-Seq file. Accepted formats: \code{"sam", "bam", "bowtie", "soap"}.
}
}
\details{
It can happen that profile will have a lot of reads starting in exactly the same position. Such reads are usually considered as artefacs resulting from "clonal" amplification during the PCR phase in ChIP-seq experiment. If this is a frequent phenomena they can have an influence on the shape of peaks. Such peaks follow a square pattern. That is why we suggest to remove so called clonal reads.
}
\value{
Produces new cleaned BED file in pointed directory. Output of \code{remove.clons} a path to the generated BED file.
}
\author{
Malgorzata Nowicka <gosia.nowicka@uzh.ch> and Wouter Van Gool under the direction of Marco Antonio Mendoza Parra <marco@igbmc.fr> and Hinrich Gronemeyer.
}
\examples{
dir <- system.file("extdata/CTCF", package="MeDiChISeq")
file <- file.path(dir, "GSM646314_GM12878_CTCF_rep1_rcd.bed")

remove.clons(file, reads.to.keep=2, output.dir=getwd(), format="bed")
}
