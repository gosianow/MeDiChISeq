\name{deconv.entire.genome.seq}
\alias{deconv.entire.genome.seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High-resolution model-based deconvolution of normalized ChIP-Seq data
}
\description{

Deconvolves an entire ChIP-Seq data set, across all chromosomes by running adjusted version of \code{chip.deconv} function from \code{MeDiChI} package on multiple, contiguous, overlapping subsets of the data and combining the results.

}
\usage{
deconv.entire.genome.seq(file.IP, file.Control=NULL, 
format="bed", genome="hg19", output.dir=NULL, output.name=NULL, 
chrom.list=NULL, limL=0, limU=Inf, potential.peaks=NULL, 
reads.elong=150, kernel, frag.length, 
quant.cutoff="q1e-5", window=20000, wig.res=10, fit.res=50, 
max.steps=100, selection.method="bic", post.proc.factor=2,
nr.boots=5, local.windows=c(1000, 2000, 5000), 
Control.corr.param=0.01, nr.cores=1, 
remove.clonal.reads=TRUE, clonal.reads.to.keep=3, 
verbose.console=TRUE, overwrite.wigs=FALSE, keep.wigs=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.IP}{
Path to the ChIP-Seq file.
}
  \item{file.Control}{
Path to the Control file.
}
  \item{format}{
Format of the ChIP-Seq file. By default \code{"bed"}. Other accepted formats: \code{"sam"}, \code{"bam"}, \code{"bowtie"}, \code{"soap"}.
}
  \item{genome}{
Chatacter specifying the genome \code{"hg18"}, \code{"hg19"}, \code{"mm8"}, \code{"mm9"}, \code{"dm2"}, \code{"dm3"}, \code{"ce4"}, \code{"ce6"}, \code{"rn3"}, \code{"rn4"}, \code{"danRer4"}, \code{"danRer6"}.
}
  \item{output.dir}{
Directory to the folder where all the output should be saved. If such a path does not exist, it will be created. By default \code{output.dir = NULL} and a folder named "MeDiChISeq_output" is created in the current working directory. 
}
  \item{output.name}{
Character. Experiment name used when output files are generated. By default \code{output.name = NULL}.
}
  \item{chrom.list}{
Character vector of chromosomes to be deconvolved. If \code{ chrom.list = NULL}, whole genome is deconvolved. 
}
  \item{limL}{
Lowest boundary of the region to deconvolve, by default \code{limL = 0}
}
  \item{limU}{
Upper boundary of the region to deconvolve, by default \code{limU = Inf}. Default values for \code{limL} and \code{limU} result in deconvolving entire chromosomes.
}
  \item{potential.peaks}{
Path or \code{data.frame} that should have columns \code{"chromosome"}, \code{"start"}, \code{"end"}. This is a list of regions in which potential peaks can be located.
}
  \item{reads.elong}{
Path or numeric. Reads are elongated by this value in order to create intensity wig files. Should correspond to the length of the sonicated DNA fragments and can be estimated with \code{\link{fit.peak.profile.seq}} function. Should be the same value which was used for reads elongation when kernel was generated.  
}
  \item{kernel}{
Path, \code{data.frame} or \code{matrix} which is one of the outputs of \code{\link{fit.peak.profile.seq}}, describing the peak shape. 
}
  \item{frag.length}{
Path or numeric value which is one of the outputs of \code{\link{fit.peak.profile.seq}}. This value is used to calculate the base of annotated peaks. 
}
  \item{quant.cutoff}{
Intensity threshold that defines the regions of potential peaks. By assigning a numeric value it will be constant for a whole genome, when defined as quantile the threshold for each window will be calculated based on a Poisson approach. Default value \code{quant.cutoff = "q1e-5"} corresponds to 1 - 1e-5 = 0.99999th quantile from Poisson distribution with lambda equal to the average intensity in that window, as if it was equally distributed.
}
  \item{window}{
Numeric in base pairs. Whole genome is divided into overlapping windows of this size in which the deconvolution takes place. When \code{quant.cutoff} is defined as quantile, the Poisson background is estimated based on the information from such windows. By default \code{window=20000}. It is suggested to use bigger window for broader profiles. See vignette.
}
  \item{wig.res}{
Resolution of the wigs that are created.
}
  \item{fit.res}{
Numeric in base pairs. Desired deconvolution resolution.
}
  \item{max.steps}{
Limit the number of LARS steps taken.
}
  \item{selection.method}{
In LARS regression the accurate model is chosen based on AIC or BIC criterion. Default is \code{"bic"}.
}
  \item{post.proc.factor}{
Post-processing filter for combining deconvolution coefficients. When equal 2, means that all the coefficients that are separated by maximum 2*\code{fit.res} will be merged. When equal 0, no merging is performed, but usually this results in profile overfitting. The value \code{post.proc.factor} * \code{fit.res} can be interpreted as minimal allowed distance between two neighbor peaks.
}
  \item{nr.boots}{
Number of bootstrap iterations to perform. If \code{nr.boots = 1}, no bootstraps are performed and no p-values will be generated. Suggested  value \code{nr.boots = 5}. 
}
  \item{local.windows}{
Numeric vector specifying the size of windows used for local assessment of confidence of putative peaks. By default \code{local.windows = c(1000, 2000, 5000)}.
}
  \item{Control.corr.param}{
Numeric value from [0,1] interval. If Control sample is available this value is the percentage of overlap that is necessary to conclude that two peaks from IP and Control are located in the same position. Consequently the confidence of IP peak is reduced. 
}
  \item{nr.cores}{
Numeric. The number of cores to be used for parallel processing. This argument is passed directly to the \code{\link{mclapply}} function.
}
  \item{remove.clonal.reads}{
Logical. Whether the clonal reads should be removed before creating WIG files.
}
  \item{clonal.reads.to.keep}{
Numeric. Number of reads that are accepted to have exactly the same position. Additional reads are removed. 
}
  \item{verbose.console}{
Logical. If \code{verbose.console = TRUE}, than the progress of fitting is displayed in console and saved in \code{output.dir}.
}
  \item{overwrite.wigs}{
Logical. Whether to overwrite WIG files with the existing ones. When \code{FALSE} one can save time on WIG construction.  
}
  \item{keep.wigs}{
Logical. Whether to keep or delete the generated WIG files. 

}
  \item{...}{
}
}
\value{
Generated output consists of the following object:

\itemize{

\item \code{All.coeffs.IP} that contains a list of all annotated peaks with their intensities, global, local and combined p-values and coefficients after control correction,

}

This object is saved in \code{output.dir}. Moreover you can find there a BED file with annotated peaks and their intensities and WIG files which can be loaded to any genome browser for visual verification of obtained results. In order to facilitate the choice of p-values cutoff a plot of intensities versus p-values is generated. The console output of \code{deconv.entire.genome.seq} is also saved. 

}

\references{
\url{http://archive.igbmc.fr/recherche/Prog_FGC/Eq_HGron/bioinfotools.html} \cr\cr

Mendoza-Parra MA, Walia M, Sankar M, Gronemeyer H. Dissecting the retinoid-induced differentiation of F9 embryonal stem cells by integrative genomics. 
Mol Syst Biol. 2011 Oct 11;7:538. doi: 10.1038/msb.2011.73. PubMed PMID: 21988834; PubMed Central PMCID: PMC3261707. 
URL: \url{http://www.nature.com/msb/journal/v7/n1/full/msb201173.html} \cr\cr

Mendoza-Parra MA, Sankar M, Walia M, Gronemeyer H. POLYPHEMUS: R package for comparative analysis of RNA polymerase II ChIP-seq profiles by non-linear normalization. 
Nucleic Acids Res. 2012 Feb;40(4):e30. doi: 10.1093/nar/gkr1205. Epub 2011 Dec 7. PubMed PMID: 22156059; PubMed Central PMCID: PMC3287170. 
URL: \url{http://nar.oxfordjournals.org/content/40/4/e30.long} \cr\cr

David J. Reiss, Marc T. Facciotti, and Nitin S. Baliga. Model-based deconvolution of genome-wide DNA binding. 
Bioinformatics, pages btm592+, 2007. doi: 10.1093/bioinformatics/btm592. 
URL: \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btm592v1}

}
\author{
Malgorzata Nowicka <gosia.nowicka@uzh.ch> and Wouter Van Gool under the direction of Marco Antonio Mendoza Parra <marco@igbmc.fr> and Hinrich Gronemeyer. Adaptation of MeDiChI by David Reiss.
}
\seealso{
\code{\link{chip.deconv.seq}}, \code{\link{fit.peak.profile.seq}}
}
\examples{
\donttest{
# if possible increase the number of cores used for parallelization
nr.cores=2
# deconvolve a small chunk of chromosome 19 from CTCF profile

dir <- system.file("extdata/CTCF", package="MeDiChISeq")
file.IP <- file.path(dir, "GSM646314_GM12878_CTCF_rep1_rcd.bed")
kernel <- file.path(dir, "MeDiChISeq_CTCF_kernel.txt")
reads.elong <- file.path(dir, "MeDiChISeq_CTCF_reads_elong.txt")
frag.length <- file.path(dir, "MeDiChISeq_CTCF_frag_length.txt")

CTCF.deconv <- deconv.entire.genome.seq(file.IP,
output.name="CTCF", genome="hg18", 
quant.cutoff="q1e-5", window=20000,
chrom.list="chr19", limL=0, limU=500000, 
reads.elong=reads.elong, kernel=kernel, frag.length=frag.length,
nr.cores=nr.cores, nr.boots=2,
remove.clonal.reads=FALSE)
}
}
