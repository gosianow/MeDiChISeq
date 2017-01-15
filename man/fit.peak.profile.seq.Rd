\name{fit.peak.profile.seq}
\alias{fit.peak.profile.seq}
\title{
Learn deconvolution kernel - peak shape - from ChIP-Seq data. 
}
\description{
Learn the peak shape, to be further used for the deconvolution via \code{deconv.entire.genome.seq}.
}
\usage{
fit.peak.profile.seq(file.IP, format="bed", genome="hg19", output.dir=NULL, 
output.name=NULL, chrom.fit=NULL, limL=0, limU=Inf, 
reads.elong=150, quant.cutoff="q1e-7", window=20000, mini.window=2000, 
wig.res=10, fit.res=50, reads.length=50, 
n.peaks = 50, n.skip = 20, re.fit=100, max.iter=500, 
selection.method="bic", post.proc.factor=2,
start.pars =c(shape =10, scale = 20), to.be.fit=c("shape", "scale"),
method = "Nelder-Mead", nr.cores=1, 
remove.clonal.reads=TRUE, clonal.reads.to.keep=3, 
write.pdf=TRUE, save.kernel=TRUE, verbose.console=TRUE, 
overwrite.wigs=FALSE, keep.wigs=TRUE, ...)
}

\arguments{
  \item{file.IP}{
Path to the ChIP-Seq file.
}
  \item{format}{
Format of the ChIP-Seq file. By default \code{"bed"}. Other accepted formats are \code{"sam"}, \code{"bam"}, \code{"bowtie"}, \code{"soap"}.
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
  \item{chrom.fit}{
Character specifying the chromosome on which the learning should be performed. If \code{fit.chrom = NULL}, the shortest chromosome is taken.
}
  \item{limL}{
Lowest boundary of the region to deconvolve, by default \code{limL = 0}
}
  \item{limU}{
Upper boundary of the region to deconvolve, by default \code{limU = Inf}. Default values for \code{limL} and \code{limU} result in entire chromosomes  deconvolution.
}
  \item{reads.elong}{
Numeric. Reads are elongated by this value in order to create intensity WIG files. Should correspond to the average length of sonicated DNA fragments. When \code{reads.elong = NULL}, average fragment length is estimated based on forward and reverse profiles. See vignette.
}
  \item{quant.cutoff}{
Intensity threshold that defines the regions of potential peaks. By assigning a numeric value it will be constant for a whole genome, when defined as quantile the threshold for each window will be calculated based on Poisson approach. Default value \code{quant.cutoff = "q1e-7"} corresponds to 1 - 1e-7 = 0.9999999th quantile from Poisson distribution with lambda equal to the average intensity in that window, as if it was equally distributed. 
}
  \item{window}{
Numeric in base pairs. Whole genome is divided into overlapping windows of this size in which the deconvolution takes place. When \code{quant.cutoff} is defined as quantile, the Poisson background is estimated based on the information from such windows. By default \code{window=20000}. It is sugested to use bigger window for broader profiles. See vignette. 
}
  \item{mini.window}{
Numeric in base pairs. Size of the window to be plotted around each of the bright peaks.
}
  \item{wig.res}{
Resolution of the wigs that are created.
}
  \item{fit.res}{
Numeric in base pairs. Desired deconvolution resolution.
}
  \item{reads.length}{
Length of alligned reads. This value is used during the estimation of \code{reads.elong}. See vignette.
}
  \item{n.peaks}{
Number of regiones with biggest peaks to learn from.
}
  \item{n.skip}{
Number of regiones with worst-fitting to skip. Enables filtering out of peaks that do not agree with the majority of isolated peaks.
}
  \item{re.fit}{
Re-run \code{deconv.entire.genome.adj} on input data using current best-fit profile, every \code{re.fit} iterations.
}
  \item{max.iter}{
Maximum number of iterations for \code{"Nelder-Mead", "BFGS", "CG"} methods.  See \code{\link{optim}} function.
}
  \item{selection.method}{
In LARS regression the accurate model is chosen based on AIC or BIC criterion. Default is \code{"bic"}.
}
  \item{post.proc.factor}{
Post-processing filter for combining deconvolution coefficients. When equal 2, means that all the coefficients that are separated by maximum 2*\code{fit.res} will be merged. When equal 0, no merging is performed, but usually this results with profile overfitting. The value \code{post.proc.factor} * \code{fit.res} can be interpreted as minimal allowed distance between two neighbor peaks.
}
  \item{start.pars}{
Starting parameters for model profile. By default \code{c(shape =10, scale = 20)}, which corresponds to the average DNA fragment length of 200bp. 
}
  \item{to.be.fit}{
Names of parameters to be learned.
}
  \item{method}{
Optimization method, by default \code{method = "Nelder-Mead"}. See \code{\link{optim}} function.
}
  \item{nr.cores}{
Numeric. The number of cores to be used for parallel processing. This argument is passed directly to the \code{\link{mclapply}} function.
}
  \item{remove.clonal.reads}{
Logical. If the clonal reads should be removed before creating WIG files. 
}
  \item{clonal.reads.to.keep}{
Numeric. Number of reads that are accepted to have exactly the same position. Additional reads are removed. 
}
  \item{write.pdf}{
Logical. Whether pdf plots of fitted kernels should be generated and saved in the \code{output.dir}.
}
  \item{save.kernel}{
Logical. Whether the files with \code{kernel} and \code{frag.length} should be saved in the \code{output.dir}.
}
  \item{verbose.console}{
Logical. If \code{verbose.console = TRUE}, than the progress of fitting is displayed in console and saved in \code{output.dir}.
}
  \item{overwrite.wigs}{
Logical. Whether to overwrite WIG files with the existing ones. When \code{FALSE} one can save time on WIGs construction.  
}
  \item{keep.wigs}{
Logical. Whether to keep or delete the generated WIG files. 
}
  \item{...}{
}
}
\value{
A list of 3 objects:
\itemize{
\item \code{reads.elong}: see above.
\item \code{kernel}: \code{data.frame} with the final peak profile.
\item \code{frag.length}: estimated fragment length. Base of the peak profile is assumed to be equal to 2*\code{frag.length}. 
}
Moreover, if \code{kernel=TRUE, verbose.console=TRUE, keep.wigs=TRUE}, than in the directory \code{output.dir} all the above outputs will be saved together with generated WIG files, console display of work progress and the PDF files with kernel plots. 


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
\code{\link{chip.deconv.seq}}, \code{\link{deconv.entire.genome.seq}}
}
\examples{
\donttest{
# if possible increase the number of cores used for parallelization
nr.cores=2

# finds kernel for CTCF profile 

dir <- system.file("extdata/CTCF", package="MeDiChISeq")
file.IP <- file.path(dir, "GSM646314_GM12878_CTCF_rep1_rcd.bed")

CTCF.peak.profile <- fit.peak.profile.seq(file.IP, genome="hg18", chrom.fit="chr19",
output.name="CTCF", window=20000, mini.window=2000, quant.cutoff=20,
limL=0, limU=1000000, re.fit=15, max.iter=20,
start.pars =c(shape = 17, scale = 10), nr.cores=nr.cores,
remove.clonal.reads=FALSE)
}
}
