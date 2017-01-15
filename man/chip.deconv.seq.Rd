\name{chip.deconv.seq}
\alias{chip.deconv.seq}
\alias{print.chip.deconv.seq}
\alias{plot.chip.deconv.seq}
\alias{coef.chip.deconv.seq}
\title{
High-resolution model-based deconvolution of ChIP-seq data in WIG format.
}
\description{
Deconvolves a subset of data on one chromosome. Used to deconvolve an entire data set across multiple chromosomes, see \code{\link{deconv.entire.genome.seq}}.
}
\usage{
chip.deconv.seq(data, center = NA, window = 20000, kernel = NA, 
quant.cutoff = "q1e-5", fit.res = 50, 
max.steps = 100, post.proc.factor = 2, selection.method = "bic", 
verbose = T, trace = F, nr.boots = 1, 
boot.sample.opt = "residual", max.peak = NA, boot.vary.res = F, 
tile.distance = NA, where = NA, 
min.npeaks = 0, max.npeaks = 99999, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
\code{data.frame} or \code{matrix} containing genome coordinates (column 1) and intensities (column 2), or a path to a WIG file.
}
  \item{center}{
Central coordinate of the data subset to deconvolve.
}
  \item{window}{
Size of the data chunk to deconvolve.
}
  \item{kernel}{
\code{data.frame} or \code{matrix} which is one of the outputs of \code{\link{fit.peak.profile.seq}}, describing the peak shape. 
}
  \item{quant.cutoff}{
Intensity or quantile cutoff for data to be processed. Limits the locations of potential sites to those that are above this threshold.
}
  \item{fit.res}{
Numeric for base pairs. Desired deconvolution resolution.
}
  \item{max.steps}{
Limit the number of LARS steps taken.
}
  \item{post.proc.factor}{
Post-processing filter for combining deconvolution coefficients. When equal 2, means that all the coefficients that are separeted by maximum 2*\code{fit.res} will be merged. 
}
  \item{selection.method}{
In LARS regression the accurate model is chosen based on AIC or BIC criterion. Default is \code{"bic"}.
}
  \item{verbose}{
If \code{TRUE}, print out status messages.
}
  \item{trace}{
If \code{TRUE}, print out LARS progress.
}
  \item{nr.boots}{
Number of bootstrap iterations to perform. Attention: in contrast to \code{chip.deconv}, this version does not produce local p-values at all, since they are calculated in different manner by \code{\link{deconv.entire.genome.seq}}. See details. 
}
  \item{boot.sample.opt}{
Bootstrap resampling option. See Details. Default is \code{"residual"}.
}
  \item{max.peak}{
Any coefficient with an intensity above this threshold is set to this value. Default is \code{NA} - no cutoff.
}
  \item{boot.vary.res}{
If \code{TRUE}, vary the resolution during bootstraps. By default \code{FALSE}.
}
  \item{tile.distance}{
Distance in bp which for ChIP-seq data corresponds to the WIGS resolution. If \code{NA} (the default) then this is computed from the data.
}
  \item{where}{
The chromosome of the subset of data to deconvolve. Inactive (\code{NA}) for ChIP-seq data.
}
  \item{min.npeaks}{
Optionally limit the minimum number of coefficients. Default is \code{0} - no lower limit.
}
  \item{max.npeaks}{
Optionally limit the maximum number of coefficients. Default is \code{99999} - no upper limit.
}
  \item{...}{
All the other parameters of original \code{chip.deconv} from \code{MeDiChI} package such as: \code{boot.sample.opt = "residual"}, \code{max.peak = NA}, \code{boot.vary.res = F}, \code{tile.distance = NA}, \code{where = NA}, \code{min.npeaks = 0}, \code{max.npeaks = 99999} and  \code{smooth=FALSE} can be also changed. Nevertheless in \code{\link{deconv.entire.genome.seq}} and \code{\link{fit.peak.profile.seq}} they are set up to the default values.
}
}
\details{
See details of \code{chip.deconv} from \code{MeDiChI} package.
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
\code{\link{deconv.entire.genome.seq}}, \code{\link{fit.peak.profile.seq}}, \code{chip.deconv}
}
\examples{
# one window deconvolution for CTCF

dir <- system.file("extdata/CTCF", package="MeDiChISeq")
data <- read.table(file.path(dir, 
"chr19_IP_GSM646314_GM12878_CTCF_rep1_rcd_res-10_dist-150_both.wig"), skip=2)
kernel <- read.table(file.path(dir, "MeDiChISeq_CTCF_kernel.txt"))

out <- chip.deconv.seq(data = data, center = 461582, window = 20000, 
kernel = kernel, quant.cutoff = 10, fit.res = 50) 

coef(out)

plot(out)


# one window deconvolution for H3K4me3

dir <- system.file("extdata/H3K4me3", package="MeDiChISeq")
data <- read.table(file.path(dir, 
"chr19_IP_GSM646424_Huvec_H3K4me3_rep1_rcd_res-10_dist-150_both.wig"), skip=2)
kernel <- read.table(file.path(dir, "MeDiChISeq_H3K4me3_kernel.txt"))

out <- chip.deconv.seq(data = data, center = 469384, window = 50000, 
kernel = kernel, quant.cutoff = 10, fit.res = 50) 

coef(out)

plot(out)

plot(out, center = 483241, window = 5000)
}

