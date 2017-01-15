### R code from vignette source 'MeDiChISeq.Rnw'

###################################################
### code chunk number 1: library
###################################################
library(MeDiChISeq)


###################################################
### code chunk number 2: fit.peak.profile.seq (eval = FALSE)
###################################################
## file.IP <- "GSM646314_GM12878_CTCF_rep1.bed"
## 
## fit.peak.profile.seq(file.IP, format="bed", genome="hg18", output.dir=NULL, 
## output.name="CTCF", chrom.fit="chr19", limL=0, limU=Inf, 
## reads.elong=150, quant.cutoff="q1e-7", window=50000, 
## mini.window=3000, wig.res=10, fit.res=50, reads.length=50, 
## n.peaks = 50, n.skip = 20, re.fit=100, max.iter=500, selection.method="bic", 
## post.proc.factor=2, start.pars = c(shape =10, scale = 30), 
## to.be.fit=c("shape", "scale"), method = "L-BFGS-B", 
## nr.cores=1, remove.clonal.reads=TRUE, clonal.reads.to.keep=3, 
## write.pdf=TRUE, save.kernel=TRUE, verbose.console=TRUE, 
## overwrite.wigs=FALSE,  keep.wigs=TRUE)


###################################################
### code chunk number 3: fit.peak.profile.seq_Short (eval = FALSE)
###################################################
## fit.peak.profile.seq(file.IP, genome="hg18", 
## output.name="CTCF", chrom.fit="chr19",
## quant.cutoff="q1e-7", window=50000, mini.window=3000,
## start.pars = c(shape =10, scale = 30),
## method = "L-BFGS-B")


###################################################
### code chunk number 4: fit.peak.profile.seq (eval = FALSE)
###################################################
## file.IP <- "GSM646424_Huvec_H3K4me3_rep1.bed"
## 
## fit.peak.profile.seq(file.IP, genome="hg18", 
## output.name="H3K4me3", chrom.fit="chr19",
## quant.cutoff="q1e-7", window=50000, mini.window=3000,
## start.pars = c(shape =10, scale = 30),
## method = "L-BFGS-B")


###################################################
### code chunk number 5: deconv.all (eval = FALSE)
###################################################
## file.IP <-  "GSM646314_GM12878_CTCF_rep1.bed"
## file.Control <-  "GSM646332_GM12878_WCE_rep1.bed"
## 
## reads.elong <- "MeDiChISeq_CTCF_reads.elong.txt"
## kernel <- "MeDiChISeq_CTCF_kernel.txt"
## frag.length <- "MeDiChISeq_CTCF_frag_length.txt"
## 
## 
## deconv.entire.genome.seq(file.IP, file.Control=file.Control, format="bed", 
## genome="hg18", output.dir=NULL, output.name="CTCF", 
## chrom.list=NULL, limL=0, limU=Inf, potential.peaks=NULL, 
## reads.elong=150, kernel=kernel, frag.length=frag.length, 
## quant.cutoff="q1e-5", window=20000, wig.res=10, 
## fit.res=50, max.steps=100, selection.method="bic", 
## post.proc.factor=2, nr.boots=5, 
## local.windows=c(1000, 2000, 5000), Control.corr.param=0.01, 
## nr.cores=1, remove.clonal.reads=F, clonal.reads.to.keep=3, 
## verbose.console=TRUE, overwrite.wigs=FALSE, keep.wigs=TRUE)


###################################################
### code chunk number 6: deconv.all.short (eval = FALSE)
###################################################
## deconv.entire.genome.seq(file.IP, file.Control=file.Control, 
## genome="hg18", output.name="CTCF",
## window=20000,
## reads.elong=reads.elong, kernel=kernel, 
## frag.length=frag.length,
## nr.cores=1, remove.clonal.reads=F)


###################################################
### code chunk number 7: deconv.all (eval = FALSE)
###################################################
## file.IP <-  "GSM646424_Huvec_H3K4me3_rep1.bed"
## file.Control <-  "GSM646430_Huvec_WCE_rep1.bed"
## 
## reads.elong <- "MeDiChISeq_H3K4me3_reads.elong.txt"
## kernel <- "MeDiChISeq_H3K4me3_kernel.txt"
## frag.length <- "MeDiChISeq_H3K4me3_frag_length.txt"
## 
## deconv.entire.genome.seq(file.IP, file.Control=file.Control, 
## genome="hg18", output.name="H3K4me3",
## window=50000, 
## reads.elong=reads.elong, kernel=kernel, 
## frag.length=frag.length,
## nr.cores=1, remove.clonal.reads=F)


###################################################
### code chunk number 8: all.coeffs (eval = FALSE)
###################################################
## out.CTCF <- read.table("MeDiChISeq_CTCF_ALL_COEFFS_IP.txt", head=T)


###################################################
### code chunk number 9: all.coeffs
###################################################
dir <- system.file("extdata/CTCF", package="MeDiChISeq")
out.CTCF <- read.table(file.path(dir, "MeDiChISeq_CTCF_ALL_COEFFS_IP_rcd.txt"), head=T)


###################################################
### code chunk number 10: all.coeffs
###################################################
head(out.CTCF)


###################################################
### code chunk number 11: deconv.piece
###################################################
dir <- system.file("extdata/CTCF", package="MeDiChISeq")
data <- read.table(file.path(dir, "chr19_IP_GSM646314_GM12878_CTCF_rep1_rcd_res-10_dist-150_both.wig"), 
                   skip=2)
kernel <- read.table(file.path(dir, "MeDiChISeq_CTCF_kernel.txt"))


###################################################
### code chunk number 12: deconv.piece (eval = FALSE)
###################################################
## data <- read.table("chr19_IP_GSM646314_GM12878_CTCF_rep1_res-10_dist-150_both.wig", 
##                    skip=2)
## kernel <- read.table("MeDiChISeq_CTCF_kernel.txt")


###################################################
### code chunk number 13: deconv.piece
###################################################
out <- chip.deconv.seq(data = data, center = 461582, window = 20000, 
kernel = kernel, quant.cutoff = "q1e-5", fit.res = 50) 

coef(out)


###################################################
### code chunk number 14: deconv.piece.fig (eval = FALSE)
###################################################
## png( "pics/medichi_deconv_CTCF.png")


###################################################
### code chunk number 15: deconv.piece.fig (eval = FALSE)
###################################################
## plot(out)


###################################################
### code chunk number 16: deconv.piece.fig (eval = FALSE)
###################################################
## dev.off()


###################################################
### code chunk number 17: deconv.piece.fig (eval = FALSE)
###################################################
## png( "pics/medichi_deconv_CTCF_zoom.png")


###################################################
### code chunk number 18: deconv.piece.fig (eval = FALSE)
###################################################
## plot(out, center=457891, window=3000)


###################################################
### code chunk number 19: deconv.piece.fig (eval = FALSE)
###################################################
## dev.off()


###################################################
### code chunk number 20: deconv.piece
###################################################
dir <- system.file("extdata/H3K4me3", package="MeDiChISeq")
data <- read.table(file.path(dir, "chr19_IP_GSM646424_Huvec_H3K4me3_rep1_rcd_res-10_dist-150_both.wig"), 
                   skip=2)
kernel <- read.table(file.path(dir, "MeDiChISeq_H3K4me3_kernel.txt"))


###################################################
### code chunk number 21: deconv.piece (eval = FALSE)
###################################################
## data <- read.table("chr19_IP_GSM646424_Huvec_H3K4me3_rep1_res-10_dist-150_both.wig", 
##                    skip=2)
## kernel <- read.table("MeDiChISeq_H3K4me3_kernel.txt")


###################################################
### code chunk number 22: deconv.piece
###################################################
out <- chip.deconv.seq(data = data, center = 469384, window = 50000, 
kernel = kernel, quant.cutoff = "q1e-5", fit.res = 50) 

coef(out)


###################################################
### code chunk number 23: deconv.piece.fig (eval = FALSE)
###################################################
## plot(out)


###################################################
### code chunk number 24: deconv.piece.fig (eval = FALSE)
###################################################
## plot(out, center=483241, window=5000)


###################################################
### code chunk number 25: deconv.piece.fig (eval = FALSE)
###################################################
## png( "pics/medichi_deconv_H3K4me3.png")
## 
## plot(out)
## 
## dev.off()


###################################################
### code chunk number 26: deconv.piece.fig (eval = FALSE)
###################################################
## png( "pics/medichi_deconv_H3K4me3_zoom.png")
## 
## plot(out, center=483241, window=5000)
## 
## dev.off()


###################################################
### code chunk number 27: threshold
###################################################
out.CTCF[out.CTCF[,"chromosome"] == "chr19" & 
out.CTCF[,"position"] <= 461791 & out.CTCF[,"position"] >= 457891, ]


