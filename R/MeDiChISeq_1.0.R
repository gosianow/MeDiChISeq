# library(Rcpp)
# library(parallel)


# .onLoad <-
#   function( libname, pkgname ) { ##.onAttach
#     cat( "Loading MeDiChISeq, version ", VERSION, " (", DATE, ")\n", sep="" )
#   }



VERSION <- "1.0.10"
DATE <- "2013-08-23"

packageStartupMessage("Loading MeDiChISeq, version ", VERSION, " (", DATE, ")\n")

#####################################################################################

# MeDiChI main functions by David Reiss - adjusted for sequencing data

#####################################################################################

chip.deconv.adj <- function (data,  center = NA, window = 20000,  kernel = NA, 
                             fit.res = 50, wig.res=10, max.steps = 100, post.proc.factor = 2,  
                             selection.method = "bic", quant.cutoff = "q1e-5", nr.boots = 1,  
                             verbose = F, trace = F, plot.status = F, potential.good.locks=NULL) 
{
  # potential.good.locks option
  # works for data without rownames
  tile.distance=wig.res
  max.peak = NA
  min.npeaks = 0
  max.npeaks = 99999
  where = NA
  boot.vary.res = F
  boot.sample.opt = "residual"
  type.lars = "lasso"
  boot.max.steps.factor <- 1.2
  matrix.return <- F
  
  #data <- load.chip.data(data, verbose = verbose)
  
  if (is.na(center) || is.na(window)) {
    center <- round(mean(data[, 1]))
    window <- diff(round(range(data[, 1])))
  }
  wind <- round(c(center - window/2, center + window/2))
  w.expand <- round(c(-0.1, 0.1) * window) + wind
  data <- data[data[, 1] >= w.expand[1], , drop = F]
  if (length(data) <= 0 || nrow(data) <= 2) 
    return(NULL)
  data <- data[data[, 1] <= w.expand[2], , drop = F]
  if (length(data) <= 0 || nrow(data) <= 2) 
    return(NULL)
  if (!is.na(where) && where %in% rownames(data)) 
    data <- data[rownames(data) == where, , drop = F]
  if (length(data) <= 0 || nrow(data) <= 2) 
    return(NULL)
  if (any(is.na(data[, 1]))) 
    data <- data[-which(is.na(data[, 1])), ]
  if (length(data) <= 0 || nrow(data) <= 2) 
    return(NULL)
  if (is.na(max.peak)) 
    max.peak <- max(data[, 2], na.rm = T) * 2
  
  # if(smooth)
  # data <- smooth.wigs(data)
  
  x <- data[, 1]
  y <- data[, 2]
  
  rm(data)
  
  if (is.na(max.steps)) {
    max.steps <- round(diff(range(x))/200)
    if (verbose) 
      cat("Using max.steps =", max.steps, "\n")
  }
  
  if(is.null(potential.good.locks)){
    
    if(is.character(quant.cutoff)){
      
      quant.poiss <- 1 - as.numeric(gsub("q", "", quant.cutoff))
      lambda <- max(1, ceiling(sum(y)/(window/wig.res) ))
      quant.cutoff <- qpois(quant.poiss, lambda)
      
    }
    
    #if (verbose) 
    cat("Using", quant.cutoff, "as data cutoff!\n")
    
    # !!!
    if (!any(y >= quant.cutoff)){
      if (verbose) 
        cat("No enrichments\n")
      return(NULL) 
    }
    
    
    
  } else if (verbose) {
    cat("Using\n")
    print(potential.good.locks)
    cat("as potential peak locuses!\n\n")
    
  }
  
  if (nr.boots < 1) 
    nr.boots <- 1
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  orig.x <- orig.orig.x <- x
  orig.y <- orig.orig.y <- y
  #   posns <- round(orig.x - min(orig.x) + 1)
  #   cnts <- as.numeric(orig.y)
  orig.fit.res <- fit.res
  best.step.1 <- NA
  all.out <- list()
  iter.boot <- 1
  
  
  while (iter.boot <= nr.boots) {
    is.boot <- FALSE
    
    
    if (iter.boot > 1) {
      is.boot <- TRUE
      if (verbose) 
        cat("*** BOOTSTRAP ITER:", iter.boot, "***\n")
      
      orig.x <- orig.orig.x
      orig.y <- orig.orig.y
      
      if (!is.na(pmatch("residual", boot.sample.opt))) {
        if (trace) 
          cat("Resampling: residual!\n")
        if (length(all.out) > 0) 
          fit <- all.out[[1]]$fit
        else fit <- matrix(c(0, 0), ncol = 2)
        dens <- density(orig.y, na.rm = T)
        level.subt <- dens$x[which.max(dens$y)]
        resids <- orig.y - level.subt - fit[, 2]
        if ("residual.1" %in% boot.sample.opt) 
          resids <- resids * sample(c(1, -1), length(resids), 
                                    replace = T)
        else resids <- resids * sample(c(-(sqrt(5) - 
                                             1)/2, (sqrt(5) + 1)/2), length(resids), prob = c((sqrt(5) - 
                                                                                                 1)/(2 * sqrt(5)), (sqrt(5) + 1)/(2 * sqrt(5))), 
                                       replace = T)
        orig.y <- resids
        quant.cutoff <- median(orig.y)
      }
      
      #       if (!is.na(boot.vary.res)) {
      #         if (is.numeric(boot.vary.res)) {
      #           fit.res <- sample(round(orig.fit.res - boot.vary.res):round(orig.fit.res + boot.vary.res), 1)
      #         }
      #         else if (is.logical(boot.vary.res) && boot.vary.res == TRUE) {
      #           fit.res <- sample(round(orig.fit.res * 4/5):round(orig.fit.res * 6/5), 1)
      #         }
      #         if (trace && fit.res != orig.fit.res) 
      #           cat("Varying the boot resolution... current fit.res =", fit.res, "\n")
      #       }
      
    }
    
    
    posns <- round(orig.x - min(orig.x) + 1)
    cnts <- as.numeric(orig.y)
    dens <- density(cnts, na.rm = T)
    level.subtract <- dens$x[which.max(dens$y)]
    cnts <- cnts - level.subtract
    
    
    if(!is.null(potential.good.locks) && iter.boot==1){
      
      potential.good.locks <- unique(as.vector(sapply(potential.good.locks, function(x) x + c(-2,-1,1,2)*fit.res)))
      
      good.locs <- unique(sapply(potential.good.locks,function(x)which.min(abs(orig.x-x))))
      
    } else{
      
      good.locs <- which(cnts > quant.cutoff - level.subtract)
      
    }
    
    if (!exists("xxx")) 
      xxx <- NULL
    
    
    # controlling if any intensities are indeed higher than quant.cutoff
    if (length(good.locs) <= 0 || length(unique(posns)) < 3) {
      
      if (iter.boot == 1){
        if(verbose)
          cat("ER1")
        return(NULL)   
      } else {
        out <- list( coeffs = matrix(nrow = 0, ncol = 2, dimnames = list(c(), c("position", "intensity"))) )
        attr(out, "class") <- "chip.deconv.seq"
        all.out[[iter.boot]] <- out
        iter.boot <- iter.boot + 1
        next
      }
      
    }
    
    
    good.posns <- unique(round(sort(posns[good.locs])))
    if (!is.null(matrix.return) && class(matrix.return) %in%  c("matrix", "dgCMatrix")) 
      xxx <- matrix.return
    tmp <- round((-max(tile.distance, fit.res) - 2):(max(tile.distance, fit.res) + 2))
    #tmp <- round((-max(tile.distance, fit.res) - 1):(max(tile.distance, fit.res) + 1)) # orig
    #tmp <- round((-max(tile.distance, fit.res*2) - 1):(max(tile.distance, fit.res*2) + 1))
    
    good.posns.hires <- unique(as.vector(sapply(good.posns, function(i) i + tmp)))
    
    good.posns.hires <- good.posns.hires[good.posns.hires %in% seq(1, diff(range(posns)) + fit.res * 10, by = fit.res)]
    
    if (plot.status) {
      if (iter.boot <= 1) 
        par(mfrow = c(2, 1))
      plot(orig.x, cnts, pch = 19, cex = 0.2)
      points(orig.x[good.locs], cnts[good.locs], pch = 19, cex = 0.5, col = "red")
      points(good.posns.hires + min(orig.x), rep(min(cnts[good.locs],  na.rm = T), length(good.posns.hires)), pch = 19,  cex = 0.3, col = "blue")
    }
    
    ### make.predictor.matrix
    
    if (!exists("xxx") || is.null(xxx) || any(!as.character(good.posns.hires) %in% colnames(xxx)) || any(!as.character(round(posns)) %in% rownames(xxx))) {
      
      tmp.posns <- round(posns)
      tmp.posns <- unique(tmp.posns)
      tmp.good.posns <- good.posns.hires
      
      if (iter.boot > 1 && !is.na(pmatch("position", boot.sample.opt))) 
        xxx <- NULL
      if (exists("xxx") && !is.null(xxx)) {
        tmp.good.posns <- tmp.good.posns[!as.character(tmp.good.posns) %in% colnames(xxx)]
        tmp.posns <- as.numeric(rownames(xxx))
      }
      if (length(tmp.good.posns) > 0) {
        if (is.null(kernel) || is.na(kernel[1])) 
          stop("No kernel provided!\n")
        kernel[, 2] <- kernel[, 2]/max(kernel[, 2])
        
        
        xxx.tmp <- try(make.predictor.matrix(tmp.posns, 
                                             kernel, fit.res = fit.res, good.posns.hires = tmp.good.posns, 
                                             sparse = T, verbose = trace), silent = !trace)
        
        
        if (class(xxx.tmp) == "try-error") {
          if(verbose)
            cat("ER2")
          return(NULL) #stop("Could not allocate predictor matrix!!!")
        }
        colnames(xxx.tmp)[1:length(tmp.good.posns)] <- tmp.good.posns
        rownames(xxx.tmp) <- tmp.posns
        if (exists("xxx") && !is.null(xxx) && nrow(xxx) == nrow(xxx.tmp)) {
          if (exists("cBind")) 
            xxx <- cBind(xxx, xxx.tmp)
          else xxx <- cbind(xxx, xxx.tmp)
        }
        else {
          xxx <- xxx.tmp
        }
      }
    }
    
    # still ???
    xx <- xxx
    if (ncol(xx) != length(good.posns.hires)) {
      if (!all(as.character(round(good.posns.hires)) %in% colnames(xx))) {
        if (trace) 
          cat("ERROR1\n")
        next
      }
      xx <- xx[, as.character(round(good.posns.hires)), drop = F]
    }
    if (nrow(xx) != length(posns)) {
      if (!all(as.character(round(posns)) %in% rownames(xx))) {
        if (trace) 
          cat("ERROR2\n")
        next
      }
      xx <- xx[as.character(round(posns)), , drop = F]
    }
    
    colnames(xx) <- as.character(1:ncol(xx))
    if (!is.na(best.step.1)) 
      max.steps <- max(20, ceiling(best.step.1 * boot.max.steps.factor))
    
    ### 1) LARS
    
    lrs <- try(lars.pos(xx, cnts, type = type.lars, max.steps = max.steps, 
                        use.Gram = F, positive = T, trace = trace), silent = !trace)
    
    
    if (class(lrs)[1] == "try-error") {
      if (trace) 
        cat("ERROR3\n")
      
      if (iter.boot == 1){
        if(verbose)
          cat("ER3")
        return(NULL)   
      } else {
        out <- list( coeffs = matrix(nrow = 0, ncol = 2, dimnames = list(c(), c("position", "intensity"))) )
        attr(out, "class") <- "chip.deconv.seq"
        all.out[[iter.boot]] <- out
        iter.boot <- iter.boot + 1
        next
      }
      
    }
    
    max.steps <- min(max.steps, length(lrs$actions))
    lrs.coeff <- predict(lrs, type = "coeff")
    
    coeff.cutoff.LARS <- 1
    
    n.coeffs <- apply(lrs.coeff$coefficients, 1, function(i) sum(i > coeff.cutoff.LARS)) + 1
    dof <- sapply(1:length(n.coeffs), function(i) max(which(n.coeffs == n.coeffs[i])))
    if (is.character(selection.method)) {
      N.x <- length(cnts)
      if (!is.null(lrs$RSS)) {
        rss <- lrs$RSS
      }
      else {
        rss <- apply(lrs.coeff$coefficients, 1, function(i) sum((xx %*% i - cnts)^2, na.rm = T))
      }
      max.row <- which.min(rss)
      hb.fit <- (xx %*% lrs.coeff$coefficients[max.row, ])[, 1]
      sigma.e.sq <- mean((hb.fit - cnts)^2, na.rm = T)
      bic <- rss/sigma.e.sq + log(N.x) * dof
      aic <- rss/sigma.e.sq + 2 * dof
      if (verbose) 
        cat("Step for min AIC:", which.min(aic), n.coeffs[which.min(aic)], "; BIC:", which.min(bic), n.coeffs[which.min(bic)], "; using:", selection.method, "\n")
      #       if (plot.status) {
      #         par(mfrow = c(2, 1))
      #         plot(rss, typ = "l")
      #         plot(log(aic), typ = "l")
      #         plot(log(bic), typ = "l")
      #       }
      aic[is.na(aic)] <- bic[is.na(bic)] <- Inf
      if (selection.method == "aic.and.bic") 
        best.step <- ceiling(mean(c(which.min(aic), which.min(bic)), na.rm = T))
      else if (selection.method == "bic") 
        best.step <- which.min(bic)
      else if (selection.method == "aic") 
        best.step <- which.min(aic)
      if (iter.boot > 1 && best.step >= max.steps) 
        warning(paste("max.steps is probably too low for model selection by", selection.method))
    }
    if (is.na(best.step) || best.step <= 1 && is.numeric(selection.method)) {
      if (selection.method %in% n.coeffs) 
        best.step <- min(which(n.coeffs == selection.method))
      else best.step <- max(which(n.coeffs <= selection.method)) - 1
    }
    if (best.step < 1) 
      best.step <- 1
    if (n.coeffs[best.step] < min.npeaks + 1 && any(n.coeffs >=  min.npeaks + 1)) 
      best.step <- min(which(n.coeffs >= min.npeaks + 1))
    if (n.coeffs[best.step] > max.npeaks + 1 && any(n.coeffs <=  max.npeaks + 1)) 
      best.step <- max(which(n.coeffs <= max.npeaks + 1))
    if (is.na(best.step) || best.step < 1) 
      best.step <- 1
    if (is.na(best.step.1)) 
      best.step.1 <- best.step
    if (!is.na(best.step.1) && best.step >= max.steps) 
      best.step.1 <- max(20, ceiling(best.step * boot.max.steps.factor))
    
    # coeffs
    coeffs <- rep(0, ncol(xx))
    names(coeffs) <- colnames(lrs.coeff$coefficients)
    coeffs[colnames(lrs.coeff$coefficients)] <- lrs.coeff$coefficients[best.step, , drop = F]
    
    if(verbose)
      cat("After LARS step: Number of coeffs:", sum(coeffs > coeff.cutoff.LARS), "\n")
    
    
    # additional part
    if(sum(coeffs > coeff.cutoff.LARS)==0){
      
      if (iter.boot == 1){
        return(NULL)   
      } else {
        out <- list( coeffs = matrix(nrow = 0, ncol = 2, dimnames = list(c(), c("position", "intensity"))) )
        attr(out, "class") <- "chip.deconv.seq"
        all.out[[iter.boot]] <- out
        iter.boot <- iter.boot + 1
        next
      }
      
    }else{
      
      good.coeffs <- coeffs[coeffs > coeff.cutoff.LARS]
      good.xx <- good.posns.hires[coeffs > coeff.cutoff.LARS] + min(orig.x)
      tmp.fit <- (xx %*% coeffs)[, 1]
      
      if(verbose)
        print(good.coeffs)
      
      ### 2) post.proc.coeffs
      
      if (!is.na(post.proc.factor)) {
        
        coeff.cutoff <- 1
        
        #         tmp.coeffs <- post.proc.coeffs(cbind(seq(along = coeffs), coeffs), fit.res = 1, factor = post.proc.factor, 
        #                                        max.coef = max.peak, mean.do = F)
        #         
        #         
        #         tmp.coeffs2 <- coeffs * 0
        #         tmp.coeffs2[round(tmp.coeffs[, 1])] <- tmp.coeffs[,2]
        #         coeffs <- tmp.coeffs2
        #         good.coeffs <- coeffs[coeffs > coeff.cutoff]
        #         good.xx <- good.posns.hires[coeffs > coeff.cutoff] + min(orig.x)
        
        
        coeffs.matrix <- cbind(good.posns.hires, coeffs)
        
        tmp.coeffs <- post.proc.coeffs.seq(coeffs.matrix, fit.res = fit.res, factor = post.proc.factor, 
                                           max.coef = max.peak, mean.do = F, coeff.cutoff=0)
        
        
        good.coeffs <- tmp.coeffs[tmp.coeffs[,2] > coeff.cutoff ,2]
        good.xx <- tmp.coeffs[tmp.coeffs[,2] > coeff.cutoff ,1] +  min(orig.x)
        
        coeffs <-  coeffs*0
        coeffs[(coeffs.matrix[,1]+  min(orig.x)) %in% good.xx] <- good.coeffs
        
        
        
        if (verbose){
          cat("After POST.PROC step: Reduced to", sum(coeffs > coeff.cutoff), "non-redundant coeffs.\n")
          print(good.coeffs)
        }
        
      }
      
      final.rss <- sum(cnts^2, na.rm = T)
      
      ### 3) solve.QP
      
      if (best.step > 1 && n.coeffs[best.step] > 0 && sum(coeffs >  coeff.cutoff) > 0) {
        
        tmp.xx <- as.matrix(xx[, coeffs > coeff.cutoff, drop = F])
        Dmat <- t(tmp.xx) %*% tmp.xx
        if (!is.positive.definite(Dmat)) 
          Dmat <- make.positive.definite(Dmat)
        
        coeffs.ls <- try(solve.QP(Dmat, t(t(cnts) %*% tmp.xx), 
                                  Amat = t(diag(ncol(tmp.xx))), bvec = rep(0, ncol(tmp.xx)), 
                                  meq = 0), silent = !trace)
        
        if (class(coeffs.ls) == "try-error") {
          if (trace) 
            cat("ERROR4\n")
          if (iter.boot == 1){
            return(NULL)   
          } else {
            out <- list( coeffs = matrix(nrow = 0, ncol = 2, dimnames = list(c(), c("position", "intensity"))) )
            attr(out, "class") <- "chip.deconv.seq"
            all.out[[iter.boot]] <- out
            iter.boot <- iter.boot + 1
            next
          }
        }
        
        coeff.cutoff <- 1
        
        good.coeffs <- coeffs.ls$solution 
        names(good.coeffs) <- colnames(tmp.xx) 
        coeffs[coeffs > coeff.cutoff] <- good.coeffs
        
        good.coeffs <- coeffs[coeffs > coeff.cutoff]
        good.xx <- good.posns.hires[coeffs > coeff.cutoff] + min(orig.x)
        
        tmp.xx <- as.matrix(xx[, coeffs > coeff.cutoff, drop = F])
        
        tmp.fit <- tmp.xx %*% good.coeffs
        final.rss <- sum((tmp.fit - cnts)^2, na.rm = T)
        
        if (verbose){
          cat("After SOLVE.QP step: Reduced to", sum(coeffs > coeff.cutoff), "coeffs.\n")
          print(good.coeffs)
        }
      }
      
      if (plot.status) {
        #par(mfrow = c(2, 1))
        plot(orig.x, cnts, pch = 19, cex = 0.5, xlim = wind, ylim = range(c(cnts, tmp.fit, coeffs)), xlab = "Genome coord.", ylab = "Chip intensity")
        lines(orig.x, tmp.fit, col = "red")
        apply(cbind(good.xx, good.coeffs), 1, function(i) lines(rep(i[1], 2), c(0, i[2]), col = "darkgreen"))
        points(good.xx, good.coeffs[1:length(good.xx)], col = "darkgreen", pch = 19, cex = 0.5)
      }
      
      #       out.info <- c(best.step = best.step, n.coeffs = n.coeffs[best.step], 
      #                     n.coeffs.nr = sum(coeffs > coeff.cutoff), dof = dof[best.step], 
      #                     sigma.e.sq = sigma.e.sq, bic = bic[best.step], aic = aic[best.step], 
      #                     wig.res = wig.res, rss = final.rss)
      
      if (iter.boot == 1){
        #         out <- list(data = cbind(orig.x, orig.y - level.subtract), 
        #                     fit = cbind(orig.x, tmp.fit), window = wind, kernel = kernel, 
        #                     coeffs = cbind(position = good.xx, intensity = good.coeffs), 
        #                     out.info = out.info)
        out <- list(fit = cbind(orig.x, tmp.fit),
                    window = wind, 
                    coeffs = cbind(position = good.xx, intensity = good.coeffs))
        attr(out, "class") <- "chip.deconv.seq"
        
      } else {
        out <- list( coeffs = cbind(position = good.xx, intensity = good.coeffs) )
      }
      
      attr(out, "class") <- "chip.deconv.seq"
      all.out[[iter.boot]] <- out
      
      iter.boot <- iter.boot + 1
    }
  }
  
  all.out[[1]]$fit <- NULL
  
  attr(all.out, "class") <- "chip.deconv.seq"
  gc()
  return(invisible(all.out))
}


chip.deconv.seq <-
  function (data, center = NA, window = 20000, kernel = NA, quant.cutoff = "q1e-5", 
            fit.res = 50, max.steps = 100, post.proc.factor = 2, selection.method = "bic",  
            verbose = T, trace = F, nr.boots = 1,
            boot.sample.opt = "residual", max.peak = NA, boot.vary.res = F, 
            tile.distance = NA,  where = NA, min.npeaks = 0, max.npeaks = 99999, ...) 
  {
    
    if (!exists("type.lars")) {
      type.lars = "lasso"
      boot.max.steps.factor <- 1.2
      lars.obj.return <- F
      matrix.return <- F
      shrink <- T
      fit.bg <- NA
      plot.status <- F
      plot.boot <- F
      interp <- T
      ls.final.do <- T
      cluster.nodes <- NA
    }
    in.args <- c(mget(names(formals()), envir = as.environment(-1)), 
                 sapply(as.list(substitute({
                   ...
                 })[-1]), deparse))
    shrink.output <- function(y) {
      y[[1]]$args$matrix.return <- y[[1]]$matrix <- y[[1]]$args$kernel <- NULL
      if (length(y) == 1) 
        return(y)
      for (j in 2:length(y)) {
        y[[j]]$args$matrix.return <- y[[j]]$args$data <- y[[j]]$fit <- y[[j]]$matrix <- NULL
        y[[j]]$args$kernel <- y[[j]]$data <- y[[j]]$kernel <- y[[j]]$all.coeffs <- y[[j]]$out.info <- NULL
        y[[j]]$args <- list(fit.res = y[[j]]$args$fit.res)
      }
      y
    }
    data <- load.chip.data(data, verbose = verbose)
    if (is.na(center) || is.na(window)) {
      center <- round(mean(data[, 1]))
      window <- diff(round(range(data[, 1])))
    }
    wind <- round(c(center - window/2, center + window/2))
    w.expand <- round(c(-0.1, 0.1) * window) + wind
    data <- data[data[, 1] >= w.expand[1], , drop = F]
    if (length(data) <= 0 || nrow(data) <= 2) 
      return(NULL)
    data <- data[data[, 1] <= w.expand[2], , drop = F]
    if (length(data) <= 0 || nrow(data) <= 2) 
      return(NULL)
    if (!is.na(where) && where %in% rownames(data)) 
      data <- data[rownames(data) == where, , drop = F]
    if (length(data) <= 0 || nrow(data) <= 2) 
      return(NULL)
    if (any(is.na(data[, 1]))) 
      data <- data[-which(is.na(data[, 1])), ]
    if (length(data) <= 0 || nrow(data) <= 2) 
      return(NULL)
    if (is.na(max.peak)) 
      max.peak <- max(data[, 2], na.rm = T) * 2
    
    # if(smooth)
    # data <- smooth.wigs(data)
    
    x <- data[, 1]
    y <- data[, 2]
    if (is.na(max.steps)) {
      max.steps <- round(diff(range(x))/200)
      if (verbose) 
        cat("Using max.steps =", max.steps, "\n")
    }
    
    if (is.na(tile.distance)) {
      if (is.null(rownames(data))) 
        rownames(data) <- rep("X1", nrow(data))
      tmp.posns <- sort(data[, 1])
      tile.distances <- numeric()
      for (i in unique(names(tmp.posns))) tile.distances <- c(tile.distances, 
                                                              diff(tmp.posns[names(tmp.posns) == i]))
      tile.distance <- median(tile.distances[tile.distances > 
                                               1])
      #         if (verbose) 
      #             cat("MEAN PROBE SPACING =", tile.distance, "\n")
      rm(tmp.posns, tile.distances)
    }
    
    #     if (is.na(quant.cutoff)) {
    #         quant.cutoff <- 0
    #     }
    #     else if (is.character(quant.cutoff)) {
    #         quant.cutoff <- as.numeric(gsub("q", "", quant.cutoff))
    #         quant.cutoff <- quantile(data[, 2], probs = quant.cutoff, 
    #             na.rm = T)
    #     }
    #     else {
    #         quant.cutoff <- quant.cutoff
    #     }
    
    if(is.character(quant.cutoff)){
      
      quant.poiss <- 1- as.numeric(gsub("q", "", quant.cutoff))
      lambda <- max(1, ceiling(sum(data[, 2])/(window/tile.distance) ))
      quant.cutoff <- qpois(quant.poiss, lambda)
      
    }
    
    if (verbose) 
      cat("Using", quant.cutoff, "as data cutoff!\n")
    if (!any(y >= quant.cutoff)) 
      return(NULL)
    
    if (nr.boots < 1) 
      nr.boots <- 1
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    orig.x <- orig.orig.x <- x
    orig.y <- orig.orig.y <- y
    posns <- round(orig.x - min(orig.x) + 1)
    cnts <- as.numeric(orig.y)
    orig.fit.res <- fit.res
    best.step.1 <- NA
    all.out <- list()
    iter.boot <- 1
    while (iter.boot <= nr.boots) {
      is.boot <- FALSE
      if (iter.boot > 1) {
        is.boot <- TRUE
        if (verbose) 
          cat("*** BOOTSTRAP ITER:", iter.boot, "***\n")
        orig.x <- orig.orig.x
        orig.y <- orig.orig.y
        if (!is.na(pmatch("wild", boot.sample.opt))) {
          if (trace) 
            cat("Resampling: wild!\n")
          if (length(all.out) > 0) 
            fit <- all.out[[1]]$fit
          else fit <- matrix(c(0, 0), ncol = 2)
          dens <- density(orig.y, na.rm = T)
          level.subt <- dens$x[which.max(dens$y)]
          resids <- orig.y - level.subt - fit[, 2]
          if ("wild.1" %in% boot.sample.opt) 
            resids <- resids * sample(c(-(sqrt(5) - 1)/2, 
                                        (sqrt(5) + 1)/2), length(resids), prob = c((sqrt(5) - 
                                                                                      1)/(2 * sqrt(5)), (sqrt(5) + 1)/(2 * sqrt(5))), 
                                      replace = T)
          else resids <- resids * sample(c(1, -1), length(resids), 
                                         replace = T)
          orig.x <- fit[, 1]
          orig.y <- fit[, 2] + resids
        }
        if (!is.na(pmatch("residual", boot.sample.opt))) {
          if (trace) 
            cat("Resampling: residual!\n")
          if (length(all.out) > 0) 
            fit <- all.out[[1]]$fit
          else fit <- matrix(c(0, 0), ncol = 2)
          dens <- density(orig.y, na.rm = T)
          level.subt <- dens$x[which.max(dens$y)]
          resids <- orig.y - level.subt - fit[, 2]
          if ("residual.1" %in% boot.sample.opt) 
            resids <- resids * sample(c(1, -1), length(resids), 
                                      replace = T)
          else resids <- resids * sample(c(-(sqrt(5) - 
                                               1)/2, (sqrt(5) + 1)/2), length(resids), prob = c((sqrt(5) - 
                                                                                                   1)/(2 * sqrt(5)), (sqrt(5) + 1)/(2 * sqrt(5))), 
                                         replace = T)
          orig.y <- resids
          quant.cutoff <- median(orig.y)
        }
        if (!is.na(pmatch("replicate", boot.sample.opt))) {
          if (trace) 
            cat("Resampling: replicate!\n")
          sds <- tapply(orig.y, orig.x, sd, na.rm = T)
          mns <- tapply(orig.y, orig.x, mean, na.rm = T)
          unq <- unique(orig.x)
          new.data <- t(apply(cbind(orig.x, orig.y), 1, 
                              function(i) {
                                ind <- which(unq == i[1])
                                c(i[1], rnorm(1, mean = mns[ind], sd = sds[ind]))
                              }))
          orig.x <- new.data[, 1]
          orig.y <- new.data[, 2]
        }
        if (!is.na(pmatch("resample", boot.sample.opt))) {
          if (trace) 
            cat("Resampling: resample!\n")
          orig.y <- sample(orig.y, replace = T)
        }
        if (!is.na(pmatch("case", boot.sample.opt))) {
          if (trace) 
            cat("Resampling: case!\n")
          samp <- sort(sample(1:length(orig.orig.y), replace = T))
          if (!1 %in% samp) 
            samp <- c(1, samp)
          orig.x <- orig.x[samp]
          orig.y <- orig.y[samp]
        }
        if (!is.na(pmatch("position", boot.sample.opt))) {
          if (trace) 
            cat("Resampling: position!\n")
          offsets <- round(rnorm(length(orig.x), mean = 0, 
                                 sd = tile.distance/4))
          sgn <- sign(offsets)
          offsets[abs(offsets) > tile.distance - 5] <- (tile.distance - 
                                                          5) * sgn[abs(offsets) > tile.distance - 5]
          orig.x <- orig.x + offsets
          ord <- order(orig.x)
          orig.y <- orig.y[ord]
          orig.x <- orig.x[ord]
        }
        if (!is.na(boot.vary.res)) {
          if (is.numeric(boot.vary.res)) {
            fit.res <- sample(round(orig.fit.res - boot.vary.res):round(orig.fit.res + 
                                                                          boot.vary.res), 1)
          }
          else if (is.logical(boot.vary.res) && boot.vary.res == 
                     TRUE) {
            fit.res <- sample(round(orig.fit.res * 4/5):round(orig.fit.res * 
                                                                6/5), 1)
          }
          if (trace && fit.res != orig.fit.res) 
            cat("Varying the boot resolution... current fit.res =", 
                fit.res, "\n")
        }
      }
      posns <- round(orig.x - min(orig.x) + 1)
      cnts <- as.numeric(orig.y)
      dens <- density(cnts, na.rm = T)
      level.subtract <- dens$x[which.max(dens$y)]
      cnts <- cnts - level.subtract
      good.locs <- which(cnts > quant.cutoff - level.subtract)
      if (!exists("xxx")) 
        xxx <- NULL
      if (length(good.locs) <= 0 || length(unique(posns)) < 
            3) {
        if (nr.boots <= 1) 
          return(NULL)
        out.info <- c(best.step = 0, n.coeffs = 0, n.coeffs.nr = 0, 
                      dof = 0, sigma.e.sq = NA, bic = NA, aic = NA, 
                      tile.distance = tile.distance, rss = NA)
        out <- list(data = cbind(orig.x, orig.y - level.subtract), 
                    args = in.args, fit = cbind(orig.x, rep(0, length(orig.x))), 
                    window = wind, kernel = kernel, coeffs = matrix(nrow = 0, 
                                                                    ncol = 2, dimnames = list(c(), c("position", 
                                                                                                     "intensity"))), out.info = out.info)
        attr(out, "class") <- "chip.deconv.seq"
        if (!(length(matrix.return) == 1 && matrix.return == 
                FALSE) && exists("xxx")) 
          out$matrix <- xxx
        if (verbose) 
          cat("Number of coeffs:", 0, "...\n")
        if (nr.boots <= 1 || iter.boot == 1) {
          if (nr.boots <= 1) 
            all.out <- out
          else {
            for (i in 1:nr.boots) all.out[[i]] <- out
            if (shrink) 
              all.out <- shrink.output(all.out)
          }
          attr(all.out, "class") <- "chip.deconv.seq"
          return(invisible(all.out))
        }
        else {
          all.out[[length(all.out) + 1]] <- out
          iter.boot <- iter.boot + 1
          next
        }
      }
      try(if (interp && tile.distance > 300 && (!is.boot || 
                                                  is.na(pmatch("position", boot.sample.opt)))) {
        require(zoo, quietly = T, warn.conflicts = F)
        orig.posns <- posns
        orig.cnts <- cnts
        tile.distances <- diff(posns)
        needs.fill <- tile.distances > 1.8 * tile.distance
        if (any(needs.fill)) {
          cnts.filled <- cnts[1]
          posns.filled <- posns[1]
          for (i in 2:length(posns)) {
            if (!needs.fill[i - 1]) {
              posns.filled <- c(posns.filled, posns[i])
              cnts.filled <- c(cnts.filled, cnts[i])
            }
            else {
              posns.filled <- c(posns.filled, rep(NA, round((posns[i] - 
                                                               posns[i - 1])/tile.distance) - 1), posns[i])
              cnts.filled <- c(cnts.filled, rep(NA, round((posns[i] - 
                                                             posns[i - 1])/tile.distance) - 1), cnts[i])
            }
          }
          posns.filled <- na.approx(posns.filled)
          if (tile.distance > 250) 
            cnts.filled[is.na(cnts.filled)] <- 0
          else cnts.filled <- na.approx(cnts.filled, posns.filled)
          posns <- posns.filled
          cnts <- cnts.filled
          orig.x <- posns + min(orig.x) - 1
          orig.y <- cnts + level.subtract
          rm(cnts.filled, posns.filled)
        }
      }, silent = T)
      good.posns <- unique(round(sort(posns[good.locs])))
      if (!is.null(matrix.return) && class(matrix.return) %in% 
            c("matrix", "dgCMatrix")) 
        xxx <- matrix.return
      tmp <- round((-max(tile.distance, fit.res) - 1):(max(tile.distance, 
                                                           fit.res) + 1))
      good.posns.hires <- unique(as.vector(sapply(good.posns, 
                                                  function(i) i + tmp)))
      good.posns.hires <- good.posns.hires[good.posns.hires %in% 
                                             seq(1, diff(range(posns)) + fit.res * 10, by = fit.res)]
      if (plot.status) {
        if (iter.boot <= 1) 
          par(mfrow = c(2, 1))
        plot(orig.x, cnts, pch = 19, cex = 0.2)
        points(orig.x[good.locs], cnts[good.locs], pch = 19, 
               cex = 0.5, col = "red")
        points(good.posns.hires + min(orig.x), rep(min(cnts[good.locs], 
                                                       na.rm = T), length(good.posns.hires)), pch = 19, 
               cex = 0.3, col = "blue")
      }
      if (!exists("xxx") || is.null(xxx) || any(!as.character(good.posns.hires) %in% 
                                                  colnames(xxx)) || any(!as.character(round(posns)) %in% 
                                                                          rownames(xxx))) {
        tmp.posns <- round(posns)
        tmp.posns <- unique(tmp.posns)
        tmp.good.posns <- good.posns.hires
        if (iter.boot > 1 && !is.na(pmatch("position", boot.sample.opt))) 
          xxx <- NULL
        if (exists("xxx") && !is.null(xxx)) {
          tmp.good.posns <- tmp.good.posns[!as.character(tmp.good.posns) %in% 
                                             colnames(xxx)]
          tmp.posns <- as.numeric(rownames(xxx))
        }
        if (length(tmp.good.posns) > 0) {
          if (is.null(kernel) || is.na(kernel[1])) 
            stop("No kernel provided!\n")
          kernel[, 2] <- kernel[, 2]/max(kernel[, 2])
          
          xxx.tmp <- try(make.predictor.matrix(tmp.posns, 
                                               kernel, fit.res = fit.res, good.posns.hires = tmp.good.posns, 
                                               sparse = T, verbose = trace), silent = !trace)
          
          if (class(xxx.tmp) == "try-error") 
            stop("Could not allocate predictor matrix!!!")
          colnames(xxx.tmp)[1:length(tmp.good.posns)] <- tmp.good.posns
          rownames(xxx.tmp) <- tmp.posns
          if (exists("xxx") && !is.null(xxx) && nrow(xxx) == 
                nrow(xxx.tmp)) {
            if (exists("cBind")) 
              xxx <- cBind(xxx, xxx.tmp)
            else xxx <- cbind(xxx, xxx.tmp)
          }
          else {
            xxx <- xxx.tmp
          }
        }
      }
      xx <- xxx
      if (ncol(xx) != length(good.posns.hires)) {
        if (!all(as.character(round(good.posns.hires)) %in% 
                   colnames(xx))) {
          if (trace) 
            cat("ERROR1\n")
          next
        }
        xx <- xx[, as.character(round(good.posns.hires)), 
                 drop = F]
      }
      if (nrow(xx) != length(posns)) {
        if (!all(as.character(round(posns)) %in% rownames(xx))) {
          if (trace) 
            cat("ERROR2\n")
          next
        }
        xx <- xx[as.character(round(posns)), , drop = F]
      }
      colnames(xx) <- as.character(1:ncol(xx))
      if (!is.na(best.step.1)) 
        max.steps <- max(20, ceiling(best.step.1 * boot.max.steps.factor))
      
      lrs <- try(lars.pos(xx, cnts, type = type.lars, max.steps = max.steps, 
                          use.Gram = F, positive = T, trace = trace, ...), 
                 silent = !trace)
      
      if (class(lrs)[1] == "try-error") {
        if (trace) 
          cat("ERROR3\n")
        iter.boot <- iter.boot + 1
        next
      }
      max.steps <- min(max.steps, length(lrs$actions))
      lrs.coeff <- predict(lrs, type = "coeff")
      coeff.cutoff <- 1
      n.coeffs <- apply(lrs.coeff$coefficients, 1, function(i) sum(i > 
                                                                     coeff.cutoff)) + 1
      dof <- sapply(1:length(n.coeffs), function(i) max(which(n.coeffs == 
                                                                n.coeffs[i])))
      if (is.character(selection.method)) {
        N.x <- length(cnts)
        if (!is.null(lrs$RSS)) {
          rss <- lrs$RSS
        }
        else {
          rss <- apply(lrs.coeff$coefficients, 1, function(i) sum((xx %*% 
                                                                     i - cnts)^2, na.rm = T))
        }
        max.row <- which.min(rss)
        hb.fit <- (xx %*% lrs.coeff$coefficients[max.row, 
                                                 ])[, 1]
        sigma.e.sq <- mean((hb.fit - cnts)^2, na.rm = T)
        bic <- rss/sigma.e.sq + log(N.x) * dof
        #aic <- rss/sigma.e.sq + 2 * dof
        aic <- rss/sigma.e.sq + 2 * dof
        if (verbose) 
          cat("Step for min AIC:", which.min(aic), n.coeffs[which.min(aic)], 
              "; BIC:", which.min(bic), n.coeffs[which.min(bic)], 
              "; using:", selection.method, "\n")
        if (plot.status) {
          par(mfrow = c(2, 1))
          plot(rss, typ = "l")
          plot(log(aic), typ = "l")
          plot(log(bic), typ = "l")
        }
        aic[is.na(aic)] <- bic[is.na(bic)] <- Inf
        if (selection.method == "aic.and.bic") 
          best.step <- ceiling(mean(c(which.min(aic), which.min(bic)), 
                                    na.rm = T))
        else if (selection.method == "bic") 
          best.step <- which.min(bic)
        else if (selection.method == "aic") 
          best.step <- which.min(aic)
        if (iter.boot > 1 && best.step >= max.steps) 
          warning(paste("max.steps is probably too low for model selection by", 
                        selection.method))
      }
      if (is.na(best.step) || best.step <= 1 && is.numeric(selection.method)) {
        if (selection.method %in% n.coeffs) 
          best.step <- min(which(n.coeffs == selection.method))
        else best.step <- max(which(n.coeffs <= selection.method)) - 
          1
      }
      if (best.step < 1) 
        best.step <- 1
      if (n.coeffs[best.step] < min.npeaks + 1 && any(n.coeffs >= 
                                                        min.npeaks + 1)) 
        best.step <- min(which(n.coeffs >= min.npeaks + 1))
      if (n.coeffs[best.step] > max.npeaks + 1 && any(n.coeffs <= 
                                                        max.npeaks + 1)) 
        best.step <- max(which(n.coeffs <= max.npeaks + 1))
      if (is.na(best.step) || best.step < 1) 
        best.step <- 1
      if (is.na(best.step.1)) 
        best.step.1 <- best.step
      if (!is.na(best.step.1) && best.step >= max.steps) 
        best.step.1 <- max(20, ceiling(best.step * boot.max.steps.factor))
      coeffs <- rep(0, ncol(xx))
      names(coeffs) <- colnames(lrs.coeff$coefficients)
      coeffs[colnames(lrs.coeff$coefficients)] <- lrs.coeff$coefficients[best.step, , drop = F]
      
      if (verbose) 
        cat("After LARS step: Number of coeffs:", sum(coeffs > coeff.cutoff), "\n")
      #             cat(sum(coeffs > coeff.cutoff), "coeffs at >", coeff.cutoff, 
      #                 "for LARS step", best.step, "\n")
      good.coeffs <- coeffs[coeffs > coeff.cutoff]
      good.xx <- good.posns.hires[coeffs > coeff.cutoff] + min(orig.x)
      #tmp.fit <- (xx %*% coeffs)[, 1]
      
      if (!is.na(post.proc.factor)) {
        #             if (verbose) 
        #                 cat("Number of coeffs:", sum(coeffs > 0), "... ")
        
        #             tmp.coeffs <- post.proc.coeffs(cbind(seq(along = coeffs), 
        #                 coeffs), fit.res = 1, factor = post.proc.factor, 
        #                 max.coef = max.peak, mean.do = F)
        #             
        #             print(tmp.coeffs)
        #             
        #             tmp.coeffs2 <- coeffs * 0
        #             tmp.coeffs2[round(tmp.coeffs[, 1])] <- tmp.coeffs[, 2]
        #             coeffs <- tmp.coeffs2
        #             good.coeffs <- coeffs[coeffs > coeff.cutoff]
        #             good.xx <- good.posns.hires[coeffs > coeff.cutoff] +  min(orig.x)
        # 
        #             
        #             print(coeffs)
        #             
        #             print(good.xx)
        
        
        coeffs.matrix <- cbind(good.posns.hires, coeffs)
        
        #print(coeffs.matrix)
        
        tmp.coeffs <- post.proc.coeffs.seq(coeffs.matrix, fit.res = fit.res, factor = post.proc.factor, 
                                           max.coef = max.peak, mean.do = F, coeff.cutoff=0)
        
        
        good.coeffs <- tmp.coeffs[tmp.coeffs[,2] > coeff.cutoff ,2]
        good.xx <- tmp.coeffs[tmp.coeffs[,2] > coeff.cutoff ,1] +  min(orig.x)
        
        coeffs <-  coeffs*0
        coeffs[(coeffs.matrix[,1]+  min(orig.x)) %in% good.xx] <- good.coeffs
        
        
        #             print(good.coeffs)
        #             print(good.xx)
        #             print(coeffs)
        
        
        if (verbose) 
          cat("After POST.PROC step: Reduced to", sum(coeffs > coeff.cutoff), "non-redundant coeffs.\n")
        #                 cat("Reduced to", sum(coeffs > coeff.cutoff), 
        #                   "non-redundant coeffs.\n")
      }
      out.info <- c(best.step = best.step, n.coeffs = n.coeffs[best.step], 
                    n.coeffs.nr = sum(coeffs > coeff.cutoff), dof = dof[best.step], 
                    sigma.e.sq = sigma.e.sq, bic = bic[best.step], aic = aic[best.step], 
                    tile.distance = tile.distance)
      final.rss <- sum(cnts^2, na.rm = T)
      
      
      if (ls.final.do && best.step > 1 && n.coeffs[best.step] > 0 && sum(coeffs > coeff.cutoff) > 0) {
        
        tmp.xx <- as.matrix(xx[, coeffs > coeff.cutoff, drop = F])
        Dmat <- t(tmp.xx) %*% tmp.xx
        if (!is.positive.definite(Dmat)) 
          Dmat <- make.positive.definite(Dmat)
        
        coeffs.ls <- try(solve.QP(Dmat, t(t(cnts) %*% tmp.xx), 
                                  Amat = t(diag(ncol(tmp.xx))), bvec = rep(0, ncol(tmp.xx)), 
                                  meq = 0), silent = !trace)
        
        if (class(coeffs.ls) == "try-error") {
          if (trace) 
            cat("ERROR4\n")
          iter.boot <- iter.boot + 1
          next
        }
        #             good.coeffs <- coeffs.ls$solution
        #             tmp.fit <- tmp.xx %*% good.coeffs
        #             coeffs[coeffs > coeff.cutoff] <- good.coeffs
        #             names(good.coeffs) <- colnames(tmp.xx)
        #             final.rss <- sum((tmp.fit - cnts)^2, na.rm = T)
        
        good.coeffs <- coeffs.ls$solution 
        names(good.coeffs) <- colnames(tmp.xx) 
        coeffs[coeffs > coeff.cutoff] <- good.coeffs
        
        good.coeffs <- coeffs[coeffs > coeff.cutoff]
        good.xx <- good.posns.hires[coeffs > coeff.cutoff] + min(orig.x)
        
        tmp.xx <- as.matrix(xx[, coeffs > coeff.cutoff, drop = F])
        
        tmp.fit <- tmp.xx %*% good.coeffs
        final.rss <- sum((tmp.fit - cnts)^2, na.rm = T)
        
        if (verbose)
          cat("After SOLVE.QP step: Reduced to", sum(coeffs > coeff.cutoff), "coeffs.\n")
        
        
      }
      
      
      out.info <- c(out.info, rss = final.rss)
      names(out.info) <- c("best.step", "n.coeffs", "n.coeffs.nr", 
                           "dof", "sigma.e.sq", "bic", "aic", "tile.distance", 
                           "rss")
      if (plot.status) {
        par(mfrow = c(2, 1))
        plot(orig.x, cnts, pch = 19, cex = 0.5, xlim = wind, 
             ylim = range(c(cnts, tmp.fit, coeffs)), xlab = "Genome coord.", 
             ylab = "Chip intensity")
        lines(orig.x, tmp.fit, col = "red")
        apply(cbind(good.xx, good.coeffs), 1, function(i) lines(rep(i[1], 
                                                                    2), c(0, i[2]), col = "darkgreen"))
        points(good.xx, good.coeffs[1:length(good.xx)], col = "darkgreen", 
               pch = 19, cex = 0.5)
      }
      out <- list(data = cbind(orig.x, orig.y - level.subtract), 
                  fit = cbind(orig.x, tmp.fit), window = wind, kernel = kernel, 
                  coeffs = cbind(position = good.xx, intensity = good.coeffs), 
                  out.info = out.info)
      attr(out, "class") <- "chip.deconv.seq"
      if (lars.obj.return) 
        out$lars.out <- list(lars.obj = lrs, lars.coeff = lrs.coeff)
      out$args <- in.args
      if (!(length(matrix.return) == 1 && matrix.return == 
              FALSE) && exists("xxx")) 
        out$matrix <- xxx
      try(rm(xx, lrs, lrs.coeff), silent = T)
      all.out[[length(all.out) + 1]] <- out
      if (plot.boot && nr.boots > 1) 
        xqz <- try(plot(out, hi.res = NA, main = paste("BOOT ITER:", 
                                                       iter.boot), ...))
      if (nr.boots <= 1) 
        all.out <- out
      iter.boot <- iter.boot + 1
    }
    #     if (nr.boots > 1) {
    #         if ((!is.na(pmatch("resample", boot.sample.opt)) || !is.na(pmatch("residual", 
    #             boot.sample.opt))) && length(all.out) > 0) {
    #             co <- all.out[[1]]$coeffs
    #             co.p <- apply(co, 1, function(i) sum(sapply(all.out, 
    #                 function(j) any(j$coeffs[, 2] >= i[2])))/length(all.out))
    #             co <- cbind(co, p.value = co.p)
    #             colnames(co)[3] <- "p.value <"
    #             all.out[[1]]$coeffs.w.p.values <- co
    #         }
    #         if (shrink) 
    #             all.out <- shrink.output(all.out)
    #     }
    attr(all.out, "class") <- "chip.deconv.seq"
    if (verbose && trace) 
      print(all.out)
    invisible(all.out)
  }



deconv.entire.genome.adj <- function (data, chroms = NA, window = 20000, step.by = 19000, centers = NULL, 
                                      quiet = F, plot = F, nr.boots = 1, fit.res = 50, quant.cutoff = "q1e-5", post.proc.factor = 3, 
                                      max.peak = NA, kernel = NA, nr.cores=1, verbose=T, ...) 
{
  # multicore changes
  
  if (!exists("fits.fin.only")) {
    fits.fin.only <- T
    in.fits <- NULL
    remove.near.edges <- T
    cluster.nodes <- NA
    save.progress <- NULL
  }
  start.time <- Sys.time()
  created.here <- FALSE
  data <- load.chip.data(data, verbose = verbose)
  fname <- attr(data, "filename")
  posns <- sort(data[, 1])
  if (any(is.na(data[, 1]))) 
    data <- data[-which(is.na(data[, 1])), ]
  if (length(data) <= 0 || nrow(data) <= 2) 
    return(NULL)
  orig.data <- data
  if (is.na(chroms) && !is.null(rownames(orig.data))) 
    chroms <- unique(rownames(orig.data))
  if (is.na(chroms)) 
    chroms <- 1
  chroms <- as.character(chroms)
  if(verbose)
    cat("Running on chromosomes:", chroms, "\n")
  
  
  
  #   if (is.na(quant.cutoff)) {
  #     quant.cutoff <- 0
  #   }
  #   else if (is.character(quant.cutoff)) {
  #     quant.cutoff <- as.numeric(gsub("q", "", quant.cutoff))
  #     quant.cutoff <- quantile(orig.data[, 2], probs = quant.cutoff, 
  #                              na.rm = T)
  #   }
  #   else {
  #     quant.cutoff <- quant.cutoff
  #   }
  
  #   if(is.character(quant.cutoff)){
  #     
  #     quant.poiss <- as.numeric(gsub("q", "", quant.cutoff))
  #     lambda <- max(1, ceiling(sum(orig.data[, 2])/((max(orig.data[,1])-min(orig.data[,1]))/tile.distance) ))
  #     quant.cutoff <- qpois(quant.poiss, lambda)
  #     
  #   }
  
  if(verbose)
    cat("Using", quant.cutoff, "as data cutoff!\n")
  
  if (is.null(kernel) || is.na(kernel[1])) 
    stop("No kernel provided!\n")
  kernel[, 2] <- kernel[, 2]/max(kernel[, 2])
  if (is.na(max.peak)) 
    max.peak <- max(orig.data[, 2]) * 2
  if (is.null(names(posns))) 
    names(posns) <- rep("XXX", length(posns))
  tmp.posns <- sort(posns)
  tile.distances <- numeric()
  for (i in unique(names(tmp.posns))) tile.distances <- c(tile.distances, 
                                                          diff(tmp.posns[names(tmp.posns) == i]))
  tile.distance <- median(tile.distances[tile.distances > 1])
  if(verbose)
    cat("MEAN PROBE SPACING =", tile.distance, "\n")
  rm(tmp.posns, tile.distances)
  fits <- list()
  if (!is.null(in.fits)) 
    fits <- in.fits
  
  
  for (where in chroms) {
    if (!is.null(rownames(orig.data))) 
      data <- orig.data[rownames(orig.data) == where, , 
                        drop = F]
    else data <- orig.data
    rng <- range(data[, 1])
    posns <- round(seq(rng[1] + round(window/2), rng[2], 
                       by = step.by))
    if (!is.null(centers)) {
      if (!is.null(names(centers))) 
        posns <- centers[which(names(centers) == where)]
      else posns <- centers
      posns <- sort(posns)
    }
    if (plot) 
      par(mfrow = c(2, 1))
    
    
    #     apply.func <- lapply
    #     is.parallel <- !no.multicore
    #     if (is.parallel) 
    #       is.parallel <- require(multicore, quietly = T, warn.conflicts = F)
    #     if (is.parallel) 
    #       is.parallel <- !multicore:::isChild()
    #     if (is.parallel) 
    #       is.parallel <- multicore:::detectCores(all.tests = TRUE) > 
    #         1
    #     if (is.parallel) 
    #       apply.func <- mclapply
    #     if (is.parallel && !quiet) 
    #       cat("Parallelizing deconvoluton of", where, "over", 
    #           multicore:::detectCores(all.tests = TRUE), "processor cores.\n")
    #     if (is.parallel) 
    #       warning("WARNING: If you are running on a Windows system, the 'multicore' option will not work.\nPlease re-start with parameter 'no.multicore=TRUE'.\n")
    
    
    out.fits <- mclapply(posns, function(pos) {
      
      tmp.fit <- chip.deconv.seq(data = data, window = window, 
                                 center = pos, where = NA, tile.distance = tile.distance, 
                                 kernel = kernel, max.peak = max.peak, quant.cutoff = quant.cutoff, 
                                 verbose = verbose, fit.res = fit.res, nr.boots = nr.boots, post.proc.factor=post.proc.factor,
                                 ...)
      
      if (is.null(tmp.fit) || length(tmp.fit) == 0 || class(tmp.fit) == 
            "try-error" || (is.list(tmp.fit) && is.character(tmp.fit[[1]]) && 
                              substr(tmp.fit[[1]], 1, 5) == "Error")) {
        if (verbose && !quiet) 
          cat("ERROR running chip.deconv.seq() on", fname, 
              where, which(posns == pos), length(posns), 
              pos, posns[length(posns)], "; DON'T PANIC - probably no data in the window...", 
              "skipping to next.\n")
        return(NULL)
      }
      if (!is.null(centers)) {
        ttmp <- which(posns >= pos - window & posns <= 
                        pos + window)
        if (length(ttmp) > 0) 
          posns <- posns[!1:length(posns) %in% ttmp]
      }
      if (plot) 
        try(plot(tmp.fit, main = paste(fname, "chip", 
                                       where, pos)))
      if (!is.null(tmp.fit$coeffs)) 
        tmp.fit <- list(`1` = tmp.fit)
      ddata <- tmp.fit[[1]]$data
      for (i in 1:length(tmp.fit)) {
        q.fit <- tmp.fit[[i]]
        q.fit$args$matrix.return <- q.fit$args$data <- q.fit$fit <- q.fit$matrix <- q.fit$args$kernel <- NULL
        q.fit$args$where <- where
        if (i > 1) {
          q.fit$data <- q.fit$kernel <- q.fit$all.coeffs <- q.fit$out.info <- NULL
          q.fit$args <- list(fit.res = q.fit$args$fit.res)
        }
        if (nrow(q.fit$coeffs) > 0 && remove.near.edges) {
          coe <- q.fit$coeffs
          tmp <- coe[, 1] <= min(ddata[, 1] + window/20) | 
            coe[, 1] >= max(ddata[, 1] - window/20)
          if (any(tmp)) 
            coe[tmp, 2] <- 0
          coe <- coe[coe[, 2] > 0, , drop = F]
          q.fit$coeffs <- coe
        }
        tmp.fit[[i]] <- q.fit
      }
      if (verbose) 
        cat(fname, where, which(posns == pos), length(posns), 
            pos, posns[length(posns)], "\t", "BEST.STEP =", 
            tmp.fit[[1]]$out.info["best.step"], "COEFFS =", 
            sum(tmp.fit[[1]]$coeffs[, 2] > 0), "\n")
      return(tmp.fit)
    }, mc.cores=nr.cores)
    
    
    
    gc()
    fits[[where]] <- list()
    for (i in 1:length(out.fits)) if (!is.null(out.fits[[i]])) 
      fits[[where]][[length(fits[[where]]) + 1]] <- out.fits[[i]]
  }
  fits.fin <- list()
  for (where in names(fits)) {
    fits.fin[[where]] <- list()
    fts <- fits[[where]]
    if (length(fts) <= 0) 
      next
    fit <- fts[[1]]
    if (!is.null(fit$coeffs)) 
      fit = list(`1` = fit)
    for (j in 1:length(fit)) {
      if (length(fts) > 1) {
        for (i in 2:length(fts)) {
          if (j > length(fts[[i]])) 
            next
          for (n in c("data", "fit", "coeffs", "all.coeffs")) if (!is.null(fts[[i]][[j]][[n]])) 
            fit[[j]][[n]] <- rbind(fit[[j]][[n]], fts[[i]][[j]][[n]])
        }
      }
      for (n in c("data", "fit", "coeffs", "all.coeffs")) {
        if (is.null(fit[[j]][[n]])) 
          next
        if (is.vector(fit[[j]][[n]])) 
          fit[[j]][[n]] <- matrix(fit[[j]][[n]], ncol = 2)
        ord <- order(fit[[j]][[n]][, 1])
        fit[[j]][[n]] <- fit[[j]][[n]][ord, , drop = F]
        fit[[j]][[n]] <- fit[[j]][[n]][!is.na(fit[[j]][[n]][, 
                                                            1]), , drop = F]
      }
      fit[[j]]$window <- range(fit[[j]]$data[, 1])
      if (verbose) 
        cat("\t", where, j, nrow(fit[[j]]$coeffs), "REDUNDANT COEFFS\n")
      attr(fit, "class") <- "chip.deconv.seq"
      fits.fin[[where]][[j]] <- post.proc.deconv(fit[[j]], 
                                                 factor = 6, fit.res = fit.res, max.coef = max.peak, 
                                                 mean.do = T)
      colnames(fits.fin[[where]][[j]]$coeffs) <- c("position", 
                                                   "intensity")
      attr(fits.fin[[where]][[j]], "class") <- "chip.deconv.seq"
      if(verbose)
        cat("\t\t", where, j, nrow(fits.fin[[where]][[j]]$coeffs), 
            "NON-REDUNDANT COEFFS\n")
    }
    attr(fits.fin[[where]], "class") <- "chip.deconv.seq"
  }
  
  
  # p-values
  
  
  if (fits.fin[[1]][[1]]$args$boot.sample.opt %in% c("resample", 
                                                     "residual")) {
    all.coeffs.boot <- do.call(rbind, lapply(fits.fin, function(i) do.call(rbind, 
                                                                           lapply(i[2:length(i)], function(j) rbind(j$coeffs, 
                                                                                                                    c(0, 0))))))
    for (i in 1:length(fits.fin)) {
      if (is.null(fits.fin[[i]]) || length(fits.fin[[i]]) <= 
            0) 
        next
      coeffs <- fits.fin[[i]][[1]]$coeffs
      pvs <- apply(coeffs, 1, function(j) (sum(all.coeffs.boot[, 
                                                               2] >= j[2]) + 1)/(nrow(all.coeffs.boot) + 1))
      tmp <- cbind(coeffs, p.value = pvs)
      colnames(tmp)[3] <- "p.value <"
      fits.fin[[i]][[1]]$coeffs.w.p.values <- tmp
    }
  }
  end.time <- Sys.time()
  if(verbose)
    cat("Completed in", difftime(end.time, start.time, units = "mins"), 
        "minutes.\n")
  out <- list(fits.fin = fits.fin)
  if (!fits.fin.only) 
    out$fits <- fits
  attr(out, "class") <- "chip.deconv.entire.genome"
  invisible(out)
}



generate.binding.profile.adj <- function (fragment.distrib = function(x, ...) dgamma(x, shape = 6, scale = 50), bs.size = 1,
                                          tile.size = 50, min.frag.size = 0, positions = seq(0, 1001, by = 50),
                                          intensity.scaling = function(x, ...) x, 
                                          hybridization.prob = function(x, ...) as.integer(x > 10), interp = T, 
                                          plot = F, verbose = F, nr.cores=1, ...) 
{
  # multicore option fixed
  
  if (!exists("kernel.method")) 
    kernel.method <- "mine"
  if (!0 %in% positions) 
    positions <- c(0, positions)
  positions <- sort(positions)
  max.dist <- max(positions)
  fragment.sizes <- fragment.distrib(positions, ...)
  if (min.frag.size > 1) 
    fragment.sizes[positions < min.frag.size] <- 0
  fragment.sizes[positions < bs.size] <- 0
  fragment.sizes[fragment.sizes < 1e-05] <- 0
  tile.start.end <- c(-ceiling(tile.size/2) + 1, ceiling(tile.size/2))
  bs.half.size <- ceiling(bs.size/2)
  tmp.pos <- 1:length(positions)
  out.distrib <- positions * 0
  
  
  
  #   apply.func <- lapply
  #   is.parallel <- !no.multicore
  #   if (is.parallel) 
  #     is.parallel <- require(multicore, quietly = T, warn.conflicts = F)
  #   if (is.parallel) 
  #     is.parallel <- !multicore:::isChild()
  #   if (is.parallel) 
  #     is.parallel <- multicore:::detectCores(all.tests = TRUE) > 
  #       1
  #   if (is.parallel) 
  #     apply.func <- mclapply
  #   if (is.parallel && verbose) 
  #     cat("Parallelizing generation of binding profile over", 
  #         multicore:::detectCores(all.tests = TRUE), "processor cores.\n")
  #   if (is.parallel) 
  #     warning("WARNING: If you are running on a Windows system, the 'multicore' option will not work.\nPlease re-start with parameter 'no.multicore=TRUE'.\n")
  #   
  
  
  tmp <- mclapply(which(fragment.sizes > 0), function(ind) {
    l <- positions[ind]
    if (verbose) 
      cat(l, "of", max(positions[which(fragment.sizes > 0)]), "\n")
    possible.fragment.locs <- seq.int(-l + bs.half.size + 1, -bs.half.size - 1)
    fragment.covers <- lapply(possible.fragment.locs, function(i) seq.int(i, i + l))
    f.d <- fragment.sizes[ind]
    factor <- max(1, l - 2 * bs.size)
    intens <- intensity.scaling(l, ...)
    for (tile.ind in tmp.pos) {
      tile.pos <- positions[tile.ind]
      tile.s.e <- tile.start.end + tile.pos
      tmp1 <- tile.s.e[1]
      if (tmp1 > l) 
        next
      tmp2 <- tile.s.e[2]
      
      
      tile.overlaps <- sapply(fragment.covers, function(i) sum(i >= tmp1 & i <= tmp2))
      tile.overlaps <- tile.overlaps[tile.overlaps > 0]
      
      
      sum.overlaps <- sum(hybridization.prob(tile.overlaps, ...))
      
      out.distrib[tile.ind] <- out.distrib[tile.ind] + sum.overlaps * intens * f.d/factor
    }
    out.distrib
  }, mc.cores = nr.cores)
  
  
  
  
  for (i in tmp) out.distrib <- out.distrib + i
  out.distrib <- out.distrib/max(out.distrib, na.rm = T)
  out.distrib[is.na(out.distrib)] <- 0
  if (all(is.na(out.distrib) | is.infinite(out.distrib) | all(out.distrib == 
                                                                0))) 
    out.distrib[1:length(out.distrib)] <- 1
  if (interp && !all(min(positions):max(positions) %in% seq(positions))) {
    ss <- predict(smooth.spline(positions, out.distrib, keep.data = F), 
                  min(positions):max(positions))
    out.distrib <- ss$y
    positions <- ss$x
  }
  out.distrib[out.distrib > 1] <- 1
  out.distrib[out.distrib < 0] <- 0
  out.distrib <- c(rev(out.distrib), out.distrib[-1])
  posns <- unique(c(-rev(positions), positions))
  if (plot) 
    plot(posns, out.distrib, typ = "l")
  invisible(cbind(posns, out.distrib))
}



fit.peak.profile.adj <- function (data, tile.size, quant.cutoff="q1e-7", n.peaks = 50, n.skip = 20, in.kernel = NA, fits = NULL, 
                                  method = "Nelder-Mead", positions = c(0, seq(20, 400+ 1, by = 20), seq(450, 1000 + 1, by = 50)), 
                                  re.fit = 100, max.iter=200, start.pars = c(shape = 10, scale = 30, bs.size = 20, h.cutoff = 0), 
                                  to.be.fit = c("shape", "scale", "bs.size", "h.cutoff"), rnd = F, 
                                  mini.window = 2000, plot = T, name = "", nr.cores=1, verbose=T, post.proc.factor=2, fit.res=50, ...) 
{
  
  default.start.pars <- c(shape = 10, scale = 10, bs.size = 20,  h.cutoff = 0, offset = 0)
  default.start.pars[names(start.pars)] <- start.pars
  in.args <- c(mget(names(formals()), envir = as.environment(-1)), 
               sapply(as.list(substitute({... })[-1]), deparse))
  data <- load.chip.data(data, verbose = verbose)
  orig.data <- data
  best.score.so.far <- 9e+09
  iteration <- 1
  best.kernel <- best.params <- NULL
  
  
  #   apply.func <- lapply
  #   is.parallel <- !no.multicore
  #   if (is.parallel) 
  #     is.parallel <- require(multicore, quietly = T, warn.conflicts = F)
  #   if (is.parallel) 
  #     is.parallel <- !multicore:::isChild()
  #   if (is.parallel) 
  #     is.parallel <- multicore:::detectCores(all.tests = TRUE) > 
  #       1
  #   if (is.parallel) 
  #     apply.func <- mclapply
  #   if (is.parallel) 
  #     cat("Parallelizing profile fitting over", multicore:::detectCores(all.tests = TRUE), 
  #         "processor cores.\n")
  #   if (is.parallel) 
  #     warning("WARNING: If you are running on a Windows system, the 'multicore' option will not work.\nPlease re-start with parameter 'no.multicore=TRUE'.\n")
  
  
  get.profile <- function(par, ...) {
    par[names(default.start.pars)[!names(default.start.pars) %in% 
                                    names(par)]] <- default.start.pars[!names(default.start.pars) %in% 
                                                                         names(par)]
    par[c("shape", "scale")] <- abs(par[c("shape", "scale")]) + 
      0.01
    hc <- round(abs(par["h.cutoff"]))
    hc <- min(tile.size - 1, hc)
    bs.size <- round(abs(par["bs.size"]))
    offset <- par["offset"]
    
    
    
    kernel <- generate.binding.profile.adj(tile.size = tile.size, interp = T, bs.size = bs.size, plot = F, 
                                           hybridization.prob = function(x, ...) as.integer(x > hc), 
                                           fragment.distrib = function(x, ...) dgamma(x - offset, shape = par["shape"],  scale = par["scale"]),
                                           positions = positions, verbose = F, nr.cores=nr.cores, ...)
    
    
    
    if (all(is.na(kernel[, 2]) | is.infinite(kernel[, 2]) | 
              kernel[, 2] == 0 | kernel[, 2] == 1)) 
      return(9e+09)
    kernel
  }
  
  
  pks <- datas <- NULL
  changed <- FALSE
  
  get.fit.score <- function(par, pks, kernel = NULL, plot = F, return.all = F)
  {
    
    if (changed && !is.na(re.fit) && iteration%%re.fit == 0) {
      if(verbose)
        cat(iteration, "--> Re-fitting current best profile to all data to find biggest peaks...\n")
      
      
      tmp.fits <- deconv.entire.genome.adj(data = orig.data, 
                                           kernel = best.kernel, plot = F, verbose = F, nr.cores=nr.cores, 
                                           quant.cutoff = quant.cutoff, fit.res=fit.res, post.proc.factor=post.proc.factor, 
                                           ...)
      
      
      ppks <- get.biggest.peaks(tmp.fits$fits.fin, n.peaks = n.peaks)
      if(verbose){
        print(ppks)
      }
      fits <<- fits
      pks <<- ppks
      datas <<- list()
      changed <<- FALSE
      for (i in 1:nrow(pks)) datas[[i]] <<- data[rownames(data) == 
                                                   rownames(pks)[i] & data[, 1] >= pks[i, 1] - mini.window * 
                                                   1.1 & data[, 1] <= pks[i, 1] + mini.window * 
                                                   1.1, , drop = F]
    }
    
    par[names(default.start.pars)[!names(default.start.pars) %in% names(par)]] <- 
      default.start.pars[!names(default.start.pars) %in% names(par)]
    
    
    if(verbose)
      cat(iteration, par, " ")
    
    
    if (is.null(kernel) || class(kernel) != "matrix") 
      kernel <- get.profile(par, ...)
    
    new.fits <- mclapply(1:nrow(pks), function(i) {
      if (i > length(datas)) 
        return(NULL)
      if (is.null(datas[[i]])) 
        return(NULL)
      dat <- datas[[i]]
      
      chip.deconv.seq(data = dat, fit.res = fit.res, window = NA, 
                      max.steps = 50, where = NA, center = NA, 
                      kernel = kernel, verbose = F, trace = F, nr.boots = 1, 
                      tile.distance = tile.size, max.npeaks = 99999, min.npeaks = 1, 
                      quant.cutoff = ifelse(is.character(quant.cutoff), "q1e-3", quant.cutoff) , post.proc.factor=post.proc.factor)
      
    }, mc.cores = nr.cores)
    
    
    pks <- t(sapply(new.fits, function(i) {
      tmp <- i$coeffs[which.max(i$coeffs[, 2]), ]
      return(if (length(tmp) > 0) tmp else rep(NA, 2))
    }))
    
    rownames(pks) <- names(new.fits) <- sapply(new.fits, 
                                               function(i) rownames(i$data)[1])
    rss <- n.data <- rep(NA, length(new.fits))
    for (f in 1:length(new.fits)) if (!is.null(new.fits[[f]])) {
      rss[f] <- new.fits[[f]]$out.info["rss"]
      n.data[f] <- nrow(new.fits[[f]]$data)
    }
    if (is.na(n.skip)) 
      n.skip <- ceiling(length(rss)/6)
    bad.pks <- rep(FALSE, length(rss))
    if (n.skip > 0) 
      bad.pks <- is.na(rss) | is.na(n.data) | n.data < 
      median(n.data, na.rm = T)/2 | rank(rss/n.data/pks[, 
                                                        2]^2, na.last = T) > length(rss) - n.skip
    out.score <- log(sum(rss[!bad.pks]/n.data[!bad.pks], 
                         na.rm = T)) - log(sum(!bad.pks))
    ding <- ""
    if (out.score <= best.score.so.far) {
      best.score.so.far <<- out.score
      best.kernel <<- kernel
      best.params <<- par
      changed <<- TRUE
      ding <- "*"
    }
    if(verbose)
      cat("| N=", sum(!bad.pks), "WORST=", max(rss[!bad.pks]), "|", out.score, ding, "\n")
    
    
    if (plot) {
      obj <- list(par = par, peaks = pks, new.fits = new.fits, 
                  is.bad = bad.pks, kernel = kernel, start = list(kernel = in.kernel, 
                                                                  par = start.pars), score = out.score, args = in.args)
      
      try(plot.fit.peak.profile(obj, n.peak.plot = 7), silent = T)
      
    }
    
    
    iteration <<- iteration + 1
    if (!return.all) 
      return(out.score)
    else return(invisible(list(score = out.score, kernel = kernel, 
                               new.fits = new.fits, is.bad = bad.pks)))
  }
  
  
  input.start.pars <- start.pars
  if (!"shape" %in% names(start.pars)) 
    start.pars["shape"] <- 7
  if (!"scale" %in% names(start.pars)) 
    start.pars["scale"] <- 50
  if (!"bs.size" %in% names(start.pars)) 
    start.pars["bs.size"] <- 20
  if (!"h.cutoff" %in% names(start.pars)) 
    start.pars["h.cutoff"] <- 15
  if (!"offset" %in% names(start.pars)) 
    start.pars["offset"] <- 0
  
  parscale <- start.pars * 0 + 1
  if ("shape" %in% names(parscale)) 
    parscale["shape"] <- 1
  if ("scale" %in% names(parscale)) 
    parscale["scale"] <- 5
  if ("bs.size" %in% names(parscale)) 
    parscale["bs.size"] <- 15
  if ("offset" %in% names(parscale)) 
    parscale["offset"] <- min(10, start.pars["offset"]/2)
  if ("h.cutoff" %in% names(parscale)) 
    parscale["h.cutoff"] <- min(50, start.pars["h.cutoff"]/2)
  if (rnd) 
    start.pars <- start.pars + (runif(length(start.pars)) - 
                                  0.5) * 2 * parscale
  
  if(verbose){
    cat(name, "Starting parameters:\n")
    print(start.pars)
    cat(name, "Allowing only these parameters to be optimized:", to.be.fit, "\n")
  }
  
  if (is.null(in.kernel) || is.na(in.kernel)) 
    in.kernel <- get.profile(start.pars[to.be.fit])
  if (is.null(fits)) 
    
    
    
    fits <- deconv.entire.genome.adj(data = orig.data, kernel = in.kernel, 
                                     plot = F, verbose = F, nr.cores=nr.cores, quant.cutoff=quant.cutoff, fit.res=fit.res, post.proc.factor=post.proc.factor, ...)
  
  
  
  pks <- get.biggest.peaks(fits$fits.fin, n.peaks = n.peaks)
  if(verbose){
    cat("Biggest peaks: \n")
    print(pks)
  }
  max.steps.2 <- 50
  if (!is.null(list(...)$max.steps)) 
    max.steps.2 <- round(list(...)$max.steps/2)
  
  
  #   quant.cutoff <- "q0.98"
  #   q.cutoff <- 0
  #   if (!is.null(list(...)$quant.cutoff)) {
  #     quant.cutoff <- list(...)$quant.cutoff
  #     q.cutoff <- as.numeric(gsub("q", "", quant.cutoff))
  #     q.cutoff <- quantile(data[, 2], probs = q.cutoff, na.rm = T)
  #   }
  
  #   if(is.character(quant.cutoff)){  
  #     quant.poiss <- as.numeric(gsub("q", "", quant.cutoff))
  #     lambda <- max(1, ceiling(sum(data[, 2])/((max(data[,1])-min(data[,1]))/tile.size) ))
  #     quant.cutoff <- qpois(quant.poiss, lambda)
  #     
  #   }
  #   
  #   if(verbose)
  #   cat("Given cutoff =", q.cutoff, ", we have", length(unique(data[data[, 
  #                                                                        2] >= q.cutoff, 1])), "possible peaks. (N.peaks =", n.peaks, 
  #       ")\n")
  
  datas <- list()
  for (i in 1:nrow(pks)) datas[[i]] <- data[rownames(data) == 
                                              rownames(pks)[i] & data[, 1] >= pks[i, 1] - mini.window * 
                                              1.1 & data[, 1] <= pks[i, 1] + mini.window * 1.1, , drop = F]
  
  
  if(verbose)
    cat("Using optimization method:", method, "\n")
  
#   if (method != "GA" && method != "SANN") {
#     
#     
#     out.params <- optim(start.pars[to.be.fit], get.fit.score, 
#                         pks = pks, plot = plot, 
#                         method = method, 
#                         control = list(maxit = max.iter, reltol = 1e-05, parscale = parscale[to.be.fit]),
#                         lower = 1, upper = 300)
#     
#     
#   }
#   else if (method == "SANN") {
#     out.params <- optim(start.pars[to.be.fit], get.fit.score, 
#                         pks = pks, plot = plot, method = "SANN", 
#                         control = list(reltol = 1e-05, temp = 0.01, maxit = max.iter, 
#                                        parscale = parscale[to.be.fit]))
#   }
  
  if (method %in% c("Nelder-Mead", "BFGS", "CG") ){
    out.params <- optim(start.pars[to.be.fit], get.fit.score, 
                        pks = pks, plot = plot, 
                        method = method, 
                        control = list(reltol = 1e-05, maxit = max.iter, parscale = parscale[to.be.fit]))
        
  }
  else if (method %in% c("L-BFGS-B", "Brent")) {
    out.params <- optim(start.pars[to.be.fit], get.fit.score, 
                        pks = pks, plot = plot, 
                        method = method, 
                        control = list(reltol = 1e-05, maxit=200, parscale = parscale[to.be.fit]),
                        lower = 1, upper = 300)
  }
  else if (method == "SANN") {
    out.params <- optim(start.pars[to.be.fit], get.fit.score, 
                        pks = pks, plot = plot, method = "SANN", 
                        control = list(reltol = 1e-05, temp = 0.01, maxit = max.iter, 
                                       parscale = parscale[to.be.fit]))
  }

  
  final.pars <- out.params$par
  iteration <- 1001
  if (!is.na(in.kernel)) {
    if(verbose)
      cat(name, "INPUT KERNEL: ")
    get.fit.score(start.pars, pks, kernel = in.kernel, plot = plot)
  }
  if(verbose)
    cat(name, "START PARAMS: ")
  start.output <- get.fit.score(start.pars, pks, plot = plot, 
                                return.all = T)
  if(verbose)
    cat(name, "FINAL PARAMS: ")
  final.output <- get.fit.score(final.pars, pks, plot = plot, 
                                return.all = T)
  final.output$par <- final.pars
  final.output$fits <- fits
  final.output$peaks <- pks
  start.output$par <- start.pars
  final.output$start <- start.output
  final.output$args <- in.args
  if (plot) {
    try(plot.fit.peak.profile(final.output, n.peak.plot = length(final.output$new.fits)), silent = T)
    try(plot.fit.peak.profile(final.output, n.peak.plot = 0), silent = T)
  }
  attr(final.output, "class") <- "fit.peak.profile"
  invisible(final.output)
}



plot.chip.deconv.seq <- function (x, boot.results = c("scaled.prob", "prob", "scale", "conf=95", "NONE")[1], where = NA, center = NA, window = NULL, 
                                  verbose = F, plot.genes = F, org = NA, hi.res = NA, quants = c(0.95, 0.5, 0.05), smooth = T, ...) 
{
  if (!exists("scale.boot.results")) 
    scale.boot.results <- NA
  in.args <- c(mget(names(formals()), envir = as.environment(-1)), 
               sapply(as.list(substitute({
                 ...
               })[-1]), deparse))
  obj <- x
  is.boot <- is.null(obj$args) && !is.null(obj[[1]]$args) && 
    obj[[1]]$args$nr.boots > 1
  if (!is.boot && !is.null(obj$args)) {
    tmp <- list(obj = obj)
    obj = tmp
  }
  if (is.null(window)) 
    window <- obj[[1]]$window
  if (length(window) == 1) {
    if (is.na(center)) 
      center <- mean(obj[[1]]$window)
    window <- round(c(-window/2, window/2)) + center
  }
  if (is.null(hi.res) || is.na(hi.res)) 
    hi.res <- obj[[1]]$args$fit.res
  fit.bg <- NA
  w.expand <- window + c(-1000, 1000)
  dat <- obj[[1]]$data
  dat <- dat[dat[, 1] >= w.expand[1] & dat[, 1] <= w.expand[2], 
             , drop = F]
  y.range <- range(c(0, dat[, 2]))
  yr <- y.range
  if (plot.genes) 
    y.range[1] <- y.range[1] - diff(y.range)/10
  x.range <- round(range(dat[, 1]))
  coeffs <- NULL
  for (ii in 1:length(obj)) coeffs <- rbind(coeffs, obj[[ii]]$coeffs)
  spline.coeffs <- coeffs[grep("spline", rownames(coeffs), 
                               fixed = T), , drop = F]
  if (nrow(coeffs) > 1) 
    coeffs <- coeffs[!is.na(coeffs[, 1]) & coeffs[, 1] >= 
                       w.expand[1] & coeffs[, 1] <= w.expand[2], , drop = F]
  if (nrow(coeffs) > 1) 
    coeffs <- coeffs[order(coeffs[, 1]), , drop = F]
  coeffs[, 1] <- coeffs[, 1] - x.range[1]
  coe <- obj[[1]]$coeffs
  coe <- coe[coe[, 1] >= w.expand[1] & coe[, 1] <= w.expand[2], 
             , drop = F]
  out.scale <- NULL
  posns <- dat[, 1]
  preds <- t(rep(0, length(posns)))
  if (!is.na(hi.res) && hi.res > 0 && !is.na(quants) && quants != 
        "NONE" && nrow(coeffs) > 0) {
    posns <- seq(1, diff(range(dat[, 1])), by = hi.res)
    xx <- make.predictor.matrix(posns, obj[[1]]$kernel, fit.res = hi.res, 
                                good.posns.hires = unique(round(coeffs[, 1])), sparse = F, 
                                verbose = verbose)
    colnames(xx) <- as.character(unique(round(coeffs[, 1])) + 
                                   x.range[1])
    preds <- matrix(0, nrow = length(obj), ncol = nrow(xx))
    for (i in 1:length(obj)) {
      ccoeffs <- obj[[i]]$coeffs
      if (nrow(ccoeffs) <= 0) 
        next
      ccoeffs <- ccoeffs[ccoeffs[, 1] >= w.expand[1] & 
                           ccoeffs[, 1] <= w.expand[2], , drop = F]
      if (nrow(ccoeffs) <= 0) 
        next
      ccoeffs <- ccoeffs[order(ccoeffs[, 1]), , drop = F]
      tmp.pos <- unique(round(ccoeffs[, 1]))
      whch <- !as.character(tmp.pos) %in% colnames(xx)
      if (any(whch)) {
        tmp.pos[whch] <- tmp.pos[whch] + 1
        whch <- !as.character(tmp.pos) %in% colnames(xx)
        if (any(whch)) {
          tmp.pos <- unique(round(ccoeffs[, 1]))
          tmp.pos[whch] <- tmp.pos[whch] - 2
          whch <- !as.character(tmp.pos) %in% colnames(xx)
        }
        if (any(whch)) {
          tmp.pos <- unique(round(ccoeffs[, 1]))
          tmp.pos <- tmp.pos[!whch]
        }
      }
      whch <- !as.character(tmp.pos) %in% colnames(xx)
      if (any(whch) || !any(!whch)) 
        next
      tmp.coe <- cbind(tmp.pos, ccoeffs[!whch, 2])
      xxx <- xx[, as.character(round(tmp.coe[, 1])), drop = F]
      preds[i, ] <- (xxx %*% tmp.coe[, 2])[, 1]
      if (i > 1 && (obj[[1]]$args$boot.sample.opt %in% 
                      c("resample", "residual"))) {
        preds <- preds[1, , drop = F]
        break
      }
    }
    if (nrow(preds) > 1) {
      preds <- t(sapply(quants, function(q) apply(preds, 
                                                  2, quantile, probs = q, na.rm = T)))
    }
  }
  else {
    if (!is.null(obj[[1]]$fit)) 
      preds <- t(obj[[1]]$fit[, 2])
  }
  boot.res.names <- c(prob = "Probability", scaled.prob = "Scaled prob.", 
                      scale = "Scale", `conf=95` = "95% Confidence")
  y.range <- range(y.range, preds, coe[, 2], 0)
  pch <- in.args$pch
  if (is.null(pch)) 
    pch <- 20
  cex <- in.args$cex
  if (is.null(cex)) 
    cex <- 0.5
  ylim <- in.args$ylim
  if (is.null(ylim)) 
    ylim <- y.range
  old.par <- par(pch = pch, cex = cex, ...)
  on.exit(par(old.par))
  
  plot(dat, type="h", xlim = window, ylim = ylim, xlab = "Genome coord.", 
       ylab = "Chip intensity", col="grey", cex.axis=1.7,  cex.lab=1.7, ...)
  
  apply(coe, 1, function(i) lines(rep(i[1], 2), c(0, i[2]), 
                                  col = "darkgreen", lwd=2))
  points(coe, col = "darkgreen", pch = 20)
  if (!is.null(preds)) 
    for (i in 1:nrow(preds)) lines(posns + x.range[1], preds[i, 
                                                             ], col = "red", lty = i + 1, lwd=2)
  if (is.boot && obj[[1]]$args$boot.sample.opt %in% c("case", 
                                                      "wild", "position", "replicate")) {
    tmp.hi.res <- hi.res
    if (is.null(tmp.hi.res) || is.na(tmp.hi.res)) 
      tmp.hi.res <- obj[[1]]$args$fit.res
    if (is.null(tmp.hi.res) || is.na(tmp.hi.res)) 
      tmp.hi.res <- 10
    for (boot.res in boot.results) {
      if (substr(boot.res, 1, 4) == "conf") {
        coeffs <- NULL
        for (ii in 1:length(obj)) coeffs <- rbind(coeffs, 
                                                  obj[[ii]]$coeffs)
        coeffs <- coeffs[!is.na(coeffs[, 1]) & coeffs[, 
                                                      1] >= w.expand[1] & coeffs[, 1] <= w.expand[2], 
                         , drop = F]
        rng <- range(coeffs[, 1])
        tmp1 <- seq(floor(rng[1]), ceiling(rng[2]), by = 1)
        tmp2 <- rep(0, length(tmp1))
        tmp.lst <- list()
        for (ii in 1:length(obj)) {
          c <- obj[[ii]]$coeffs
          ttmp <- tmp2
          ttmp[tmp1 %in% round(c[, 1])] <- c[, 2]
          tmp.lst[[ii]] <- cbind(tmp1, ttmp)
        }
        coeffs <- do.call(rbind, tmp.lst)
        coeffs <- coeffs[coeffs[, 1] >= w.expand[1] & 
                           coeffs[, 1] <= w.expand[2], , drop = F]
        conf.int <- 0.95 * length(obj)
        if (substr(boot.res, 5, 5) == "=") 
          conf.int <- as.integer(strsplit(boot.res, "=")[[1]][2])/100 * 
          length(obj)
        ccoeffs <- coeffs[coeffs[, 2] > 0, , drop = F]
        tile.dist <- max(obj[[1]]$out.info["tile.distance"], 
                         100)
        for (i in 1:nrow(coe)) {
          coord <- coe[i, 1]
          dists <- abs(ccoeffs[, 1] - coord)
          for (dist in seq(1, tile.dist, by = tmp.hi.res)) {
            frac <- sum(dists <= dist)
            if (frac >= conf.int) 
              break
          }
          if (frac >= conf.int) {
            hts <- ccoeffs[dists < dist, 2]
            hts.low.hi <- quantile(hts, c(1 - conf.int/length(obj), 
                                          conf.int/length(obj)), na.rm = T)
            lines(rep(coord - dist, 2), c(0, hts.low.hi[2]), 
                  col = "green")
            lines(rep(coord + dist, 2), c(0, hts.low.hi[2]), 
                  col = "green")
            if (dist < tile.dist) 
              dist <- tile.dist
            lines(rep(coord, 2), rep(hts.low.hi[2], 2), 
                  col = "green")
            lines(c(coord - dist, coord + dist), rep(hts.low.hi[1], 
                                                     2), col = "green")
            lines(c(coord - dist, coord + dist), rep(hts.low.hi[2], 
                                                     2), col = "green")
          }
        }
      }
      else {
        tmp <- get.chip.boot.probs(obj, boot.results = boot.res, 
                                   window = w.expand, hi.res = tmp.hi.res, smooth = smooth)
        ccoeffs <- tmp[, c(1, 3), drop = F]
        coef.rng <- range(ccoeffs[, 2])
        if (is.na(scale.boot.results)) 
          scale.boot.results <- y.range[2]
        out.scale <- scale.boot.results/max(ccoeffs[, 
                                                    2])
        ccoeffs[, 2] <- ccoeffs[, 2] * out.scale
        lines(rbind(cbind(c(-999, min(ccoeffs[, 1]) - 
                              tmp.hi.res), c(0, 0)), ccoeffs, cbind(c(max(ccoeffs[, 
                                                                                  1]) + tmp.hi.res, 9e+09), c(0, 0))), col = "green")
        lines(ccoeffs, col = "green")
        tmp1 <- seq(0, max(coef.rng), length = 4)
        tmp2 <- seq(0, scale.boot.results, length = length(tmp1))
        axis(4, tmp2, labels = sprintf("%.1f", tmp1), 
             ...)
        m.cex <- 1
        mtext(boot.res.names[boot.res], 4, line = 0, 
              padj = 1.2, outer = F, cex = m.cex)
      }
    }
  }
  if (is.na(where)) 
    where <- obj[[1]]$args$where
  if (plot.genes) 
    plot.genes.in.window(mean(w.expand), where, diff(range(w.expand)), 
                         new.plot = F, yoff = y.range[1] + 0.4, org = org, 
                         yscale = abs(yr[1] - y.range[1])/2, ...)
  if (!is.null(out.scale)) 
    return(invisible(out.scale))
}


plot.fit.peak.profile.adj <- function (x) 
{
  plot.spline = F
  obj <- x
  pks <- obj$peaks
  new.fits <- obj$new.fits
  
  rang <- 1:max(x$args$positions)
  
  par <- abs(obj$par)
  par2 <- abs(obj$start$par)
  plot(rang, dgamma(rang, shape = par["shape"], scale = par["scale"]), 
       typ = "l", xlab = "Fragment length", ylab = "Fragment length distrib.", 
       main = sprintf("FINAL: mean = %.3f\nSTART: mean = %.3f", 
                      par["shape"] * par["scale"], par2["shape"] * 
                        par2["scale"]))
  lines(rang, dgamma(rang, shape = par2["shape"], scale = par2["scale"]), 
        col = "green", lty = 3)
  
  wheres <- unlist(sapply(new.fits, function(i) rownames(i$data)[1]))
  tmp <- list()
  for (w in unique(wheres)) {
    tmp[[w]] <- list(list(data = NULL, coeffs = NULL))
    ind <- 1
    for (ww in which(wheres == w)) {
      if (obj$is.bad[ww]) 
        next
      tmp[[w]][[1]]$data <- rbind(tmp[[w]][[1]]$data, new.fits[[ww]]$data)
      tmp[[w]][[1]]$coeffs <- rbind(tmp[[w]][[1]]$coeffs, 
                                    new.fits[[ww]]$coeffs)
      ind <- ind + 1
    }
  }
  tmp <- get.peak.profile(tmp, n.peaks = nrow(pks), out.to = x$args$mini.window, 
                          filter = F)
  
  plot(tmp, pch = 20, cex = 0.5, main = c(paste(sprintf("%.3f", abs(obj$par)), collapse = " ")), 
       xlab = "Offset (bp)", ylab = "Scaled intensity", ylim=c(0,1.2), xlim=c(-1000, 1000))
  lines(obj$kernel, col = "red", lwd = 3, xlim=c(-1000, 1000))
  if (plot.spline) {
    lines(obj$start$kernel, col = "green", lty = 3)
    ss <- predict(smooth.spline(tmp, keep.data = F), min(tmp[, 
                                                             1]):max(tmp[, 1]))
    lines(ss$x, ss$y, col = "blue", lty = 3)
  }
}



post.proc.coeffs.seq <-
  function (coeffs, fit.res=50, max.coef = NA, factor = 2, mean.do = F, coeff.cutoff=0) 
  {
    
    coeffs.org <- coeffs
    coeffs[,1] <- (coeffs[,1]-1)/fit.res
    coeffs.org2 <- coeffs
    coeffs <- coeffs[coeffs[, 2] > coeff.cutoff, , drop = F]
    
    if (nrow(coeffs) <= 1) 
      return(coeffs.org[coeffs.org[, 2] > coeff.cutoff, , drop = F])
    if (is.na(max.coef)) 
      max.coef <- 100
    coeffs[coeffs[, 2] > max.coef, 2] <- max.coef
    ind <- 1
    old.row <- coeffs[1, ]
    tmp <- c(old.row, ind)
    for (i in 2:nrow(coeffs)) {
      row <- coeffs[i, ]
      
      #if (row[1] - old.row[1] > fit.res * factor || (row[1] - old.row[1] > fit.res * factor * 2 && (row[2] > old.row[2] * 5 || old.row[2] >= row[2] * 5)) ) 
      
      if (row[1] - old.row[1] >  factor ) 
        
        #if (row[1] - old.row[1] >= fit.res * factor || ( row[1] - old.row[1] >= fit.res * (factor+1) && (row[2] < old.row[2] * 2 || old.row[2] < row[2] * 2)) ) 
        
        ind <- ind + 1
      
      
      tmp <- rbind(tmp, c(row, ind))
      old.row <- row
    }
    
    coeffs.out <- NULL
    for (i in unique(tmp[, 3])) {
      rows <- tmp[tmp[, 3] == i, , drop = F]
      
      if (!mean.do) 
        out.row <- c(coeffs.org[which.min(abs(coeffs.org2[,1]-weighted.mean(rows[, 1], rows[, 2]))), 1], 
                     sum(rows[, 2]))
      else out.row <- c(coeffs.org[which.min(abs(coeffs.org2[,1]-weighted.mean(rows[, 1], rows[, 2]))), 1], 
                        mean(rows[, 2]))
      coeffs.out <- rbind(coeffs.out, out.row)
    }
    
    rownames(coeffs.out) <- NULL
    coeffs.out
  }


#####################################################################################

# MedichiSeq functions 

#####################################################################################


# # smooth.wigs <- function(data, wig.res=10, smooth.window=50){

# require(TTR)

# min.pos <- min(data[,1])
# max.pos <- max(data[,1])

# new.pos <- seq(min.pos, max.pos, by=wig.res)

# new.data <- cbind(new.pos, rep(0, length(new.pos)))

# new.data[new.data[,1] %in%  data[,1], 2] <- data[,2]

# n <- max(3, round(smooth.window/wig.res))

# new.data[,2] <- round(EMA(new.data[,2], n=n))

# new.data <- as.matrix(new.data[!is.na(new.data[,2]) & new.data[,2] != 0,])

# return(new.data)

# }



time.formatting <- function(secs){
  
  h <- secs %/% 3600
  
  min <- (secs %% 3600) %/% 60
  
  sec <- round(secs %% 60)
  
  return(paste(h, "h", min, "min", sec, "sec")) 
}


getChromosomes <- function(genome){
  
  # chromosome sizes for Homo Sapiens gene build 18
  hg18 <- c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,146274826,140273252,135374737,134452384,132349534,114142980,106368585,100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432,154913754,57772954)
  
  # add chromosome names to chromosome sizes
  names(hg18) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
  
  # chromosome sizes for Homo Sapiens gene build 19
  hg19 <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753, 81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
  
  # add chromosome names to chromosome sizes
  names(hg19) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
  
  # chromosome lengths for mouse genome assembly 8
  mm8 <- c(197069962,181976762,159872112,155029701,152003063,149525685,145134094,132085098,124000669,129959148,121798632,120463159,120614378,123978870,103492577,98252459,95177420,90736837,61321190,165556469,165556469)
  # add chromosome names to chromosome sizes
  names(mm8) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")
  
  # chromosome lengths for mouse genome assembly 9
  mm9 <- c(197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430,166650296,15902555)
  # add chromosome names to chromosome sizes
  names(mm9) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")
  
  dm2 <- c(27905053,23771897,22407834,22224390,20766785,8724946,2955737,1694122,1281640,396896,359526,88110,19517 )
  names(dm2) <- c("chr3R","chr3L","chr2L","chrX","chr2R","chrU","chr3h","chr2h","chr4","chrYh","chrXh","chr4h","chrM")
  
  dm3 <- c(29004656,27905053,24543557,23011544,22422827,21146708,10049037,3288761,2555491,2517507,1351857,368872,347038,204112,19517)
  names(dm3) <- c("chrUextra","chr3R","chr3L","chr2L","chrX","chr2R","chrU","chr2RHet","chr3LHet","chr3RHet","chr4","chr2LHet","chrYHet","chrXHet","chrM")
  
  ce4 <- c(20919398,17718852,17493784,15279316,15072419,13783681,13794)
  names(ce4) <- c("chrV","chrX","chrIV","chrII","chrI","chrIII","chrM") 
  
  ce6 <- c(20919568,17718854,17493785,15279323,15072421,13783681,13794)
  names(ce6) <- c("chrV","chrX","chrIV","chrII","chrI","chrIII","chrM") 
  
  rn3 <- c(268121971,258222147,187371129,173106704,170969371,160775580,147642806,143082968,129061546,113649943,112220682,111348958,110733352,109774626,97307196,90224819,87800381,87338544,75822765,59223525,55296979,46649226)
  names(rn3) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr14","chr13","chr10","chr15","chr17","chr16","chr11","chr18","chrUn","chr19","chr20","chr12")
  
  rn4 <- c(267910886,258207540,187126005,173096209,171063335,160699376,147636619,143002779,129041809,113440463,112194335,111154910,110718848,109758846,97296363,90238779,87759784,87265094,75822765,59218465,55268282,46782294)
  names(rn4) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr14","chr13","chr10","chr15","chr17","chr16","chr11","chr18","chrUn","chr19","chr20","chr12")
  
  danRer4 <- c(91717235,87691871,84656180,77179095,70589895,69554819,66798501,65489547,64258675,63653707,63411520,61889685,59765243,58719258,57214918,56255777,55712184,54070595,53215897,52342180,51715404,47751166,47249802,46081529,40315040,16596)
  names(danRer4) <- c("chr14","chr7","chr5","chr3","chr1","chr6","chr8","chr16","chr13","chr20","chr17","chr2","chr18","chr12","chr15","chr21","chr9","chr10","chr23","chr11","chr19","chr22","chr4","chr24","chr25","chrM")
  
  danRer6 <- c(76918211,74451498,71658100,61647013,60907308,59305620,58009534,55568185,54736511,52930158,51890894,51884995,50748729,49469313,49271716,48708673,47572505,47237297,46853116,44714728,44116856,43467561,41415389,40403431,38768535)
  names(danRer6) <- c("chr7","chr5","chr4","chr6","chr3","chr1","chr2","chr8","chr9","chr14","chr16","chr20","chr13","chr17","chr18","chr19","chr21","chr15","chr12","chr23","chr11","chr10","chr22","chr24","chr25") 
  
  
  # chromosome sizes for Homo Sapiens gene build 18
  test3 <- c(95177420,90736837,61321190)
  
  # add chromosome names to chromosome sizes
  names(test3) <- c("chr17","chr18","chr19")
  
  
  # chromosome sizes for Homo Sapiens gene build 18
  test2 <- c(90736837,61321190)
  
  # add chromosome names to chromosome sizes
  names(test2) <- c("chr18","chr19")
  
  
  # chromosome sizes for Homo Sapiens gene build 18
  test1 <- c(95177420)
  
  # add chromosome names to chromosome sizes
  names(test1) <- c("chr17")
  
  
  chromosomeVec <- vector(mode="character")
  if(genome=="hg18"){
    chromosomeVec <- hg18  
  }
  
  if(genome=="hg19"){
    chromosomeVec <- hg19  
  }
  
  if(genome=="mm8"){
    chromosomeVec <- mm8  
  }
  
  if(genome=="mm9"){
    chromosomeVec <- mm9  
  }
  
  if(genome=="dm2"){
    chromosomeVec <- dm2  
  }
  
  if(genome=="dm3"){
    chromosomeVec <- dm3  
  }
  
  if(genome=="ce4"){
    chromosomeVec <- ce4  
  }
  
  if(genome=="ce6"){
    chromosomeVec <- ce6  
  }
  
  if(genome=="rn3"){
    chromosomeVec <- rn3  
  }
  
  if(genome=="rn4"){
    chromosomeVec <- rn4  
  }
  
  if(genome=="danRer4"){
    chromosomeVec <- danRer4  
  }
  
  if(genome=="danRer6"){
    chromosomeVec <- danRer6  
  }
  if(genome=="test3"){
    chromosomeVec <- test3
  }
  if(genome=="test2"){
    chromosomeVec <- test2
  }
  if(genome=="test1"){
    chromosomeVec <- test1
  }
  return(chromosomeVec)
}


removeClonalReads <- function(fileIn, readsToKeep, outputDir, format){
  r <- .Call("removeClonalReads", fileIn, readsToKeep, outputDir, format, PACKAGE = "MeDiChISeq")
  return(r[1])
}


remove.clons <- function(file, reads.to.keep, output.dir, format){
  
  extString <- ""
  if(format == "sam"){
    extString = "sam"
  }else if(format == "bam"){
    extString = "bam"
  }else if(format == "bowtie"){
    extString = "map"
  }
  #  else if(format == "eland"){
  #    extString = "map"
  #  }
  else if(format == "soap"){
    extString = "sop"
  }else{
    extString = "bed"
  }
  
  suffixLength <- nchar(extString) + 1
  
  file.name.base <- substr(basename(file),0,(nchar(as.character(basename(file)))-suffixLength))
  
  cleaned.file <- paste(output.dir, "/", file.name.base, "_", reads.to.keep, "-duplicate_cleaned.bed", sep="")
  
  if(file.exists(cleaned.file)){
    
    cat(paste("\nFile ", cleaned.file,  " exist already!!!\n",sep=""))
    
    return(cleaned.file)
    
  }else{
    
    r <- removeClonalReads(fileIn=file, readsToKeep=reads.to.keep, outputDir=output.dir, format)
    return(r)
    
  }
}



writeWigs <- function(chromName, fileName, sampleType, outputDir, windowSize, fragLength, split, genomeBuild, zeros, format){
  .Call( "wigWriter", fileName, sampleType, outputDir, chromName, windowSize, fragLength, split, genomeBuild, zeros, format, PACKAGE = "MeDiChISeq" )
}



wig.file <- function(output.dir, chr, sample.type, file, wig.res, reads.elong, strand, format){
  
  extString <- ""
  if(format == "sam"){
    extString = "sam"
  }else if(format == "bam"){
    extString = "bam"
  }else if(format == "bowtie"){
    extString = "map"
  }
  #  else if(format == "eland"){
  #    extString = "map"
  #  }
  else if(format == "soap"){
    extString = "sop"
  }else{
    extString = "bed"
  }
  
  suffixLength <- nchar(extString) + 1
  
  file.name.base <- substr(basename(file),0,(nchar(as.character(basename(file)))-suffixLength))
  wig.path <- paste(output.dir, "/", chr, "_", sample.type, "_", file.name.base, "_res-", wig.res, "_dist-", reads.elong,"_",strand, ".wig", sep="")
  return(wig.path)
}



write.wigs.parallel <- function(file, output.dir, chromosomes, sample.type,  wig.res, reads.elong, split=FALSE, genome, format, zeros=F,
                                nr.cores=1, overwrite.wigs=FALSE, verbose=T)
{
  
  
  if(split){
    
    wig.file.fw <- wig.file(output.dir, chr=chromosomes[1], sample.type, file, wig.res, reads.elong, strand="fw", format)
    wig.file.rev <- wig.file(output.dir, chr=chromosomes[1], sample.type, file, wig.res, reads.elong, strand="rev", format)
    
    
    if(file.exists(wig.file.fw) && file.exists(wig.file.rev) && !overwrite.wigs){
      if(verbose)
        cat(paste("\nFiles ", wig.file.fw, "   and  ", wig.file.rev, " exist already!!!\n",sep=""))
      
      wigs.vector <- c(wig.file.fw, wig.file.rev)
      
    }else{
      
      wigs.vector <- writeWigs(chromName=chromosomes[1], fileName=file, sampleType=sample.type, outputDir=output.dir, 
                               windowSize=wig.res, fragLength=reads.elong, split=split, genomeBuild=genome, zeros=zeros, format=format)
    }
    
  }else{
    
    
    wigs.vector <- mclapply(chromosomes, function(chr){
      
      wig.file.both <- wig.file(output.dir, chr=chr, sample.type, file, wig.res, reads.elong, strand="both", format)
      
      if(file.exists(wig.file.both) && !overwrite.wigs){
        if(verbose)
          cat(paste("\nFile ", wig.file.both, " exists already!!!\n",sep=""))
        
        
      }else{
        
        wig.file.both <- writeWigs(chromName=chr, fileName=file, sampleType=sample.type, outputDir=output.dir, 
                                   windowSize=wig.res, fragLength=reads.elong, split=split, genomeBuild=genome, zeros=zeros, format=format)
        
        wig.file.both <- wig.file.both[1]
        
      }
      
      return(wig.file.both)
      
    }, mc.cores=nr.cores)
    
    
    wigs.vector <- unlist(wigs.vector)
    names(wigs.vector) <- chromosomes
    
    
  }
  
  return(invisible(wigs.vector))
  
}


countReads <- function(bedIn){  
  r <- .Call("countReads", bedIn);
  return(r[1]) 
}


prepare.data <- function(data, chrom, limL=0, limU=Inf){
  data <- as.matrix(data)
  attr(data,"dimnames") <- list(rep(chrom, nrow(data)), c("positions","intensities"))
  data <- data[data[,1]<=limU & data[,1]>=limL,]
  return(data)
}


load.prepare.wigs.data <- function(path.wigs, chr=NULL, limL=0, limU=Inf, verbose=T, optimal=T){
  
  if(verbose)
    cat("\nLoading and preparing ", path.wigs, " data...\n") 
  
  if(path.wigs=="ILLFORMED"){
    cat("Input file format seems to be ill-formed. No WIG files could be created...\n")
    cat("Please inspect your input file and try again...\n")
    data.wigs <- NULL
    return(data.wigs)
  }
  
  if(!file.exists(path.wigs)){
    cat("WIG file ", path.wigs, "doesn't exist! Skipping...\n")
    data.wigs <- NULL
    return(data.wigs)
  }
  
  data.wigs <- tryCatch({ data.wigs <- read.table(path.wigs, skip=2)}, 
                        error = function(e) {
                          cat("WIGS file is empty! Skipping", chr, " and removing empty file ", path.wigs,"\n")
                          unlink(path.wigs)
                          data.wigs <- NULL  } )
  
  if(is.null(data.wigs))
    return(data.wigs)
  
  data.wigs <- as.matrix(data.wigs)
  if(optimal){
    colnames(data.wigs) <- c("positions", "intensities")
    rownames(data.wigs) <- NULL
  }else{
    attr(data.wigs,"dimnames") <- list(rep(chr, nrow(data.wigs)), c("positions","intensities"))
  }
  data.wigs <- data.wigs[data.wigs[,1]<=limU & data.wigs[,1]>=limL,]
  return(data.wigs)
}



fragment.length <- function(output.fit.peak.profile){ 
  par <- abs(output.fit.peak.profile$par)
  d <- par["shape"] * par["scale"]
  return(as.vector(round(d,digits=0)))
}



cleaning <- function(tmp.fit){ 
  
  if(!is.null(tmp.fit)){ 
    new <- list()
    new[[1]] <- list()
    new[[1]]$coeffs <- tmp.fit[[1]]$coeffs
    
    if(length(tmp.fit)>=2){
      for(i in 2:length(tmp.fit)){
        new[[i]] <- list()
        new[[i]]$coeffs <- tmp.fit[[i]]$coeffs    
      }
    }
    
    gc()
    return(new) 
  }
  return(tmp.fit)
}



parallel.chip.deconv.adj <- function(wigs.IP, wigs.Control, chrom.name, kernel, max.steps, selection.method, post.proc.factor,
                                     fit.res, window, limL, frag.length, potential.peaks, nr.boots, quant.cutoff, nr.cores, verbose=T,
                                     output.dir, output.name, wig.res){
  # performing multicore chip.deconv.adj
  # cleaning of outputs from overlaps between windows 
  # ( I. when default deconvolutions and II. when potential.peaks are defined )
  # saving bootstrap information, if interested in recalculating p-values
  # no original David's p-values!
  # dynamic quant.cutoff per chromosome
  # potential.good.locs option
  
  if.potential.good.locs = T
  
  if(is.null(wigs.IP))
    return(NULL)
  
  if(verbose)
    cat("\nDeconvolving" ,chrom.name, "...\n")
  
  quant.cutoff.IP <- quant.cutoff.Control <- quant.cutoff
  
  #   if (is.character(quant.cutoff)) {
  #     
  #     quant.cutoff <- as.numeric(gsub("q", "", quant.cutoff))
  #     quant.cutoff.IP <- quantile(wigs.IP[,2], probs = quant.cutoff, na.rm = T)
  #     
  #     if (verbose) 
  #       cat("Using", quant.cutoff.IP, "as data cutoff for IP!\n")
  # 
  #     if(!is.null(wigs.Control) && !if.potential.good.locs){
  #       
  #       quant.cutoff.Control <- quantile(wigs.Control[,2], probs = quant.cutoff, na.rm = T)
  #       
  #       if (verbose) 
  #         cat("Using", quant.cutoff.Control, "as data cutoff for Control!\n")   
  #     }
  #   }
  
  
  wind.overlap <- 2*frag.length
  stepSize <- window-wind.overlap
  
  if(is.null(potential.peaks)){
    
    max <- max(wigs.IP[,1])
    # min <- min(wigs.IP[,1])
    min <- limL
    windows.centers <- seq(min+window/2, max, by=stepSize)
    additional.cleaning <- FALSE
    
  }else{    
    
    which.chrom <- potential.peaks[,1]==chrom.name
    windows.centers <- sort(unique(round((potential.peaks[which.chrom,3] + potential.peaks[which.chrom,2])/2)))
    additional.cleaning <- TRUE
    
  }
  
  n.windows.centers <- length(windows.centers)
  
  #write.table(windows.centers, paste(output.dir, "/MeDiChISeq", output.name , "_Centers_",chrom.name,".txt", sep=""), sep="\t", row.names = FALSE, col.names = T, quote = FALSE)
  
  time.chr <- system.time ( out.fits <- mclapply(1:n.windows.centers, function(i){
    #i=21
    
    tmp.fitIP <- NULL
    
    tmp.fitIP <- chip.deconv.adj(data=wigs.IP, 
                                 max.steps = max.steps,
                                 fit.res = fit.res, 
                                 wig.res=wig.res,
                                 nr.boots = nr.boots,
                                 kernel = kernel,
                                 verbose = F, 
                                 window=window, 
                                 center=windows.centers[i], 
                                 quant.cutoff = quant.cutoff.IP,
                                 selection.method=selection.method, 
                                 post.proc.factor=post.proc.factor)
    
    
    if(class(tmp.fitIP)=="chip.deconv.seq"){
      
      nwindL <- tmp.fitIP[[1]]$window[1] + wind.overlap/2
      nwindU <- tmp.fitIP[[1]]$window[2] - wind.overlap/2
      
      coeffs <- coeffs.org <- tmp.fitIP[[1]]$coeffs  
      coeffs <- coeffs[nwindL < coeffs[,"position"] & coeffs[,"position"] <= nwindU, , drop=FALSE]
      
      this.wind.additional.cleaning <- additional.cleaning && i < n.windows.centers && windows.centers[i+1]-windows.centers[i] < stepSize
      
      if(this.wind.additional.cleaning)
        coeffs <- coeffs[coeffs[,"position"] <= (windows.centers[i+1] - window/2 + wind.overlap/2), , drop=FALSE]
      
      
      if(nrow(coeffs)==0){    
        if(verbose)
          cat(chrom.name , ": window ", i," out of", n.windows.centers , " deconvolved!  0  peaks found! \n")   
        gc()
        return(list(IP=NULL, Control=NULL))
      }
      
      
      tmp.fitIP[[1]]$coeffs <- data.frame(chromosome=rep(chrom.name, nrow(coeffs)), start = floor(coeffs[,"position"]- frag.length),
                                          end = ceiling(coeffs[,"position"]+ frag.length), position=coeffs[,"position"], 
                                          intensity=coeffs[,"intensity"], row.names = NULL)
      
      
      if(nr.boots > 1)
        for(b in 2:nr.boots){
          #b=2
          
          coeffs <- tmp.fitIP[[b]]$coeffs
          coeffs <- coeffs[nwindL < coeffs[,"position"] & coeffs[,"position"] <= nwindU, , drop=FALSE]
          
          if(this.wind.additional.cleaning)
            coeffs <- coeffs[coeffs[,"position"] <= (windows.centers[i+1] - window/2 + wind.overlap/2), , drop=FALSE]
          
          if(nrow(coeffs)==0)
            tmp.fitIP[[b]]$coeffs <- NULL else
              tmp.fitIP[[b]]$coeffs <- data.frame(chromosome=rep(chrom.name, nrow(coeffs)), position=coeffs[,"position"], 
                                                  intensity=coeffs[,"intensity"], bootstrap=rep(b, nrow(coeffs)), row.names = NULL)
          
        }
      
      
      tmp.fitControl <- NULL
      
      if(!is.null(wigs.Control)){
        # i=2758
        
        potential.good.locks <- NULL
        if(if.potential.good.locs)
          potential.good.locks <- coeffs.org[,1]
        
        
        tmp.fitControl <- chip.deconv.adj(data=wigs.Control, 
                                          max.steps = max.steps,
                                          fit.res = fit.res,
                                          wig.res=wig.res,
                                          nr.boots = nr.boots,
                                          kernel = kernel,
                                          verbose = F, 
                                          window=window, 
                                          center=windows.centers[i], 
                                          quant.cutoff = quant.cutoff.Control,
                                          potential.good.locks=potential.good.locks,
                                          selection.method=selection.method, 
                                          post.proc.factor=post.proc.factor)
        
        
        if(class(tmp.fitControl)=="chip.deconv.seq"){
          
          coeffs <- tmp.fitControl[[1]]$coeffs
          coeffs <- coeffs[nwindL < coeffs[,"position"] & coeffs[,"position"] <= nwindU, , drop=FALSE]
          
          if(this.wind.additional.cleaning)
            coeffs <- coeffs[coeffs[,"position"] <= (windows.centers[i+1] - window/2 + wind.overlap/2), , drop=FALSE]
          
          if(nrow(coeffs)==0){
            if(verbose)
              cat(chrom.name , ": window ", i," out of", n.windows.centers , " deconvolved! ", nrow(tmp.fitIP[[1]]$coeffs) ," peaks found! \n")   
            gc()
            return(list(IP=cleaning(tmp.fitIP), Control=NULL))
          }
          
          tmp.fitControl[[1]]$coeffs <- data.frame(chromosome=rep(chrom.name, nrow(coeffs)), start = floor(coeffs[,"position"]- frag.length),
                                                   end = ceiling(coeffs[,"position"]+ frag.length), position=coeffs[,"position"], 
                                                   intensity=coeffs[,"intensity"], row.names = NULL) 
          
          if(nr.boots > 1)
            for(b in 2:nr.boots){
              #b=2
              
              coeffs <- tmp.fitControl[[b]]$coeffs
              coeffs <- coeffs[ nwindL < coeffs[,"position"] & coeffs[,"position"] <= nwindU, , drop=FALSE]
              
              if(this.wind.additional.cleaning)
                coeffs <- coeffs[coeffs[,"position"] <= (windows.centers[i+1] - window/2 + wind.overlap/2), , drop=FALSE]
              
              if(nrow(coeffs)==0)
                tmp.fitControl[[b]]$coeffs <- NULL else
                  tmp.fitControl[[b]]$coeffs <- data.frame(chromosome=rep(chrom.name, nrow(coeffs)), position=coeffs[,"position"], 
                                                           intensity=coeffs[,"intensity"], bootstrap=rep(b, nrow(coeffs)), row.names = NULL)
              
            }
          
        }else     
          tmp.fitControl <- NULL
      }
      
    } else
      tmp.fitIP <- tmp.fitControl <- NULL
    
    if(verbose)
      cat(chrom.name , ": window ", i," out of", n.windows.centers , " deconvolved! ", ifelse(is.null(tmp.fitIP), 0, nrow(tmp.fitIP[[1]]$coeffs))," peaks found! \n")   
    gc()
    return(list(IP=cleaning(tmp.fitIP), Control=cleaning(tmp.fitControl)))
    
  }, mc.cores=nr.cores) )
  
  if(verbose)
    cat("\nChromosome ", chrom.name, " done in ", time.chr["elapsed"], " sec!\n")
  
  gc()
  return(out.fits)
}



merging.output <- function(results, output.dir, output.name, nr.boots, if.Control=FALSE){
  # for data.frame
  # merging genome parallel.chip.deconv.adj output
  
  co.IP <- do.call(rbind, lapply(results, function(chr) # for chromosomes
    do.call(rbind, lapply(chr, function(w)# for windows
      rbind(w[["IP"]][[1]]$coeffs, NULL)))))
  
  row.names(co.IP) <- NULL
  
  write.table(co.IP[,c("chromosome", "start", "end", "intensity")], paste(output.dir, "/MeDiChISeq", output.name , "_ALL_COEFFS_", "IP", ".txt", sep=""), sep="\t", row.names = FALSE, col.names = T, quote = FALSE)
  
  write.table(co.IP[,c("chromosome", "start", "end", "intensity")], paste(output.dir, "/MeDiChISeq", output.name , "_intensities_", "IP", ".bed", sep=""), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  cob.IP <- NULL
  
  if(nr.boots>1){
    
    cob.IP <- do.call(rbind, lapply(results, function(chr) # for chromosomes
      do.call(rbind, lapply(chr, function(w)# for windows
        do.call(rbind, lapply(w[["IP"]][2:nr.boots], function(b)
          rbind(b$coeffs, NULL)))))))
    
    row.names(cob.IP) <- NULL
    
    #write.table(cob.IP, paste(output.dir, "/MeDiChISeq", output.name, "_Coeffs_boot_", "IP", ".txt", sep=""), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    
  }
  
  co.Control <- NULL
  cob.Control <- NULL
  
  if(if.Control){
    
    co.Control <- do.call(rbind, lapply(results, function(chr) # for chromosomes
      do.call(rbind, lapply(chr, function(w)# for windows
        rbind(w[["Control"]][[1]]$coeffs, NULL)))))
    
    row.names(co.Control) <- NULL
    
    #write.table(co.Control[,c("chromosome", "start", "end", "intensity")], paste(output.dir, "/MeDiChISeq", output.name , "_ALL_COEFFS_", "Control", ".txt", sep=""), sep="\t", row.names = FALSE, col.names = T, quote = FALSE)
    
    #write.table(co.Control[,c("chromosome", "start", "end", "intensity")], paste(output.dir, "/MeDiChISeq", output.name , "_intensities_", "Control", ".bed", sep=""), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    if(nr.boots>1){
      
      cob.Control <- do.call(rbind, lapply(results, function(chr) # for chromosomes
        do.call(rbind, lapply(chr, function(w)# for windows
          do.call(rbind, lapply(w[["Control"]][2:nr.boots], function(b)
            rbind(b$coeffs, NULL)))))))
      
      row.names(cob.Control) <- NULL
      
      #write.table(cob.Control, paste(output.dir, "/MeDiChISeq", output.name, "_Coeffs_boot_", "Control", ".txt", sep=""), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
    }
    
  }
  
  
  gc()
  return(list(All.coeffs.IP = co.IP, All.coeffs.boot.IP = cob.IP, All.coeffs.Control = co.Control, All.coeffs.boot.Control = cob.Control))
  
}



combining.p.values <- function(pvs){
  
  k<-ncol(pvs)
  chi.sq2 <- -2*rowSums(log(pvs))
  comb.pvs <- pchisq(q=chi.sq2, df=k, lower.tail = FALSE)
  return(comb.pvs)
  
}



local.global.p.value <- function(co, all.coeffs.boot, experiment=c("IP", "Control")[1], local.windows=c(1000, 2000, 5000), output.dir, output.name, 
                                 nr.cores, nr.boots, verbose=T){
  # calculates global p-values
  # calculates local p-values for each of local.windows (using David's approach)
  # no David's p-values!!!
  
  if(!is.null(co)){
    
    if(verbose)
      cat("Global p-values for", experiment, "... ")
    
    n.co <- nrow(co)    
    n.all.coeffs.boot <- nrow(all.coeffs.boot)
    
    gpv.time <- system.time( gpv <- mclapply(1:n.co, function(j) {
      (sum(all.coeffs.boot[, "intensity"] >= co[j,"intensity"]) + 1)/(n.all.coeffs.boot + 1)
    }, mc.cores=nr.cores) )
    
    final.gpv <- do.call("rbind", gpv)
    colnames(final.gpv) <- "global.p.value"
    if(verbose)
      cat("done in" , gpv.time["elapsed"],"sec! \n\n")
    
    
    if(!is.null(local.windows)){
      
      if(verbose)
        cat("Local p-values for", experiment, "... ")
      
      chromosomes <- as.vector(unique(co[, "chromosome"]))
      
      lpv.time <- system.time( lpv.genome <- mclapply(chromosomes, function(chr){
        #chr="chr19"
        
        all.coeffs.boot.chr <- all.coeffs.boot[all.coeffs.boot[, "chromosome"]==chr, , drop=FALSE]  
        co.chr <- co[co[, "chromosome"]==chr, , drop=FALSE]
        
        if(!is.null(co.chr)){
          
          lpv <- matrix(0, nrow(co.chr), length(local.windows))
          
          local.pvs <- paste("local.p.value.", local.windows, sep="")
          colnames(lpv) <- local.pvs  
          
          for(n in 1:nrow(co.chr)){  
            #n=1
            
            for(p in 1:length(local.windows)){
              #p=3
              
              if.overlap.wind <- all.coeffs.boot.chr[, "position"] <= co.chr[n,"position"] + local.windows[p]/2 &   all.coeffs.boot.chr[, "position"] >= co.chr[n,"position"] - local.windows[p]/2
              
              all.coeffs.boot.chr.overlap <- all.coeffs.boot.chr[if.overlap.wind,, drop=FALSE]
              
              lpv[n, p] <- ( sum(sapply(2:nr.boots, function(b) any(all.coeffs.boot.chr.overlap[all.coeffs.boot.chr.overlap[, "bootstrap"]==b, "intensity"] >= co.chr[n,"intensity"]))) + 1) / nr.boots
              
            } 
            
          }
          
          gc()
          return(lpv)
          
        }else
          
          return(NULL)
        
        
      }, mc.cores=nr.cores) )
      if(verbose)
        cat("done in" , lpv.time["elapsed"],"sec! \n\n")
      
      final.lpv.genome <- do.call("rbind", lpv.genome)              
      
      comb.lgpvs <- combining.p.values(cbind(final.lpv.genome, final.gpv))
      
      out <- cbind(co[c("chromosome", "start", "end", "position", "intensity")], final.lpv.genome, final.gpv, combined.local.global.p.values=comb.lgpvs)
      
      if(experiment=="IP"){
        
        write.table(out, paste(output.dir, "/MeDiChISeq", output.name , "_ALL_COEFFS_", experiment, ".txt", sep=""), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        #write.table(out[c("chromosome", "start", "end", "combined.local.global.p.values")], paste(output.dir, "/MeDiChISeq", output.name , "_Combined_local_global_p-values_", experiment, ".bed", sep=""), sep="\t", row.names = FALSE, col.names = F, quote = FALSE)
        
        pdf(paste(output.dir,"/MeDiChISeq", output.name , "_combined_local_global_p-values_VS_intensity_", experiment, ".pdf",sep=""))
        
        plot(log10(out[, "intensity"]), -log10(out[, "combined.local.global.p.values"]), 
             xlab="log10(intensity)", ylab="-log10(combined.local.global.p.values)")
        
        dev.off()
        
      }
      
      return(out)
      
    }
    
    out <- cbind(co[c("chromosome", "start", "end", "position", "intensity")],final.gpv, combined.local.global.p.values=final.gpv)
    
    if(experiment=="IP"){
      
      write.table(out, paste(output.dir, "/MeDiChISeq", output.name , "_ALL_COEFFS_", experiment, ".txt", sep=""), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      #write.table(out[c("chromosome", "start", "end", "combined.local.global.p.values")], paste(output.dir, "/MeDiChISeq", output.name , "_Combined_local_global_p-values_", experiment, ".bed", sep=""), sep="\t", row.names = FALSE, col.names = F, quote = FALSE)
      
      pdf(paste(output.dir,"/MeDiChISeq", output.name , "_combined_local_global_p-values_VS_intensity_", experiment, ".pdf",sep=""))
      
      plot(log10(out[, "intensity"]), -log10(out[, "combined.local.global.p.values"]), 
           xlab="log10(intensity)", ylab="-log10(combined.local.global.p.values)")
      
      dev.off()
      
    }
    
    return(out)
    
    
  }else
    
    return(NULL)
  
}



Control.corrections <- function(P.values.IP, P.values.Control=NULL, overlap=0.01, output.dir, output.name, frag.length, nr.cores, verbose=T)
{
  # v3: tc = -log(pvIP)*IntIP - (-log(pvControl)*IntControl)  
  
  start.time <- proc.time()
  
  if(verbose)
    cat("Control corrections...")
  
  chromosomes <- as.vector(unique(P.values.IP[,"chromosome"]))
  
  corr.coeffs.genome <- mclapply(chromosomes, function(chr){
    #chr="chr19"
    
    coeffsIP <- P.values.IP[P.values.IP[,"chromosome"]==chr, , drop=FALSE]  
    coeffsControl <- P.values.Control[P.values.Control[,"chromosome"]==chr, , drop=FALSE]
    
    corrected.coeffsIP <- cbind(if.control.overlap=rep(0, nrow(coeffsIP)), -log10(coeffsIP[, "combined.local.global.p.values"])* coeffsIP[, "intensity"])
    
    #names.pv <- colnames(coeffsIP[, "combined.local.global.p.values"])
    #names.tc <- paste("tc.", names.pv, sep="")
    names.pv <- "combined.local.global.p.values"
    names.tc <- "control.correction"
    
    colnames(corrected.coeffsIP) <- c("if.control.overlap", names.tc)
    
    
    if(nrow(coeffsControl)!=0 && !is.null(P.values.Control)){
      
      
      for(n in 1:nrow(coeffsIP)){ 
        #n=17
        
        WhichOverlap <- coeffsControl[,"position"] <= coeffsIP[n, "position"]+(1-overlap)*2*frag.length & coeffsControl[,"position"] >= coeffsIP[n, "position"]-(1-overlap)*2*frag.length
        
        if(sum(WhichOverlap)!=0){
          
          corrected.coeffsIP[n, "if.control.overlap"] <- 1
          
          for(p in 1:length(names.pv)){
            #p=1
            
            which.pv <- names.pv[p]
            
            min.pv <- min(coeffsControl[WhichOverlap,which.pv])     
            max.int.pv <- max(coeffsControl[WhichOverlap & coeffsControl[,which.pv]==min.pv, "intensity"])
            position.pv <- coeffsControl[WhichOverlap & coeffsControl[,which.pv]==min.pv & coeffsControl[,"intensity"]==max.int.pv, "position"]
            
            numerator <- -log10(coeffsIP[n,which.pv]) * coeffsIP[n,"intensity"]
            
            denominator <- -log10(min.pv) * max.int.pv
            
            tc <- numerator-denominator
            
            corrected.coeffsIP[n, names.tc[p]] <- tc
            
          }
          
        }
        
      }
      
    }
    
    return(cbind(coeffsIP, corrected.coeffsIP))
    
  }, mc.cores=nr.cores)
  
  
  final.corr.coeffs <- do.call("rbind", corr.coeffs.genome)
  
  end.time <- proc.time()
  
  if(verbose)
    cat("done in", end.time["elapsed"]- start.time["elapsed"], "seconds. \n")
  
  write.table(final.corr.coeffs, paste(output.dir, "/MeDiChISeq", output.name , "_ALL_COEFFS_IP.txt", sep=""), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  #write.table(final.corr.coeffs[,c("chromosome", "start", "end","control.correction")], paste(output.dir, "/MeDiChISeq", output.name , "_control_correction_IP.bed", sep=""), sep="\t", row.names = FALSE, col.names = F, quote = FALSE)
  
  
  gc()
  return(final.corr.coeffs)
  
}



write.kernel <- function(kernel.obj, file){
  write.table(kernel.obj, file, quote=FALSE)
}



read.kernel <- function(file,skip.header=FALSE){
  
  k <- read.table(file)
  if(skip.header){
    k <- k[3:nrow(k),]
  }
  return(k)
}



fit.peak.profile.seq <- function(file.IP, format="bed", genome="hg19", output.dir=NULL, output.name=NULL, chrom.fit=NULL, limL=0, limU=Inf, reads.elong=150,
                                 quant.cutoff="q1e-7", window=20000, mini.window=2000, wig.res=10, fit.res=50, reads.length=50, 
                                 n.peaks = 50, n.skip = 20, re.fit=100, max.iter=500, selection.method="bic", post.proc.factor=2,
                                 start.pars =c(shape =10, scale = 20), to.be.fit=c("shape", "scale"),
                                 method = "Nelder-Mead", nr.cores=1, remove.clonal.reads=TRUE, clonal.reads.to.keep=3, 
                                 write.pdf=TRUE, save.kernel=TRUE, verbose.console=TRUE, overwrite.wigs=FALSE,  keep.wigs=TRUE,...)
{
  # additionally returns frag.length
  
  start.time <- proc.time()
  
  reads.elong0 <- reads.elong
  step.by <- window - 500
  start.pars <- c(start.pars, h.cutoff = 0)
  
  format <- tolower(format)
  
  if(is.null(output.dir))
    output.dir <- file.path(getwd(), "MeDiChISeq_output")
  
  output.dir <- paste(output.dir, "/", sep="")
  
  if(!file.exists(output.dir))
    dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
  
  output.name <- ifelse(is.null(output.name), "", paste("_", output.name, sep=""))
  
  if(verbose.console){
    sink(file=paste(output.dir,"MeDiChISeq",output.name , "_output_console_fit_peak_profile_seq.txt", sep=""), split=TRUE)
    cat("\n*** fit.peak.profile.seq *** version", VERSION, "*** date", DATE, "***")
    cat("\n\n", date(), "\n\n")
    
    #memo.start <- gc(reset=T)
    #print(memo.start)
    in.args <- c(mget(names(formals()), envir = as.environment(-1)), 
                 sapply(as.list(substitute({
                   ...
                 })[-1]), deparse))
    print(in.args)
    rm(in.args)
  }
  
  chromosomes <- getChromosomes(genome)
  
  if(is.null(chrom.fit)){
    
    tmpChroms <- chromosomes[(names(chromosomes)!="chrX" & names(chromosomes)!="chrY")]
    chrom.fit <- names(tmpChroms[tmpChroms==min(tmpChroms)])
    
  }
  
  output.dir.wigs <- paste(output.dir,"MeDiChISeq", "" , "_WIGS/", sep="" )
  dir.create(output.dir.wigs,showWarnings = FALSE)
  
  if(is.null(reads.elong)){
    
    if(verbose.console)
      cat("\n* Estimation of fragment length...\n\n")
    
    if(remove.clonal.reads){
      #       if(verbose.console)
      #         cat("Removing clonal reads from file... \n")
      file.IP <- remove.clons(file=file.IP, format=format, reads.to.keep=clonal.reads.to.keep, output.dir=output.dir.wigs)
      format="bed"
      remove.clonal.reads <- FALSE
    }
    
    if(verbose.console)
      cat("Converting",file.IP, "into forward and reverse WIG files (reads elongated to",reads.length, "nucleotides)...\n")
    
    wigsOut <- write.wigs.parallel(file=file.IP, output.dir=output.dir.wigs, chromosomes=chrom.fit, sample.type="IP", 
                                   wig.res=wig.res, reads.elong=reads.length, split=TRUE, genome=genome, format=format, zeros=F,
                                   nr.cores=nr.cores, overwrite.wigs=overwrite.wigs, verbose=verbose.console)
    
    
    wig.fw <- load.prepare.wigs.data(path.wigs=wigsOut[1], chr=chrom.fit, limL=limL, limU=limU, verbose=verbose.console, optimal=F)
    if(is.null(wig.fw))
      stop(paste("WIG file for forward strand contains no intensity scores, please use another chromosome for estimating the elongation distance!\n",sep=""))
    
    wig.rev <- load.prepare.wigs.data(path.wigs=wigsOut[2], chr=chrom.fit, limL=limL, limU=limU, verbose=verbose.console, optimal=F)
    if(is.null(wig.rev))
      stop(paste("WIG file for reverse strand contains no intensity scores, please use another chromosome for estimating the elongation distance!\n",sep=""))
    
    if(write.pdf)
      pdf(paste(output.dir,"/MeDiChISeq", output.name, "_details_fitted_kernels_fw.pdf",sep=""))
    
    if(verbose.console)
      cat("\nEstimating fragment length for forward strand from file ", wigsOut[1],"...\n") 
    
    output.fit.peak.profile.fw <- fit.peak.profile.adj(wig.fw, tile.size=wig.res, quant.cutoff=quant.cutoff,fit.res=fit.res, 
                                                       re.fit=re.fit, max.iter=max.iter, plot=write.pdf, 
                                                       start.pars = start.pars, method=method,
                                                       to.be.fit = to.be.fit,window=window, mini.window=mini.window, step.by=step.by,nr.cores=nr.cores, verbose=verbose.console,
                                                       n.peaks = n.peaks, n.skip = n.skip,
                                                       selection.method=selection.method, post.proc.factor=post.proc.factor) 
    
    if(write.pdf)
      dev.off()
    
    if(write.pdf)
      pdf(paste(output.dir,"/MeDiChISeq", output.name, "_details_fitted_kernels_rev.pdf",sep=""))
    
    if(verbose.console)
      cat("\nEstimating fragment length for reverse strand from file ", wigsOut[2],"...\n") 
    
    output.fit.peak.profile.rev <- fit.peak.profile.adj(wig.rev, tile.size=wig.res, quant.cutoff=quant.cutoff,fit.res=fit.res, 
                                                        re.fit=re.fit,max.iter=max.iter, plot=write.pdf, 
                                                        start.pars = start.pars, method=method,
                                                        to.be.fit=to.be.fit,window=window, mini.window=mini.window, step.by=step.by,nr.cores=nr.cores,verbose=verbose.console,
                                                        n.peaks = n.peaks, n.skip = n.skip,
                                                        selection.method=selection.method, post.proc.factor=post.proc.factor) 
    
    if(write.pdf)  
      dev.off()
    
    estimated.fl.fw <- fragment.length(output.fit.peak.profile.fw)
    estimated.fl.rev <- fragment.length(output.fit.peak.profile.rev)
    
    estimated.fl <- (estimated.fl.fw + estimated.fl.rev) - reads.length
    reads.elong <- estimated.fl
    
    if(verbose.console){
      cat("Estimated fragment length forward strand: ", estimated.fl.fw, "\n")
      cat("Estimated fragment length reverse strand: ", estimated.fl.rev, "\n\n")
      cat("Estimated fragment length: ", estimated.fl, "\n\n")
    }
    
  } # end of if(is.null(reads.elong))
  
  
  
  if(remove.clonal.reads){
    #     if(verbose.console)
    #       cat("Removing clonal reads from  file... \n")
    file.IP <- remove.clons(file=file.IP, format=format, reads.to.keep=clonal.reads.to.keep, output.dir=output.dir.wigs) 
    format <- "bed"
  }
  if(verbose.console)
    cat("Converting", file.IP, "into WIG file (reads elongated to",reads.elong, "nucleotides)...\n")
  
  
  elongated.wig <- write.wigs.parallel(file=file.IP, output.dir=output.dir.wigs, chromosomes=chrom.fit, sample.type="IP", 
                                       wig.res=wig.res, reads.elong=reads.elong, split=FALSE, genome=genome, format=format, zeros=F,
                                       nr.cores=nr.cores, overwrite.wigs=overwrite.wigs, verbose=verbose.console)
  
  
  elongated.data <- load.prepare.wigs.data(path.wigs=elongated.wig[1], chr=chrom.fit, limL=limL, limU=limU, verbose=verbose.console, optimal=F)
  if(is.null(elongated.data))
    stop(paste("WIG file contains no intensity scores, please use another chromosome for estimating the elongation distance!\n",sep=""))
  
  
  if(write.pdf)
    pdf(paste(output.dir,"/MeDiChISeq", output.name , "_details_fitted_kernels_final.pdf",sep=""))
  
  if(verbose.console)
    cat("\n\n* Fitting kernel...\n\n")
  
  output.fit.peak.profile.elongated <- fit.peak.profile.adj(data=elongated.data, tile.size=wig.res, quant.cutoff=quant.cutoff,
                                                            fit.res=fit.res, re.fit=re.fit, max.iter=max.iter,plot=write.pdf, 
                                                            start.pars = start.pars, method=method,
                                                            to.be.fit=to.be.fit, window=window, mini.window=mini.window, step.by=step.by, nr.cores=nr.cores,verbose=verbose.console,
                                                            n.peaks = n.peaks, n.skip = n.skip,
                                                            selection.method=selection.method, post.proc.factor=post.proc.factor) 
  
  if(write.pdf)  
    dev.off()
  
  
  estimated.kernel.both <- output.fit.peak.profile.elongated$kernel
  estimated.fl.both <- fragment.length(output.fit.peak.profile.elongated) 
  
  
  if(save.kernel){
    
    if(is.null(reads.elong0)){      
      pdf( paste(output.dir,"/MeDiChISeq", output.name , "_fitted_kernels_rev_fw_final.pdf",sep=""))
      par(mfrow = c(3, 2))
      try(plot.fit.peak.profile.adj(output.fit.peak.profile.rev), silent = T)
      try(plot.fit.peak.profile.adj(output.fit.peak.profile.fw), silent = T)
      try(plot.fit.peak.profile.adj(output.fit.peak.profile.elongated), silent = T)
      dev.off()
    }
    
    pdf( paste(output.dir,"/MeDiChISeq", output.name , "_fitted_kernel_final.pdf",sep=""))
    try(plot.fit.peak.profile.adj(output.fit.peak.profile.elongated), silent = T)
    dev.off()
    
    write.table(reads.elong, paste(output.dir,"/MeDiChISeq", output.name , "_reads_elong.txt",sep=""), row.names = F, col.names = F, quote = F)
    write.table(estimated.kernel.both, paste(output.dir,"/MeDiChISeq", output.name , "_kernel.txt",sep=""), quote=FALSE)
    write.table(estimated.fl.both, paste(output.dir,"/MeDiChISeq", output.name , "_frag_length.txt",sep=""), row.names = F, col.names = F, quote = F)
    
  }
  
  if(verbose.console){
    end.time <- proc.time()
    #memo.end <- gc()
    cat("\n\n* Total time of running fit.peak.profile.seq function:" , time.formatting(end.time["elapsed"]- start.time["elapsed"]) , "!")
    if(save.kernel)
      cat("\n* Output is saved in:", output.dir)
    cat("\n\n", date(), "\n\n")
    sink(file=NULL)
  }
  
  if(!keep.wigs)
    unlink(output.dir.wigs, recursive=T)
  
  gc()
  
  invisible(list(reads.elong=as.numeric(reads.elong), kernel=estimated.kernel.both, frag.length=as.numeric(estimated.fl.both)))  
  
}    



deconv.entire.genome.seq <- function(file.IP, file.Control=NULL, format="bed", genome="hg19", output.dir=NULL, output.name=NULL, chrom.list=NULL, 
                                     limL=0, limU=Inf, potential.peaks=NULL, reads.elong=150, kernel, frag.length, 
                                     quant.cutoff="q1e-5", window=20000, wig.res=10, fit.res=50, 
                                     max.steps=100, selection.method="bic", post.proc.factor=2,                                    
                                     nr.boots=5, local.windows=c(1000, 2000, 5000), Control.corr.param=0.01, 
                                     nr.cores=1, remove.clonal.reads=TRUE, clonal.reads.to.keep=3, 
                                     verbose.console=TRUE, overwrite.wigs=FALSE, keep.wigs=TRUE, ...)
{
  
  start.time <- proc.time()  
  
  # CHECKING THE INPUT ARGUMENTS
  
  if(is.null(output.dir))
    output.dir <- file.path(getwd(), "MeDiChISeq_output")
  
  output.dir <- paste(output.dir, "/", sep="")
  
  if(!file.exists(output.dir))
    dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
  
  output.name <- ifelse(is.null(output.name), "", paste("_", output.name, sep=""))
  
  if(verbose.console){
    sink(file=paste(output.dir,"/MeDiChISeq",output.name , "_output_console_deconv_entire_genome_seq.txt", sep=""), split=TRUE)
    cat("\n*** deconv.entire.genome.seq *** version", VERSION, "*** date", DATE, "***")
    cat("\n\n", date(), "\n\n")  
    #memo.start <- gc(reset=T)
    #print(memo.start)
    in.args <- c(mget(names(formals()), envir = as.environment(-1)), 
                 sapply(as.list(substitute({
                   ...
                 })[-1]), deparse))
    
    print(in.args)
    rm(in.args)
  }
  
  if(!is.data.frame(kernel) && !is.matrix(kernel)){
    kernel <- tryCatch({t <- read.table(kernel)}, 
                       error = function(e) {
                         t <- NULL
                         stop(paste("WRONG kernel! \n Must be either a data frame or a matrix or a path. \n\n",sep=""))
                       })
    
  }
  
  if(!is.numeric(frag.length)){
    frag.length <- tryCatch({t <- as.numeric(read.table(frag.length))}, 
                            error = function(e) {
                              t <- NULL
                              stop(paste("WRONG frag.length! \n Must be either numeric or a path. \n\n",sep=""))
                            })
  }
  
  if(!is.numeric(reads.elong)){
    reads.elong <- tryCatch({t <- as.numeric(read.table(reads.elong))}, 
                            error = function(e) {
                              t <- NULL
                              stop(paste("WRONG reads.elong! \n Must be either numeric or a path. \n\n",sep=""))
                            })
  }
  
  
  if(!is.null(potential.peaks)){
    if(!is.data.frame(potential.peaks)){
      potential.peaks <- tryCatch({t <- read.table(potential.peaks)}, 
                                  error = function(e) {
                                    t <- NULL
                                    stop(paste("WRONG potential.peaks! \n Must be either a data frame or a path. \n\n",sep=""))
                                  })
      
    }
  }
  
  
  all.chromosomes <- names(getChromosomes(genome))  
  
  if(!is.null(potential.peaks)){  
    chrom.list <- as.vector(unique(potential.peaks[,1]))  
  } 
  
  if(is.null(chrom.list)){
    chromosomes <- all.chromosomes
  }else
    chromosomes <- intersect(all.chromosomes, chrom.list)  
  
  if(length(chromosomes)==0){
    cat("Error in chrom.list parameter! Chromosomes must come from: ", all.chromosomes)
    return(NULL)
  }  
  
  # PREPARING IP WIGS
  
  output.dir.wigs.IP <- paste(output.dir,"MeDiChISeq","", "_WIGS/", sep="" )
  dir.create(output.dir.wigs.IP,showWarnings = FALSE)   
  
  if(remove.clonal.reads){
    file.IP <- remove.clons(file=file.IP, reads.to.keep=clonal.reads.to.keep, output.dir=output.dir.wigs.IP, format=format)
    format="bed"
  }
  if(verbose.console)
    cat("Converting",file.IP, "to WIG files (tags elongated to",reads.elong, "nucleotides)...\n")
  
  wigsVecIP <-  write.wigs.parallel(file=file.IP, output.dir=output.dir.wigs.IP, chromosomes=chromosomes, sample.type="IP", 
                                    wig.res=wig.res, reads.elong=reads.elong, split=FALSE, genome=genome, format=format, zeros=F,
                                    nr.cores=nr.cores, overwrite.wigs=overwrite.wigs, verbose=verbose.console)
  
  
  # PREPARING CONTROL WIGS
  
  if(!is.null(file.Control)){
    
    output.dir.wigs.Control <- paste(output.dir,"MeDiChISeq", "", "_WIGS/", sep="" )
    dir.create(output.dir.wigs.Control,showWarnings = FALSE) 
    
    if(remove.clonal.reads)   {
      file.Control <- remove.clons(file=file.Control, format=format, reads.to.keep=clonal.reads.to.keep, output.dir=output.dir.wigs.Control)
      format="bed"
    }
    if(verbose.console)
      cat("Converting",file.Control, "to WIG files (tags elongated to",reads.elong, "nucleotides)...\n")
    
    wigsVecControl <- write.wigs.parallel(file=file.Control, output.dir=output.dir.wigs.Control, chromosomes=chromosomes, sample.type="Control", 
                                          wig.res=wig.res, reads.elong=reads.elong, split=FALSE, genome=genome, format=format, zeros=F,
                                          nr.cores=nr.cores, overwrite.wigs=overwrite.wigs, verbose=verbose.console)
    
  }
  
  
  
  # DECONVOLUTION
  
  
  time.genome <- system.time( results <- lapply(chromosomes, function(chr) {
    
    #chr="chr19"
    wigs.IP <- load.prepare.wigs.data(path.wigs=wigsVecIP[chr], chr, limL=limL, limU=limU, verbose=verbose.console)
    wigs.Control <- NULL
    
    if(!is.null(file.Control))
      wigs.Control <- load.prepare.wigs.data(path.wigs=wigsVecControl[chr], chr, limL=limL, limU=limU, verbose=verbose.console)
    
    output <- parallel.chip.deconv.adj(wigs.IP=wigs.IP, wigs.Control=wigs.Control, chrom.name=chr, kernel=kernel, 
                                       max.steps=max.steps, selection.method=selection.method, post.proc.factor=post.proc.factor,
                                       fit.res=fit.res, window=window, limL=limL,
                                       frag.length=frag.length, potential.peaks=potential.peaks, nr.boots=nr.boots, 
                                       quant.cutoff=quant.cutoff, nr.cores=nr.cores, verbose=verbose.console,
                                       output.dir=output.dir, output.name=output.name, wig.res=wig.res)
    
    gc()
    return(output) 
    
  } ) )
  
  names(results) <- chromosomes
  
  # END OF DECONVOLUTION --> CALCULATING P-VALUES, CONTROL CORRECTION
  
  #   if(verbose.console)
  #     cat("\nTotal time for parallel.chip.deconv.adj() with",nr.cores ,"cores: ", time.genome["elapsed"], "seconds. \n\n")
  
  Merged.coeffs.all <- merging.output(results, output.dir, output.name, nr.boots, if.Control=!is.null(file.Control))
  
  All.coeffs.IP = Merged.coeffs.all$All.coeffs.IP
  All.coeffs.Control = Merged.coeffs.all$All.coeffs.Control
  
  
  if(nr.boots>1){
    All.coeffs.IP <- local.global.p.value(co = Merged.coeffs.all$All.coeffs.IP, all.coeffs.boot = Merged.coeffs.all$All.coeffs.boot.IP, experiment="IP", local.windows=local.windows, 
                                          output.dir=output.dir,output.name=output.name, nr.cores=nr.cores, nr.boots=nr.boots, verbose=verbose.console)
    
    if(!is.null(file.Control)){
      
      All.coeffs.Control <- local.global.p.value(co=Merged.coeffs.all$All.coeffs.Control, all.coeffs.boot=Merged.coeffs.all$All.coeffs.boot.Control, experiment="Control", local.windows=local.windows, 
                                                 output.dir=output.dir,output.name=output.name, nr.cores=nr.cores, nr.boots=nr.boots, verbose=verbose.console)
      
    }
    
    All.coeffs.IP <- Control.corrections(P.values.IP=All.coeffs.IP, P.values.Control=All.coeffs.Control, overlap=Control.corr.param, frag.length=frag.length, output.dir=output.dir, output.name=output.name, nr.cores=nr.cores, verbose=verbose.console)
    
    
  }
  
  
  if(verbose.console){
    end.time <- proc.time()
    cat("\n\n* Total time of deconvolution:" , time.formatting(end.time["elapsed"]- start.time["elapsed"]), "!")
    cat("\n* Output is saved in:", output.dir)
    cat("\n\n", date(), "\n\n")
    sink(file=NULL)
  }
  
  if(!keep.wigs)
    unlink(output.dir.wigs.IP, recursive=T)
  
  
  invisible(list(All.coeffs.IP = All.coeffs.IP))
  
  
}



