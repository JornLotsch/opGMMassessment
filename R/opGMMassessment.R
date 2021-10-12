# Analysis of a Gaussian mixture structure in the data Statistical
# justification with likelihood ratio tests, AIC or BIC, or GAP criterion
#' @importFrom ClusterR GMM
#' @importFrom mclust densityMclust
#' @importFrom mixtools normalmixEM
#' @importFrom parallel detectCores makeCluster clusterExport stopCluster parLapply
#' @importFrom NbClust NbClust
#' @importFrom AdaptGauss InformationCriteria4GMM LikelihoodRatio4Mixtures
#' @importFrom methods hasArg
#' @importFrom utils head tail sessionInfo
#' @importFrom stats ks.test rnorm sd
#' @importFrom cluster pam maxSE
#' @importFrom multimode modetest
#' @importFrom DistributionOptimization DistributionOptimization
#' @importFrom mixAK NMixMCMC
opGMMassessment <- function(Data, FitAlg = "MCMC", Criterion = "LR", MaxModes = 8,
  MaxCores = 2048, PlotIt = FALSE, KS = FALSE, Seed) {

  DIM <- function(...) {
    args <- list(...)
    unlist(lapply(args, function(x) {
      if (is.null(dim(x)))
        return(1)
      dim(x)[2]
    }))
  }

  if (!hasArg("Data"))
    stop("GMMassessment: No data.")
  if (length(Data) < 2)
    stop("GMMassessment: Too few data.")
  if (DIM(Data) != 1 | is.numeric(Data) == FALSE)
    stop("GMMassessment: Data must be a one-dimensional numerical vector.")
  list.of.FitAlgs <- c("ClusterRGMM", "densityMclust", "DO", "MCMC", "normalmixEM")
  if (!FitAlg %in% list.of.FitAlgs)
    stop("GMMassessment: Fit algorithm not implemented. Use ClusterRGMM, DO, densityMclust or normalmixEM")
  list.of.Criteria <- c("AIC", "BIC", "FM", "GAP", "LR", "NbClust", "SI")
  if (!Criterion %in% list.of.Criteria)
    stop("GMMassessment: Criterion not implemented. Use AIC, BIC, FM, GAP, LR, NbClust, or SI.")

  requireNamespace("parallel")
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    num_workers <- 2L
  } else {
    num_workers <- parallel::detectCores()
  }
  nProc <- min(num_workers - 1, MaxCores, MaxModes)

  # mclapply hack for Windows
  mclapply.hack <- function(...) {
    size.of.list <- length(list(...)[[1]])
    cl <- makeCluster(min(size.of.list, detectCores()))
    loaded.package.names <- c(sessionInfo()$basePkgs, names(sessionInfo()$otherPkgs))
    tryCatch({
      this.env <- environment()
      while (identical(this.env, globalenv()) == FALSE) {
        clusterExport(cl, ls(all.names = TRUE, envir = this.env), envir = this.env)
        this.env <- parent.env(environment())
      }
      clusterExport(cl, ls(all.names = TRUE, envir = globalenv()), envir = globalenv())
      parLapply(cl, 1:length(cl), function(xx) {
        lapply(loaded.package.names, function(yy) {
          require(yy, character.only = TRUE)
        })
      })
      return(parLapply(cl, ...))
    }, finally = {
      stopCluster(cl)
    })
  }

  mclapply <- switch(Sys.info()[["sysname"]], Windows = {
    mclapply.hack
  }, Linux = {
    mclapply
  }, Darwin = {
    mclapply
  })

  DataOrigLength <- length(as.vector(Data))
  DataOrignoNA <- which(!is.na(Data) & !is.infinite(Data))
  Data <- as.vector(Data[DataOrignoNA])
  GMMdata <- Data
  list.of.Modes <- 1:MaxModes
  BestGMM <- 1
  Means <- as.vector(mean(GMMdata, na.rm = TRUE))
  SDs <- as.vector(sd(GMMdata, na.rm = TRUE))
  Weights <- 1
  Mixtures <- cbind(Means, SDs, Weights)

  if (!missing(Seed)) {
    ActualSeed <- Seed
  } else {
    ActualSeed <- tail(get(".Random.seed", envir = globalenv()), 1)
  }

  is.integer0 <- function(x) {
    is.integer(x) && length(x) == 0L
  }

  idBestGMM_AICBIC <- function(GMMdata, GMMfit, Criterion, ActualSeed) {
    AICBIC <- lapply(1:length(GMMfit), function(x) {
      AICBICi <- AdaptGauss::InformationCriteria4GMM(Data = GMMdata, Means = lapply(GMMfit,
        "[[", 2)[[x]][, 1], SDs = lapply(GMMfit, "[[", 2)[[x]][, 2], Weights = lapply(GMMfit,
        "[[", 2)[[x]][, 3])
      return(AICBICi)
    })
    AICBIC <- unlist(lapply(1:length(GMMfit), function(x) {
      switch(Criterion, AIC = {
        AICBIC <- AICBIC[[x]]$AIC
      }, BIC = {
        AICBIC <- AICBIC[[x]]$BIC
      })
      return(AICBIC)
    }))
    BestGMM <- which.min(AICBIC)
    return(BestGMM)
  }

  idBestGMM_FM <- function(GMMdata, MaxModes, nProc) {
    list.of.Modes <- 1:(MaxModes - 1)
    if (.Platform$OS.type != "windows" & MaxCores > 1) {
      ExcessMassi <- mclapply(list.of.Modes, function(x) {
        pExM <- multimode::modetest(Data, method = "FM", mod0 = x, B = 60)$p.value
        return(pExM)
      }, mc.cores = nProc)
    } else {
      ExcessMassi <- lapply(list.of.Modes, function(x) {
        pExM <- multimode::modetest(Data, method = "FM", mod0 = x, B = 60)$p.value
        return(pExM)
      })
    }
    BestGMM <- 1
    for (i in 2:length(ExcessMassi)) {
      if (ExcessMassi[i] < 0.05)
        firstBestGMM <- i else break
    }
    return(BestGMM)
  }

  idBestGMM_GAP <- function(GMMdata, MaxModes, ActualSeed, nProc) {
    pam1 <- function(x, k) {
      list(cluster = cluster::pam(x, 3, metric = "euclidean", stand = FALSE,
        cluster.only = TRUE))
    }
    gsPam1 <- clusGapP(x = cbind(GMMdata, GMMdata), FUNcluster = pam1, K.max = MaxModes,
      B = 60, spaceH0 = "original", MaxCores = nProc)
    BestGMM <- with(gsPam1, maxSE(Tab[, "gap"], Tab[, "SE.sim"], method = "globalSEmax"))
    return(BestGMM)
  }

  idBestGMM_LR <- function(GMMdata, GMMfit, ActualSeed) {
    LRi <- c(1, unlist(lapply(2:MaxModes, function(x) {
      AdaptGauss::LikelihoodRatio4Mixtures(Data = GMMdata, NullMixture = lapply(GMMfit,
        "[[", 2)[[x - 1]], OneMixture = lapply(GMMfit, "[[", 2)[[x]], PlotIt = FALSE)$Pvalue
    })))
    firstBestGMM <- 1
    for (i in 2:length(LRi)) {
      if (LRi[i] < 0.05)
        firstBestGMM <- i else break
    }
    BestGMM <- firstBestGMM
    return(BestGMM)
  }

  idBestGMM_NbClust <- function(GMMdata, MaxModes) {
    NBres <- suppressWarnings(NbClust::NbClust(GMMdata, method = "kmeans", distance = "euclidean",
      max.nc = MaxModes))
    BestGMM <- length(unique(NBres$Best.partition))
    return(BestGMM)
  }

  idBestGMM_SI <- function(GMMdata, MaxModes, nProc) {
    list.of.Modes <- 1:(MaxModes - 1)
    if (.Platform$OS.type != "windows" & MaxCores > 1) {
      ExcessMassi <- mclapply(list.of.Modes, function(x) {
        pExM <- multimode::modetest(Data, method = "SI", mod0 = x, B = 60)$p.value
        return(pExM)
      }, mc.cores = nProc)
    } else {
      ExcessMassi <- lapply(list.of.Modes, function(x) {
        pExM <- multimode::modetest(Data, method = "SI", mod0 = x, B = 60)$p.value
        return(pExM)
      })
    }
    BestGMM <- 1
    for (i in 2:length(ExcessMassi)) {
      if (ExcessMassi[i] < 0.05)
        firstBestGMM <- i else break
    }
    return(BestGMM)
  }

  # Start of GMM fit Do the GMM fit
  nProc <- min(num_workers - 1, MaxModes, MaxCores)

  switch(FitAlg, ClusterRGMM = {
    GMMfit <- lapply(list.of.Modes, function(x, Mixtures = Mixtures) {
      set.seed(ActualSeed)
      GMMfit_Mode <- try(ClusterR::GMM(data = data.frame(GMMdata), gaussian_comps = list.of.Modes[x],
        dist_mode = "eucl_dist"), TRUE)
      if (class(GMMfit_Mode) != "try-error") {
        Mixtures <- cbind(GMMfit_Mode$centroids, sqrt(GMMfit_Mode$covariance_matrices),
          GMMfit_Mode$weights)
      } else Mixtures <- Mixtures
      return(list(GMMfit_Mode, Mixtures))
    })
  }, densityMclust = {
    GMMfit <- mclapply(list.of.Modes, function(x, Mixtures = Mixtures) {
      set.seed(ActualSeed)
      GMMfit_Mode <- try(mclust::densityMclust(data = GMMdata, G = x), TRUE)
      if (class(GMMfit_Mode) != "try-error") {
        res = GMMfit_Mode$parameters
        Mixtures <- cbind(unname(res$mean), sqrt(unname(res$variance$sigmasq)),
          unname(res$pro))
      } else Mixtures <- Mixtures
      return(list(GMMfit_Mode, Mixtures))
    }, mc.cores = nProc)
  }, DO = {
    GMMfit <- mclapply(list.of.Modes, function(x, Mixtures = Mixtures) {
      set.seed(ActualSeed)
      GMMfit_Mode <- try(DistributionOptimization::DistributionOptimization(Data = GMMdata,
        Modes = list.of.Modes[x], Monitor = 0, CrossoverRate = 0.9, ErrorMethod = "chisquare",
        Seed = ActualSeed), TRUE)
      if (class(GMMfit_Mode) != "try-error") {
        Mixtures <- cbind(GMMfit_Mode$Means, GMMfit_Mode$SDs, GMMfit_Mode$Weights)
      } else Mixtures <- Mixtures
      return(list(GMMfit_Mode, Mixtures))
    }, mc.cores = nProc)
  }, MCMC = {
    GMMfit <- mclapply(list.of.Modes, function(x, Mixtures = Mixtures) {
      set.seed(ActualSeed)
      Prior <- list(priorK = "fixed", Kmax = x)
      nMCMC <- c(burn = 5000, keep = 10000, thin = 5, info = 1000)
      GMMfit_Mode <- try(mixAK::NMixMCMC(y0 = GMMdata, prior = Prior, nMCMC = nMCMC),
        TRUE)
      if (class(GMMfit_Mode) != "try-error") {
        MeansMCMC <- GMMfit_Mode[[1]]$poster.mean.mu * GMMfit_Mode[[1]]$scale$scale +
          GMMfit_Mode[[1]]$scale$shift
        SDsMCMC <- sqrt(GMMfit_Mode[[1]]$scale$scale^2 * as.numeric(GMMfit_Mode[[1]]$poster.mean.Sigma))
        WeightsMCMC <- GMMfit_Mode[[1]]$poster.mean.w[seq(0, (GMMfit_Mode[[1]]$nx_w -
          1) * GMMfit_Mode[[1]]$prior$Kmax, by = GMMfit_Mode[[1]]$prior$Kmax) +
          c(1:x)]
        Mixtures <- cbind(MeansMCMC, SDsMCMC, WeightsMCMC)
      } else Mixtures <- Mixtures
      return(list(GMMfit_Mode, Mixtures))
    }, mc.cores = nProc)
  }, normalmixEM = {
    GMMfit <- mclapply(list.of.Modes, function(x, Mixtures = Mixtures) {
      set.seed(ActualSeed)
      GMMfit_Mode <- try(mixtools::normalmixEM(GMMdata, mu = kmeans(GMMdata,
        x)$centers, ECM = TRUE, maxrestarts = 1e+05), TRUE)
      if (class(GMMfit_Mode) != "try-error") {
        Mixtures <- cbind(GMMfit_Mode$mu, GMMfit_Mode$sigma, GMMfit_Mode$lambda)
      } else Mixtures <- Mixtures
      return(list(GMMfit_Mode, Mixtures))
    }, mc.cores = nProc)
  })

  # Identify best fit based on selected criterion
  switch(Criterion, AIC = {
    BestGMM <- idBestGMM_AICBIC(GMMdata = GMMdata, GMMfit = GMMfit, Criterion = Criterion,
      ActualSeed = ActualSeed)
  }, BIC = {
    BestGMM <- idBestGMM_AICBIC(GMMdata = GMMdata, GMMfit = GMMfit, Criterion = Criterion,
      ActualSeed = ActualSeed)
  }, FM = {
    BestGMM <- idBestGMM_FM(GMMdata = GMMdata, MaxModes = MaxModes, nProc = nProc)
  }, GAP = {
    BestGMM <- idBestGMM_GAP(GMMdata = GMMdata, MaxModes = MaxModes, nProc = nProc,
      ActualSeed = ActualSeed)
  }, LR = {
    BestGMM <- idBestGMM_LR(GMMdata = GMMdata, GMMfit = GMMfit, ActualSeed = ActualSeed)
  }, NbClust = {
    BestGMM <- idBestGMM_NbClust(GMMdata = GMMdata, MaxModes = MaxModes)
  }, SI = {
    BestGMM <- idBestGMM_SI(GMMdata = GMMdata, MaxModes = MaxModes, nProc = nProc)
  }, BestGMM <- BestGMM)

  # Extract GMM parameters
  switch(FitAlg, ClusterRGMM = {
    Means <- as.vector(GMMfit[[BestGMM]][[1]]$centroids)
    SDs <- sqrt(as.vector(GMMfit[[BestGMM]][[1]]$covariance_matrices))
    Weights <- as.vector(GMMfit[[BestGMM]][[1]]$weights)
  }, densityMclust = {
    Means <- as.vector(GMMfit[[BestGMM]][[2]][, 1])
    SDs <- as.vector(GMMfit[[BestGMM]][[2]][, 2])
    Weights <- as.vector(GMMfit[[BestGMM]][[2]][, 3])
  }, DO = {
    Means <- as.vector(GMMfit[[BestGMM]][[1]]$Means)
    SDs <- as.vector(GMMfit[[BestGMM]][[1]]$SDs)
    Weights <- as.vector(GMMfit[[BestGMM]][[1]]$Weights)
  }, MCMC = {
    Means <- as.vector(GMMfit[[BestGMM]][[2]][, 1])
    SDs <- as.vector(GMMfit[[BestGMM]][[2]][, 2])
    Weights <- as.vector(GMMfit[[BestGMM]][[2]][, 3])
  }, normalmixEM = {
    Means <- as.vector(GMMfit[[BestGMM]][[2]][, 1])
    SDs <- as.vector(GMMfit[[BestGMM]][[2]][, 2])
    Weights <- as.vector(GMMfit[[BestGMM]][[2]][, 3])
  })

  # Calculate Bayes boundaries
  Boundaries <- c()
  ClassesB <- rep(1, length(GMMdata))
  if (BestGMM > 1) {
    Boundaries <- AdaptGauss::BayesDecisionBoundaries(Means = Means, SDs = SDs,
      Weights = Weights)
    if (is.integer0(Boundaries) == FALSE)
      Boundaries <- Boundaries[Boundaries >= min(Means) & Boundaries <= max(Means)]
    if (length(Boundaries) > 0)
      ClassesB <- cutGMM(x = GMMdata, breaks = Boundaries)
  }

  Classes <- rep(NA, length(DataOrigLength))
  Classes[DataOrignoNA] <- ClassesB

  # Do Kolmogorov-Smirnov test
  if (KS == TRUE) {
    set.seed(ActualSeed)
    Pred <- CreateGMM(Means = Means, SDs = SDs, Weights = Weights, n = 1000)$Data
    KStest <- suppressWarnings(ks.test(x = GMMdata, y = Pred))
  } else KStest <- NA

  # Prepare plot
  p1 <- GMMplotGG(Data = GMMdata, Means = Means, SDs = SDs, Weights = Weights,
    Hist = TRUE)
  if (PlotIt == TRUE)
    print(p1)

  # Return results
  return(list(Cls = Classes, Means = Means, SDs = SDs, Weights = Weights, Boundaries = Boundaries,
    Plot = p1, KS = KStest))
}
