MLEGeneCount2 <- function(tr, geneCountData, mMax = NULL, geomMean = NULL,
                          dirac = NULL, useRootStateMLE = FALSE,
                          conditioning = c("oneOrMore", "twoOrMore",
                                           "oneInBothClades", "none"),
                          equalBDrates = FALSE, fixedRetentionRates = FALSE, 
                          startingBDrates = c(0.01, 0.02), startingQ = NULL) {
  if (isTRUE(all.equal(geneCountData, round(geneCountData)))) 
    geneCountData <- round(geneCountData)
  else stop("data are not of integer type")
  if (dim(geneCountData)[2] != length(tipLabels(tr))) 
    stop("the number of species in the tree does not match the number of species in the data file")
  if (sum(names(geneCountData) %in% tipLabels(tr)) != length(tipLabels(tr))) 
    stop("the column names do not match the species names")
  if (is.null(mMax)) 
    mMax = max(apply(geneCountData, 1, sum))
  if ((mMax != floor(mMax)) || mMax <= 0) 
    stop("mMax must be a positive integer. Was ", mMax)
  if (mMax < max(apply(geneCountData, 1, sum))) 
    warning("the value of the parameter mMax should be larger to avoid approximations.", 
            immediate. = TRUE)
  if (mMax < max(geneCountData)) 
    stop(paste("the value of mMax should be no smaller than the largest gene count:", 
               max(geneCountData)))
  isNotNullGeomMean = !is.null(geomMean)
  isNotNullDirac = !is.null(dirac)
  if ((isNotNullGeomMean + isNotNullDirac + useRootStateMLE) != 
      1) 
    stop("Use exactly one of these three arguments: geomMean or dirac or useRootStateMLE")
  if (isNotNullDirac) {
    isDiracInteger = (dirac == floor(dirac))
    if (!isDiracInteger) {
      dirac <- floor(dirac)
      warning("the dirac value has to be an integer. It has been rounded down", 
              immediate. = TRUE)
    }
    if (dirac < 1) 
      stop("the dirac value needs to be positive.")
  }
  if (isNotNullGeomMean) {
    if (geomMean < 1) 
      stop("the mean of the prior geometric distribution has to be greater or equal to 1")
  }
  if (useRootStateMLE) {
    if (conditioning == "oneInBothClades" | conditioning == 
        "twoOrMore") 
      stop(paste0("The 'rootStateMLE' is not implemented with the conditioning '", 
                  conditioning, "'"))
  }
  input = processInput(tr, equalBDrates, fixedRetentionRates, 
                       startingBDrates, startingQ)
  nWGD = sum(input$phyloMat$type == "WGD")
  nWGT = sum(input$phyloMat$type == "WGT")
  nWG = nWGD + nWGT
  if ((nWG == 0) & (!fixedRetentionRates)) {
    fixedRetentionRates = TRUE
    warning("The tree has no WGD or WGT events, setting fixedRetentionRates to TRUE", 
            immediate. = TRUE)
  }
  para = input$para
  lower = input$lower
  upper = input$upper
  if (!is.null(geomMean)) {
    geomProb = 1/geomMean
  }
  else {
    geomProb = NULL
  }
  conditioning = match.arg(conditioning)
  if (conditioning == "twoOrMore") {
    geneNumberInEachFamily <- rowSums(geneCountData, dims = 1)
    familiesNotAppropriate = which(geneNumberInEachFamily == 
                                     1)
    if (length(familiesNotAppropriate) > 0) {
      print("The following families are responsible for the error:\n  ")
      print(familiesNotAppropriate)
      stop("conditioning: we cannot use the type twoOrMore\n since at least one family has only 1 gene copy")
    }
  }
  else if (conditioning == "oneInBothClades") {
    WGDgc:::.checkClades(input$phyloMat, geneCountData, input$nLeaf)
  }
  result <- optim(para, getLikGeneCount, input = input,
                  geneCountData = geneCountData, 
                  mMax = mMax, geomProb = geomProb, dirac = dirac,
                  useRootStateMLE = useRootStateMLE,
                  conditioning = conditioning,
                  fixedRetentionRates = fixedRetentionRates, 
                  equalBDrates = equalBDrates, method = "L-BFGS-B",
                  lower = lower, upper = upper)
  para = result$par
  npara <- length(para)
  loglik = -result$value
  if (equalBDrates) {
    lambda = exp(para[1])
    lenlammu = 1
    mu = "same as birth rate"
  }
  else {
    lambda = exp(para[1])
    mu = exp(para[2])
    lenlammu = 2
  }
  if (!fixedRetentionRates) {
    input$wgdTab$retain2 = para[(lenlammu + 1):length(para)]
    input$wgdTab$retain1 = 1 - input$wgdTab$retain2
    if (nWGT) {
      iT = input$wgdTab$type == "WGT"
      input$wgdTab$retain3[iT] = input$wgdTab$retain2[iT]^2
      input$wgdTab$retain2[iT] = 2 * input$wgdTab$retain2[iT] * 
        input$wgdTab$retain1[iT]
      input$wgdTab$retain1[iT] = input$wgdTab$retain1[iT]^2
    }
  }
  aic_val <- -2 * loglik + 2*npara
  return(list(birthrate = lambda, deathrate = mu, loglikelihood = loglik, 
              WGDtable = input$wgdTab, phyloMat = input$phyloMat,
              call = match.call(), convergence = result$convergence,
              mMax = mMax, aic = aic_val))
}

aicw_calc <- function(printable_res) {
  # https://stat.ethz.ch/pipermail/r-help/2008-December/183087.html
  # printable_res <- printable_res[order(printable_res$AIC), ]
  # for (i in seq_len(nrow(printable_res))) {
  #   printable_res$diff[i] <- printable_res$AIC[1] - printable_res$AIC[i]
  # }
  # printable_res$wi <- exp(1)^(-0.5 * printable_res$diff)
  # printable_res$aic.weights <- printable_res$wi / sum(printable_res$wi)
  # printable_res
  cbind(printable_res, qpcR::akaike.weights(printable_res$AIC))
}
