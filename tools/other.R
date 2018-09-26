akaike.weights <- function(x) {
  # From qpcR
  x <- x[!is.na(x)]
  delta.aic <- x - min(x, na.rm = TRUE)
  rel.LL <- exp(-0.5 * delta.aic)
  sum.LL <- sum(rel.LL, na.rm = TRUE)
  weights.aic <- rel.LL/sum.LL
  return(list(deltaAIC = delta.aic, rel.LL = rel.LL, weights = weights.aic))
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
  cbind(printable_res, akaike.weights(printable_res$AIC))
}

extract_slot <- function(model_res, slt_nm) {
  vapply(X = model_res, FUN = function(x) {
    x[[slt_nm]]}, FUN.VALUE = numeric(1))
}
