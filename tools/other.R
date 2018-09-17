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

extract_slot <- function(model_res, slt_nm) {
  vapply(X = model_res, FUN = function(x) {
    x[[slt_nm]]}, FUN.VALUE = numeric(1))
}