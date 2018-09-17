aicw_calc <- MLEGeneCount2 <- NULL

# Libs
suppressMessages(suppressWarnings(library(WGDgc)))
source('mldgenecount2.R')

# Functions
extract_slot <- function(model_res, slt_nm) {
  vapply(X = model_res, FUN = function(x) {
    x[[slt_nm]]}, FUN.VALUE = numeric(1))
}

# Gene counts
gene_counts <- read.csv('gene_counts.csv', stringsAsFactors = FALSE)

# Tree strings
trstrs <- readLines('scenario_treestrings.txt')
trstrs <- strsplit(x = trstrs, split = ' = ')
trstrs_scenarios <- vapply(X = trstrs, FUN = '[[', i = 1, FUN.VALUE = character(1))
trstrs <- vapply(X = trstrs, FUN = '[[', i = 2, FUN.VALUE = character(1))
names(trstrs) <- trstrs_scenarios

# Testing
model_res <- vector(mode = 'list', length(trstrs))
names(model_res) <- names(trstrs)
for (nm in names(trstrs)) {
  cat('.... ', nm, '\n')
  tree <- read.simmap(text = trstrs[[nm]])
  tmp <- MLEGeneCount2(tr = tree, geneCountData = gene_counts, dirac = 1,
                       conditioning = "twoOrMore")
  model_res[[nm]] <- tmp
}
saveRDS(object = model_res, file = 'raw_model_output.RData')

# Printable results
printable_res <- data.frame(scenario = names(model_res), loglikelihood =
                              extract_slot(model_res = model_res,
                                           slt_nm = 'loglikelihood'),
                            AIC = extract_slot(model_res = model_res,
                                               slt_nm = 'aic'))
rownames(printable_res) <- NULL
# add weighted AIC
printable_res <- aicw_calc(printable_res)
write.csv(x = printable_res, file = 'model_results_no_bonus.csv',
          row.names = FALSE, quote = FALSE)
