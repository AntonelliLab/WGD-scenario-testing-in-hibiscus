# Libs ----
cat('Loading libraries and functions ....', '\n')
extract_slot <- aicw_calc <- MLEGeneCount2 <- NULL
suppressMessages(suppressWarnings(library(WGDgc)))
source(file.path('tools', 'mldgenecount2.R'))
source(file.path('tools', 'other.R'))

# Data ----
cat('Importing data ....', '\n')
# Gene counts
gene_counts <- read.csv(file.path('data', 'gene_counts.csv'),
                        stringsAsFactors = FALSE)
# Tree strings
trstrs <- readLines(file.path('data', 'scenario_treestrings.txt'))
trstrs <- strsplit(x = trstrs, split = ' = ')
trstrs_scenarios <- vapply(X = trstrs, FUN = '[[', i = 1,
                           FUN.VALUE = character(1))
trstrs <- vapply(X = trstrs, FUN = '[[', i = 2, FUN.VALUE = character(1))
names(trstrs) <- trstrs_scenarios

# Testing ----
cat('Running scenario tests ....', '\n')
model_res <- vector(mode = 'list', length(trstrs))
names(model_res) <- names(trstrs)
for (nm in names(trstrs)) {
  cat('.... ', nm, '\n')
  tree <- read.simmap(text = trstrs[[nm]])
  tmp <- MLEGeneCount2(tr = tree, geneCountData = gene_counts, dirac = 1,
                       conditioning = "twoOrMore")
  model_res[[nm]] <- tmp
}

# Output ----
cat('Writing out results ....', '\n')
saveRDS(object = model_res, file = file.path('results',
                                             'raw_model_output.RData'))
# model_res <- readRDS(file = file.path('results', 'raw_model_output.RData'))
printable_res <- data.frame(scenario = names(model_res), loglikelihood =
                              extract_slot(model_res = model_res,
                                           slt_nm = 'loglikelihood'),
                            AIC = extract_slot(model_res = model_res,
                                               slt_nm = 'aic'))
rownames(printable_res) <- NULL
# add weighted AIC
printable_res <- aicw_calc(printable_res)
write.csv(x = printable_res, row.names = FALSE, quote = FALSE,
          file = file.path('results', 'model_results.csv'))
cat('Done!\n')
