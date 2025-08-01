% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunctionPeakGrouping.R
\name{convertSemiTargetedResult2LongTable}
\alias{convertSemiTargetedResult2LongTable}
\title{convertSemiTargetedResult2LongTable}
\usage{
convertSemiTargetedResult2LongTable(result_semi_table)
}
\description{
convert wide semi targeted clustering result to long table
}
\examples{
result_semi_table <- readr::read_csv('I:/00_projects/03_MetDNA2/00_data/20200805_targeted_annotation_evaluation_200STD/01_targeted_annotation/200STD_neg_reverse/annot_credential/semi_targeted_annotation_result.csv')
test <- convertSemiTargetedResult2LongTable(result_semi_table)
}
\author{
Zhiwei Zhou
}
