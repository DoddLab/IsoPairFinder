% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunctionOthers.R
\name{readMSP}
\alias{readMSP}
\title{readMSP}
\usage{
readMSP(
  file,
  mode = c("all", "standard"),
  source = c("MetAnalyzer", "MSDIAL", "Other")
)
}
\arguments{
\item{file}{the file name}

\item{mode}{standard: extract name, mz and RT from the file, which fit for MSP data exported from QI software; all: extract all information in the file}
}
\value{
A list object. info: MS1 information, including ms1, rt etc.; spec: the MS/MS spectrum
}
\description{
read MSP spectra files
}
\examples{
test <- readMSP(file = 'F:/01 MetIMMS/00 data processing/190515 external validation msms data extraction/zhumetlib_validation_pos_20v_190520.msp', mode = 'all')
}
\author{
Zhiwei Zhou
\email{zhouzw@stanford.edu}
}
