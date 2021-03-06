% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findCorrelation.R
\name{findCorrelation}
\alias{findCorrelation}
\title{Function to perfrom correlation analysis in TADs}
\usage{
findCorrelation(object, method = "pearson", adj.pval = FALSE,
  plot.proportions = FALSE)
}
\arguments{
\item{object}{InTADSig object with signals and genes combined in TADS}

\item{method}{Correlation method: "pearson" (default), "kendall", "spearman"}

\item{adj.pval}{Perform p-value adjsutment and include q-values in result}

\item{plot.proportions}{Plot proportions of signals and genes in correlation}
}
\value{
A table with correlation values for signal-gene pairs including
correlation p-value, euclidian distance and rank.
}
\description{
This function combines genes and signals in inside of TADs
}
\examples{
## perform analysis on test data
inTadSig <- newSigInTAD(enhSel, enhSelGR, rpkmCountsSel, txsSel)
inTadSig <- filterGeneExpr(inTadSig, geneType = "protein_coding")
inTadSig <- combineInTAD(inTadSig, tadGR)
corData <- findCorrelation(inTadSig, method="pearson")

}
