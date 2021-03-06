% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoreMethods.R
\name{decorate_confusion}
\alias{decorate_confusion}
\title{decorate_confusion}
\usage{
decorate_confusion(confusion, y, cellnames, sname)
}
\arguments{
\item{confusion}{confusion matrix from "transfer_map_consensus" function in Julia}

\item{y}{table of count for each cell type (query)\cr
- y <- table(field(query_signature, "celltypes"))}

\item{cellnames}{cell type of query data}

\item{sname}{cell type of side-information}
}
\value{
a confusion matrix that has unassign count
}
\description{
add unassign count from confusion matrix (internal function of predict_cell)
}
\examples{
NA

}
