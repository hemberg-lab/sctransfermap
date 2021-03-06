% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_code.R
\name{tissue_merge}
\alias{tissue_merge}
\title{tissue_merge}
\usage{
tissue_merge(tissue)
}
\arguments{
\item{tissue}{input tissue name divided by "_"
\itemize{
\item tissue_name: Heart_1, Heart_2, Lung_1, Lung_2 ...
\item BrainNonMyeloid: Brain
\item BrainMyeloid: Brain
}}
}
\value{
collapsed tissue name
}
\description{
merge same tissue type from tissue name (internal function for where tissue name should be merged)
}
\examples{
new_tissue <- tissue_merge("Heart_1")

}
