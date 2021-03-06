% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoreMethods.R
\name{make_w2v}
\alias{make_w2v}
\title{make_w2v}
\usage{
make_w2v(w2v_file)
}
\arguments{
\item{w2v_file}{Pubmed index file for w2v\cr
- ex: http://bio.nlplab.org/ -> Word vectors -> http://evexdb.org/pmresources/vec-space-models/ -> PubMed-w2v.bin}
}
\value{
an embedded format of word2vec file in julia object
}
\description{
Embedding Pubmed index file to make word2vec object
}
\examples{
w2v <- make_w2v("PubMed-w2v.bin")

}
