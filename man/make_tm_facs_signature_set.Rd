% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_code.R
\name{make_tm_facs_signature_set}
\alias{make_tm_facs_signature_set}
\title{make_tm_facs_signature_set}
\usage{
make_tm_facs_signature_set(tm_facs_input, s, annot, data_name, sdatatype)
}
\arguments{
\item{tm_facs_input}{result of make_input_tm_facs}

\item{s}{side-infomration (julia object)}

\item{annot}{cell_ontology (julia object)}

\item{data_name}{name for the signature}

\item{sdatatype}{type of side-information (cell_ontology, ATAC-seq, or user-defined term)}
}
\value{
signature set
}
\description{
creating signature set from result of make_input_tm_facs
}
\examples{
tm_facs_input <- make_input_tm_facs(tm_facs)
w2v <- make_w2v("Pubmed.index")
cell_ontology <- make_cellontology("cell_ontogloy_file")
stopword <- make_stopwords(stopwords)
side_info <- make_w2v_sideinfo(w2v, cell_ontology, stopword)
sig_set <- make_tm_facs_signature_set(tm_facs_input, side_info, cell_ontology, "user", "cell_ontology")

}
