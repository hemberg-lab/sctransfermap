#sctransfermap: CoreMethods.R
#
 
#' Description of importdata
#'
#' import_data will import 3 types of single-cell data format to create gene-by-cell matrix
#'
#' @param file file of the single-cell data (ex: scRNA-seq, scATAC)
#' @param opt data format of the single-cell data ("raw": tab-delimited text file, "SingleCellExperiment": SingleCellExperiment rds file, "sparse": MatrixMarket format (written by writeMM)
#' @param datatype "RNAseq" or "atac"
#' @param colname_file colname_file for sparse matrix (usually cell name), when opt is "sparse". Otherwise, it is FALSE
#' @param rowname_file rowname_file for sparse matrix (gene name: official gene symbol or peak location), when opt is "sparse". Otherwise it is FALSE
#'
#' @return a matrix class composed of gene/peak_location(row)-by-cell(column)
#'
#' @examples
#' data <- importdata("filename") (default: opt = "raw", datatype = "RNAseq", colname_file = FALSE, rowname_file = FALSE)
#'
#' @export
importdata <- function(file, opt = "raw", datatype = "RNAseq", colname_file = FALSE, rowname_file = FALSE){
#for SingleCellExperiment -> datatype: RNAseq or atac
	if (opt == "SingleCellExperiment"){
	#	require(scater)
	#	require(SingleCellExperiment)
	#	require(Matrix)
		tmp_data <- readRDS(file)
		e <- as.matrix(counts(tmp_data))
		if (datatype == "RNAseq"){
			colnames(e) <- as.character(colData(tmp_data)$cell_type1)
		}
		else if (datatype == "atac"){
			colnames(e) <- as.character(colData(tmp_data)$cell_label)
		}
		rownames(e) <- as.character(rownames(tmp_data))
		res <- as.matrix(e)
	}
	else if (opt == "raw"){
		data <- read.table(file, sep = "\t", check.name = FALSE)
		res <- as.matrix(data)
	}
	else if (opt == "sparse"){
		tmp_cellnames <- readLines(colname_file)
		cellnames <- tmp_cellnames[2:length(tmp_cellnames)]

		tmp_Rownames <- readlines(rowname_file)
		Rownames <- readLines[2:length(rowname_file)]

		sparse_mat <- readMM(file)
		colnames(sparse_mat) <- cellnames
		rownames(sparse_mat) <- Rownames
		res <- as.matrix(sparse_mat)
	}
	else {
		res <- print("Wrong data type\n")
	}
	
	return(res)
}

#convert_suerat: convert seurat object into matrix
#' Title
#'
#' covert seurat object into matrix format (same result format of importdata function)
#'
#' @param seurtat_object R variable which is a seurat object
#' @param slot_option slot parameter in GetAssayData function in Seurat package (default: "counts") 
#'
#' @return a matrix class composed of gene(row)-by-cell(column)
#'
#' @examples
#'
#' @export
convert_seurat <- function(seurat_object, slot_option = "counts"){
	#require(Seurat)
	res <- GetAssayData(seurat_object, slot = slot_option)
	res<-as.matrix(res)
	return(res)
}

#' Title
#'
#' setup the environment for calling Julia packages and Julia code (CoreFunction.jl) (importing pre-existing julia object is optional)
#'
#' @param julia_obj_index it should be TRUE if importing pre-existing julia object (default FALSE)
#'
#' @return return the julia_object if julia_obj_index is TRUE. Otherwise, there is no output
#'
#' @examples
#'
#' @export
julia_environment <- function(julia_obj_index = FALSE){
	julia_setup()
	cat("\nInstallation of Julia packages if necessary ...\n")
	julia_install_package_if_needed("DelimitedFiles")
	julia_install_package_if_needed("StatsBase")
	julia_install_package_if_needed("Statistics")
#	julia_install_package_if_needed("LinearAlgebra")
	julia_install_package_if_needed("Random")
	julia_install_package_if_needed("Hadamard")
	julia_install_package_if_needed("Embeddings")
#	julia_install_package_if_needed("DataFrames")
#	julia_install_package_if_needed("GLM")
	julia_install_package_if_needed("Combinatorics")
	julia_install_package_if_needed("SparseArrays")
#	julia_install_package_if_needed("DataStructures")
	julia_install_package_if_needed("Compose")
	julia_install_package_if_needed("Gadfly")
#	julia_install_package_if_needed("MultivariateStats")
#	julia_install_package_if_needed("ChemometricsTools")
#	julia_install_package_if_needed("StatsModels")
	julia_install_package_if_needed("JLD")

	cat("\nimporting Julia packages ... \n")
	julia_library("DelimitedFiles")
	julia_library("StatsBase")
	julia_library("Statistics")
#	julia_library("LinearAlgebra")
	julia_library("Random")
	julia_library("Hadamard")
	julia_library("Embeddings")
#	julia_library("DataFrames")
#	julia_library("GLM")
	julia_library("Combinatorics")
	julia_library("SparseArrays")
#	julia_library("DataStructures")
	julia_library("Compose")
	julia_library("Gadfly")
#	julia_library("MultivariateStats")
#	julia_library("ChemometricsTools")
#	julia_library("StatsModels")
	julia_library("JLD")

	cat("\nloading core functions from Julia\n")
	julia_source(system.file("CoreFunction.jl", package = "sctransfermap"))
	if (julia_obj_index) {
		cat("\nImporting Julia object\n")
		julia_file <- system.file("julia_file.jld", package = "sctransfermap")
		julia_obj <- julia_call("load", julia_file, need_return = "R")
		return (julia_obj)
	}	
}

#julia_file.jld
#$side_infonrmation (word2vec)
#$$annotation (information about cell type for side information
#attribute_set ('V' matrix)
#signature_set (for logistic regression modeling)
#meta_model (logistic regression model)



#' Title
#'
#' user-friendly quick prediction for single-cell data
#' side-information is cell ontology
#'
#' @param input_var matrix class composed of gene/peak_location(row)-by-cell(column) (same as result of importdata or convert_seurat)
#' @param name any name that user wants for the result
#' @param opt prediction mode whethere "consensus" or "logistic" (default: "consensus")
#'
#' @return a result of predicted cell types for each query cell cluster 
#'
#' @examples
#'
#' @export
quick_prediction <- function(input_var, name, opt = "consensus"){
#opt: consensus, single (or tissue name), logistic

	cat("\nBefore running this code, user should make matrix object for single-cell data (col: cell, row: gene).\n")

	cat("\nImporting Julia object\n")
	julia_file <- system.file("data/julia_file.jld", package = "sctransfermap")
	julia_obj <- julia_call("load", julia_file, need_return = "R")


	sct_input <- make_input(input_var)
	#julia_obj$side_information <- imported side information (by celltype2vec)
	#julia_obj$annotation <- imported information of celltype in side information
	query_sig <- create_signature(name, "query_tissue", julia_obj$side_information, sct_input[[1]], sct_input[[2]], sct_input[[3]],  julia_obj$annotation, "cellontology", Stopwords = JuliaObject(""))

	#julia_obj$attribute_set <- attribute set ('V' matrix) from each tissue
	if (opt == "consensus"){
		query_res <- predict_cell(julia_obj$attribute_set, query_sig)
	}
	else if (opt == "logistic"){
		query_res <-  logistic_prediction_cell(julia_obj$attribute_set, julia_obj$meta_model, query_sig, tissue_opt = FALSE, celltype_specific = FALSE)
	}	

	query_res <- estimate_celltype(query_res[[1]])
	return (query_res)
}

#' Title
#'
#' quick traing of zero-shot-learing (ZSL)
#'
#' @param input_var matrix class composed of gene/peak_location(row)-by-cell(column) (same as result of importdata or convert_seurat)
#' @param name any name that user wants for the result
#' @param tissue tissue of the sample (please make sure the word of tissue must be same as pre-existing tissue list. If it is new tissue type, it doesn't matter) 
#' Mouse: Aorta, Brain, Diaphragm, Fat, Heart, Kidney, LargeIntestine, LimbMuscle, Liver, Lung, MammaryGland, Marrow, Pancreas, Skin, Spleen, Thymus, Tongue, Trachea
#'
#' @return an attribute for ZSL (trained data) (julia object)
#'
#' @examples
#'
#' @export
quick_training <- function(input_var, name, tissue){
	cat("\nBefore running this code, user should make matrix object for single-cell data (col: cell, row: gene).\n")

	cat("\nImporting Julia object\n")
	julia_file <- system.file("data/julia_file.jld", package = "sctransfermap")
	julia_obj <- julia_call("load", julia_file, need_return = "R")
	sct_input <- make_input(input_var)
	#julia_obj$side_information <- imported side information (by celltype2vec)
	#julia_obj$annotation <- imported information of celltype in side information

	train_att <- create_attribute(name, tissue, julia_obj$side_information, sct_input[[1]], sct_input[[2]], sct_input[[3]], julia_obj$annotation, "cellontology", Stopwords = JuliaObject(""))
	return (train_att)


}

#' Title
#'
#' generating side information of word2vec
#'
#' @param w2v_obj the result of make_w2v
#' @param cell_ontology_obj te result of make_cellontology
#' @param stopword the result of make_stopwords
#'
#' @return side information of cell_ontology in julia object
#'
#' @examples
#'
#' @export
make_w2v_sideinfo <- function(w2v_obj, cell_ontology_obj, stopword, R_output = FALSE){
	cat("\nGenerating side information or word2vec. Use side information object with w2v_index = FALSE during create_signature or create_attribute\n")
	if (R_output) {
		s <- julia_call("make_sideinfo", w2v, cell_ontology_obj, stopword, need_return = "R")
	}
	else {
		s <- julia_call("make_sideinfo", w2v, cell_ontology_obj, stopword, need_return = "Julia")
	}
	return (s)
}


#' Title
#'
#' Embedding Pubmed index file to make word2vec object
#'
#' @param w2v_file Pubmed index file for w2v (ex: http://bio.nlplab.org/ -> Word vectors -> http://evexdb.org/pmresources/vec-space-models/ -> PubMed-w2v.bin)
#'
#' @return an embedded format of word2vec file in julia object 
#'
#' @examples
#'
#' @export
make_w2v <- function(w2v_file){
	cat("\nEmbedding word2vec file from Pubmed Index\n")
	w2v <- julia_call("load_emb", w2v_file, need_return = "Julia")
	return (w2v)
}


#' Title
#'
#' make julia object of stopwords for word2vec
#'
#' @param stopwords stopword file for word2vec (text file: line by line)
#'
#' @return a julia object of stopword file
#'
#' @examples
#'
#' @export
make_stopwords <- function(stopwords){
	stopwords <- JuliaObject(stopwords)
	return (stopwords)
}

#' Title
#'
#' save MatrixMarket format of sparse matrix
#'
#' @param input_var matrix class composed of gene/peak_location(row)-by-cell(column) (same as result of importdata or convert_seurat)
#' @param filename filename to be saved
#'
#' @return no return value. But, it saves sparse matrix in MatrixMarket format with cell annotation file (for column name) and gene name or peak location, etc... (for rowname)
#'
#' @examples
#'
#' @export
save_sparse_data <- function(input_var, filename){
#input_var: output of importdata (=> col: cell, row: peak position or gene)
	cells <- colnames(input_var)
	peakloc <- rownames(input_var)
	mat <- as(input_var, "sparseMatrix")
	
	write.table(cells, paste0(filename, "_cells.csv"), quote = FALSE, row.names = FALSE)
	write.table(peakloc, paste0(filename, "_peaks.csv"), quote = FALSE, row.names = FALSE)
	sparse_mat <- as(counts(mat), "sparseMatrix")
	writeMM(sparse_mat, paste0(filename, "_count.tsv")) 
}


#' Title
#'
#' make input format for sctransfer map from matrix file  
#'
#' @param input_var matrix class composed of gene/peak_location(row)-by-cell(column) (same as result of importdata or convert_seurat)
#' @param min_thr minimum threshold of expressing cell ratio for a given gene/peak_location (default: 0.1)
#' @param max_thr maximum threshold of expressing cell ratio for a given gene/peak_location (default: 0.9)
#' @param pseudo_cnt pseudo count for count matrix for log-transformation (default: 1)
#' @param train TRUE if it is data for training (default: FALSE) 
#' train -> TRUE: adjust min_thr and max_thr. In addition, make the gene/peak_location count for hadamard matrix
#' train -> FALSE: no gene/peak_location filtering
#' @param gene_usage option for gene/peak_location filtering to make hadamard matrix (default: FALSE; train parameter must be TRUE)
#' gene_usage -> FALSE: gene/peak_location filtering by discarding low expressing cell ratio genes/peak_locations
#'
#' @return return a list variable composed of 4 objects (d, genename, cellname, celltype) \cr
#' d: expression or peak matrix (row: gene, col: cell) \cr
#' genename: gene or peak location \cr
#' cellname: cell annotation \cr
#' celltype: celltype == sort(unique(cellname)) \cr
#'
#' @examples 
#' data1 <- make_input(input_var) [query case]
#' data1 <- make_input(input_var, train = TRUE) [training case]
#' usage: d -> data1[[1]], genename -> data1[[2]], cellname -> data1[[3]], celltype -> data1[[4]]
#'
#' @export
make_input <- function(input_var, min_thr = 0.1, max_thr = 0.9, pseudo_cnt = 1.0, train = FALSE, gene_usage = FALSE){
#gene_usage: FALSE: usual hadamard gene filtering, other: 2048, direct gene count
	cat("\nMaking input format for sctransfermap\n")
	d <- input_var[,which(colnames(input_var) != "NA" & colnames(input_var) != "unknown" & colnames(input_var) != "Unknown")]
	cellname <- colnames(input_var)[which(colnames(input_var) != "NA" & colnames(input_var) != "unknown")]

	genesize <- dim(d)[2]
	genes <- rowSums(d!=0) / genesize
#hadamard matrix is only needed in trainset
	if (train) {
		genes_order <- order(genes, decreasing = TRUE)
		d <- d[genes_order,]
		genes <- genes[genes_order]
		genes_cond <- which(genes > min_thr & genes < max_thr)
		d <- d[genes_cond,]
	
		tmp_cnt = dim(d)[1]
		hadamard_number <- 512 * round(floor(tmp_cnt/512))
	
		if (hadamard_number < 512){
			cat("\nNot enough genes are remaind for further analysis. It should be more than at least 512 genes. Change the min_thr or max_thr. Otherwise, it is not suitable for given data\n")
			return(FALSE)
		}
		if (gene_usage){
			d <- d[1:gene_usage,]
		}
		else{
			d <- d[1:hadamard_number,]
		}
	}
	else {
#		genes_cond <- which(genes > min_thr & genes < max_thr)
#		d <- d[genes_cond,]	
	}

	d <- log2(pseudo_cnt + d)
	genename <- rownames(d)
	#cellname <- gsub(" ", "_", cellname)
	celltype <- sort(unique(cellname))

	d<-as.matrix(d)
	rownames(d) <- NULL
	colnames(d) <- NULL
	return (list(d, genename, cellname, celltype))
}

#' Title
#'
#' make cell ontology variable (celltype | description)
#'
#' @param cell_ontology cell_ontology file which is composed of celltype and description of corresponding celltype separated by |
#' @param celltype extract the cell_ontology with a given celltype list or use all the celltypes (default: c())
#' celltype -> c(): use all the cell types in cell_ontology file
#' celltype -> celltype: use only the cell types in cellname variable
#'
#' @return a julia object of cell ontology
#'
#' @examples
#' cell_ontology <- make_cellontology("cell_ontology_file")
#' -> using all the cell types in cell_ontology file
#'
#' d <- make_input(input_var)
#' cell_ontology <- make_cellontology("cell_ontology_file", celltype = d[[4]])
#' @export
make_cellontology <- function(cell_ontology, celltype = c()){
#default usage: cell_ontology = cell_ontology
	cat("\nImporting cellontology file for side information\n")
	cell_ontology <- read.table(cell_ontology, sep = "|", quote = "")
	#cell_names <- as.character(cell_ontology$V1)
	#cell_names <- gsub(" ", "_", cell_names)
	#cell_ontology$V1 <- cell_names


	cell_ontology <- as.matrix(cell_ontology)
	if (length(celltype) > 0){
		head = TRUE
		for (i in c(1:dim(cell_ontology)[1])){
			cell <- cell_ontology[i,1]
			if (is.element(cell, celltype)){
				if (head){
					res <- cell_ontology[i,]
					head = FALSE
					next
				}
				res <- rbind(res, cell_ontology[i,])
				
			}
      		}
	}
	else {
		res <- cell_ontology
	}	
	res <- data.frame(res)
	cell_names <- as.character(res$V1)

	julia_order <- julia_call("sort", JuliaObject(cell_names), need_return = "R")
	new_index <- r2julia_sort(cell_names, julia_order)

	res <- res[new_index,]
	res <- as.matrix(res)
	rownames(res) <- NULL		
	res <- JuliaObject(res)
	return (res)
}

#' Title
#'
#' Generation of side-information from expression or peak matrix
#'
#' @param Side_mat matrix class composed of gene/peak_location(row)-by-cell(column) (same as result of importdata or convert_seurat)
#' @param cellnames cell name of Side_mat (column name)
#' @param Sdim the size of side-information for each cell type (default: 200)
#'
#' @return a side-information which could be directly used in create_siganture or create_attribute with w2v_index = FALSE (julia_object)
#'
#' @examples
#' d <- make_input(input_var)
#' s <- generate_sideinformation(d[[1]], d[[3]])
#' @export
generate_sideinformation <- function(Side_mat, Cellnames, Sdim = 200){
	side_mat <- JuliaObject(Side_mat)
#Side_mat: ex: output of make_input(atac-seq)[[1]]
	cellname <- JuliaObject(Cellname) 
#output of make_input(atac-seq)[[3]]
	sdim <- julia_call("Int64", Sdim, need_return = "Julia")
	res <- julia_call("generate_sideinformation", side_mat, cellname, sdim, need_return = "Julia")
	return (res)
}
#result should be directly used with w2v_index = FALSE


#' Title
#'
#' make attribute object for ZSL (trained data)
#'
#' @param Name user-defined name for the output
#' @param Tissue tissue of the sample (please make sure the word of tissue must be same as pre-existing tissue list. If it is new tissue type, it doesn't matter)
#' Mouse: Aorta, Brain, Diaphragm, Fat, Heart, Kidney, LargeIntestine, LimbMuscle, Liver, Lung, MammaryGland, Marrow, Pancreas, Skin, Spleen, Thymus, Tongue, Trachea
#' @param w2v word2vector object
#' output of make_w2v
#' pre-existing side-information object (julia object) (ex: julia_obj$side_information) 
#' output of generate_sideinformation
#' @param Expr expression/peak matrix (first index from the result of make_input: d <- make_input(input_var), d[[1]])
#' @param Genename gene/peak_location of Expr (second index from the result of make_input: d<- make_input(input_var), d[[2]]) 
#' @param Cellname cellname of Expr (third index from the result of make_input: d<- make_input(input_var), d[[3]])
#' @param cell_ontology julia object of cell ontology (the result of make_cellontology)
#' @param Sdatatype type of side-information (cell_ontology, ATAC-seq, or user-defined term)
#' @param w2v_index TRUE if user needs to make side-information from cell ontology. FALSE if julia object for side-information exists (default: FALSE)
#' @param Xdatatype type of Expr (default: scRNA-seq; or scATAC-seq)
#' @param Ngene a number of vector size after JL transformation (default: 64)
#' @param Stopwords julia object of stopwords (the result of make_stopwords) (default: pre-existing variable: stopwords)
#' @param Fdrthr fdr threshold for calculating similarity score during training (default: 0.1)
#' @param Gamma paramter for ZSL (default: 1)
#' @param Lambda parameter for ZSL (default: 1)
#' @param Xpeaks TRUE if data type of x (or Expr) is ATAC-seq (default: FALSE)
#' @param Sep delimiter for peak_location (default: "_")
#' ex: chr1_3084739_3085712 (make sure delimiter for peak_location of data must be "_") 
#' @param Mincellcount minimum cell count for a given cell type for training (default: 10)
#' @param Mincellfrac minimum fraction of cell type for training (default: 0.05)
#' @param Simthres initial similarity threshold (default: -Inf)
#' @param R_output TRUE ifthe  result of the function should be R object. FALSE if the result of the function should be julia object (default: FALSE)
#'
#' @return an attribute for ZSL (trained data) (julia object)
#'
#' @examples
#'
#' @export
create_attribute <- function(Name, Tissue, w2v, Expr, Genename, Cellname, cell_ontology, Sdatatype, w2v_index = FALSE, Xdatatype = "scRNA-seq", Ngene = 64, Stopwords = stopwords, Fdrthr = 0.1, Gamma = 1.0, Lambda = 1.0, Xpeaks = FALSE, Sep = "_", Mincellcount = 10, Mincellfrac = 0.05, Simthres = -Inf, R_output = FALSE){
#Sdatatype: cell ontology filenmae or ATAC
	name <- JuliaObject(Name)
	tissue <- JuliaObject(Tissue)
	#w2v
	expr <- JuliaObject(Expr)
	genename <- JuliaObject(Genename)
	cellname <- JuliaObject(Cellname)
	#celltype <- JuliaObject(Celltype)
	#cell_ontology
	sdatatype <- JuliaObject(Sdatatype)
	w2v_index <- JuliaObject(w2v_index)
	xdatatype <- JuliaObject(Xdatatype)
	ngene <- julia_call("Int64", Ngene, need_return = "Julia")
	#stopwords
	fdrthr <- julia_call("Float64", Fdrthr, need_return = "Julia")
	gamma <- julia_call("Float64", Gamma, need_return = "Julia")
	lambda <- julia_call("Float64", Lambda, need_return = "Julia")
	xpeaks <- JuliaObject(Xpeaks)
	sep <- JuliaObject(Sep)
	simthres <- JuliaObject(Simthres)
	mincellcount <- julia_call("Int64", Mincellcount, need_return = "Julia")
	mincellfrac <- julia_call("Float64", Mincellfrac, need_return = "Julia")

	cat("\nMaking model 'V' matrix by training data\n")
	if (R_output){
		Index_res <- julia_call("create_attribute", name, tissue, w2v, w2v_index, expr, genename, cellname, cell_ontology, sdatatype, xdatatype, ngene, Stopwords, fdrthr, gamma, lambda, xpeaks, sep, mincellfrac, mincellcount, simthres,  need_return = "R") 
	}	
	else {
		Index_res <- julia_call("create_attribute", name, tissue, w2v, w2v_index, expr, genename, cellname, cell_ontology, sdatatype, xdatatype, ngene, Stopwords, fdrthr, gamma, lambda, xpeaks, sep, mincellfrac, mincellcount, simthres, need_return = "Julia")
	}
	return (Index_res)
}


#' Title
#'
#' make signature object for ZSL (input for prediction)
#'
#' @param Name user-defined name for the output
#' @param Tissue tissue of the sample (please make sure the word of tissue must be same as pre-existing tissue list. If it is new tissue type, it doesn't matter)
#' Mouse: Aorta, Brain, Diaphragm, Fat, Heart, Kidney, LargeIntestine, LimbMuscle, Liver, Lung, MammaryGland, Marrow, Pancreas, Skin, Spleen, Thymus, Tongue, Trachea
#' @param w2v word2vector object
#' output of make_w2v
#' pre-existing side-information object (julia object) (ex: julia_obj$side_information)
#' output of generate_sideinformation
#' @param Expr expression/peak matrix (first index from the result of make_input: d <- make_input(input_var), d[[1]])
#' @param Genename gene/peak_location of Expr (second index from the result of make_input: d<- make_input(input_var), d[[2]])
#' @param Cellname cellname of Expr (third index from the result of make_input: d<- make_input(input_var), d[[3]])
#' @param cell_ontology julia object of cell ontology (the result of make_cellontology)
#' @param Sdatatype type of side-information (cell_ontology, ATAC-seq, or user-defined term)
#' @param w2v_index TRUE if user needs to make side-information from cell ontology. FALSE if julia object for side-information exists (default: FALSE)
#' @param Stopwords julia object of stopwords (the result of make_stopwords) (default: pre-existing variable: stopwords)
#' @param Xdatatype type of Expr (default: scRNA-seq; or scATAC-seq)
#' @param R_output TRUE ifthe  result of the function should be R object. FALSE if the result of the function should be julia object (default: FALSE)
#' @param celltype_cnt_thr if the number for cell type exceeds, warning message will be sent and return FALSE. To ignore it, change the parameter into FALSE
#'
#' @return an attribute for ZSL (trained data) (julia object)
#'
#' @examples
#'
#' @export
create_signature <- function(Name, Tissue, w2v, Expr, Genename, Cellname, cell_ontology,Sdatatype, w2v_index = FALSE, Stopwords = stopwords, Xdatatype = "genes",  R_output = FALSE, celltype_cnt_thr = TRUE){
	celltype_cnt = length(unique(Cellname))
	if (celltype_cnt_thr) {
		if (celltype_cnt > 1000){
			cat("\n Too many celltypes to predict (more than 1000). It will take lots of time to predict. If you still want to run, make parameter 'celltype_cnt_thr = FALSE'\n")
			return (FALSE)
		}
	}
	name <- JuliaObject(Name)
	tissue <- JuliaObject(Tissue)
	#w2v
	expr <- JuliaObject(Expr)
	genename <- JuliaObject(Genename)
	cellname <- JuliaObject(Cellname)
	#celltype <- JuliaObject(Celltype)
	#cell_ontology
	sdatatype <- JuliaObject(Sdatatype)
	#stopwords
	xdatatype <- JuliaObject(Xdatatype)
	w2v_index <- JuliaObject(w2v_index)

	cat("\nMaking signature format of query (or test set) for prediction\n")
	if (R_output){
		Map_res <- julia_call("create_signature", name, tissue, w2v, w2v_index, expr, genename, cellname, cell_ontology, sdatatype, Stopwords, xdatatype, need_return = "R")
	}
	else {
		Map_res <- julia_call("create_signature", name, tissue, w2v, w2v_index, expr, genename, cellname, cell_ontology, sdatatype, Stopwords, xdatatype, need_return = "Julia")
	}
	return (Map_res)
}


#' Title
#'
#' adding attribute object into attribute set
#'
#' @param attribute julia object of attribute (the result of create_attribute or quick_training)
#' @param iscell pre-existing attribute set (julia object)
#'
#' @return merged attribute set 
#'
#' @examples
#'
#' @export
add_attribute_set <- function(attribute, iscell){
	iscell <- julia_call("IndexList", attribute, iscell, need_return = "Julia")
	return (iscell)

}	


#' Title
#'
#' making attribute set
#'
#' @param attribute julia object of attribute (the result of create_attribute or quick_training)
#'
#' @return initiallized attribute set
#'
#' @examples
#'
#' @export
make_attribute_set <- function(attribute){
	iscell <- julia_call("IndexList", attribute, need_return = "Julia")
	return (iscell)
}

#' Title
#'
#' adding signature object into signature set
#'
#' @param signature julia object of signature (the result of create_signature)
#' @param mscell per-existing signature set
#'
#' @return merged signature set
#'
#' @examples
#'
#' @export
add_signature_set <- function(signature, mscell){
	mscell <- julia_call("MappingList", signature, mscell, need_return = "Julia")
	return (mscell)
}

#' Title
#'
#' making signature set 
#'
#' @param signature julia object of signature (the result of create_signature)
#'
#' @return initiallized signature set
#'
#' @examples
#'
#' @export
make_signature_set <- function(signature){
	mscell <- julia_call("MappingList", signature, need_return = "Julia")
	return (mscell)
}

#' Title
#'
#' make a confusion matrix of query data by ZSL with consensus manner
#'
#' @param iscell attribute set
#' @param query signature object of query (the result of create_signature)
#' @param simthr default simthreshold for prediction (default: Inf -> considers only simthr in attribute set)
#' @param gn_overlap_thr minimum gene overlap ratio of attribute in attribute set between query data (default: 0.5)
#' @param R_output TRUE ifthe  result of the function should be R object. FALSE if the result of the function should be julia object (default: TRUE)
#'
#' @return a confusion matrix of prediction result
#' row: query cell type or cluster name
#' column: cell types in side-information (from query)
#' first column: unassign count
#'
#' @examples
#'
#' @export
predict_cell <- function(iscell, query, simthr = Inf, gn_overlap_thr = 0.5, R_output = TRUE){
	simthr <- JuliaObject(simthr)
	gn_overlap_thr <- JuliaObject(gn_overlap_thr)
	cat("\nPredicting celltype by model from ZSL in consensus manner with various tissue\n")
	if (R_output){
		confusion <- julia_call("transfer_map_consensus", iscell, query, simthr, gn_overlap_thr, need_return = "R")
		

		y <- table(field(query, "celltypes"))
		tmp <- JuliaObject(field(query, "celltypes"))
		tmp <- julia_call("unique", tmp, need_return = "Julia")
		julia_cell_list <- julia_call("sort", tmp, need_return = "R")
		y_cell_list <- rownames(y)
		new_index <- r2julia_sort(y_cell_list, julia_cell_list)
		y <- y[new_index]
		cellnames <- y_cell_list[new_index]

		confusion <- decorate_confusion(confusion, y, cellnames, field(query, "snames"))
	}
	else {
		confusion <- julia_call("transfer_map_consensus", iscell, query, simthr, gn_overlap_thr, need_return = "Julia")
	}
	return (confusion)
}

#' Title
#'
#' make a index for a given list order 
#'
#' @param r_list list that want to change
#' @param julia_list target list
#' @param Index TRUE if the result should be an index for r_list. FALSE if the result should be the changed list or r_list
#'
#' @return index or changed list of r_list
#'
#' @examples
#'
#' @export
r2julia_sort <- function(r_list, julia_list, Index = TRUE){
#sort r list into the order of other list (ex: julia sorted list)
#Index True: return only new index for r list, Index False: return reordered r list (== julia list)
	new_index <- c()
	for (i in julia_list){
		new_index <- c(new_index, match(i, r_list))
	}

	if (Index){
		return(new_index)
	}
	else {
		return(r_list[new_index])
	}
}




#' Title
#'
#' add unassign count from confusion matrix
#'
#' @param confusion confusion matrix from "transfer_map_consensus" function in Julia
#' @param y table of count for each cell type (query)
#' y <- table(field(query_signature, "celltypes"))
#' @param cellnames cell type of query data
#' @param sname cell type of side-information
#'
#' @return a confusion matrix that has unassign count
#'
#' @examples
#'
#' @export
decorate_confusion <- function(confusion, y, cellnames, sname){
	confusion <- as.matrix(confusion)
	assign_sum <- apply(confusion, 1, sum)
	unassign <- y - assign_sum
	unassign <- as.matrix(unassign)
	unassign <- data.frame(unassign)
	rownames(unassign ) <- cellnames
	colnames(unassign) <- "unassign"
	colnames(confusion) <- sname
	confusion <- cbind(unassign, confusion)
	return (confusion)
}

#' Title
#'
#' convert count confusion matrix into ratio of each query cell type, respectively
#'
#' @param confusion confusion matrix from predict_cell
#'
#' @return ratio of each query cell type 
#'
#' @examples
#'
#' @export
confusion_ratio <- function(confusion) {
	total <- apply(confusion, 1, sum)
	ratio_confusion <- confusion / total
	
	return (ratio_confusion)
}

#' Title
#'
#' validation of known data
#'
#' @param con_ratio confusion ratio matrix from confusion_ratio 
#'
#' @return True postivie rate and unassign ratio for each cell type, respectively 
#'
#' @examples
#'
#' @export
validation <- function(conf_ratio){
	tp = c()
	unassign = conf_ratio[,1]
	for (i in c(1:dim(conf_ratio)[1])){
		tp_index <- which(colnames(conf_ratio)==rownames(conf_ratio)[i])
		tp <- c(tp, conf_ratio[i,tp_index])

	}
	tp <- as.matrix(tp)
	unassign <- as.matrix(unassign)
	res<-cbind(tp, unassign)
	colnames(res) <- c("TP", "Unassign")
	rownames(res) <- rownames(conf_ratio)
	return(res)

}
#during creating logistic regression model, it should be stacked with Indexdata that uses all the celltypes without mincell thresholding. 
#' Title
#'
#' creating logistic regression model 
#'
#' @param iscell attribute set
#' @param mscell signature set
#' @param gn_overlap_thr minimum gene overlap ratio of attribute in attribute set between query data (default: 0.5)
#' @param ncore number of cores for multi-threading (default = 1)
#'
#' @return return list composed of (model_list, sname, zsl_tissue, meta_column_name)
#' model_list: logistic regression model (composed of predicting cell and corresponding coefficients)
#' sname: cell types that are used in each coefficient (T cell, B cell, fibroblast ....)
#' zsl_tissue: tissue types that are used in logistic regression (Lung, Heart ...)
#' meta_column_name: pseudo name of variables in logistic regression (Heart_cell1, Heart_cell2...)
#'
#' @examples
#'
#' @export
create_meta_index_av <- function(iscell, mscell, gn_overlap_thr = 0.5, ncore = 1){
	cat("\nCreating logistic regression model")
	cat("\nGenerating ZSL results for stacking from each tissue\n")

	registerDoParallel(ncore)
	start_time <- Sys.time()

	gn_overlap_thr <- JuliaObject(gn_overlap_thr)

	zsl_tissue <- c()
	sname <- field(mscell[1], "snames")
#	model_list <- list()
	first = TRUE
	for (q_index in 1:length(mscell)){
		query = mscell[q_index]
		comment <- paste("\nTraining set: ", as.character(q_index), "/", as.character(length(mscell)), sep = "")
		cat(comment)
		meta_info <- julia_call("create_meta_index", iscell, query, gn_overlap_thr, need_return = "R")
		meta_table_tmp <- data.frame(meta_info[[1]])
		if (first){
		#make ff (formula), tissue & tissue size, meta_table_column	
		#	first = FALSE (after meta_table)
			head = TRUE
			meta_column_name <- c()
			for (tissue in meta_info[[2]]){
				if (tissue %in% zsl_tissue){
					next
				}
				zsl_tissue <- c(zsl_tissue, tissue)
				#zsl_tissue: tissue list used in ZSL prediction (=tissue of attribute sets)
				for (cell_index in 1:length(sname)){
					cellname <- paste("cell", as.character(cell_index), sep = "")
					tissue_cell <- paste(tissue, cellname, sep = "_")
					meta_column_name <- c(meta_column_name, tissue_cell)
					if (head){
						head <- FALSE
						ff <- paste("F~", tissue_cell, sep = "")
						next
					}
					ff <- paste(ff, tissue_cell, sep = "+")
				}
			}
			train_tissue <- as.character(meta_info[[2]])
			train_tissue_unique <- unique(train_tissue)
		#unique doesn't change the order	
			train_tissue_overlap_index <- list()
			for (w in 1:length(train_tissue_unique)){
				temp_tissue <- train_tissue_unique[w]
				train_tissue_overlap_index[[w]] <- which(train_tissue == temp_tissue) 
			}	

			train_tissue_size <- length(train_tissue)	

		}
		block_size <- dim(meta_table_tmp)[2] / train_tissue_size
		for (w in 1:length(train_tissue_overlap_index)){
			temp_index <- train_tissue_overlap_index[[w]]
			temp_first = TRUE
			
			for (i in temp_index){
				pre_tmp_table <- meta_table_tmp[,(1+(i-1)*block_size):(i*block_size)]
				pre_tmp_table <- t(pre_tmp_table)
				if (temp_first){
					tmp_table <- pre_tmp_table
					temp_first = FALSE
					next
				}		
				tmp_table <- tmp_table + pre_tmp_table		

			}
			tmp_table <- tmp_table / length(temp_index)
			if (i == 1){
				meta_table <- tmp_table
				next
			}
			meta_table <- cbind(meta_table, tmp_table)
			
		}

		meta_table <- data.frame(meta_table)
		cells <- field(query, "celltypes")
		if (first){
			first = FALSE
			total_meta_table <- meta_table
			total_cells <- cells
			next
		}
	
		total_meta_table <- rbind(total_meta_table, meta_table)
		total_cells <- c(total_cells, cells)

	}
	total_meta_table <- data.frame(total_meta_table)
	colnames(total_meta_table) <- meta_column_name
	rownames(total_meta_table) <- c(1:dim(total_meta_table)[1])

#just for reducing long rownames
	total_celltype <- sort(unique(total_cells))
	cat("\nGenerating logistic regression model ...\n\n")
#	for (j in 1:length(total_celltype)){

	model_list <- foreach (j = 1:length(total_celltype)) %dopar% {
		comment <- paste("\nLogistic regression modeling: ", as.character(j), "/", as.character(length(total_celltype)), "\n", sep = "")
		cat(comment)
		cell <- total_celltype[j]
		tmp <- total_cells
		tmp[which(total_cells == cell)] <- 1
		tmp[which(total_cells != cell)] <- 0

		total_meta_table <- as.matrix(total_meta_table)
		tmp <- as.numeric(tmp)
		fit <- fastglm(total_meta_table, tmp, family = binomial, method = 3)
#		tmp <- as.factor(tmp)
#		ff <- as.formula(ff)	
#		F <- tmp
#		fit <- glm(ff, data = total_meta_table, family = binomial)
	}	
	
	names(model_list) <- total_celltype
	end_time <- Sys.time()
	proc_time <- end_time - start_time
	print(proc_time)	
	registerDoParallel(1)
	return (list(model_list, sname, zsl_tissue, meta_column_name))

}	



#tissue_opt: show tissue origin
#celltype_specific: during matching the tissue, whether uses all the celltype or only the matched celltype
#' Title
#'
#' predict cell type with logistic regression model
#'
#' @param iscell attribute set
#' @param meta_model_list result from create_meta_index
#' @param query signature object of query
#' @param gn_overlap_thr minimum gene overlap ratio of attribute in attribute set between query data (default: 0.5)
#' @param logistic_thr threshold for each cell whether predicted score from logistic regression model is confident (default: 0.5)
#' @param tissue_opt TRUE if user wants to find tissue-origin (default: FALSE)
#' @param celltype_specific TRUE if user wants to calculate the weight fo each tissue with only matched cell type. Otherwise, weight for each tissue will be calculated with whole celltypes (default: FALSE; tissue_opt should be TRUE to use this parameter)
#'
#' @return return 2 confusion matrix 
#' res[[1]]: a confusion matrix of prediction result
#' row: query cell type or cluster name
#' column: cell types in side-information (from query)
#' first column: unassign count
#' res[[2]]: a confusion matrix of prediction result in each cell
#' row: query cell
#'
#' @examples
#'
#' @export
logistic_prediction_cell <- function(iscell, meta_model_list, query, gn_overlap_thr = 0.5, logistic_thr = 0.5, tissue_opt = FALSE, celltype_specific = FALSE){
	cat("\nPredicting celltype by logistic regression model\n")
	meta_cellname <- as.character(meta_model_list[[2]])
	gn_overlap_thr <- JuliaObject(gn_overlap_thr)
	predict_info <- julia_call("create_meta_index", iscell, query, gn_overlap_thr, need_return = "R")
	
	query_model_cellname <- as.character(predict_info[[3]])
	query_model_tissue <- as.character(predict_info[[2]])
	query_model_tissue_uniq <- unique(query_model_tissue)

	model_list <- meta_model_list[[1]]
	
	fullmeta_cellname <- meta_model_list[[4]]
#	tmp_glm_tissue_cell <- meta_model_list[[4]]
	glm_tissue <- meta_model_list[[3]]

#make real tissue + cell name (tissue cellname for glm is tissue_cell1, tissue_cell2 ...)
	
	real_tissue_cellname <- c()
	for (tmp_tissue in glm_tissue){
		for(tmp_cell in meta_cellname){
			tmp_tissue_cell <- paste(tmp_tissue, tmp_cell, sep = "_")
			real_tissue_cellname <- c(real_tissue_cellname, tmp_tissue_cell)
		}
	}
	##make query_table with (x*v*s score) and column: new cellname (full_meta_cellname) / row: cells
	query_table_tmp <- data.frame(predict_info[[1]])

	tissue_overlap_index <- list()	
	for (w in 1:length(query_model_tissue_uniq)){
		temp_tissue <- query_model_tissue_uniq[w]
		tissue_overlap_index[[w]] <- which(query_model_tissue == temp_tissue)
	}


	new_cellnames <- c()
	query_overlap_tissue <- c()
	for (i in query_model_tissue){
		if (i %in% query_overlap_tissue){
			next
		}
		query_overlap_tissue <- c(query_overlap_tissue, i)
		for (j in query_model_cellname){

			new_cellname <- paste(i, j, sep = "_")
			new_cellnames <- c(new_cellnames, new_cellname)
		}
	}
	block_size <- dim(query_table_tmp)[2] / length(query_model_tissue)
	for (w in 1:length(tissue_overlap_index)){
		temp_index <- tissue_overlap_index[[w]]
		temp_first = TRUE

		for (i in temp_index){
			pre_tmp_table <- query_table_tmp[,(1+(i-1)*block_size):(i*block_size)]
			pre_tmp_table <- t(pre_tmp_table)
			if (temp_first){
				tmp_table <- pre_tmp_table
				temp_first = FALSE
				next
			}
		}
			tmp_table <- tmp_table / length(temp_index)
		 if (i == 1){
			 query_table <- tmp_table
			 next
		}
		query_table <- cbind(query_table, tmp_table)
	}
#if query stypenames (celltype) exceed celltype in logistic regression model, it automatically discards missing exceeded celltype
#if query has missing celltype, zero matrix will be given for logistic regression model
	new_cell_index <- r2julia_sort(new_cellnames, real_tissue_cellname)
	query_final_cellname <- new_cellnames[new_cell_index]
	
	query_table <- query_table[,new_cell_index]
	query_table[is.na(query_table)] <- 0
	colnames(query_table) <- meta_model_list[[4]]
	query_table <- data.frame(query_table)
	
	##


	##pre-generating empty all_confusion (all cells) and confusion (each celltype)
	cells <- field(query, "celltypes")
	celltype <- sort(unique(cells))
	tmp_confusion <- matrix(nrow = length(cells), ncol = length(names(model_list)))
		

	if (tissue_opt){
		model_celltype <- real_tissue_cellname
		all_confusion <- matrix(0, nrow = length(cells), ncol = length(model_celltype))

	}
	else{
		model_celltype <- names(model_list)
		all_confusion <- matrix(0, nrow = length(cells), ncol = length(model_celltype))
	}

	confusion <- matrix(0,nrow = length(celltype), ncol = length(model_celltype))
	confusion <- data.frame(confusion)
	all_confusion <- data.frame(all_confusion)
	##

	##prediction by logistic regression model
	for (j in 1:length(model_list)){
		comment <- paste("\nPrediction by logistic regression modeling: ", as.character(j), "/", as.character(length(model_list)), "celltype\n", sep = "")
		cat(comment)
		predict_all <- predict(model_list[[j]], as.matrix(query_table), type = "response")
#		predict_all <- predict(model_list[[j]], query_table, type = "response")
		tmp_confusion[,j] <- predict_all
	}	
	##

	##searching best matched celltype (max score)
	max_predict_cell <- apply(tmp_confusion, 1, which.max)
	max_values <- apply(tmp_confusion, 1, max)
	#thresholding
	thr_index <- max_values >= logistic_thr
	##

	##build confusion matrix
	if (tissue_opt){
		cat("\nTissue option will take long time\n")
		coef_model <- list()
		
		for (j in 1:length(model_list)){
			tmp_coef <- coef(model_list[[j]])
			tmp_coef[is.na(tmp_coef)] <- 0
			tmp_coef <- tmp_coef[2:length(tmp_coef)]
			coef_model[[j]] <- tmp_coef
		}
		if (celltype_specific){
			
		}
		else {
			block_size <- length(coef_model) / length(glm_tissue)
		}
	}
	cat("\nMaking confusion matrix...\n")
	for (i in 1:length(celltype)) {
		cell = celltype[i]
		assign_index <- which(cells == cell & thr_index)

		if (length(assign_index) == 0){
			next
		}

		tmp_index <- max_predict_cell[which(cells == cell & thr_index)]

		if (tissue_opt){
			for (j in 1:length(assign_index)){
				tmp_assign_index <- assign_index[j]
				matched_coef <- coef_model[[tmp_index[j]]]
				score_coef <- matched_coef * query_table[j,]
				tmp_tissue_conf <- matrix(0, nrow = 1, ncol = length(glm_tissue))
				tmp_tissue_conf <- data.frame(tmp_tissue_conf)
				predicted_celltype <- as.character(names(model_list)[tmp_index[j]])

				if (celltype_specific){
					celltype_specific_score_coef <- score_coef[which(meta_cellname == predicted_celltype)]
					tmp_tissue_conf[1,] <- celltype_specific_score_coef
					

				}
				else{
					for (tissue_index in 1:length(glm_tissue)){
						tmp_tissue_conf[1,tissue_index] <- sum(score_coef[(1+(tissue_index-1)*block_size):(tissue_index*block_size)])
					}
				}
				max_tissue <- apply(tmp_tissue_conf, 1, which.max)
				matched_tissue <- glm_tissue[max_tissue]
				matched_tissue_cell <- paste(matched_tissue, cells[tmp_assign_index], sep = "_") 	
				new_index <- which(model_celltype == matched_tissue_cell)


				all_confusion[tmp_assign_index, new_index] <- all_confusion[tmp_assign_index, new_index] + 1
				confusion[i, new_index] <- confusion[i, new_index] + 1
			}
		}
		else {
			for (j in tmp_index){
				confusion[i,j] <- confusion[i,j] + 1
			}
			for (q in 1:length(assign_index)){
				j = tmp_index[q]
				j1 = assign_index[q]
				all_confusion[j1,j] <- all_confusion[j1,j] + 1
			}
			
		}
	}
	##
	##make row, col names (same with julia sort)
	y <- table(cells)
	
	julia_y_list <- julia_call("sort", JuliaObject(celltype), need_return = "R")
	y_cell_list <- rownames(y)
	new_index <- r2julia_sort(y_cell_list, julia_y_list)
	y <- y[new_index]

	new_index_row <- r2julia_sort(celltype, julia_y_list)
	
	
	julia_model_names <- julia_call("sort", JuliaObject(model_celltype), need_return = "R")
	new_index_col <- r2julia_sort(model_celltype, julia_model_names)
	reordered_model_celltype <- model_celltype[new_index_col]
	
	confusion <- confusion[new_index_row, new_index_col]

	if (tissue_opt){	
		query_tissue <- field(query, "tissue")
		new_celltype <- c()
		for (i in celltype){
			new_tmp <- paste(query_tissue, i, sep = "_")
			new_celltype <- c(new_celltype, new_tmp)
		}
	}
	else{
		new_celltype<- celltype
	}

	confusion <- decorate_confusion(confusion, y, new_celltype, reordered_model_celltype)
	all_confusion <- all_confusion[,new_index_col]
	colnames(all_confusion) <- reordered_model_celltype
	return (list(confusion, all_confusion))
	##
}



#'Title
#'
#' prediction of cell type by a given confusion matrix from "predict_cell". Assume binomial distribution: n: total cell count for each query celltype / p: 1/celltypes from cell ontology / X (for z score): cell count for each celltype from cell ontology. Only consider max cell count (draw -> mixed)
#'
#' @param conf confusion matrix from "predict_cell"
#'
#' @return suggests matched cell type for corresponding query cell type with q-value, respectively  
#'
#' @examples
#'
#' @export
estimate_celltype <- function(conf){

	tmp1 <- max.col(conf, ties.method = "first")
	tmp2 <- max.col(conf, ties.method = "last")
	available <- tmp1==tmp2
	filtered_conf <- conf[available,]
	mixed_conf <- conf[!available,]
	
	nonzero_set <- rowSums(filtered_conf != 0)
	p <- 1/nonzero_set
	n_set <- apply(filtered_conf, 1, sum)
	#nonzero for q value correction
	max_values <- apply(filtered_conf, 1, max)
	max_index <- max.col(filtered_conf)
	predicted_celltype <- colnames(conf)[max_index]

	mean_set <- n_set * p
	sigma_set <- sqrt(n_set * p * (1-p))

	z_set <- (max_values - mean_set) / sigma_set
	pvalues <- pnorm(-abs(z_set))
	
	qvalues <- pvalues * nonzero_set
	tmp_mat <- matrix(nrow = dim(conf)[1], ncol = 2)
	rownames(tmp_mat) <- rownames(conf)
	colnames(tmp_mat) <- c("celltype", "q-value")
	tmp_mat[available, 1] <- colnames(conf)[max_index]
	tmp_mat[available, 2] <- qvalues
	mixed_annot <- rep("mixed", dim(mixed_conf)[1])
	na_annot <- rep("NA", dim(mixed_conf)[1])

	tmp_mat[!available, 1] <- mixed_annot
	tmp_mat[!available, 2] <- na_annot	
	res <- cbind(tmp_mat, conf)
	return(res)


}	

#saving julia object
#julia_call("save, "data/julia_storage.jld", "variable name", variable name, ...)
#loading julia object
#julia_data <- julia_call("load", "data/julia_storage.jld", need_return = "R")
#usage: julia_data$ ...


