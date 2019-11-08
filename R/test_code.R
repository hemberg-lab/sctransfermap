#sctransfermap: test_code.R
#
#start_time <- Sys.time()
#function
#end_time <- Sys.time()
#proc_time <- end_time - start_time
#proc_time <- as.numeric(proc_time)

#object.size(variable_name)

#' import_tm_facs
#'
#' "importdata" for TM (or multiple "raw" matrix files)
#'
#' @param address scdata address list (format: "raw")\cr
#' - filename: tissue_1_cnt, tissue_2_cnt ...
#'
#' @return list of matrix (names: each tissue or filename)
#'
#' @examples
#' tm_facs <- import_tm_facs("address_file")
#'
#' @export
import_tm_facs <- function(address){
	start_time <- Sys.time()
	address_list <- readLines(address)
	data_list <- list()
	for (i in 1:length(address_list)){
		address_tissue <- address_list[i]
		import_tissue <- importdata(address_tissue)
		tmp <- strsplit(address_tissue, "_cnt")
		tmp1 <- strsplit(tmp[[1]], "/")
		tissue_name <- tmp1[[1]][length(tmp1[[1]])]
		print(tissue_name)
		data_list[[i]] <- import_tissue
		names(data_list)[[i]] <- tissue_name
	}
	end_time <- Sys.time()
	proc_time <- end_time - start_time
	print(proc_time)
	return (data_list)
}

#' import_sparse
#'
#' "importdata" for TM (or multiple "sparse" matrix files)
#'
#' @param address scdata address list (contains only the prefix of the file: Heart, Lung ...)
#' @param suffix_file suffix for rowname_file, cellname_file, and sparse matrix
#' \itemize{
#' \item format of the suffix_file (should be same order as shown below)
#' \item _cells.tsv
#' \item _row.tsv
#' \item _count.mtx
#' }
#' @param row_opt indicates the column number that sholud be used for the gene(peak) name
#'
#' @return list of matrix (names: each tissue or filename)
#'
#' @examples
#' tm_sparse <- import_sparse("address_file", "suffix_file")
#'
#' @export
import_sparse <- function(address, suffix_file, row_opt = 1){
	start_time <- Sys.time()
	suffix_list <- readLines(suffix_file)
	address_list <- readLines(address)
	data_list <- list()
	for (i in 1:length(address_list)){
		address_tissue <- address_list[i]
		tmp_cell <- paste(address_tissue, suffix_list[1], sep = "")
		tmp_row <- paste(address_tissue, suffix_list[2], sep = "")
		tmp_mat <- paste(address_tissue, suffix_list[3], sep = "")
		import_tissue <- importdata(tmp_mat, opt = "sparse", colname_file = tmp_cell, rowname_file = tmp_row, row_opt = row_opt)
		tmp <- strsplit(address_tissue, "/")
		tissue_name <- tmp[[1]][length(tmp[[1]])]
		print(tissue_name)
		data_list[[i]] <- import_tissue
		names(data_list)[[i]] <- tissue_name
	}
	end_time <- Sys.time()
	proc_time <- end_time - start_time
	print(proc_time)
	return(data_list)
}

#' sparse_tm_facs
#'
#' "save_sparse_data" for multiple "raw" matrix object 
#'
#' @param tm_facs list variable of importdata (same as import_tm_facs or import_sparse or list of importdata)
#'
#' @return save files in sparse format (see save_sparse_data)
#'
#' @examples
#' tm_facs <- import_tm_facs("address_file")
#' sparse_tm_facs(tm_facs)
#'
#' @export
sparse_tm_facs <- function(tm_facs){
	start_time <- Sys.time()
	for (i in 1:length(tm_facs)){
		tissue_name <- names(tm_facs)[i]
		print(tissue_name)
		save_sparse_data(tm_facs[[i]], tissue_name)
	}
	end_time <- Sys.time()
	proc_time <- end_time - start_time
	print(proc_time)	

}

#' make_input_tm_facs
#'
#' "make_input" for result of import_tm_facs
#'
#' @param data_list result of import_tm_facs or import_sparse
#' @param gene_usage option for gene/peak_location filtering to make hadamard matrix
#' \itemize{
#' \item FALSE: gene/peak_location filtering by discarding low expressing cell ratio genes/peak_locations
#' \item certain number (ex: 2048): the number of gene to use for training (must be suitable number to generate hadamard matrix)
#' }
#' @param Train TRUE if it is data for training
#' \itemize{
#' \item TRUE: adjust min_thr and max_thr. In addition, make the gene/peak_location count for hadamard matrix
#' \item FALSE: no gene/peak_location filtering
#' }
#' @param ncore multi-thread
#'
#' @return list of "make_input" results from each tissue (or sample)
#'
#' @examples
#' tm_facs <- import_tm_facs("address")
#' tm_facs_input <- make_input_tm_facs(tm_facs)
#'
#' @export
make_input_tm_facs <- function(data_list, gene_usage = FALSE, Train = TRUE, ncore = 1){
	registerDoParallel(ncore)
	start_time <- Sys.time()
	tissue_name <- names(data_list)
	tm_facs_input <- foreach (i = 1:length(data_list)) %dopar% {
		tissue_name <- names(data_list)[i]
		print(tissue_name)
		data <- data_list[[i]]
		d <- make_input(data, train = Train, gene_usage = gene_usage)
	}

	names(tm_facs_input) <- tissue_name
	end_time <- Sys.time()
	proc_time <- end_time - start_time
	print(proc_time)
	registerDoParallel(1)
	return (tm_facs_input)

}

#' tissue_merge
#'
#' merge same tissue type from tissue name (internal function for where tissue name should be merged)
#'
#' @param tissue input tissue name divided by "_"
#' \itemize{
#' \item tissue_name: Heart_1, Heart_2, Lung_1, Lung_2 ...
#' \item BrainNonMyeloid: Brain
#' \item BrainMyeloid: Brain
#' }
#'
#' @return collapsed tissue name
#'
#' @examples
#' new_tissue <- tissue_merge("Heart_1")
#'
#' @export
tissue_merge <- function(tissue){
	if (tissue == "BrainNonMyeloid"){
		tissue <- "Brain"
	}
	else if (tissue == "BrainMyeloid"){
		tissue <- "Brain"
	}
	tissue <- strsplit(tissue, "_")[[1]][1]
	return (tissue)
}


#' make_tm_facs_signature_set
#'
#' creating signature set from result of make_input_tm_facs
#'
#' @param tm_facs_input result of make_input_tm_facs
#' @param s side-infomration (julia object)
#' @param annot cell_ontology (julia object)
#' @param data_name name for the signature
#' @param sdatatype type of side-information (cell_ontology, ATAC-seq, or user-defined term)
#' @return signature set
#'
#' @examples
#' tm_facs_input <- make_input_tm_facs(tm_facs)
#' w2v <- make_w2v("Pubmed.index")
#' cell_ontology <- make_cellontology("cell_ontogloy_file")
#' stopword <- make_stopwords(stopwords)
#' side_info <- make_w2v_sideinfo(w2v, cell_ontology, stopword)
#' sig_set <- make_tm_facs_signature_set(tm_facs_input, side_info, cell_ontology, "user", "cell_ontology")
#'
#' @export
make_tm_facs_signature_set <- function(tm_facs_input, s, annot, data_name, sdatatype){
	start_time <- Sys.time()
	head = TRUE
	for(i in 1:length(tm_facs_input)){
		tissue_name <- names(tm_facs_input)[i]
		print(tissue_name)
		tissue_name <- tissue_merge(tissue_name)		
		d <- tm_facs_input[[i]]
		d_sig <- create_signature(data_name, tissue_name, s, d[[1]], d[[2]], d[[3]], annot, sdatatype, Stopwords = JuliaObject(""))
		if (head){
			tm_facs_signature_set <- make_signature_set(d_sig)
			head = FALSE
			next
		}
		tm_facs_signature_set <- add_signature_set(d_sig, tm_facs_signature_set)
	}
	end_time <- Sys.time()
	proc_time <- end_time - start_time
	print(proc_time)
	return (tm_facs_signature_set)

}

#' make_tm_facs_attribute_set
#'
#' creating attribute set from the result of make_input_tm_facs
#'
#' @param tm_facs_input the result of make_input_tm_facs
#' @param s side-information (jula object)
#' @param annot cell ontology (julia object)
#' @param mincell TRUE if user wants to give minimum cell type thresholding
#' @param data_name name for the attribute
#' @param sdatatype type of side-information (cell_ontology, ATAC-seq, or user-defined term)
#' @param ngene a number of vector size after JL transformation
#' @param fdrthr fdr threshold for calculating similarity score during training
#' @param gamma paramter for ZSL
#' @param lambda parameter for ZSL
#'
#' @return attribute set
#'
#' @examples
#' tm_facs_input <- make_input_tm_facs(tm_facs)
#' w2v <- make_w2v("Pubmed.index")
#' cell_ontology <- make_cellontology("cell_ontogloy_file")
#' stopword <- make_stopwords(stopwords)
#' side_info <- make_w2v_sideinfo(w2v, cell_ontology, stopword)
#' sig_set <- make_tm_facs_attribute_set(tm_facs_input, side_info, cell_ontology, FALSE, "user", "cell_ontology")
#'
#' @export
make_tm_facs_attribute_set <- function(tm_facs_input, s, annot, mincell, data_name, sdatatype, ngene = 64, fdrthr = 0.1, lambda = 1, gamma = 1){
	start_time <- Sys.time()
	head = TRUE
	for (i in 1:length(tm_facs_input)){
		tissue_name <- names(tm_facs_input)[i]
		print(tissue_name)
		tissue_name <- tissue_merge(tissue_name)
		d <- tm_facs_input[[i]]
		
		if (mincell){
			d_att <- create_attribute(data_name, tissue_name, s, d[[1]], d[[2]], d[[3]], annot, sdatatype, Stopwords = JuliaObject(""), Ngene = ngene, Fdrthr = fdrthr,Gamma = gamma, Lambda = lambda)		
		}
		else {
			d_att <- create_attribute(data_name, tissue_name, s, d[[1]], d[[2]], d[[3]], annot, sdatatype, Stopwords = JuliaObject(""), Mincellcount = 0, Mincellfrac = 0, Ngene = ngene, Fdrthr = fdrthr,Gamma = gamma, Lambda = lambda)	
		}
		if (head){
			tm_facs_att_set <- make_attribute_set(d_att)
			head = FALSE
			next
		}
		tm_facs_att_set <- add_attribute_set(d_att, tm_facs_att_set)
	}
	end_time <- Sys.time()
	proc_time <- end_time - start_time
	print(proc_time)
	return (tm_facs_att_set)
}

#' leave_one_val
#'
#' leave one tissue validation of TM (the parameter for side-information of att_set and sig_set must be same)
#'
#' @param att_set attribute set
#' @param sig_set signature set
#'
#' @return 6 results: conf_list, conf_ratio_list, val_list, tissue_list, f1_list, unassign_list
#' \itemize{
#' \item conf_list: confusion matrix
#' \item conf_ratio_list: ratio_list (confusion_ratio_matrix)
#' \item val_list: validation matrix
#' \item tissue_list
#' \item f1_list: f1 score
#' \item unassign_list
#' }
#' - each result is composed of sample in TM (ex: res$conf_list[1]: confusion matrix of Aorta tissue)
#'
#' @examples
#' res <- leave_one_val(att_set, sig_set)
#'
#' @export
leave_one_val <-function(att_set, sig_set){
	start_time <- Sys.time()
	res_list <- list()
	ratio_list <- list()
	val_list <- list()
	tissue_list <- c()
	f1_list <- c()
	unassign_list <- c()
	for (i in 1:length(att_set)){
		i1 <- julia_call("Int64", 1, need_return = "Julia")
		new_iscell <- julia_call("leave_one", att_set, i1, need_return = "Julia")

		tmp_tissue <- field(sig_set[i], "tissue")
		print(tmp_tissue)
		tmp_res <- predict_cell(new_iscell, sig_set[i])
		
		tmp_ratio <- confusion_ratio(tmp_res)
		tmp_val <- validation(tmp_ratio)		

		res_list[[i]] <- tmp_res
		ratio_list[[i]] <- tmp_ratio
		val_list[[i]] <- tmp_val
		tissue_list <- c(tissue_list, tmp_tissue)

		tmp_val_sample <- validation_sample(tmp_res)
		f1_score <- tmp_val_sample[[1]]
		unassign <- tmp_val_sample[[2]]

		f1_list <- c(f1_list, f1_score)
		unassign_list <- c(unassign_list, unassign)

	}

	names(f1_list) <- rep("F1", length(f1_list))

	end_time <- Sys.time()
	proc_time <- end_time - start_time
	print(proc_time)
	total_res_list <- list(res_list, ratio_list, val_list, tissue_list, f1_list, unassign_list)
	names(total_res_list) <- c("conf_list", "conf_ratio_list", "val_list", "tissue_list", "f1_list", "unassign_list")
	return (total_res_list)
}

#' validation_sample
#'
#' calculate F1 scores and unassign ratio
#'
#' @param conf confusion matrix from predict_cell
#'
#' @return f1_score, unassign
#' \itemize{
#' \item f1_score: list of F1 scores for each sample
#' \item unassign: list of unassign ratio for each sample
#' }
#'
#' @examples
#' conf <- predict_cell(att_set, query_signature)
#' res <- validation_sample(conf)
#'
#' @export
validation_sample <- function(conf){
	tmp_ratio <- confusion_ratio(conf)
	tmp_val <- validation(tmp_ratio)
	
	total_true <- sum(conf)
	total_positive <- sum(conf[,2:dim(conf)[2]])
	tmp_tp <- colSums(validation(conf))[1]
	
	precision <- tmp_tp / total_positive
	recall <- tmp_tp / total_true
	if (precision == 0 & recall == 0){
		f1_score <- 0
	}
	else {
		f1_score <- (2*precision*recall) / (precision + recall)
	}

	tmp_unassign <- colSums(validation(conf))[2]
	unassign <- tmp_unassign / total_true
	
	res_list <- list(f1_score, unassign)
	names(res_list) <- c("f1_score", "unassign")	
	return (list(f1_score, unassign))

}


#' tm_facs_logistic
#'
#' prediction by logistic regression for multiple data
#'
#' @param att_set attribute set
#' @param meta_model logistic regression model from create_meta_index_av
#' @param sig_set signature set
#' @param Tissue_opt TRUE if user wants to find tissue-origin
#' @param celltype_specific TRUE if user wants to calculate the weight fo each tissue with only matched cell type. Otherwise, weight for each tissue will be calculated with whole celltypes (deault: FALSE; tissue_opt should be TRUE to use this parameter)
#'
#' @return 6 results: conf_list, conf_ratio_list, val_list, tissue_list, f1_list, unassign_list
#' \itemize{
#' \item conf_list: confusion matrix
#' \item conf_ratio_list: ratio_list (confusion_ratio_matrix)
#' \item val_list: validation matrix
#' \item tissue_list
#' \item f1_list: f1 score
#' \item unassign_list
#' }
#' - each result is composed of sample in TM (ex: res$conf_list[1]: confusion matrix of Aorta tissue)
#'
#' @examples
#' meta_model <- create_meta_index_av(att_set1, sig_set1)
#' res <- tm_facs_logistic(att_set1, meta_model, sig_set2)
#'
#' @export
tm_facs_logistic <- function(att_set, meta_model, sig_set, Tissue_opt = FALSE, celltype_specific = FALSE){
	start_time <- Sys.time()
	res_list = list()
	ratio_list = list()
	val_list = list()
	tissue_list = list()
	f1_list = c()
	unassign_list = c()



	for (i in 1:length(sig_set)){
		tmp_tissue <- field(sig_set[i], "tissue")
		print(tmp_tissue)

		tmp_res <- logistic_prediction_cell(att_set, meta_model, sig_set[i], tissue_opt = Tissue_opt, celltype_specific = celltype_specific)
		tmp_res <- as.matrix(tmp_res[[1]])
		tmp_ratio <- confusion_ratio(tmp_res)
		tmp_val <- validation(tmp_ratio)	
	
		res_list[[i]] <- tmp_res
		ratio_list[[i]] <- tmp_ratio
		val_list[[i]] <- tmp_val
		tissue_list <- c(tissue_list, tmp_tissue)

		tmp_val_sample <- validation_sample(tmp_res)
		f1_score <- tmp_val_sample[[1]]
		unassign <- tmp_val_sample[[2]]

		f1_list <- c(f1_list, f1_score)
		unassign_list <- c(unassign_list, unassign)
	}

	f1_name <- rep("F1", length(f1_list))
	names(f1_list) <- f1_name
	end_time <- Sys.time()
	proc_time <- end_time - start_time
	print(proc_time)
	total_res_list <- list(res_list, ratio_list, val_list, tissue_list, f1_list, unassign_list)
    names(total_res_list) <- c("conf_list", "conf_ratio_list", "val_list", "tissue_list", "f1_list", "unassign_list")
    return (total_res_list)
}

#tm_facs: result of data_list from import_tm_facs
#k: ratio between ZSL(1) and logistic regression modeling (k) (k=1 -> ZSL:glm = 1:1, k=2 -> ZSL:glm = 1:2)
#during logistic regression modeling and test step -> data won't have any gene filtering
#' Title
#'
#' Description
#'
#' @param tm_facs result of import_tm_facs
#'
#' @return att_set, meta_set, logistic_result
#' \itemize{
#' \item att_set: attribute_set
#' \item meta_set: result of logistic regression (same as create_meta_index_av)
#' \item logistic_result: result of logistic regression prediction (same a leave_one_val)
#' \item n_fold: all the results needs n-fold index (ex: res$att_set$1_fold: att_set)
#' }
#'
#' @examples
#' tm_facs <- input_tm_facs("address")
#' w2v <- make_w2v("Pubmed.index")
#' cell_ontology <- make_cellontology("cell_ontogloy_file")
#' stopword <- make_stopwords(stopwords)
#' side_info <- make_w2v_sideinfo(w2v, cell_ontology, stopword)
#' res <- n_fold_logistic(tm_facs, 2048, side_info, cell_ontology, FALSE, 64, "user", "cell_ontology")
#' 
#' @export
n_fold_logistic <- function(tm_facs, Gene_usage, s, annot, mincell, Ngene, data_name, sdatatype, fdrthr = 0.1, n = 5, k = 1, tissue_opt = FALSE, Celltype_specific = FALSE, ncore = 1){
	tm_att_result<- list()
	tm_meta_result <- list()
	tm_logistic_result <- list()
	for (n_index in 5:n){
		paste("\n\n", n_index, "/",  n, "-fold change test round\n\n", sep = "") 
	
		start_time <- Sys.time()
		val_list <- list()
		zsl_list <- list()
		glm_list <- list()

		for (i in 1:length(tm_facs)){
		
			tissue_name <- names(tm_facs)[i]
			print(tissue_name)
			tmp_raw <- tm_facs[[i]]
			tmp_cells <- colnames(tmp_raw)
			tmp_raw <- tmp_raw[,which(colnames(tmp_raw) != "NA" & colnames(tmp_raw) != "unknown")]


			cell_number <- length(colnames(tmp_raw))
			random_cell_index <- sample(1:cell_number, cell_number, replace = FALSE)

			tmp_raw <- tmp_raw[,random_cell_index]
			tmp_cells <- tmp_cells[which(tmp_cells != "NA" & tmp_cells != "unknown")]
			tmp_cells <- tmp_cells[random_cell_index]
			tmp_celltype <- sort(unique(tmp_cells))
			val_cell_index_total <- c()
			zsl_train_index_total <- c()
			glm_train_index_total <- c()

			for (j in 1:length(tmp_celltype)){
				j_tmp_celltype <- tmp_celltype[j]
				tmp_cell_index <- which(tmp_cells == j_tmp_celltype)

				
				tmp_size <- round(length(tmp_cell_index) / n)


				if (n_index == n){
					 val_cell_index <- tmp_cell_index[(1 + (tmp_size * (n_index-1))):length(tmp_cell_index)]
				}
				else{
					val_cell_index <- tmp_cell_index[(1 + (tmp_size * (n_index-1))):(tmp_size * (n_index))]
				}
				total_train_index <-setdiff(tmp_cell_index, val_cell_index)
				tmp_size_zsl <- round(length(total_train_index) / (k+1))
				zsl_train_index <- total_train_index[1: tmp_size_zsl]
				glm_train_index <- total_train_index[(tmp_size_zsl + 1):length(total_train_index)]
				
				val_cell_index_total <- c(val_cell_index_total, val_cell_index)
				zsl_train_index_total <- c(zsl_train_index_total, zsl_train_index)
				glm_train_index_total <- c(glm_train_index_total, glm_train_index)

			}
			val_raw <- tmp_raw[,val_cell_index_total]
			zsl_raw <- tmp_raw[,zsl_train_index_total]
			glm_raw <- tmp_raw[,glm_train_index_total]
		
			cat("Make input for 'validation set', 'zsl set', and 'glm set'\n")				
			val_input <- make_input(val_raw)
			zsl_input <- make_input(zsl_raw, train = TRUE, gene_usage = Gene_usage)
			glm_input <- make_input(glm_raw)
			
			val_list[[i]] <- val_input
			names(val_list)[i] <- tissue_name
			zsl_list[[i]] <- zsl_input
			names(zsl_list)[i] <- tissue_name
			glm_list[[i]] <- glm_input
			names(glm_list)[i] <- tissue_name	
								
		}
		cat("\n\nSignature for GLM\n")
		tm_sig_glm_set <- make_tm_facs_signature_set(glm_list, s, annot, data_name, sdatatype)
		cat("\nSignature for Valdiation")
		tm_sig_val_set <- make_tm_facs_signature_set(val_list, s, annot, data_name, sdatatype)
		tm_att_set <- make_tm_facs_attribute_set(zsl_list, s, annot, mincell, data_name, sdatatype, ngene = Ngene, fdrthr = fdrthr)

		meta <- create_meta_index_av(tm_att_set, tm_sig_glm_set, ncore = ncore)
		logistic_res <- tm_facs_logistic(tm_att_set, meta, tm_sig_val_set, Tissue_opt = tissue_opt, celltype_specific = Celltype_specific)
		
		tm_att_result[[n_index]] <- tm_att_set
		names(tm_att_result)[n_index] <- paste(n_index, "fold", sep = "_")
		tm_meta_result[[n_index]] <- meta
		names(tm_meta_result)[n_index] <- paste(n_index, "fold", sep = "_")
		tm_logistic_result[[n_index]] <- logistic_res
		names(tm_logistic_result)[n_index] <- paste(n_index, "fold", sep = "_")

		end_time <- Sys.time()
		proc_time <- end_time - start_time
		print(proc_time)
	}
	
	total_res_list <- list(tm_att_result, tm_meta_result, tm_logistic_result)
	names(total_res_list) <- c("att_set", "meta_set", "logistic_result")	
	return (total_res_list)
}



