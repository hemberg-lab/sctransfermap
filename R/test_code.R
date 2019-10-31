#sctransfermap: test_code.R
#
#start_time <- Sys.time()
#function
#end_time <- Sys.time()
#proc_time <- end_time - start_time
#proc_time <- as.numeric(proc_time)

#object.size(variable_name)

#' Title
#'
#' "importdata" for TM
#'
#' @param address scdata address list
#'
#' @return 
#'
#' @examples
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
#brainmyeloid, brainnonmyeloid -> brain


# gene_usage: 2048 or total
#' Title
#'
#' "make_input" for result of import_tm_facs
#'
#' @param data_list result of import_tm_facs
#'
#' @return 
#'
#' @examples
#'
#' @export
make_input_tm_facs <- function(data_list, gene_usage = FALSE, Train = TRUE, ncore = 1){
	registerDoParallel(ncore)
	start_time <- Sys.time()
	tissue_name <- names(data_list)
#	tm_facs_input <- list()
#	for (i in 1:length(data_list)){
	tm_facs_input <- foreach (i = 1:length(data_list)) %dopar% {
		tissue_name <- names(data_list)[i]
		print(tissue_name)
		data <- data_list[[i]]
		d <- make_input(data, train = Train, gene_usage = gene_usage)
#		tm_facs_input[[i]] <- d
#		names(tm_facs_input)[[i]] <- tissue_name
	}

	names(tm_facs_input) <- tissue_name
	end_time <- Sys.time()
	proc_time <- end_time - start_time
	print(proc_time)
	registerDoParallel(1)
	return (tm_facs_input)

}

#' Title
#'
#' creating signature set from result of make_input_tm_facs
#'
#' @param tm_facs_input result of make_input_tm_facs
#' @param s side-infomration (julia object)
#' @param annot cell_ontology (julia object)
#'
#' @return signature set
#'
#' @examples
#'
#' @export
make_tm_facs_signature_set <- function(tm_facs_input, s, annot){
	start_time <- Sys.time()
	head = TRUE
	for(i in 1:length(tm_facs_input)){
		tissue_name <- names(tm_facs_input)[i]
		print(tissue_name)
		if (tissue_name == "BrainNonMyeloid"){
			tissue_name <- "Brain"
		}
		else if (tissue_name == "BrainMyeloid"){
			tissue_name <- "Brain"
		}
		
		d <- tm_facs_input[[i]]
		d_sig <- create_signature("Tabula_Muris_Facs", tissue_name, s, d[[1]], d[[2]], d[[3]], annot, "cellontology", Stopwords = JuliaObject(""))
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

#mincell: TRUE -> do minimum cell thresholding
#' Title
#'
#' creating attribute set from the result of make_input_tm_facs
#'
#' @param tm_facs_input the result of make_input_tm_facs
#' @param s side-information (jula object)
#' @param annot cell ontology (julia object)
#' @param mincell TRUE if user wants to give minimum cell type thresholding
#'
#' @return attribute set
#'
#' @examples
#'
#' @export
make_tm_facs_attribute_set <- function(tm_facs_input, s, annot, mincell, ngene = 64, fdrthr = 0.1, lambda = 1, gamma = 1){
	start_time <- Sys.time()
	head = TRUE
	for (i in 1:length(tm_facs_input)){
		tissue_name <- names(tm_facs_input)[i]
		print(tissue_name)
		if (tissue_name == "BrainNonMyeloid"){
			tissue_name <- "Brain"
		}
		else if (tissue_name == "BrainMyeloid"){
			tissue_name <- "Brain"
		}
		d <- tm_facs_input[[i]]
		
		if (mincell){
			d_att <- create_attribute("Tabula_Muris_Facs", tissue_name, s, d[[1]], d[[2]], d[[3]], annot, "cellontology", Stopwords = JuliaObject(""), Ngene = ngene, Fdrthr = fdrthr,Gamma = gamma, Lambda = lambda)		
		}
		else {
			d_att <- create_attribute("Tabula_Muris_Facs", tissue_name, s, d[[1]], d[[2]], d[[3]], annot, "cellontology", Stopwords = JuliaObject(""), Mincellcount = 0, Mincellfrac = 0, Ngene = ngene, Fdrthr = fdrthr,Gamma = gamma, Lambda = lambda)	
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

#' Title
#'
#' leave one tissue validation of TM (the parameter for side-information of att_set and sig_set must be same)
#'
#' @param att_set attribute set
#' @param sig_set signature set
#'
#' @return res_list, ratio_list, val_list, tissue_list, f1_list, unassign_list
#' res[[1]]: res_list (confusion matrix)
#' res[[2]]: ratio_list (confusion_ratio_matrix)
#' res[[3]]: val_list (validation matrix)
#' res[[4]]: tissue_list
#' res[[5]]: f1_list (f1 score)
#' res[[6]]: unassign_list
#' each result is composed of sample in TM (ex: res[[1]][1]: confusion matrix of Aorta tissue)
#' @examples
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
	return (list(res_list, ratio_list, val_list, tissue_list, f1_list, unassign_list))
}

#' Title
#'
#' Description
#'
#' @param description of parameters
#'
#' @return 
#'
#' @examples
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

	return (list(f1_score, unassign))

}


#modeling logistic regression is enough from CoreMethod.R
#' Title
#'
#' Description
#'
#' @param description of parameters
#'
#' @return 
#'
#' @examples
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
	return (list(res_list, ratio_list, val_list, tissue_list, f1_list, unassign_list))


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
#' @return tm_att_result, tm_meta_result, tm_logistic_result
#' res[[1]]: tm_att_result -> attribute_set
#' res[[2]]: tm_meta_result -> result of logistic regression (same as create_meta_index)
#' res[[3]]: tm_logistic_result -> result of logistic regression prediction (same a leave_one_val)
#' all the result should decide n-fold index (ex: res[[1]][[1]]: tm_att_result, res[[2]][[5]]: tm_meta_result)
#' @examples
#'
#' @export
n_fold_logsitic <- function(tm_facs, Gene_usage, s, annot, mincell, Ngene, fdrthr = 0.1, n = 5, k = 1, tissue_opt = FALSE, Celltype_specific = FALSE, ncore = 1){
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
		tm_sig_glm_set <- make_tm_facs_signature_set(glm_list, s, annot)
		cat("\nSignature for Valdiation")
		tm_sig_val_set <- make_tm_facs_signature_set(val_list, s, annot)
		tm_att_set <- make_tm_facs_attribute_set(zsl_list, s, annot, mincell, ngene = Ngene, fdrthr = fdrthr)

		meta <- create_meta_index_av(tm_att_set, tm_sig_glm_set, ncore = ncore)
		logistic_res <- tm_facs_logistic(tm_att_set, meta, tm_sig_val_set, Tissue_opt = tissue_opt, celltype_specific = Celltype_specific)
		
		tm_att_result[[n_index]] <- tm_att_set
		tm_meta_result[[n_index]] <- meta
		tm_logistic_result[[n_index]] <- logistic_res
		

		end_time <- Sys.time()
		proc_time <- end_time - start_time
		print(proc_time)
	}
	return (list(tm_att_result, tm_meta_result, tm_logistic_result))
}



