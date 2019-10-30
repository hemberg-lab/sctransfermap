
mutable struct IndexData
	name::String
	tissue::String
	xdatatype::String
	sdatatype::String
	celltypes::Array{String, 1}
	xnames::Array{String, 1}
	snames::Array{String, 1}
	v::Array{Float64, 2}
	fjltx::Array{Float64, 2}
	fjlts::Array{Float64, 2}
	simthres::Float64
end

mutable struct MappingData
	name::String
	tissue::String
	xdatatype::String
	sdatatype::String
	x::Array{Float64, 2}
	celltypes::Array{String, 1}
	xnames::Array{String, 1}
	snames::Array{String, 1}
	s::Array{Float64, 2}
	stypenames::Array{String, 1}
end

##############################################################
function load_emb(WordIndexFile)
	  res = load_embeddings(Word2Vec, WordIndexFile);
	  return res;
end

function IndexList(index1, iscell = IndexData[])
	push!(iscell, index1);
	return iscell;
end

function MappingList(map1, mscell = MappingData[])
	push!(mscell, map1);
	return mscell;
end

function leave_one(iscell, index)
	tmp_iscell = copy(iscell)
	deleteat!(tmp_iscell, index)
	return tmp_iscell
end

function Geneoverlap(overlap_cnt::Int64, train_cnt::Int64, overlap_thr::Float64)
	#overlap_cnt: length(indexgeneinds)
	#train_cnt: length(indexes[k].xnames)
	#overlap_thr: overlap_thr
	flag = overlap_cnt / train_cnt;
	if flag >= overlap_thr
		return false
	else
		return true
	end
	#false: pass the threshold: just for code readability
end

function make_sideinfo(w2v, annotation, stopwords, CALC)
	#CALC: true -> mean, false -> max for w2v_sentences score
	s = word2vec_tissue(w2v, vec(annotation[:,2],), stopwords, CALC);
	return s
end

###############################################################
# ATAC-seq data processing #

function parse_peak_locations(peaklocs::Array{String, 1}, sep::String="_")
	peaklocs = hcat(convert.(Array{String, 1}, split.(peaklocs, sep))...);
	pos = parse.(Int, peaklocs[2:3,:]);
	return peaklocs[1,:], pos;
end

function find_overlapping_peak_inds(querychrs::Array{String, 1}, querypeaks::Array{Int64, 2}, refchrs::Array{String, 1}, refpeaks::Array{Int64, 2})
#Use this method to identify the peaks that are found in both and that can be used for the projection
	queryoverlapinds = Int64[];
	refoverlapinds = Int64[];
	querychrsunique = sort(unique(querychrs));
	for i in 1:length(querychrsunique)
		querychrinds = findall(querychrsunique[i].==querychrs);
		refchrinds = findall(querychrsunique[i].==refchrs);
		if (length(querychrinds)>0) & (length(refchrinds)>0)
			#sort both query and ref so that we only need to traverse once
			queryposorder = sortperm(querypeaks[1,querychrinds]);
			querypossorted = querypeaks[:,querychrinds[queryposorder]];
			refposorder = sortperm(refpeaks[1,refchrinds]);
			refpossorted = refpeaks[:,refchrinds[refposorder]];
			iq = 1;
			ir = 1;
			while iq<=length(queryposorder)
				while (ir<length(refchrinds)) & (refpossorted[1,ir]<=querypossorted[2,iq])
					if ((refpossorted[1,ir]<=querypossorted[1,iq]) & (refpossorted[2,ir]>=querypossorted[1,iq])) | ((refpossorted[1,ir]<=querypossorted[2,iq]) & (refpossorted[2,ir]>=querypossorted[2,iq])) | ((refpossorted[1,ir]>=querypossorted[1,iq]) & (refpossorted[1,ir]<=querypossorted[1,iq])) | ((refpossorted[2,ir]>=querypossorted[2,iq]) & (refpossorted[2,ir]<=querypossorted[2,iq]))
						push!(queryoverlapinds, querychrinds[queryposorder[iq]])
						push!(refoverlapinds, refchrinds[refposorder[ir]])
					end
					ir += 1;
				end
				iq += 1;
			end
		end
	end
	return queryoverlapinds, refoverlapinds;
end

#####################################################################

#make [cell, celltype] matrix with (cell == celltype: 1, cell != cetype: -1)
function get_assignments(cells::Array{String, 1}; celltypesunique::Array{String, 1}=String[])
	if isempty(celltypesunique)
		celltypesunique = sort(unique(cells));
	end
	ncells = length(cells);
	assignment = -ones(Int, ncells, length(celltypesunique));
	for i in 1:length(celltypesunique)
		cellinds = vec(findall(celltypesunique[i].==cells[1:ncells]));
		assignment[cellinds, i] = vec(ones(length(cellinds), 1));
	end
	return assignment;
end


# extract normalized centroids (c-min/max-min: min, max: celltype, respectively) for celltypes from JL-transformed matrix by collecting median of values matched with target celltype
# making new side information from matrix (cell and gene exp or atac seq peak count ...)
function generate_sideinformation(s::Array{Float64, 2}, cells::Array{String, 1}, sdim::Int64)
	celltypesunique = sort(unique(cells));
	
	fjlts = FastJLTransformation(size(s,1), sdim);

	centroids = zeros(length(celltypesunique), size(fjlts, 1));
	for i in 1:length(celltypesunique)
		centroids[i,:] = median(fjlts*s[:, findall(cells.==celltypesunique[i])], dims=2);
		#scale the centroids so that they are in the interval [0, 1]
		centroids[i, :] = (centroids[i, :] .- minimum(centroids[i, :]))./(maximum(centroids[i, :]) - minimum(centroids[i, :]));
	end
	return centroids;
end

# excludeinds: exclude combination (celltype A,B): extract index from total cell, keepinds, index for remaining celltypes
function get_excluded_celltypeinds(celltypes::Array{String, 1}, excludetypeinds; celltypesunique::Array{String, 1}=String[])
	if isempty(celltypesunique)
		celltypesunique = sort(unique(celltypes));
	end
	excludeinds = union(findall(celltypes.==celltypesunique[excludetypeinds[1]]), findall(celltypes.==celltypesunique[excludetypeinds[2]]));
	for i in 3:length(excludetypeinds)
		excludeinds = union(excludeinds, findall(celltypes.==celltypesunique[excludetypeinds[i]]))
	end
	keepinds = setdiff(1:length(celltypes), excludeinds);
	return excludeinds, keepinds
end

################################################################

function create_signature(name::String, tissue::String, w2v, w2v_index, x::Array{Float64, 2}, xnames::Array{String, 1}, cells::Array{String, 1}, annotation::Array{String, 2}, sdatatype::String, stopwords, xdatatype::String)
	if w2v_index
		s = word2vec_tissue(w2v, vec(annotation[:,2]), stopwords);
	else
		s = w2v
	end
	snames = sort(annotation[:,1]);
	#sdatatype: cell_ontology filename, or atac
	return MappingData(name, tissue, xdatatype, sdatatype, x, cells, xnames, snames, s, snames)
end

function create_attribute(name::String, tissue::String, w2v, w2v_index, x::Array{Float64, 2}, xnames::Array{String, 1}, cells::Array{String, 1}, annotation::Array{String,2}, sdatatype::String, xdatatype::String, ngenes::Int64, stopwords, fdrthreshold::Float64, gamma::Float64, lambda::Float64, xpeaks::Bool, sep::String, mincellfrac::Float64, mincellcount::Int64, simthres)
	if w2v_index
		s = word2vec_tissue(w2v, vec(annotation[:,2]), stopwords);
	else
		s = w2v
	end

	snames = sort(annotation[:,1]);
	scelltypeunique = sort(unique(snames));
	celltypesunique = sort(unique(cells));
	
	#remove the cell types that are at low frequency in the reference
	minc = maximum([mincellfrac*size(x,2), mincellcount]);	
	tmp_count = Int64[];
	for i in 1:length(celltypesunique)
		push!(tmp_count, length(findall(celltypesunique[i].==cells)));
	end
	celltypesunique = celltypesunique[findall(tmp_count.>minc)];

#intersection of celltypes between x matrix and s matrix
	tmpcelltype_index = indexin(scelltypeunique, celltypesunique);
	tmpcelltype_index = tmpcelltype_index[findall(tmpcelltype_index.!=nothing)];
	celltypesunique = celltypesunique[tmpcelltype_index]

	#celltypesunique: filtered celltypesunique;

	xindskeep = Int64[];
	for i in 1:length(celltypesunique)
		append!(xindskeep, findall(celltypesunique[i].==cells));
	end
	sort!(xindskeep);
	x = x[:,xindskeep];
	cells = cells[xindskeep];
	#x: filtered x, cells: filtered cells
	
#filtering celltype from s matrix
	sind = Int64[];
	for i in 1:length(celltypesunique)
		append!(sind, findall(celltypesunique[i].==snames));
	end
	sort!(sind);
	snames = snames[sind];
	s = s[sind,:];
	#snames: filtered snames, s: filtered s

	assignment = get_assignments(cells, celltypesunique = celltypesunique);

	fjltx = FastJLTransformation(size(x, 1), ngenes);
	println("Creating index ")
	if (fdrthreshold>0) & (length(celltypesunique)>4) #find out what threshold is required 
		simthres = estimate_similarity_threshold(s, x, fjltx, assignment, cells, celltypesunique, fdrthreshold, xnames, snames, xpeaks=xpeaks, sep=sep, fjlts=Matrix{Float64}(I, 0, 0))
	elseif (fdrthreshold>0) & (length(celltypesunique)<=4)
		println("Unable to calculate threshold as there are only " * string(length(celltypesunique)) * " cell-types present.")
	end
	v = create_index_cluster(s, fjltx*x, assignment);
	return IndexData(name, tissue, xdatatype, sdatatype, celltypesunique, xnames, snames, v, fjltx, Matrix{Float64}(I, 0, 0), simthres);
end


function create_index_cluster(s::Array{Float64, 2}, x::Array{Float64, 2}, y::Array{Int64, 2}; gamma::Float64=1.0, lambda::Float64=1.0)
	#Use this for the ESZSL (the equivalent of scmap-cluster. x is the training data, for example ATAC-seq data while s is the data from the other mode (e.g. RNA-seq data) while y are the cell-types. Here, s only the median of each gene for the clusters. s is a gene by cluster matrix while x is a peak by cell matrix and y is cell by cluster matrix
	X = x*x'+gamma*Matrix{Float64}(I, size(x, 1), size(x, 1));
	S = s*s' + lambda*Matrix{Float64}(I, size(s, 1), size(s, 1));
	Y = X \ x;
	Z = S' \ s;
	v = Y*y*Z;
	return v
end

				
function create_meta_index(indexes::Array{IndexData, 1}, mappingdata::MappingData, gn_overlap_thr::Float64)
	#Build a logistic regression model to combine the outputs of the several classifiers. This method is intended to be used for cross-modality problems. The idea is to train the model using the side information from a different training dataset to emulate the situation where we get entirely new cell type descriptions. This should allow us to be more confident about which of the indexes will produce the most reliable outcome

	scelltypes = mappingdata.snames;



	ms = [] #Array{Array{String, 1}, 1};
	meta_tissue = []
	for j in 1:length(indexes) #find the prediction for each index
		m = [];
		mapgeneinds = indexin(indexes[j].xnames, mappingdata.xnames);
		mapgeneinds = unique(mapgeneinds[findall(mapgeneinds.!=nothing)]);
		indexgeneinds = indexin(mappingdata.xnames, indexes[j].xnames);
		indexgeneinds = unique(indexgeneinds[findall(indexgeneinds.!=nothing)]);
		gn_overlap_flag = Geneoverlap(length(indexgeneinds), length(indexes[j].xnames), gn_overlap_thr)
		if gn_overlap_flag
			continue
		#false: overlap
		end

		

		#need to add an additional category to allow for unassigned as well when using the meta index
		for i in 1:size(mappingdata.x, 2)
			push!(m, prediction_zsl(indexes[j].fjltx[:,indexgeneinds]*mappingdata.x[mapgeneinds, i], indexes[j].v, convert(Array{Float64, 2}, mappingdata.s')));
		end
		push!(meta_tissue, indexes[j].tissue)
		push!(ms, m);
	end
	return ms, meta_tissue, scelltypes

end

# difference btw estimate_similarity_threshold_centroid: calculate centroid with only keepinds of s matrix
function estimate_similarity_threshold(s::Array{Float64, 2}, x::Array{Float64, 2}, fjltx::Array{Float64, 2}, assignment::Array{Int64, 2}, celltypes::Array{String, 1}, celltypesunique::Array{String, 1}, fdrthreshold::Float64, xnames::Array{String, 1}, snames::Array{String, 1}; nexcludetypes::Int64=2, xpeaks::Bool=false, sep::String="_", fjlts::Array{Float64, 2}=Matrix{Float64}(I, 0, 0))

							 
	#enumerate the combinations
	testsets = collect(combinations(1:length(celltypesunique), nexcludetypes));
	thres = [];
	thres_index = [];

	for i in 1:length(testsets)
		excludetypeinds = testsets[i];
		#println("Estimating threshold while excluding types " * celltypesunique[excludetypeinds[1]] * " and " * celltypesunique[excludetypeinds[2]])
		keeptypeinds = setdiff(1:length(celltypesunique), excludetypeinds);
		excludeinds, keepinds = get_excluded_celltypeinds(celltypes, excludetypeinds; celltypesunique=celltypesunique);
		#v = create_index_cluster(centroids[keeptypeinds,:], fjltx*x[:,keepinds], assignment[keepinds, keeptypeinds]);
		centroids = s[keeptypeinds,:]
		v = create_index_cluster(centroids, fjltx*x[:,keepinds], assignment[keepinds, keeptypeinds]);

		#index = create_index("", "", assignment[keepinds, keeptypeinds], x[:,keepinds], 0, centroids, celltypes[keepinds], xnames, snames, fjltx=fjltx, fjlts=fjlts);
		index = IndexData("", "", "", "", celltypesunique[keeptypeinds], xnames, [""], v, fjltx, fjlts, -Inf);
		mappingdata = MappingData("", "", "", "", x[:,excludeinds], celltypes[excludeinds], xnames, snames, s[excludetypeinds,:], celltypesunique[excludetypeinds]);
		confusion, scores, correctscores = transfer_map_for_thr(index, mappingdata, xpeaks=xpeaks, sepindex=sep, sepmap=sep)
		#find out what threshold is required to have no more than fdrthreshold false positives

		tmp_thr_info = find_similarity_threshold(vec(scores), correctscores, fdrthreshold);
		push!(thres, tmp_thr_info[1]);
		push!(thres_index, tmp_thr_info[2]);
	end

### opt_thr: true(1) -> passed fdr_thr, false(0) -> doesn't pass fdr_thr, opt_thr : the highest accuracy (opt_thr > thr: false, opt_thr < thr: true)	
	
	total_true_cnt = sum(thres_index);
	
	n = length(thres);
	total_false_cnt = n - total_true_cnt; 
		
	
	accuracy = Float64[];
	thres_unique = unique(thres);
	#avoid draw score
	for j in 0:length(thres_unique)
		if j == 0
			acc = total_true_cnt
			push!(accuracy, acc)
			continue
		end
		tmp_indexes = findall(thres.>thres_unique[j]);
		tmp_length = length(tmp_indexes);
		if tmp_length != 0
			tp_acc = sum(thres_index[tmp_indexes]);
			fp_acc = total_false_cnt - tmp_length + tp_acc;
			acc = tp_acc + fp_acc;
		else
			acc = total_false_cnt - tmp_length
		end	
		push!(accuracy, acc);
	end
	max_acc_index =	findmax(accuracy)[2]
	if max_acc_index == 1
	#	return median(thres)
		return -Inf
	else
		return thres_unique[max_acc_index - 1]
	#	return median(thres)
	end
end




function find_similarity_threshold(scores::Array{Float64, 1}, correctscores::Array{Float64, 1}, fdrthreshold::Float64)
	n = length(scores);
	initial_t = length(correctscores)

	initial_fpr = 1- (initial_t / n)
	if initial_fpr < fdrthreshold
		return [-Inf, 1.0]
	end
	sort!(scores);
	sort!(correctscores);
	ind = 1;
	n = length(scores);

	while ind<n
		tp = length(findall(correctscores.>scores[ind]))
		fpr = (n-ind-tp)/(n-ind); #false positive rate
		if fpr.<fdrthreshold
			break;
		end
		ind += 1;
	end
	if ind == n
		return [scores[ind], 0.0];
	else
		return [scores[ind], 1.0];
	end
end
				
##################################################################

function transfer_map_for_thr(index::IndexData, mappingdata::MappingData; xpeaks::Bool=false, sepindex::String="_", sepmap::String="_", simthres::Float64=Inf)
	celltypesuniquemapping = sort(unique(mappingdata.snames));
	celltypesuniquedata = sort(unique(mappingdata.celltypes));
	confusion = zeros(Int, length(celltypesuniquedata), length(celltypesuniquemapping));
	#find the genes that are required for the index
	indexxinds = Int[];
	mapxinds = Int[];
	if (xpeaks) 
		ci, pi = parse_peak_locations(index.xnames, sep=sepindex);
		cm, pm = parse_peak_locations(mappingdata.xnames, sep=sepmap);
		indexxinds, mapxinds = find_overlapping_peak_inds(ci, pi, cm, pm);
		println("Found " * string(length(mapxinds)) * " peaks in common for the data.")
	else
		indexxinds = indexin(mappingdata.xnames, index.xnames);
		indexxinds = (indexxinds[findall(indexxinds.!=nothing)]);
		mapxinds = indexin(index.xnames, mappingdata.xnames);
		mapxinds = (mapxinds[findall(mapxinds.!=nothing)]);
	end

	centroids = convert(Array{Float64, 2}, mappingdata.s);


	scores = zeros(1, size(mappingdata.x, 2));
	correctscores = Float64[]
	for i in 1:length(celltypesuniquedata)
		exinds = findall(mappingdata.celltypes.==celltypesuniquedata[i]);
		for j in 1:length(exinds)
			tmp = predict_max_celltype(index.fjltx[:,indexxinds]*mappingdata.x[mapxinds, exinds[j]], index.v, convert(Array{Float64, 2}, centroids'); returnscore=true);
			if (tmp[1]>index.simthres) | (tmp[1]>simthres)
				confusion[i,tmp[2]] += 1;
				if tmp[2]==i
					push!(correctscores, tmp[1]);
				end
			end
			scores[exinds[j]] = tmp[1];
		end
	end
	return confusion, scores, correctscores
end

							
							
							
							
function transfer_map_consensus(indexes::Array{IndexData, 1}, mappingdata::MappingData, simthres::Float64, gn_overlap_thr::Float64)
	celltypesuniquemapping = mappingdata.snames;
	celltypesuniquedata = sort(unique(mappingdata.celltypes));
	confusion = zeros(Int, length(celltypesuniquedata), length(celltypesuniquemapping));
	#find the genes that are required for the index
	predictions = zeros(Int, size(mappingdata.x, 2), length(celltypesuniquemapping));
	for k in 1:length(indexes)
		println("Tissue from ZSL: " * string(k) * "/" * string(length(indexes))) 
		if indexes[k].simthres < Inf
			mapgeneinds = indexin(indexes[k].xnames, mappingdata.xnames);
			mapgeneinds = unique(mapgeneinds[findall(mapgeneinds.!=nothing)]);

			indexgeneinds = indexin(mappingdata.xnames, indexes[k].xnames);
			indexgeneinds = unique(indexgeneinds[findall(indexgeneinds.!=nothing)]);
			gn_overlap_flag = Geneoverlap(length(indexgeneinds), length(indexes[k].xnames), gn_overlap_thr)
			if gn_overlap_flag
				continue
			#false: overlap	
			end
			

			for i in 1:length(celltypesuniquedata)
				exinds = findall(mappingdata.celltypes.==celltypesuniquedata[i]);
				for j in 1:length(exinds)
					tmp = predict_max_celltype(indexes[k].fjltx[:,indexgeneinds]*mappingdata.x[mapgeneinds, exinds[j]], indexes[k].v, convert(Array{Float64, 2}, mappingdata.s'); returnscore=true);
					if (tmp[1]>indexes[k].simthres) | (tmp[1]>simthres)
						predictions[exinds[j], tmp[2]] += 1;
					end
				end	
			end
		end	
	end
	println("consensus mapping...")
	for i in 1:length(celltypesuniquedata)
		exinds = findall(mappingdata.celltypes.==celltypesuniquedata[i]);
		for j in 1:length(exinds)
			if sum(predictions[exinds[j],:])>0
				confusion[i,findmax(predictions[exinds[j],:])[2]] += 1;
			end
		end
	end
	return confusion;
end


function predict_max_celltype(x::Array{Float64, 1}, v::Array{Float64, 2}, s::Array{Float64, 2}; returnscore::Bool=false)
	if returnscore
		return findmax(vec((reshape(x, 1, length(x))*v*s)));
	end
	return findmax(vec((reshape(x, 1, length(x))*v*s)))[2];
end

function prediction_zsl(x::Array{Float64, 1}, v::Array{Float64, 2}, s::Array{Float64, 2})
	return vec((reshape(x, 1, length(x))*v*s))


end

##################################################################

function parse_annotation_file(filename::String)
	annotation = convert(Array{String, 2}, readdlm(filename, '|'));
	return annotation[sortperm(vec(annotation[:,1])),:];
end

function normalize_w2v_scores(v::Array{Float32, 1})
	#The ESZSL algorithm prefers attributes in the range [0, 1], so we need to re-scale
	return (v .- minimum(v))./(maximum(v) - minimum(v));
end

function normalize_w2v_scores(v::Array{Float64, 1})
	return (v .- minimum(v))./(maximum(v) - minimum(v));
end

function word2vec_tissue(w2v::Embeddings.EmbeddingTable, annotations::Array{String, 1}, stopwords::Array{String, 1}, calc::Bool=true)
	#calc: true -> mean value or w2v_sentence, false -> max value
	#calculate scores for all cell type descriptions as the mean of their word vectors. Re-scale each vector to the interval [0, 1]
	d = size(w2v.embeddings, 1);
	vecs = zeros(length(annotations), d);
	for i in 1:length(annotations)
		#vecs[i,:] = mean(word2vec_sentence(w2v, convert(String, annotations[i]), stopwords)[1], dims=1);
		if calc
			vecs[i,:] = normalize_w2v_scores(vec(mean(word2vec_sentence(w2v, convert(String, annotations[i]), stopwords)[1], dims=1)));
		else
			vecs[i,:] = normalize_w2v_scores(vec(maximum(word2vec_sentence(w2v, convert(String, annotations[i]), stopwords)[1], dims=1)));
		end
	end
	return vecs;
end

function word2vec_sentence(w2v::Embeddings.EmbeddingTable, desc::String, stopwords::Array{String, 1})
	#remove dots and the like
	desc = replace(replace(replace(replace(replace(desc, "." => ""), ";" => ""), "," => ""), "(" => " "), ")" => " ")
	#tokenize the sentence
	tokens = split(desc);
	#remove the stop words
	tokens = tokens[findall(indexin(map(lowercasefirst, tokens), stopwords).==nothing)];
	#find the word2vec representation for each remaining word
	vecs = Array{Float64, 2}(undef, 0, 200);
	for i in 1:length(tokens)
		ind = findfirst(w2v.vocab.=="\n" * tokens[i]);
		if ind!=nothing
			#vecs = vcat(vecs, normalize_w2v_scores(w2v.embeddings[:,ind])');
			vecs = vcat(vecs, w2v.embeddings[:,ind]');
		else #try with lower case
			ind = findfirst(w2v.vocab.=="\n" * lowercasefirst(tokens[i]));
			if ind!=nothing
				#vecs = vcat(vecs, normalize_w2v_scores(w2v.embeddings[:,ind])');
				vecs = vcat(vecs, w2v.embeddings[:,ind]');
			end
		end
	end
	return vecs, tokens
end

##################################################################


function FastJLTransformation(n::Int64, d::Int64) 
	#Calculate the fast Johnson-Lindenstrauss transform and project onto d dimensions. We ues this to speed up the distance calculation. Assume that x is a gene x cell matrix
	#Find a Hadamard matrix of suitable size. The requirement is that the factors that are not 2, form a product that is <64
	h = n;
	H = map(Float64, hadamard(h))/sqrt(h);
	D = Diagonal(vec(sign.(randn(1, h))));
	#R = eye(h)[randperm(h)[1:d], :] #
	R = SparseSubGaussianMatrix(d, h); #randn(d, h)
	return R*H*D;
end

function SparseSubGaussianMatrix(n::Int64, m::Int64)
	x = zeros(n, m)
	for i in 1:n
		for j in 1:m
			r = rand()
			if r>5/6
				x[i, j] = 1
			elseif r<1/6
				x[i, j] = -1
			end
		end
	end
	return x*3/n
end

