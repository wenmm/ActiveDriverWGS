

#### Regression Test

get_obs_exp = function(hyp, select_positions, dfr) {
	obs_mut = sum(dfr[select_positions, "muts"])
	predict_mut = stats::predict(hyp, type="response")[select_positions]
	exp_mut = sum(predict_mut)
	exp_boot = sapply(1:100, function(x) sum(na.rm=T, sample(predict_mut, replace=T)))
	obs_enriched = obs_mut > stats::quantile(exp_boot, 0.95)
	list(obs_mut, exp_mut, obs_enriched)
}

get_dfr_by_3n = function(this_3n, this_muts, this_tag, cols_3n) {

	# first make sure that only valid trinucleotides are included (exclude those with N)	
	this_muts = this_muts[this_muts$mut_tag %in% cols_3n,,drop=F]	
	this_3n = this_3n[names(this_3n) %in% cols_3n]
	
	# subtract mutated trinucleotides from all trinucleotides
	# we subtract positions in trinucleotide context that have been mutated at least once
	# as opposed to all counts of trinucleotides where recurrent mutations would be counted several times
	mutated_3n = table(unique(this_muts[,c("mut_pos1", "mut_tag")])[,"mut_tag"])
	this_3n[names(mutated_3n)] = this_3n[names(mutated_3n)] - mutated_3n
	# some indel situations may hit a small site (one trinucleotide) but overlap it
	# thus subtraction will cause negative numbers; these need to be flattened
	this_3n[this_3n<0] = 0
	
	# first get mutated positions by recurrence
	mut_count_by_3n = table(this_muts$mut_pos1, this_muts$mut_tag)
	trinucs_used = colnames(mut_count_by_3n)
	mut_count_by_3n = lapply(colnames(mut_count_by_3n), function(x) mut_count_by_3n[,x])
	mut_count_by_3n = lapply(mut_count_by_3n, function(x) unname(x[x>0]))
	names(mut_count_by_3n) = trinucs_used
	
	mut_dfr = do.call(rbind, lapply(names(mut_count_by_3n), 
			function(x) do.call(rbind, lapply(mut_count_by_3n[[x]], 
					function(y)  cbind(muts=y,trinuc=x) )) ))
	# there could be a site that is all mutated, no silent trinucleotides
	quiet_dfr = NULL
	if (sum(this_3n)>0) {
		quiet_dfr = cbind(muts=0, trinuc=unlist(lapply(names(this_3n), function(x) rep(x,this_3n[x]) )))
	}
	
	dfr = data.frame(rbind(mut_dfr, quiet_dfr), element=this_tag, stringsAsFactors=F)
	dfr$muts = as.numeric(dfr$muts)
	dfr
}

get_site_dfr = function(frag_id, site_id, this_site_3n, site_muts, cols_3n) {
	single_site_3n = unlist(this_site_3n[this_site_3n$frag_id==frag_id & this_site_3n$site_id==site_id,cols_3n])
	single_site_muts = site_muts[site_muts$reg_frag_id==frag_id & site_muts$reg_site_id==site_id,]

	get_dfr_by_3n(single_site_3n, single_site_muts, paste("site", frag_id, site_id, sep="_"), cols_3n)
}

regress_test = function(id, mutations_in_sites, mutations_in_elements, mutations_in_windows, 
							site_3n, element_3n, window_3n) {

	blank_result = data.frame(id, pp_site=NA, pp_element=NA, element_muts_obs=NA, element_muts_exp=NA, 
					element_enriched=NA, site_muts_obs=NA, site_muts_exp=NA, site_enriched=NA,
					selected_sites=NA, h0_df=NA, h1_df=NA, stringsAsFactors=F)

	# select mutations in this element, sites and surrounding window
	site_muts = mutations_in_sites[mutations_in_sites$reg_id==id,, drop=F]
	element_muts = mutations_in_elements[mutations_in_elements$reg_id==id,, drop=F]
	window_muts = mutations_in_windows[mutations_in_windows$reg_id==id,, drop=F]
	
	# return NA earlier if no mutations found in element
	if (nrow(element_muts)==0) {
		return(blank_result)
	}
	
	if (!is.null(site_muts[[1]])) {
		site_muts$mut_tag = toupper(gsub(">.$", "", site_muts$mut_tag))
	}
	element_muts$mut_tag = toupper(gsub(">.$", "", element_muts$mut_tag))
	window_muts$mut_tag = toupper(gsub(">.$", "", window_muts$mut_tag))

	# select trinucleotides in this element, sites, and surrounding window	
	cols_3n = c(names(Biostrings::trinucleotideFrequency(Biostrings::DNAString("AA"))), "INDEL")
	
	this_site_3n = site_3n[site_3n$id==id,,drop=F]
	this_site_3n_total = 0
	mut_site_frag_id = ""
	if (!is.null(this_site_3n[[1]])) {
		# keep only mutated sites
		this_site_3n = site_3n[site_3n$id==id,,drop=F]
		this_site_3n$site_frag_id = paste(sep=":", this_site_3n$frag_id, this_site_3n$site_id)
		mut_site_frag_id = unique(paste(sep=":", site_muts$reg_frag_id, site_muts$reg_site_id))
		this_site_3n = this_site_3n[this_site_3n$site_frag_id %in% mut_site_frag_id,,drop=F]
		this_site_3n_total = apply(this_site_3n[,cols_3n, drop=F], 2, sum)		
	}

	this_element_3n = unlist(element_3n[element_3n$id==id,cols_3n, drop=T])
	this_window_3n = unlist(window_3n[window_3n$id==id,cols_3n, drop=T])
	
	# exclude narrower muts from wider
	window_muts = window_muts[!window_muts$mut_pos1 %in% element_muts$mut_pos1,]
	element_muts = element_muts[!element_muts$mut_pos1 %in% site_muts$mut_pos1,]

	# exclude narrower trinucleotide counts from wider
	this_window_3n = this_window_3n-this_element_3n
	this_element_3n = this_element_3n-this_site_3n_total
	
	# indels might hit multiple sites or elements !!!
	# therefore keep only one copy of every individual mutation
	muts_cols = c("mut_chr", "mut_pos1", "mut_pos2", "mut_ref", 
			"mut_alt", "mut_patient", "mut_tag")
	site_muts = site_muts[!duplicated(site_muts[,muts_cols]),,drop=F]
	element_muts = element_muts[!duplicated(element_muts[,muts_cols]),,drop=F]
	window_muts = window_muts[!duplicated(window_muts[,muts_cols]),,drop=F]
	
	# remove kataegis-style multi-mutations in the same one patient
	site_muts = site_muts[!duplicated(site_muts[,"mut_patient"]),,drop=F]
	element_muts = element_muts[!duplicated(element_muts[,"mut_patient"]),,drop=F]
	element_muts = element_muts[!element_muts[,"mut_patient"] %in% site_muts[,"mut_patient"],, drop=F]
	
	# some promoters have a lot of NNNN nucleotides, skipping		
	if (sum(this_element_3n[-length(this_element_3n)]) < 25) {
		return(blank_result)
	}
	
	element_dfr = cbind(get_dfr_by_3n(this_element_3n, element_muts, "element", cols_3n), any_site=FALSE)
	window_dfr = cbind(get_dfr_by_3n(this_window_3n, window_muts, "Window", cols_3n), any_site=FALSE)

	site_dfr = NULL
	
	if (!is.null(nrow(site_muts)) && nrow(site_muts)>0) {
		site_dfr = cbind(get_dfr_by_3n(this_site_3n_total, site_muts, "element", cols_3n), any_site=TRUE)
	}
	
	merged_dfr = rbind(site_dfr, element_dfr, window_dfr)
	
	rm(site_dfr, element_dfr, window_dfr)
	gc()
	
	# now remove trinucleotides that are never mutated - their weights in the model won't matter
	never_mutated_trinucs = names(which(by(merged_dfr$muts, merged_dfr$trinuc, sum)==0))
	merged_dfr = merged_dfr[!merged_dfr$trinuc %in% never_mutated_trinucs,]	

	# if element not mutated, don't calculate further	
	obs_mut_element = sum(merged_dfr[merged_dfr$element=="element", "muts"])
	if (obs_mut_element==0) {
		return(blank_result)		
	}
	
	# if only one type of trinucleotides present then trinucleotide cannot be included as predictor
	h0_formula = stats::as.formula("muts~trinuc")
	h1_formula = stats::as.formula("muts~trinuc+element")
	if (length(unique(merged_dfr$trinuc))==1) {
		h0_formula = stats::as.formula("muts~1")
		h1_formula = stats::as.formula("muts~element")
	}
	
	h0 = stats::glm(h0_formula, family=stats::poisson, data=merged_dfr)
	h1 = stats::glm(h1_formula, family=stats::poisson, data=merged_dfr)
	chi_sq = stats::anova(h0, h1, test="Chisq")
	pp_element = chi_sq[2,5]
	h0_df = chi_sq[1,3]
	h1_df = chi_sq[2,3]

	pp_site = site_muts_obs = site_muts_exp = site_enriched = selected_sites = NA

	if (!is.null(this_site_3n[[1]]) && !is.null(site_muts[[1]]) && nrow(this_site_3n)>0 && nrow(site_muts)>0 && length(which(merged_dfr[,1]>0))>1) {

		h2_formula = stats::as.formula(ifelse(length(unique(merged_dfr$trinuc))>1, 
				"muts ~ trinuc + element + any_site", 
				"muts ~ element + any_site"))
		
		h2 = stats::glm(h2_formula, data=merged_dfr, family=stats::poisson)
		pp_site = stats::anova(h1, h2, test="Chisq")[2,5]

		site_stats = get_obs_exp(h1, merged_dfr[,"any_site"], merged_dfr) 	
		site_muts_obs = site_stats[[1]]
		site_muts_exp = site_stats[[2]]
		site_enriched = site_stats[[3]]
		
		if (!site_enriched & !is.na(pp_site) & pp_site<0.5) {
			pp_site = 1 - pp_site
		}
		
		rm(h2)
	}
	
	element_stats = get_obs_exp(h0, merged_dfr$element=="element", merged_dfr) 
	element_muts_obs = element_stats[[1]]
	element_muts_exp = element_stats[[2]]
	element_enriched = element_stats[[3]]
	
	if (!element_enriched & !is.na(pp_element) & pp_element<0.5) {
		pp_element = 1 - pp_element
	}
	
	rm(h0, h1, merged_dfr)

	data.frame(id, pp_site, pp_element, 
			element_muts_obs, element_muts_exp, element_enriched, 
			site_muts_obs, site_muts_exp, site_enriched,
			selected_sites = paste(mut_site_frag_id, collapse=";"), h0_df, h1_df,
			stringsAsFactors=F)
}


#### Prepare Elements

split_coord_fragments_in_BED = function(i, coords) {
	lens = as.numeric(strsplit(coords[i, "V11"], ',')[[1]])
	starts = as.numeric(strsplit(coords[i, "V12"], ',')[[1]]) + coords[i, "V2"]
	ends = starts+lens
	
	dfr = cbind(chr=coords[i,"V1"], starts, ends, id=coords[i,"V4"], frag_id=1:length(lens))
	dfr	
}

prepare_element_coords_from_BED = function(fname) {
	input = utils::read.delim(fname, stringsAsFactors=F, header=F)
	colnames(input) = paste0("V", 1:ncol(input))
	input$V1 = gsub("chr", "", input$V1, ignore.case = TRUE)
	input$V1 = paste0("chr", input$V1)
	coords = do.call(rbind, lapply(1:nrow(input), split_coord_fragments_in_BED, input))
	coords = data.frame(coords, stringsAsFactors=F)
	coords$starts = as.numeric(coords$starts)
	coords$ends = as.numeric(coords$ends)
	coords$frag_id = as.numeric(coords$frag_id)
	
	chr = by(coords$chr, coords$id, function(x) unique(as.character(x)), simplify=F)
	exclude_XY = names(which(sapply(chr, length)==2))
	legal_chr = paste0("chr", c(1:22, "M", "X", "Y"))
	coords = coords[!(coords$id %in% exclude_XY) & coords$chr %in% legal_chr,]
	# remove repeated elements
	coords = unique(coords)
	coords
}

coords_sanity = function(coords) {
	if (is.null(coords[[1]])) {
		return(NULL)
	}
	max_end = GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[as.character(coords$chr)]
	# to accommodate -1..+1 for trinucleotide computation
	# maximum coordinate is max-1
	coords$ends[coords$ends>(max_end-1)] = (max_end-1)[coords$ends>(max_end-1)]
	# minimum coordinate is 2
	coords$starts[coords$starts<2] = 2
	coords	
}

get_3n = function(coords) {
	
	if (is.null(coords[[1]])) {
		return(NULL)
	}

	gr_3n = GenomicRanges::GRanges(coords$chr, IRanges::IRanges(coords$starts-1, coords$ends+1), mcol=coords[,c("id", "frag_id")])	
	seqs = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, gr_3n)
	trinucs = trinucs1 = Biostrings::trinucleotideFrequency(seqs)

	# CT in main signature; translate to alternative strand if GA
	new_tag = colnames(trinucs)
	mid_tag = gsub("(.)(.)(.)", "\\2", new_tag)
	which_to_complement = which(mid_tag %in% c("G", "A"))
	which_to_keep = which(!mid_tag %in% c("G", "A"))

	trinucs_keep = trinucs[,which_to_keep,drop=F]
	trinucs_complement = trinucs[,which_to_complement,drop=F]
	# complement names of dfr to add up to primary trinucleotides
	colnames(trinucs_complement) = as.character(Biostrings::complement(Biostrings::DNAStringSet(colnames(trinucs_complement))))
	trinucs_keep[,colnames(trinucs_keep)] = 
		trinucs_keep[,colnames(trinucs_keep)] + trinucs_complement[,colnames(trinucs_keep)]
	# set secondary trinucleotides to zero
	trinucs_complement[,] = 0
	# complement names of dfr again to get back to original trinucleotides (but these are zero for now)
	colnames(trinucs_complement) = as.character(Biostrings::complement(Biostrings::DNAStringSet(colnames(trinucs_complement))))
	trinucs = cbind(trinucs_keep, trinucs_complement)

	trinucs = data.frame(trinucs, INDEL=apply(trinucs, 1, sum), id=coords$id, frag_id=coords$frag_id, stringsAsFactors=F)
	if (length(coords$site_id)==nrow(trinucs)) {
		trinucs$site_id = coords$site_id
	} else {
		trinucs$site_id = NA
	}
	
	trinucs
}

add_trinucs_by_elements = function(this_3n) {

	cols_3n = c(names(Biostrings::trinucleotideFrequency(Biostrings::DNAString("AA"))), "INDEL")
	trinucs_dfr = do.call(rbind, by(this_3n, this_3n[, "id"], function(x) apply(x[, cols_3n], 2, sum) ))
	trinucs_dfr = data.frame(trinucs_dfr, id=rownames(trinucs_dfr), stringsAsFactors=F)
	rownames(trinucs_dfr) = NULL
	trinucs_dfr
}

get_window_coords = function(coords, flank_window) {

	min_start = c(by(coords$starts, coords$id, min) - flank_window)
	max_end = c(by(coords$ends, coords$id, max) + flank_window)

	chr = by(coords$chr, coords$id, function(x) unique(as.character(x)))
	if (!all(sapply(c(min_start, max_end, chr), length)==1)) {
		stop("Error getting window coords")
	}

	data.frame(chr=c(chr), starts=min_start, ends=max_end, id=names(chr), frag_id=1, 
			stringsAsFactors=F)		
}

get_sites_per_element = function(i, element_coords, site_gr) {
	if (stats::runif(1)<0.01) gc()
	
	element_gr = GenomicRanges::GRanges(element_coords[i, "chr"], 
			IRanges::IRanges(start=element_coords[i, "starts"], end=element_coords[i, "ends"]), 
			mcol=data.frame(id=element_coords[i,"id"], frag_id=element_coords[i, "frag_id"], 
			stringsAsFactors=F))

	overlaps = GenomicRanges::findOverlaps(site_gr, element_gr)
	site_matched = site_gr[S4Vectors::queryHits(overlaps)]
	if (length(site_matched)==0) {
		return(NULL)
	}

	covr = GenomicRanges::coverage(site_matched)
	rm(overlaps, site_matched)
	
	site_islands = NA
	if (length(covr)>0) {
		site_islands = IRanges::slice(covr, lower=1)[[element_coords[i, "chr"]]]
	}
	
	# starts and ends of sites need to be adjusted to not exceed element borders
	ends = sapply(BiocGenerics::end(site_islands), min, element_coords[i, "ends"])
	starts = sapply(BiocGenerics::start(site_islands), max, element_coords[i, "starts"])

	data.frame(chr=element_coords[i, "chr"], starts, ends, id=element_coords[i,"id"], 
			frag_id=element_coords[i, "frag_id"], site_id=1:length(site_islands), 
			stringsAsFactors=F)
}

#' Preparation of elements and active sites for ActiveDriverWGS analysis with the \code{\link{find_drivers}} function
#'
#' @param elements_file path to a BED12 file of the elements to analyze
#' @param sites_file path to a BED4 file of active sites to use in the analysis
#' @param window_size numeric value for the size of the background sequence window to use on each side of the given elements. Default value is 50000.
#' @param recovery_dir path to a directory for ActiveDriverWGS to store intermediate computations in.
#'             In the case that prepare_elements does not complete, run it again with the same value for this parameter,
#'             and it will read the intermediate computations rather than compute them again. Please specify a
#'             different directory, or delete the ActiveDriverWGS data in a previously used directory, if running with
#'             different data. Default value is NA, and prepare_elements will not save any intermediate information with this value.
#' @param mc.cores numeric indicating the number of computing cores to use. Default value is 1.
#'
#' @return list containing element and active site information needed for ActiveDriverWGS analysis
#'
#' @export
prepare_elements = function(elements_file, sites_file, window_size = 50000, mc.cores = 1, recovery_dir = NA) {
	try_error = "try-error"

	element_coords_filename = paste0(recovery_dir, "/ADWGS_elements_coords_recovery_file.rsav")
	element_coords_1 = NULL
	load_result = suppressWarnings(try(load(element_coords_filename), silent = TRUE))
	if (class(load_result) == try_error) {
		element_coords_1 = prepare_element_coords_from_BED(elements_file)
		if (!is.na(recovery_dir)) {
			save(element_coords_1, file = element_coords_filename)
		}
		cat("Successfully read elements\n")
	} else {
		cat("Recovered elements\n")
	}

	element_coords = coords_sanity(element_coords_1)

	element_3n_filename = paste0(recovery_dir, "/ADWGS_element_3n_recovery_file.rsav")
	element_3n = NULL
	load_result = suppressWarnings(try(load(element_3n_filename), silent = TRUE))
	if (class(load_result) == try_error) {
		element_3n = add_trinucs_by_elements(get_3n(element_coords))
		if (!is.na(recovery_dir)) {
			save(element_3n, file = element_3n_filename)
		}
	}

	window_coords_filename = paste0(recovery_dir, "/ADWGS_window_coords_recovery_file.rsav")
	window_coords = NULL
	load_result = suppressWarnings(try(load(window_coords_filename), silent = TRUE))
	if (class(load_result) == try_error) {
		window_coords = coords_sanity(get_window_coords(element_coords, window_size))
		if (!is.na(recovery_dir)) {
			save(window_coords, file = window_coords_filename)
		}
	}

	window_3n_filename = paste0(recovery_dir, "/ADWGS_window_3n_recovery_file.rsav")
	window_3n = NULL
	load_result = suppressWarnings(try(load(window_3n_filename), silent = TRUE))
	if (class(load_result) == try_error) {
		window_3n = add_trinucs_by_elements(get_3n(window_coords))
		if (!is.na(recovery_dir)) {
			save(window_3n, file = window_3n_filename)
		}
	}
	
	raw_site_coords_filename = paste0(recovery_dir, "/ADWGS_raw_site_coords_recovery_file.rsav")
	site_coords_filename = paste0(recovery_dir, "/ADWGS_site_coords_recovery_file.rsav")
	raw_site_coords = NULL
	site_coords = NULL
	sites_load_result = suppressWarnings(try(load(site_coords_filename), silent = TRUE))
	if (class(sites_load_result) == try_error) {
		load_result = suppressWarnings(try(load(raw_site_coords_filename), silent = TRUE))
		if (class(load_result) == try_error) {
			raw_site_coords = utils::read.delim(sites_file, stringsAsFactors = FALSE, header = FALSE)
			colnames(raw_site_coords) = c("chr", "starts", "ends", "site_info")
			raw_site_coords$chr = gsub("chr", "", raw_site_coords$chr, ignore.case = TRUE)
			raw_site_coords$chr = paste0("chr", raw_site_coords$chr)
			if (!is.na(recovery_dir)) {
				save(raw_site_coords, file = raw_site_coords_filename)
			}
			cat("Successfully read active sites\n")
		} else {
			cat("Recovered active sites\n")
		}

		if (!is.null(raw_site_coords[[1]])) {
			# first only keep sites that overlap with element coordinates
			site_gr = GenomicRanges::GRanges(raw_site_coords$chr,
				IRanges::IRanges(start=raw_site_coords$starts, end=raw_site_coords$ends))
			element_gr = GenomicRanges::GRanges(element_coords$chr, 
				IRanges::IRanges(start=element_coords$starts, end=element_coords$ends))

			rm(raw_site_coords)

			overlaps_found = GenomicRanges::findOverlaps(element_gr, site_gr)
			site_here_gr = site_gr[S4Vectors::subjectHits(overlaps_found)]

			rm(site_gr, element_gr, overlaps_found)
	
			site_coords = do.call(rbind, parallel::mclapply(1:nrow(element_coords), 
				get_sites_per_element, element_coords, site_here_gr, mc.cores=mc.cores))
			rm(site_here_gr)
			site_coords = coords_sanity(site_coords)
		}
		if (!is.na(recovery_dir)) {
			save(site_coords, file = site_coords_filename)
		}
	} else {
		cat("Recovered active sites\n")
	}
	rm(sites_load_result, element_coords)

	site_3n_filename = paste0(recovery_dir, "/ADWGS_site_3n_recovery_file.rsav")
	site_3n = NULL
	load_result = suppressWarnings(try(load(site_3n_filename), silent = TRUE))
	if (class(load_result) == try_error) {
		site_3n = get_3n(site_coords)
		if (!is.na(recovery_dir)) {
			save(site_3n, file = site_3n_filename)
		}
	}
	
	list(element_coords=element_coords_1, element_3n=element_3n, window_coords=window_coords, window_3n=window_3n, site_coords=site_coords, site_3n=site_3n)

}


#### Find Drivers

get_3n_context_of_mutations = function(mutations) {
	
	legal_dna = c("A", "C", "G", "T")
	
	mutations_snv = mutations[mutations$ref %in% legal_dna & mutations$alt %in% legal_dna,]
	mutations_mnv = mutations[!(mutations$ref %in% legal_dna & mutations$alt %in% legal_dna),]
	
	# snvs can have flanks
	flank_ranges = GenomicRanges::GRanges(mutations_snv$chr, 
			IRanges::IRanges(start=mutations_snv$pos1-1, end=mutations_snv$pos2+1), strand="*")
	triples = as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, flank_ranges))
	
	# complement trinucleotide where necessary to force into one signature space
	new_triples = triples
	new_alt = mutations_snv$alt
	which_to_complement = which(mutations_snv$ref %in% c("G", "A"))

	new_triples[which_to_complement] = as.character(Biostrings::complement(Biostrings::DNAStringSet(new_triples[which_to_complement])))
	new_alt[which_to_complement] = as.character(Biostrings::complement(Biostrings::DNAStringSet(new_alt[which_to_complement])))
	mutations_snv$tag = paste0(new_triples, ">", new_alt)
	
	if (nrow(mutations_mnv)==0) {
		return(mutations_snv)
	}
	
	mutations_mnv$tag = "indel>X"
	rbind(mutations_snv, mutations_mnv)
}

format_muts = function(muts_filename, filter_hyper_MB) {

	muts = utils::read.delim(muts_filename, stringsAsFactors = FALSE, header = FALSE)
	
	if (ncol(muts) != 6) {
		stop("Error: mutation file is not in the correct format")
	}

	colnames(muts) = c("chr", "pos1", "pos2", "ref", "alt", "patient")
	
	# remove hypermutated samples, according to muts/megabase rate defined
	if (!is.na(filter_hyper_MB) & filter_hyper_MB>0) {
		total_muts_filter = 3000*filter_hyper_MB
		sample_mut_count = table(muts$patient)
		hyper_tab = sample_mut_count[sample_mut_count>total_muts_filter]
		spl_rm = names(hyper_tab)
		no_mut_rm = sum(hyper_tab)
		muts = muts[!muts$patient %in% spl_rm,, drop=F]
	}

	# keep only relevant chrs, make sure CHR is present in address
	muts$chr = gsub("chr", "", muts$chr, ignore.case=TRUE)
	muts = muts[muts$chr %in% c(1:22, "Y", "X", "M"),]	
	muts$chr = paste0("chr", muts$chr)
	
	muts$pos1 = as.numeric(muts$pos1)
	muts$pos2 = as.numeric(muts$pos2)
	
	# reverse start/end coordinates of deletions
	rev_coords = which(muts$pos2-muts$pos1<0)
	if (length(rev_coords)>0) {
		pos1 = muts[rev_coords, "pos1"]
		pos2 = muts[rev_coords, "pos2"]
		muts[rev_coords, "pos1"] = pos2
		muts[rev_coords, "pos2"] = pos1
		rm(pos1, pos2)
	}
	muts = get_3n_context_of_mutations(muts)
	muts
}

merge_granges_overlaps = function(muts_gr, muts_tab, regs_gr, regs_tab, compute_overlap=F) {
	overlaps_found = GenomicRanges::findOverlaps(muts_gr, regs_gr)
	muts_in_regs = muts_tab[S4Vectors::queryHits(overlaps_found),]
	regs_with_muts = regs_tab[S4Vectors::subjectHits(overlaps_found),]
	
	colnames(muts_in_regs) = paste("mut", colnames(muts_in_regs), sep="_")
	colnames(regs_with_muts) = paste("reg", colnames(regs_with_muts), sep="_")
	merged_ranges = cbind(muts_in_regs, regs_with_muts)
	
	if (compute_overlap) {
		gr_intersects = GenomicRanges::pintersect(muts_gr[S4Vectors::queryHits(overlaps_found)], regs_gr[S4Vectors::subjectHits(overlaps_found)])
		gr_unions = GenomicRanges::punion(muts_gr[S4Vectors::queryHits(overlaps_found)], regs_gr[S4Vectors::subjectHits(overlaps_found)])
		percentIntersect = BiocGenerics::width(gr_intersects) / BiocGenerics::width(muts_gr[S4Vectors::queryHits(overlaps_found)])
		jaccard = BiocGenerics::width(gr_intersects)/BiocGenerics::width(gr_unions)
		merged_ranges = cbind(merged_ranges, percentIntersect, jaccard)
	}
	merged_ranges
}

merge_elements_snvs = function(coords, snvs) {
	if (is.null(coords[[1]])) {
		return(NULL)
	}
	e_gr = GenomicRanges::GRanges(coords$chr, IRanges::IRanges(start=coords$starts, end=coords$ends), mcols=coords$element_id)
	m_gr = GenomicRanges::GRanges(snvs$chr, IRanges::IRanges(start=snvs$pos1, end=snvs$pos2), mcols=snvs$tag)
	m_e = merge_granges_overlaps(m_gr, snvs, e_gr, coords)
	rownames(m_e) = NULL
	m_e
}

get_todo_for_elements = function(element_coords, mutations_in_elements) {
	elements_with_muts = intersect(element_coords$id, mutations_in_elements$reg_id)
}

get_unmutated_elements = function(element_coords, mutations_in_elements) {
	elements_without_muts = setdiff(element_coords$id, mutations_in_elements$reg_id)
	
	if (length(elements_without_muts) == 0) {
		return(NULL)
	}

	blank_result = data.frame(id=elements_without_muts, pp_site=NA, pp_element=NA, element_muts_obs=NA, 
					element_muts_exp=NA, element_enriched=NA, site_muts_obs=NA, site_muts_exp=NA, site_enriched=NA,
					selected_sites=NA, h0_df=NA, h1_df=NA, stringsAsFactors=F)
	as.matrix(blank_result)
}

fix_all_results = function(all_results) {
	
	resi = data.frame(all_results, stringsAsFactors=F)
	resi$pp_site = as.numeric(resi$pp_site)
	resi$pp_element = as.numeric(resi$pp_element)
	resi$element_muts_obs = as.numeric(resi$element_muts_obs)
	resi$element_muts_exp = as.numeric(resi$element_muts_exp)
	resi$element_enriched = as.logical(gsub("\\s+", "", resi$element_enriched))
	resi$site_muts_obs = as.numeric(resi$site_muts_obs)
	resi$site_muts_exp = as.numeric(resi$site_muts_exp)
	resi$site_enriched = as.logical(gsub("\\s+", "", resi$site_enriched))

	# DO NOT set NA's to ones here - deal with downstream with significance filtering	
	resi[!is.na(resi$pp_element) & resi$pp_element==0,"pp_element"] = 1e-300
	resi[!is.na(resi$pp_site) & resi$pp_site==0,"pp_site"] = 1e-300
	
	## remove duplicated entries that come from computing the same element/mut in multiple rounds
	resi_tag = resi[,"id"]
	resi = resi[!duplicated(resi_tag),]
	resi
}

get_signf_results = function(all_res) {
	this_results = all_res
	if (nrow(this_results)==0) {
		return(NULL)
	}
	# this is FDR treating element-level NAs as 1s
	this_results$fdr_element = stats::p.adjust(this_results$pp_element, method="fdr", n=nrow(this_results))
	
	this_results = this_results[order(this_results$fdr_element),]
	
	filtered_results = this_results[!is.na(this_results$fdr_element) & this_results$fdr_element<0.05,]
	unsignf_results = this_results[is.na(this_results$fdr_element) | this_results$fdr_element>=0.05,]

	if (nrow(filtered_results) + nrow(unsignf_results) != nrow(this_results)) {
		stop("Error: Something unexpected happened when formatting results")
	}

	if (nrow(filtered_results)!=0) {
		# site-level FDR perform only on elements with pre-selection of FDR<0.05
		filtered_results$fdr_site = stats::p.adjust(filtered_results$pp_site, method="fdr", n=nrow(filtered_results))
		filtered_results$has_site_mutations = !is.na(filtered_results$fdr_site) & filtered_results$fdr_site<0.05
		filtered_results$has_site_mutations = c("","V")[1+c(filtered_results$has_site_mutations)]
	}
	
	if (nrow(unsignf_results) !=0) {
		unsignf_results$fdr_site = NA
		unsignf_results$has_site_mutations = NA
	}
	
	
	final_results = rbind(filtered_results, unsignf_results)
	final_results$pp_element = replace(final_results$pp_element, which(is.na(final_results$pp_element)), 1)
	final_results$pp_site = replace(final_results$pp_site, which(is.na(final_results$pp_site)), 1)
	final_results$fdr_element = replace(final_results$fdr_element, which(is.na(final_results$fdr_element)), 1)
	final_results$fdr_site = replace(final_results$fdr_site, which(is.na(final_results$fdr_site)), 1)
	final_results
}

#' ActiveDriverWGS analysis of cancer mutation information
#'
#' @param prepared_elements the elements and active sites to analyze, as prepared by \code{\link{prepare_elements}}
#' @param mutations_file path to a tab deliminated text file with 6 columns
#'    (chromosome, start position, end position, reference allele, alternate allele, patient identifier)
#' @param filter_hyper_MB numeric specifying the mutation frequency to filter. Default value is 30
#' @param mc.cores numeric indicating the number of computing cores to use. Default value is 1.
#' @param recovery_dir path to a directory for ActiveDriverWGS to store intermediate computations in.
#'             In the case that find_drivers does not complete, run it again with the same value for this parameter,
#'             and it will read the intermediate computations rather than compute them again. Please specify a
#'             different directory, or delete the ActiveDriverWGS data in previously used directory, if running with
#'             different data. Default value is NA, and find_drivers will not save any intermediate information.
#'
#' @return data.frame containing all results from the ActiveDriverWGS analysis
#'
#' @export
find_drivers = function(prepared_elements, mutations_file, filter_hyper_MB = 30, mc.cores = 1, recovery_dir = NA) {
	try_error = "try-error"

	element_coords = prepared_elements$element_coords
	element_3n = prepared_elements$element_3n
	window_coords = prepared_elements$window_coords
	window_3n = prepared_elements$window_3n
	site_coords = prepared_elements$site_coords
	site_3n = prepared_elements$site_3n
	rm(prepared_elements)

	muts_filename = paste0(recovery_dir, "/ADWGS_mutations_recovery_file.rsav")
	mutations_in_windows_filename = paste0(recovery_dir, "/ADWGS_mutations_in_windows_recovery_file.rsav")
	muts = NULL
	mutations_in_windows = NULL
	in_windows_load_result = suppressWarnings(try(load(mutations_in_windows_filename), silent = TRUE))
	if (class(in_windows_load_result) == try_error) {
		load_result = suppressWarnings(try(load(muts_filename), silent = TRUE))
		if (class(load_result) == try_error) {
			muts = format_muts(mutations_file, filter_hyper_MB)
			if (!is.na(recovery_dir)) {
				save(muts, file = muts_filename)
			}
			cat("Successfully read mutations\n")
		} else {
			cat("Recovered mutations\n")
		}
	}
	rm(muts_filename, mutations_file)

	mutations_in_sites = NULL
	mutations_in_sites_filename = paste0(recovery_dir, "/ADWGS_mutations_in_sites_recovery_file.rsav")
	load_result = suppressWarnings(try(load(mutations_in_sites_filename), silent = TRUE))
	if (class(load_result) == try_error) {
		mutations_in_sites = merge_elements_snvs(site_coords, muts)
		if (!is.na(recovery_dir)) {
			save(mutations_in_sites, file = mutations_in_sites_filename)
		}
	}
	rm(mutations_in_sites_filename)

	mutations_in_elements = NULL
	mutations_in_elements_filename = paste0(recovery_dir, "/ADWGS_mutations_in_elements_recovery_file.rsav")
	load_result = suppressWarnings(try(load(mutations_in_elements_filename), silent = TRUE))
	if(class(load_result) == try_error) {
		mutations_in_elements = merge_elements_snvs(element_coords, muts)
		if (!is.na(recovery_dir)) {
			save(mutations_in_elements, file = mutations_in_elements_filename)
		}
	}
	rm(mutations_in_elements_filename)

	if (class(in_windows_load_result) == try_error) {
		mutations_in_windows = merge_elements_snvs(window_coords, muts)
		if (!is.na(recovery_dir)) {
			save(mutations_in_windows, file = mutations_in_windows_filename)
		}
		cat("Successfully prepared mutations\n")
	} else {
		cat("Recovered prepared mutations\n")
	}
	rm(mutations_in_windows_filename, muts, in_windows_load_result)

	not_done = NULL
	not_done_filename = paste0(recovery_dir, "/ADWGS_TODO_recovery_file.rsav")
	load_result = suppressWarnings(try(load(not_done_filename), silent = TRUE))
	if (class(load_result) == try_error) {
		not_done = get_todo_for_elements(element_coords, mutations_in_elements)
		if (!is.na(recovery_dir)) {
			save(not_done, file = not_done_filename)
		}
	}

	recovered_results = NULL
	recovered_result_numbers = c()
	if (!is.na(recovery_dir)) {
		results_filenames = list.files(recovery_dir, pattern = "ADWGS_result[0123456789]+_recovery_file.rsav")
		if (length(results_filenames) > 0) results_filenames = paste0("/", results_filenames)
		recovered_results = do.call(rbind, lapply(results_filenames, function(filename) {
				load_result = suppressWarnings(try(load(paste0(recovery_dir, filename)), silent = TRUE))
				if (class(load_result) == try_error) return(NULL)
				if (ncol(result) != 13) return(NULL)
				result
			}))
		recovered_result_numbers = recovered_results$result_number
	}

	cat("Tests to do: ", length(not_done), "\n")
	if (length(recovered_result_numbers) > 0) cat("Tests recovered: ", length(recovered_result_numbers), "\n")

	mutated_results = do.call(rbind, parallel::mclapply(1:length(not_done), function(i) {
		if (i %% 100 == 0) cat(i, "\n")
		if (i %in% recovered_result_numbers) return(NULL)
		element_id = not_done[i]
		result = regress_test(element_id, mutations_in_sites, mutations_in_elements,
			mutations_in_windows, site_3n, element_3n, window_3n)
		if (!is.na(recovery_dir)) {
			result$result_number = i
			save(result, file = paste0(recovery_dir, "/ADWGS_result", i, "_recovery_file.rsav"))
		}
		result
	}, mc.cores = mc.cores))
	rm(element_3n, window_coords, window_3n, site_coords, site_3n, mutations_in_sites, mutations_in_windows, not_done)

	mutated_results = rbind(recovered_results, mutated_results)
	rm(recovered_results, recovered_result_numbers)

	if (!is.na(recovery_dir)) {
		mutated_results = mutated_results[, 1:12]
	}

	mutated_results = as.matrix(mutated_results)
	unmutated_results =  get_unmutated_elements(element_coords, mutations_in_elements)
	rm(element_coords, mutations_in_elements)

	all_results = rbind(mutated_results, unmutated_results)
	rm(mutated_results, unmutated_results)

	all_results = fix_all_results(all_results)
	all_results = get_signf_results(all_results)
	all_results
}
