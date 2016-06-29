suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(optparse))

help_text = "\nThis R script builds a coverage table from bed files that described hicoverage interval within bmr references (gene intervals)."  
option_list <- list(
	make_option(c("-T", "--threads"), type="integer", default=1,
							help="Number of threads for gene BMR reference creation. BMR reference creation is time intensive. Running in parallel is suggested."),
	
	make_option(c("-G","--geneSubset"),type="character", default=NULL,
							help = "Required parameter used to provide a file specifying a list of genes for which to build the coverage table."),

	make_option(c("-B","--bmrBeds"),type="character", default=NULL,
							help = "Directory of beds overlapping gene intervals as defined by BMR references."),
	
	make_option(c("-S","--sampleSubset"),type="character", default=NULL,
							help = "Optional parameter used to provide a file specifying subset of samples for which to create coverage table. Otherwise will use all samples provided in bmrBeds."),
	
	make_option(c("-R","--bmrRefs"),type="character", default=NULL,
							help = "Directory of BMR references used to create coverage table."),
	make_option(c("-C","--combineAmbigs"),type="logical", action="store_true",default=F,
							help = "Option to combine coverage data for nonoverlapping gene intervals. Default action is to process and report data for \"ambigs\" separately.")
	
)


args <- parse_args(OptionParser(description=help_text, option=option_list))
if(args$threads > 1){
	registerDoParallel(cores=args$threads)
	print(paste("Running coverage table creation in parallel with ", args$threads, " threads.", sep=""))	
}

if(is.null(args$geneSubset)){
	print("File containing list of genes to create coverage table was not provided. Aborting")
	quit()
}else{
	print(paste("Creating coverage table using the genes specified in ", args$geneSubset,".", sep=""))
	genes_considered = read.delim(args$geneSubset, header=F, stringsAsFactors = F)[,1]
}

if(is.null(args$bmrRefs)){
	print("Directory containing BMR references was not provided. Aborting")
	quit()
}

if(is.null(args$bmrBeds)){
# 	print("Directory containing bed files used to create coverage table was not provided. Aborting")
# 	quit()
	
	print("No BMR bed file provided. Aborting")
	quit()
}else{
	print(paste("Creating coverage table using the bmr beds specified in ", args$bmrBeds,".", sep=""))
# 	all_samps = list.files(path = args$bmrBeds, pattern = "*.hg19_covered.bps.bed")
#   
# 	if(!is.null(args$sampleSubset)){
# 		print(paste("Creating coverage table using samples specified in", args$geneSubset,".", sep=""))
# 		samps = read.delim(args$sampleSubset, header=F, stringsAsFactors = F)[,1]
# 		samps = samps[which(paste(samps,".hg19_covered.bps.bed",sep="") %in% all_samps)]
# 	}else{samps = substr(all_samps, 1, nchar(all_samps)-21)}
# 	
# 	if(length(samps) == 0){
# 		print('No bmr bed files ending with ".hg19_covered.bps.bed" were provided. Aborting')
# 		quit()
# 	}
# 	samp_covs = list()
# 	samp_index = 1
# 	for(samp in samps){
# 		print(paste("Reading bmr bed file: ", samp_index, " of ", length(samps),sep=''))
# 		samp_index = samp_index + 1
# 		samp_covs[[length(samp_covs)+1]] = read.delim(paste(args$bmrBeds,"/",samp,".hg19_covered.bps.bed",sep=""), header=F, stringsAsFactors = F)
# 	}
	
# 	covs <- foreach (samp_index=1:length(samps), .combine = rbind) %dopar% {
# 	
# 		samp = samps[samp_index]
# 		print(paste("Reading bmr bed file: ", samp_index, " of ", length(samps),sep=''))
# 		samp_index = samp_index + 1
# 		intervals = read.delim(paste(args$bmrBeds,"/",samp,".hg19_covered.bps.bed",sep=""), header=F, stringsAsFactors = F)
# 		cbind(rep(samp, nrow(intervals)), intervals, intervals[,4])
# 	}
	write("Started reading.", file="pull_coverage.log", append=T)

	covs = read.delim(args$bmrBeds,header=F, stringsAsFactors = F)
	samps = unique(covs[,5])
	covs = cbind(covs, covs[,4])
	print("Finished reading bmr beds. Processing by gene.")
	write("Finished reading bmr beds. Processing by gene.", file="pull_coverage.log", append=T)

	covs = split(data.frame(covs), covs[,6])
	covs = covs[genes_considered[which(genes_considered %in% names(covs))]]


# 	samp_covs <- foreach (samp_index=1:length(samps)) %dopar% {
# 		
# 		samp = samps[samp_index]
# 		print(paste("Reading bmr bed file: ", samp_index, " of ", length(samps),sep=''))
# 		samp_index = samp_index + 1
# 		intervals = read.delim(paste(args$bmrBeds,"/",samp,".hg19_covered.bps.bed",sep=""), header=F, stringsAsFactors = F)
# 		intervals
# 	}
# 	names(samp_covs) =  samps
	
# 	samp_covs = list()
# 	for(samp in samps){
# 		samp_covs[[length(samp_covs)+1]] = covs[which(covs[,1] == samp),2:5]
# 	}
# 	samp_covs = split(data.frame(covs), covs[,1])
# 	covs = NULL
	
	#print("Finished reading bed files")
}

#### Order of rows from BMRs
## silent:1-7
## nonsilent:1-7
## noncoding:1-7
# 

categs = c(paste(rep("silent",7), 1:7,sep=":"), paste(rep("nonsilent",7), 1:7,sep=":"), paste(rep("noncoding",7), 1:7,sep=":"))
if(args$combineAmbigs){
  print("Combining coverage data for nonoverlapping gene intervals")	
}else{
	print("Treating nonoverlapping gene intervals (amigs) separately. Default action.")
}

print("Generating coverage table.")
write("Generating coverage table.", file="pull_coverage.log", append=T)

#print(names(samp_covs))
#for(gene in genes_considered){

gene_names = names(covs)
not_found = genes_considered[which(!(genes_considered %in% gene_names))]
print(paste("# of genes not found in BMR beds: ", length(not_found),sep=""))
print("Printing missing genes to missing.genes.txt")
write.table(not_found, file="missing.genes.txt", quote=F, col.names=F, row.names=F)

#for(gene_cov in covs){
results <- foreach (gene_cov=covs, .combine = rbind) %dopar% {
	#index = which(genes_considered == gene)
	gene = unique(gene_cov[,4])
	gene_cov = gene_cov[,-4]

	if(is.na(gene)){print(gene_cov)}
	gene_index = which(gene == gene_names)
	if(gene_index %% 3 == 0){print(paste("Starting gene: ", gene_index, " of ", length(gene_names)," at time ", Sys.time(), ".",sep=""))
		write(paste("Starting gene: ", gene_index, " of ", length(gene_names)," at time ", Sys.time(), ".",sep=""), file="pull_coverage.log", append=T)}
	#print(gene)
	#write(gene, file="pull_coverage.log", append=T)
	samp_covs = split(data.frame(gene_cov), gene_cov[,4])
	gene_bmrFile <- tryCatch({readLines(paste(args$bmrRefs,"/",gene,".bmr",sep=""))},
													 error=function(cond){
													 	message(paste("BMR not present: ",gene, sep=""))
													 	return(NULL)
													 	}, 
													 warning=function(cond){
													 	message(paste("BMR not present: ",gene, sep=""))
													 	return(NULL)
													 	})
	if(is.na(gene_bmrFile[1])){ return(NULL)}

	cover_per_tumor = matrix(0, nrow=21, ncol=length(samps))
	rownames(cover_per_tumor) = paste(rep(gene,21), categs, sep=":")
	colnames(cover_per_tumor) = samps
	
	gene_info = strsplit(gene_bmrFile[1],":")
	gene_name = gene_info[[1]][1]
	gene_chrom = gene_info[[1]][2]
	gene_interval = as.numeric(strsplit(gene_info[[1]][3],"-")[[1]])
	gene_bmr = matrix(as.numeric(unlist(lapply(gene_bmrFile[-1], function(.line) strsplit(.line, '')[[1]]) )), nrow=21, byrow =T)
	gene_bmr = gene_bmr / 3

	for(tmp in 1:length(samp_covs)){
		samp = names(samp_covs)[tmp]
		sub_cov = samp_covs[[tmp]]
		if(is.null(sub_cov)){
			next
		}
		starts = as.numeric(sub_cov[,2])-gene_interval[1]+1
		ends = as.numeric(sub_cov[,3])-gene_interval[1]
		cover_per_tumor[,samp] = rowSums(gene_bmr[,unlist(mapply(seq, starts,ends)), drop=F])
	}
   cover_per_tumor
}
colnames(results) = samps
cover_out = round(results)
info = matrix(unlist(strsplit(rownames(cover_out), ":")), nrow=length(rownames(cover_out)), byrow=T)
old_names = colnames(cover_out)
cover_out = cbind(info, cover_out)
colnames(cover_out) = c("gene", "effect", "categ", old_names)
write.table(cover_out, file="MutSig_coverage_data_table_mTCGA.txt", col.names=T, row.names=F, quote=F, sep="\t")
print("Finished generation.")


												 	
												 	