suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))

args = commandArgs(trailingOnly=TRUE)


print(args[1])
cl<-makeCluster(as.numeric(args[1]),outfile="")
registerDoParallel(cl)

#### Order of rows from BMRs
## silent:1-7
## nonsilent:1-7
## noncoding:1-7
# 
# gene_bmrFile = readLines("bmrs/finished_bmrs/Adamts6.bmr.csv")
# gene_info = strsplit(gene_bmrFile[1],":")
# gene_name = gene_info[[1]][1]
# gene_chrom = gene_info[[1]][2]
# gene_interval = as.numeric(strsplit(gene_info[[1]][3],"-")[[1]])
# gene_bmr = matrix(as.numeric(unlist(lapply(gene_bmrFile[-1], function(.line) strsplit(.line, '')[[1]]) )), nrow=21, byrow =T)
# gene_bmr = gene_bmr / 3
# 
# sub_cov = covered[which(covered[,4] == gene_name),]
categs = c(paste(rep("silent",7), 1:7,sep=":"), paste(rep("nonsilent",7), 1:7,sep=":"), paste(rep("noncoding",7), 1:7,sep=":"))

samps = read.delim("samps.txt", header=F, stringsAsFactors = F)[,1]
samp_covs = list()
samp_header = samps
for(samp in samps){
	samp_covs[[length(samp_covs)+1]] = read.delim(paste("finished_covs/",samp,".covered_bps.bed",sep=""), header=F, stringsAsFactors = F)
}


genes_considered = read.delim("genes_considered.txt", header=F, stringsAsFactors = F)[,1]
cover_per_tumor = matrix(NA, nrow=length(genes_considered)*21, ncol=length(samps))
cover_per_tumor_double = matrix(NA, nrow=length(genes_considered)*21, ncol=length(samps))

rownames(cover_per_tumor) = paste(rep(genes_considered,each=21), categs, sep=":")
rownames(cover_per_tumor_double) = paste(rep(genes_considered,each=21), categs,sep=":")
colnames(cover_per_tumor) = samps
colnames(cover_per_tumor_double) = samps

print(Sys.time())


#for(gene in genes_considered[1:20]){
results <- foreach (index=1:length(genes_considered), .combine = rbind) %dopar% {
	#index = which(genes_considered == gene)
	gene = genes_considered[index]

	write(gene, file="pull_coverage.log", append=T)
	gene_bmrFile = readLines(paste("finished_bmrs/",gene,".bmr.csv",sep=""))
	gene_info = strsplit(gene_bmrFile[1],":")
	gene_name = gene_info[[1]][1]
	gene_chrom = gene_info[[1]][2]
	gene_interval = as.numeric(strsplit(gene_info[[1]][3],"-")[[1]])
	gene_bmr = matrix(as.numeric(unlist(lapply(gene_bmrFile[-1], function(.line) strsplit(.line, '')[[1]]) )), nrow=21, byrow =T)
	gene_bmr = gene_bmr / 3
	start_index = ((index -1) * 21 + 1)
	stop_index = start_index+20
	
	for(tmp in 1:length(samp_covs)){
		
		covered = samp_covs[[tmp]]
		sub_cov = covered[which(covered[,4] == gene_name),]
		if(nrow(sub_cov) == 0){
			cover_per_tumor[start_index:stop_index,tmp] = rep(0,7)
			cover_per_tumor_double[start_index:stop_index,tmp] = rep(0,7)
			next
		}
	
		samp_chrom = unique(sub_cov[,1])[1]
		if(samp_chrom == gene_chrom){
			
			starts = as.numeric(sub_cov[,2])-gene_interval[1]+1
			ends = as.numeric(sub_cov[,3])-gene_interval[1]
			
		}else{
			cover_per_tumor[start_index:stop_index,tmp] = rep(0,7)
		  cover_per_tumor_double[start_index:stop_index,tmp] = rep(0,7)
		next}
		
	
		cover_per_tumor[start_index:stop_index,tmp] = rowSums(gene_bmr[,unlist(mapply(seq, starts,ends))])
		cover_per_tumor_double[start_index:stop_index,tmp] = rowSums(ceiling(gene_bmr[,unlist(mapply(seq, starts,ends))]))
	}

	cbind(cover_per_tumor[start_index:stop_index,], cover_per_tumor_double[start_index:stop_index,])
}
stopCluster(cl)
cover_out = results[,1:length(samps)]
cover_out_double = results[,(length(samps)+1):ncol(results)]

write.table(round(cover_out), file="MutSig_coverage_data_table_mTCGA.txt", col.names=T, row.names=T, quote=F, sep="\t")
write.table(cover_out_double, file="MutSig_coverage_data_table_mTCGA_double_count.txt", col.names=T, row.names=T, quote=F, sep="\t")




												 	
												 	