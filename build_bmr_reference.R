suppressPackageStartupMessages(library(optparse))

help_text = "\nThis R script builds a BMR reference to be used for MutSigCV coverage table generation."  
option_list <- list(
	make_option(c("-T", "--threads"), type="integer", default=1,
							help="Number of threads for gene BMR reference creation. BMR reference creation is time instensive. Running in parallel is suggested."),
	           
	make_option(c("-S","--geneSubset"),type="character", default=NULL,
							help = "Optimal parameter to provide a file that contains a list of genes for which to build BMR reference files."),
	
	make_option(c("-F","--fasta"),type="character", default=NULL,
							help = "The reference genome in fasta format to be used for BMR reference creation"),
	
	make_option(c("-R","--refGene"),type="character", default=NULL,
							help = "The transcript annotation file to be used in the majority of cases during BMR reference creation in refFlat format (gene names as IDs). refGene is suggested."),
	
	make_option(c("-E","--ensl"),type="character", default=NULL,
							help = "The transcript annotation file to be used in the minor cases where gene annotation is nonsensical (for example, the number of coding basepairs is not divisible by three or translation is nonsense). This happens in genes with micro-introns of bps<3. Ensbl is suggested."))

#Since mutation annotations are prioritized s by consequence (as per snpEff), 
#I am going to the same with the coverage data (that is if the bp could be coding or noncoding, it will be specified as coding)
args <- parse_args(OptionParser(description=help_text, option=option_list))
#args = commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library("bedr"))
suppressPackageStartupMessages(library("Biostrings"))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))

if(is.null(args$fasta)){
	print("Reference fasta file not given. Aborting")
	quit()
}else{
#	fastaFile = "~/Raw/NGS/genome.builds/mm10/mm10.norandom.fa"
	fastaFile = args$fasta
}

if(is.null(args$refGene)){
	print("Reference major gene annotation file not given. Aborting")
	quit()
}else{
	
	refGene = read.delim(args$refGene, header=T, stringsAsFactors = F)
	print(paste("Using ", args$refGene, " as major gene annotation reference.", sep=""))
	colnames(refGene)  = c("Gene", "Isoform", "Chrom", "Strand", "TxStart", "TxEnd","CodingStart", 
												 "CodingEnd", "NumExons", "ExonStarts", "ExonEnds")
	
	refGene$Chrom = substr(refGene$Chrom, 4, nchar(refGene$Chrom))
	##Remove non-canonical chromosomes
	refGene = refGene[which(refGene$Chrom %in% c(1:22, "X", "Y", "M")),]
	##Only keep transcript that specify coding regions
}

if(is.null(args$ensl)){
	print("Reference minor gene annotation file not given. Aborting")
	quit()
}else{
	##Used "mm10_genes_Ensembl2.bed" in first generation
  ensl = read.delim(args$ensl, header=T, stringsAsFactors = F)	
  print(paste("Using ", args$ensl, " as minor gene annotation reference.", sep=""))
  #colnames(ensl)[4] = "Strand"
  ensl$Chrom = substr(ensl$Chrom, 4, nchar(ensl$Chrom))
  ##Remove non-canonical chromosomes
  ensl = ensl[which(ensl$Chrom %in% c(1:22, "X", "Y", "M")),]
  colnames(ensl) = c("Gene", "Isoform", "Chrom", "Strand", "TxStart", "TxEnd","CodingStart", 
                        "CodingEnd", "NumExons", "ExonStarts", "ExonEnds")
}

if(args$threads > 1){
	registerDoParallel(cores=args$threads)
	print(paste("Running BMR reference creation in parallel with ", args$threads, " threads.", sep=""))	
}

if(!is.null(args$geneSubset)){
	print(paste("Running BMR reference creation on the subset of genes provided in ", args$geneSubset, sep=""))	
	genes_to_run = read.delim(args$geneSubset,header=F,stringsAsFactors = F)[,1]
}else{
	genes_to_run = unique(refGene[,2])
	print("Creating BMR reference for each gene in major gene annotation")
}

#####annotations
options(scipen=999)

nts = c("A", "T","C","G")
codon_table = rbind(c("ATT","I"), c("ATC","I"),c("ATA","I"),c("CTT","L"),c("CTC","L"),c("CTA","L"),c("CTG","L"),
										c("TTA","L"),c("TTG","L"), c("GTT","V"),c("GTC","V"),c("GTA","V"),c("GTG","V"),c("TTT","F"),
										c("TTC","F"),c("ATG","M"),c("TGT","C"), c("TGC","C"),c("GCT","A"),c("GCC","A"),c("GCA","A"),
										c("GCG","A"),c("GGT","G"),c("GGC","G"),c("GGA","G"),c("GGG","G"),c("CCT","P"),c("CCC","P"),
										c("CCA","P"),c("CCG","P"),c("ACT","T"),c("ACC","T"),c("ACA","T"),c("ACG","T"),c("TCT","S"),
										c("TCC","S"),c("TCA","S"),c("TCG","S"),c("AGT","S"),c("AGC","S"),c("TAT","Y"),c("TAC","Y"),
										c("TGG","W"),c("CAA","Q"),c("CAG","Q"),c("AAT","N"),c("AAC","N"),c("CAT","H"),c("CAC","H"),
										c("GAA","E"),c("GAG","E"),c("GAT","D"),c("GAC","D"),c("AAA","K"),c("AAG","K"),c("CGT","R"),
										c("CGC","R"),c("CGA","R"),c("CGG","R"),c("AGA","R"),c("AGG","R"),c("TAA","Stop"),c("TAG","Stop"),
										c("TGA","Stop"))

rownames(codon_table) = codon_table[,1]
transversions = c("A:T", "T:A", "A:C","C:A", "T:G", "G:T", "G:C", "C:G")
transitions = c("A:G", "G:A", "C:T", "T:C")
categs = 1:7
names(categs) = c("CpGtransition", "CpGtransversion", "C:Gtransition","C:Gtransversion", "A:Ttransition", "A:Ttransversion","indels")
categ_rows = c(paste(rep("silent",7),c(1:7), sep=""), paste(rep("nonsilent",7),c(1:7), sep=""), paste(rep("noncoding",7),c(1:7), sep=""))

logFile = "bmr.reference.log"
failedLog = "bmr.failed.log"
protein_fail = "protein_noncanonical.log"

createBMR <- function(transcriptSubset, ambiguousGeneNumber=NULL){

	subGenes = transcriptSubset
	gene = subGenes[1,1]
	##Get total interval of gene
	intervalStart = as.numeric(min(subGenes[,"TxStart"]))-5000
	intervalEnd = as.numeric(max(subGenes[,"TxEnd"]))+5000
	intervalLength = intervalEnd - intervalStart 
	chrom = subGenes[1,"Chrom"]
	fastaSeq = capture.output(get.fasta(x=paste("chr",subGenes[1,"Chrom"],":",intervalStart-1,"-",intervalEnd+1,sep=""),fasta=fastaFile))
	fastaSeq = substr(fastaSeq[20], regexpr(" ", fastaSeq[20])+1, nchar(fastaSeq[20]))

	if(subGenes[1,"Strand"] == "-"){
		fastaSeq = toString(complement(DNAString(fastaSeq)))
	}
  
	coding_bools = matrix(0, nrow(subGenes), intervalLength)
	colnames(coding_bools) = (intervalStart+1):intervalEnd
	firstPass = TRUE
	numGenes = nrow(subGenes)
	switchEnsl = rep(F, numGenes)
	for(row in 1:numGenes){
	
		subGene = subGenes[row,]
		###Define coding and noncoding regions
		codingStart = as.numeric(subGene[,"CodingStart"])+1
		starts = as.numeric(strsplit(subGene[,"ExonStarts"],",")[[1]])+1
		ends = as.numeric(strsplit(subGene[,"ExonEnds"],",")[[1]])
		coding_bool = matrix(0, 1, intervalLength)
		names(coding_bool) = (intervalStart+1):intervalEnd
	
		for(i in 1:as.numeric(subGene[,"NumExons"])){
			
			bps = starts[i]:ends[i]
			bps = bps[which(bps > as.numeric(subGene[,"CodingStart"]) & bps <= as.numeric(subGene[,"CodingEnd"]))]
			
			coding_bool[as.character(bps)] = 1
		} 
	 
		coding_bools[row,] = coding_bool

		if(length(which(coding_bool == 1)) %% 3 != 0){
		  switchEnsl[row] = T
		}
	  
	}
	if(all(switchEnsl)){
		correctCoding = rep(F, numGenes)
	  subGenes = ensl[which(ensl$Gene == gene),]
    subGenes = subGenes[which(((subGenes$TxStart >= intervalStart & subGenes$TxStart <= intervalEnd) | (subGenes$TxEnd >= intervalStart & subGenes$TxEnd <= intervalEnd))  & subGenes$Chrom == chrom),]
           
	  if(nrow(subGenes) == 0){
	    write(paste("Nonsense transcripts in main, no transcripts in minor: ", gene, ambiguousGeneNumber, sep=" "), file=failedLog, append=T) 
	    print("yup")
	  	return(NULL)
	  }
	  ##Get total interval of gene
	  intervalStart = as.numeric(min(subGenes[,"TxStart"]))-5000
	  intervalEnd = as.numeric(max(subGenes[,"TxEnd"]))+5000
	  intervalLength = intervalEnd - intervalStart 
	  fastaSeq = capture.output(get.fasta(x=paste("chr",subGenes[1,"Chrom"],":",intervalStart-1,"-",intervalEnd+1,sep=""),fasta=fastaFile))
	  fastaSeq = substr(fastaSeq[20], regexpr(" ", fastaSeq[20])+1, nchar(fastaSeq[20]))
	  if(subGenes[1,"Strand"] == "-"){
	    fastaSeq = toString(complement(DNAString(fastaSeq)))
	  }
	  
	  
	  coding_bools = matrix(0, nrow(subGenes), intervalLength)
	  colnames(coding_bools) = (intervalStart+1):intervalEnd
	  firstPass = TRUE
	  numGenes = nrow(subGenes)
	  for(row in 1:numGenes){
	    subGene = subGenes[row,]
	    ###Define coding and noncoding regions

	    starts = as.numeric(strsplit(subGene[,"ExonStarts"],",")[[1]])+1
	    ends = as.numeric(strsplit(subGene[,"ExonEnds"],",")[[1]])
	    
	    
	    coding_bool = matrix(0, 1, intervalLength)
	    names(coding_bool) = (intervalStart+1):intervalEnd
	    sumbp = 0
	    for(i in 1:as.numeric(subGene[,"NumExons"])){
	      
	      bps = starts[i]:ends[i]
	      bps = bps[which(bps > as.numeric(subGene[,"CodingStart"]) & bps <= as.numeric(subGene[,"CodingEnd"]))]
	      sumbp = sumbp + length(bps)
	      coding_bool[as.character(bps)] = 1
	    } 
	    coding_bools[row,] = coding_bool
	    if(length(which(coding_bool == 1)) %% 3 == 0){
	    	correctCoding[row] = T
	    }
	  }
	  if(all(!correctCoding)){
	  	write(paste("No transcripts with coding basepairs divisible by 3:", gene,ambiguousGeneNumber, sep=" "), file=failedLog, append=T) 	
	  	return(NULL)
	  }else{
	  	
	    subGenes = subGenes[which(correctCoding),,drop=F]	
	  	coding_bools = coding_bools[which(correctCoding),,drop=F]
	  }
	  
	}else{
		
		subGenes = subGenes[which(!switchEnsl),,drop=F]	
		coding_bools = coding_bools[which(!switchEnsl),,drop=F]
	}

	pos_sums = colSums(coding_bools)
	ambiguous_positions = colnames(coding_bools)[which(pos_sums != 0 & pos_sums != nrow(coding_bools))]
	
	###Get background mutation possibilites for each transcripts
	
  bmrs = list()

	for(row in 1:nrow(subGenes)){
	 
		subGene = subGenes[row,]
    coding_bool = coding_bools[row,]

		coding_seq = paste(strsplit(fastaSeq, "")[[1]][which(coding_bool == 1)+1],collapse="")
		coding_seq_index = toupper(strsplit(fastaSeq, "")[[1]][which(coding_bool == 1)+1])
		names(coding_seq_index) = names(coding_bool)[which(coding_bool == 1)]
		
		nucs = toupper(strsplit(fastaSeq, "")[[1]])
		names(nucs) = (intervalStart):(intervalEnd+1)
		categ_and_effects = matrix(0, 21, intervalLength)
		rownames(categ_and_effects) = categ_rows
		colnames(categ_and_effects) = (intervalStart+1):(intervalEnd)
		codon_positions = rep(c(1,2,3),length(which(coding_bool == 1)) / 3)
		protein = c()
     #fastaSeq2 = c()

		if(subGene[,"Strand"] == "+"){
			coding_index = 1
			if(row ==  1){to_compute = 1:(intervalLength)}else{to_compute = which(names(coding_bool) %in% ambiguous_positions)}
			if(all(coding_bool[to_compute] == 0 & length(bmrs) > 0)){
				next}
			
			for(index in to_compute){
				locus = ""
				nuc = nucs[index+1]
				if(nuc == "C"){
					if(nucs[index+2] == "G"){
						locus = "CpG"
					}else{locus = "C:G"}
				}else if(nuc == "G"){
					if(nucs[index] == "C"){
						locus = "CpG"
					}else{locus = "C:G"}
					
				}else{locus = "A:T"}

				if(coding_bool[index] == 1){### This means that the effect is either silent or nonsilent (not noncoding)
					##Index of nucleotide along coding sequence
					#nuc_index = which(names(coding_seq_index) == pos)
					if(row != 1){
						nuc_index = which(names(coding_seq_index) == names(coding_bool)[index])
						}else{
							nuc_index = coding_index
					   coding_index = coding_index + 1}
					#fastaSeq2 = c(fastaSeq2, nuc)
					##Nucleotide position within codon
					codon_pos = codon_positions[nuc_index]
			
					##Get refernce codon and amino acid
					ref_codon = coding_seq_index[(nuc_index-codon_pos+1):(nuc_index+2-codon_pos+1)]
					ref_aa = codon_table[paste(ref_codon,collapse=""),2]
					protein = c(protein, ref_aa)
					alts = nts[which(nts != nuc)]
					alt_codons = matrix("",3,2)
					
					for(alt in alts){
						
						alt_codon = ref_codon
						alt_codon[codon_pos] = alt
						alt_aa = codon_table[paste(alt_codon, collapse=""),2]
						
						mut =  paste(nuc,":",alt,sep="")
						if(mut %in% transitions){
							type = "transition"
						}else{
							type = "transversion"	
						}
						paste(locus,type,sep="")
						categ = categs[paste(locus,type,sep="")]
						if(alt_aa != ref_aa){ ##nonsilent mutation
							effect = "nonsilent"
						}else{effect = "silent"}
						categ_and_effects[paste(effect,categ,sep=""),index] = categ_and_effects[paste(effect,categ,sep=""),index] + 1
						
					}
					
				}else{
					
					### Tabulate over noncoding basepairs
					effect = "noncoding"
					alts = nts[which(nts != nuc)]
					
					for(alt in alts){
						mut =  paste(nuc,":",alt,sep="")
						if(mut %in% transitions){
							type = "transition"
						}else{
							type = "transversion"	
						}
						categ = categs[paste(locus,type,sep="")]
						categ_and_effects[paste(effect,categ,sep=""),index] = categ_and_effects[paste(effect,categ,sep=""),index] + 1
					}
					
					
				}
		
			}
			
		}else{ ##Sequence is - strand
		
			coding_index = 1
			if(row ==  1){to_compute = 1:(intervalLength)}else{to_compute = which(names(coding_bool) %in% ambiguous_positions)}

			if(all(coding_bool[to_compute] == 0 & length(bmrs) > 0)){
				next}

			for(index in to_compute){
			
				
				###Get Sequence context to detemine if CpG
				locus = ""
				nuc = nucs[index+1]
				if(nuc == "C"){
					if(nucs[index] == "G"){
						locus = "CpG"
					}else{locus = "C:G"}
				}else if(nuc == "G"){
					if(nucs[index+2] == "C"){
						locus = "CpG"
					}else{locus = "C:G"}
					
				}else{locus = "A:T"}
				
				
				if(coding_bool[index] == 1){### This means that the effect is either silent or nonsilent (not noncoding)
					##Index of nucleotide along coding sequence
					if(row != 1){nuc_index = which(names(coding_seq_index) == names(coding_bool)[index])}else{nuc_index = coding_index
					     coding_index = coding_index + 1}
					
					#fastaSeq2 = c(fastaSeq2, nuc)
					##Nucleotide position within codon
					codon_pos = codon_positions[nuc_index]
			
					##Get refernce codon and amino acid
					ref_codon = rev(coding_seq_index[(nuc_index-codon_pos+1):(nuc_index+2-codon_pos+1)])
					ref_aa = codon_table[paste(ref_codon,collapse=""),2]
					protein = c(protein, ref_aa)
					alts = nts[which(nts != nuc)]
					alt_codons = matrix("",3,2)
					
					
					for(alt in alts){
						
						alt_codon = rev(ref_codon)
						alt_codon[codon_pos] = alt
						alt_aa = codon_table[rev(paste(alt_codon, collapse="")),2]
						
						mut =  paste(nuc,":",alt,sep="")
						if(mut %in% transitions){
							type = "transition"
						}else{
							type = "transversion"	
						}
						paste(locus,type,sep="")
						categ = categs[paste(locus,type,sep="")]
						if(alt_aa != ref_aa){ ##nonsilent mutation
							effect = "nonsilent"
						}else{effect = "silent"}
						categ_and_effects[paste(effect,categ,sep=""),index] = categ_and_effects[paste(effect,categ,sep=""),index] + 1
						
					}
					
				}else{
					
					### Tabulate over noncoding basepairs
					effect = "noncoding"
					alts = nts[which(nts != nuc)]
					
					for(alt in alts){
						mut =  paste(nuc,":",alt,sep="")
						if(mut %in% transitions){
							type = "transition"
						}else{
							type = "transversion"	
						}
						categ = categs[paste(locus,type,sep="")]
						categ_and_effects[paste(effect,categ,sep=""),index] = categ_and_effects[paste(effect,categ,sep=""),index] + 1
					}
				}
			}
		}
	  if(row == 1 & length(which(coding_bool == 1)) > 0){
	  	if( !(protein[1] %in% c("Stop","M")) | !(protein[length(protein)] %in% c("Stop","M")) ){
	  		write(paste(gene, subGene[,"Isoform"],subGene[,"Strand"], paste(protein[seq(1,length(protein),3)],collapse=""), paste("ambig",ambiguousGeneNumber,sep=""), sep="\t"), file=protein_fail, append=T)
	  	}
	  }
 
		categ_and_effects["silent7",] = colSums(categ_and_effects[1:6,])
		categ_and_effects["nonsilent7",] = colSums(categ_and_effects[8:13,])
		categ_and_effects["noncoding7",] = colSums(categ_and_effects[15:20,])
		bmrs[[length(bmrs)+1]] = categ_and_effects
	
	}	

final = bmrs[[1]]
for(ambi_pos in ambiguous_positions){
	to_convert = min(which(coding_bools[,ambi_pos] == 1))
	final[,ambi_pos] = bmrs[[to_convert]][,ambi_pos]

}
if(all(colSums(final[1:14,]) == 0)){
	write(paste("No sense protein regions found: ", gene," ",ambiguousGeneNumber, sep=""), file=protein_fail, append=T)
}
if(is.null(ambiguousGeneNumber)){
	write.table(final, file=paste(gene,".bmr.csv", sep=""), row.names=F,col.names=c(paste(gene,":chr",subGenes[1,"Chrom"],":",intervalStart+1,"-",intervalEnd,sep=""),rep("",ncol(final)-1)), sep="",quote=F)
}else{
	write.table(final, file=paste(gene,".ambig",ambiguousGeneNumber,".bmr.csv", sep=""), row.names=F,col.names=c(paste(gene,":chr",subGenes[1,"Chrom"],":",intervalStart+1,"-",intervalEnd,sep=""),rep("",ncol(final)-1)), sep="",quote=F)
}

}


#for(gene_index in 1:length(genes_to_run)){

results <- foreach (gene_index=1:length(genes_to_run), .packages=c('bedr', 'Biostrings'), .errorhandling="remove") %dopar% {
	gene = genes_to_run[gene_index]
	
	write(gene, file=logFile, append=T)
	subGenes = refGene[which(refGene$Gene == gene),]
	if(nrow(subGenes) == 0){
		write(paste("Gene not present in gene annotation (possibly due to noncanonical chromosome): ", gene, sep=""), file=failedLog, append=T)	
		#next
	  return(NULL)
	}
	subGenes = subGenes[which(subGenes$CodingEnd - subGenes$CodingStart > 0),]
	if(nrow(subGenes) == 0){
		write(paste("No coding regions detected: ", gene, sep=""), file=failedLog, append=T)	
		#next
		return(NULL)
	}
	

	###See if there is ambiguous gene annotations, that is that transcripts do not overlap 
	###and are located on different parts of the genoem
	positions = list()
	positions[[1]] = subGenes[1,]
  if(nrow(subGenes) > 1){
  	for(x in 2:nrow(subGenes)){
  		transcript = subGenes[x,]
  		found = F
  		for(y in 1:length(positions)){
  			uniqueSet = positions[[y]]
  			beginSet = unique(min(uniqueSet$TxStart))
  			endSet = unique(max(uniqueSet$TxEnd))
  			if( ((transcript$TxStart >= beginSet && transcript$TxStart <= endSet) || (transcript$TxEnd >= beginSet && transcript$TxEnd <= endSet))  && transcript$Chrom == uniqueSet$Chrom){
  					positions[[y]] = rbind(uniqueSet, transcript)
  					found = T
  					break
  			}
  		}
  		
  		if(!found){
  			positions[[length(positions)+1]] = transcript
  			
  		}
  	}
	}


	if(length(positions) == 1){
		createBMR(subGenes)
	}else{
		for(ambigNum in 1:length(positions)){
			
	  	createBMR(positions[[ambigNum]], ambigNum)
		}
	}
	
}	

