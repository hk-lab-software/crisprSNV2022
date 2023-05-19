#=======================================================================#
# Method_to_find_passenger_SNVs_for_CRISPRCas9_mediated_cancer_therapy
#=======================================================================#
# !!Prior requisites!!
# 1) Cas-OFFinder
# Download the offline version of Cas-OFFinder from
# http://www.rgenome.net/cas-offinder/portable
# according to the publishers requirements

# 2) Reference genome fasta file
# download an appropriate reference genome fasta file
# this would be an input for Cas-OFFinder.

# Set working directory
# setwd("your directory")

# Install required packages
if(!require("data.table", quietly = T))
  install.packages("data.table")

if(!require("BiocManager", quietly = T))
  install.packages("BiocManager")

if(!require("crisprScore", quietly = T)){
  BiocManager::install(version = "devel")
  BiocManager::install(version = "crisprScore")
}

# Load packages
library(data.table)
library(crisprScore)

# Loading a Genome from "BSgenome" package in R for sequence search
# We currently are available for human and mouse genomes only
# This step needs to be checked by each researcher!!
# and is a different step from Cas-OFFinder!!
# Load the identical genome version performed during mutation analysis
# Below are some examples
# Ex) Human GRCh38
#if(!require("BSgenome.Hsapiens.NCBI.GRCh38", quietly = T))
#  BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
#library(BSgenome.Hsapiens.NCBI.GRCh38)

# Ex) Mouse mm10
#if(!require("BSgenome.Mmusculus.UCSC.mm10", quietly = T))
#  BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
#library(BSgenome.Mmusculus.UCSC.mm10)

#==================================#
# Function to get reverse sequence
#==================================#
rev.comp<-function(x,rev=TRUE)
{
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb]=="A") y[bbb]<-"T"        
    if(xx[bbb]=="C") y[bbb]<-"G"        
    if(xx[bbb]=="G") y[bbb]<-"C"        
    if(xx[bbb]=="T") y[bbb]<-"A"
  }
  if(rev==FALSE) 
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T)
  {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(yy)    
}


#===================================================#
# Input
#===================================================#
# 1. Mutation file (containing only SNVs)
# Modifiy your output file obtained from tools such as Mutect2 etc.
# A data.frame where each row is a SNV and
# columns should include the below headings
# (please use identical names for each column)
# Gene: any identifier is acceptable as long as it matches with the expression file and gene set file
# Sample: should be included in the expression file
# Chromosome: Chromosome of the SNV (ex. "chr10"), should match the reference genome
# Position: Genomic coordinate of the SNV
# Reference_Allele: a reference allele (from the reference genome)
# Tumor_Allele: the SNV allele
# Tumor_count: number of reads with SNV (numeric)
# Total_depth: total sequencing depth (total number of read) (numeric)
# Genotype: heterozygous (het) or homozygous (hom)
mut_tb<-fread("CRC_cell_lines_mutation.txt", sep = "\t", data.table = F)
mut_tb$Chromosome<-sapply(1:nrow(mut_tb), function(x) paste0("chr", mut_tb$Chromosome[x]))

# 2. Expression file
# A data.frame with gene expression values where
# Samples in colnames
# Genes in rownames
# Should be in log2 scale
exp_tb<-read.delim("CRC_cell_lines_expression.txt", sep = "\t", header = T, stringsAsFactors = F)

#===================================================#
# Parameters
#===================================================#
# 1. cutoff
# choose an allele frequency theshold to use.
# This will select SNVs with allele frequencies above the calculated threshold.
# default: Het_Median_plus_MAD
# Mean: Mean of all SNVs
# Median: Median of all SNVs
# Het_Mean: Mean of heterozygous SNVs
# Het_Median: Median of heterozogous SNVs
# Het_Median_plus_MAD: Sum of median of heterozygous SNVs and median absolute deviation of heterozygous SNVs
# Het_Median_minus_MAD: Subtraction of medain absolute deviation of heterozygous SNVs from median of heterozygous SNVs
# NULL: apply no allele frequency threshold
cutoff<-"Het_Median_plus_MAD"

# 2. Gene set file
# Use our collected gene set file
# or provide a customized tab-delimited file which includes
# 2-1. To use the default gene set file
gene_set_dir<-NULL
# and select a gene set to analyze
# Essential / Haploinsufficient / Oncogenes
gene_set<-"Essential"
# !!Should use HUGO gene symbols as identifiers!!

# 2-2. To use a customized gene set file
# File should be a tab-delimited file,
# Containing one gene in each row and
# two columns:
# Gene: identifier of gene
# Gene_set: name of gene set 
gene_set_dir<-"File directory"
# You should use the Genotype parameter
# where "het" selects heterozygous SNVs,
# and "hom" selects homozygous SNVs.
# If you wish to analyze all SNVs, use c("het", "hom")
genotypes<-c("het", "hom")

# 3. Directory where reference file is located
# fasta file ex. hg38.fa
reference_path<-'E:\\KHH\\HSK Lab\\Data\\UCSC\\Reference\\mouse\\'

# 4. Directory of CasOFFinder
cas_offinder<-'"E:\\KHH\\HSK Lab\\Tools\\cas-offinder.exe"'

# 5. Ontarget score cutoff
# An on-target score cutoff (calculated with Azimuth 2.0 method),
# tool will select sgRNAs with on-target scores above this cutoff
# Default is 0.5
ontarget_cutoff<-0.5

# 6. Offtarget score cutoff
# An off-target score cutoff (calculated with CFD method)
# tool will filter out sgRNAs with off-target scores above this cutoff
# Default is 0.175
offtarget_cutoff<-0.175

# 7. Species
# Select a species
# Hsapiens(default), Mmusculus (mouse)
species<-"Hsapiens"

# The tool will produce one output tab-delimited file: 
# Final_candidate_table.txt
# with candidate SNVs and corresponding sgRNAs.
CRISPR_Tool(mutation_table = mut_tb, expression_table = exp_tb,
            cutoff = cutoff, gene_set_dir = NULL, gene_set = gene_set, genotypes = NULL,
            reference_path = reference_path, cas_offinder = cas_offinder,
            ontarget_cutoff = ontarget_cutoff, offtarget_cutoff = offtarget_cutoff, species = species)

#===================================================#
# Tool
#===================================================#
CRISPR_Tool<-function(mutation_table, expression_table,
                      cutoff = "Het_Median_plus_MAD",
                      gene_set_dir = NULL, gene_set = "Essential", genotypes = NULL,
                      reference_path, cas_offinder,
                      ontarget_cutoff, offtarget_cutoff, species = "Hsapiens"){
  mut_tb<-mutation_table
  exp_tb<-expression_table
  
  #=================================================#
  # 1. Calculate Allele Frequency of each SNV
  #=================================================#
  mut_tb$AF<-as.numeric(mut_tb$Tumor_count)/as.numeric(mut_tb$Total_depth)
  
  #=====================================================#
  # 2. Calculate allele frequency cutoff for each sample
  #=====================================================#
  sample<-unique(mut_tb$Sample)
  cutoffs<-data.frame("Sample"=sample, "Cutoff"=NA, stringsAsFactors = F)
  if(is.null(cutoff)){
    cutoffs$Cutoff<-0
  }else{
    if(cutoff=="Mean"){
      cutoffs$Cutoff<-sapply(1:nrow(cutoffs), function(x) mean(subset(mut_tb, mut_tb$Sample==cutoffs$Sample[x])$AF))
    }else if(cutoff=="Median"){
      cutoffs$Cutoff<-sapply(1:nrow(cutoffs), function(x) median(subset(mut_tb, mut_tb$Sample==cutoffs$Sample[x])$AF))
    }else if(cutoff=="Het_Mean"){
      cutoffs$Cutoff<-sapply(1:nrow(cutoffs),
                             function(x) mean(subset(mut_tb, mut_tb$Sample==cutoffs$Sample[x] & mut_tb$Genotype=="het")$AF))
    }else if(cutoff=="Het_Median"){
      cutoffs$Cutoff<-sapply(1:nrow(cutoffs),
                             function(x) median(subset(mut_tb, mut_tb$Sample==cutoffs$Sample[x] & mut_tb$Genotype=="het")$AF))
    }else if(cutoff=="Het_Median_plus_MAD"){
      cutoffs$Cutoff<-sapply(1:nrow(cutoffs), function(x){
        tmp<-subset(mut_tb, mut_tb$Sample==cutoffs$Sample[x] & mut_tb$Genotype=="het")$AF
        return(median(tmp) + mad(tmp, constant = 1))
      }
      )
    }else if(cutoff=="Het_Median_minus_MAD"){
      cutoffs$Cutoff<-sapply(1:nrow(cutoffs), function(x){
        tmp<-subset(mut_tb, mut_tb$Sample==cutoffs$Sample[x] & mut_tb$Genotype=="het")$AF
        lowerThan_median<-tmp[which(tmp<median(tmp))]
        if(length(lowerThan_median)==0){
          return(median(tmp)-mad(tmp, constant = 1))
        }else{
          return(median(tmp)-mad(lowerThan_median, constant = 1))
        }
      }
      )
    }
  }

  #=================================================#
  # 3. Retrieve gene
  #=================================================#
  if(is.null(gene_set_dir)){
    gene_set_tb<-read.delim("Gene_set_merged.txt",
                            sep = "\t", header = T, stringsAsFactors = F)
    gene_set_tb<-subset(gene_set_tb, gene_set_tb$Gene_set==gene_set)
    
    if(species=="Hsapiens"){
      genes<-unique(gene_set_tb$Human_gene)
    }else if(species=="Mmusculus"){
      genes<-unique(gene_set_tb$Mouse_gene)
    }
    
    if(gene_set=="Essential"){
      genotypes<-c("hom")
    }else{
      genotypes<-c("het", "hom")
    }
  }else{
    gene_set_tb<-read.delim(gene_set_dir, sep = "\t", header = F, stringsAsFactors = F)
    genes<-gene_set_tb$Gene
  }
  
  #==============================================================#
  # 4. Filter with gene (+ genotype) and allele frequency cutoff
  #==============================================================#
  # Gene
  mut_tb<-subset(mut_tb, mut_tb$Gene%in%genes & mut_tb$Genotype%in%genotypes)
  
  idx<-sapply(1:nrow(mut_tb), function(x){
    if(mut_tb$AF[x]>=subset(cutoffs, cutoffs$Sample==mut_tb$Sample[x])$Cutoff){
      return(x)
    }else{
      return(NA)
    }
  })
  if(any(is.na(idx))){
    idx<-idx[-which(is.na(idx))]
  }
  
  
  if(length(idx)==0){
    return(NA)
  }else{
    mut_tb<-mut_tb[idx,]
  }
  
  #=================================================#
  # 5. Filter with gene expression
  #=================================================#
  # Log2(FPKM+1)>=1
  idx<-sapply(1:nrow(mut_tb), function(x){
    if(mut_tb$Gene[x]%in%rownames(exp_tb)){
      tmp<-exp_tb[mut_tb$Gene[x], mut_tb$Sample[x]]
      if(tmp>=1){
        return(x)
      }else{
        return(NA)
      }
    }else{
      return(NA)
    }
    
  })
  
  if(any(is.na(idx))){
    idx<-idx[-which(is.na(idx))]
  }
  
  if(length(idx)==0){
    return(NA)
  }else{
    mut_tb<-mut_tb[idx,]
  }
  
  #=================================================#
  # 6. Get sequence surrounding SNV site
  #=================================================#
  # 32bp upstream
  # 32bp downstream
  # Upstream sequence
  mut_tb$Sequence_UP<-sapply(1:nrow(mut_tb), function(x){
    tmp_Chr<-mut_tb$Chromosome[x]
    
    if(species == "Hsapiens"){
      tmp_Chr<-gsub("chr", "", tmp_Chr)
    }
    
    tmp_start<-mut_tb$Position[x]-31
    tmp_end<-mut_tb$Position[x]
    
    tryCatch(expr = {
      return(as.character(getSeq(get(species), tmp_Chr, start = tmp_start, end = tmp_end)))
    },
    error = function(e){return(NA)})
  })
  
  if(any(is.na(mut_tb$Sequence_UP))){
    mut_tb<-mut_tb[-which(is.na(mut_tb$Sequence_UP)),]
  }
  
  # Downstream sequence
  mut_tb$Sequence_DW<-sapply(1:nrow(mut_tb), function(x){
    tmp_Chr<-mut_tb$Chromosome[x]
    
    if(species == "Hsapiens"){
      tmp_Chr<-gsub("chr", "", tmp_Chr)
    }
    
    tmp_start<-mut_tb$Position[x]
    tmp_end<-mut_tb$Position[x]+31
    
    tryCatch(expr = {
      return(as.character(getSeq(get(species), tmp_Chr, start = tmp_start, end = tmp_end)))
    },
    error = function(e){return(NA)})
  })
  
  if(any(is.na(mut_tb$Sequence_DW))){
    mut_tb<-mut_tb[-which(is.na(mut_tb$Sequence_DW)),]
  }
  
  # Convert sequence to containing SNV
  # Upstream sequence
  mut_tb$SequenceUP_snv<-sapply(1:nrow(mut_tb), function(x){
    tmp<-strsplit(mut_tb$Sequence_UP[x], split = "")[[1]]
    tmp[32]<-mut_tb$Tumor_Allele[x]
    tmp<-paste0(tmp, collapse = "")
    return(tmp)
  })
  
  # Downstream sequence
  mut_tb$SequenceDW_snv<-sapply(1:nrow(mut_tb), function(x){
    tmp<-strsplit(mut_tb$Sequence_DW[x], split = "")[[1]]
    tmp[1]<-mut_tb$Tumor_Allele[x]
    tmp<-paste0(tmp, collapse = "")
    return(tmp)
  })
  
  #=================================================#
  # 7. Search for PAM sites from sequence
  #=================================================#
  # SpCas9 PAM: NGG
  Get_PAM_UpReverse<-function(SequenceUp_snv, First_Index, Last_Index, PAM, Cut_Index){
    SequenceUp_snv<-strsplit(SequenceUp_snv, split = "", fixed = T)[[1]]
    PAM_sites<-"."
    non <- 0
    for(i in First_Index:Last_Index){
      search<-paste0(SequenceUp_snv[i:(i+Cut_Index-1)], collapse = "")
      if(search == PAM){
        if(non==0){
          PAM_sites<-paste0(as.character(33-(i+Cut_Index-1)), "-", as.character(33-(i)))
        }
        if(non==1){
          PAM_sites<-paste0(PAM_sites, ":", as.character(33-(i+Cut_Index-1)), "-", as.character(33-(i)))
        }
        non<-1
      }
    }
    return(PAM_sites)
  }
  Get_PAM_DwForward<-function(SequenceDw_snv, Last_Index, PAM, Cut_Index){
    SequenceDw_snv<-strsplit(SequenceDw_snv, split = "", fixed = T)[[1]]
    PAM_sites<-"."
    non <- 0
    for(i in 1:Last_Index){
      search<-paste0(SequenceDw_snv[i:(i+Cut_Index-1)], collapse = "")
      if(search == PAM){
        if(non==0){
          PAM_sites<-paste0(as.character(i), "-", as.character(i+Cut_Index-1))
        }
        if(non==1){
          PAM_sites<-paste0(PAM_sites, ":", as.character(i), "-", as.character(i+Cut_Index-1))
        }
        non<-1
      }
    }
    return(PAM_sites)
  }
  
  mut_tb$PAM_Up<-NA
  mut_tb$PAM_Dw<-NA
  
  mut_tb$PAM_Up<-sapply(1:nrow(mut_tb), function(x){
    sequence_UP<-mut_tb$SequenceUP_snv[x]
    return(Get_PAM_UpReverse(sequence_UP, First_Index = 9, Last_Index = 31, PAM = "CC", Cut_Index = 2))
  })
  
  mut_tb$PAM_Dw<-sapply(1:nrow(mut_tb), function(x){
    sequence_DW<-mut_tb$SequenceDW_snv[x]
    return(Get_PAM_DwForward(sequence_DW, Last_Index = 22, PAM = "GG", Cut_Index = 2))
  })
  
  #=================================================#
  # 8. On-target scoring with Azimuth and filtering
  #=================================================#
  forward_down_tb<-data.frame("Sample"=NA, "Gene" = NA, "Chromosome"=NA, "Position"=NA,
                              "Reference_Allele"=NA, "Tumor_Allele"=NA, "Total_depth"=NA, "Tumor_count"=NA,
                              "PAM_site_first"=NA, "PAM_site_second"=NA,
                              "Position_Ontarget_start"=NA, "Position_Ontarget_End"=NA,
                              "SNVSequence5prime3prime" = NA, "SNVSequence5prime3Prime_30mer" = NA,
                              "CutsitePosition"=NA, "CutsiteSequence"=NA,
                              stringsAsFactors=F)
  
  reverse_up_tb<-data.frame("Sample"=NA, "Gene" = NA, "Chromosome"=NA, "Position"=NA,
                            "Reference_Allele"=NA, "Tumor_Allele"=NA, "Total_depth"=NA, "Tumor_count"=NA,
                            "PAM_site_first"=NA, "PAM_site_second"=NA,
                            "Position_Ontarget_start"=NA, "Position_Ontarget_End"=NA,
                            "SNVSequence5prime3prime" = NA, "SNVSequence5prime3Prime_30mer" = NA,
                            "CutsitePosition"=NA, "CutsiteSequence"=NA,
                            stringsAsFactors=F)
  # Forward Down
  for(k in 1:nrow(mut_tb)){
    tmp<-strsplit(mut_tb$PAM_Up[k], split = ":", fixed = T)[[1]]
    if(!"."%in%tmp){
      for(x in 1:length(tmp)){
        line<-c(mut_tb$Sample[k],
                mut_tb$Gene[k],
                mut_tb$Chromosome[k],
                mut_tb$Position[k],
                mut_tb$Reference_Allele[k],
                mut_tb$Tumor_Allele[k],
                mut_tb$Total_depth[k],
                mut_tb$Tumor_count[k])
        
        idx<-as.numeric(strsplit(tmp[x], split = "-", fixed = T)[[1]])
        
        start<-as.character(mut_tb$Position[k] - (idx[1]))
        end<-as.character(mut_tb$Position[k] + (23 - idx[2]))
        
        line<-c(line, idx[1], idx[2], start, end)
        
        # Get sgRNA sequence
        # 5'>3'
        sequence<-c(mut_tb$SequenceUP_snv[k], mut_tb$SequenceDW_snv[k])
        sequence<-c(strsplit(sequence[1], split = "")[[1]],
                    strsplit(sequence[2], split = "")[[1]][-1])
        names(sequence)<-as.character(seq(from = mut_tb$Position[k]-31, to = mut_tb$Position[k]+31))
        sequence<-paste0(sequence[as.character(seq(from = as.numeric(start), to = as.numeric(end)))], collapse = "")
        sequence<-rev.comp(sequence)
        
        line<-c(line, sequence)
        
        # Get sequence 30-mer target site plus context
        # 4bp before, 3bp after
        sequence<-c(mut_tb$SequenceUP_snv[k], mut_tb$SequenceDW_snv[k])
        sequence<-c(strsplit(sequence[1], split = "")[[1]],
                    strsplit(sequence[2], split = "")[[1]][-1])
        names(sequence)<-as.character(seq(from = mut_tb$Position[k]-31, to = mut_tb$Position[k]+31))
        
        context<-paste0(sequence[as.character(seq(from = as.numeric(start)-3, to = as.numeric(end)+4))], collapse = "")
        context<-rev.comp(context)
        line<-c(line, context)
        
        # Cutsite
        # 5 base after start site
        cutsite<-as.character(as.numeric(start)+5)
        line<-c(line, cutsite)
        
        # Cutsite base
        cutsitebase<-sequence[cutsite]
        cutsitebase<-rev.comp(cutsitebase)
        line<-c(line, cutsitebase)
        
        forward_down_tb<-rbind(forward_down_tb, line)
      }
    }
  }
  forward_down_tb<-forward_down_tb[-1,]
  
  # Reverse Up
  for(k in 1:nrow(mut_tb)){
    tmp<-strsplit(mut_tb$PAM_Dw[k], split = ":", fixed = T)[[1]]
    if(!"."%in%tmp){
      for(x in 1:length(tmp)){
        line<-c(mut_tb$Sample[k],
                mut_tb$Gene[k],
                mut_tb$Chromosome[k],
                mut_tb$Position[k],
                mut_tb$Reference_Allele[k],
                mut_tb$Tumor_Allele[k],
                mut_tb$Total_depth[k],
                mut_tb$Tumor_count[k])
        
        idx<-as.numeric(strsplit(tmp[x], split = "-", fixed = T)[[1]])
        
        start<-as.character(mut_tb$Position[k] - (23 - idx[2]))
        end<-as.character(mut_tb$Position[k] + (idx[1]))
        
        line<-c(line, idx[1], idx[2], start, end)
        
        # Get sgRNA sequence
        # 5'>3'
        sequence<-c(mut_tb$SequenceUP_snv[k], mut_tb$SequenceDW_snv[k])
        sequence<-c(strsplit(sequence[1], split = "")[[1]],
                    strsplit(sequence[2], split = "")[[1]][-1])
        names(sequence)<-as.character(seq(from = mut_tb$Position[k]-31, to = mut_tb$Position[k]+31))
        sequence<-paste0(sequence[as.character(seq(from = as.numeric(start), to = as.numeric(end)))], collapse = "")
        
        line<-c(line, sequence)
        
        # Get sequence 30-mer target site plus context
        # 3bp before, 4bp after
        sequence<-c(mut_tb$SequenceUP_snv[k], mut_tb$SequenceDW_snv[k])
        sequence<-c(strsplit(sequence[1], split = "")[[1]],
                    strsplit(sequence[2], split = "")[[1]][-1])
        names(sequence)<-as.character(seq(from = mut_tb$Position[k]-31, to = mut_tb$Position[k]+31))
        
        context<-paste0(sequence[as.character(seq(from = as.numeric(start)-4, to = as.numeric(end)+3))], collapse = "")
        line<-c(line, context)
        
        # Cutsite
        # 5 base after start site
        cutsite<-as.character(as.numeric(end)-5)
        line<-c(line, cutsite)
        
        # Cutsite base
        cutsitebase<-sequence[cutsite]
        cutsitebase<-rev.comp(cutsitebase)
        line<-c(line, cutsitebase)
        
        reverse_up_tb<-rbind(reverse_up_tb, line)
      }
    }
  }
  reverse_up_tb<-reverse_up_tb[-1,]
  
  forward_down_tb$PAM_site_first<-as.numeric(forward_down_tb$PAM_site_first)
  forward_down_tb$PAM_site_second<-as.numeric(forward_down_tb$PAM_site_second)
  
  reverse_up_tb$PAM_site_first<-as.numeric(reverse_up_tb$PAM_site_first)
  reverse_up_tb$PAM_site_second<-as.numeric(reverse_up_tb$PAM_site_second)
  
  # Filter out SNVs far from PAM site
  # Filter with PAM_site_first<=14
  forward_down_tb<-subset(forward_down_tb, forward_down_tb$PAM_site_first<=14)
  reverse_up_tb<-subset(reverse_up_tb, reverse_up_tb$PAM_site_first<=14)
  
  # Calculate Azimuth score
  # Forward Down
  if(nrow(forward_down_tb)!=0){
    forward_down_tb$AzimuthScore<-NA
    if(!all(is.na(forward_down_tb))){
      predictions<-getAzimuthScores(forward_down_tb$SNVSequence5prime3Prime_30mer)
      forward_down_tb$AzimuthScore<-as.numeric(predictions$score)
    }
  }
  
  # Reverse Up
  if(nrow(reverse_up_tb)!=0){
    reverse_up_tb$AzimuthScore<-NA
    if(!all(is.na(forward_down_tb))){
      predictions<-getAzimuthScores(reverse_up_tb$SNVSequence5prime3Prime_30mer)
      reverse_up_tb$AzimuthScore<-as.numeric(predictions$score)
    }
  }
  
  # Filter with ontarget score
  forward_down_tb<-subset(forward_down_tb, forward_down_tb$AzimuthScore>=ontarget_cutoff)
  reverse_up_tb<-subset(reverse_up_tb, reverse_up_tb$AzimuthScore>=ontarget_cutoff)
  
  #====================================================#
  # 9. Run CasOFFinder to search for off-target sites
  # Will search for off-target sites under 3 mismatches
  #====================================================#
  # Make input file for Cas-OFFinder
  pattern<-"NNNNNNNNNNNNNNNNNNNNNGG"
  
  # Two files
  forward_input_file<-""
  reverse_input_file<-""
  # Forward Down
  if(nrow(forward_down_tb)!=0){
    forward_input_file<-sapply(1:nrow(forward_down_tb),
                               function(x) paste0(forward_down_tb$SNVSequence5prime3prime[x], "\t", 3))
    forward_input_file<-c(reference_path, pattern, forward_input_file)
  }
  
  # Reverse Up
  if(nrow(reverse_up_tb)!=0){
    reverse_input_file<-sapply(1:nrow(reverse_up_tb),
                               function(x) paste0(reverse_up_tb$SNVSequence5prime3prime[x], "\t", 3))
    reverse_input_file<-c(reference_path, pattern, reverse_input_file)
  }
  
  write.table(forward_input_file, paste0(gene_set, "_CasOFFinder_Forward_Down.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(reverse_input_file, paste0(gene_set, "_CasOFFinder_Reverse_Up.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  
  # Run CasOFFinder
  # Forward up
  shell(paste0(cas_offinder, " ", paste0(gene_set, "_CasOFFinder_Forward_Down.txt"), " G ", paste0(gene_set, "_CasOFFinder_Forward_Down.txt")))
  # Reverse down
  shell(paste0(cas_offinder, " ", paste0(gene_set, "_CasOFFinder_Reverse_Up.txt"), " G ", paste0(gene_set, "_CasOFFinder_Reverse_Up.txt")))
  
  #==========================================================================#
  # 10. Run CFD (off-target scoring) for sgRNAs containing off-target sites
  #==========================================================================#
  # Two files
  CasOFFinder_out_forward_down_tb<-fread(paste0(gene_set, "_CasOFFinder_Forward_Down.txt"), data.table = F)
  CasOFFinder_out_reverse_up_tb<-fread(paste0(gene_set, "_CasOFFinder_Reverse_Up.txt"), data.table = F)
  
  colnames(CasOFFinder_out_forward_down_tb)<-c("sgRNA_PAM","Chromosome", "Position", "DNA", "Direction", "Mismatches")
  colnames(CasOFFinder_out_reverse_up_tb)<-c("sgRNA_PAM","Chromosome", "Position", "DNA", "Direction", "Mismatches")
  
  # Forward down
  CasOFFinder_out_forward_down_tb$spacer<-NA
  CasOFFinder_out_forward_down_tb$pam<-NA
  CasOFFinder_out_forward_down_tb$protospacer<-NA
  for(i in 1:nrow(CasOFFinder_out_forward_down_tb)){
    tmp<-strsplit(CasOFFinder_out_forward_down_tb$sgRNA_PAM[i], split = "", fixed = T)[[1]]
    CasOFFinder_out_forward_down_tb$spacer[i]<-paste0(tmp[1:20], collapse = "")
    CasOFFinder_out_forward_down_tb$pam[i]<-paste0(tmp[21:23], collapse = "")
    CasOFFinder_out_forward_down_tb$protospacer[i]<-toupper(paste0(strsplit(CasOFFinder_out_forward_down_tb$DNA[i],
                                                                            split = "", fixed = T)[[1]][1:20],
                                                                   collapse = ""))
  }
  CasOFFinder_out_forward_down_tb$CFD_score<-getCFDScores(spacers = CasOFFinder_out_forward_down_tb$spacer,
                                                          protospacers = CasOFFinder_out_forward_down_tb$protospacer,
                                                          pams = CasOFFinder_out_forward_down_tb$pam)$score
  
  
  # Reverse up
  CasOFFinder_out_reverse_up_tb$spacer<-NA
  CasOFFinder_out_reverse_up_tb$pam<-NA
  CasOFFinder_out_reverse_up_tb$protospacer<-NA
  for(i in 1:nrow(CasOFFinder_out_reverse_up_tb)){
    tmp<-strsplit(CasOFFinder_out_reverse_up_tb$sgRNA_PAM[i], split = "", fixed = T)[[1]]
    CasOFFinder_out_reverse_up_tb$spacer[i]<-paste0(tmp[1:20], collapse = "")
    CasOFFinder_out_reverse_up_tb$pam[i]<-paste0(tmp[21:23], collapse = "")
    CasOFFinder_out_reverse_up_tb$protospacer[i]<-toupper(paste0(strsplit(CasOFFinder_out_reverse_up_tb$DNA[i],
                                                                            split = "", fixed = T)[[1]][1:20],
                                                                   collapse = ""))
  }
  CasOFFinder_out_reverse_up_tb$CFD_score<-getCFDScores(spacers = CasOFFinder_out_reverse_up_tb$spacer,
                                                          protospacers = CasOFFinder_out_reverse_up_tb$protospacer,
                                                          pams = CasOFFinder_out_reverse_up_tb$pam)$score
  
  # Filtering with off-target score
  # sgRNAs with off-targets
  CasOFFinder_out_forward_down_filtered_tb<-subset(CasOFFinder_out_forward_down_tb,
                                                   CasOFFinder_out_forward_down_tb$CFD_score>=offtarget_cutoff)
  CasOFFinder_out_reverse_up_filtered_tb<-subset(CasOFFinder_out_reverse_up_tb,
                                                 CasOFFinder_out_reverse_up_tb$CFD_score>=offtarget_cutoff)
  
  # Forward down
  idx<-c()
  for(i in 1:nrow(forward_down_tb)){
    sgRNA<-forward_down_tb$SNVSequence5prime3prime[i]
    chromosome<-forward_down_tb$Chromosome[i]
    position<-forward_down_tb$Position_Ontarget_start[i]
    
    tmp<-subset(CasOFFinder_out_forward_down_filtered_tb, CasOFFinder_out_forward_down_filtered_tb$sgRNA_PAM==sgRNA)
    if(nrow(tmp)==0){
      idx<-c(idx, i)
    }else{
      # Remove original SNV site off-targets
      tmp1<-intersect(which(tmp$Chromosome==chromosome), which(tmp$Position==(as.numeric(position)-1)))
      if(length(tmp1)!=0){
        tmp<-tmp[-tmp1,]
      }
      if(nrow(tmp)==0){
        idx<-c(idx, i)
      }
    }
  }
  forward_down_tb<-forward_down_tb[idx,]
  
  # Reverse up
  idx<-c()
  for(i in 1:nrow(reverse_up_tb)){
    sgRNA<-reverse_up_tb$SNVSequence5prime3prime[i]
    chromosome<-reverse_up_tb$Chromosome[i]
    position<-reverse_up_tb$Position_Ontarget_start[i]
    
    tmp<-subset(CasOFFinder_out_reverse_up_filtered_tb, CasOFFinder_out_reverse_up_filtered_tb$sgRNA_PAM==sgRNA)
    if(nrow(tmp)==0){
      idx<-c(idx, i)
    }else{
      # Remove original SNV site off-targets
      tmp1<-intersect(which(tmp$Chromosome==chromosome), which(tmp$Position==(as.numeric(position)-1)))
      if(length(tmp1)!=0){
        tmp<-tmp[-tmp1,]
      }
      
      if(nrow(tmp)==0){
        idx<-c(idx, i)
      }
      
    }
  }
  reverse_up_tb<-reverse_up_tb[idx,]
  
  # Add CFD score
  if(nrow(forward_down_tb)!=0){
    forward_down_tb$Number_of_offtargets<-NA
    forward_down_tb$CFD_scores<-NA
    for(i in 1:nrow(forward_down_tb)){
      sgRNA<-forward_down_tb$SNVSequence5prime3prime[i]
      chromosome<-forward_down_tb$Chromosome[i]
      position<-forward_down_tb$Position_Ontarget_start[i]
      
      tmp<-subset(CasOFFinder_out_forward_down_tb, CasOFFinder_out_forward_down_tb$sgRNA_PAM==sgRNA)
	  tmp<-unique(tmp)
      if(nrow(tmp)!=0){
		# Remove original SNV site off-targets
		tmp1<-intersect(which(tmp$Chromosome==chromosome), which(tmp$Position==(as.numeric(position)-1)))
		if(length(tmp1)!=0){
		tmp<-tmp[-tmp1,]
		}
		forward_down_tb$Number_of_offtargets[i]<-nrow(tmp)
		forward_down_tb$CFD_scores[i]<-paste0(sort(tmp$CFD_score), collapse = ";")
      }
    }
  }else{
    forward_down_tb<-cbind(forward_down_tb, data.frame(matrix(ncol = 2, nrow = 0)))
    colnames(forward_down_tb)[(ncol(forward_down_tb)-1):ncol(forward_down_tb)]<-c("Number_of_offtargets", "CFD_scores")
  }
  
  
  if(nrow(reverse_up_tb)!=0){
    reverse_up_tb$Number_of_offtargets<-NA
    reverse_up_tb$CFD_scores<-NA
    for(i in 1:nrow(reverse_up_tb)){
      sgRNA<-reverse_up_tb$SNVSequence5prime3prime[i]
      chromosome<-reverse_up_tb$Chromosome[i]
      position<-reverse_up_tb$Position_Ontarget_start[i]
      
      tmp<-subset(CasOFFinder_out_reverse_up_tb, CasOFFinder_out_reverse_up_tb$sgRNA_PAM==sgRNA)
	  tmp<-unique(tmp)
      if(nrow(tmp)!=0){
		# Remove original SNV site off-targets
		tmp1<-intersect(which(tmp$Chromosome==chromosome), which(tmp$Position==(as.numeric(position)-1)))
		if(length(tmp1)!=0){
			tmp<-tmp[-tmp1,]
		}
		reverse_up_tb$Number_of_offtargets[i]<-nrow(tmp)
        reverse_up_tb$CFD_scores[i]<-paste0(sort(tmp$CFD_score), collapse = ";")
      }
    }
  }else{
    reverse_up_tb<-cbind(reverse_up_tb, data.frame(matrix(ncol = 2, nrow = 0)))
    colnames(reverse_up_tb)[(ncol(reverse_up_tb)-1):ncol(reverse_up_tb)]<-c("Number_of_offtargets", "CFD_scores")
  }

  #=================================================#
  # Final target table
  #=================================================#
  final_tb<-rbind(forward_down_tb, reverse_up_tb)
  if(nrow(final_tb)!=0){
    # Add Allele frequency, gene and gene set
  	final_tb$AF<-as.numeric(final_tb$Tumor_count)/as.numeric(final_tb$Total_depth)
	final_tb$Gene_set<-gene_set
    final_tb$Gene<-NA
    for(k in 1:nrow(final_tb)){
      tmp<-subset(mut_tb, mut_tb$Sample==final_tb$Sample[k] & mut_tb$Chromosome==final_tb$Chromosome[k] & mut_tb$Position==final_tb$Position[k])
      final_tb$Gene[k]<-tmp$Gene
    }
  }
  
  write.table(final_tb, "Final_candidate_table.txt", sep = "\t", quote = F, row.names = F, col.names = T)
}

