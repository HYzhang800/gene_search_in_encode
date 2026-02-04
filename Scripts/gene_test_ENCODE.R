args <- commandArgs(TRUE)

defaultW <- getOption("warn")
options(warn = -1)

devtools::install_github("dzhang32/ggtranscript")
#devtools::install_github("fursham-h/factR")


suppressPackageStartupMessages({
  library(ORFik)
  library(GenomicFeatures)
  library(data.table)
  library(Biostrings)
  library(GenomicRanges)
  library(tidyverse)
  library(tidyr)
  library(rtracklayer)
  library(dplyr)
  library(ggplot2)
  library(ggtranscript)
  library(factR)
  library(RColorBrewer)
  library(reshape2)
})

options(warn = defaultW, show.error.messages = T, keep.source = T)

# options(error = quote({
#   setwd('Results/Rlog/'); # Set working directory where you want the dump to go, since dump.frames() doesn't seem to accept absolute file paths.
#   dump.frames("errorDump", to.file=TRUE, include.GlobalEnv=TRUE); # First dump to file; this dump is not accessible by the R session.
#   sink(file="error.log"); # Specify sink file to redirect all output.
#   dump.frames(); # Dump again to be able to retrieve error message and write to error log; this dump is accessible by the R session since not dumped to file.
#   cat(attr(last.dump,"error.message")); # Print error message to file, along with simplified stack trace.
#   cat('\nTraceback:');
#   cat('\n');
#   traceback(2); # Print full traceback of function calls with all parameters. The 2 passed to traceback omits the outermost two function calls.
#   sink();
#   q()}))



filtered_gtf <- args[1]
gene_info <- args[2]
metadata_dir <- args[3]
quanti_file_path <- args[4]
reference_genome <- args[5]
reference_annotation <- args[6]
organ_info_path <- args[7]
wkdir <- args[8]


#sink("Results/Rlog/Rlog.log")


cat("Input argument 1 is a ", class(filtered_gtf))

stage <- 'read input files'
# filtered_gtf <- "C:/Users/Haoyu/Desktop/snake_test/Results/filtered_gtf_list.txt"
# gene_info <- "C:/Users/Haoyu/Desktop/snake_test/Results/gencode_genes.csv"
# metadata_dir <- "E:/Encode_tsv/Metadata/lung.tsv"
# quanti_file_path <- "E:/Encode_tsv/all/"
# reference_genome <- "E:/reference_for_APTARS/hg38.fa"
# reference_annotation <- "E:/reference_for_APTARS/gencode.v29.short.gtf"
# organ_info_path <- "E:/Encode_annotation/all/organ_info_full.csv"
# wkdir <- "C:/Users/Haoyu/Desktop/snake_test/"


# get metadata


metadata_tsv <- read.table(file = metadata_dir, header = T, sep = "\t")
#metadata_tsv <- read.table("E:/Encode_annotation/Metadata/metadata_lung.tsv", header = T, sep = "\t")
organ_info <- read.csv(file = organ_info_path, header = T)
#head(organ_info)

#organ_info <- read.csv("E:/Encode_annotation/all/organ_info_full.csv")

metadata_tsv$organ <- organ_info$organ[match(metadata_tsv$Biosample.term.name, organ_info$Biosample.term.name)]

sample_size <- length(unique(metadata_tsv$File.accession))

quanti_tissue <- metadata_tsv %>%
  dplyr::select(File.accession, Experiment.accession, Biosample.term.name, organ)

#quanti_tissue$organ <- gsub(".tsv", "", quanti_tissue$organ)


# read gene to be processed


genes_in_gencode <- read.csv(gene_info, header = T)

gene_of_interest <- unique(genes_in_gencode$gene_name)


cat("Genes to be processed:", genes_in_gencode$gene_name)



# Set the reference DNA sequence
ref <- readDNAStringSet(reference_genome)   #change Add another config

# Read GENCODE V29
gencode <- rtracklayer::import(reference_annotation) %>%
  as_tibble()

# Read ENCODE GTF files
file_list <- read.table(filtered_gtf, header = F, sep = "\t")


annotation <- grep("ENCFF", file_list[,1], value = TRUE)


annotation <- paste("Merged_gtf/", annotation, sep = "")

cat("Read file_list has ", nrow(file_list), "files, ", length(annotation), "will be used for next steps.")


#setwd("F:/Data/Encode_annotation/all/encode_raw/Selected_gtf/")


f <- list()

cat("Start to read filtered ENCODE annotation files:", annotation[1], "...")





for (i in 1:length(annotation)) {
  gtf <- rtracklayer::import(annotation[i])

  gtf <- as_tibble(gtf)
  f[[i]] <- gtf
}



gtf_for_all_genes <- bind_rows(f)

gtf_transcripts <- gtf_for_all_genes %>% filter(type == "transcript")

export(gtf_transcripts, "Results/encode_annotation_transcripts.gtf", format = "gtf")

gtf_transcripts <- import("Results/encode_annotation_transcripts.gtf")

## ---- snakemake/HPC: per-gene tryCatch + debug capture (ADD) ----
dir.create("Results", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Rlog", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/debug_failed_genes", showWarnings = FALSE, recursive = TRUE)

.failed_rows <- list()

.log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%F %T"), " | ", paste(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = "Results/Rlog/per_gene_status.log", append = TRUE)
}

.safe_saveRDS <- function(obj, filepath) {
  try(saveRDS(obj, filepath), silent = TRUE)
}

.record_failure <- function(gene, stage, err_msg) {
  .failed_rows[[length(.failed_rows) + 1]] <<- data.frame(
    gene = gene,
    stage = stage,
    error_message = err_msg,
    time = format(Sys.time(), "%F %T"),
    stringsAsFactors = FALSE
  )
}

.write_placeholder_outputs <- function(gene, stage, err_msg) {
  # Your Snakemake expects Results/transcript_usage_<gene>.csv to exist
  out_csv <- paste0("Results/transcript_usage_", gene, ".csv")
  if (!file.exists(out_csv)) {
    write.csv(
      data.frame(
        gene = gene,
        status = "FAILED",
        stage = stage,
        error_message = err_msg,
        stringsAsFactors = FALSE
      ),
      out_csv,
      row.names = FALSE
    )
  }
}

.save_debug_bundle <- function(gene, stage, err_msg) {
  dbg_dir <- file.path("Results", "debug_failed_genes", gene)
  dir.create(dbg_dir, showWarnings = FALSE, recursive = TRUE)

  # Always save a meta file
  meta <- list(
    gene = gene,
    stage = stage,
    error_message = err_msg,
    time = format(Sys.time(), "%F %T"),
    sessionInfo = capture.output(sessionInfo())
  )
  .safe_saveRDS(meta, file.path(dbg_dir, "meta.rds"))

  # Save key objects IF they exist (best-effort, won’t crash)
  if (exists("processing_gene_info", inherits = FALSE)) .safe_saveRDS(processing_gene_info, file.path(dbg_dir, "processing_gene_info.rds"))
  if (exists("gene_region", inherits = FALSE))         .safe_saveRDS(gene_region, file.path(dbg_dir, "gene_region.rds"))
  if (exists("overlapping_trans", inherits = FALSE))   .safe_saveRDS(overlapping_trans, file.path(dbg_dir, "overlapping_trans.rds"))
  if (exists("gtf_tibble", inherits = FALSE))          .safe_saveRDS(gtf_tibble, file.path(dbg_dir, "gtf_tibble.rds"))

  # This is the most important for your predictNMD crash:
  if (exists("my_id_whole_gtf", inherits = FALSE))     .safe_saveRDS(my_id_whole_gtf, file.path(dbg_dir, "my_id_whole_gtf.rds"))
  if (exists("my_id_grange", inherits = FALSE))        .safe_saveRDS(my_id_grange, file.path(dbg_dir, "my_id_grange.rds"))

  # If you already exported the per-gene gtf (you do just before predictNMD), copy that path into meta:
  gtf_path <- paste0("Results/my_id_whold_gtf_", gene, ".gtf")
  if (file.exists(gtf_path)) {
    # Keep a copy under debug dir as well (so you don’t lose it if Results gets cleaned)
    file.copy(gtf_path, file.path(dbg_dir, paste0("my_id_whold_gtf_", gene, ".gtf")), overwrite = TRUE)
  }

  # Traceback log (useful on HPC)
  tb <- paste(capture.output(traceback(2)), collapse = "\n")
  cat(tb, "\n", file = file.path("Results/Rlog", paste0("traceback_", gene, ".log")), append = FALSE)
}


for (g in 1:length(gene_of_interest)) {
  
processing_gene <- gene_of_interest[g]

cat("Processing gene", processing_gene)

 stage <- "init"
  .log_msg("START ", processing_gene)

tryCatch({ 
processing_gene_info <- genes_in_gencode %>% filter(gene_name == processing_gene)

gene_name <- processing_gene
gene_id <- processing_gene_info$gene_id ##############CHECK ID FIELD

gene_strand <- processing_gene_info$strand
gene_start <- processing_gene_info$start
gene_end <- processing_gene_info$end
gene_chrom <- processing_gene_info$seqname

gene_region <- GRanges(
  seqnames = gene_chrom,
  ranges = IRanges(start = as.numeric(gene_start), end = as.numeric(gene_end)),
  strand = Rle(strand(gene_strand))
)


## If target gene is an overlapping locus (with another gene), only reads fully contained in the locus will be kept
## Otherwise, all reads with at least 50 bps overlapping will be kept
stage <- 'check overlapping'
if(is.na(processing_gene_info$tag)){
  overlapping_trans <- subsetByOverlaps(gtf_transcripts, gene_region, minoverlap = 50) 

}else if(processing_gene_info$tag == "overlapping_locus"){
  overlapping_trans <- subsetByOverlaps(gtf_transcripts, gene_region, type = "within")
} else{
  overlapping_trans <- subsetByOverlaps(gtf_transcripts, gene_region, minoverlap = 50)
}

overlapping_trans <- unique(overlapping_trans$transcript_id)

# If no transcripts detected for this gene, move to the next gene 
if(length(overlapping_trans) == 0){
  
  next
}


  
# keep unique 'exons'
stage <- 'keep unique exons'
gtf <- gtf_for_all_genes %>% filter(type == "exon") %>%
  filter(transcript_id %in% overlapping_trans) %>%
  unique()


# rank exons for each transcripts from 5' to 3'



gtf <- gtf %>% group_by(transcript_id) %>% group_split()


if(gene_strand == "-"){
for (i in 1:length(gtf)){
  gtf[[i]] <- gtf[[i]][order(-gtf[[i]]$start),] 
}
} else {
  for (i in 1:length(gtf)){
    gtf[[i]] <- gtf[[i]][order(gtf[[i]]$start),] 
  
  }
}

gtf_tibble <- bind_rows(gtf)

gtf_tibble <- gtf_tibble %>% group_by(transcript_id) %>% mutate(exon_number = row_number())



###############



# Create GRangesList from gtf_tibble
grl <- makeGRangesListFromDataFrame(gtf_tibble, split.field = "transcript_id", names.field = "transcript_id")

stage <- 'get transcript/ORF sequences'
# Extract transcript sequences
seqs <- extractTranscriptSeqs(ref, grl)
#writeXStringSet(seqs, "/mnt/e/Encode_tsv/all/transcript_seq.fa")


# Convert the sequences to a data frame

gene_seq <- data.frame(
  annot_transcript_id = names(seqs),
  transcript_sequence = paste(seqs)
)

orf_orfik <- findMapORFs(grl, seqs,longestORF = T,minimumLength = 0, startCodon = startDefinition(6), stopCodon = stopDefinition(1),
                         groupByTx = T) #GRangesList and DNAStringSet have to be in the same order 

orftibble <- as_tibble(orf_orfik)
orftibble <- unique(orftibble)

orftibble <- orftibble %>% group_by(names) %>% mutate(total_length = rowsum(width, group)) %>% ungroup()
orftibble <- orftibble %>% separate(names, c('transcript_id', 'n_orf'), sep = '_')
orftibble$total_length <- as.integer(orftibble$total_length)


#select longest ORF for each transcript
orf_longest <- orftibble %>% group_by(transcript_id) %>% top_n(1,total_length)


orf_longest <- unique(orf_longest)

#make gtf structure
all_cds_gtf <- orf_longest %>% dplyr::select(seqnames:transcript_id)


#Used quantified transcripts, so unquantified suspicous transcripts were not included 


all_cds_gtf$type <- "CDS"

#gene_cds_gtf <- gene_cds_gtf[,c(1,2,3,4,5,7,6)]



all_cds_gtf <- unique(all_cds_gtf)




whole_gtf <- dplyr::bind_rows(gtf_tibble, all_cds_gtf)
whole_gtf$gene_name <- gene_name
whole_gtf$gene_id <- gene_id


talon_id <- whole_gtf$transcript_id %>% unique()


#####
files <- list.files(path = quanti_file_path)
quantification <- grep("tsv*", files, value = TRUE)

qv <- metadata_tsv$File.accession
qv <- paste(qv, ".tsv", sep = "")
quantification <- intersect(quantification, qv)


if (!endsWith(quanti_file_path, "/")) {
  # If it doesn't end with "/", add "/" to the end
  quanti_file_path <- paste(quanti_file_path, "/", sep = "")
}

quantification_path <- paste(quanti_file_path, quantification, sep = "")

f <- list()



#select previously identified transcripts and calculate their usage
stage <- 'calculate TPM'
for (i in 1:length(quantification)) {
  f[[i]] <- read.table(quantification_path[i], sep = '\t', header = TRUE)
  f[[i]] <- f[[i]] %>% filter(annot_transcript_id %in% talon_id) %>%
    filter(transcript_novelty != "Genomic")  #remoce Genomic transcripts
  
  if (nrow(f[[i]]) != 0) {
    f[[i]]$sample <- colnames(f[[i]])[12]
    colnames(f[[i]])[12] <- 'read_counts'
    f[[i]]$proportion <- f[[i]]$read_counts / sum(f[[i]]$read_counts)
    colnames(f[[i]])[14] <- paste('proportion', f[[i]]$sample[1])
    f[[i]]$file <- quantification[i]
  }
}


f <- f[sapply(f, function(df) nrow(df) > 0)]


sample_detected <- length(f)


quanti_info <- bind_rows(f)
quanti_info[is.na(quanti_info)] <- 0
quanti_info$file <- gsub(".tsv", "", as.character(quanti_info$file))



# Add tissue information from quanti_tissue data frame
quanti_info$tissue <- quanti_tissue$Biosample.term.name[match(quanti_info$file, quanti_tissue$File.accession)]

# Add organ information from quanti_tissue data frame
quanti_info$organ <- quanti_tissue$organ[match(quanti_info$file, quanti_tissue$File.accession)]
########

#get gene length



gencode_gene <- gencode %>% filter(type == "gene")

#get TPM for genes
#spike-ins removed
#'Antisense' and 'Intergenic' genes removed
f <- list()


for (i in 1:length(quantification)) {
  df <- read.table(quantification_path[i], sep = '\t', header = TRUE)
  
  df <- df %>% filter(transcript_novelty != "Genomic")
  colnames(df)[12] <- "read_counts"
  
  # For CLN3, deal with overlapping readthrough gene.
  # df$annot_gene_id <- ifelse(df$annot_transcript_id %in% cln3_id, gene_id, 
  #                            ifelse(df$annot_transcript_id %in% fusion_id, gene_id, df$annot_gene_id))
  
  
  df$annot_gene_id <- ifelse(df$annot_transcript_id %in% talon_id, gene_id, df$annot_gene_id)
  df <- df %>% 
    filter(gene_novelty == "Known")
  
  df_gencode <- df %>% 
    filter(annot_gene_id %in% gencode_gene$gene_id)
  
  df_gencode$gene_length <- gencode_gene$width[match(df_gencode$annot_gene_id, gencode_gene$gene_id)]
  
  df_mis <- df %>% filter(!annot_gene_id %in% gencode_gene$gene_id)
  
  #remove spike-ins
  df_mis <- df_mis %>%
    filter(!grepl("Spikein", annot_gene_id))
  
  dup_id <- df_mis$annot_gene_id[duplicated(df_mis$annot_gene_id)]
  dup <- df_mis %>% 
    filter(annot_gene_id %in% dup_id) %>% 
    group_by(annot_gene_id) %>%
    summarise(gene_length = max(length))
  
  df_mis$gene_length <- dup$gene_length[match(df_mis$annot_gene_id, dup$annot_gene_id)]
  df_mis$gene_length <- ifelse(df_mis$annot_gene_id %in% dup$annot_gene_id, df_mis$gene_length, df_mis$length)
  
  df <- bind_rows(df_gencode, df_mis)
  df <- df %>% 
    dplyr::select(annot_gene_id, read_counts, gene_length)
  df <- df %>%
    group_by(annot_gene_id, gene_length) %>%
    summarise(read_counts = sum(read_counts))
  
  
  df$length_kb <- df$gene_length / 1000
  df$rpk <- df$read_counts / df$length_kb
  total <- sum(df$rpk)
  df$TPM <- (df$rpk / total) * 1e6
  df <- df %>% dplyr::select(-c(length_kb, rpk))
  df$file <- gsub(".tsv", "", as.character(quantification[[i]]))
  
  f[[i]] <- df
}

gene_quanti <- bind_rows(f)

sample_depth <- gene_quanti %>% filter(annot_gene_id == gene_id)


sample_depth$Biosample.term.name <- quanti_tissue$Biosample.term.name[match(sample_depth$file, quanti_tissue$File.accession)]

organs <- unique(quanti_tissue$organ)
organ_depth <- list()


#make separate rows for same tissues in different organs
for (i in 1:length(organs)) {
  tissues <- quanti_tissue %>% filter(organ == organs[i])
  df <- sample_depth %>% dplyr::filter(Biosample.term.name %in% tissues$Biosample.term.name)
  df$organ <- organs[i]
  organ_depth[[i]] <- df
}

sample_depth <- bind_rows(organ_depth)

sample_depth_gene <- sample_depth %>% 
  filter(annot_gene_id == gene_id)

sample_depth_gene$gene <- gene_name

#summarise CLN3 TPM for each organ
sample_info <- sample_depth_gene %>%
  group_by(organ) %>%
  summarize(
    TPM_mean = mean(TPM),
    TPM_median = median(TPM),
    TPM_max = max(TPM),
    TPM_min = min(TPM),
    reads_mean = mean(read_counts),
    reads_median = median(read_counts),
    reads_max = max(read_counts),
    reads_min = min(read_counts),
    num_samples = n(),
    tissues = paste(unique(Biosample.term.name), collapse = ', ')
  )



whole_gtf <- whole_gtf %>% filter(transcript_id %in% quanti_info$annot_transcript_id) %>% unique()



# Get CDS coordinates for each transcript
cds_coord <- whole_gtf %>% filter(type == "CDS") %>% GRanges()

cds_coord <- split(cds_coord,cds_coord$transcript_id)

# Extract CDS sequences
cds <- extractTranscriptSeqs(ref, cds_coord)

# Create a data frame with annotated transcript ID and CDS sequence
gene_cds <- data.frame(annot_transcript_id = names(cds), cds_sequence = paste(cds))
#3 TALON transcripts do not have CDSs


# Merge quanti_cds with gene_cds based on annot_transcript_id
quanti_cds <- quanti_info %>% dplyr::select(-c('read_counts', 'sample'))
#quanti_cds <- merge(quanti_cds, gene_cds, by = 'annot_transcript_id')

quanti_cds <- left_join(quanti_cds, gene_cds, by = 'annot_transcript_id')

quanti_cds$cds_sequence <- ifelse(is.na(quanti_cds$cds_sequence), "no", quanti_cds$cds_sequence)

# Retrieve unique CDS sequences and calculate rank and number of amino acids
gene_cds <- quanti_cds %>%
  dplyr::select(cds_sequence) %>%
  unique() %>%
  arrange(desc(str_length(cds_sequence))) %>%
  mutate(rank = 1:n()) %>%
  mutate(num_aa = ifelse(cds_sequence != "no", (str_length(cds_sequence) - 3) / 3, 0)) %>%
  mutate(orf_id = ifelse( num_aa != 0, paste0(gene_name,"_", rank, '_', num_aa, 'aa'),paste0(gene_name,"_", rank, '_No_ORF') )) %>%
  dplyr::select(-c('rank', 'num_aa'))

# Merge quanti_cds with gene_cds based on cds_sequence
quanti_cds <- merge(quanti_cds, gene_cds, by = "cds_sequence")

#####determine transcript novelty

annot_gencode <- gencode %>% 
  filter(gene_name %in% quanti_info$annot_gene_name) 
  

rtracklayer::export(annot_gencode, paste("Results/gencode_gene_",processing_gene,".gtf", sep = ""), format = "gtf")


annot_gencode <- annot_gencode %>% filter(type == "exon")

encode_lr <- gtf_tibble

encode_lr$exon_number <- as.character(encode_lr$exon_number)

# Create coordinates column
encode_lr$coordinates <- paste0(encode_lr$start, "-", encode_lr$end)
annot_gencode$coordinates <- paste0(annot_gencode$start, "-", annot_gencode$end)

# Merge encode and gencode transcripts
transcripts_whole <- bind_rows(encode_lr, annot_gencode) %>%
  group_by(coordinates) %>%
  mutate(exon_id = cur_group_id())

# Assign exon_id to gencode and encode
annot_gencode$exon_id <- transcripts_whole$exon_id[match(annot_gencode$coordinates, transcripts_whole$coordinates)]
encode_lr$exon_id <- transcripts_whole$exon_id[match(encode_lr$coordinates, transcripts_whole$coordinates)]


encode_list <- encode_lr %>%
  group_by(transcript_id) %>%
  group_split()

# Assign IDs to encode transcripts using unique exon IDs
encode_lr <- bind_rows(lapply(encode_list, function(df) {
  df$transcript_num <- paste(df$exon_id, collapse = "-")
  df
}))


annot_gencode_list <- annot_gencode %>%
  group_by(transcript_id) %>%
  group_split()


annot_gencode <- bind_rows(lapply(annot_gencode_list, function(df) {
  df$transcript_num <- paste(df$exon_id, collapse = "-")
  df
}))


encode_lr$transcript_novelty <- ifelse(encode_lr$transcript_num %in% annot_gencode$transcript_num, "Known", "Novel")

# Assign transcript_novelty to quanti_cds
quanti_cds$transcript_novelty <- encode_lr$transcript_novelty[match(quanti_cds$annot_transcript_id, encode_lr$transcript_id)]


##UTRs novelty check
stage <- 'UTR novelty check'
txdb_known <- makeTxDbFromGFF(paste("Results/gencode_gene_",processing_gene,".gtf", sep = ""), format="gtf") #this will take some time
five_known <- fiveUTRsByTranscript(txdb_known, use.names = T)
three_known <- threeUTRsByTranscript(txdb_known, use.names = T)


five_known <- as_tibble(five_known)
three_known <- as_tibble(three_known)


##exon_name in encode dataset and gencode data are not the same !


#####make transcript database using un-merged annotation file with 'CDS'
export(whole_gtf, paste("Results/gene_encode_",processing_gene,".gtf", sep = ""), format = "gtf")

txdb <- makeTxDbFromGFF(paste("Results/gene_encode_",processing_gene,".gtf", sep = ""), format = "gtf")
####

stage <- 'check UTR novelty'
five_txdb_encode <- fiveUTRsByTranscript(txdb, use.names = T)
three_txdb_encode <- threeUTRsByTranscript(txdb, use.names = T)

five_UTR <- as_tibble(five_txdb_encode)
three_UTR <- as_tibble(three_txdb_encode)


##check 5utr novelty
five_UTR$coordinates <- paste0(five_UTR$start,"-",five_UTR$end, sep = "" )
five_known$coordinates <- paste0(five_known$start,"-",five_known$end, sep= "")


five_whole <- bind_rows(five_UTR, five_known)

#assign numeric IDs for each exons
five_whole <- five_whole %>% 
  group_by(coordinates) %>% 
  mutate(exon_id = cur_group_id())

#check which IDs are from known/long-read data
five_known$exon_id <- five_whole$exon_id[match(five_known$coordinates, five_whole$coordinates)]
five_UTR$exon_id <- five_whole$exon_id[match(five_UTR$coordinates, five_whole$coordinates)]


five_UTR <- five_UTR %>% group_by(group_name) %>% group_split()

#make numeric IDs for each UTR using exon IDs
five <- list()
for (i in 1:length(five_UTR)){
  
  df <- five_UTR[[i]]
  df$five_id <- paste(df$exon_id, collapse = "-")
  
  five[[i]] <- df
  
}


five <- bind_rows(five)


five_known <- five_known %>% group_by(group_name) %>% group_split()
gencode_five <- list()

for (i in 1:length(five_known)){
  df <- five_known[[i]]
  df$five_id <- paste(df$exon_id, collapse = "-")
  
  gencode_five[[i]] <- df
  
}

gencode_five <- bind_rows(gencode_five)

five$five_novelty <- ifelse(five$five_id %in% gencode_five$five_id, "Known", "Novel") 


#check novelty for 3utr
##########################################

three_UTR$coordinates <- paste0(three_UTR$start,"-",three_UTR$end, sep = "" )
three_known$coordinates <- paste0(three_known$start,"-",three_known$end, sep= "")


three_whole <- bind_rows(three_UTR, three_known)



three_whole <- three_whole %>% group_by(coordinates) %>% mutate(exon_id = cur_group_id())

three_known$exon_id <- three_whole$exon_id[match(three_known$coordinates, three_whole$coordinates)]
three_UTR$exon_id <- three_whole$exon_id[match(three_UTR$coordinates, three_whole$coordinates)]

three_UTR <- three_UTR %>% group_by(group_name) %>% group_split()

three <- list()
for (i in 1:length(three_UTR)){
  
  df <- three_UTR[[i]]
  df$three_id <- paste(df$exon_id, collapse = "-")
  
  three[[i]] <- df
  
}


three <- bind_rows(three)


three_known <- three_known %>% group_by(group_name) %>% group_split()
gencode_three <- list()
for (i in 1:length(three_known)){
  
  df <- three_known[[i]]
  df$three_id <- paste(df$exon_id, collapse = "-")
  
  gencode_three[[i]] <- df
  
}

gencode_three <- bind_rows(gencode_three)

three$three_novelty <- ifelse(three$three_id %in% gencode_three$three_id, "Known", "Novel") 




five_novelty <- five %>% dplyr::select(group_name,five_novelty) %>% unique()
three_novelty <- three %>% dplyr::select(group_name,three_novelty) %>% unique()

##transcript merging
stage <- 'transcript merging'

five_UTR <- as_tibble(five_txdb_encode)
three_UTR <- as_tibble(three_txdb_encode)


info <- quanti_cds


five_UTR$orf_id <- info$orf_id[match(five_UTR$group_name, info$annot_transcript_id)]

# Count the number of exons per UTR
five_UTR <- five_UTR %>% 
  group_by(group_name) %>% 
  mutate(n_5prime_exons = n())


five_UTR$coordinates <- paste0(five_UTR$start,'-',five_UTR$end)

five_UTR_list <- five_UTR %>% 
  group_by(group_name) %>% 
  group_split()

#Extract first exons and rest exons for each UTR
first_five <- list()
rest_five <- list()
for (i in 1:length(five_UTR_list)){
  first_five[[i]] <- five_UTR_list[[i]][five_UTR_list[[i]]$exon_rank == min(five_UTR_list[[i]]$exon_rank),]
  rest_five[[i]] <- five_UTR_list[[i]][five_UTR_list[[i]]$exon_rank != min(five_UTR_list[[i]]$exon_rank),]
  
}




first_five_df <- bind_rows(first_five)
rest_five_df <- bind_rows(rest_five)


rest_five_df <- rest_five_df %>% 
  group_by(coordinates) %>% 
  mutate(my_exon_id = cur_group_id())

rest_five <- rest_five_df %>% group_by(group_name) %>% group_split()

# Order rest exons by exon_rank and assign rest_id
five_prime_rest_l <- lapply(rest_five, function(df) {
  df <- df[order(df$exon_rank), ]
  my_combined_id <- df$my_exon_id
  df$rest_id <- paste(my_combined_id, collapse = '_')
  df
})

five_prime_rest <- bind_rows(five_prime_rest_l)

five_prime_id <- five_prime_rest %>% ungroup() %>% dplyr::select(group_name, rest_id) %>% unique()

five_UTR$rest_id <- five_prime_id$rest_id[match(five_UTR$group_name, five_prime_id$group_name)]

####
#assign id for first exon in 5 prime

#add id for 3 prime of first exon


five_prime_first <- first_five_df
five_prime_first$rest_id <- five_prime_id$rest_id[match(five_prime_first$group_name, five_prime_id$group_name)]
five_prime_first[is.na(five_prime_first)] <- 'no'

five_prime_first <- five_prime_first %>% 
  group_by(start) %>% 
  mutate(start_first = cur_group_id())

five_prime_first$first_start_id_rest <- paste0(five_prime_first$start_first,'_',five_prime_first$rest_id) 



five_prime_first_l <- five_prime_first %>% 
  group_by(first_start_id_rest) %>% 
  group_split()

k = 1

start_merged_l <- list()

if(gene_strand == "-") {

for (m in 1:length(five_prime_first_l)){
  df <- five_prime_first_l[[m]]
  df <- df[order(df$end),]
  df <- data.frame(df)
  
  
  for (i in 1:nrow(df)){
    if(i == 1){
      df$exon_1id[i] <- k
    } else if (i >1 & df$end[i] < df[df$exon_1id == df$exon_1id[i-1],]$end[1] + 40){
      df$exon_1id[i] <- k
      
    } else if (i >1 & df$end[i] > df[df$exon_1id == df$exon_1id[i-1],]$end[1] + 40 ){   
      df$exon_1id[i] <- k+1
      k = k+1
    } 
  }
  start_merged_l[m]<- list(df)
  k= k+1
}
} else {
  
  for (m in 1:length(five_prime_first_l)){
    df <- five_prime_first_l[[m]]
    df <- df[order(df$end),]
    df <- data.frame(df)
    
    
    for (i in 1:nrow(df)){
      if(i == 1){
        df$exon_1id[i] <- k
      } else if (i >1 & df$start[i] < df[df$exon_1id == df$exon_1id[i-1],]$start[1] + 40){
        df$exon_1id[i] <- k
        
      } else if (i >1 & df$start[i] > df[df$exon_1id == df$exon_1id[i-1],]$start[1] + 40 ){   
        df$exon_1id[i] <- k+1
        k = k+1
      } 
    }
    start_merged_l[m]<- list(df)
    k= k+1  
  
}
}
  
start_merged <- bind_rows(start_merged_l)

start_merged$five_id <- paste0(start_merged$exon_1id,'_', start_merged$first_start_id_rest)

start_merged <- start_merged %>% 
  group_by(five_id) %>% 
  mutate(final_five_id = cur_group_id())

five_id <- start_merged %>% 
  ungroup() %>% 
  dplyr::select(group_name, final_five_id) %>%
  unique() %>%
  mutate(final_five_id = paste0('5UTR_', final_five_id))



###########Do the same for 3'UTRs


three_UTR$orf_id <- info$orf_id[match(three_UTR$group_name, info$annot_transcript_id)]

three_UTR <- three_UTR %>% 
  group_by(group_name) %>% 
  mutate(n_3prime_exons = n())



three_UTR$coordinates <- paste0(three_UTR$start,'-',three_UTR$end)


three_UTR_list <- three_UTR %>% 
  group_by(group_name) %>% 
  group_split()


#separate the last exons and rest exons
last_three <- list()
rest_three <- list()
for (i in 1:length(three_UTR_list)){
  last_three[[i]] <- three_UTR_list[[i]][three_UTR_list[[i]]$exon_rank == max(three_UTR_list[[i]]$exon_rank),]
  rest_three[[i]] <- three_UTR_list[[i]][three_UTR_list[[i]]$exon_rank != max(three_UTR_list[[i]]$exon_rank),]
  
}




last_three_df <- bind_rows(last_three)
rest_three_df <- bind_rows(rest_three)


##assign id for exons which are not the 'last' exons then arrange id using exon id, for each transcript
rest_three_df <- rest_three_df %>% group_by(coordinates) %>% mutate(my_exon_id = cur_group_id())

rest_three <- rest_three_df %>% group_by(group_name) %>% group_split()


three_prime_rest_l <- lapply(rest_three, function(df) {
  df <- df[order(df$exon_rank), ]
  my_combined_id <- df$my_exon_id
  df$rest_id <- paste(my_combined_id, collapse = '_')
  df
})



three_prime_rest <- bind_rows(three_prime_rest_l)

three_prime_id <- three_prime_rest %>%
  ungroup() %>% 
  dplyr::select(group_name, rest_id) %>% 
  unique()

three_UTR$rest_id <- three_prime_id$rest_id[match(three_UTR$group_name, three_prime_id$group_name)]

####

#add id for 3 prime of first exon


three_prime_last <- last_three_df
three_prime_last$rest_id <- three_prime_id$rest_id[match(three_prime_last$group_name, three_prime_id$group_name)]
three_prime_last[is.na(three_prime_last)] <- 'no'

three_prime_last <- three_prime_last %>% 
  group_by(end) %>% 
  mutate(end_last = cur_group_id())

three_prime_last$last_end_id_rest <- paste0(three_prime_last$end_last,'_',three_prime_last$rest_id) 




end_merged_l <- list()

k=1 

three_prime_last_l <- three_prime_last %>% group_by(last_end_id_rest) %>% group_split()

if(gene_strand == "-"){

for (m in 1:length(three_prime_last_l)){
  df <- three_prime_last_l[[m]]
  df <- df[order(df$start),]
  df <- data.frame(df)
  
  
  for (i in 1:nrow(df)){
    if(i == 1){
      df$exon_1id[i] <- k
    } else if (i >1 & df$start[i] < df[df$exon_1id == df$exon_1id[i-1],]$start[1] + 40){                  
      df$exon_1id[i] <- k                                            
      
    } else if (i >1 & df$start[i] > df[df$exon_1id == df$exon_1id[i-1],]$start[1] + 40 ){                    
      df$exon_1id[i] <- k+1
      k = k+1
    } 
    
  }
  end_merged_l[m]<- list(df)
  k= k+1
}
} else{
  
  for (m in 1:length(three_prime_last_l)){
    df <- three_prime_last_l[[m]]
    df <- df[order(df$start),]
    df <- data.frame(df)
    
    
    for (i in 1:nrow(df)){
      if(i == 1){
        df$exon_1id[i] <- k
      } else if (i >1 & df$end[i] < df[df$exon_1id == df$exon_1id[i-1],]$end[1] + 40){                  
        df$exon_1id[i] <- k                                            
        
      } else if (i >1 & df$end[i] > df[df$exon_1id == df$exon_1id[i-1],]$end[1] + 40 ){                    
        df$exon_1id[i] <- k+1
        k = k+1
      } 
      
    }
    end_merged_l[m]<- list(df)
    k= k+1
  }
}


end_merged <- bind_rows(end_merged_l)



end_merged$three_id <- paste0(end_merged$last_end_id_rest,'_',end_merged$exon_1id)

end_merged <- end_merged %>% 
  group_by(three_id) %>% 
  mutate(final_three_id = cur_group_id())

three_id <- end_merged %>% 
  ungroup() %>% 
  dplyr::select(group_name, final_three_id) %>%
  unique() %>% 
  mutate(final_three_id = paste0("3UTR_",final_three_id))



info$three_id <- three_id$final_three_id[match(info$annot_transcript_id, three_id$group_name)]
info$five_id <- five_id$final_five_id[match(info$annot_transcript_id, five_id$group_name)]    


info$five_id[is.na(info$five_id)] <- '5UTR_no'
info$three_id[is.na(info$three_id)] <- '3UTR_no'
info$new_id <- paste(info$orf_id,'_',info$five_id,'_',info$three_id, sep = "")

#for transcripts with no ORFs, use TALON ID as well to distinguish

if("no_ORF" %in% info$orf_id){
info$new_id <- ifelse(grep(pattern = "no_ORF", info$orf_id), 
                      paste0(info$new_id, "_", info$annot_transcript_id, sep= ""), 
                      info$new_id)
}

# Assign new UTR IDs
five_novelty$five_id <- info$five_id[match(five_novelty$group_name, info$annot_transcript_id)]
three_novelty$three_id <- info$three_id[match(three_novelty$group_name, info$annot_transcript_id)]

# Group and split 5' UTR novelty
five_novelty <- five_novelty %>% 
  group_by(five_id) %>% 
  group_split()

# Update 5' UTR novelty status
lt <- list()
for (i in 1:length(five_novelty)) {
  df_novelty <- five_novelty[[i]]
  merged_status <- ifelse('Known' %in% df_novelty$five_novelty, 'Known', 'Novel')
  df_novelty$five_novelty <- merged_status
  lt[[i]] <- as_tibble(df_novelty)
}

five_novelty <- bind_rows(lt) %>% 
  dplyr::select(five_id, five_novelty) %>% 
  unique()

# Group and split 3' UTR novelty
three_novelty <- three_novelty %>% 
  group_by(three_id) %>% 
  group_split()

# Update 3' UTR novelty status
lt <- list()
for (i in 1:length(three_novelty)) {
  df_novelty <- three_novelty[[i]]
  merged_status <- ifelse('Known' %in% df_novelty$three_novelty, 'Known', 'Novel')
  df_novelty$three_novelty <- merged_status
  lt[[i]] <- as_tibble(df_novelty)
}
three_novelty <- bind_rows(lt) %>% 
  dplyr::select(three_id, three_novelty) %>% 
  unique()

# Update ORF UTR novelty
info$five_novelty <- five_novelty$five_novelty[match(info$five_id, five_novelty$five_id)]
info$three_novelty <- three_novelty$three_novelty[match(info$three_id, three_novelty$three_id)]

# Check merged transcript novelty
merged_transcript_novelty <- info %>% 
  dplyr::select(new_id, annot_transcript_id,transcript_novelty)

merged_transcript_novelty <- merged_transcript_novelty %>% 
  group_by(new_id) %>% 
  group_split()

# Update merged transcript novelty status
lt <- list()
for (i in 1:length(merged_transcript_novelty)) {
  df_novelty <- merged_transcript_novelty[[i]]
  merged_status <- ifelse('Known' %in% df_novelty$transcript_novelty, 'Known', 'Novel')
  df_novelty$merged_novelty <- merged_status
  lt[[i]] <- as_tibble(df_novelty)
}
merged_transcript_novelty <- bind_rows(lt) %>% 
  dplyr::select(new_id, merged_novelty)

#add merged novelty to info
info$transcript_novelty <- merged_transcript_novelty$merged_novelty[match(info$new_id, merged_transcript_novelty$new_id)]



trans_merged_quanti <- info %>% 
  dplyr::select(new_id, cds_sequence, colnames(info)[grep('rep', colnames(info))])
trans_merged_quanti[is.na(trans_merged_quanti)] <- 0

# Merge rows and calculate occurrence
trans_merged_quanti <- trans_merged_quanti %>% 
  group_by(new_id, cds_sequence) %>% 
  dplyr::summarise_each(list(sum))

for (i in 1:nrow(trans_merged_quanti)) {
  trans_merged_quanti[i, 'occurrence'] <- length(which(trans_merged_quanti[i, 3:(2 + length(grep('rep', colnames(trans_merged_quanti))))] != 0))
}


info$occurrence <- trans_merged_quanti$occurrence[match(info$new_id, trans_merged_quanti$new_id)]

# Filter merged ORF data based on occurrence
trans_merged_filtered <- trans_merged_quanti %>% 
  filter(occurrence >= 3)

# Update transcript novelty in filtered ORF data
trans_merged_filtered$transcript_novelty <- merged_transcript_novelty$merged_novelty[match(trans_merged_filtered$new_id, merged_transcript_novelty$new_id)]


##get transcript level quantification
stage <- 'get transcript level quantification'

# Replace 0 with NA
transcript_info <- trans_merged_quanti
transcript_info[transcript_info == 0] <- NA

# Calculate mean, min, median, and max proportions
transcript_info$mean_prop <- apply(transcript_info[, grep('rep', colnames(transcript_info))], 1, function(x) mean(na.omit(as.numeric(x))))
transcript_info$min_prop <- apply(transcript_info[, grep('rep', colnames(transcript_info))], 1, function(x) min(na.omit(as.numeric(x))))
transcript_info$median_prop <- apply(transcript_info[, grep('rep', colnames(transcript_info))], 1, function(x) median(na.omit(as.numeric(x))))
transcript_info$max_prop <- apply(transcript_info[, grep('rep', colnames(transcript_info))], 1, function(x) max(na.omit(as.numeric(x))))

# Select columns and order by occurrence in descending order
transcript_info <- transcript_info %>% dplyr::select(new_id, mean_prop, min_prop, median_prop, max_prop, occurrence)
transcript_info <- transcript_info[order(-transcript_info$occurrence), ]

# Extract ENST_ID for known transcripts
enstid <- info %>% filter(transcript_novelty == "Known") %>% dplyr::select(new_id, annot_transcript_id) %>% unique()
transcript_info$ENST_ID <- enstid$annot_transcript_id[match(transcript_info$new_id, enstid$new_id)]

####


whole_gtf$transcript_novelty <- info$transcript_novelty[match(whole_gtf$transcript_id, info$annot_transcript_id)]



gene_cds_gtf <- whole_gtf %>% dplyr::filter(type == 'CDS')
gene_exons <- whole_gtf %>% dplyr::filter(type == 'exon')


my_final_id <- unique(info$new_id)
final_five_id <- unique(five_id$final_five_id)
final_three_id <- unique(three_id$final_three_id)

#As CLN3 is on - strand, "end" of the 5' most exon is 5' end
five_end_coord <- sapply(final_five_id, function(id) {
  orf_5 <- info %>% filter(five_id == id)
  gene_exons_5 <- gene_exons %>% filter(transcript_id %in% orf_5$annot_transcript_id)
  five_start <- round(mean(gene_exons_5[gene_exons_5$exon_number == min(as.integer(gene_exons_5$exon_number)),]$end), digits = 0)
  five_start
})

three_end_coord <- sapply(final_three_id, function(id) {
  orf_3 <- info %>% filter(three_id == id)
  gene_exons_3 <- gene_exons %>% filter(transcript_id %in% orf_3$annot_transcript_id)
  three_end <- round(mean(gene_exons_3[gene_exons_3$exon_number == max(as.integer(gene_exons_3$exon_number)),]$start), digits = 0)
  three_end
})

five_coord <- data.frame(final_five_id, five_end_coord, row.names = NULL)
three_coord <- data.frame(final_three_id, three_end_coord, row.names = NULL)


result_gtf <- list()
result_cds <- list()



for (i in 1:length(my_final_id)) {
  orf_example <- info %>% filter(new_id == my_final_id[i])
  
  
  gene_exons_example <- gene_exons %>% filter(transcript_id %in% orf_example$annot_transcript_id )
  
  gene_cds_example <- gene_cds_gtf %>% filter(transcript_id %in% orf_example$annot_transcript_id )
  
  #check merged transcript novelty
  transcript_status <- gene_exons_example$transcript_novelty
  merged_status <- ifelse('Known' %in% transcript_status, 'Known', 'Novel')
  
  
  gene_exons_example <- gene_exons_example %>% group_by(transcript_id) %>% group_split()
  merge_gtf <- gene_exons_example[[1]]
  merge_gtf <- merge_gtf[order(as.integer(merge_gtf$exon_number)),]
  
  
  five_start <- five_coord[five_coord$final_five_id == unique(orf_example$five_id),]$five_end_coord
  three_end <- three_coord[three_coord$final_three_id == unique(orf_example$three_id),]$three_end_coord
  
  merge_gtf$end[1] <- ifelse(length(five_start) > 0, as.integer(five_start), merge_gtf$end[1])
  
  merge_gtf$start[nrow(merge_gtf)] <- ifelse(length(three_end) > 0, as.integer(three_end), merge_gtf$start[nrow(merge_gtf)])
  
  merged_id <- info[info$annot_transcript_id == merge_gtf$transcript_id[1],]$new_id
  merge_gtf$transcript_id <- my_final_id[i]
  
  merge_gtf$transcript_status <- merged_status
  
  
  
  gene_cds_example <- gene_cds_example %>% group_by(transcript_id) %>% group_split()
  if(length(gene_cds_example) != 0){
    merged_cds <- gene_cds_example[[1]]
    merged_cds$transcript_status <- merged_status
    merged_cds$transcript_id <- my_final_id[i]
  } else {
    merged_cds <- list()
  }
  
  result_gtf[[i]] <- merge_gtf
  result_cds[[i]] <- merged_cds
  
  
}

my_id_gtf <- unique(bind_rows(result_gtf))
my_id_cds <- unique(bind_rows(result_cds))

my_id_whole_gtf <- rbind(my_id_gtf,my_id_cds)




my_id_filtered <- my_id_whole_gtf %>% filter(transcript_id %in% trans_merged_filtered$new_id)




##NMD prediction
stage <- 'NMD prediction'

                                  
export(my_id_whole_gtf, paste("Results/my_id_whold_gtf_",processing_gene,".gtf", sep = ""), format = "gtf")
my_id_grange <- import(paste("Results/my_id_whold_gtf_",processing_gene,".gtf", sep = ""))

bad <- my_id_grange[end(my_id_grange) < start(my_id_grange)]

if (length(bad) > 0) {
  # save bad ranges for debugging

  saveRDS(bad, file.path(dbg_dir, "bad_ranges_granges.rds"))
  stop("Invalid ranges detected before predictNMD (end < start).")
}

                                  
nmd_info <- predictNMD(my_id_grange)



nmd_info <- as.data.frame(nmd_info)

##ORF usage
stage <- 'get ORF usage'
trans_merged_filtered <- trans_merged_filtered %>% ungroup()
trans_merged_filtered$orf_id <- info$orf_id[match(trans_merged_filtered$new_id, info$new_id)]


orf_prop <- trans_merged_filtered %>% dplyr::select('orf_id', colnames(trans_merged_filtered)[grep('rep', colnames(trans_merged_filtered))])
#orf_prop <- info %>% dplyr::select('orf_id', colnames(info)[grep('rep', colnames(info))])


column_sums <- colSums(orf_prop[,2:(sample_detected+1)], na.rm = TRUE)

# Create a new data frame with column names as variables and sums as values
sum_df <- data.frame(variable = names(column_sums), sum_value = column_sums)



orf_prop <- orf_prop %>% group_by(orf_id) %>% summarise_all(sum)

# Check occurrence of ORFs
orf_prop$occurrence <- rowSums(orf_prop[, 2:(1+length(grep('rep', colnames(trans_merged_quanti))))] != 0)

orf_prop_long <- reshape2::melt(orf_prop)
orf_prop_long <- orf_prop_long %>% dplyr::filter(variable != 'occurrence')
orf_prop_long$occurrence <- orf_prop$occurrence[match(orf_prop_long$orf_id, orf_prop$orf_id)]

orf_prop_long[orf_prop_long == 0] <- NA
orf_prop_long <- na.omit(orf_prop_long)

#ORF novelty and usage
grl_gencode <- makeGRangesListFromDataFrame(annot_gencode, split.field = "transcript_id", names.field = "transcript_id")

seqs_gencode <- extractTranscriptSeqs(ref, grl_gencode)

gencode_gene_seq <- data.frame(annot_transcript_id = names(seqs_gencode), transcript_sequence = paste(seqs_gencode))

orf_gencode <- findMapORFs(grl_gencode, seqs_gencode, longestORF = TRUE, minimumLength = 0, startCodon = startDefinition(6), stopCodon = stopDefinition(1), groupByTx = TRUE) #GRangesList and DNAStringSet have to be in the same order 

orftibble <- as_tibble(orf_gencode) %>% unique()

orftibble <- orftibble %>%
  group_by(names) %>%
  mutate(total_length = sum(width)) %>%
  ungroup() %>%
  separate(names, c('transcript_id', 'n_orf'), sep = '_') %>%
  mutate(total_length = as.integer(total_length))

orf_longest <- orftibble %>%
  group_by(transcript_id) %>%
  top_n(1, total_length) %>%
  unique()

grl_gencode_cds <- makeGRangesListFromDataFrame(orf_longest, split.field = "transcript_id", names.field = "transcript_id")

gencode_cds <- extractTranscriptSeqs(ref, grl_gencode_cds)

gencode_gene_cds <- data.frame(gencode_transcript_id = names(gencode_cds), cds_sequence = paste(gencode_cds))

novel_orf <- gene_cds %>%
  filter(!cds_sequence %in% gencode_gene_cds$cds_sequence)


# Calculate mean, min, median, and max proportions for each ORF
orf_prop[orf_prop == 0] <- NA

for (i in 1:nrow(orf_prop)) {
  prop <- as.numeric(orf_prop[i, grep('rep', colnames(orf_prop))])
  prop <- na.omit(prop)
  
  orf_prop[i, 'mean_prop'] <- mean(prop)
  orf_prop[i, 'min_prop'] <- min(prop)
  orf_prop[i, 'median_prop'] <- median(prop)
  orf_prop[i, 'max_prop'] <- max(prop)
}


#Summary for ORF of transcripts with occurrence >=3

orf_prop_sum <- orf_prop %>%
  dplyr::select(orf_id, mean_prop, min_prop, median_prop, max_prop, occurrence) %>%
  arrange(desc(occurrence))



orf_prop_sum$novelty <- ifelse(orf_prop_sum$orf_id %in% novel_orf$orf_id, "Novel", "Known")



#transcript class determination
stage <- 'transcript class determination'                               

info$five_novelty <- five_novelty$five_novelty[match(info$five_id, five_novelty$five_id)]
info$three_novelty <- three_novelty$three_novelty[match(info$three_id, three_novelty$three_id)]

info$nmd <- nmd_info$is_NMD[match(info$new_id, nmd_info$transcript)]


info$orf_novelty <- ifelse(info$cds_sequence %in% gencode_gene_cds$cds_sequence, "Known", "Novel")


gene_cds$protein_coding <- ifelse(str_length(gene_cds$cds_sequence) >= 303, T, F)

info$protein_coding <- gene_cds$protein_coding[match(info$orf_id, gene_cds$orf_id)]

info$protein_coding <- ifelse(info$orf_id == "no_ORF", T, info$protein_coding)

info$coding <- ifelse(info$protein_coding == TRUE, "coding", "non_coding")

#Both Non-coding (<150aa) and NMD transcripts (By Sqanti) don't have ORF_seq predicted by Sqanti

novelty <- info %>%
  #filter(type == "Coding") %>%
  dplyr::select(new_id, five_novelty, three_novelty, orf_id, orf_novelty, nmd, coding, transcript_novelty)





novelty$transcript_class <- NA  # Create a new column for transcript class


# Non-coding transcripts
non_coding_known <- novelty$coding == "non_coding" & novelty$transcript_novelty == "Known"
non_coding_novel <- novelty$coding == "non_coding" & novelty$transcript_novelty == "Novel"

#NMD
nmd_known <- novelty$nmd == TRUE & novelty$coding == "coding" & novelty$transcript_novelty == "Known"
nmd_novel <- novelty$nmd == TRUE & novelty$coding == "coding" & novelty$transcript_novelty == "Novel"

# Coding transcripts
coding_known_match <- novelty$nmd == FALSE & novelty$coding == "coding" & novelty$transcript_novelty == "Known" & (novelty$five_novelty == "Known" | is.na(novelty$five_novelty)) & (novelty$three_novelty == "Known" | is.na(novelty$three_novelty))
coding_known_alt <- novelty$nmd == FALSE & novelty$coding == "coding" & novelty$transcript_novelty == "Known" & (novelty$five_novelty == "Novel" | novelty$three_novelty == "Novel")
coding_novel <- novelty$nmd == FALSE & novelty$coding == "coding" & novelty$transcript_novelty == "Novel"


test <- novelty %>% filter(novelty$nmd == FALSE & novelty$coding == "coding" & novelty$transcript_novelty == "Known" & novelty$five_novelty == "Novel" & novelty$three_novelty == "Novel")


# Classify coding_novel transcripts


# novelty$transcript_class[coding_novel & novelty$orf_novelty == "Known" & ((novelty$five_novelty == "Novel" | novelty$three_novelty == "Novel") | (novelty$five_novelty == "Novel" & novelty$three_novelty == "Novel"))] <- "Coding novel(UTR only)"
# 
# novelty$transcript_class[coding_novel & novelty$orf_novelty == "Novel" & ((novelty$five_novelty == "Novel" | novelty$three_novelty == "Novel") | (novelty$five_novelty == "Novel" & novelty$three_novelty == "Novel"))] <- "Coding novel(ORF and UTR)"
# novelty$transcript_class[coding_novel & novelty$orf_novelty == "Novel" & novelty$five_novelty == "Known" & novelty$three_novelty == "Known"] <- "Coding novel(ORF only)"
# novelty$transcript_class[coding_novel & novelty$orf_novelty == "Known" & novelty$five_novelty == "Known" & novelty$three_novelty == "Known"] <- "Coding novel(combination)"


novelty$transcript_class[coding_novel] <- "Coding novel"

# Classify non-coding and NMD transcripts
novelty$transcript_class[non_coding_known] <- "Non-coding known"
novelty$transcript_class[non_coding_novel] <- "Non-coding novel"
novelty$transcript_class[nmd_known] <- "NMD known"
novelty$transcript_class[nmd_novel] <- "NMD novel"
# Classify coding_known transcripts
novelty$transcript_class[coding_known_match] <- "Coding known(complete match)"
novelty$transcript_class[coding_known_alt] <- "Coding known(alternate 5'/3' end)"

trans_merged_filtered$transcript_class <- novelty$transcript_class[match(trans_merged_filtered$new_id, novelty$new_id)]
info$transcript_class <- novelty$transcript_class[match(info$new_id, novelty$new_id)]


#
#####################
#make output files

stage <- 'make output files'

cat("Start generating output files.")

#####Determine how many transcripts to show

#Top transcripts by occurrence
trans_rank <- transcript_info #%>% filter(occurrence >= 3) 
trans_rank <- trans_rank[with(trans_rank, order(-trans_rank$occurrence)),]

if(nrow(trans_rank) >= 30){
  trans_rank <- trans_rank[1:30,]
} else{
  
  trans_rank <- trans_rank
}



###generate n most distinctive colors
info_show <- info %>% filter(new_id %in% trans_rank$new_id)


#trans_merged_filtered$orf_id <- info$orf_id[match(trans_merged_filtered$new_id, info$new_id)]

n <- length(unique(info_show$orf_id))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))





if(n >= 15){
  n = n #generate colours for all ORFs from top 30 transcripts
  orf_show <- unique(info_show$orf_id)
} else {
  
  n = 15 #if top 30 transcripts don't contain top 15 ORFs, generate colours for top 15 ORF (by occurrence) to cover
  orf_show <- orf_prop_sum[1:15,]$orf_id
}


color_used <- col_vector[1:n]

#pie(rep(1,n), col=sample(color_used, n))


for (i in 1:length(color_used)){
  
  names(color_used)[i] <- orf_show[i]
  
}



# determine colors for transcript classes
n_types <- length(unique(info_show$transcript_class))

types <- unique(info_show$transcript_class)
color_types <- col_vector[60:(60 + n_types - 1)]


for (i in 1:length(color_types)) {
  names(color_types)[i] <- types[i]
}


#ORF usage

#select top 15 ORFs by occurrence




top_15_orf <- orf_prop_sum[1:15,]$orf_id

orf_prop_top <- orf_prop_long %>% filter(orf_id %in% top_15_orf)

orf_usage_full <- ggplot(orf_prop_top, aes(reorder(x=orf_id, occurrence), y= value, fill = orf_id)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.9, alpha=0.9, width = 0.25) +
  scale_y_continuous(n.breaks = 9, labels = scales::percent, limits = c(0, 1)) +
  ggtitle("All Tissues") +
  ylab("Proportion") +
  xlab("ORFs") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        title = element_text(face = "bold"),
        plot.title = element_text(size=24), 
        legend.position = "bottom") +
  coord_flip() +
  scale_fill_manual(values = color_used) +
  guides(fill = guide_legend(title = "ORFs", nrow = 5))


#Differential ORF usage

quanti_organ <- info %>% dplyr::select(orf_id, tissue, grep('rep', colnames(info)))

organ_info <- read.csv(organ_info_path, sep = ",")


existing_organ <- organ_info %>% filter(Biosample.term.name %in% quanti_organ$tissue)

existing_organ <- unique(existing_organ$organ)

#for same tissue types within different organs
orf_organ <- list()

for (i in 1:length(existing_organ)) {
  print(paste("Iteration:", i))
  tissues <- quanti_tissue %>% filter(organ == existing_organ[i])
  if(nrow(tissues) == 0){
    next()
  }
  df <- quanti_organ %>% dplyr::filter(tissue %in% tissues$Biosample.term.name)
  df$organ <- existing_organ[i]  
  #add organ names used in existing_organ from the new organ_info,
  #not old organ definitions from metadata_tsv (vector 'organs')
  orf_organ[[i]] <- df
}

quanti_organ <- bind_rows(orf_organ)



quanti_organ <- quanti_organ %>% dplyr::select(-tissue) %>% group_by(orf_id, organ) %>% summarise_each(list(sum))






quanti_organ_top <- quanti_organ %>% 
  ungroup() %>% 
  filter(quanti_organ$orf_id %in% top_15_orf)



quanti_organ_top <- reshape2::melt(quanti_organ_top) %>%
  filter(value != 0)

quanti_organ_top <- quanti_organ_top %>% 
  dplyr::select(-variable) %>% 
  group_by(orf_id, organ) %>% 
  summarise_each(list(sum))  

#reorder
quanti_organ_top <- quanti_organ_top %>%
  group_by(orf_id) %>%
  mutate(mean_usage = mean(value)) %>%
  arrange(desc(-mean_usage))


#rank organ by the most abundant orf

top_orf <- orf_prop_sum$orf_id[1]

organ_rank <- quanti_organ_top %>% 
  group_by(organ) %>%
  mutate(scaled_usage = value / sum(value)) %>%
  filter(orf_id == top_orf) %>%
  arrange(desc(scaled_usage)) %>%
  ungroup() %>%
  mutate(organ_rank = row_number())

quanti_organ_top$organ_rank <- organ_rank$organ_rank[match(quanti_organ_top$organ, organ_rank$organ)]

quanti_organ_top$orf_id <- factor(quanti_organ_top$orf_id, levels = unique(quanti_organ_top$orf_id[order(quanti_organ_top$value)]))


quanti_organ_top$organ <- gsub("_", " ", quanti_organ_top$organ)

#quanti_organ_top$orf_id <- factor(quanti_organ_top$orf_id, levels = rev(levels(quanti_organ_top$orf_id)))

dou_plot <- ggplot(data = quanti_organ_top, aes(x=reorder(organ, -organ_rank), y=value, fill= orf_id)) + 
  geom_bar(position = "fill", stat = "identity")+
  xlab("Organs") + 
  ylab("ORFs usage") + 
  coord_flip()+
  theme_bw()+
  scale_fill_manual(values = color_used) + 
  scale_y_continuous(labels = scales::percent) + 
  guides(fill=guide_legend(title="ORF IDs")) + 
  theme(axis.title = element_text(face = "bold"))


#Select top transcripts to show

my_id_cds <- my_id_whole_gtf %>% filter(type == "CDS") %>% as.data.frame() %>% unique()

#only keep exons for gtf
my_id_gtf <- my_id_whole_gtf %>% filter(type == "exon")


my_id_gtf$occurrence <- trans_merged_quanti$occurrence[match(my_id_gtf$transcript_id, trans_merged_quanti$new_id)]
my_id_cds$occurrence <- trans_merged_quanti$occurrence[match(my_id_cds$transcript_id, trans_merged_quanti$new_id)]


my_id_gtf$orf_id <- info$orf_id[match(my_id_gtf$transcript_id, info$new_id)]
my_id_cds$orf_id <- info$orf_id[match(my_id_cds$transcript_id, info$new_id)]





my_id_cds_top15 <- my_id_cds %>% filter(transcript_id %in% trans_rank$new_id)
my_id_gtf_top15 <- my_id_gtf %>% filter(transcript_id %in% trans_rank$new_id)




get_locus <- function(seqnames, start, end, strand, gene_id) {
  
  locus_subset <- 
    GRanges(seqnames = seqnames,
            ranges = IRanges(start = start, 
                             end = end),
            strand = strand)
}



#zoomed in for top 15
gene_locus <- get_locus(gene_chrom, as.numeric(gene_start),as.numeric(gene_end),gene_strand, gene_id )

structure_filtered_top <- my_id_gtf_top15 %>% 
  ggplot(aes(
    xstart = start,
    xend = end,
    y = reorder(transcript_id, occurrence),
    
  )) +
  geom_range(
    # fill = "white",
    height = 0.25
  ) +
  geom_range(
    data = my_id_cds_top15,
    aes(fill = orf_id)
  ) +
  geom_intron(
    data = to_intron(my_id_gtf_top15, "transcript_id"),
    aes(strand = strand),
    arrow.min.intron.length = 500,
  ) + 
  theme_bw()+ 
  xlab("Coordinates on CHR") + 
  ylab(NULL) + 
  scale_fill_manual(values = color_used)+
  guides(fill = guide_legend(ncol = 1))+
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))+ 
  xlim(start(gene_locus), 
       end(gene_locus))





filtered_prop <- info %>% 
  filter(new_id %in% trans_rank$new_id) %>%
  dplyr::select(new_id, colnames(trans_merged_filtered)[grep('rep', colnames(trans_merged_filtered))])


prop_filtered_long <- reshape2::melt(filtered_prop)
prop_filtered_long <- prop_filtered_long %>% dplyr::filter(variable != 'occurrence')
prop_filtered_long$occurrence <- transcript_info$occurrence[match(prop_filtered_long$new_id, transcript_info$new_id)]



prop_filtered_long[prop_filtered_long == 0] <- NA
prop_filtered_long <- prop_filtered_long[complete.cases(prop_filtered_long),]

orf_long_top <- prop_filtered_long %>% filter(prop_filtered_long$new_id %in% trans_rank$new_id)


orf_long_top$orf_id <- info$orf_id[match(orf_long_top$new_id, info$new_id)]

usage_top <- ggplot(orf_long_top, aes(reorder(x=new_id, occurrence), y= value, fill = orf_id)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.25) +
  scale_y_continuous(n.breaks = 9, labels = scales::percent, limits = c(0,1) )+
  ylab("Proportion") + 
  xlab("Transcripts") +
  theme(plot.title = element_text(size=11),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 1),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  coord_flip() + 
  scale_fill_manual(values = color_used) + 
  guides(fill = guide_legend(ncol = 1))



occur_top <- info %>% 
  filter(new_id %in% trans_rank$new_id)

#occur_top$occurrence <- transcript_info$occurrence[match(occur_top$new_id, transcript_info$new_id)]


occurrence_top <-  ggplot(occur_top,
                          aes(x=reorder(new_id, occurrence), 
                              y=occurrence, fill = transcript_class)) +
  geom_col( position = "identity") + 
  theme(axis.text.x = element_text(angle = 90))+ coord_flip() + 
  xlab("Transcript ID") + 
  ylab('Occurrence') + 
  scale_y_continuous(breaks = c(0,20,40,60,80,sample_size)) + 
  theme_bw() + 
  scale_fill_manual(values = color_types) + 
  theme(axis.title.y = element_text(size = rel(1.5), face = "bold"),
        axis.title.x = element_text(size = rel(1), face = "bold"),
        axis.text.x = element_text(size = 11))



fig2_top <- cowplot::plot_grid(occurrence_top + 
                                 theme(legend.position = "left"), 
                               structure_filtered_top + 
                                 theme(axis.text.y = element_blank(),
                                       axis.ticks.y = element_blank(),
                                       axis.title.y = element_blank(), 
                                       legend.position = "none"), 
                               usage_top + 
                                 theme(axis.title.y = element_blank(),                                                    
                                       axis.text.y = element_blank(),
                                       legend.justification = c(0,0.5),
                                       legend.position = "right"),
                               nrow=1,
                               labels=NA, 
                               rel_widths = c(1.2, 2.5, 1.2), 
                               align = "h")



#########output files

                              
#change gene specific name for output files
#change output file names in the snakemake rule, or no output
stage <- 'generating output files'
                                  
ggsave(paste("Results/top_transcripts_sum_", gene_name, ".png", sep = ""), plot = fig2_top ,width = 30, height = 12)
ggsave(paste("Results/dou_", gene_name, ".png", sep = ""), plot = dou_plot ,width = 7, height = 4.5)
ggsave(paste("Results/orf_usage_",gene_name, ".png", sep = ""), plot = orf_usage_full ,width = 10, height = 20)


write.csv(transcript_info, file = paste("Results/transcript_usage_",gene_name,".csv", sep = ""), row.names = F)
write.csv(info, file = paste("Results/transcript_summary_",gene_name,".csv", sep = ""), row.names = F)
write.csv(sample_info, file = paste("Results/gene_quanti_",gene_name,".csv", sep = ""), row.names = F)
write.csv(orf_prop_sum, file = paste("Results/orf_usage_sum_", gene_name, ".csv", sep = ""), row.names = F)

                                  
                                  .log_msg("DONE  ", processing_gene)

  }, error = function(e) {

    err_msg <- conditionMessage(e)
    .log_msg("FAIL  ", processing_gene, " | stage=", stage, " | ", err_msg)

    .record_failure(processing_gene, stage, err_msg)
    .save_debug_bundle(processing_gene, stage, err_msg)
    .write_placeholder_outputs(processing_gene, stage, err_msg)
    
if (grepl("each range must have an end", err_msg, fixed = TRUE) && exists("my_id_whole_gtf", inherits = FALSE)) {
  bad <- my_id_whole_gtf %>% dplyr::filter(is.na(start) | is.na(end) | start > end)
  try(saveRDS(bad, file.path("Results/debug_failed_genes", processing_gene, "bad_ranges.rds")), silent = TRUE)
}
    # continue to next gene
    NULL
  })
}


## ---- write failed gene summary (ADD) ----
failed_df <- if (length(.failed_rows) == 0) {
  data.frame(gene=character(), stage=character(), error_message=character(), time=character(),
             stringsAsFactors = FALSE)
} else {
  do.call(rbind, .failed_rows)
}

write.table(failed_df, file = "Results/failed_genes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(failed_df$gene, con = "Results/failed_genes.txt")

.log_msg("FAILED_GENES_COUNT ", nrow(failed_df))                                  

