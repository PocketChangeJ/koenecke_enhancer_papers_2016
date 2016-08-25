library(plyr)

message("Reading in RNA lines...")
tx <- read.delim(pipe("gzcat dmel-all-no-analysis-r5.57.gff.gz | awk '{ if ($3 == \"mRNA\" || $3 == \"pre_miRNA\" || $3 == \"ncRNA\" || $3 == \"rRNA\" || $3 == \"snRNA\" || $3 == \"snoRNA\" || $3 == \"tRNA\") print }'"), stringsAsFactors=F, header=F)

tx$transcript_id <- gsub("^ID=(.*?);.*$",        "\\1", tx$V9)
tx$gene_id       <- gsub("^.*;Parent=(.*?);.*$", "\\1", tx$V9)

fb.transcripts <- data.frame(stringsAsFactors=F, fb_gene_id = tx$gene_id, 
                                                 fb_tx_id   = tx$transcript_id, 
                                                 chr        = paste("chr", tx$V1, sep=""),
                                                 start      = tx$V4,
                                                 end        = tx$V5,
                                                 strand     = ifelse(tx$V7 == "+", 1, -1),
                                                 type       = tx$V3)

# select one transcript among those with the same coordinates
selected.tx <- ddply(fb.transcripts, .var=.(fb_gene_id, chr, start, end, strand, type), summarize, selected_transcript_id=fb_tx_id[1], .progress="text")

fb.transcripts <- subset(fb.transcripts, fb_tx_id %in% selected.tx$selected_transcript_id)

message("Reading annotation ID map file...")
fbnames <- read.delim("fbgn_annotation_ID_fb_2014_03.tsv.gz", stringsAsFactors=F, header=F, skip=5)[, c(1,2,4)]
names(fbnames) <- c("fb_symbol", "fb_gene_id", "fb_cg_id")

fb.transcripts.df <- merge(fbnames, fb.transcripts, all.y=T)

saveRDS(fb.transcripts.df, file="fb.transcripts.r5.57.rds")
