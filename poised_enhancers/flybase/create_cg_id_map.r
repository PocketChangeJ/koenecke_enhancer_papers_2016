
ids <- read.delim("fbgn_annotation_ID_fb_2014_03.tsv.gz", stringsAsFactors=F, header=F, skip=5)
names(ids) <- c("fb_symbol", "fb_gene_id", "prev_id", "fb_cg_id", "prev_cg_id")
ids$fb_symbol <- NULL
ids$prev_id <- NULL

ids.expand <- subset(ids, prev_cg_id != "")
ids.expand <- transform(ids.expand[rep(seq(nrow(ids.expand)), sapply(v.s <- strsplit(ids.expand$prev_cg_id, split=","), length)),], prev_cg_id=unlist(v.s))

ids.same <- unique(data.frame(stringsAsFactors=F, 
                              fb_gene_id = as.character(ids$fb_gene_id), 
                              fb_cg_id   = as.character(ids$fb_cg_id),
                              prev_cg_id = as.character(ids$fb_cg_id)))

ids.all <- rbind(ids.same, ids.expand)

fbidmap <- ids.all
saveRDS(fbidmap, file="fbcgidmap.rds")
