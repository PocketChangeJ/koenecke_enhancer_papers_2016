library(GenomicRanges)
library(parallel)
library(dplyr)
library(compiler)

load_motif <- function(motif_name) {
  updateObject(readRDS(file.path(motif_granges_path(), paste0(motif_name, ".granges.rds"))))
}

match_motifs_to_fb_genes <- function(motif_files) {
  m2g <- data_frame(motif=gsub(".granges.rds", "", motif_files))

  m2g$motif_fb_id <- gsub(".*_(FBgn.......)$", "\\1", m2g$motif)
  m2g$motif_fb_id[m2g$motif == m2g$motif_fb_id] <- NA

  fbidmap <- readRDS("flybase/fbidmap.rds")
  fb.genes <- unique(flybase_txs()[, c("fb_gene_id", "fb_symbol")])

  m2g$matched_id <- fbidmap$fb_gene_id[match(m2g$motif_fb_id, fbidmap$prev_id)]

  jaspar_i <- grep("^JASPAR", m2g$motif)
  if(length(jaspar_i) > 0) {
    jaspar_symbol <- gsub("^JASPAR_....-(.*)-M.*$", "\\1", m2g$motif[jaspar_i])
    jaspar_symbol <- gsub("^(.*)_.*$", "\\1", jaspar_symbol)
    jaspar_matched_id <- fb.genes$fb_gene_id[match(jaspar_symbol, fb.genes$fb_symbol)]
    jaspar_matched_id <- ifelse(is.na(jaspar_matched_id), jaspar_symbol, jaspar_matched_id)
    m2g$matched_id[jaspar_i] <- jaspar_matched_id
  }

  m2g$motif_fb_id <- m2g$matched_id
  m2g$matched_id <- NULL
  m2g$motif_symbol <- fb.genes$fb_symbol[match(m2g$motif_fb_id, fb.genes$fb_gene_id)]

  m2g
}

single_motif_counts <- function(motif_gr_rds, peaks.gr) {
  stopifnot(file.exists(motif_gr_rds))
  motifs.gr <- updateObject(readRDS(motif_gr_rds))

  motif_name <- motifs.gr$motif[1]

  ol <- findOverlaps(query=motifs.gr, subject=peaks.gr, type="within", ignore.strand=TRUE)
  with_motif <- length(unique(subjectHits(ol)))
  without_motif <- length(peaks.gr) - with_motif
  stopifnot(with_motif + without_motif == length(peaks.gr))
  data.frame(stringsAsFactors=FALSE, 
             motif         = motif_name, 
             with_motif    = with_motif,
             without_motif = without_motif)
}
single_motif_counts %<>% cmpfun

motif_counts_rds <- function(peaks.gr, motif_files, cores=3) {
  message("Counting ", pn(length(motif_files)), " motifs in ", pn(length(peaks.gr)), " peaks")
  counts.df <- motif_files %>%
               mclapply(single_motif_counts, peaks.gr, mc.cores=cores, mc.preschedule=FALSE) %>%
               bind_rows
  counts.df
}

proportion_test <- function(values, enrichment=TRUE) {
  stopifnot(is.logical(enrichment))
  m <- matrix(values, nrow=2, byrow=TRUE)
  alt.test <- ifelse(enrichment, "greater", "less")

  data_frame(test_type = ifelse(enrichment, "enrichment", "depletion"),
             pvalue    = prop.test(m, alternative=alt.test)$p.value)
}
proportion_test %<>% cmpfun

single_motif_test <- function(i, set.df) {
  expected_proportion <- set.df$s2_W[i] / (set.df$s2_W[i] + set.df$s2_WO[i])
  actual_proportion   <- set.df$s1_W[i] / (set.df$s1_W[i] + set.df$s1_WO[i])
  
  enrichment <- actual_proportion >= expected_proportion
  
  proportion_test(as.numeric(set.df[i, c("s1_W", "s1_WO", "s2_W", "s2_WO")]), enrichment)
}
single_motif_test %<>% cmpfun

# tests enrichment or depletion
motif_count_comparison <- function(set1.df, set2.df) {
  stopifnot(nrow(set1.df) == nrow(set2.df))
  
  names(set1.df)[2:3] <- c("s1_W", "s1_WO")
  names(set2.df)[2:3] <- c("s2_W", "s2_WO")
  
  set.df <- merge(set1.df, set2.df)
  stopifnot(nrow(set.df) == nrow(set1.df))
  
  set_test.df   <- subset(set.df, s1_W > 0 | s2_W > 0)
  
  test_results.df <- bind_rows(lapply(1:nrow(set_test.df), single_motif_test, set_test.df))

  set_test.df$test_type <- test_results.df$test_type
  set_test.df$pvalue    <- test_results.df$pvalue

  set_notest.df <- subset(set.df, s1_W == 0 & s2_W == 0)
  if(nrow(set_notest.df) > 0) {
    set_notest.df$pvalue <- 1
    set_notest.df$test_type <- "enrichment"
    set.df <- rbind(set_test.df, set_notest.df)
  } else {
    set.df <- set_test.df
  }
  
  stopifnot(nrow(set.df) == nrow(set1.df))
  set.df <- transform(set.df, enrichment = (s1_W/(s1_W+s1_WO)) / (s2_W/(s2_W+s2_WO)) )
  set.df <- transform(set.df, enrichment = ifelse(test_type == "enrichment", enrichment, -1 / enrichment))
  set.df
}

clustered_motif_heatmap <- function(df, title, motif.order=c(), group.order=c(), enrichment.max=10) {
  df <- as.data.frame(df)
  e.df <- reshape(df[, c("motif", "group_label", "enrichment")],
                 idvar="motif", timevar="group_label", v.names="enrichment",
                 direction="wide")
  e.m <- as.matrix(e.df[, -1])
  rownames(e.m) <- e.df$motif
  colnames(e.m) <- gsub("enrichment.", "", colnames(e.m))

  motifs.d <- as.dendrogram(hclust(dist(e.m)))
  factors.d <- as.dendrogram(hclust(dist(t(e.m))))
  
  if(length(motif.order) == 0) {
    motif.order <- rownames(e.m)[order.dendrogram(motifs.d)]
  }
  
  if(length(group.order) == 0) {
    group.order <- rownames(t(e.m))[order.dendrogram(factors.d)]
  }
  
  message("Motif order: ", paste0(motif.order, collapse=", "))
  message("Group order: ", paste0(group.order, collapse=", "))
  
  df$motif <- factor(df$motif, levels=motif.order)
  df$group_label <- factor(df$group_label, levels=group.order)
  
  stars <- pmin(5, -1 * ceiling(log10(sig.df$adj_pv)))
  
  df$sig_label <- stars %>%
                  lapply(function(count) {
                    paste0(rep("*", times=count), collapse="")
                  }) %>%
                  unlist
  
  df$sig_label[df$adj_pv >= pv_cutoff()] <- ""
  
  e.max <- pmin(enrichment.max, ceiling(max(abs(df$enrichment))))
  df$enrichment <- pmin(e.max, df$enrichment)
  df$enrichment <- pmax(-e.max, df$enrichment)
  
  g <- ggplot(df, aes(x=group_label, y=motif, fill=enrichment)) +
       geom_tile(color=NA) +
       geom_text(aes(label=sig_label), color="gray10", size=3) +
       scale_y_discrete(expand=c(0, 0)) +
       scale_x_discrete(expand=c(0, 0)) +
       scale_fill_gradientn(name="Enrichment", space="Lab", 
                            values=c(-e.max, -1, 1, e.max), 
                            colours=c("darkblue", "white", "white", "darkred"), 
                            rescaler=function(x,...) x, oob=identity,
                            limits=c(-e.max, e.max),
                            guide=guide_colorbar()) +
       theme_manuscript(base_size=10) +
       labs(x="", y="", title=title) +
       theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
             axis.text.y=element_text(hjust=1),
             panel.grid.minor=element_blank(),
             panel.grid.major=element_blank(),
             axis.line=element_blank(),
             axis.ticks=element_line(size = 0.5, color="gray50"))
  g
}


