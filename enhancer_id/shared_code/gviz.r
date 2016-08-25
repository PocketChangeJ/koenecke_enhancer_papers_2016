library(Gviz)

gviz_data_track <- function(bw, title, track_color, max_value) {
    
  reg_track <- DataTrack(range=bw,
                         genome="dm3", name=title, 
                         type="polygon",
                         ylim=c(0, max_value),
                         col=track_color,
                         fill.mountain=c(track_color, track_color))
  reg_track
}

max_value <- function(bw, region.gr) {
  stopifnot(length(bw) == 1)
  stopifnot(length(region.gr) == 1)
  regionMaxs(region.gr, bw)
}

gviz_dual_data_track <- function(m.bw, de.bw, region.gr, title, colors, max_value) {
  
  colors %<>% rev
  colors <- c(colors, "black")
  
  data_matrix <- c(unlist(import(m.bw,  which=region.gr, as="NumericList")),
                   unlist(import(de.bw, which=region.gr, as="NumericList")),
                   rep(0, times=width(region.gr))) %>%
                 matrix(nrow=3, byrow=TRUE)
  
  data_matrix[, c(1, ncol(data_matrix))] <- 0
  
  reg_track <- DataTrack(data=data_matrix,
                         start=start(region.gr):end(region.gr), width=0, #ncol(data_matrix), 
                         chromosome=as.character(seqnames(region.gr)),
                         genome="dm3", name=title, 
                         groups=c("Toll10b", "gd7", "zero"),
                         type="a",#alpha=0.5,
                         ylim=c(0, max_value),
                         legend=TRUE,
                         col.mountain=colors, fill.mountain=colors, col.line=colors, col=colors)
  reg_track
}
