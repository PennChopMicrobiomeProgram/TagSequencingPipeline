library(magrittr)
library(ape)
library(dplyr)
library(pheatmap)

###================================================
###
### Some useful R functions are collected here
###
###================================================

###=====
###  adapted from simplify_assignments()
###  long_assignments <- function(adf) 
###  adf: taxonomic assignment table with 7 columns
###=====

long_assignments <- function(adf) {
  apply(adf, 1, function(x) {
    x <- na.omit(as.character(x))
    paste(x, collapse=" ")
  })
}

###=====
###  thanks to Ceylan
###  make_pcoa_plot <- function(uu, s, shape_by, color_by, title)
###  uu: distance, s: mapping file, shape_by: variable used for shape, color_by: variable used for color, title: title
###=====

make_pcoa_plot <- function(uu, s, shape_by, color_by, title) {
  uu_pcoa <- pcoa(uu)
  uu_df <- merge(s, uu_pcoa$vectors[, 1:5], by.x="SampleIDorig", by.y="row.names")
  uu_pct <- round(uu_pcoa$values$Relative_eig * 100)
  
  g_uu = ggplot(uu_df, aes(x=Axis.1, y=Axis.2)) +
    #coord_equal() +
    theme_bw() +
    xlab(paste0("PCoA axis 1 (", uu_pct[1], "%)")) +
    ylab(paste0("PCoA axis 2 (", uu_pct[2], "%)")) +
    scale_shape_discrete(name=sub("_", " ", shape_by)) + 
    scale_colour_discrete(name=sub("_", " ", color_by)) +
    ggtitle(title)
  
  if (is.null(shape_by) & !is.null(color_by)) {
    g_uu <- g_uu + geom_point(aes(colour=factor(get(color_by))))
  } else if (!is.null(shape_by) & !is.null(color_by)) {
    g_uu <- g_uu + geom_point(aes(colour=factor(get(color_by)), shape=factor(get(shape_by))))
  } else {
    g_uu <- g_uu + geom_point()
  }
  return(g_uu)
}

###=====
###  thanks to Ceylan
###  filter_low_coverage <- function(props, frac_cutoff=0.6) 
###  get taxa that are found in at least "frac_cutoff" proportion of the samples 
###=====
 
filter_low_coverage <- function(props, frac_cutoff=0.6){
  frac_nonzero <- function (x) sum(x > 0) / length(x)
  apply(props, 1, frac_nonzero) >= frac_cutoff
}

###=====
###  thanks to Chunyu
###  heatmap_grouped <- function(genus_props, heatmap_s, grps = c("study_group", "study_day"), fname=NULL, thre=0.8, option=1)
###  option=1: rows_to_keep <- filter_low_coverage(heatmap_props, perc_cutoff=thre) ## taxa found in at least 80% of samples
###  option=2: rows_to_keep <- apply(heatmap_props,1,max) >= 0.01 ## taxa with abundance in any sample exceeding 1%
###=====

heatmap_grouped <- function(summed_props, heatmap_s, grps = c("study_group", "study_day"), fname=NULL, thre=0.8, option=1, prop_cut=0.01, satu_limit=0.4){
  
  #color = saturated_rainbow(101)
  color = saturated_rainbow(101, saturation_limit=satu_limit)
  breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))
  
  heatmap_props <- summed_props[,heatmap_s$SampleID]
  
  if (option == 1) {
    rows_to_keep <- filter_low_coverage(heatmap_props, frac_cutoff=thre) 
  } else if (option == 2) {
    rows_to_keep <- apply(heatmap_props,1,max) >= prop_cut 
  }
  heatmap_props <- heatmap_props[rows_to_keep,]
  
  ## group the SampleIDs
  heatmap_s %<>% arrange_(.dots=grps)
  heatmap_props <- heatmap_props[, heatmap_s$SampleID]
  
  ## update the annotation
  annc <- heatmap_s[,grps] %>% as.data.frame()
  rownames(annc) <- heatmap_s$SampleID
  colnames(annc) <- grps
  
  ## heatmap time
  if (!is.null(fname))
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, filename = fname, 
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
  else
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, 
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
}

###=====
###  adapted from Chunyu
###  heatmap_row_gap <- function(genus_props, heatmap_s, BSH_table = b, grps = c("study_group", "study_day"), fname=NULL, thre=0.8, option=1)
###  BSH_table: Bile Salt Hydrolase genes information with two columns: "taxa_id", "HasBSH" 
###  option=1: rows_to_keep <- filter_low_coverage(heatmap_props, perc_cutoff=thre) ## taxa found in at least 80% of samples
###  option=2: rows_to_keep <- apply(heatmap_props,1,max) >= 0.01 ## taxa with abundance in any sample exceeding 1%
###=====

heatmap_row_gap <- function(genus_props, heatmap_s, BSH_table = b, grps = c("study_group", "study_day"), fname=NULL, thre=0.8, option=1){
  
  color = saturated_rainbow(101)
  breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))
  
  heatmap_props <- genus_props[,heatmap_s$SampleID]
  
  if (option == 1) {
    rows_to_keep <- filter_low_coverage(heatmap_props, perc_cutoff=thre) 
  } else {
    rows_to_keep <- apply(heatmap_props,1,max) >= 0.01 
  }
  
  heatmap_props <- heatmap_props[rows_to_keep,] 
  all_heatmap_taxa <- rownames(heatmap_props) ## list of taxa in the heatmap
  
  tt <- names(a)[a %in% all_heatmap_taxa] ## match this with taxa_id 
  bs <- b %>% filter(taxa_id %in% tt)
  aa <- as.data.frame(a[names(a) %in% bs$taxa_id]) 
  aa$taxa_id <- rownames(aa)
  rownames(aa) <- NULL
  colnames(aa) <- c("organism", "taxa_id")
  m <- merge(aa, b, by="taxa_id", all.x=T) %>%
    select(organism, HasBSH)
  m <- m[!duplicated(m),]
  
  ### determine the order of taxa in the heatmap
  taxa_order <- data.frame(organism=all_heatmap_taxa) %>% merge(m, by="organism", all.x=T) 
  taxa_order$order[!is.na(taxa_order$HasBSH) & taxa_order$HasBSH==TRUE] <- 1
  taxa_order$order[!is.na(taxa_order$HasBSH) & taxa_order$HasBSH==FALSE] <- 2
  taxa_order$order[is.na(taxa_order$HasBSH)] <- 3
  taxa_order %<>% arrange(order, organism)
  
  ### position of the row gaps
  gaps_row <- c(max(which(taxa_order$order==1)), max(which(taxa_order$order==2)))
  
  ## group the SampleIDs and sort rows by taxa order
  heatmap_s %<>% arrange_(.dots=grps)
  heatmap_props <- heatmap_props[taxa_order$organism, heatmap_s$SampleID]
  
  ## update the annotation
  annc <- heatmap_s[,grps] %>% as.data.frame()
  rownames(annc) <- heatmap_s$SampleID
  colnames(annc) <- grps
  
  ## heatmap time
  if (!is.null(fname))
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, filename = fname, 
             fontsize_col = 8, fontsize_row = 8, gaps_row = gaps_row, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
  else
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, 
             fontsize_col = 8, fontsize_row = 8, gaps_row = gaps_row, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
}

###=====
###  adapted from Chunyu
###  tidy_lm <- function(summed_cts, grp="Age")
###  BSH_table: Bile Salt Hydrolase genes information with two columns: "taxa_id", "HasBSH" 
###=====

tidy_lm <- function(genus_cts_top_df, grp="Age") {
  #genus_cts_top_df <- genus_props_top_df
  formula <- paste0("LogProp ~ ", grp)
  expr <- lazyeval::interp(quote(! is.na(x)), x = as.name(grp)) 
  genus_cts_top_df %<>%  filter_(expr) %>% droplevels()
  
  lm_models <- genus_cts_top_df %>%
    group_by(Taxon) %>%
    do(mod = lm(as.formula(formula), data=.)) %>%
    ungroup()
  
  #lapply(lm_models$mod, summary)
  
  summaries <- lapply(1:length(lm_models$mod), function(x) data.frame(tidy(lm_models$mod[[x]]), taxa=lm_models$Taxon[[x]]))
  
  summaries_df <- do.call(rbind, summaries) %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value <= sig.pval.cut) %>%
    dplyr::select(taxa, everything())
  
  if (dim(summaries_df)[1] > 0){
    taxa <- summaries_df$taxa %>% as.character()
    nrow = ceiling(length(taxa) / 4)
    
    fig <- genus_cts_top_df %>%
      filter(Taxon %in% taxa) %>%
      droplevels() %>%
      ggplot(aes(x=get(grp), y=Proportion, color=get(grp))) + 
      geom_boxplot(coef = 10) +
      geom_quasirandom(group = 1) +
      scale_color_brewer(palette = "Set2", name=grp) +
      scale_y_log10() + 
      facet_wrap(~ TaxonLabel, scales = "free_y", nrow = nrow) +
      labs(x = "", y = "Relative Abundance") + 
      theme_bw() +
      ggtitle(grp) +
      theme(plot.title = element_text(hjust = 0.5))
    
    plot(fig)
  }
  pander(as.data.frame(summaries_df)) ### without as.data.frame, somehow pander does not work 
}

###=====
###  Thanks to Chunyu
###  a function used to read the kegg assignment file
###  tidy_lm <- function(summed_cts, grp="Age")
###  BSH_table: Bile Salt Hydrolase genes information with two columns: "taxa_id", "HasBSH" 
###=====

read_ko_assign <- function(filePath){
  
  ko <- read.table(filePath, sep='\t', header=TRUE, stringsAsFactors = FALSE, quote="", comment='',row.names = 1)
  ko <- sweep(ko, 2, colSums(ko), "/")
  colnames(ko) %<>%
    sub("^PCMP_", "", ., perl=T)
  
  return(ko)
}

