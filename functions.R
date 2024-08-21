plot_rich_reads_samlenames_lm <- function(physeq, group = "Site", label = "Repeat"){
  rish <- estimate_richness(physeq, measures = "Observed")
  reads.sum <- as.data.frame(sample_sums(physeq))
  reads.summary <- cbind(rish, reads.sum)
  colnames(reads.summary) <- c("otus","reads")
  reads.summary["Repeat"] <-unlist(purrr::map(stringr::str_split(rownames(physeq@sam_data), "\\_", 2), function(x) x[[2]]))
  reads.summary["Site"] <- physeq@sam_data[[group]]
  library(ggrepel)
  require(ggforce)
  p1 <- ggplot(data=reads.summary) + 
    geom_point(aes(y=otus, x=log2(reads), color=Site),size=3) + 
    geom_text_repel(aes(y=otus, x=log2(reads), label=paste0(Repeat))) + 
    theme_bw() +
    geom_smooth(aes(y=otus, x=log2(reads), fill=Site, color=Site),method=lm, se=FALSE, ymin = 1) + 
    scale_x_continuous(sec.axis = sec_axis(sec.axis ~ 2**.)) 

  return(p1)
}


beta_custom_norm_NMDS_elli_w <- function(ps, seed = 7888, normtype="vst", Color="What", Group="Repeat"){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  library(ggforce)
  
  unicode_minus <- function(x) sub('^-', '\U2212', format(x))
  ps@otu_table[ps@otu_table < 0] <- 0
  ordination.b <- ordinate(ps, "NMDS", "bray")
  mds <- as.data.frame(ordination.b$points)
  p  <-  plot_ordination(ps,
                         ordination.b,
                         type="sample",
                         color = Color,
                         title="NMDS - Bray-Curtis",
                         axes = c(1,2) ) + 
    theme_bw() + 
    theme(text = element_text(size = 10)) + 
    geom_point(size = 3) +
    annotate("text",
    x=min(mds$MDS1) + abs(min(mds$MDS1))/4,
    y=max(mds$MDS2),
    label=paste0("Stress -- ", round(ordination.b$stress, 3))) +
    geom_mark_ellipse(aes_string(group = Group, label = Group),
                      label.fontsize = 10,
                      label.buffer = unit(2, "mm"),
                      label.minwidth = unit(5, "mm"),
                      con.cap = unit(0.1, "mm"),
                      con.colour='gray') +
    scale_y_continuous(labels = unicode_minus) +
    scale_x_continuous(labels = unicode_minus) +
    theme(legend.position = "none") +
    scale_colour_viridis_d(option = "magma", 
                       aesthetics = "color", 
                       begin = 0, 
                       end = 0.8)

  
  return(p)
}

plot_alpha_w_toc_2 <- function(ps, group, metric) {
  
  require(phyloseq)
  require(ggplot2)
  
  ps_a <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  er <- estimate_richness(ps_a)
  df_er <- cbind(ps_a@sam_data, er)
  df_er <- df_er %>% dplyr::select(c(group, metric))
  stat.test <- aov(as.formula(paste0(metric, "~", group)), data = df_er) %>%
    rstatix::tukey_hsd()
  y <- seq(max(er[[metric]]), length=length(stat.test$p.adj.signif[stat.test$p.adj.signif != "ns"]), by=max(er[[metric]]/20))
  
  plot_richness(ps_a, x=group, measures=metric) + 
    geom_boxplot() +
    geom_point(size=1.2, alpha=0.3) +
    ggpubr::stat_pvalue_manual(
      stat.test, 
      label = "p.adj.signif", 
      y.position = y,
      hide.ns=TRUE) +
    theme_light() + 
    scale_color_brewer(palette="Dark2") +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x=element_blank()) +
    labs(y=paste(metric, "index")) 
}

phyloseq_to_ampvis2 <- function(physeq) {
  #check object for class
  if(!any(class(physeq) %in% "phyloseq"))
    stop("physeq object must be of class \"phyloseq\"", call. = FALSE)
  
  #ampvis2 requires taxonomy and abundance table, phyloseq checks for the latter
  if(is.null(physeq@tax_table))
    stop("No taxonomy found in the phyloseq object and is required for ampvis2", call. = FALSE)
  
  #OTUs must be in rows, not columns
  if(phyloseq::taxa_are_rows(physeq))
    abund <- as.data.frame(phyloseq::otu_table(physeq)@.Data)
  else
    abund <- as.data.frame(t(phyloseq::otu_table(physeq)@.Data))
  
  #tax_table is assumed to have OTUs in rows too
  tax <- phyloseq::tax_table(physeq)@.Data
  
  #merge by rownames (OTUs)
  otutable <- merge(
    abund,
    tax,
    by = 0,
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
  )
  colnames(otutable)[1] <- "OTU"
  
  #extract sample_data (metadata)
  if(!is.null(physeq@sam_data)) {
    metadata <- data.frame(
      phyloseq::sample_data(physeq),
      row.names = phyloseq::sample_names(physeq), 
      stringsAsFactors = FALSE, 
      check.names = FALSE
    )
    
    #check if any columns match exactly with rownames
    #if none matched assume row names are sample identifiers
    samplesCol <- unlist(lapply(metadata, function(x) {
      identical(x, rownames(metadata))}))
    
    if(any(samplesCol)) {
      #error if a column matched and it's not the first
      if(!samplesCol[[1]])
        stop("Sample ID's must be in the first column in the sample metadata, please reorder", call. = FALSE)
    } else {
      #assume rownames are sample identifiers, merge at the end with name "SampleID"
      if(any(colnames(metadata) %in% "SampleID"))
        stop("A column in the sample metadata is already named \"SampleID\" but does not seem to contain sample ID's", call. = FALSE)
      metadata$SampleID <- rownames(metadata)
      
      #reorder columns so SampleID is the first
      metadata <- metadata[, c(which(colnames(metadata) %in% "SampleID"), 1:(ncol(metadata)-1L)), drop = FALSE]
    }
  } else
    metadata <- NULL
  
  #extract phylogenetic tree, assumed to be of class "phylo"
  if(!is.null(physeq@phy_tree)) {
    tree <- phyloseq::phy_tree(physeq)
  } else
    tree <- NULL
  
  #extract OTU DNA sequences, assumed to be of class "XStringSet"
  if(!is.null(physeq@refseq)) {
    #convert XStringSet to DNAbin using a temporary file (easiest)
    fastaTempFile <- tempfile(pattern = "ampvis2_", fileext = ".fa")
    Biostrings::writeXStringSet(physeq@refseq, filepath = fastaTempFile)
  } else
    fastaTempFile <- NULL
  
  #load as normally with amp_load
  ampvis2::amp_load(
    otutable = otutable,
    metadata = metadata,
    tree = tree,
    fasta = fastaTempFile
  )
}

norm_anc <- function(ps, anc){
  samp_frac <-  anc$samp_frac
  samp_frac[is.na(samp_frac)] <-  0 
  log_obs_abn <-  log(anc$feature_table + 1)
  log_corr_abn <-  t(t(log_obs_abn) - samp_frac)
  ps.out <- ps
  otu_table(ps.out) <- otu_table(t(log_corr_abn), taxa_are_rows = FALSE)
  return(ps.out)
}
