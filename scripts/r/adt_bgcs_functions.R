# function for comparing diversity metrics across different sites

compare_diversity <- function(metric) {
  metric_val <- deparse(substitute(metric))
  metric_var <- enquo(metric)
  test <- as.formula(paste0(metric_val,"~site"))
  kw <- kruskal.test(test, data = shannon)
  dt <- dunnTest(test, data = shannon, method="bh")
  cld <- cldList(P.adj ~ Comparison,
                 data = dt$res,
                 threshold = 0.05)
  cld$site <- unique(sort(shannon$site))
  max <- max(shannon[[metric_val]])
  max_groups <- shannon %>%
    group_by(site) %>%
    summarize(!! metric_var := max(!!metric_var) + 0.05 * max)
  max_groups <- left_join(max_groups, cld)
  e2 <- epsilonSquared(x = shannon[[metric_val]], g=shannon$site)
  return_list <- list("stats" = max_groups, "e2" = e2)
}

# functions for using DESeq to identify site-specific enriched BGCs

## threshold cutoffs for differential abundance analysis
pvalue_cutoff <- 1e-8
fc_cutoff <- 2

## function for comparing each site to all others
dds_by_site <- function(site, countMatrix, metaData) {
  body_site <- deparse(substitute(site))
  design_val <- paste("~", body_site)
  dds_site <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = metaData,
                              design = as.formula(design_val))
  dds_site <- DESeq(dds_site)
  res_site <- results(dds_site, contrast = c(body_site,body_site,"other"))
}

get_enriched <- function(dds){
  # convert dds object to tibble and set rowname as a column called "cluster"
  dds_df <- as.data.frame(dds)
  dds_df <- rownames_to_column(dds_df, var = "cluster") %>% as_tibble()
  dds_df <- dds_df %>%
    filter(padj <= pvalue_cutoff & log2FoldChange >= fc_cutoff) %>%
    select(cluster)
}

## calculate theoretical max y-intercept for volcano plots if P-value = 0
max_intercept <- function(dds){
  padj = "padj"
  p_val_no_na <- (dds[[padj]][which(!is.na(dds[[padj]]))])
  p_val_min <- min(p_val_no_na[which(p_val_no_na > 0)])
  y_int <- -log10(10^-1 * p_val_min)
}

## function for generating volcano plots for each site
volcano_plot_site <- function(res, plot_title) {EnhancedVolcano(res,
                                                                lab = rownames(res),
                                                                x = 'log2FoldChange',
                                                                y = 'padj',
                                                                title = plot_title,
                                                                subtitle = NULL,
                                                                caption = NULL,
                                                                xlim = c(-12,12),
                                                                ylim = c(2e-03,320),
                                                                xlab = bquote(~log[2]~ 'Fold Difference'),
                                                                ylab = bquote(~-log[10]~adjusted~italic(P)),
                                                                pCutoff = pvalue_cutoff,
                                                                FCcutoff = fc_cutoff,
																                                cutoffLineWidth = 0.25,
                                                                labSize = 4,
                                                                colAlpha = 0.25,
                                                                legendPosition = "bottom",
																                                legendLabels = c('NS', bquote(log[2] ~ 'Fold Difference'), bquote(italic(P)), bquote(italic(P)~ and ~ log[2] ~ 'Fold Difference')),
                                                                selectLab = c("clusterx"), # ensures that no clusters are labeled, which would otherwise clutter the plot
)
}

# function to extract median abundance counts for each cluster in specific enriched site from function get_enriched

cluster_med_counts <- function(site, enriched_object) {
  med_counts <- hmp_adt_kallisto_counts %>%
    filter(cluster %in% enriched_object$cluster & site == !!site) %>%
    group_by(cluster) %>%
    summarize(cluster_med = median(log_normalized_counts))
}

# function for Jaccard similarity
## adapted from: jsb at https://stats.stackexchange.com/questions/176613/jaccard-similarity-in-r

jaccard <- function(x, m) {
  if (m == 1 | m == 2) {
    M_00 <- apply(x, m, sum) == 0
    M_11 <- apply(x, m, sum) == 2
    if (m == 1) {
      x <- x[!M_00, ]
      JSim <- sum(M_11) / nrow(x)
    } else {
      x <- x[, !M_00]
      JSim <- sum(M_11) / length(x)
    }
    JDist <- 1 - JSim
    return(c(JSim = JSim, JDist = JDist))
  } else break
}
