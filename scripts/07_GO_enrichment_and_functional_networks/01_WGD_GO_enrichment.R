#!/usr/bin/env Rscript

# Required R packages:
#   tidyverse
#   clusterProfiler
#   enrichplot

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(enrichplot)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8) {

  stop(
    paste(
      "\nUsage:\n",
      "Rscript 01_WGD_GO_enrichment.R",
      "<sample>",
      "<ks_genelevel.tsv>",
      "<MCScanX.collinearity>",
      "<InterProScan.tsv>",
      "<WGD1_min> <WGD1_max>",
      "<WGD2_min> <WGD2_max>\n"
    )
  )
}


sample_id <- args[1]

ks_file  <- args[2]
col_file <- args[3]
ipr_file <- args[4]

wgd1_min <- as.numeric(args[5])
wgd1_max <- as.numeric(args[6])

wgd2_min <- as.numeric(args[7])
wgd2_max <- as.numeric(args[8])


outdir <- paste0(sample_id, "_WGD_functional_analysis")

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)


normalize_id <- function(x) {

  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^>", "", x)
  x <- sub("\\s.*$", "", x)
  x <- gsub("X1$", ".1", x)

  x
}


ks <- read.table(
  ks_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


if (ncol(ks) < 3) {
  stop("Ks file must contain at least three columns: gene1, gene2, Ks")
}


colnames(ks)[1:3] <- c(
  "gene1",
  "gene2",
  "Ks"
)


ks <- ks %>%

  mutate(
    gene1 = normalize_id(gene1),
    gene2 = normalize_id(gene2),
    Ks = as.numeric(Ks)
  ) %>%

  filter(
    is.finite(Ks),
    Ks > 0,
    Ks < 5
  )


cat(
  "[INFO] Gene pairs with valid Ks:",
  nrow(ks),
  "\n"
)

pdf(
  file.path(
    outdir,
    paste0(
      sample_id,
      "_syntenic_gene_pair_Ks_distribution.pdf"
    )
  ),
  width = 7,
  height = 5
)


hist(
  ks$Ks,
  breaks = 100,
  main = "Ks distribution of syntenic gene pairs",
  xlab = "Ks"
)


dev.off()


lines <- readLines(col_file)


block_df <- data.frame(
  block = character(),
  gene1 = character(),
  gene2 = character(),
  stringsAsFactors = FALSE
)


for (l in lines) {

  if (grepl("^[[:space:]]*[0-9]+-", l)) {

    m <- regmatches(

      l,

      regexec(
        "^[[:space:]]*([0-9]+)-\\s+[0-9]+:\\s+(\\S+)\\s+(\\S+)",
        l
      )

    )[[1]]


    if (length(m) == 4) {

      block_df <- rbind(

        block_df,

        data.frame(
          block = paste0("block_", m[2]),
          gene1 = normalize_id(m[3]),
          gene2 = normalize_id(m[4]),
          stringsAsFactors = FALSE
        )
      )
    }
  }
}


cat(
  "[INFO] Gene pairs extracted from MCScanX:",
  nrow(block_df),
  "\n"
)


make_pair_key <- function(g1, g2) {

  paste(
    pmin(g1, g2),
    pmax(g1, g2),
    sep = "__"
  )
}


ks$key <- make_pair_key(
  ks$gene1,
  ks$gene2
)


block_df$key <- make_pair_key(
  block_df$gene1,
  block_df$gene2
)


block_ks <- merge(
  block_df,
  ks[, c("key", "Ks")],
  by = "key"
)


block_ks_clean <- block_ks %>%

  select(
    block,
    gene1,
    gene2,
    Ks
  ) %>%

  distinct()


cat(
  "[INFO] MCScanX pairs with gene-level Ks:",
  nrow(block_ks_clean),
  "\n"
)


write.table(
  block_ks_clean,
  file.path(
    outdir,
    paste0(
      sample_id,
      "_collinear_gene_pairs_with_Ks.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

block_summary <- block_ks_clean %>%

  group_by(block) %>%

  summarise(

    median_Ks = median(
      Ks,
      na.rm = TRUE
    ),

    n_pairs = n(),

    .groups = "drop"
  )


write.table(
  block_summary,
  file.path(
    outdir,
    paste0(
      sample_id,
      "_syntenic_block_median_Ks.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  "[INFO] Number of syntenic blocks:",
  nrow(block_summary),
  "\n"
)

pdf(
  file.path(
    outdir,
    paste0(
      sample_id,
      "_block_median_Ks_distribution.pdf"
    )
  ),
  width = 7,
  height = 5
)


hist(
  block_summary$median_Ks,
  breaks = 40,
  main = "Block median Ks distribution",
  xlab = "Median Ks"
)


dev.off()


wgd1_blocks <- block_summary %>%

  filter(
    median_Ks >= wgd1_min,
    median_Ks <= wgd1_max
  )


wgd2_blocks <- block_summary %>%

  filter(
    median_Ks >= wgd2_min,
    median_Ks <= wgd2_max
  )


write.table(
  wgd1_blocks,
  file.path(
    outdir,
    paste0(
      sample_id,
      "_WGD1_blocks.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  wgd2_blocks,
  file.path(
    outdir,
    paste0(
      sample_id,
      "_WGD2_blocks.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  "[INFO] WGD1 blocks:",
  nrow(wgd1_blocks),
  "\n"
)

cat(
  "[INFO] WGD2 blocks:",
  nrow(wgd2_blocks),
  "\n"
)

wgd1_genes <- block_ks_clean %>%

  filter(
    block %in% wgd1_blocks$block
  ) %>%

  select(
    gene1,
    gene2
  ) %>%

  unlist(
    use.names = FALSE
  ) %>%

  unique()


wgd2_genes <- block_ks_clean %>%

  filter(
    block %in% wgd2_blocks$block
  ) %>%

  select(
    gene1,
    gene2
  ) %>%

  unlist(
    use.names = FALSE
  ) %>%

  unique()


background_genes <- unique(

  c(
    block_ks_clean$gene1,
    block_ks_clean$gene2
  )

)


wgd1_genes <- normalize_id(
  wgd1_genes
)

wgd2_genes <- normalize_id(
  wgd2_genes
)

background_genes <- normalize_id(
  background_genes
)


cat(
  "[INFO] WGD1 genes:",
  length(wgd1_genes),
  "\n"
)

cat(
  "[INFO] WGD2 genes:",
  length(wgd2_genes),
  "\n"
)

cat(
  "[INFO] Syntenic background genes:",
  length(background_genes),
  "\n"
)


write.table(
  wgd1_genes,
  file.path(
    outdir,
    paste0(
      sample_id,
      "_WGD1_genes.txt"
    )
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


write.table(
  wgd2_genes,
  file.path(
    outdir,
    paste0(
      sample_id,
      "_WGD2_genes.txt"
    )
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


ipr <- read.delim(
  ipr_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  fill = TRUE
)


if (ncol(ipr) < 14) {

  stop(
    "InterProScan TSV must contain at least 14 columns."
  )
}


# Normalize protein/gene identifiers before GO mapping

ipr$V1 <- normalize_id(
  ipr$V1
)


ipr_go <- ipr %>%

  filter(
    !is.na(V14),
    V14 != "-",
    V14 != ""
  ) %>%

  select(
    gene = V1,
    GO = V14
  ) %>%

  separate_rows(
    GO,
    sep = "\\||,"
  ) %>%

  mutate(

    GO = sub(
      "\\(InterPro\\)",
      "",
      GO
    ),

    GO = trimws(GO)

  ) %>%

  filter(
    grepl(
      "^GO:[0-9]+$",
      GO
    )
  ) %>%

  distinct()


term2gene <- ipr_go %>%

  select(
    term = GO,
    gene
  )


cat(
  "[INFO] Unique gene-GO associations:",
  nrow(term2gene),
  "\n"
)


go_annotated_genes <- unique(
  term2gene$gene
)


background_go <- intersect(
  background_genes,
  go_annotated_genes
)


wgd1_go_genes <- intersect(
  wgd1_genes,
  background_go
)


wgd2_go_genes <- intersect(
  wgd2_genes,
  background_go
)


cat(
  "[INFO] GO-annotated WGD1 genes:",
  length(wgd1_go_genes),
  "\n"
)

cat(
  "[INFO] GO-annotated WGD2 genes:",
  length(wgd2_go_genes),
  "\n"
)


run_enrichment <- function(
  genes,
  background,
  term2gene
) {

  if (length(genes) == 0) {
    return(NULL)
  }


  enricher(
    gene = genes,
    universe = background,
    TERM2GENE = term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
}


ego1 <- run_enrichment(
  genes = wgd1_go_genes,
  background = background_go,
  term2gene = term2gene
)


ego2 <- run_enrichment(
  genes = wgd2_go_genes,
  background = background_go,
  term2gene = term2gene
)

export_enrichment <- function(
  ego,
  filename
) {

  if (
    is.null(ego) ||
    nrow(as.data.frame(ego)) == 0
  ) {

    return(
      tibble()
    )
  }


  result <- as.data.frame(
    ego
  ) %>%

    as_tibble() %>%

    filter(
      !is.na(p.adjust),
      p.adjust <= 0.05
    ) %>%

    arrange(
      p.adjust
    )


  write.table(
    result,
    filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )


  result
}


wgd1_go_results <- export_enrichment(

  ego1,

  file.path(
    outdir,
    paste0(
      sample_id,
      "_WGD1_GO_enrichment.tsv"
    )
  )
)


wgd2_go_results <- export_enrichment(

  ego2,

  file.path(
    outdir,
    paste0(
      sample_id,
      "_WGD2_GO_enrichment.tsv"
    )
  )
)


if (nrow(wgd1_go_results) > 0) {

  p1 <- dotplot(
    ego1,
    showCategory = min(
      20,
      nrow(wgd1_go_results)
    )
  ) +

    ggtitle(
      paste0(
        sample_id,
        " WGD1 GO enrichment"
      )
    )


  ggsave(
    file.path(
      outdir,
      paste0(
        sample_id,
        "_WGD1_GO_dotplot.pdf"
      )
    ),
    p1,
    width = 8,
    height = 6
  )
}


if (nrow(wgd2_go_results) > 0) {

  p2 <- dotplot(
    ego2,
    showCategory = min(
      20,
      nrow(wgd2_go_results)
    )
  ) +

    ggtitle(
      paste0(
        sample_id,
        " WGD2 GO enrichment"
      )
    )


  ggsave(
    file.path(
      outdir,
      paste0(
        sample_id,
        "_WGD2_GO_dotplot.pdf"
      )
    ),
    p2,
    width = 8,
    height = 6
  )
}


wgd1_domains <- ipr %>%

  filter(
    V1 %in% wgd1_genes,
    !is.na(V13),
    V13 != "-",
    V13 != ""
  ) %>%

  count(
    V13,
    sort = TRUE,
    name = "count"
  ) %>%

  rename(
    InterPro_description = V13
  )


wgd2_domains <- ipr %>%

  filter(
    V1 %in% wgd2_genes,
    !is.na(V13),
    V13 != "-",
    V13 != ""
  ) %>%

  count(
    V13,
    sort = TRUE,
    name = "count"
  ) %>%

  rename(
    InterPro_description = V13
  )


write.table(
  wgd1_domains,
  file.path(
    outdir,
    paste0(
      sample_id,
      "_WGD1_top_InterPro_domains.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  wgd2_domains,
  file.path(
    outdir,
    paste0(
      sample_id,
      "_WGD2_top_InterPro_domains.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("Analysis completed:", sample_id, "\n")
cat("============================================================\n")

cat(
  "WGD1 Ks window:",
  wgd1_min,
  "-",
  wgd1_max,
  "\n"
)

cat(
  "WGD2 Ks window:",
  wgd2_min,
  "-",
  wgd2_max,
  "\n"
)

cat(
  "WGD1 syntenic blocks:",
  nrow(wgd1_blocks),
  "\n"
)

cat(
  "WGD2 syntenic blocks:",
  nrow(wgd2_blocks),
  "\n"
)

cat(
  "WGD1 genes:",
  length(wgd1_genes),
  "\n"
)

cat(
  "WGD2 genes:",
  length(wgd2_genes),
  "\n"
)

cat(
  "Background syntenic genes:",
  length(background_genes),
  "\n"
)

cat(
  "Significant WGD1 GO terms:",
  nrow(wgd1_go_results),
  "\n"
)

cat(
  "Significant WGD2 GO terms:",
  nrow(wgd2_go_results),
  "\n"
)

cat(
  "Output directory:",
  outdir,
  "\n"
)

cat("============================================================\n")
