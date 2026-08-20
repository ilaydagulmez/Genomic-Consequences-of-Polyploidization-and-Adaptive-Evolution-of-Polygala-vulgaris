#!/usr/bin/env Rscript

# Required R packages:
#   tidyverse
#   clusterProfiler
#   GO.db
#   AnnotationDbi
#   scales
#   ggrepel

suppressPackageStartupMessages({

  library(tidyverse)
  library(clusterProfiler)
  library(GO.db)
  library(AnnotationDbi)
  library(scales)
  library(ggrepel)

  # Reload after AnnotationDbi to avoid namespace conflicts
  library(dplyr)

})

args <- commandArgs(
  trailingOnly = TRUE
)


if (length(args) < 1) {

  stop(
    paste0(
      "\nUsage:\n",
      "Rscript 02_combined_functional_meta_network.R ",
      "<input_dir> [output_dir]\n"
    )
  )
}


input_dir <- args[1]


output_dir <- ifelse(
  length(args) >= 2,
  args[2],
  "combined_meta_network_results"
)


dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

samples <- c(
  "20628_1",
  "20628_8",
  "20628_9",
  "20896_1",
  "20896_2"
)

find_existing_file <- function(paths) {

  existing <- paths[
    file.exists(paths)
  ]


  if (length(existing) == 0) {

    stop(
      paste0(
        "None of the expected input files were found:\n",
        paste(
          paths,
          collapse = "\n"
        )
      )
    )
  }


  existing[1]
}


sample_config <- tibble(
  sample_name = samples
) %>%

  mutate(

    ks_file = map_chr(
      sample_name,
      function(x) {

        find_existing_file(
          c(
            file.path(
              input_dir,
              paste0(
                x,
                "_ks_genelevel.tsv"
              )
            )
          )
        )
      }
    ),


    ipr_file = map_chr(
      sample_name,
      function(x) {

        find_existing_file(
          c(

            file.path(
              input_dir,
              paste0(
                x,
                "_interproscan.tsv"
              )
            ),

            file.path(
              input_dir,
              paste0(
                "AGAT_",
                x,
                "_protein.fa.tsv"
              )
            ),

            file.path(
              input_dir,
              paste0(
                x,
                "_protein.fa.tsv"
              )
            )

          )
        )
      }
    )

  )


print(
  sample_config
)

wgd_windows <- tibble(

  WGD = c(
    "WGD1",
    "WGD2"
  ),

  age_label = c(
    "55 MYA",
    "33 MYA"
  ),

  ks_min = c(
    0.45,
    0.24
  ),

  ks_max = c(
    0.70,
    0.40
  )

)

out_prefix <- "Polygala_5individual_combined_meta_GO_network"

padj_method <- "BH"

padj_cutoff <- 0.05

include_other_module <- FALSE

terms_per_module_per_set <- 2

max_terms_per_set <- 8

term_wrap_width <- 18


set.seed(
  123
)

module_fill_cols <- c(

  "Stress / defense" =
    "#7DA7E8",

  "Transport / ion homeostasis" =
    "#4FC3B3",

  "Regulation / development" =
    "#E78670",

  "Carbohydrate / cell wall" =
    "#F2C14E",

  "Genome / proteostasis" =
    "#A38AD7",

  "Other" =
    "#C8D0DE"

)

module_edge_cols <- c(

  "Stress / defense" =
    "#3E6FB1",

  "Transport / ion homeostasis" =
    "#1D9484",

  "Regulation / development" =
    "#CC5C45",

  "Carbohydrate / cell wall" =
    "#D59A00",

  "Genome / proteostasis" =
    "#6E58B0",

  "Other" =
    "#98A3B5"

)

sample_cols <- c(

  "20628_1" =
    "#1F77B4",

  "20628_8" =
    "#2CA02C",

  "20628_9" =
    "#D6278B",

  "20896_1" =
    "#FF7F0E",

  "20896_2" =
    "#8C564B"

)

preferred_module_order <- c(

  "Transport / ion homeostasis",

  "Stress / defense",

  "Regulation / development",

  "Carbohydrate / cell wall",

  "Genome / proteostasis",

  "Other"

)


out_file <- function(suffix) {

  file.path(
    output_dir,
    paste0(
      out_prefix,
      suffix
    )
  )
}

clean_gene_id <- function(x) {

  x %>%

    as.character() %>%

    str_replace(
      "\\s.*$",
      ""
    ) %>%

    str_replace(
      "^>",
      ""
    ) %>%

    str_replace(
      "\\.[0-9]+$",
      ""
    ) %>%

    str_replace(
      "X[0-9]+$",
      ""
    )
}


safe_go_term <- function(go_id) {

  obj <- GO.db::GOTERM[
    [go_id]
  ]


  if (is.null(obj)) {

    return(
      go_id
    )
  }


  out <- AnnotationDbi::Term(
    obj
  )


  if (
    length(out) == 0 ||
    is.na(out)
  ) {

    return(
      go_id
    )
  }


  out
}


assign_module <- function(desc) {

  case_when(

    str_detect(
      desc,
      regex(
        paste0(
          "ion|cation|potassium|calcium|metal|",
          "transmembrane|transport|antiporter|channel|",
          "ABC|nitrate|NRT|PTR|membrane|carrier"
        ),
        ignore_case = TRUE
      )
    ) ~

      "Transport / ion homeostasis",


    str_detect(
      desc,
      regex(
        paste0(
          "stress|salt|drought|water deprivation|",
          "osmotic|oxidative|detox|defense|immune|",
          "pathogen|response to|jasmonate response"
        ),
        ignore_case = TRUE
      )
    ) ~

      "Stress / defense",


    str_detect(
      desc,
      regex(
        paste0(
          "hormone|auxin|abscisic|ABA|jasmonate|",
          "ethylene|gibberellin|cytokinin|signaling|",
          "transcription|zinc finger|WRKY|MYB|Dof|B3|",
          "root|leaf|epidermis|meristem|growth|development|",
          "gravitropic|protein kinase|phosphatase|receptor"
        ),
        ignore_case = TRUE
      )
    ) ~

      "Regulation / development",


    str_detect(
      desc,
      regex(
        paste0(
          "starch|sucrose|glucan|glucose|fructose|mannose|",
          "carbohydrate|carbon|cellulose|pectin|xylan|lignin|",
          "phenylpropanoid|cell wall|cuticle|lipid|fatty acid"
        ),
        ignore_case = TRUE
      )
    ) ~

      "Carbohydrate / cell wall",


    str_detect(
      desc,
      regex(
        paste0(
          "DNA|chromatin|recombination|repair|replication|",
          "transpos|ubiquitin|proteasome|chaperone|Clp|",
          "translation|ribosom|RNA modification|",
          "protein folding|proteolysis"
        ),
        ignore_case = TRUE
      )
    ) ~

      "Genome / proteostasis",


    TRUE ~

      "Other"

  )
}

build_module_centers <- function(
  modules_present
) {

  modules_present <- intersect(
    preferred_module_order,
    unique(
      modules_present
    )
  )


  n <- length(
    modules_present
  )


  coords <- switch(

    as.character(
      n
    ),


    "1" = tibble(
      cx = 0,
      cy = 0
    ),


    "2" = tibble(
      cx = c(
        -4.2,
        4.2
      ),
      cy = c(
        0,
        0
      )
    ),


    "3" = tibble(
      cx = c(
        -5.0,
        0,
        5.0
      ),
      cy = c(
        0,
        3.2,
        0
      )
    ),


    "4" = tibble(
      cx = c(
        -5.2,
        5.2,
        -5.2,
        5.2
      ),
      cy = c(
        2.8,
        2.8,
        -2.8,
        -2.8
      )
    ),


    "5" = tibble(
      cx = c(
        -5.6,
        0,
        5.6,
        -3.4,
        3.4
      ),
      cy = c(
        3.0,
        3.8,
        3.0,
        -2.8,
        -2.8
      )
    ),


    tibble(
      cx = c(
        -5.6,
        0,
        5.6,
        -3.4,
        3.4,
        0
      ),
      cy = c(
        3.0,
        3.8,
        3.0,
        -2.8,
        -2.8,
        -4.6
      )
    )

  )


  tibble(
    module = modules_present
  ) %>%

    bind_cols(
      coords[
        seq_len(n),
      ]
    )
}

safe_rescale <- function(
  x,
  to
) {

  finite_x <- x[
    is.finite(x)
  ]


  if (
    length(finite_x) == 0
  ) {

    return(
      rep(
        mean(to),
        length(x)
      )
    )
  }


  if (
    length(
      unique(finite_x)
    ) == 1
  ) {

    return(
      rep(
        mean(to),
        length(x)
      )
    )
  }


  scales::rescale(
    x,
    to = to,
    from = range(
      finite_x,
      na.rm = TRUE
    )
  )
}

run_one_sample <- function(
  sample_name,
  ks_file,
  ipr_file
) {

  message(
    "Running ",
    sample_name
  )

  ks <- read_tsv(
    ks_file,
    show_col_types = FALSE
  ) %>%

    mutate(

      gene1 = clean_gene_id(
        gene1
      ),

      gene2 = clean_gene_id(
        gene2
      ),

      Ks = as.numeric(
        Ks
      )

    ) %>%

    filter(
      !is.na(Ks),
      is.finite(Ks)
    )

  ipr <- read_tsv(
    ipr_file,
    col_names = FALSE,
    comment = "#",
    show_col_types = FALSE
  )

  if (
    ncol(ipr) < 14
  ) {

    stop(
      "InterProScan TSV must contain at least 14 columns: ",
      ipr_file
    )
  }


  ipr <- ipr[
    ,
    1:min(
      ncol(ipr),
      15
    )
  ]

  colnames(ipr)[1:14] <- c(

    "gene_raw",

    "md5",

    "seq_length",

    "analysis",

    "signature_accession",

    "signature_description",

    "start",

    "stop",

    "score",

    "status",

    "date",

    "ipr_accession",

    "ipr_description",

    "go_raw"

  )


  if (
    ncol(ipr) >= 15
  ) {

    colnames(ipr)[15] <-
      "pathway_raw"
  }


  ipr2 <- ipr %>%

    mutate(
      gene = clean_gene_id(
        gene_raw
      )
    )

  gene2go <- ipr2 %>%

    filter(
      !is.na(go_raw),
      go_raw != "-",
      go_raw != ""
    ) %>%

    separate_rows(
      go_raw,
      sep = "\\|"
    ) %>%

    mutate(
      GO = str_extract(
        go_raw,
        "GO:\\d+"
      )
    ) %>%

    filter(
      !is.na(GO)
    ) %>%

    distinct(
      GO,
      gene
    )

  term2name <- tibble(
    GO = unique(
      gene2go$GO
    )
  ) %>%

    mutate(
      Description = map_chr(
        GO,
        safe_go_term
      )
    ) %>%

    distinct(
      GO,
      Description
    )

  run_go_enrichment <- function(
    genes,
    wgd_name
  ) {

    g <- intersect(
      unique(
        genes
      ),
      background_genes
    )


    if (
      length(g) == 0
    ) {

      return(
        tibble()
      )
    }


    ego <- enricher(

      gene = g,

      universe = background_genes,

      TERM2GENE = gene2go,

      TERM2NAME = term2name,

      pAdjustMethod =
        padj_method,

      pvalueCutoff =
        1,

      qvalueCutoff =
        1

    )


    if (
      is.null(ego) ||
      nrow(
        as.data.frame(
          ego
        )
      ) == 0
    ) {

      return(
        tibble()
      )
    }


    as.data.frame(
      ego
    ) %>%

      as_tibble() %>%

      mutate(

        sample_name =
          sample_name,

        WGD =
          wgd_name,

        set_id =
          paste(
            sample_name,
            wgd_name,
            sep = " | "
          ),

        p.adjust.safe =
          pmax(
            p.adjust,
            .Machine$double.xmin
          ),

        module =
          assign_module(
            Description
          ),

        minus_log10_padj =
          -log10(
            p.adjust.safe
          )

      ) %>%

      filter(
        !is.na(
          p.adjust
        ),
        p.adjust <=
          padj_cutoff
      )
  }


  map_dfr(

    seq_len(
      nrow(
        wgd_windows
      )
    ),

    function(i) {

      wgd_name <-
        wgd_windows$WGD[i]

      ks_min <-
        wgd_windows$ks_min[i]

      ks_max <-
        wgd_windows$ks_max[i]


      genes <- ks %>%

        filter(
          Ks >= ks_min,
          Ks <= ks_max
        ) %>%

        dplyr::select(
          gene1,
          gene2
        ) %>%

        pivot_longer(
          cols = everything(),
          values_to = "gene"
        ) %>%

        distinct(
          gene
        ) %>%

        pull(
          gene
        )


      run_go_enrichment(
        genes,
        wgd_name
      )
    }
  )
}

all_results <- pmap_dfr(

  sample_config,

  function(
    sample_name,
    ks_file,
    ipr_file
  ) {

    run_one_sample(
      sample_name,
      ks_file,
      ipr_file
    )
  }

)


if (
  nrow(
    all_results
  ) == 0
) {

  stop(
    "No significant GO terms were detected across all samples."
  )
}


write_csv(

  all_results,

  out_file(
    "_all_GO_enrichment_results.csv"
  )

)


write_tsv(

  wgd_windows,

  out_file(
    "_WGD_windows.tsv"
  )

)

selected_per_set <- all_results %>%

  {

    if (
      !include_other_module
    ) {

      filter(
        .,
        module != "Other"
      )

    } else {

      .

    }
  } %>%

  group_by(
    set_id,
    module
  ) %>%

  slice_min(

    order_by =
      p.adjust,

    n =
      terms_per_module_per_set,

    with_ties =
      FALSE

  ) %>%

  ungroup() %>%

  group_by(
    set_id
  ) %>%

  slice_min(

    order_by =
      p.adjust,

    n =
      max_terms_per_set,

    with_ties =
      FALSE

  ) %>%

  ungroup()


if (
  nrow(
    selected_per_set
  ) == 0
) {

  stop(
    "No GO terms remained after representative-term selection."
  )
}


write_csv(

  selected_per_set,

  out_file(
    "_selected_terms_per_set.csv"
  )

)


term_summary <- selected_per_set %>%

  group_by(
    ID,
    Description,
    module
  ) %>%

  summarise(

    support_n =
      n_distinct(
        set_id
      ),

    mean_sig =
      mean(
        minus_log10_padj,
        na.rm = TRUE
      ),

    min_padj =
      min(
        p.adjust,
        na.rm = TRUE
      ),

    .groups =
      "drop"

  )


term_summary$scaled_sig <- safe_rescale(

  term_summary$mean_sig,

  to = c(
    0.6,
    3.2
  )

)


term_summary <- term_summary %>%

  mutate(

    node_key =
      paste0(
        "TERM__",
        ID
      ),

    label =
      str_wrap(
        str_to_title(
          Description
        ),
        width =
          term_wrap_width
      ),

    score =
      support_n +
      scaled_sig

  )

set_nodes <- selected_per_set %>%

  distinct(
    set_id,
    sample_name,
    WGD
  ) %>%

  mutate(

    node_key =
      paste0(
        "SET__",
        set_id
      ),

    sample_col =
      sample_cols[
        sample_name
      ]

  )

edge_tbl <- selected_per_set %>%

  transmute(

    from =
      paste0(
        "TERM__",
        ID
      ),

    to =
      paste0(
        "SET__",
        set_id
      ),

    module =
      module,

    p.adjust =
      p.adjust,

    set_id =
      set_id,

    sample_name =
      sample_name,

    WGD =
      WGD

  ) %>%

  distinct()

term_nodes <- term_summary %>%

  transmute(

    node_key,

    node_type =
      "term",

    label,

    module,

    support_n,

    mean_sig,

    score

  )


set_nodes2 <- set_nodes %>%

  transmute(

    node_key,

    node_type =
      "set",

    label =
      set_id,

    module =
      "set",

    sample_name,

    sample_col,

    support_n =
      1,

    mean_sig =
      NA_real_,

    score =
      1

  )


nodes <- bind_rows(

  term_nodes,

  set_nodes2

)

centers <- build_module_centers(

  unique(
    term_nodes$module
  )

)


term_pos <- term_nodes %>%

  left_join(
    centers,
    by = "module"
  ) %>%

  group_by(
    module
  ) %>%

  arrange(

    desc(
      support_n
    ),

    desc(
      mean_sig
    ),

    .by_group =
      TRUE

  ) %>%

  mutate(

    idx =
      row_number(),

    n_terms =
      n(),

    angle =
      if_else(
        n_terms == 1,
        pi / 2,
        2 * pi *
          (idx - 1) /
          n_terms
      ),

    radius =
      case_when(

        n_terms == 1 ~
          0,

        n_terms == 2 ~
          1.7,

        n_terms <= 4 ~
          2.3,

        TRUE ~
          2.9

      ),

    x =
      cx +
      radius *
      cos(
        angle
      ),

    y =
      cy +
      radius *
      sin(
        angle
      )

  ) %>%

  ungroup() %>%

  dplyr::select(
    node_key,
    x,
    y,
    module
  )

set_pos <- edge_tbl %>%

  left_join(

    term_pos %>%

      dplyr::select(
        from = node_key,
        tx = x,
        ty = y
      ),

    by = "from"

  ) %>%

  group_by(
    node_key = to
  ) %>%

  summarise(

    x =
      mean(
        tx,
        na.rm = TRUE
      ),

    y =
      mean(
        ty,
        na.rm = TRUE
      ),

    .groups =
      "drop"

  ) %>%

  mutate(

    r =
      sqrt(
        x^2 +
        y^2
      ),

    x =
      x +
      ifelse(
        r == 0,
        0,
        0.9 * x / r
      ) +
      rnorm(
        n(),
        0,
        0.55
      ),

    y =
      y +
      ifelse(
        r == 0,
        0,
        0.9 * y / r
      ) +
      rnorm(
        n(),
        0,
        0.55
      )

  ) %>%

  dplyr::select(
    node_key,
    x,
    y
  )

node_pos <- bind_rows(

  term_pos %>%

    dplyr::select(
      node_key,
      x,
      y
    ),

  set_pos %>%

    dplyr::select(
      node_key,
      x,
      y
    )

) %>%

  right_join(

    nodes %>%

      dplyr::select(
        node_key
      ),

    by =
      "node_key"

  ) %>%

  mutate(

    x =
      ifelse(
        is.na(x),
        rnorm(
          n(),
          0,
          0.5
        ),
        x
      ),

    y =
      ifelse(
        is.na(y),
        rnorm(
          n(),
          0,
          0.5
        ),
        y
      )

  )


nodes_for_plot <- nodes %>%

  left_join(
    node_pos,
    by = "node_key"
  )


term_index <- which(
  nodes_for_plot$node_type ==
    "term"
)


nodes_for_plot$size_plot <-
  3.9


nodes_for_plot$size_plot[
  term_index
] <- safe_rescale(

  nodes_for_plot$score[
    term_index
  ],

  to = c(
    15,
    31
  )

)

edge_plot <- edge_tbl %>%

  left_join(

    nodes_for_plot %>%

      dplyr::select(
        from = node_key,
        x_from = x,
        y_from = y
      ),

    by =
      "from"

  ) %>%

  left_join(

    nodes_for_plot %>%

      dplyr::select(
        to = node_key,
        x_to = x,
        y_to = y
      ),

    by =
      "to"

  )


write_tsv(

  nodes_for_plot,

  out_file(
    "_nodes.tsv"
  )

)


write_tsv(

  edge_tbl,

  out_file(
    "_edges.tsv"
  )

)

sample_legend_df <- tibble(

  sample_name =
    names(
      sample_cols
    ),

  col =
    unname(
      sample_cols
    ),

  x =
    12.4,

  y =
    seq(
      4.4,
      0.4,
      length.out =
        length(
          sample_cols
        )
    )

)

term_label_df <- nodes_for_plot %>%

  filter(
    node_type ==
      "term"
  ) %>%

  mutate(

    raw_label =
      str_replace_all(
        label,
        "\\n",
        " "
      ),

    force_inside =
      str_detect(

        raw_label,

        regex(
          paste0(
            "DNA.?Binding.*Transcription.*Factor|",
            "Response.*Auxin|",
            "Pectinesterase|",
            "Cell Wall.*Modification"
          ),
          ignore_case =
            TRUE
        )

      ),


    label =
      dplyr::case_when(

        str_detect(
          raw_label,
          regex(
            "DNA.?Binding.*Transcription.*Factor",
            ignore_case = TRUE
          )
        ) ~

          "DNA-Binding\nTranscription\nFactor\nActivity",


        str_detect(
          raw_label,
          regex(
            "Response.*Auxin",
            ignore_case = TRUE
          )
        ) ~

          "Response To\nAuxin",


        str_detect(
          raw_label,
          regex(
            "Pectinesterase",
            ignore_case = TRUE
          )
        ) ~

          "Pectinesterase\nActivity",


        str_detect(
          raw_label,
          regex(
            "Cell Wall.*Modification",
            ignore_case = TRUE
          )
        ) ~

          "Cell Wall\nModification",


        str_detect(
          raw_label,
          regex(
            "Sequence.?Specific.*DNA",
            ignore_case = TRUE
          )
        ) ~

          "Sequence-Specific\nDNA Binding",


        TRUE ~
          label
      ),


    module_col =
      module_edge_cols[
        module
      ],


    label_nchar =
      nchar(
        gsub(
          "\\n",
          " ",
          label
        )
      ),


    label_nlines =
      stringr::str_count(
        label,
        "\\n"
      ) +
      1,

    label_mode =
      dplyr::case_when(

        force_inside ~
          "inside",

        str_detect(
          raw_label,
          regex(
            "Sequence.?Specific.*DNA",
            ignore_case = TRUE
          )
        ) ~
          "outside",

        size_plot >= 23 &
          label_nchar <= 24 &
          label_nlines <= 3 ~
          "inside",

        size_plot >= 20 &
          label_nchar <= 16 &
          label_nlines <= 2 ~
          "inside",

        TRUE ~
          "outside"

      ),


    inside_size =
      dplyr::case_when(

        str_detect(
          label,
          regex(
            "DNA-Binding",
            ignore_case = TRUE
          )
        ) ~
          2.55,

        str_detect(
          label,
          regex(
            "Response To\\nAuxin",
            ignore_case = TRUE
          )
        ) ~
          2.75,

        label_nlines >= 3 ~
          2.65,

        TRUE ~
          3.05

      ),


    r =
      sqrt(
        x^2 +
        y^2
      ),


    r =
      ifelse(
        r == 0,
        1,
        r
      ),


    nudge_x =
      dplyr::case_when(

        str_detect(
          raw_label,
          regex(
            "Sequence.?Specific.*DNA",
            ignore_case = TRUE
          )
        ) ~
          1.70,

        TRUE ~
          1.55 *
          x /
          r

      ),


    nudge_y =
      dplyr::case_when(

        str_detect(
          raw_label,
          regex(
            "Sequence.?Specific.*DNA",
            ignore_case = TRUE
          )
        ) ~
          -0.18,

        TRUE ~
          1.30 *
          y /
          r

      ),


    nudge_y =
      ifelse(

        abs(y) < 1.0 &

          !str_detect(
            raw_label,
            regex(
              "Sequence.?Specific.*DNA",
              ignore_case = TRUE
            )
          ),

        nudge_y +
          0.55 *
          sign(
            ifelse(
              y == 0,
              1,
              y
            )
          ),

        nudge_y

      )

  )

inside_term_labels <-
  term_label_df %>%

  dplyr::filter(
    label_mode ==
      "inside"
  ) %>%

  dplyr::mutate(

    y_label =
      dplyr::case_when(

        str_detect(
          label,
          regex(
            "Sequence-Specific",
            ignore_case = TRUE
          )
        ) ~

          y - 0.38,

        TRUE ~
          y

      )

  )


inside_term_labels_long <-
  inside_term_labels %>%

  dplyr::filter(
    inside_size <
      2.9
  )


inside_term_labels_regular <-
  inside_term_labels %>%

  dplyr::filter(
    inside_size >=
      2.9
  )


outside_term_labels <-
  term_label_df %>%

  dplyr::filter(
    label_mode ==
      "outside"
  )

combined_meta_network <-

  ggplot() +

  geom_curve(

    data =
      edge_plot,

    aes(

      x =
        x_from,

      y =
        y_from,

      xend =
        x_to,

      yend =
        y_to,

      colour =
        module

    ),

    curvature =
      0.12,

    linewidth =
      0.62,

    alpha =
      0.52,

    lineend =
      "round"

  ) +


  geom_point(

    data =
      nodes_for_plot %>%
      filter(
        node_type ==
          "set"
      ),

    aes(
      x =
        x,
      y =
        y
    ),

    shape =
      21,

    size =
      3.9,

    fill =
      set_nodes2$sample_col[

        match(

          (
            nodes_for_plot %>%
            filter(
              node_type ==
                "set"
            )
          )$node_key,

          set_nodes2$node_key

        )

      ],

    colour =
      "white",

    stroke =
      0.55,

    alpha =
      0.98

  ) +

  geom_point(

    data =
      nodes_for_plot %>%
      filter(
        node_type ==
          "term"
      ),

    aes(

      x =
        x,

      y =
        y,

      size =
        size_plot,

      fill =
        module

    ),

    shape =
      21,

    colour =
      "#F5F7FA",

    stroke =
      1.45,

    alpha =
      0.99

  ) +

  geom_text(

    data =
      inside_term_labels_regular,

    aes(

      x =
        x,

      y =
        y_label,

      label =
        label

    ),

    family =
      "serif",

    fontface =
      "bold",

    size =
      3.05,

    lineheight =
      0.88,

    colour =
      "#1F2937"

  ) +

  geom_text(

    data =
      inside_term_labels_long,

    aes(

      x =
        x,

      y =
        y_label,

      label =
        label

    ),

    family =
      "serif",

    fontface =
      "bold",

    size =
      2.55,

    lineheight =
      0.82,

    colour =
      "#1F2937"

  ) +

  ggrepel::geom_text_repel(

    data =
      outside_term_labels,

    aes(

      x =
        x,

      y =
        y,

      label =
        label

    ),

    nudge_x =
      outside_term_labels$nudge_x,

    nudge_y =
      outside_term_labels$nudge_y,

    family =
      "serif",

    fontface =
      "bold",

    size =
      3.05,

    colour =
      "#1F2937",

    box.padding =
      0.68,

    point.padding =
      2.65,

    segment.colour =
      alpha(
        "#AAB3C2",
        0.95
      ),

    segment.size =
      0.34,

    segment.curvature =
      0.06,

    segment.ncp =
      2,

    min.segment.length =
      0,

    max.overlaps =
      Inf,

    seed =
      123,

    force =
      2.65,

    force_pull =
      0.06,

    direction =
      "both",

    bg.color =
      alpha(
        "white",
        0.86
      ),

    bg.r =
      0.10

  ) +


  geom_text(

    data =
      tibble(
        x = 12.4,
        y = 5.25,
        lab = "Samples"
      ),

    aes(
      x =
        x,
      y =
        y,
      label =
        lab
    ),

    inherit.aes =
      FALSE,

    family =
      "serif",

    fontface =
      "bold",

    size =
      4.1,

    hjust =
      0,

    colour =
      "#111827"

  ) +

  geom_point(

    data =
      sample_legend_df,

    aes(
      x =
        x,
      y =
        y
    ),

    inherit.aes =
      FALSE,

    shape =
      21,

    size =
      4.0,

    fill =
      sample_legend_df$col,

    colour =
      "white",

    stroke =
      0.55

  ) +


  # ----------------------------------------------------------
  # Sample legend labels
  # ----------------------------------------------------------

  geom_text(

    data =
      sample_legend_df,

    aes(

      x =
        x + 0.42,

      y =
        y,

      label =
        sample_name

    ),

    inherit.aes =
      FALSE,

    family =
      "sans",

    size =
      3.2,

    hjust =
      0,

    colour =
      "#374151"

  ) +

  scale_fill_manual(

    values =
      module_fill_cols,

    name =
      "Functional Module",

    drop =
      TRUE

  ) +


  scale_colour_manual(

    values =
      module_edge_cols,

    guide =
      "none",

    drop =
      TRUE

  ) +


  scale_size_identity() +

  coord_equal(

    xlim =
      c(
        -10.8,
        15.0
      ),

    ylim =
      c(
        -6.6,
        6.7
      ),

    clip =
      "off"

  ) +


  theme_void(
    base_family =
      "serif"
  ) +


  theme(

    legend.position =
      "bottom",

    legend.title =
      element_text(
        size = 8.8,
        face = "plain"
      ),

    legend.text =
      element_text(
        size = 8.8
      ),

    plot.title =
      element_text(
        hjust = 0.5,
        face = "bold",
        size = 18,
        colour = "#111827"
      ),

    plot.subtitle =
      element_text(
        hjust = 0.5,
        size = 10.2,
        colour = "#4B5563"
      ),

    plot.caption =
      element_text(
        hjust = 0.5,
        size = 8.8,
        colour = "#6B7280"
      ),

    plot.margin =
      margin(
        10,
        12,
        10,
        10
      )

  ) +


  guides(

    fill =
      guide_legend(

        override.aes =
          list(

            size =
              4.8,

            colour =
              "#F5F7FA"

          ),

        nrow =
          1

      )

  ) +


  labs(

    title =
      expression(
        italic(
          "Polygala vulgaris"
        ) *
          " Combined Functional Meta-Network of WGD-Retained Duplicated Genes"
      ),

    subtitle =
      paste0(
        "Large nodes are enriched GO terms pooled across five individuals ",
        "and two WGD windows; small colored nodes are sample-WGD sets."
      ),

    caption =
      paste0(
        "Node size reflects recurrence and mean enrichment strength ",
        "across sample-WGD sets."
      )

  )

print(
  combined_meta_network
)


