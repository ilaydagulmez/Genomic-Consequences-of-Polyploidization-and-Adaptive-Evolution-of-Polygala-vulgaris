# Genomic Consequences of Polyploidization and Adaptive Evolution of *Polygala vulgaris*

This repository contains genomic resources, analysis scripts, and reproducible computational workflows associated with the comparative and evolutionary genomic analysis of tetraploid *Polygala vulgaris*.

The workflow integrates genome assembly and annotation, whole-genome duplication (WGD) inference, synteny and gene-level Ks analyses, orthology and gene-family evolution, functional annotation, and GO enrichment analyses to investigate the genomic consequences of polyploidization and the functional retention of duplicated genes.

---

## Overview

Five *Polygala vulgaris* individuals were analyzed:

| Sample  |
| ------- |
| 20628_1 |
| 20628_8 |
| 20628_9 |
| 20896_1 |
| 20896_2 |

Whole-genome Illumina paired-end sequencing data were used for genome assembly and downstream comparative genomic analyses.

Raw sequencing reads are publicly available through the **NCBI Sequence Read Archive (SRA)** under **BioProject PRJNA1060933**.

**NCBI SRA:**
https://www.ncbi.nlm.nih.gov/sra/PRJNA1060933

---

## Repository Structure

```text
.
├── data/
│   ├── proteins/
│   ├── cds/
│   ├── gff3/
│
├── scripts/
│   ├── 01_genome_assembly/
│   ├── 02_gene_annotation/
│   ├── 03_WGD_and_Ks/
│   ├── 04_collinearity_and_gene_level_Ks/
│   ├── 05_orthology_and_gene_family_evolution/
│   ├── 06_functional_annotation/
│   └── 07_GO_enrichment_and_functional_networks/
│
└── README.md
```

---

## Genomic Resources

The `data/` directory contains the principal genomic resources generated and used in this study.


Because of their large file sizes, genome assemblies are provided through external permanent download links rather than being stored directly in the GitHub repository.

| Sample  | Assembly                  |
| ------- | ------------------------- |
| 20628_1 | Download link to be added |
| 20628_8 | Download link to be added |
| 20628_9 | Download link to be added |
| 20896_1 | Download link to be added |
| 20896_2 | Download link to be added |

### Protein sequences

Predicted protein sequences are available in:

```text
data/proteins/
```

### Coding sequences

Predicted coding sequences (CDS) are available in:

```text
data/cds/
```

### Structural genome annotations

Final structural genome annotations in GFF3 format are available in:

```text
data/gff3/
```

---

## Raw Sequencing Data

Raw Illumina whole-genome sequencing reads are publicly available through the **NCBI Sequence Read Archive (SRA)**.

**NCBI BioProject: PRJNA1060933**

https://www.ncbi.nlm.nih.gov/sra/PRJNA1060933

| Sample  | NCBI SRA                |
| ------- | ----------------------- |
| 20628_1 | SRR32576552             |
| 20628_8 | BioProject PRJNA1060933 |
| 20628_9 | BioProject PRJNA1060933 |
| 20896_1 | BioProject PRJNA1060933 |
| 20896_2 | BioProject PRJNA1060933 |

Individual SRA Run accessions are linked through the corresponding BioProject record.

---

# Computational Workflow

The computational workflow implemented in this study is organized into seven major analytical stages.

```text
Genome assembly
      │
      ▼
Gene annotation
      │
      ▼
Whole-genome duplication and Ks analysis
      │
      ▼
Collinearity and gene-level Ks estimation
      │
      ▼
Orthology and gene-family evolution
      │
      ▼
Functional annotation
      │
      ▼
GO enrichment and functional meta-network analysis
```

---

## 01. Genome Assembly

Directory:

```text
scripts/01_genome_assembly/
```

Genome assemblies were generated and refined using a multi-step workflow involving:

* **MEGAHIT**
* **SOAPdenovo2**
* **Redundans**
* **Racon**
* **RagTag**

Assembly completeness was evaluated using **BUSCO**.

The complete assembly workflow is provided in:

```text
run_assembly.sh
```

---

## 02. Gene Annotation

Directory:

```text
scripts/02_gene_annotation/
```

Structural gene prediction was primarily performed using **Helixer** with the land-plant model.

Predicted CDS, protein, and transcript sequences were extracted from the resulting genome annotations.

**AUGUSTUS** is additionally included in the workflow as an optional alternative gene-prediction approach.

---

## 03. Whole-Genome Duplication and Ks Analysis

Directory:

```text
scripts/03_WGD_and_Ks/
```

Whole-genome duplication signals and synonymous substitution-rate distributions were investigated using **wgd v2**.

The workflow includes:

* sequence-similarity searches,
* Ks distribution estimation,
* syntenic anchor detection,
* and identification of WGD-associated Ks peaks.

Two major WGD signals were investigated:

| Event | Approximate age | Gene-level Ks interval used for meta-network analysis |
| ----- | --------------: | ----------------------------------------------------: |
| WGD1  |         ~55 MYA |                                             0.45–0.70 |
| WGD2  |         ~33 MYA |                                             0.24–0.40 |

---

## 04. Collinearity and Gene-Level Ks

Directory:

```text
scripts/04_collinearity_and_gene_level_Ks/
```

Gene-level synteny and Ks estimation were performed using the following workflow:

```text
DIAMOND
   ↓
MCScanX
   ↓
Collinear gene pairs
   ↓
MAFFT
   ↓
PAL2NAL
   ↓
KaKs_Calculator
   ↓
Gene-level Ks estimates
```

Protein similarity searches were performed with **DIAMOND**, and collinear blocks were detected using **MCScanX**.

For each collinear gene pair:

1. protein sequences were aligned using **MAFFT**,
2. codon-aware alignments were generated using **PAL2NAL**,
3. synonymous substitution rates were estimated with **KaKs_Calculator** using the **Yang–Nielsen (YN)** model.

The resulting gene-level Ks tables were used for downstream WGD-associated gene analyses.

---

## 05. Orthology and Gene-Family Evolution

Directory:

```text
scripts/05_orthology_and_gene_family_evolution/
```

Orthologous gene families were inferred using **OrthoFinder**.

Gene-family expansion and contraction analyses were subsequently performed using **CAFE5** with an ultrametric species tree.

The workflow includes:

* orthogroup inference,
* gene-family count generation,
* filtering of extreme gene families,
* and estimation of lineage-specific gene-family expansions and contractions.

---

## 06. Functional Annotation

Directory:

```text
scripts/06_functional_annotation/
```

Predicted proteins were functionally annotated using **InterProScan**.

The workflow retrieves:

* Pfam signatures,
* InterPro annotations,
* Gene Ontology (GO) terms,
* and associated functional annotations.

Protein-to-GO mappings generated from InterProScan output were used for subsequent enrichment analyses.

---

## 07. GO Enrichment and Functional Networks

Directory:

```text
scripts/07_GO_enrichment_and_functional_networks/
```

This directory contains two complementary analyses.

### WGD-associated GO enrichment

```text
01_WGD_GO_enrichment.R
```

This script:

1. integrates gene-level Ks estimates with MCScanX syntenic blocks,
2. calculates median Ks values for syntenic blocks,
3. identifies blocks associated with WGD Ks intervals,
4. extracts WGD-retained duplicated genes,
5. performs GO enrichment using `clusterProfiler`,
6. summarizes InterPro functional domains.

### Combined functional meta-network

```text
02_combined_functional_meta_network.R
```

This script integrates results from all five *P. vulgaris* individuals and both WGD intervals.

Significantly enriched GO terms are grouped into major functional modules:

* Stress / defense
* Transport / ion homeostasis
* Regulation / development
* Carbohydrate / cell wall
* Genome / proteostasis

Representative enriched GO terms are subsequently combined into a meta-network connecting functional terms with individual × WGD datasets.

The script generates:

* complete GO enrichment tables,
* selected representative GO-term tables,
* network node and edge tables,
* publication-quality PDF figures,
* 600-dpi PNG figures,
* and 600-dpi TIFF figures.

---

## Software

Major software and packages used in the workflow include:

| Analysis                               | Software                          |
| -------------------------------------- | --------------------------------- |
| Genome assembly                        | MEGAHIT, SOAPdenovo2              |
| Assembly refinement                    | Redundans, Racon, RagTag          |
| Assembly assessment                    | BUSCO                             |
| Gene prediction                        | Helixer                           |
| Optional gene prediction               | AUGUSTUS                          |
| WGD analysis                           | wgd v2                            |
| Similarity searches                    | DIAMOND                           |
| Collinearity                           | MCScanX                           |
| Protein alignment                      | MAFFT                             |
| Codon alignment                        | PAL2NAL                           |
| Ka/Ks estimation                       | KaKs_Calculator                   |
| Orthology inference                    | OrthoFinder                       |
| Gene-family evolution                  | CAFE5                             |
| Functional annotation                  | InterProScan                      |
| GO enrichment                          | clusterProfiler                   |
| Statistical analysis and visualization | R / tidyverse / ggplot2 / ggrepel |

Detailed software versions and parameters are provided within the corresponding scripts and the Methods section of the associated manuscript.

---

## Reproducibility

The scripts in this repository are organized according to the analytical workflow used in the associated study.

Where possible, sample-specific and computing-system-specific paths were removed so that the workflows can be adapted to other systems.

Computational requirements vary considerably among analyses. Genome assembly, InterProScan, OrthoFinder, and large-scale pairwise Ks estimation are recommended to be executed in a high-performance computing environment.

---

## Data Availability

Raw whole-genome sequencing reads have been deposited in the **NCBI Sequence Read Archive (SRA)** under **BioProject PRJNA1060933**.

The repository additionally provides:

* final structural genome annotations (GFF3),
* predicted protein sequences,
* predicted coding sequences (CDS),
* computational workflows,
* and scripts required to reproduce the principal comparative and evolutionary genomic analyses.

Final genome assemblies are available through the permanent download links provided in the `data/assemblies/` directory.

---

## Contact

For questions regarding the genomic resources or analytical workflows, please contact:

**İlayda Gülmez**

e-mail: ilaydagulmez@hacettepe.edu.tr

Molecular Plant Systematic Laboratory (MOBIS)
Department of Biology, Faculty of Science
Hacettepe University
Ankara, Türkiye
