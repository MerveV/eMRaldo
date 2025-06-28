# eMRaldo: MR Signaling Exploration Platform for Transcriptomic Data

eMRaldo is an interactive Shiny-based web application developed to enable exploration of transcriptomic datasets involving mineralocorticoid receptor (MR/NR3C2) overexpression and aldosterone treatments. The app is designed to support both RNA-seq and microarray datasets sourced from publicly available experiments in NCBI GEO ([link](https://www.ncbi.nlm.nih.gov/geo/)).

## 🔬 Key Features

- Analyze datasets using **limma** (Ritchie et al., 2015) for two-group contrasts
- Define your own gene selection thresholds or use pre-integrated gene lists (e.g., KEGG)
- Interactive visualization of results with:
  - PCA plots and volcano plots
  - Expression heatmaps of selected genes
  - Gene set enrichment analysis (MSigDB)
  - DisGeNet disease association networks
  - Transcription factor (TF) activity heatmaps and networks (Dorothea + mSigDB)
- Custom scatterplots between contrasts for exploratory gene comparison
- Data quality summaries and QC plots in the Data Summary tab

## 🧪 Analysis Workflow Overview

After selecting a dataset:
1. Choose from four contrast comparisons (e.g., control vs aldosterone)
2. Click **"Analyze"** to compute differential expression
3. Use the result table to:
   - Apply user-defined gene filtering thresholds
   - Choose from curated gene sets (e.g., KEGG)
4. Visualize selected genes via:
   - Heatmaps (logFC values)
   - Scatterplots
   - Gene set enrichment tables → interactive heatmaps of MSigDB terms
   - DisGeNet and TF network views


## 📂 Dataset Summary

- **Cell lines**: MCF7 and ZR75-1 breast cancer models
- **Treatments**: Aldosterone (100nM) and Spironolactone (10µM)
- **Design**: Transfection with MR (OV) or empty vector (EV), followed by hormone treatments
- **Platform**: Illumina NextSeq500 (75bp single-end RNA-seq)

## 📦 Installation (Local use)

This app depends on the following R packages:

```r
# Required packages
packages <- c(
  "shiny", "shinydashboard", "shinydashboardPlus", "rintrojs", "fresh",
  "shinyWidgets", "DT", "heatmaply", "plotly", "limma", "edgeR", "DESeq2",
  "visNetwork", "shinyBS", "shinycssloaders", "EnrichmentBrowser",
  "stringr", "markdown", "msigdbr", "hypeR", "enrichplot", "clusterProfiler",
  "DOSE", "gprofiler2", "RColorBrewer", "EnhancedVolcano", "shinyalert"
)

# Install if not already installed
install.packages(setdiff(packages, rownames(installed.packages())))

# Load the app
shiny::runApp(".")
# Install required packages
install.packages(c("shiny", "ggplot2", "dplyr", "readr", "DT", "plotly"))

# Run the app
shiny::runApp(".")
```
Or open `app.R` in RStudio and click **Run App**.

## 🧬 Experimental Protocols

- **Growth Protocol**: MCF7 and ZR75-1 cells cultured at 37°C, 5% CO₂ in DMEM or RPMI with FBS, NEAA, Na-pyruvate
- **Treatment Protocol**: Transfection with MR or EV vectors, followed by exposure to aldosterone and/or spironolactone for 24 hours
- **RNA Extraction**: Total RNA isolated using QIAZOL, RIN > 7
- **Library Construction**: Prepared and sequenced on Illumina NextSeq500 (EMBL Genomics Core)

## 🖼️ Example Outputs

- PCA and volcano plots for overall expression patterns
- Interactive heatmaps of filtered gene sets
- Enrichment tables and annotated networks from DisGeNet and MSigDB
- Scatterplots comparing logFC across contrasts

## 📊 Data Summary Tab

The app includes a dedicated section for examining:
- Quality control metrics (read counts, mapping rates, RIN values)
- Normalized expression distributions
- Sample similarity and clustering
- Differential expression statistics

## 📨 Contact

If you have questions, feedback, or suggestions, please contact:  
[mervevuralozdeniz@gmail.com](mailto:mervevuralozdeniz@gmail.com)

## 📄 License

This project is released under the [MIT License](LICENSE).

---

*eMRaldo was developed to facilitate MR-related transcriptomic analyses and enable hypothesis generation through integrative visualization and analysis tools.*
