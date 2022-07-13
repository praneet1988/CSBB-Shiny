[![Github All Releases](https://img.shields.io/github/downloads/praneet1988/CSBB-Shiny/total.svg)]()

# CSBB-Shiny
Computational Suite for Bioinformaticians and Biologists (CSBB), is a RShiny application developed with an intention to empower researchers from wet and dry lab to perform downstream Bioinformatics analysis. CSBB powered by RShiny is packed with 8 modules Visualization, Normalization, Basic Stats, Differential Expression, Correlation Profiles, Function/Pathway Enrichment, ChIP-ATAC Seq and Single Cell RNA-Seq analysis. These modules are designed in order to help researchers design a hypothesis or answer research questions with little or no expertise in Bioinformatics. CSBB is also available as a command line application and has Next generation sequencing data processing capabilities. New modules and functionalities will be added periodically.

## CSBB-Shiny featured in Top 7 Shiny dashboard examples in Life Sciences by R-bloggers 
[https://www.r-bloggers.com/2022/03/r-shiny-in-life-sciences-top-7-dashboard-examples/]

# Run CSBB-Shiny from your R console or RStudio with one command
Open your R console or RStudio and paste the commands provided below. CSBB-Shiny will automatically install all required dependencies (R packages).

```
library(shiny)
runGitHub("CSBB-Shiny", "praneet1988", launch.browser = TRUE)

```
In case R package shiny is not installed please run the following command.

```
install.packages("shiny", repos="http://cran.us.r-project.org")
library(shiny)

```

# End to End Pipeline Workflow with CSBB
![Graph](CSBB.png)

# Video Tutorial
https://youtu.be/c0P7TMu_IyY

# Single Cell Transcriptomics Analysis Tutorial
https://www.youtube.com/watch?v=s8Q4o1e-f1E

# CSBB Single Cell Browser Tutorial
https://www.dropbox.com/s/9po853gc5gzdxnl/CSBB-Shiny_scRNA_Browser.mov?dl=0

