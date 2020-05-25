options(repos = BiocManager::repositories())
library(shiny)
library(servr)
library(ggplot2)
library(pheatmap)
library(M3C)
library(RUVSeq)
library(scales)
library(dtwclust)
library(dplyr)
library(ggcorrplot)
library(tibble)
library(ReactomePA)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(EnhancedVolcano)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(clusterProfiler)

options(shiny.maxRequestSize=100*1024^2)
ui <- fluidPage(
  
  tags$head(includeHTML(("GoogleAnalytics.html"))),
  # App title ----
  titlePanel("Computational Suite For Bioinformaticians and Biologists"),

  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of Dataset ----
      selectInput("Module",
                  label = "Module",
                  choices = c("Normalization", "Basic Stats", "Visualization", "Differential Expression", "Correlation Profiles", "Functional and Pathway Enrichment", "ChIP-ATAC Seq Analysis"),
                  selected = "Visualization"),

      conditionalPanel(
        condition = "input.Module == 'Normalization' || input.Module == 'Visualization' || input.Module == 'Basic Stats' || input.Module == 'Correlation Profiles'",
        fileInput("File",
                  label = "Upload Expression File (Accepted Format: tab delimited text)",
                  accept = c("text", "text", ".txt"))),

      conditionalPanel(
        condition = "input.Module == 'Functional and Pathway Enrichment'",
        fileInput("GeneList_FP",
                  label = "Upload a Gene List (Accepted Format: .txt)",
                  accept = c("text", "text", ".txt"))),

      conditionalPanel(
        condition = "input.Module == 'Differential Expression'",
        fileInput("Counts",
                  label = "Upload Counts File (Accepted Format: tab delimited text)",
                  accept = c("text", "text", ".txt"))),

      conditionalPanel(
        condition = "input.Module == 'Normalization'",
        selectInput("Norm",
                    label = "Choose Normalization",
                    choices = c("upper quantile", "median", "full", "log2", "zScore", "none"),
                    selected = "upper quantile")),

      conditionalPanel(
        condition = "input.Module == 'Differential Expression'",
        selectInput("Controls", 
                    label = "Select Control Samples",
                    choices = NULL,
                    selected = NULL,
                    multiple = TRUE)),

      conditionalPanel(
        condition = "input.Module == 'Differential Expression'",
        selectInput("Treatments", 
                    label = "Select Treatment Samples",
                    choices = NULL,
                    selected = NULL,
                    multiple = TRUE)),

      conditionalPanel(
        condition = "input.Module == 'Differential Expression'",
        numericInput("DEFilterLog", 
                    label = "Log2 Fold Change Cutoff",
                    value = 0.5)),

      conditionalPanel(
        condition = "input.Module == 'Differential Expression'",
        numericInput("DEFilterFDR", 
                    label = "False Discovery Rate (FDR) Cutoff",
                    value = 0.1)),
      
      conditionalPanel(
        condition = "input.Module == 'Correlation Profiles'",
        selectInput("Correlation", 
                    label = "Get Correlation among",
                    choices = c("Genes","Samples"),
                    selected = "Samples",
                    multiple = FALSE)),

      conditionalPanel(
        condition = "input.Module == 'Correlation Profiles' && input.Correlation == 'Genes'",
        fileInput("GeneList",
                  label = "Upload Gene List",
                  accept = c("text", "text", ".txt"))),
      
      conditionalPanel(
        condition = "input.Module == 'Correlation Profiles'",
        selectInput("CorrelationMethod",
                  label = "Select Correlation Method",
                  choices = c("pearson", "kendall", "spearman"),
                  selected = "pearson")),

      conditionalPanel(
        condition = "input.Module == 'Normalization' || input.Module == 'Visualization' || input.Module == 'Differential Expression'",
        selectInput("PlotType", 
                    label = "Visualization",
                    choices = c("pca", "tsne", "heatmap"),
                    selected = "pca")),

      conditionalPanel(
        condition = "input.PlotType == 'tsne'",
        sliderInput("perplexity", 
                    label = "Choose a perplexity value (If number of samples <= 50 : recommended value 10",
                    min = 2,
                    max = 30,
                    value = 10)),

      conditionalPanel(
        condition = "input.PlotType == 'tsne'",
        fileInput("GroupFile",
                  label = "Upload a Sample Group File",
                  accept = c("text", "text", ".txt"))),

      conditionalPanel(
        condition = "input.PlotType == 'heatmap'",
        selectInput("Cluster", 
                    label = "Perform Clustering on",
                    choices = c("Rows", "Columns", "Rows and Columns", "None"),
                    selected = "Rows and Columns")),

      conditionalPanel(
        condition = "input.PlotType == 'heatmap'",
        selectInput("Scaling", 
                    label = "Perform Scaling on",
                    choices = c("row", "column", "none"),
                    selected = "row")),

      conditionalPanel(
        condition = "input.PlotType == 'Functional and Pathway Enrichment'",
        selectInput("SpeciesUse", 
                    label = "Choose Species",
                    choices = c("human", "mouse"),
                    selected = "human")),

      conditionalPanel(
        condition = "input.Module == 'Functional and Pathway Enrichment'",
        selectInput("Species", 
                    label = "Choose Species",
                    choices = c("human", "mouse"),
                    selected = "human")),

      conditionalPanel(
        condition = "input.Module == 'ChIP-ATAC Seq Analysis'",
        fileInput("PeakFile",
                  label = "Upload Peak File (Peak file in .narrowPeak, .broadPeak or .bed  and in .gz format are accepted",
                  accept = c(".narrowPeak", ".broadPeak", ".bed", ".narrowPeak.gz", ".broadPeak.gz", ".bed.gz"))),

      conditionalPanel(
        condition = "input.Module == 'ChIP-ATAC Seq Analysis'",
        selectInput("SpeciesChIP", 
                    label = "Choose Species and genome version",
                    choices = c("human (hg19)", "human (hg38)", "mouse (mm10)", "mouse (mm9)"),
                    selected = "human (hg19)")),

      conditionalPanel(
        condition = "input.Module == 'ChIP-ATAC Seq Analysis'",
        selectInput("PlotChIP", 
                    label = "Select Analysis for Peaks",
                    choices = c("Coverage Plot (Visualize coverage of peaks across chromosomes)", "Heatmap showing profile of peaks binding to TSS regions", "Average Profile of peaks binding to TSS regions", "Peak Annotation (finding closest genes, postion of Peaks wrt TSS etc.)"),
                    selected = "Peak Annotation (finding closest genes, postion of Peaks wrt TSS etc.)")),

      conditionalPanel(
        condition = "input.Module == 'ChIP-ATAC Seq Analysis' && input.PlotChIP == 'Peak Annotation (finding closest genes, postion of Peaks wrt TSS etc.)'",
        selectInput("PlotAnnotation", 
                    label = "Select Visualization",
                    choices = c("Peak genomic annotation", "Functional enrichemnt"),
                    selected = "Peak genomic annotation"))
    ),
    mainPanel(
          tabsetPanel(type = "tabs",
              tabPanel("About", fluidRow(
                p(strong("Computational Suite for Bioinformaticians and Biologists (CSBB)"), "is a RShiny application developed with an intention to empower researchers from wet and dry lab to perform downstream Bioinformatics analysis. CSBB powered by RShiny is packed with 7 modules", strong("Visualization, Normalization, Basic Stats, Differential Expression, Correlation Profiles, Function/Pathway Enrichment and ChIP/ATAC Seq analysis"), ". These modules are designed in order to help researchers design a hypothesis or answer research questions with little or no expertise in Bioinformatics. CSBB is also available as a command line application and has Next generation sequencing data processing capabilities. New modules and functionalities will be added periodically.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p("CSBB RShiny is avaibale on", a("GitHub", href = "https://github.com/praneet1988/CSBB-Shiny", target = "_blank"), ", if interested in hosting on your own servers.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p("Please post issues, suggestions and improvements using", a("Issues/suggestions", href = "https://github.com/praneet1988/CSBB-Shiny", target = "_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p("To use CSBB command line application please access", a("CSBB CMD", href = "https://github.com/praneet1988/Computational-Suite-For-Bioinformaticians-and-Biologists", target = "_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p("If using CSBB RShiny in your research please cite the GitHub page", a("Cite", href = "https://github.com/praneet1988/CSBB-Shiny", target = "_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p("Developed and maintained by Praneet Chaturvedi. To view other tools and contributions please visit", a("GitHub", href = "https://github.com/praneet1988/", target = "_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px")), imageOutput('Pipeline')),
              tabPanel("Result Window", textOutput('DisplayText'), DT::dataTableOutput("result"), downloadButton('downloadResult', 'Download Results')), 
              tabPanel("Visualization Window", textOutput('DisplayText1'), downloadButton('downloadPlot', 'Save Plot'), plotOutput("Plot")),
              tabPanel("Getting Started", fluidRow(
                p(strong("CSBB"), "is easy to use and is packed with some very powerful modules to help you analyze your data. Results generated from the modules are loaded on the Result window whereas the Visualization plots are displyed on Visualization windows. Now let's see what each module helps with and what are the options users can explore.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Tutorial Video"), a("YouTube", href = "https://youtu.be/c0P7TMu_IyY", target = "_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Normalization module"), "can help users perform normalization on their data using following methods: upper quantile, median, full, log2 and zScore. Normalized data can also visualized using Principal component analysis (pca), t-stochastic neighbor embedding (tSNE) and heatmap. For tSNE visualization a group file is required. Group file should provide group name for each sample in the data. PCA is linear dimension reduction technique and tSNE is non-linear dimension reduction technique.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Visualization module"), "lets user visualize their data using pca, tSNE and heatmap", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Basic Stats module"), "is very helpful for estimating mean, median, standard deviation, median adjusted deviation, sum, min and max expression per gene in the expression matrix.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Differential Expression module"), "helps user perform differential expression analysis on raw counts of genes across samples using RUVSeq. Please cite RUVSeq if using Differential Expression module in your research using", a("Cite", href = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4404308/", target = "_blank"), ". Differentially expressed (DE) genes are reported in result tab and results can be filtered using logFC and FDR filters. Users can visualize their data based on DE genes using pca, tSNE, heatmap, volcano plots and perform functional/pathway enrichemnt on DE genes", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Correlation Profiles module"), "is developed to help users analyze Correlation among the samples in the data or see how a gene set is correlated based on the expression across samples", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Functional/Pathway Enrichment module"), "is designed to compute and visualize top enriched functions and pathways based on user provided gene list. ReactomePA R package is used to perform and visualize enrichment. Please cite ReactomePA when using the module in your research", a("Cite", href = "https://pubs.rsc.org/en/content/articlelanding/2016/MB/C5MB00663E#!divAbstract", target = "_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("ChIP-ATAC Seq Analysis module"), "is designed to perform downstream analysis with peaks like obtaining coverage plot across chromosomes, quantifying profile of peaks binding to transcription start site (TSS) and performing peak annotation. TSS is defined as +- 3kb regions of flanking sequences of the TSS sites. ChIPseeker package in R is used for performing downsream analysis on peaks. Please cite ChIPseeker when using this module.", a("Cite", href = "https://academic.oup.com/bioinformatics/article/31/14/2382/255379", target = "_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p("Please access sample files for each module using", a("Sample Files", href = "https://github.com/praneet1988/CSBB-Shiny", target = "_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"))),
			        tabPanel("What's New in CSBB Shiny", fluidRow(
			  	      p(strong("Version 1.1 Log:"), "Brushed off some known bugs which entered the system and added new Module ChIP-ATAC Seq Analysis for users. CSBB Shiny will be updated periodically. Some known bugs kicked out include result files and plots file naming, threshold changes, text changes", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Version 1.0 Log:"), "CSBB powered by Shiny has 6 modules and is built with an idea to empower researchers in performing bioinformatics analysis with little or no bioinformatics expertise", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px")))
          )
      )
  )
)
 

server <- function(input, output, session) {

 output$Pipeline <- renderImage({
  list(src = 'www/CSBB.png',
         contentType = 'image/png',
         width = 800,
         height = 800,
         alt = "This is alternate text")
  }, deleteFile = F)
 
 observe({
     if(input$Module == "Differential Expression"){
       inFile <- input$Counts
       if (is.null(inFile))
        return(NULL)
       data <- as.matrix(read.table(inFile$datapath, sep="\t", header=T, row.names=1, check.names=TRUE))
       samples <- colnames(data)
       samples <- data.frame(samples)
       samples <- samples$samples
       updateSelectInput(session, "Controls", 
                         label = "Select Control Samples",
                         choices = samples)
    }
 })

 observe({
     if(input$Module == "Differential Expression"){
       inFile <- input$Counts
       if (is.null(inFile))
        return(NULL)
       data <- as.matrix(read.table(inFile$datapath, sep="\t", header=T, row.names=1, check.names=TRUE))
       samples <- colnames(data)
       samples <- data.frame(samples)
       samples <- samples$samples
       updateSelectInput(session, "Treatments", 
                         label = "Select Treatment Samples",
                         choices = samples)
    }
 })

 observe({
     if(input$Module == "Differential Expression"){
       updateSelectInput(session, "PlotType", 
                         label = "Visualization",
                         choices = c("pca", "tsne", "heatmap", "Functional and Pathway Enrichment", "Volcano Plots"),
                         selected = "heatmap")
    }
 })

 NormalizationData <- reactive({
    if(input$Module == "Normalization"){
      inFile <- input$File
      if (is.null(inFile))
      	return(NULL)
      if(input$Norm == "upper quantile"){
        source("app_bin/UQ_Norm.r", local = TRUE)
        return(UData)
      }
      else if(input$Norm == "median"){
        source("app_bin/Median_Norm.r", local = TRUE)
        return(UData)
      }
      else if(input$Norm == "full"){
        source("app_bin/Full_Norm.r", local = TRUE)
        return(UData)
      }
      else if(input$Norm == "none"){
        source("app_bin/No_Norm.r", local = TRUE)
        return(UData)
      }
      else if(input$Norm == "log2"){
        source("app_bin/log2_Norm.r", local = TRUE)
        return(UData)
      }
      else if(input$Norm == "zScore"){
        source("app_bin/Zscore_Norm.r", local = TRUE)
        return(UData)
      }
    }
})

BasicStatsData <- reactive({
    if(input$Module == "Basic Stats"){
      inFile <- input$File
      if (is.null(inFile))
        return(NULL)
      source("app_bin/BasicStats.r", local = TRUE)
      return(output)
    }
})

VisualizationData <- reactive({
    if(input$Module == "Visualization"){
      inFile <- input$File
      if (is.null(inFile))
      	return(NULL)
      data <- read.table(inFile$datapath, sep="\t", header=T, row.names=1, check.names=F)
      return(data)
    }
})

DEData <- reactive({
    if(input$Module == "Differential Expression"){
      inFile <- input$Counts
      if (is.null(inFile))
      	return(NULL)
      data <- as.matrix(read.table(inFile$datapath, sep="\t", header=T, row.names=1, check.names=TRUE))
      data <- data.frame(data)
      samplelist <- c(input$Controls, input$Treatments)
      data_temp <- data[,samplelist]
      data_temp <- as.matrix(data_temp)
      lengthcontrol <- length(input$Controls)
      lengthtreatment <- length(input$Treatments)
      if((lengthcontrol == 1)&(lengthtreatment == 1)){
        source("app_bin/RUVseq_NoReps.r", local = TRUE, echo=FALSE)
        DEresult_filter <-  DEresult %>% rownames_to_column('gene') %>% filter(logFC >= input$DEFilterLog | logFC <= -1*(input$DEFilterLog), FDR <= input$DEFilterFDR) %>% column_to_rownames('gene')
        if(is.null(DEresult_filter))
        	return(NULL)
        return(DEresult_filter)
      }
      else if((lengthcontrol > 1)&(lengthtreatment >= 1)){
        cutoff <- round((lengthcontrol + lengthtreatment)/2)
        source("app_bin/RUVseq_replicates.r", local = TRUE, echo=FALSE)
        DEresult_filter <-  DEresult %>% rownames_to_column('gene') %>% filter(logFC >= input$DEFilterLog | logFC <= -1*(input$DEFilterLog), FDR <= input$DEFilterFDR) %>% column_to_rownames('gene')
        if(is.null(DEresult_filter))
        	return(NULL)
        return(DEresult_filter)
      }
      else if((lengthcontrol >= 1)&(lengthtreatment > 1)){
        cutoff <- round((lengthcontrol + lengthtreatment)/2)
        source("app_bin/RUVseq_replicates.r", local = TRUE, echo=FALSE)
        DEresult_filter <-  DEresult %>% rownames_to_column('gene') %>% filter(logFC >= input$DEFilterLog | logFC <= -1*(input$DEFilterLog), FDR <= input$DEFilterFDR) %>% column_to_rownames('gene')
        if(is.null(DEresult_filter))
        	return(NULL)
        return(DEresult_filter)
      }
      else if((lengthcontrol == 0)&(lengthtreatment == 0)){
      	return(NULL)
      }
    }
})

CorrelationData <- reactive({
      if(input$Module == "Correlation Profiles"){
        if(input$Correlation == "Samples"){
        inFile <- input$File
        if (is.null(inFile))
        	return(NULL)
        data <- as.matrix(read.table(inFile$datapath, sep="\t", header=T, row.names=1, check.names=F))
        cormat <- cor(data, method=input$CorrelationMethod, use = "na.or.complete")
        return(cormat)
      }
      else if(input$Correlation == "Genes"){
        inFile <- input$File
        if (is.null(inFile))
        	return(NULL)
        data <- as.matrix(read.table(inFile$datapath, sep="\t", header=T, row.names=1, check.names=F))
        Genes_List <- input$GeneList
        if (is.null(Genes_List))
        	return(NULL)
        GenesUpload <- as.matrix(read.table(Genes_List$datapath, sep="\n", header=F))
        GenesUpload <- data.frame(GenesUpload)
        GenesUpload <- unique(GenesUpload$V1)
        datause <- subset(data, rownames(data) %in% GenesUpload)
        datause <- t(datause)
        cormat <- cor(datause, method=input$CorrelationMethod, use = "na.or.complete")
        return(cormat)
      }
    }
})

FPEnrichmentData <- reactive({
      if(input$Module == "Functional and Pathway Enrichment"){
        if(input$Species == "human"){
          inFile <- input$GeneList_FP
          if (is.null(inFile))
          	return(NULL)
          GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
          GenesUpload <- data.frame(GenesUpload)
          GenesUpload <- unique(GenesUpload$V1)
          genes_ids <- mapIds(org.Hs.eg.db, as.character(GenesUpload), 'ENTREZID', 'SYMBOL')
          enrichemnt <- enrichPathway(gene = genes_ids, pvalueCutoff = 0.05, readable=T, organism = "human")
          return(enrichemnt)
        }
        else if(input$Species == "mouse"){
          inFile <- input$GeneList_FP
          if (is.null(inFile))
          	return(NULL)
          GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
          GenesUpload <- data.frame(GenesUpload)
          GenesUpload <- unique(GenesUpload$V1)
          genes_ids <- mapIds(org.Mm.eg.db, as.character(GenesUpload), 'ENTREZID', 'SYMBOL')
          enrichemnt <- enrichPathway(gene = genes_ids, pvalueCutoff = 0.05, readable=T, organism = "mouse")
          return(enrichemnt)
        }
      }
})

CHIPData <- reactive({
      if(input$Module == "ChIP-ATAC Seq Analysis"){
        if(input$SpeciesChIP == "human (hg19)"){
          if((input$PlotChIP == "Coverage Plot (Visualize coverage of peaks across chromosomes)")|(input$PlotChIP == "Heatmap showing profile of peaks binding to TSS regions")|(input$PlotChIP == "Average Profile of peaks binding to TSS regions")){
          	inFile <- input$PeakFile
            if (is.null(inFile))
            	return(NULL)
            peak <- readPeakFile(inFile$datapath)
            return(peak)
          }
          else if(input$PlotChIP == "Peak Annotation (finding closest genes, postion of Peaks wrt TSS etc.)"){
            inFile <- input$PeakFile
            if (is.null(inFile))
            	return(NULL)
            peak <- readPeakFile(inFile$datapath)
            txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
            if(input$PlotAnnotation == "Peak genomic annotation"){
            	peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
            	return(peakAnno)
            }
            else if(input$PlotAnnotation == "Functional enrichemnt"){
            	peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
            	pathway <- enrichPathway(as.data.frame(peakAnno)$geneId, pvalueCutoff = 0.05, readable=T, organism = "human")
              	return(pathway)
            }
          }
        }
        else if(input$SpeciesChIP == "human (hg38)"){
          if((input$PlotChIP == "Coverage Plot (Visualize coverage of peaks across chromosomes)")|(input$PlotChIP == "Heatmap showing profile of peaks binding to TSS regions")|(input$PlotChIP == "Average Profile of peaks binding to TSS regions")){
            inFile <- input$PeakFile
            if (is.null(inFile))
            	return(NULL)
            peak <- readPeakFile(inFile$datapath)
            return(peak)
          }
          else if(input$PlotChIP == "Peak Annotation (finding closest genes, postion of Peaks wrt TSS etc.)"){
            inFile <- input$PeakFile
            if (is.null(inFile))
            	return(NULL)
            peak <- readPeakFile(inFile$datapath)
            txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
            if(input$PlotAnnotation == "Peak genomic annotation"){
            	peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
            	return(peakAnno)
            }
            else if(input$PlotAnnotation == "Functional enrichemnt"){
            	peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
            	pathway <- enrichPathway(as.data.frame(peakAnno)$geneId, pvalueCutoff = 0.05, readable=T, organism = "human")
              	return(pathway)
            }
          }
        }
        else if(input$SpeciesChIP == "mouse (mm10)"){
          if((input$PlotChIP == "Coverage Plot (Visualize coverage of peaks across chromosomes)")|(input$PlotChIP == "Heatmap showing profile of peaks binding to TSS regions")|(input$PlotChIP == "Average Profile of peaks binding to TSS regions")){
          	inFile <- input$PeakFile
            if (is.null(inFile))
            	return(NULL)
            peak <- readPeakFile(inFile$datapath)
            return(peak)
          }
          else if(input$PlotChIP == "Peak Annotation (finding closest genes, postion of Peaks wrt TSS etc.)"){
            inFile <- input$PeakFile
            if (is.null(inFile))
            	return(NULL)
            peak <- readPeakFile(inFile$datapath)
            txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
            if(input$PlotAnnotation == "Peak genomic annotation"){
            	peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
            	return(peakAnno)
            }
            else if(input$PlotAnnotation == "Functional enrichemnt"){
            	peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
            	pathway <- enrichPathway(as.data.frame(peakAnno)$geneId, pvalueCutoff = 0.05, readable=T, organism = "mouse")
            	return(pathway)
            }
          }
        }
        else if(input$SpeciesChIP == "mouse (mm9)"){
          if((input$PlotChIP == "Coverage Plot (Visualize coverage of peaks across chromosomes)")|(input$PlotChIP == "Heatmap showing profile of peaks binding to TSS regions")|(input$PlotChIP == "Average Profile of peaks binding to TSS regions")){
            inFile <- input$PeakFile
            if (is.null(inFile))
            	return(NULL)
            peak <- readPeakFile(inFile$datapath)
            return(peak)
          }
          else if(input$PlotChIP == "Peak Annotation (finding closest genes, postion of Peaks wrt TSS etc.)"){
            inFile <- input$PeakFile
            if (is.null(inFile))
            	return(NULL)
            peak <- readPeakFile(inFile$datapath)
            txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
            if(input$PlotAnnotation == "Peak genomic annotation"){
            	peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
            	return(peakAnno)
            }
            else if(input$PlotAnnotation == "Functional enrichemnt"){
            	peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
            	pathway <- enrichPathway(as.data.frame(peakAnno)$geneId, pvalueCutoff = 0.05, readable=T, organism = "mouse")
              	return(pathway)
            }
          }
        }
      }
})

 PCAplot <- reactive({
     inFile <- input$File
     if (is.null(inFile))
     	return(NULL)
     data <- c()
     if(input$Module == "Normalization"){
     	data <- as.matrix(NormalizationData())
     }
     else if(input$Module == "Visualization"){
     	data <- as.matrix(VisualizationData())
     }
     data <- data[apply(data[,-1], 1, function(x) !all(x==0)),]
     data.t <- t(data)
     pca <- prcomp(data.t, center=T, scale. = T)
     pc1 <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
     pc2 <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
     PC1_use <- paste0("PC1", "(", pc1, "%)")
     PC2_use <- paste0("PC2", "(", pc2, "%)")
     Samples_temp <- rownames(data.t)
     Samples <- factor(Samples_temp)
     scores <- data.frame(Samples_temp, pca$x[,1:3])
     MIN_X <- min(scores$PC1)
     Max_X <- max(scores$PC1)
     header <- "Principal Component Analysis"
     qplot(x=PC1, y=PC2, data=scores, colour=Samples, xlim=c(MIN_X-75,Max_X+75)) + xlab(PC1_use) + ylab(PC2_use) + geom_point(shape=1) + geom_text(aes(label=Samples_temp), hjust=0, vjust=0) + scale_size_area() + theme(axis.text = element_text(size = 14),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"),legend.key = element_rect(fill = "white"),legend.background = element_rect(fill = "white"),panel.grid.major = element_line(),panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white")) + ggtitle(header) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
 })

 TSNEplot <- reactive({
     inFile <- input$File
     if (is.null(inFile))
     	return(NULL)
     inGroupFile <- input$GroupFile
     if (is.null(inGroupFile))
     	return(NULL)
     data <- c()
     if(input$Module == "Normalization"){
     	data <- as.matrix(NormalizationData())
     }
     else if(input$Module == "Visualization"){
     	data <- as.matrix(VisualizationData())
     }
     Group <- as.matrix(read.table(inGroupFile$datapath, header=T, sep="\t", row.names=1, check.names=F))
     Group <- data.frame(Group)
     GroupUse <- as.factor(Group$group)
     tsne(data, labels=GroupUse, perplex=input$perplexity) + xlab("tSNE_1") + ylab("tSNE_2") + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
 })

 HEATMAPplot <- reactive({
     inFile <- input$File
     if (is.null(inFile))
     	return(NULL)
     data <- c()
     if(input$Module == "Normalization"){
     	data <- as.matrix(NormalizationData())
     }
     else if(input$Module == "Visualization"){
     	data <- as.matrix(VisualizationData())
     }
     data <- data[apply(data, MARGIN = 1, FUN = function(x) sd(x) != 0),]
     if(input$Cluster == "Rows"){
     	pheatmap(data, scale=input$Scaling, cluster_rows=TRUE, cluster_cols=FALSE, main="Heatmap", border_color = "NA")
     }
     else if(input$Cluster == "Columns"){
      	pheatmap(data, scale=input$Scaling, cluster_rows=FALSE, cluster_cols=TRUE, main="Heatmap", border_color = "NA")
     }
     else if(input$Cluster == "Rows and Columns"){
      	pheatmap(data, scale=input$Scaling, cluster_rows=TRUE, cluster_cols=TRUE, main="Heatmap", border_color = "NA")
     }
     else if(input$Cluster == "None"){
      	pheatmap(data, scale=input$Scaling, cluster_rows=FALSE, cluster_cols=FALSE, main="Heatmap", border_color = "NA")
     }
 })

 CorrelationPlot <- reactive({
     if(input$Correlation == "Samples"){
     	inFile <- input$File
      	if (is.null(inFile))
      		return(NULL)
      	titleuse <- paste0("Displaying Correlation Plot of ", input$Correlation)
      	data <- CorrelationData()
      	ggcorrplot(data, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation") + labs(title = titleuse)
     }
     else if(input$Correlation == "Genes"){
      	inFile <- input$File
      	if (is.null(inFile))
        	return(NULL)
      	Genes_List <- input$GeneList
      	if (is.null(Genes_List))
        	return(NULL)
      	titleuse <- paste0("Displaying Correlation Plot of ", input$Correlation)
      	data <- 	CorrelationData()
      	ggcorrplot(data, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation") + labs(title = titleuse)
     }
 })

 DifferentialExpressionPlot <- reactive({
     inFile <- input$Counts
     if (is.null(inFile))
     	return(NULL)
     DEresult <- DEData()
     DEgenes <- rownames(DEresult)
     data <- as.matrix(read.table(inFile$datapath, sep="\t", header=T, row.names=1, check.names=F))
     datause <- subset(data, rownames(data) %in% DEgenes)
     datause <- log2(datause+1)
     if(input$PlotType == "pca"){
     	data <- as.matrix(datause)
      	data <- data[apply(data[,-1], 1, function(x) !all(x==0)),]
      	data.t <- t(data)
      	pca <- prcomp(data.t, center=T, scale. = T)
      	pc1 <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
     	pc2 <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
     	PC1_use <- paste0("PC1", "(", pc1, "%)")
      	PC2_use <- paste0("PC2", "(", pc2, "%)")
      	Samples_temp <- rownames(data.t)
      	Samples <- factor(Samples_temp)
      	scores <- data.frame(Samples_temp, pca$x[,1:3])
      	MIN_X <- min(scores$PC1)
      	Max_X <- max(scores$PC1)
      	header <- "Principal Component Analysis"
      	qplot(x=PC1, y=PC2, data=scores, colour=Samples, xlim=c(MIN_X-75,Max_X+75)) + xlab(PC1_use) + ylab(PC2_use) + geom_point(shape=1) + geom_text(aes(label=Samples_temp), hjust=0, vjust=0) + scale_size_area() + theme(axis.text = element_text(size = 14),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"),legend.key = element_rect(fill = "white"),legend.background = element_rect(fill = "white"),panel.grid.major = element_line(),panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white")) + ggtitle(header) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
    }
    else if(input$PlotType == "tsne"){
      inGroupFile <- input$GroupFile
      if (is.null(inGroupFile))
      	return(NULL)
      data <- as.matrix(datause)
      Group <- as.matrix(read.table(inGroupFile$datapath, header=T, sep="\t", row.names=1, check.names=F))
      Group <- data.frame(Group)
      GroupUse <- as.factor(Group$group)
      tsne(data, labels=GroupUse, perplex=input$perplexity) + xlab("tSNE_1") + ylab("tSNE_2") + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
    }
    else if(input$PlotType == "heatmap"){
      data <- as.matrix(datause)
      data <- data[apply(data, MARGIN = 1, FUN = function(x) sd(x) != 0),]
      if(input$Cluster == "Rows"){
      	pheatmap(data, scale=input$Scaling, cluster_rows=TRUE, cluster_cols=FALSE, main="Heatmap", border_color = "NA")
      }
      else if(input$Cluster == "Columns"){
       pheatmap(data, scale=input$Scaling, cluster_rows=FALSE, cluster_cols=TRUE, main="Heatmap", border_color = "NA")
      }
      else if(input$Cluster == "Rows and Columns"){
       pheatmap(data, scale=input$Scaling, cluster_rows=TRUE, cluster_cols=TRUE, main="Heatmap", border_color = "NA")
      }
      else if(input$Cluster == "None"){
       pheatmap(data, scale=input$Scaling, cluster_rows=FALSE, cluster_cols=FALSE, main="Heatmap", border_color = "NA")
      }
    }
    else if(input$PlotType == "Functional and Pathway Enrichment"){
      if(input$SpeciesUse == "human"){
       data <- as.matrix(datause)
       genes_de <- rownames(data)
       genes_de_id <- mapIds(org.Hs.eg.db, genes_de, 'ENTREZID', 'SYMBOL')
       enrichemnt <- enrichPathway(gene = genes_de_id, pvalueCutoff = 0.05, readable=T, organism = "human", maxGSSize = 5000)
       emapplot(enrichemnt)
      }
      else if(input$SpeciesUse == "mouse"){
       data <- as.matrix(datause)
       genes_de <- rownames(data)
       genes_de_id <- mapIds(org.Mm.eg.db, genes_de, 'ENTREZID', 'SYMBOL')
       enrichemnt <- enrichPathway(gene = genes_de_id, pvalueCutoff = 0.05, readable=T, organism = "mouse", maxGSSize = 5000)
       emapplot(enrichemnt)
      }
    }
    else if(input$PlotType == "Volcano Plots"){
       data <- DEresult
       EnhancedVolcano(data, lab = rownames(data), x = 'logFC', y = 'FDR')
    }
 })

 FPEnrichmentPlot <- reactive({
 	inFile <- input$GeneList_FP
    if (is.null(inFile))
    	return(NULL)
    enrichemnt <- FPEnrichmentData()
    emapplot(enrichemnt)
 })

 CHIPSeqPlot <- reactive({
    inFile <- input$PeakFile
    if (is.null(inFile))
    	return(NULL)
    peak <- readPeakFile(inFile$datapath)
    if(input$SpeciesChIP == "human (hg19)"){
    	txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
      	if(input$PlotChIP == "Coverage Plot (Visualize coverage of peaks across chromosomes)"){
        	covplot(peak, weightCol="V5")
      	}
      	else if(input$PlotChIP == "Heatmap showing profile of peaks binding to TSS regions"){
        	peakHeatmap(peak, TxDb=txdb, upstream=3000, downstream=3000, color="red")
      	}
      	else if(input$PlotChIP == "Average Profile of peaks binding to TSS regions"){
        	plotAvgProf2(peak, TxDb=txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
      	}
      	else if(input$PlotChIP == "Peak Annotation (finding closest genes, postion of Peaks wrt TSS etc.)"){
        	if(input$PlotAnnotation == "Peak genomic annotation"){
          		plotAnnoPie(CHIPData())
        	}
        	else if(input$PlotAnnotation == "Functional enrichemnt"){
          		emapplot(CHIPData())
        	}
      	}
    }
    else if(input$SpeciesChIP == "human (hg38)"){
    	txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
      	if(input$PlotChIP == "Coverage Plot (Visualize coverage of peaks across chromosomes)"){
        	covplot(peak, weightCol="V5")
      	}
      	else if(input$PlotChIP == "Heatmap showing profile of peaks binding to TSS regions"){
        	peakHeatmap(peak, TxDb=txdb, upstream=3000, downstream=3000, color="red")
      	}
      	else if(input$PlotChIP == "Average Profile of peaks binding to TSS regions"){
        	plotAvgProf2(peak, TxDb=txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
      	}
      	else if(input$PlotChIP == "Peak Annotation (finding closest genes, postion of Peaks wrt TSS etc.)"){
        	if(input$PlotAnnotation == "Peak genomic annotation"){
          		plotAnnoPie(CHIPData())
        	}
        	else if(input$PlotAnnotation == "Functional enrichemnt"){
          		emapplot(CHIPData())
        	}
     	}
    }
    else if(input$SpeciesChIP == "mouse (mm10)"){
    	txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
      	if(input$PlotChIP == "Coverage Plot (Visualize coverage of peaks across chromosomes)"){
        	covplot(peak, weightCol="V5")
      	}
      	else if(input$PlotChIP == "Heatmap showing profile of peaks binding to TSS regions"){
        	peakHeatmap(peak, TxDb=txdb, upstream=3000, downstream=3000, color="red")
      	}
      	else if(input$PlotChIP == "Average Profile of peaks binding to TSS regions"){
        	plotAvgProf2(peak, TxDb=txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
      	}
      	else if(input$PlotChIP == "Peak Annotation (finding closest genes, postion of Peaks wrt TSS etc.)"){
        	if(input$PlotAnnotation == "Peak genomic annotation"){
          		plotAnnoPie(CHIPData())
        	}
        	else if(input$PlotAnnotation == "Functional enrichemnt"){
          		emapplot(CHIPData())
        	}
      	}
    }
    else if(input$SpeciesChIP == "mouse (mm9)"){
    	txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
      	if(input$PlotChIP == "Coverage Plot (Visualize coverage of peaks across chromosomes)"){
        	covplot(peak, weightCol="V5")
      	}
      	else if(input$PlotChIP == "Heatmap showing profile of peaks binding to TSS regions"){
        	peakHeatmap(peak, TxDb=txdb, upstream=3000, downstream=3000, color="red")
      	}
      	else if(input$PlotChIP == "Average Profile of peaks binding to TSS regions"){
        	plotAvgProf2(peak, TxDb=txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
      	}
      	else if(input$PlotChIP == "Peak Annotation (finding closest genes, postion of Peaks wrt TSS etc.)"){
        	if(input$PlotAnnotation == "Peak genomic annotation"){
          		plotAnnoPie(CHIPData())
        	}
        	else if(input$PlotAnnotation == "Functional enrichemnt"){
          		emapplot(CHIPData())
        	}
      	}
    }
 })

 output$result <- DT::renderDataTable({
    if(input$Module == "Normalization") {
    	DT::datatable(NormalizationData())
    }
    else if(input$Module == "Visualization") {
    	DT::datatable(VisualizationData())
    }
    else if(input$Module == "Correlation Profiles") {
    	DT::datatable(CorrelationData())
    }
    else if(input$Module == "Differential Expression") {
    	DT::datatable(DEData())
    }
    else if(input$Module == "Functional and Pathway Enrichment"){
    	DT::datatable(data.frame(FPEnrichmentData()))
    }
    else if(input$Module == "ChIP-ATAC Seq Analysis"){
      	DT::datatable(data.frame(CHIPData()))
    }
    else if(input$Module == "Basic Stats"){
      	DT::datatable(BasicStatsData())
    }
 })

 output$Plot <- renderPlot({
    if(input$Module == "Normalization") {
      if(input$PlotType == "pca") {
        PCAplot()
      }
      else if(input$PlotType == "tsne") {
        TSNEplot()
      }
      else if(input$PlotType == "heatmap") {
        HEATMAPplot()
      }
    }
    else if(input$Module == "Visualization") {
      if(input$PlotType == "pca") {
        PCAplot()
      }
      else if(input$PlotType == "tsne") {
        TSNEplot()
      }
      else if(input$PlotType == "heatmap") {
        HEATMAPplot()
      }
    }
    else if(input$Module == "Correlation Profiles") {
      CorrelationPlot()
    }
    else if(input$Module == "Differential Expression") {
      DifferentialExpressionPlot()
    }
    else if(input$Module == "Functional and Pathway Enrichment"){
      FPEnrichmentPlot()
    }
    else if(input$Module == "ChIP-ATAC Seq Analysis"){
      CHIPSeqPlot()
    }
  }, width=1200, height=1000)

 output$downloadResult <- downloadHandler(
      filename = function() {
        paste0(input$Module, "_Analysis_Result", "-", Sys.Date(), ".txt")
      },
      content = function(file) {
        if(input$Module == "Normalization") {
          NormalizationData()
          write.table(NormalizationData(), file, sep="\t", quote=F)
        }
        else if(input$Module == "Visualization") {
          VisualizationData()
          write.table(VisualizationData(), file, sep="\t", quote=F)
        }
        else if(input$Module == "Basic Stats") {
          BasicStatsData()
          write.table(BasicStatsData(), file, sep="\t", quote=F)
        }
        else if(input$Module == "Differential Expression") {
          DEData()
          write.table(DEData(), file, sep="\t", quote=F)
        }
        else if(input$Module == "Functional and Pathway Enrichment"){
          data.frame(FPEnrichmentData())
          write.table(data.frame(FPEnrichmentData()), file, sep="\t", quote=F)
        }
        else if(input$Module == "ChIP-ATAC Seq Analysis"){
          data.frame(CHIPData())
          write.table(data.frame(CHIPData()), file, sep="\t", quote=F)
        }
        else if(input$Module == "Correlation Profiles"){
          CorrelationData()
          write.table(CorrelationData(), file, sep="\t", quote=F)
        }
      }
 )

  output$DisplayText <- renderText({"Please Upload Input Data"})
  output$DisplayText1 <- renderText({"Please Upload Input Data"})

  observe({
  	if(input$Module == "Normalization"){
  		output$DisplayText <- renderText({"Please Upload Input Data"})
  		inFile <- input$File
  		if((!is.null(inFile))&&(is.null(NormalizationData()))){
  			output$DisplayText <- renderText({"CSBB is processing your data"})
  			output$DisplayText1 <- renderText({"CSBB is generating selected visualization plot"})
  		}
  		else if((!is.null(inFile))&&(!is.null(NormalizationData()))){
  			output$DisplayText <- renderText({"Processing Complete"})
  			output$DisplayText1 <- renderText({"Plot Generated"})
  		}
  	}
  	else if(input$Module == "Visualization"){
  		output$DisplayText <- renderText({"Please Upload Input Data"})
  		inFile <- input$File
  		if((!is.null(inFile))&&(is.null(VisualizationData()))){
  			output$DisplayText <- renderText({"CSBB is processing your data"})
  			output$DisplayText1 <- renderText({"CSBB is generating selected visualization plot"})
  		}
  		else if((!is.null(inFile))&&(!is.null(VisualizationData()))){
  			output$DisplayText <- renderText({"Processing Complete"})
  			output$DisplayText1 <- renderText({"Plot Generated"})
  		}
  	}
  	else if(input$Module == "Basic Stats"){
  		output$DisplayText <- renderText({"Please Upload Input Data"})
  		inFile <- input$File
  		if((!is.null(inFile))&&(is.null(BasicStatsData()))){
  			output$DisplayText <- renderText({"CSBB is processing your data"})
  			output$DisplayText1 <- renderText({"No Visualization avaibale"})
  		}
  		else if((!is.null(inFile))&&(!is.null(BasicStatsData()))){
  			output$DisplayText <- renderText({"Processing Complete"})
  			output$DisplayText1 <- renderText({"No Visualization avaibale"})
  		}
  	}
  	else if(input$Module == "Differential Expression"){
  		output$DisplayText <- renderText({"Please Upload Input Data"})
  		inFile <- input$Counts
  		if((!is.null(inFile))&&(is.null(DEData()))){
  			output$DisplayText <- renderText({"CSBB is processing your data"})
  			output$DisplayText1 <- renderText({"CSBB is generating selected visualization plot"})
  		}
  		else if((!is.null(inFile))&&(!is.null(DEData()))){
  			output$DisplayText <- renderText({"Processing Complete"})
  			output$DisplayText1 <- renderText({"Plot Generated"})
  		}
  	}
  	else if(input$Module == "Correlation Profiles"){
  		output$DisplayText <- renderText({"Please Upload Input Data"})
  		inFile <- input$File
  		if((!is.null(inFile))&&(is.null(CorrelationData()))){
  			output$DisplayText <- renderText({"CSBB is processing your data"})
  			output$DisplayText1 <- renderText({"CSBB is generating selected visualization plot"})
  		}
  		else if((!is.null(inFile))&&(!is.null(CorrelationData()))){
  			output$DisplayText <- renderText({"Processing Complete"})
  			output$DisplayText1 <- renderText({"Plot Generated"})
  		}
  	}
  	else if(input$Module == "Functional and Pathway Enrichment"){
  		output$DisplayText <- renderText({"Please Upload Input Data"})
  		inFile <- input$GeneList_FP
  		if((!is.null(inFile))&&(is.null(FPEnrichmentData()))){
  			output$DisplayText <- renderText({"CSBB is processing your data"})
  			output$DisplayText1 <- renderText({"CSBB is generating selected visualization plot"})
  		}
  		else if((!is.null(inFile))&&(!is.null(FPEnrichmentData()))){
  			output$DisplayText <- renderText({"Processing Complete"})
  			output$DisplayText1 <- renderText({"Plot Generated"})
  		}
  	}
  	else if(input$Module == "ChIP-ATAC Seq Analysis"){
  		output$DisplayText <- renderText({"Please Upload Input Data"})
  		inFile <- input$PeakFile
  		if((!is.null(inFile))&&(is.null(CHIPData()))){
  			output$DisplayText <- renderText({"CSBB is processing your data"})
  			output$DisplayText1 <- renderText({"CSBB is generating selected visualization plot"})
  		}
  		else if((!is.null(inFile))&&(!is.null(CHIPData()))){
  			output$DisplayText <- renderText({"Processing Complete"})
  			output$DisplayText1 <- renderText({"Plot Generated"})

  		}
  	}
  })

  output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0(input$Module, "_Analysis_Plot", "-", Sys.Date(), ".png")
      },
      content = function(file) {
        if(input$Module == "Normalization") {
          if(input$PlotType == "pca") {
            ###png(file, width = 600, height = 600, res = 600)
            PCAplot()
            ggsave(file, width = 15, height = 15)
          }
          else if(input$PlotType == "tsne") {
            TSNEplot()
            ggsave(file, width = 15, height = 15)
          }
          else if(input$PlotType == "heatmap") {
            png(file, width = 1200, height = 1000)
            print(HEATMAPplot(), useSource=TRUE)
            dev.off()
          }
        }
        else if(input$Module == "Visualization") {
          if(input$PlotType == "pca") {
            PCAplot()
            ggsave(file, width = 15, height = 15)
          }
          else if(input$PlotType == "tsne") {
            TSNEplot()
            ggsave(file, width = 15, height = 15)
          }
          else if(input$PlotType == "heatmap") {
            png(file, width = 1200, height = 1000)
            print(HEATMAPplot(), useSource=TRUE)
            dev.off()
          }
       }
       else if(input$Module == "Correlation Profiles") {
          CorrelationPlot()
          ggsave(file, width = 15, height = 15)
       }
       else if(input$Module == "Differential Expression") {
          png(file, width = 1200, height = 1000)
          print(DifferentialExpressionPlot(), useSource=TRUE)
          dev.off()
       }
       else if(input$Module == "Functional and Pathway Enrichment"){
          FPEnrichmentPlot()
          ggsave(file, width = 15, height = 15)
       }
       else if(input$Module == "ChIP-ATAC Seq Analysis"){
          CHIPSeqPlot()
          ggsave(file, width = 15, height = 15)
       }
    },
    contentType = 'image/png'
  )
}
shinyApp(ui = ui, server = server)