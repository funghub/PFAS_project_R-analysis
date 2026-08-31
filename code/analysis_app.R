
# knitr::opts_knit$set(root.dir = here::here()) # left from Rmd
# knitr::opts_chunk$set(echo = TRUE)
# set the root dir as [1] "/home/lawrence/Documents/School/thesis/PFAS_R-analysis"

# Load additional libraries that might be useful (these may need to be installed)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(gProfileR)

# convert gene id to symbols
library(org.Hs.eg.db)
library(annotate)

library(dplyr) # this is to load piping %>%
library(tibble)  # for `rownames_to_column` and `column_to_rownames`
select <- dplyr::select # this is to help with the conflict of select()
count <- dplyr::count

library(DESeq2)
library(ggplot2)
library(pheatmap)

library(shiny)
library(zip)

setwd(here::here())

REF_CHEMICALS <- c("acetaminophen","aflatoxin B1","cyclosporin A","WY-14643",
                   "rosiglitazone","troglitazone","caffeine","sucrose","potassium chloride")
PFAS_CHEMICALS <- c("PFOA","PFBS","PFOS","PFDS")
CONTROLS <- c("DMSO_ref", "DMSO_pfas", "water")
DURATION_CHOICES <- c(1, 10, 4, 14)
DEFAULT_DESIGN <- "~ chemical + duration"




# CREATE master metadata + counts (everything above the old user-choice / DESEQ section)
# FIRST function
build_masters <- function() {


# CREATE the Master Metadata File (PFAS dataset + Ref Dataset)

## Retrieve the SRA Metadata csv: Reference dataset
ref_metadata <- read.csv("./data/metadata/ref/SraRunTable.csv") %>% #, row.names = 1)
  rename('SRR' = 'Run')
ref_metadata$timepoint <- as.numeric(sub("d", "", sub("24h", "1", ref_metadata$timepoint)))

SRR_table <- ref_metadata %>% select('SRR','chemical','timepoint','dose','replicate')
SRR_table_ref <- SRR_table %>%
  rename('chemical.y' = 'chemical', 'duration' = 'timepoint')

## Export the entire SRR table to csv to combine with the PFAS one in the future
dir.create("./data/data_input_deseq/metadata_ref", recursive = TRUE, showWarnings = FALSE)
write.csv(SRR_table_ref,"./data/data_input_deseq/metadata_ref/ref_metadata_all.csv",row.names = FALSE)

## Retrieve the SRA Metadata csv: PFAS dataset
PFAS_metadata <- read.csv("data/metadata/pfas/SraRunTable.csv", row.names = 1)

# SRR → GSM → title to recover concentration/dose and timepoint/duration
# Documentation of GEOquery: https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html
library(GEOquery)

gse <- getGEO("GSE144775", GSEMatrix = FALSE)
gsm_list <- GSMList(gse)

# In gsm_titles, pull out duration and concentration in the title into two columns append
gsm_titles <- data.frame(
  GSM = names(gsm_list),
  title = sapply(gsm_list, function(x) Meta(x)$title),
  chemical = sub("^chemical: ","",sapply(gsm_list, function(x) Meta(x)$characteristics_ch1[1])),
  duration = as.numeric(sub("^time point \\(days\\): ","",sapply(gsm_list, function(x) Meta(x)$characteristics_ch1[3]))), # days
  # duration_full = sapply(gsm_list, function(x) Meta(x)$characteristics_ch1[3]),
  dose = as.numeric(sub("^dose\\(um\\): ","",sapply(gsm_list, function(x) Meta(x)$characteristics_ch1[2]))), # (um)
  # dose_full = sapply(gsm_list, function(x) Meta(x)$characteristics_ch1[2]),
  stringsAsFactors = FALSE
)

# Meta(gsm_list[["GSM4295460"]])$title
# Meta(gsm_list[[5]])$characteristics_ch1

# Find the unique chemical, dose, and duration
options(scipen = 999) # tells R to strongly avoid scientific notation whenever a fixed-notation alternative is reasonably close
unique_combos <- unique(gsm_titles[, c("chemical", "dose", "duration")])
unique_combos[order(unique_combos$chemical, unique_combos$duration, unique_combos$dose), ]
# order computes the row positions in order for those columns chem dur dose, then [takes those rows, no columns entered]

# Merge PFAS_metadata\$GEO_Accession..exp. file with gsm_titles\$GSM on the equivalent GSM#s --> reveal which SRR read belongs to which GSM on metadata file
# rename PFAS_metadata\$GEO_Accession..exp. to PFAS_metadata\$GSM
PFAS_metadata <- PFAS_metadata %>%
  rename('GSM' = 'GEO_Accession..exp.') %>%
  rownames_to_column('SRR') # rownames to be a column so don't disappear after merge

# put all df into list / merge all df in list
df_GSMtoSRR_merged <- full_join(PFAS_metadata, gsm_titles, by = 'GSM')

# Validation, row count check before and after merge
c(nrow(PFAS_metadata), nrow(gsm_titles))  # rows in each individual df
nrow(df_GSMtoSRR_merged)

# Validation, find if anything in dataframe is NA
df_GSMtoSRR_merged %>% filter(if_any(everything(), is.na))
which(is.na(df_GSMtoSRR_merged))
sum(is.na(df_GSMtoSRR_merged))

# This is now a central dataframe for SRR conversion to which chemical, duration, and dose for each sample
GSMtoSRR <- df_GSMtoSRR_merged %>% select('SRR','chemical.y','duration','dose')

# Export the entire SRR table to csv to combine with the ref one in the future
# - this is so that right before running DESeq, we can add chemical, duration, and dosage to the sample name for readability
dir.create("./data/data_input_deseq/metadata_PFAS", recursive = TRUE, showWarnings = FALSE)
write.csv(GSMtoSRR,"./data/data_input_deseq/metadata_PFAS/PFAS_metadata_all.csv",row.names = FALSE)

# Combine the SRR conversion tables from PFAS and Ref into a single one!
# Append the two tables together to make a master index (currently just PFAS)
SRR_all <- rbind(select(SRR_table_ref, -'replicate'), GSMtoSRR)

# Make a master metadata dataframe adding the information for duration and dose
# - but first, need to filter out those SRR in metadata_PFAS and metadata_ref for only duration == 1,10
# - then we combine those two metadata files but need to configure that metadata_PFAS have the same columns - will need to add back duration and dose to it to continue

# we need to drop some col or enter the missing data into the metadata from the index
metadata_ref_final <- ref_metadata %>%
  select(-c('cell_line','Collection_Date','geo_loc_name_country','geo_loc_name_country_continent','geo_loc_name','replicate')) %>%
  rename('duration' = 'timepoint') %>%
  mutate(chemical = ifelse(chemical == "DMSO", "DMSO_ref", chemical))

# rename metadata_PFAS column so that it matches columns of metadata_ref
metadata_PFAS_final <- PFAS_metadata %>%
  rename('Library.Name' = 'GSM')

# put dose,timepoint into metadata_PFAS_final (from SRR_all)
# merge dose and timepoint by SRR
metadata_PFAS_final_merged <- merge(metadata_PFAS_final, GSMtoSRR, by='SRR', all=TRUE) %>%
  mutate(chemical = chemical.y) %>% # this so that we only have one chemical col (will be used by DESeq)
  select(!chemical.y) %>%
  mutate(chemical = ifelse(chemical == "DMSO", "DMSO_pfas", chemical))

# make the master metadata dataframe combine metadata_PFAS_final with metadata_ref_final
metadata_all <- rbind(metadata_ref_final, metadata_PFAS_final_merged)

# validate that the sum of the 2 df is the merged
stopifnot(nrow(metadata_all) == nrow(metadata_ref_final) + nrow(metadata_PFAS_final_merged))

# validate on terminal
cat('\nMetadata all\n')
print(head(metadata_all))
cat('\nMetadata end\n')
print(colnames(metadata_all))

cat("\nValidate Unique Chemical\n")
print(unique(metadata_all$chemical))

cat("\nValidate Unique Dose\n")
print(unique(metadata_all$dose))

cat("\nValidate Unique Duration\n")
print(unique(metadata_all$duration))

# Durations and Ref Chemical we care about
# durations <- c(1,10)
# ref <- c("acetaminophen","aflatoxin B1","cyclosporin A","WY-14643","rosiglitazone","troglitazone","caffeine","sucrose","potassium chloride")


dir.create("./data/data_input_deseq", recursive = TRUE, showWarnings = FALSE)
write.csv(metadata_all, "./data/data_input_deseq/metadata_all.csv", row.names = FALSE)



# CREATE the Master Counts File (PFAS dataset + Ref Dataset)

## Ref dataset counts
markdups_counts_ref <- read.table("./data/featureCounts/ref/markdups_counts_ref.txt", header=TRUE, sep = "\t", row.names = 1)

## Replace everything after _ and the _ with nothing in all column names
markdups_counts_ref_cleancols <- markdups_counts_ref %>%
  rename_with(~gsub("_.*", "", .x)) %>%
  select(-c("Chr", "Start", "End", "Strand", "Length"))

## PFAS dataset counts
# Lead the reference counts files for PFAS dataset
markdups_counts_pfas <- read.table("./data/featureCounts/pfas/markdups_counts_pfas.txt", header=TRUE, sep = "\t", row.names = 1)
# raw_counts_pfas <- read.table("./data/featureCounts/pfas/raw_counts_pfas.txt", header=TRUE, sep ="\t", row.names = 1)

# Replace everything after _ and the _ with nothing in all column names
markdups_counts_pfas_cleancols <- markdups_counts_pfas %>%
  rename_with(~gsub("_.*", "", .x)) %>%
  select(-c("Chr", "Start", "End", "Strand", "Length"))

# raw_counts_pfas_cleancols <- raw_counts_pfas %>%
#   rename_with(~gsub("_.*", "", .x)) %>%
#   select(-c("Chr", "Start", "End", "Strand", "Length"))

# check that the row names are the same so that means genes are the same
# then combine the columns if true
if (identical(rownames(markdups_counts_pfas_cleancols), rownames(markdups_counts_ref_cleancols))) {
  markdups_counts_all <- cbind(markdups_counts_pfas_cleancols, markdups_counts_ref_cleancols)
} else {
  stop("Row names do NOT match between markdups_counts_pfas_cleancols and markdups_counts_ref_cleancols — check gene order/IDs before combining.")
}

# align metadata rows to counts columns, then rename SRR to chemical_duration_dose
# (replicates share a combo, so make.unique appends _1, _2, ...)
metadata_all <- metadata_all[match(colnames(markdups_counts_all), metadata_all$SRR), ]
stopifnot(identical(colnames(markdups_counts_all), metadata_all$SRR))

sample_id <- make.unique(
  paste(gsub(" ", "_", metadata_all$chemical), 
        paste0(metadata_all$duration, "d"), 
        paste0(metadata_all$dose,"uM"), #\u00b5M
        sep = "-"),
  sep = "-"
)
metadata_all$sample_id <- sample_id
rownames(metadata_all) <- sample_id
colnames(markdups_counts_all) <- sample_id

print("metadata_all validate sample_id START")
print(head(metadata_all$sample_id))
print(head(metadata_all$chemical))
print(head(metadata_all$duration))
print(head(metadata_all$dose))
print("metadata_all validate sample_id END")


write.csv(metadata_all, "./data/data_input_deseq/metadata_all.csv", row.names = TRUE)
write.csv(tibble::rownames_to_column(markdups_counts_all, "gene"),
          "./data/data_input_deseq/markdups_counts_all.csv", row.names = FALSE)

# validate on terminal
cat('\nmarkdups_counts_all #cols == metadata_all #rows\n')
print(ncol(markdups_counts_all) == nrow(metadata_all))
cat('\ncomparison validation ncols==nrows end\n')

print("\n\nWARNING: data includes ref chemicals not of focus and durations 4 and 14 not of focus")
print("DMSO for ref and for pfas dataset has been separated into DMSO_ref and DMSO_pfas\n\n")

print("Chemicals to choose from:")
print("-Reference chemicals")
print(REF_CHEMICALS)
print("\n-PFAS chemicals")
print(PFAS_CHEMICALS)
print("Controls")
print(CONTROLS)
print("Durations to choose from:")
print("1 day, 10 day (default); 4 and 14 also in the data")
print("Default DESeq design: ~ chemical + duration")

  list(metadata_all = metadata_all, markdups_counts_all = markdups_counts_all)
}
# end of the build_master function()


# -----------------------------------------------------------------------------


# CHOOSE the chemicals, controls, and durations to choose from
# chemicals_chosen <- c("acetaminophen","aflatoxin B1","cyclosporin A","WY-14643","rosiglitazone","troglitazone","caffeine","sucrose","potassium chloride","PFOA","PFBS","PFOS","PFDS")
# controls_chosen <- c("DMSO_ref", "DMSO_pfas", "water")

# SECOND function
run_deseq_analysis <- function(chemicals_chosen, controls_chosen, duration_chosen,
                               metadata_all, markdups_counts_all, design_txt) {
  
  # combine chemicals_chosen and controls_chosen
  chemicals_and_controls_chosen <- c(chemicals_chosen, controls_chosen)
  
  # filter metadata_all to match chemicals_and_controls_chosen and duration_chosen
  metadata_all_filtered <- metadata_all %>%
    filter(chemical %in% chemicals_and_controls_chosen) %>%
    filter(duration %in% duration_chosen)
  rownames(metadata_all_filtered) <- metadata_all_filtered$sample_id
  
  # filter markdups_counts_all to match metadata_all_filtered
  # counts col samples == metadata rownames
  markdups_counts_all_filtered <- markdups_counts_all[, rownames(metadata_all_filtered)]
  
  # create dds
  dds <- DESeqDataSetFromMatrix(
    countData = markdups_counts_all_filtered, # sample gene counts unnormalized
    colData = metadata_all_filtered, # metadata info about samples
    # design = ~ chemical
    # design = ~ chemical + duration # make this a user choice in the UI, and this as default
    design = as.formula(design_txt)
  )
  
  # keep genes with at least smallestGroupSize reads in at least one sample
  smallestGroupSize <- min(count(metadata_all_filtered, chemical, dose)$n)
  keep <- rowSums(counts(dds) > 0) >= smallestGroupSize # lowest gene count to keep is 1
  # keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize # original
  dds <- dds[keep, ]
  
  # run DESeq: running the differential expression pipeline
  dds <- DESeq(dds)
  
  # create results table
  res <- results(dds)
  res <- res[order(res$padj), ]
  res$significant <- res$padj < 0.05 & abs(res$log2FoldChange) > 1
  
  vsd <- vst(dds, blind = TRUE) # Variance-stabilizing transformation
  # is this for the PCA plot? rplot can be used?
  
  # return the dds, results, vsd, and metadata_all_filtered
  list(dds = dds, results = res, vsd = vsd, metadata = metadata_all_filtered)
}


# RUN the FIRST function
masters <- build_masters()
# masters is a list containing metadata_all and markdups_counts_all
# master metadata and master counts for ref and pfas datasets

# retrieve the metadata and counts from the masters list after running above func()
metadata_all <- masters$metadata_all
# metadata_all is a dataframe containing the metadata for all samples
markdups_counts_all <- masters$markdups_counts_all
# markdups_counts_all is a dataframe containing the counts for all samples


# -------------------------------------------------------------
# CREATING the UI

ui <- fluidPage(
  titlePanel("PFAS / Reference Chemical DESeq2 Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("chemicals", "Chemical(s):",
                  choices = c(REF_CHEMICALS, PFAS_CHEMICALS),
                  multiple = TRUE, selected = "PFBS"),
      selectInput("controls", "Control:",
                  choices = CONTROLS, selected = "DMSO_pfas"),
      checkboxGroupInput("durations", "Duration(s):",
                          choices = DURATION_CHOICES, selected = c(1, 10)),
      checkboxGroupInput("design_terms",
                         paste("DESeq design (default: chemical + duration);",
                         "include duration if you selected all durations"),
                         choices = c("chemical", "duration", "dose"),
                         selected = c("chemical", "duration")
                         ),
      actionButton("run", "Run DESeq2"),
      hr(), # Horizontal dividing line
      downloadButton("download_data", "Download All Results & Plots (.zip)")
    ),
    mainPanel(
      plotOutput("pca_plot"),
      plotOutput("volcano_plot"),
      tableOutput("top_genes"),
      plotOutput("heatmap")
    )
  )
)

server <- function(input, output, session) {

  # terms <- input$design_terms # moved inside reactive
  # put chemical last so results() is the chemical contrast
  # terms <- c(setdiff(terms, "chemical"), "chemical") # moved this down to the run_deseq_analysis function
  # because it needs to be done after the user has selected the design terms

  # design_txt <- paste("~", paste(terms, collapse = " + "))
  # design_use <- as.formula(design_txt)

  # Run the DESeq2 analysis, but make sure user has selected all the necessary inputs in eventReactive()
  analysis <- eventReactive(input$run, {
    req(input$chemicals, input$controls, input$durations, input$design_terms)
    terms <- input$design_terms
    terms <- c(setdiff(terms, "chemical"), "chemical") # must put chemical last so results() is the chemical contrast
    # Note: In order to benefit from the default settings of the package, you should put the variable of interest at the end of the formula and make sure the control level is the first level
    run_deseq_analysis(
      chemicals_chosen   = input$chemicals,
      controls_chosen    = input$controls,
      duration_chosen    = as.numeric(input$durations),
      metadata_all        = metadata_all,
      markdups_counts_all = markdups_counts_all,
      design_txt = paste("~", paste(terms, collapse = " + "))
    )
  })

  # 3. PCA Plot
  output$pca_plot <- renderPlot({
    req(analysis()$vsd) # make sure analysis() is not NULL and has run
    plotPCA(analysis()$vsd, intgroup = "chemical")
  })

  #  4. Volcano Plot
  output$volcano_plot <- renderPlot({
    req(analysis()$results) # make sure analysis() is not NULL and has run
    res_df <- as.data.frame(analysis()$results)
    ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
      geom_point() +
      theme_minimal() +
      labs(x = "Log2 Fold Change", y = "-Log10 P-value")
  })

  # 1. Full DESeq2 Results -> 2. Top 20 Genes
  output$top_genes <- renderTable({
    req(analysis()$results) # make sure analysis() is not NULL and has run
    as.data.frame(head(analysis()$results, 20))
  }, rownames = TRUE)
  
  # 5. Distance Heatmap
  output$heatmap <- renderPlot({
    res <- analysis()
    req(res$vsd)
    # Heatmap of Sample Distances
    sampleDists <- dist(t(assay(res$vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    pheatmap(sampleDistMatrix, main = "Sample Distance Heatmap")
  })
  
  
  # Download handler for ZIP bundle
  output$download_data <- downloadHandler(
    filename = function() {
      paste0(paste("DESeq2_analysis_results", 
            Sys.Date(), 
            paste(input$chemicals, collapse = "_"), 
            paste(input$controls, collapse = "_"), 
            paste(input$durations, collapse = "_"), 
            paste(input$design_terms, collapse = "_"), 
            sep = "-"), ".zip")
    },
    content = function(file) {
      res <- analysis() # takes results from earlier deseq analysis
      req(res$vsd, res$results) # make sure we have analysis results first
      
      tmpdir <- tempdir() # path to temp folder by R
      owd <- setwd(tmpdir) # switches wd to that temporary for all generated things
                            # ALSO returns the prev wd = app dir
      on.exit(setwd(owd)) # when exit, set wd back to original working folder
                          # takes that app dir returned before
      
      files_to_zip <- c() # compile list of files to zip
      
      # 1. Full DESeq2 Results (CSV)
      res_df <- as.data.frame(res$results)
      csv_file <- "deseq2_results.csv"
      write.csv(res_df, csv_file, row.names = TRUE)
      files_to_zip <- c(files_to_zip, csv_file)
      
      # 2. Top 20 Genes (CSV)
      top_csv_file <- "top_20_genes.csv"
      write.csv(head(res_df, 20), top_csv_file, row.names = TRUE)
      files_to_zip <- c(files_to_zip, top_csv_file)
      
      # 3. PCA Plot (PNG)
      pca_file <- "pca_plot.png"
      png(pca_file, width = 800, height = 600, res = 120)
      print(plotPCA(res$vsd, intgroup = "chemical"))
      dev.off() # write the image file
      files_to_zip <- c(files_to_zip, pca_file)
      
      # 4. Volcano Plot (PNG)
      if (!"significant" %in% colnames(res_df)) {
        # Evaluates each gene's adjusted p-value (FDR)
        # return TRUE for genes stat threshold FDR<0.05 & magnitude of expression
        res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
      }
      p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
        geom_point(alpha = 0.6) +
        theme_minimal() +
        labs(x = "Log2 Fold Change", y = "-Log10 P-value")
      
      volcano_file <- "volcano_plot.png"
      ggsave(volcano_file, plot = p_volcano, width = 7, height = 5, dpi = 300)
      files_to_zip <- c(files_to_zip, volcano_file)
      
      # 5. Distance Heatmap (PNG)
      heatmap_file <- "sample_distance_heatmap.png"
      sampleDists <- dist(t(assay(res$vsd)))
      sampleDistMatrix <- as.matrix(sampleDists)
      
      png(heatmap_file, width = 800, height = 800, res = 120)
      pheatmap(sampleDistMatrix, main = "Sample Distance Heatmap")
      dev.off()
      files_to_zip <- c(files_to_zip, heatmap_file)
      
      # Zip everything together
      zip::zipr(zipfile = file, files = files_to_zip)
    }
  )
}



shinyApp(ui, server) 
# in console, run:
# setwd("/home/lawrence/Documents/School/thesis/PFAS_project_R-analysis")
# shiny::runApp("code/analysis_app.R")