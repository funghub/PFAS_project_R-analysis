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

library(DESeq2)
library(ggplot2)
library(pheatmap)

library(shiny)
library(zip)

# for PCA
library(plotly)

library(fgsea)
library(msigdbr)
library(data.table)

# DESeq2/S4Vectors and AnnotationDbi mask these; keep dplyr verbs
select <- dplyr::select
count <- dplyr::count
rename <- dplyr::rename

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
  dir.create("./results_from_analysis/data_input_deseq/metadata_ref", recursive = TRUE, showWarnings = FALSE)
  write.csv(SRR_table_ref,"./results_from_analysis/data_input_deseq/metadata_ref/ref_metadata_all.csv",row.names = FALSE)
  
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
  dir.create("./results_from_analysis/data_input_deseq/metadata_PFAS", recursive = TRUE, showWarnings = FALSE)
  write.csv(GSMtoSRR,"./results_from_analysis/data_input_deseq/metadata_PFAS/PFAS_metadata_all.csv",row.names = FALSE)
  
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
  
  
  dir.create("./results_from_analysis/data_input_deseq", recursive = TRUE, showWarnings = FALSE)
  write.csv(metadata_all, "./results_from_analysis/data_input_deseq/metadata_all.csv", row.names = FALSE)
  
  
  
  # CREATE the Master Counts File (PFAS dataset + Ref Dataset)
  
  ## Ref dataset counts (featureCounts raw, not markdups)
  raw_counts_ref <- read.table("./data/featureCounts/bowtie/ref/raw_counts.txt", header=TRUE, sep = "\t", row.names = 1)
  
  ## Replace everything after _ and the _ with nothing in all column names
  raw_counts_ref_cleancols <- raw_counts_ref %>%
    rename_with(~gsub("_.*", "", .x)) %>%
    select(-c("Chr", "Start", "End", "Strand", "Length"))
  
  ## PFAS dataset counts (featureCounts raw, not markdups)
  raw_counts_pfas <- read.table("./data/featureCounts/bowtie/pfas/raw_counts.txt", header=TRUE, sep = "\t", row.names = 1)
  
  # Replace everything after _ and the _ with nothing in all column names
  raw_counts_pfas_cleancols <- raw_counts_pfas %>%
    rename_with(~gsub("_.*", "", .x)) %>%
    select(-c("Chr", "Start", "End", "Strand", "Length"))
  
  # check that the row names are the same so that means genes are the same
  # then combine the columns if true
  if (identical(rownames(raw_counts_pfas_cleancols), rownames(raw_counts_ref_cleancols))) {
    raw_counts_all <- cbind(raw_counts_pfas_cleancols, raw_counts_ref_cleancols)
  } else {
    stop("Row names do NOT match between raw_counts_pfas_cleancols and raw_counts_ref_cleancols — check gene order/IDs before combining.")
  }
  
  # align metadata rows to counts columns, then rename SRR to chemical_duration_dose
  # (replicates share a combo, so make.unique appends _1, _2, ...)
  metadata_all <- metadata_all[match(colnames(raw_counts_all), metadata_all$SRR), ]
  stopifnot(identical(colnames(raw_counts_all), metadata_all$SRR))
  
  sample_id <- make.unique(
    paste(gsub(" ", "_", metadata_all$chemical), 
          paste0(metadata_all$duration, "d"), 
          paste0(metadata_all$dose,"uM"), #\u00b5M
          sep = "-"),
    sep = "-"
  )
  metadata_all$sample_id <- sample_id
  rownames(metadata_all) <- sample_id
  colnames(raw_counts_all) <- sample_id
  
  print("metadata_all validate sample_id START")
  print(head(metadata_all$sample_id))
  print(head(metadata_all$chemical))
  print(head(metadata_all$duration))
  print(head(metadata_all$dose))
  print("metadata_all validate sample_id END")
  
  
  write.csv(metadata_all, "./results_from_analysis/data_input_deseq/metadata_all.csv", row.names = TRUE)
  write.csv(tibble::rownames_to_column(raw_counts_all, "gene"),
            "./results_from_analysis/data_input_deseq/raw_counts_all.csv", row.names = FALSE)
  
  # validate on terminal
  cat('\nraw_counts_all #cols == metadata_all #rows\n')
  print(ncol(raw_counts_all) == nrow(metadata_all))
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
  
  list(metadata_all = metadata_all, raw_counts_all = raw_counts_all)
}
# end of the build_master function()


# -----------------------------------------------------------------------------


# CHOOSE the chemicals, controls, and durations to choose from
# chemicals_chosen <- c("acetaminophen","aflatoxin B1","cyclosporin A","WY-14643","rosiglitazone","troglitazone","caffeine","sucrose","potassium chloride","PFOA","PFBS","PFOS","PFDS")
# controls_chosen <- c("DMSO_ref", "DMSO_pfas", "water")

# SECOND function
run_deseq_analysis <- function(chemicals_chosen,controls_chosen,duration_chosen,
                               dosages_chosen,padj_value,log2fold_absolute,
                               sort_deseq_res, sort_decreasing,
                               metadata_all, raw_counts_all, 
                               design_txt) {
  
  # combine chemicals_chosen and controls_chosen
  # (in case a control was selected in chemicals_chosen)
  all_controls <- unique(c(controls_chosen, CONTROLS))
  chemicals_and_controls_chosen <- c(chemicals_chosen, controls_chosen)
  
  # filter metadata_all to match chemicals_and_controls_chosen and duration_chosen
  metadata_all_filtered <- metadata_all %>%
    filter(chemical %in% chemicals_and_controls_chosen) %>%
    filter(duration %in% duration_chosen) %>%
    # run each row in metadata through these statements, keep if one of these true
    filter(  
      chemical %in% all_controls |               # Keep controls regardless of dose, even if control is in chemicals_chosen
        # is.na(dose) |                                 # Safeguard: Keep any sample with a missing/NA dose value
        "all" %in% dosages_chosen |                    # test chemicals_chosen, Keep all if "All" is checked
        is.null(dosages_chosen) |                      # Keep all if input isn't rendered yet in UI
        dose %in% as.numeric(dosages_chosen)           # if dosage selected, Match selected doses for those chemicals_chosen
    )
  # if comparing control with control, MUST select ALL in the dosages chosen, or else FMSO_pfas would be dropped
  # b/c the dose is not NA, but 0, will fail second condition
  
  rownames(metadata_all_filtered) <- metadata_all_filtered$sample_id
  
  # filter raw_counts_all to match metadata_all_filtered
  # counts col samples == metadata rownames
  raw_counts_all_filtered <- raw_counts_all[, rownames(metadata_all_filtered)]
  
  # create dds
  dds <- DESeqDataSetFromMatrix(
    countData = raw_counts_all_filtered, # sample gene counts unnormalized
    colData = metadata_all_filtered, # metadata info about samples
    # design = ~ chemical
    # design = ~ chemical + duration # make this a user choice in the UI, and this as default
    design = as.formula(design_txt)
  )
  
  # specifying the reference level
  #   so that comparisons will NOT be based on the alpha order of the levels
  dds$chemical <- relevel(dds$chemical, ref = controls_chosen)
  cat("Controls Chosen", controls_chosen)
  
  # keep genes with at least smallestGroupSize reads in at least one sample
  smallestGroupSize <- min(count(metadata_all_filtered, chemical, dose)$n)
  keep <- rowSums(counts(dds) > 0) >= smallestGroupSize # lowest gene count to keep is 1
  # keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize # original
  dds <- dds[keep, ]
  
  # run DESeq: running the differential expression pipeline
  dds <- DESeq(dds)
  
  # create results table
  res <- results(dds)
  res <- res[order(res[[sort_deseq_res]], decreasing = sort_decreasing), ]
  # compute the significance column T/F
  # res$significant <- res$padj < 0.05 & abs(res$log2FoldChange) > 1 # can change threshold
  res$significant <- res$padj < padj_value & abs(res$log2FoldChange) > log2fold_absolute # can change threshold
  # THIS A possibility on shiny to edit these values interface slider bar
  
  # T/F those with small pval & logfold>1 == significant
  print("Summary Result Table after deseq")
  summary(res)
  
  vsd <- vst(dds, blind = TRUE) # Variance-stabilizing transformation
  # is this for the PCA plot? rplot can be used?
  
  # Prep for LFC: Log fold change shrinkage for visualization and ranking
  cat("resultsNames\n")
  print(resultsNames(dds))
  
  # results table with shrunken LFC
  # get the last item in vector for the comparison name
  coef_name <- resultsNames(dds)[length(resultsNames(dds))]
  resLFC <- lfcShrink(dds, coef=coef_name, type="apeglm")
  # apeglm has no `stat`; skip sort if the radio column has stat included
  if (sort_deseq_res %in% names(resLFC)) {
    resLFC <- resLFC[order(resLFC[[sort_deseq_res]], decreasing = sort_decreasing), ]
  }
  # recompute the significance column T/F as did with res
  # resLFC$significant <- resLFC$padj < 0.05 & abs(resLFC$log2FoldChange) > 1
  resLFC$significant <- resLFC$padj < padj_value & abs(resLFC$log2FoldChange) > log2fold_absolute
  # T/F those with small pval & logfold>1 == significant
  
  
  # Create Gene List for gsea based on results (not results lfc shrink)
  # gene_list <- res[[sort_deseq_res]] # extract the $stat column
  # Checks if the variable sort_deseq_res is set to either "stat"
  metric <- if (sort_deseq_res %in% c("stat", "log2FoldChange")) {
    sort_deseq_res
  }  else {
    "stat"
  }
  gene_list <- res[[metric]] # Pulls out the column corresponding to metric
  names(gene_list) <- rownames(res) # add row names into this df gene_list
  gene_list <- gene_list[!is.na(gene_list)] # remove na stat rows
  rnk_df <- data.frame("SYMBOL" = names(gene_list), "METRIC" = gene_list) %>%
    arrange(desc(METRIC))
  
  print("Number of Duplicate Gene Symbols:")
  print(sum(duplicated(rnk_df$SYMBOL)))# check if there are duplicates gene symbols
  
  
  # GSEA https://www.youtube.com/watch?v=Mi6u4r0lJvo
  # https://www.youtube.com/watch?v=Tm0LhciYxk8
  # with gene list, you can run GSEA
  # gse <- gseGO
  
  
  
  # return the dds, results, vsd, and metadata_all_filtered
  list(dds = dds, results = res, resultsLFC = resLFC, vsd = vsd,
       metadata = metadata_all_filtered, rank_genes = rnk_df)
}

# FGSEA function
run_fgsea_function <- function(pathways, ranks, minSize, maxSize) {
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    eps      = 0.0,
                    minSize  = minSize,
                    maxSize  = maxSize)
}

# RUN the FIRST function
masters <- build_masters()
# masters is a list containing metadata_all and raw_counts_all
# master metadata and master counts for ref and pfas datasets

# retrieve the metadata and counts from the masters list after running above func()
metadata_all <- masters$metadata_all
# metadata_all is a dataframe containing the metadata for all samples
raw_counts_all <- masters$raw_counts_all
# raw_counts_all is a dataframe containing the counts for all samples

# print(head(rownames(raw_counts_all))) #test

# -------------------------------------------------------------
# CREATING the UI

ui <- fluidPage(
  titlePanel("PFAS / Reference Chemical DESeq2 Explorer"),
  p("This interface only works for 1v1 comparisons. For example, adding more 
    than 1 chemical will render the dosage selection useless because the chemicals
    do not share the same range of dosages."),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("chemicals", "Chemical(s):",
                  choices = c(REF_CHEMICALS, PFAS_CHEMICALS,
                              CONTROLS # add controls here so that we can compare DMSO vs water
                  ),
                  multiple = TRUE, selected = "PFBS"),
      selectInput("controls", "Control:",
                  choices = CONTROLS, selected = "DMSO_pfas"),
      checkboxGroupInput("durations", "Duration(s):",
                         choices = DURATION_CHOICES, selected = c(1, 10)),
      p("If you choose one duration, you must unselect duration\" from DESeq 
      design b/c one duration would not explain the variation in the data."),
      
      uiOutput("dose_selector"), # sidebar for dosage
      
      checkboxGroupInput("design_terms",
                         paste("DESeq design (default: chemical + duration);",
                               "include duration if you selected >1 durations"),
                         choices = c("chemical", "duration", "dose"),
                         selected = c("chemical", "duration")
      ),
      
      p("Select the threshold value to label genes TRUE for stat threshold 
        , for example, FDR less than 0.05 & magnitude of expression greater than 1"),
      
      sliderInput("padj_val", "Adjusted P-values:",
                  min = 0, max = 1,
                  value = 0.05, step = 0.05),
      
      sliderInput("log2fold_abs", "Absolute Value Log2FoldChange:",
                  min = 0, max = 20,
                  value = 1, step = 0.50),
      
      # https://docs.gsea-msigdb.org/#GSEA/GSEA_User_Guide/#gseapreranked-page
      # First sort
      h3("First Sort"),
      h4("*requires to click Run DESeq Again if choose a different sort*"),
      radioButtons(inputId="sort_res", label="Sort the DESeq2 results by:", 
                   choices=c("log2FoldChange","lfcSE","stat","pvalue","padj"), selected = "stat"),
      radioButtons(inputId="sort_dir", label="Sort order:",
                   choices=c("Descending" = "decreasing",
                             "Ascending" = "ascending"),selected = "decreasing"),
      # Second sort
      h3("Second Sort"),
      h4("*DOES NOT require to click Run DESeq Again if choose a different sort*"),
      p("If the first sort is padj (adjusted p-value, not raw pvalue), genes with
        adjusted p-value above the Adjusted P-values slider are removed, then this sort runs."),
      radioButtons(inputId="sort_res2", label="Sort AGAIN by (not for GSEA rank genes):",
                   choices=c("None" = "none", "log2FoldChange","lfcSE","stat","pvalue","padj"),
                   selected = "none"),
      radioButtons(inputId="sort_dir2", label="Sort order for second sort (not for GSEA rank genes):",
                   choices=c("Descending" = "decreasing", "Ascending" = "ascending"),
                   selected = "decreasing"),
      
      actionButton("run", "Run DESeq2"),
      hr(), # Horizontal dividing line
      downloadButton("download_data", "Download All Results & Plots (.zip)"),
      p("File name will be: DESeq2_analysis_results-{date}-{chemicals}-{controls}-{durations}-{designterms}-{dosages}-{sortresultsby}-{sortorder}.zip")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Plots & Tables",
          h3("PCA Plot"),
          plotOutput("pca_plot"),
          hr(),
          
          h3("PCA Plot Interactive"),
          p("Sample Names are: {chemical}-{duration}-{dose}-{unique replicate}"),
          plotlyOutput("pca_plot_interactive"),
          hr(),
          
          h3("Volcano Plot"),
          plotOutput("volcano_plot"),
          hr(),
          
          h3("Volcano Plot Interactive"),
          plotlyOutput("volcano_plot_interactive"),
          hr(),
          
          h3("Top 20 Genes"),
          tableOutput("top_genes"),
          hr(),
          
          h3("Top 20 Genes (Shrunken LFC)"),
          p("If you second-sorted by 'stat', this table stays unsorted (apeglm drops stat)."),
          tableOutput("top_genes_lfcShrink"),
          hr(),
          
          h3("Gene Expression Heatmap (selected genes vs samples)"),
          textInput(
            inputId = "genes_interested", 
            label = "Enter Interested Genes from top pathways from GSEA or other (with comma with separator):", 
            value = "ALL"
          ),
          textInput(
            inputId = "samples_interested", 
            label = "Enter Samples Interested if you have specific ones other for further filtering after deseq: (enter as DMSO_pfas-1d-0uM with comma sep)", 
            value = "ALL"
          ),
          plotOutput("heatmap", height = "700px"),
          hr(),
          
          h3("Gene Expression Heatmap: Significant Genes Only (genes vs samples)"),          textOutput("heatmap_significant_text"),
          plotOutput("heatmap_significant",height = "700px"),
          hr(),
          
          h3("MA Plot"),
          plotOutput("ma_plot"),
          
          h3("Dispersion Plot"),
          plotOutput("dispersion_plot"),
          
          h3("Ranked Gene List"),
          tableOutput("ranked_genes")
        ),
        tabPanel(
          "FGSEA Analysis",
          h3("FGSEA Results"),
          p("Runs after DESeq2. Ranked genes use stat or log2FoldChange (same as the .rnk file)."),
          radioButtons(inputId="species", label="Species:",
                       choices=c("Human" = "Homo sapiens",
                                 "Mouse" = "Mus musculus"),
                       selected = "Homo sapiens"),
          radioButtons(inputId="collection", label="Gene Collections",
                       choices=c("Hallmark" = "H",
                                 "KEGG Medicus Pathways" = "CP:KEGG_MEDICUS",
                                 "Reactome Pathways" = "CP:REACTOME",
                                 "WikiPathways" = "CP:WIKIPATHWAYS",
                                 "GO Biological Process" = "GO:BP",
                                 "Cell Type Signature" = "C8"),
                       selected = "H"),
          sliderInput("minsize", "minSize gene overlap:",
                      min = 5, max = 100,
                      value = 15, step = 5),
          sliderInput("maxsize", "maxSize gene overlap:",
                      min = 50, max = 2000,
                      value = 500, step = 25),
          actionButton("run_fgsea", "Run FGSEA"),
          hr(),
          h4("Top pathways"),
          tableOutput("fgsea_table"),
          hr(),
          h4("NES per pathway"),
          plotOutput("fgsea_enrichment_plot"),
          hr(),
          h4("Top pathways with enrichment curves"),
          plotOutput("fgsea_enrichment_plot_labels"),
          hr(),
          selectInput("fgsea_pathway", "Enrichment plot for pathway:", choices = NULL),
          plotOutput("fgsea_plot")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # dose selector slider bar dynamic with the selected chemicals
  output$dose_selector <- renderUI({
    req(input$chemicals)
    
    # Filter metadata for selected chemicals (excluding controls)
    sub_meta <- metadata_all %>%
      filter(chemical %in% input$chemicals, !is.na(dose), dose > 0)
    
    avail_doses <- sort(unique(sub_meta$dose)) # get the dosages for chemicals chosen
    
    if (length(avail_doses) == 0) {
      return(NULL) # No dose variation present
    }
    
    checkboxGroupInput(
      inputId  = "dosages",
      label    = "Select Dosage(s) in uM (Optional):",
      choices  = c("All" = "all", avail_doses),
      selected = "all" # Defaults to selecting all doses
    )
  })
  
  
  # Run the DESeq2 analysis, but make sure user has selected all the necessary inputs
  analysis <- eventReactive(input$run, {
    req(input$chemicals, input$controls, input$durations, input$design_terms,
        input$padj_val, input$log2fold_abs, input$sort_res, input$sort_dir,
        !any(input$chemicals %in% input$controls) # Ensures no overlap between chemical and control selections
    )
    terms <- input$design_terms
    terms <- c(setdiff(terms, "chemical"), "chemical") # must put chemical last so results() is the chemical contrast
    # Note: In order to benefit from the default settings of the package, you should put the variable of interest at the end of the formula and make sure the control level is the first level
    stopifnot("chemical" %in% terms)
    
    run_deseq_analysis(
      chemicals_chosen   = input$chemicals,
      controls_chosen    = input$controls,
      duration_chosen    = as.numeric(input$durations),
      dosages_chosen     = input$dosages,
      padj_value         = input$padj_val,
      log2fold_absolute  = input$log2fold_abs,
      sort_deseq_res     = input$sort_res,
      sort_decreasing    = input$sort_dir == "decreasing",
      metadata_all        = metadata_all,
      raw_counts_all = raw_counts_all,
      design_txt = paste("~", paste(terms, collapse = " + "))
    )
  })
  
  # 3. PCA Plot
  output$pca_plot <- renderPlot({
    req(analysis()$vsd) # make sure analysis() is not NULL and has run
    plotPCA(analysis()$vsd, intgroup = "chemical")
  })
  
  # 3.2 PCA Plot Interactive
  output$pca_plot_interactive <- renderPlotly({
    req(analysis()$vsd) # make sure analysis() is not NULL and has run
    
    # just return the data.frame
    pca_df <- plotPCA(analysis()$vsd, intgroup = "chemical", returnData=TRUE)
    percent_var <- round(attributes(pca_df)$percentVar*100)
    
    # put the dataframe and use ggplot2 to plot the PCA
    df_in_ggplot <- ggplot(data = pca_df, 
                           aes(PC1, PC2, color = chemical, text = name)) +
      geom_point(size = 2) + # adding data points
      coord_fixed() + # fixing coordinates
      xlab(paste0("PC1: ", percent_var[1],"% variance")) +
      ylab(paste0("PC2: ", percent_var[2],"% variance"))
    
    # put the ggplot in a plotly wrapper to make interactive
    ggplotly(df_in_ggplot, tooltip = "text")
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
  
  # 4.1 Volcano Plot Interactive
  output$volcano_plot_interactive <- renderPlotly({
    req(analysis()$results)
    res_df <- as.data.frame(analysis()$results)
    res_df$gene <- rownames(res_df) # gene symbol for hover tooltip
    
    df_in_ggplot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue),
                                       color = significant, text = gene)) +
      geom_point() +
      theme_minimal() +
      labs(x = "Log2 Fold Change", y = "-Log10 P-value")
    
    ggplotly(df_in_ggplot, tooltip = "text")
  })
  
  # 1. Full DESeq2 Results -> 2. Top 20 Genes
  output$top_genes <- renderTable({
    req(analysis()$results) # make sure analysis() is not NULL and has run
    res_table <- analysis()$results
    
    # second sort (live, no need to re-run DESeq2)
    if (input$sort_res2 != "none") {
      # first sort was adjusted p-value: keep only padj below the slider, then sort
      if (input$sort_res == "padj") {
        res_table <- res_table[!is.na(res_table$padj) & res_table$padj < input$padj_val, ]
      }
      res_table <- res_table[
        order(res_table[[input$sort_res2]],
              decreasing = input$sort_dir2 == "decreasing"), ]
    }
    as.data.frame(head(res_table, 20))
  }, rownames = TRUE)
  
  
  # 1.1. Full DESeq2 Results -> 2. Top 20 Genes - lfcShrink
  output$top_genes_lfcShrink <- renderTable({
    req(analysis()$resultsLFC) # make sure analysis() is not NULL and has run
    LFCshrink_table <- analysis()$resultsLFC
    
    if (input$sort_res2 != "none") { 
      validate(need(input$sort_res2 %in% names(LFCshrink_table),
                    # the %in% is safety for if "stat" is still in the columns of table
                    message = "Shrunken LFC results have no 'stat' column. Choose another sort by."))
      if (input$sort_res == "padj") {
        LFCshrink_table <- LFCshrink_table[!is.na(LFCshrink_table$padj) & LFCshrink_table$padj < input$padj_val, ]
      }
      LFCshrink_table <- LFCshrink_table[
        order(LFCshrink_table[[input$sort_res2]],
              decreasing = input$sort_dir2 == "decreasing"), ]
    }
    as.data.frame(head(LFCshrink_table, 20))
  }, rownames = TRUE)
  
  
  # create the filtered heatmap reactive function (no filtering significant genes)
  heatmap_vsd_expression <- reactive({
    
    req(input$genes_interested, analysis()$vsd, input$samples_interested)
    
    #TEST
    v <- assay(analysis()$vsd)
    print(head(colnames(v), 50))
    print(setdiff(trimws(strsplit(input$genes_interested, ",")[[1]]), rownames(v)))
    print(setdiff(trimws(strsplit(input$samples_interested, ",")[[1]]), colnames(v)))
    
    
    # parse the genes out from the text (account whitespaces, and comma as sep)
    genes_no_white_spaces <-gsub(" ","",input$genes_interested)
    samples_no_white_spaces <-gsub(" ","",input$samples_interested)
    
    vsd_expression <- assay(analysis()$vsd)
    
    if (genes_no_white_spaces != "ALL") {
      # convert genes to vector
      genes_vector <- strsplit(genes_no_white_spaces, ",")[[1]]
      
      # filtering the genes I want
      vsd_expression <- vsd_expression[rownames(vsd_expression) %in% genes_vector, , drop = FALSE]
      validate(need(nrow(vsd_expression) > 0, "None of the entered genes were found."))
      
      #TEST
      print(dim(vsd_expression))
      print(setdiff(genes_vector, rownames(assay(analysis()$vsd))))
    }
    
    if (samples_no_white_spaces != "ALL") {
      # filtering the samples I want, based on the entered by user
      # convert genes to vector
      samples_vector <- strsplit(samples_no_white_spaces, ",")[[1]]
      sample_names_no_replicate <- sub("-[0-9]+$", "", colnames(vsd_expression))
      vsd_expression <- vsd_expression[, colnames(vsd_expression) %in% samples_vector |
                                         sample_names_no_replicate %in% samples_vector, drop = FALSE]
      validate(need(ncol(vsd_expression) > 1, "Need at least 2 matching samples."))
    }
    print(dim(vsd_expression))
    vsd_expression
  })
  
  # 5. Distance Heatmap (no filtering significant)
  output$heatmap <- renderPlot({
    vsd_expression <- heatmap_vsd_expression()

    vsd_expression <- vsd_expression[apply(vsd_expression, 1, var) > 0, , drop = FALSE]
    validate(need(nrow(vsd_expression) >= 1, "Selected genes have no variance across these samples."))
    
    pheatmap(vsd_expression,
             scale = if (nrow(vsd_expression) > 1) "row" else "none", # if > 1 gene = centered/scaled row direction, else none
             cluster_rows = nrow(vsd_expression) > 1, # >1 gene = cluster row
             cluster_cols = ncol(vsd_expression) > 2, # >2 samples = cluster col
             show_rownames = TRUE, show_colnames = TRUE,
             main = "Selected genes (vst)")
  })
  
  # 5.1 Distance Heatmap (Filtered vsd to remove )
  # 2. Filter for significant genes (e.g., padj < 0.05 and absolute log2 fold change > 1)
  # sig_genes <- rownames(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))
  # 3. Subset your VST matrix to ONLY include those significant genes
  # mat <- assay(res$vsd)[sig_genes, ]
  # reactive: vsd expression of significant genes only (all samples)
  heatmap_vsd_expression_significant <- reactive({
    req(analysis()$results, analysis()$vsd, input$genes_interested, input$samples_interested)  
    res <- analysis()$results
    vsd_expression <- assay(analysis()$vsd)
    
    # Filter for interested samples
    samples_no_white_spaces <- gsub(" ", "", input$samples_interested)
    if (samples_no_white_spaces != "ALL") {
      samples_vector <- strsplit(samples_no_white_spaces, ",")[[1]]
      
      sample_names_no_replicate <- sub("-[0-9]+$", "", colnames(vsd_expression)) #combine all replictes
      vsd_expression <- vsd_expression[, colnames(vsd_expression) %in% samples_vector | # match exact specific replicate
                                         sample_names_no_replicate %in% samples_vector, drop = FALSE] # Match group basename and all its replicates
      validate(need(ncol(vsd_expression) > 1,
                    paste("Only", ncol(vsd_expression), "matching sample(s).")))
    }
    
    # Filtering interested genes
    sig_genes <- rownames(res)[which(res$significant)] # which() drops NA padj rows
    validate(need(length(sig_genes) > 1, "Fewer than 2 significant genes at the chosen thresholds."))
    
    # Step 1: Keep only significant genes
    vsd_expression <- vsd_expression[rownames(vsd_expression) %in% sig_genes, , drop = FALSE]
    
    # Step 2: FURTHER filter by user-typed genes (if not "ALL")
    genes_no_white_spaces <- gsub(" ", "", input$genes_interested)
    if (genes_no_white_spaces != "ALL") {
      genes_vector <- strsplit(genes_no_white_spaces, ",")[[1]]
      vsd_expression <- vsd_expression[rownames(vsd_expression) %in% genes_vector, , drop = FALSE]
      validate(need(nrow(vsd_expression) > 0, "None of the entered genes were found among the significant genes."))
    }
    
    
    print(dim(vsd_expression))
    vsd_expression
  })
  
  
  output$heatmap_significant_text <- renderText({
    vsd_expression <- heatmap_vsd_expression_significant()
    paste0("Filtered to ", nrow(vsd_expression), " significant genes (padj < ",
           input$padj_val, " and |log2FoldChange| > ", input$log2fold_abs,
           "), using the thresholds set when Run DESeq2 was clicked.")
  })
  
  output$heatmap_significant <- renderPlot({
    vsd_expression <- heatmap_vsd_expression_significant()
    
    vsd_expression <- vsd_expression[apply(vsd_expression, 1, var) > 0, , drop = FALSE]
    validate(need(nrow(vsd_expression) >= 2, "Fewer than 2 significant genes with variance in these samples."))
    
    pheatmap(vsd_expression, scale = "row", # for each gene calculate z-score (mean expression gene used in z-score)
             cluster_cols = ncol(vsd_expression) > 2,
             show_rownames = nrow(vsd_expression) <= 60,
             show_colnames = ncol(vsd_expression) <= 40,
             main = paste0("Significant genes (vst, row z-score): ", nrow(vsd_expression), " genes"))
    })
  
  
  # 6. MA-Plot
  output$ma_plot <- renderPlot({
    req(analysis())
    plotMA(analysis()$resultsLFC, ylim = c(-2, 2))
  })
  
  # 6.1 Dispersion Plot
  output$dispersion_plot <- renderPlot({
    req(analysis())
    plotDispEsts(analysis()$dds)
    # plotDispEsts(analysis()$dds,fitType = "local")
  })
  
  #7. Gene List for GSEA
  output$ranked_genes <- renderTable({
    req(analysis()$rank_genes) # make sure analysis() is not NULL and has run
    as.data.frame(head(analysis()$rank_genes, 20))
  }, rownames = TRUE)
  
  # FGSEA new tab
  # retrieve gene set (only after Run FGSEA is clicked)
  genesets <- reactive({
    if (input$collection %in% c("H", "C8")) {
      msigdbr(species = input$species, collection = input$collection)
    } else {
      msigdbr(species = input$species, subcollection = input$collection)
    }
  }) %>% bindEvent(input$run_fgsea)
  
  # running fgsea
  fgsea_run <- reactive({
    req(analysis()$rank_genes)
    
    gs <- genesets() # retrieve geneset from MSigDB
    
    # divide the geneset into gene symbols per geneset name
    pathways <- split(gs$gene_symbol, gs$gs_name)
    # https://stephenturner.github.io/deseq-to-fgsea/#load_some_results_from_deseq2
    ranks <- tibble::deframe(analysis()$rank_genes) # df -> named vector (SYMBOL -> METRIC)
    
    table <- run_fgsea_function(
      pathways = pathways,
      ranks    = ranks,
      minSize  = input$minsize,
      maxSize  = input$maxsize
    )
    
    list(table = table, pathways = pathways, stats = ranks)
  }) %>% bindEvent(input$run_fgsea)
  
  # populate pathway dropdown after fgsea runs
  observe({
    updateSelectInput(session, "fgsea_pathway",
                      choices = names(fgsea_run()$pathways))
  }) %>% bindEvent(fgsea_run())
  
  output$fgsea_table <- renderTable({
    req(fgsea_run()$table)
    tbl <- as.data.frame(fgsea_run()$table) # convert table to data frame
    tbl <- tbl[order(tbl$padj), ] # most significant first
    tbl$leadingEdge <- vapply(tbl$leadingEdge, function(x) paste(x, collapse = ", "), character(1)) # show the genes nicely
    head(tbl[, c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")], 30)
  })
  
  # plot normalized enrichment scores of pathways
  output$fgsea_enrichment_plot <- renderPlot({
    req(fgsea_run()$table)
    tbl <- as.data.frame(fgsea_run()$table)
    tbl <- tbl[order(tbl$padj), ]
    
    ggplot(tbl, aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill = padj < 0.05)) +
      coord_flip() +
      labs(x = "Pathway", y = "Normalized Enrichment Score",
           title = "Pathways NES from GSEA") + 
      theme_minimal()
  }, height = function() {
    max(400, nrow(fgsea_run()$table) * 20) # 20px per pathway, minimum 400px
  })
  
  # top up/down pathways with enrichment curves (data.table syntax, table kept as data.table)
  output$fgsea_enrichment_plot_labels <- renderPlot({
    req(fgsea_run()$table)
    tbl <- fgsea_run()$table
    tbl <- tbl[order(padj)]
    
    topPathwaysUp <- tbl[ES > 0][head(order(pval), n = 10), pathway]
    topPathwaysDown <- tbl[ES < 0][head(order(pval), n = 10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    fgsea::plotGseaTable(fgsea_run()$pathways[topPathways],
                         fgsea_run()$stats,
                         tbl,
                         gseaParam = 0.5)
  })
  
  output$fgsea_plot <- renderPlot({
    req(fgsea_run(), input$fgsea_pathway)
    pw <- input$fgsea_pathway
    req(pw %in% names(fgsea_run()$pathways))
    fgsea::plotEnrichment(fgsea_run()$pathways[[pw]], fgsea_run()$stats) +
      ggplot2::labs(title = pw)
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
                   paste(input$dosages, collapse = "_"),
                   paste(input$sort_res, collapse = "_"),
                   paste(input$sort_dir, collapse = "_"),
                   sep = "-"), ".zip")
    },
    content = function(file) {
      res <- analysis() # takes results from earlier deseq analysis
      if (!isTruthy(input$run_fgsea)) {
        showNotification("Run FGSEA first.", type = "error")
      }
      req(res$results, input$genes_interested, analysis()$vsd, input$samples_interested, fgsea_run()$table) # make sure we have analysis results first
      
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
      
      
      # VST results
      vsd_df <- as.data.frame(assay(res$vsd))
      csv_file <- "vsd_results.csv"
      write.csv(vsd_df, csv_file, row.names = TRUE)
      files_to_zip <- c(files_to_zip, csv_file)
      
      
      # 2. Top 20 Genes (CSV)
      top_csv_file <- "top_20_genes.csv"
      write.csv(head(res_df, 20), top_csv_file, row.names = TRUE)
      files_to_zip <- c(files_to_zip, top_csv_file)
      
      # lfcShrink
      # 1.1. Full DESeq2 Results (CSV) - lfcShrink
      res_df <- as.data.frame(res$resultsLFC)
      csv_file <- "deseq2_results-LFC.csv"
      write.csv(res_df, csv_file, row.names = TRUE)
      files_to_zip <- c(files_to_zip, csv_file)
      
      # 2.1. Top 20 Genes (CSV) - lfcShrink
      top_csv_file <- "top_20_genes-LFC.csv"
      write.csv(head(res_df, 20), top_csv_file, row.names = TRUE)
      files_to_zip <- c(files_to_zip, top_csv_file)
      
      # 3. PCA Plot (PNG)
      pca_file <- "pca_plot.png"
      png(pca_file, width = 800, height = 600, res = 120)
      print(plotPCA(res$vsd, intgroup = "chemical"))
      dev.off() # write the image file
      files_to_zip <- c(files_to_zip, pca_file)
      
      # 4. Volcano Plot (PNG)
      # if (!"significant" %in% colnames(res_df)) {
      #   # Evaluates each gene's adjusted p-value (FDR)
      #   # return TRUE for genes stat threshold FDR<0.05 & magnitude of expression
      #   res_df$significant <- res_df$padj < input$padj_val & abs(res_df$log2FoldChange) > input$log2fold_abs
      # }
      
      res_df_volcano <- as.data.frame(res$results)
      p_volcano <- ggplot(res_df_volcano, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
        geom_point(alpha = 0.6) +
        theme_minimal() +
        labs(x = "Log2 Fold Change", y = "-Log10 P-value")
      
      volcano_file <- "volcano_plot.png"
      ggsave(volcano_file, plot = p_volcano, width = 7, height = 5, dpi = 300)
      files_to_zip <- c(files_to_zip, volcano_file)
      
      # 5. Distance Heatmap ALL genes and samples (PNG)
      heatmap_file <- "sample_all_gene_all_distance_heatmap.png"
      sampleDists <- dist(t(assay(res$vsd)), method="euclidean")
      sampleDistMatrix <- as.matrix(sampleDists)
      
      png(heatmap_file, width = 800, height = 800, res = 120)
      pheatmap(sampleDistMatrix, main = "Sample Distance All and Genes All Heatmap")
      dev.off()
      files_to_zip <- c(files_to_zip, heatmap_file)
      
      
      # 5.1 Gene Expression Heatmap specific genes and samples (PNG)
      heatmap_file_filtered <- "gene_sample_filtered_heatmap.png"
      
      vsd_expression <- heatmap_vsd_expression()
      vsd_expression_filtered <- vsd_expression[apply(vsd_expression, 1, var) > 0, , drop = FALSE]
      
      png(heatmap_file_filtered, width = 800, height = 800, res = 120)
      pheatmap(vsd_expression_filtered,
               scale = if (nrow(vsd_expression_filtered) > 1) "row" else "none",
               cluster_rows = nrow(vsd_expression_filtered) > 1,
               cluster_cols = ncol(vsd_expression_filtered) > 2,
               main = "Selected genes (vst)")
      dev.off()
      files_to_zip <- c(files_to_zip, heatmap_file_filtered)
      
      # 5.1 Gene Expression Heatmap specific genes and samples .csv (EXPORT the vst of filtered samples and genes)
      csv_file <- "gene_sample_filtered_heatmap.csv"
      write.csv(vsd_expression, csv_file, row.names = TRUE)
      files_to_zip <- c(files_to_zip, csv_file)
      
      
      # 5.3 Gene Expression Heatmap significant genes only (PNG + CSV)
      vsd_expression_sig <- tryCatch(heatmap_vsd_expression_significant(), error = function(e) NULL)
      
      if (!is.null(vsd_expression_sig)) {
        heatmap_file_significant <- "gene_sample_significant_heatmap.png"
        vsd_expression_sig_filtered <- vsd_expression_sig[apply(vsd_expression_sig, 1, var) > 0, , drop = FALSE]
        
        png(heatmap_file_significant, width = 1000, height = 1000, res = 120)
        pheatmap(vsd_expression_sig_filtered, scale = "row",
                 cluster_cols = ncol(vsd_expression_sig_filtered) > 2,
                 show_rownames = nrow(vsd_expression_sig_filtered) <= 60,
                 show_colnames = ncol(vsd_expression_sig_filtered) <= 40,
                 main = paste0("Significant genes (vst, row z-score): ", nrow(vsd_expression_sig_filtered), " genes"))
        dev.off()
        files_to_zip <- c(files_to_zip, heatmap_file_significant)
        
        csv_file <- "gene_sample_significant_heatmap.csv"
        write.csv(vsd_expression_sig, csv_file, row.names = TRUE)
        files_to_zip <- c(files_to_zip, csv_file)
      }
      
      
      # 6. MA Plot (PNG)
      ma_file <- "ma_plot.png"
      png(ma_file, width = 800, height = 600, res = 120)
      plotMA(res$resultsLFC, ylim = c(-2, 2))
      dev.off() # write the image file
      files_to_zip <- c(files_to_zip, ma_file)
      
      
      #7. Gene List for GSEA
      ranked_genes_file <- "deseq_ranked.rnk"
      write.table(res$rank_genes, ranked_genes_file,
                  sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      files_to_zip <- c(files_to_zip, ranked_genes_file)
      
      # 8. FGSEA results (only if FGSEA has been run)
      if (!is.null(fgsea_run()) && nrow(fgsea_run()$table) > 0) {
        fgsea_df <- as.data.frame(fgsea_run()$table)
        fgsea_df <- fgsea_df[order(fgsea_df$padj), ]
        fgsea_df$leadingEdge <- vapply(
          fgsea_df$leadingEdge, function(x) paste(x, collapse = ","), character(1)
        )
        fgsea_file <- paste0("fgsea_",
                             gsub(" ", "-", input$species), "_",
                             gsub(":", "-", input$collection), ".csv")
        write.csv(fgsea_df, fgsea_file, row.names = FALSE)
        files_to_zip <- c(files_to_zip, fgsea_file)
      }
      
      # Zip everything together
      zip::zipr(zipfile = file, files = files_to_zip)
    }
  )
}



shinyApp(ui, server) 
# in console, run:
# setwd("/home/lawrence/Documents/School/thesis/PFAS_project_R-analysis")
# shiny::runApp("code/analysis_app.R")