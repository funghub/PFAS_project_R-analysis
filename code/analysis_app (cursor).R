
knitr::opts_knit$set(root.dir = here::here())
knitr::opts_chunk$set(echo = TRUE)
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

setwd(here::here())

REF_CHEMICALS <- c("acetaminophen","aflatoxin B1","cyclosporin A","WY-14643",
                   "rosiglitazone","troglitazone","caffeine","sucrose","potassium chloride")
PFAS_CHEMICALS <- c("PFOA","PFBS","PFOS","PFDS")
CONTROLS <- c("DMSO_ref", "DMSO_pfas", "water")
DURATION_CHOICES <- c(1, 10, 4, 14)
DEFAULT_DESIGN <- "~ chemical + duration"

# CREATE master metadata + counts (everything above the old user-choice / DESEQ section)
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

# one comparison: filter masters, DESeq, graphs
# design default ~ chemical + duration:
#   is there a difference between chemical and control, with duration adjusted out.
#   Dose stays out — DMSO does not span a dose range.
run_one_comparison <- function(masters, chemicals_chosen, controls_chosen, duration_chosen, design_txt) {
  chemicals_and_controls_chosen <- c(chemicals_chosen, controls_chosen)

  metadata_all_filtered <- masters$metadata_all %>%
    filter(chemical %in% chemicals_and_controls_chosen) %>%
    filter(duration %in% duration_chosen)
  rownames(metadata_all_filtered) <- metadata_all_filtered$sample_id

  markdups_counts_all_filtered <- masters$markdups_counts_all[, rownames(metadata_all_filtered), drop = FALSE]

  cat('\nmarkdups_counts_all_filtered #cols == metadata_all_filtered #rows\n')
  print(ncol(markdups_counts_all_filtered) == nrow(metadata_all_filtered))
  print(unique(as.character(metadata_all_filtered$chemical)))
  print(unique(metadata_all_filtered$duration))

  stopifnot(identical(colnames(markdups_counts_all_filtered), rownames(metadata_all_filtered)))
  stopifnot(ncol(markdups_counts_all_filtered) == nrow(metadata_all_filtered))

  metadata_all_filtered$chemical <- factor(metadata_all_filtered$chemical)
  metadata_all_filtered$duration <- factor(metadata_all_filtered$duration)
  if (!controls_chosen %in% levels(metadata_all_filtered$chemical)) {
    stop("Control '", controls_chosen, "' is not in the filtered samples. Check chemical/control/duration.")
  }
  metadata_all_filtered$chemical <- relevel(metadata_all_filtered$chemical, ref = controls_chosen)

  design_use <- as.formula(design_txt)
  if ("duration" %in% all.vars(design_use) && nlevels(metadata_all_filtered$duration) < 2) {
    cat("Only one duration in this subset; using ~ chemical instead of ", design_txt, "\n", sep = "")
    design_use <- ~ chemical
  }

  dds <- DESeqDataSetFromMatrix(
    countData = markdups_counts_all_filtered,
    colData = metadata_all_filtered,
    design = design_use
  )

  smallestGroupSize <- min(count(metadata_all_filtered, chemical, dose)$n)
  keep <- rowSums(counts(dds) > 0) >= smallestGroupSize
  dds <- dds[keep, ]

  dds <- DESeq(dds)

  vsd <- vst(dds, blind = TRUE)
  pca_plot <- plotPCA(vsd, intgroup = "chemical")
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)

  res_by_chem <- lapply(chemicals_chosen, function(ch) {
    res <- results(dds, contrast = c("chemical", ch, controls_chosen))
    res <- res[order(res$padj), ]
    res$significant <- res$padj < 0.05 & abs(res$log2FoldChange) > 1
    res
  })
  names(res_by_chem) <- chemicals_chosen

  list(
    label = paste0(paste(chemicals_chosen, collapse = "+"), " vs ", controls_chosen,
                   " | days ", paste(duration_chosen, collapse = ","),
                   " | ", design_txt),
    chemicals = chemicals_chosen,
    control = controls_chosen,
    dds = dds,
    res_by_chem = res_by_chem,
    pca_plot = pca_plot,
    sampleDistMatrix = sampleDistMatrix
  )
}

run_logged <- function(fun) {
  msgs <- character()
  res <- NULL
  out <- capture.output({
    withCallingHandlers(
      tryCatch(
        res <<- fun(),
        error = function(e) {
          msgs <<- c(msgs, paste("ERROR:", conditionMessage(e)))
          NULL
        }
      ),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
  }, split = TRUE)
  list(result = res, log = paste(c(out, msgs), collapse = "\n"))
}

ui <- fluidPage(
  titlePanel("PFAS DESeq2"),
  p("Defaults: duration 1 and 10 days. Design ", code("~ chemical + duration"),
    " (chemical contrast, duration adjusted out). Dose is not in the default design."),
  tabsetPanel(
    id = "tabs",
    tabPanel(
      "Comparisons",
      br(),
      fluidRow(
        column(
          4,
          checkboxGroupInput("chems_ref", "Reference chemicals", REF_CHEMICALS),
          checkboxGroupInput("chems_pfas", "PFAS chemicals", PFAS_CHEMICALS, selected = "PFOA"),
          selectInput("control", "Control", CONTROLS, selected = "DMSO_pfas"),
          checkboxGroupInput("duration", "Duration (days)", DURATION_CHOICES, selected = c(1, 10)),
          textInput("design", "DESeq design", value = DEFAULT_DESIGN),
          helpText("Default is ~ chemical + duration. Examples: ~ chemical  or  ~ chemical + duration"),
          actionButton("add_comp", "Add comparison"),
          actionButton("remove_last", "Remove last"),
          actionButton("clear_comps", "Clear all")
        ),
        column(
          8,
          h4("Queued comparisons"),
          tableOutput("comp_table"),
          br(),
          actionButton("run", "Run setup + DESeq", class = "btn-primary"),
          helpText("Runs master metadata/counts first (Log tab), then DESeq for each comparison.")
        )
      )
    ),
    tabPanel(
      "Log",
      br(),
      p("Terminal output from setup (everything above DESEQ) and from each comparison:"),
      verbatimTextOutput("log", placeholder = TRUE)
    ),
    tabPanel(
      "Results",
      br(),
      selectInput("which_comp", "Comparison", choices = NULL),
      selectInput("which_chem", "Chemical contrast (vs control)", choices = NULL),
      h4("Results (first 50 rows)"),
      tableOutput("res_head"),
      h4("PCA"),
      plotOutput("pca", height = "400px"),
      h4("Sample distance heatmap"),
      plotOutput("heatmap", height = "500px"),
      h4("Volcano"),
      plotOutput("volcano", height = "400px")
    )
  )
)

server <- function(input, output, session) {
  comps <- reactiveVal(data.frame(
    chemicals = character(),
    control = character(),
    duration = character(),
    design = character(),
    stringsAsFactors = FALSE
  ))
  log_rv <- reactiveVal("Add comparisons, then click Run. Setup log will appear here.\n")
  results_rv <- reactiveVal(list())

  append_log <- function(txt) {
    log_rv(paste0(log_rv(), txt, "\n"))
  }

  output$comp_table <- renderTable(comps())
  output$log <- renderText(log_rv())

  observeEvent(input$add_comp, {
    chems <- c(input$chems_ref, input$chems_pfas)
    if (length(chems) == 0 || length(input$control) == 0 || length(input$duration) == 0) {
      showNotification("Pick at least one chemical, a control, and a duration.", type = "error")
      return()
    }
    if (!nzchar(trimws(input$design))) {
      showNotification("Design cannot be empty (default is ~ chemical + duration).", type = "error")
      return()
    }
    comps(rbind(comps(), data.frame(
      chemicals = paste(chems, collapse = ", "),
      control = input$control,
      duration = paste(input$duration, collapse = ", "),
      design = input$design,
      stringsAsFactors = FALSE
    )))
  })

  observeEvent(input$remove_last, {
    df <- comps()
    if (nrow(df) > 0) comps(df[-nrow(df), , drop = FALSE])
  })

  observeEvent(input$clear_comps, {
    comps(comps()[0, ])
  })

  observeEvent(input$run, {
    if (nrow(comps()) == 0) {
      showNotification("Add at least one comparison first.", type = "error")
      return()
    }
    updateTabsetPanel(session, "tabs", selected = "Log")
    log_rv("=== Running setup (master metadata + counts) ===\n")
    results_rv(list())

    withProgress(message = "Setup + DESeq", value = 0, {
      incProgress(0.1, detail = "Building master files (GEO + counts)...")
      setup <- run_logged(build_masters)
      append_log(setup$log)
      if (is.null(setup$result)) {
        append_log("Setup failed.")
        return()
      }
      masters <- setup$result

      out <- list()
      n <- nrow(comps())
      for (i in seq_len(n)) {
        row <- comps()[i, ]
        chems <- trimws(strsplit(row$chemicals, ",", fixed = TRUE)[[1]])
        durs <- as.numeric(trimws(strsplit(row$duration, ",", fixed = TRUE)[[1]]))
        incProgress(0.1 + 0.8 * i / n, detail = paste("DESeq", i, "of", n))
        append_log(paste0("\n=== Comparison ", i, "/", n, ": ", row$chemicals, " vs ", row$control, " ===\n"))
        one <- run_logged(function() {
          run_one_comparison(masters, chems, row$control, durs, row$design)
        })
        append_log(one$log)
        if (!is.null(one$result)) out[[one$result$label]] <- one$result
      }
      results_rv(out)
      if (length(out) > 0) {
        updateSelectInput(session, "which_comp", choices = names(out), selected = names(out)[[1]])
        append_log("\nDone. Open the Results tab.")
      }
    })
  })

  observe({
    req(input$which_comp)
    res <- results_rv()
    if (!input$which_comp %in% names(res)) return()
    updateSelectInput(session, "which_chem",
                      choices = res[[input$which_comp]]$chemicals,
                      selected = res[[input$which_comp]]$chemicals[[1]])
  })

  current <- reactive({
    req(input$which_comp)
    results_rv()[[input$which_comp]]
  })

  current_res <- reactive({
    cmp <- current()
    req(cmp, input$which_chem)
    cmp$res_by_chem[[input$which_chem]]
  })

  output$res_head <- renderTable({
    res <- current_res()
    df <- as.data.frame(res)
    df$gene <- rownames(df)
    head(df[, c("gene", "baseMean", "log2FoldChange", "pvalue", "padj", "significant")], 50)
  })

  output$pca <- renderPlot({
    print(current()$pca_plot)
  })

  output$heatmap <- renderPlot({
    pheatmap(current()$sampleDistMatrix, main = "Sample Distance Heatmap")
  })

  output$volcano <- renderPlot({
    res <- as.data.frame(current_res())
    ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
      geom_point() +
      theme_minimal() +
      labs(title = paste(input$which_chem, "vs", current()$control),
           x = "Log2 Fold Change", y = "-Log10 P-value")
  })
}

shinyApp(ui, server)
