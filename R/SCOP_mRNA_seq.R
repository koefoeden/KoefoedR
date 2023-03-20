# Tidyverse imports -----------------------------------------------------------------
#' @import ggplot2
#' @import dplyr
#' @import edgeR
#' @import purrr
#' @import readr
#' @import tidyr
#' @import tibble
#' @import stringr
#' @import forcats
NULL


# Data processing -------------------------------------------------------------------
#' Generate external gene names for ENSEMBL-IDS.
#'
#' @param dataset Choose human or mouse biomart ensembl dataset 
#' @param gene_counts_raw Raw salmon output with ensemble IDS as rownames.
#'
#' @return The gene info list
#' @export
generate_gene_info_list <- function(dataset="hsapiens_gene_ensembl", gene_counts_raw) {
  
  gene_info_list <- biomaRt::getBM(filters = "ensembl_gene_id",
                                   attributes = c("ensembl_gene_id","external_gene_name", "description"),
                                   values = rownames(gene_counts_raw),
                                   mart = biomaRt::useDataset(dataset, biomaRt::useMart("ensembl"))) %>%
    dplyr::rename("Name"=external_gene_name,
                  "ENSEMBL_ID"=ensembl_gene_id) %>%
    select(ENSEMBL_ID, Name, description)
  
  return(gene_info_list)
}


#' Get a named vector that can be used to translate ENSEMBL-ids into gene synmbols
#'
#' @param dataset The biomart dataset 
#' @param gene_counts_raw the gne counts raw object to extract ENSEMBL_IDS from
#'
#' @return A named vector that can be used for translation
#' @export
get_translate_vec <- function(dataset="hsapiens_gene_ensembl", gene_counts_raw) {
  gene_info_list <- generate_gene_info_list(dataset, gene_counts_raw)
  
  symbol_translate_vec <- gene_info_list %>% 
    select(ENSEMBL_ID, Name) %>% deframe()
  
  return(symbol_translate_vec)
}

#' Get MD data for all samples in long tibble
#'
#' @param object The DGElist object with transcript/gene counts and sample information
#' encoded in the y$sample object
#' @return A long tibble with three columns: sample label, mean and diff values
#' @export


get_MD_data_for_all_samples_updated <- function(object) {
  sample_indices_w_names <- colnames(object) %>%
    purrr::set_names()
  
  MA_plot_data_long <- purrr::map(sample_indices_w_names,
                                  ~cbmr::get_MD_data_for_sample(object = object,
                                                                column =.)) %>%
    dplyr::bind_rows(.id="Sample") %>%
    dplyr::mutate(Sample=factor(Sample,
                                levels=gtools::mixedsort(colnames(object))))
}


#' Get number of significant genes/GO-terms per contrast
#' @param results_DF_list An list of either DE genes or GO-terms results,
#' with an element for each contrast.
#' @return A facetted ggplot
#' @export
n_sig_genes_pr_contrast <- function(results_DF_list) {
  
  purrr::map(results_DF_list,
             ~.x %>%
               filter(if_any(any_of(c("FDR", "adj.P.Val")),
                             ~.x < 0.05)) %>%
               nrow()) %>%
    unlist()
}



#' Get number of significant DEGs or GO-terms from dataframe
#' @param df Data frame as produced by get_GO/DE_datatable
#' @param genes_type character vector, NCBI or ENSEMBL
#' @param type DE or GO
#'
#' @return character vector: gene names or ensembl IDs with FDR<0.5
#' @export
get_sig_entities_from_df <- function(df, genes_type="NCBI", type="DE") {
  
  sig_entities <- df %>%
    filter(FDR<0.05) %>%
    select(any_of(c(case_when(genes_type=="ENSEMBL" & type == "DE"  ~ c("ENSEMBL_ID"),
                              genes_type=="NCBI"& type == "DE" ~ c("Name"),
                              TRUE ~ "TERM")))) %>%
    pull()
  return(sig_entities)
}

#' Get number of significant DEGs or GO-terms from tsv file as produced by get_DE_datatable
#' @param path character vector: Path to the .tsv file
#' @param type character vector: DE or GO
#' @param genes_type character vector: DE or GO
#' @return character vector: gene names or ensembl IDs with FDR<0.5
#' @export
get_sig_entities_from_path <- function(path, type, genes_type) {
  
  sig_entities <- path %>%
    read_tsv %>%
    get_sig_entities_from_df(type=type, genes_type = genes_type)
  return(sig_entities)
  
}

#' Gets a nested list of significant entities (genes or GO-terms) between contrasts across tissues.
#' Assumes there is a data dir with each tissue described in its name, with GO and DE-tables in them
#' with corresponding names.
#' @param type Either "DE" or "GO" depending on the type of entity
#' @param tissue_str Character vector: Comma-seperated list of tissues as present in the project root, i.e. Hyp,NAc
#' @return A nested list, where the upper levels corresponds to supplied  contrasts,
#' and the lower levels each contain a vector of significant gene names or GO-terms for the tissue
#' @export
get_sig_entities_between_conds_across_tissues <- function(type, tissue_str, cond_str, genes_type) {
  # Get tissues in vector
  tissues <- tissue_str %>%
    str_split(pattern=",") %>%
    unlist()
  
  # Get intermediat tissue pattern
  tissue_pattern <- tissues %>%
    paste(collapse = "|")
  
  # get final tissue pattern
  tissue_pattern_final <- str_glue(".*analysis.*[{tissue_pattern}].*_results")
  
  contrasts_pattern <- cond_str %>%
    str_split(pattern=",") %>%
    unlist() %>%
    paste(collapse = "|")
  print("Conditions pattern")
  print(contrasts_pattern)
  
  
  list_files_contrasts_pattern <- str_glue("{type}_({contrasts_pattern}).tsv") %>% as.character()
  
  extract_pattern <- str_glue("{type}_({contrasts_pattern}).tsv") %>% as.character()
  
  data_dirs <- list.files(pattern=tissue_pattern_final) %>%
    set_names(tissues)
  
  
  print(data_dirs)
  tables_per_tissue_list <- map(data_dirs,
                                
                                ~{
                                  print("list_files_contrasts_pattern:")
                                  print(list_files_contrasts_pattern)
                                  
                                  result_files <- list.files(path = .x,
                                                             pattern=list_files_contrasts_pattern,
                                                             full.names = T,)
                                  
                                  unformatted_names <-list.files(path = .x,
                                                                 pattern=list_files_contrasts_pattern,
                                                                 full.names = F)
                                  
                                  print("Unformatted names:")
                                  print(unformatted_names)
                                  
                                  print("Extract pattern:")
                                  print(extract_pattern)
                                  
                                  formatted_names<- map(unformatted_names,
                                                        ~str_match(string = .x,
                                                                   pattern =extract_pattern )[[2]])
                                  
                                  print("Contrasts:")
                                  print(contrasts)
                                  
                                  
                                  tables <- map(.x = result_files,
                                                .f = ~get_sig_entities_from_path(path = .x, type = type, genes_type =genes_type))
                                  
                                  
                                  names(tables) <- formatted_names
                                  
                                  return(tables)}
  )
  return(tables_per_tissue_list)
}
#' Gets a nested list of significant entities (genes or GO-terms) between tissues across contrasts.
#' Assumes there is a data dir with each tissue described in its name, with GO and DE-tables in them
#' with corresponding names.
#' @param type Either "DE" or "GO" depending on the type of entity
#' @param cond_str Character vector: Comma-separated list of contrasts as present in the data directory, i.e. Group2,Group3
#' @param tissue_str Character vector: Comma-seperated list of tissues as present in the project root, i.e. Hyp,NAc
#' @param genes_type Either ENSEMBL or NCBI
#' @return A nested list, where the upper levels corresponds to supplied tissues,
#' and the lower levels each contain a vector of significant gene names or GO-terms for the given condition
#' @export
get_sig_entities_between_tissues_across_conds <- function(type, cond_str, tissue_str, genes_type) {
  
  
  contrasts <- cond_str %>%
    str_split(pattern=",") %>%
    unlist()
  
  tissue_pattern <- tissue_str %>%
    str_split(pattern=",") %>%
    unlist() %>%
    paste(collapse = "|")
  
  
  condition_pattern <-  str_glue("{type}_{contrasts}.tsv") %>%
    as.vector() %>%
    set_names(contrasts)
  
  all_result_files <- map(condition_pattern, ~list.files(pattern=.x, recursive = T))
  
  tables_per_condition_list <- map(all_result_files,
                                   ~{tables <- map(.x, ~get_sig_entities_from_path(path = .x, type = type, genes_type = genes_type))
                                   
                                   tissue_names <- map(.x,
                                                       ~str_match(.x,
                                                                  pattern = str_glue("analysis_({tissue_pattern})_.*"))[[2]])
                                   
                                   names(tables) <- tissue_names
                                   return(tables)
                                   })
  return(tables_per_condition_list)
}

#' Function to get the intersection values correpsodning to a UpsetR-plot
#'
#' @param The binary membership df, can be geneated from named list via fromList()
#' @param ... The unique combination of intersections to get the result for
#'
#' @return
#' @export
get_intersect_members <- function (x, ...){
  require(dplyr)
  require(tibble)
  x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
  n <- names(x)
  x %>% rownames_to_column() -> x
  l <- c(...)
  a <- intersect(names(x), l)
  ar <- vector('list',length(n)+1)
  ar[[1]] <- x
  i=2
  for (item in n) {
    if (item %in% a){
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '>= 1')
        i <- i + 1
      }
    } else {
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '== 0')
        i <- i + 1
      }
    }
  }
  do.call(filter_, ar) %>% column_to_rownames() -> x
  return(x)
}

#' Small re-write of the UpsetR::fromList function.
#'
#' @param input
#'
#' @return
#' @export
fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

#' Get intersection values from a named list, corresponding to the Upset-plot
#'
#' @param listInput Named list which can also be used in Upset fromList() function
#' @param sort whether to sort the output,
#'
#' @return Returns the labels of the intersecting values
#' @export
overlapGroups <- function (listInput, sort = TRUE) {
  
  listInputmat    <- fromList(listInput) == 1
  listInputunique <- unique(listInputmat)
  grouplist <- list()
  # going through all unique combinations and collect elements for each in a list
  for (i in 1:nrow(listInputunique)) {
    currentRow <- listInputunique[i,]
    myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
    attr(myelements, "groups") <- currentRow
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
  }
  if (sort) {
    grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  # return(map(grouplist, names))
  return(grouplist)
  # save element list to facilitate access using an index in case rownames are not named
}

#' Get specific intersection
#'
#' @param named_list A named list of entities
#' @param extract_vector A list of extrcation vectors, i.e. a list of named boolean vectors,
#' where the names correspond to the given set, and the boolean whether it should be in that set or not
#' @return
#' @export
extract_specific_intersection <- function(named_list, extract_vectors) {
  
  overlapping_groups <- overlapGroups(named_list)
  
  across_vecs <- map(extract_vectors,
                     function(bool_vector) {
                       idx <- map(overlapping_groups,
                                  ~identical(attr(.x, which = "groups"), bool_vector)) %>%
                         unlist() %>%
                         which()
                       
                       return(overlapping_groups[[idx]] %>% names())
                     }
  )
  
  return(across_vecs)
}



# Plots -------------------------------------------------------------------

#' Draws a Venn diagram of significant entities (genes or GO-terms) between tissues across contrasts,
#' or between contrasts across tissues
#'
#' @param between_contrasts boolean: Whether to plot overlaps between contrasts (TRUE) or between tissues (FALSE)
#' @param type character vector: Either "DE" or "GO" depending on the type of entity
#' @param cond_str character vector: Comma-separated list of contrasts as present in the data directory, i.e. Group2,Group3
#' @param tissue_str character vector: Comma-seperated list of tissues as present in the project root, i.e. Hyp, NAc
#'
#' @return A Venn Diagram of the ggplot-kind.
#' @export
make_venn <- function(between_contrasts, type, tissue_str, cond_str, genes_type, data_dir) {
  
  entity_type <- case_when(type=="DE"~"genes",
                           type=="GO"~ "GO-terms")
  
  if (between_contrasts) {
    tables_pr_tissue_or_cond_list <- get_sig_entities_between_conds_across_tissues(type = type , tissue_str = tissue_str, cond_str = cond_str, genes_type = genes_type)
    title <- "Significant {entity_type} (FDR<0.05) between contrasts in {.y}"
  }
  else {
    tables_pr_tissue_or_cond_list <-  get_sig_entities_between_tissues_across_conds(type = type, tissue_str=tissue_str, cond_str=cond_str, genes_type=genes_type)
    title <- "Significant {entity_type} (FDR<0.05) between tissues in {.y}"
  }
  print(tables_pr_tissue_or_cond_list)
  
  plots <- imap(tables_pr_tissue_or_cond_list,
                ~{plot <- ggVennDiagram(.x, label_geom="label") +
                  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
                  theme(legend.position = "none") +
                  theme(plot.title = element_text(face = "bold")) +
                  ggtitle(str_glue(title)) +
                  scale_x_continuous(expand = expansion(mult = .2))
                
                process_region_data(Venn(.x)) %>%
                  mutate(item=map(item, ~paste(.x, collapse=",")) %>% unlist()) %>%
                  select(-component) %>%
                  rename(shared_items=item,
                         intersection=name) %>%
                  relocate(intersection) %>%
                  write_tsv(str_glue("{data_dir}/{type}_across_{.y}.tsv"))
                return(plot)})
  return(plots)
}


# **** CHANGED: Added documentation & export, and IS NOW INTERACTIVE  ****
#' Make a volcano plot using ggplot
#'
#' @param df gene results from EdgeR or voom
#' @param genes character vector, "all" for plotting all genes,
#' otherwise only plot vector of genes
#' @param only_sig Whether or not only to print FDR significant genes
#' @param interactive_plot Whether or not to print a plotly interactive plot
#' @return The ggplot object
#' @export
ggplot_volcano_updated <- function(df,
                                   plot_genes="all",
                                   highlight_genes=c("interesting_gene_name_1", "interesting_gene_name_2"),
                                   only_FDR_sig = FALSE,
                                   only_nom_sig = FALSE,
                                   interactive_plot=TRUE,
                                   title="") {
  if(!identical(plot_genes,"all")) {
    df <- df %>% filter(Name %in% plot_genes)
  }
  if (only_nom_sig) {
    df <- df %>% filter(if_any(any_of(c("PValue", "P.Value")), ~.x<0.05))
  }
  
  if (only_FDR_sig) {
    df <- df %>% filter(if_any(any_of(c("FDR","adj.P.Val")), ~.x<0.05))
  }
  
  allLogFc <- df %>% pull("logFC")
  minLogFc <- allLogFc %>% min()
  maxLogFc <- allLogFc %>% max()
  
  minPval <- df %>%
    select(any_of(c("PValue", "P.Value"))) %>%
    pull() %>%
    log10() %>%
    `*`(-1) %>%
    max()
  
  maxPval <- 0
  
  df_formatted <- df %>%
    mutate(across(any_of(c("PValue", "P.Value")),
                  .fns = ~-log10(.x) %>% signif(3),
                  .names="log10Pval")) %>%
    mutate(across(any_of(c("FDR","adj.P.Val")),
                  .fns = ~.x < 0.05,
                  .names="FDR<0.05")) %>%
    mutate(highlighted_genes = Name %in% highlight_genes,
           FC=signif(2^(logFC),3),
           logFC=signif(logFC, 3)) %>%
    mutate(custom_col=case_when(highlighted_genes & `FDR<0.05` ~ "dark blue",
                                highlighted_genes & (!`FDR<0.05`) ~ "light blue",
                                !highlighted_genes & `FDR<0.05` ~ "red",
                                !highlighted_genes & (!`FDR<0.05`) ~ "orange"))
  
  plot <- ggplot(df_formatted, aes_string(x = "logFC",
                                          y = "log10Pval",
                                          colour = "custom_col",#"FDR<0.05",
                                          text="Name",
                                          alpha="highlighted_genes",
                                          label="FC")) +
    geom_point(size = 0.1) +
    scale_color_manual(values = c("dark blue"="dark blue",
                                  "light blue"="light blue",
                                  "red"="red",
                                  "orange"="orange"))+
    # `FALSE` = "black"),
    # name = "Significance",
    # labels = c(`TRUE` = "adj.P.Val < 0.05",
    # `FALSE` = "adj.P.Val ≥ 0.05"),
    # breaks = c("TRUE", "FALSE")) +
    scale_alpha_manual(values = c(`TRUE` = 0.7,
                                  `FALSE` = 0.7)) +
    ylab("-log10 P-value") +
    xlab("Log2 Fold Change") +
    scale_y_continuous(expand = expansion(c(0,0.5)),
                       limits =c(0,minPval)) +
    geom_hline(yintercept = -log10(0.05), lty="dashed", col="grey") +
    geom_vline(xintercept = 0, lty="dashed", col="grey") +
    theme(legend.position = "none") +
    ggtitle(title)
  
  if(interactive_plot) {
    plot <- ggplotly(plot, tooltip = c("text", "x", "y", "colour", "label"))
  }
  
  return(plot)
}

# *** CORRECTED: Can now take multiple color cols...***
#' Plot samples in two-dimensional space using MDS
#'
#' @param y The DGE list object
#' @param dims The dimensions to plot as a vector
#' @param color_by The vector to colour the plots by
#' @return Returns a list of ggplots coloured according to color_by
#' @export
ggplot_mds_repel_updated <- function (y, dims, color_by) {
  
  plotMDS_obj <- edgeR::plotMDS.DGEList(y, dim.plot = dims,
                                        plot = FALSE)
  x_y_data <- plotMDS_obj$eigen.vectors[, dims] %>% as.data.frame()
  x_y_data <- cbind(x_y_data, y$samples)
  
  var_explained_per_dim <- plotMDS_obj$var.explained[dims] %>%
    signif(2) %>% `*`(100)
  
  axis_labels <- plotMDS_obj$axislabel
  
  
  plots <- map(color_by,
               ~x_y_data %>%
                 ggplot(aes_string(x = colnames(x_y_data)[1],
                                   y = colnames(x_y_data)[2],
                                   text="Sample.ID",
                                   colour = .x)) +
                 ggplot2::geom_point() +
                 xlab(str_glue("{axis_labels}. {dims[1]} ({var_explained_per_dim[1]} % var. explained)")) +
                 ylab(str_glue("{axis_labels}. {dims[2]} ({var_explained_per_dim[2]} % var. explained)")) +
                 ggtitle(str_glue("MDS-plot colored by {.x}. Dimensions: {dims[1]} & {dims[2]}")))
  
  interactive_plots <- map(plots, ~ggplotly(.x, tooltip = c("text", "colour","x","y")))
  
  return(interactive_plots)
}

# *** ADDED: Now with added cols option ****
#' Plot MD-figures for all samples in DGElist.
#'
#' @param object The DGElist object with transcript/gene counts and sample information
#' encoded in the y$sample object
#' @param samples The samples to plot as a character vector
#' @param ncol Number of columns in the combined plot as integer vector
#' @return A facetted ggplot
#' @export
ggplot_MD_updated <-function(object, samples="all", ncol=3) {
  MD_data <- get_MD_data_for_all_samples_updated(object)
  
  if (!identical(samples, "all")) {
    MD_data <- MD_data %>%
      dplyr::filter(Sample %in% as.character(samples))
  }
  
  MD_data %>%
    ggplot(aes(x=Mean, y=Diff, col=Diff>0)) +
    geom_point(size=1, alpha=0.5) +
    facet_wrap(~Sample,ncol = ncol) +
    xlab("Average log CPM (this sample and others)") +
    ylab("log-ratio (this sample vs others)") +
    theme(legend.position = "none")
}

# Tables ------------------------------------------------------------------
# *** CHANGED: NOW EXPECTS BOTH GENE NAME AND ENSEMBL_ID, AND PASSES FURTHER ARGUMENTS TO DT***
#' Return a nicely formatted DT::datatable object of the DE results
#'
#' @param df The dataframe with the DE results
#' @param type Whether the results is edgeR or voom
#'
#' @return A DT::datatable
#' @export
get_DE_datatable_updated <- function(df, ...) {
  my_datatable <- df %>%
    select(any_of(c("Name", "ENSEMBL_ID", "description", "logFC", "PValue","FDR")))  %>%
    mutate(across(any_of(c("logFC", "PValue","FDR")),
                  ~signif(.x,2))) %>%
    DT::datatable(extensions = 'Buttons',
                  # filter="top",
                  rownames = FALSE,
                  options = list(dom = 'Bfrtip',
                                 buttons = c('copy', 'csv', 'excel', 'colvis'),
                                 pageLength = 10),
                  height = 700,
                  width = "100%",
                  ...)
  
  return(my_datatable)
}

#' *** CANGED: NOW PASSES FURTHER ARGUMENTS TO DATATABLE, AND
#' USES SET PIXEL HEIGHT TO AVOID PROBLEMS IN LOOPS ****
#' Return a nicely formatted DT::datatable of the GO-results
#'
#' @param df The GO-results
#' @param ... extra arguments passed to DT::datatable()
#'
#' @return A DT::datatable
#' @export
#'
#' @examples
#' 
#' get_GO_datatable()
get_GO_datatable_updated <- function(df, ...) {
  df %>%
    transmute(ID,
              Direction=as.factor(Direction),
              TERM,
              `#genes`=NGenes,
              PValue=signif(PValue,2),
              FDR=signif(FDR, 2)) %>%
    DT::datatable(#extensions = 'Buttons',
      # filter="top",
      rownames = FALSE,
      options = list(dom = 'Bfrtip',
                     buttons = c('copy', 'csv', 'excel', 'colvis'),
                     pageLength = 10),
      height = 700,
      width = "100%",
      ...)
}


# Render functions --------------------------------------------------------

#' Render a single tissue
#' @param markdown_path The absolute path to the markdown fle
#' @param tissue Character vector of length 1, indicating the tissue to keep,
#' as present in the metadata sheet for the project
#' @param exclude Character vector of length 1, indicating the Sample IDS to
#' remove seperated by commas, as present in the metadata sheet for the project
#' @return Returns nothing, but creates the corresponding .html files with descriptive names
#' @export
render_markdown_file_single_tissue_updated <- function (markdown_path, tissue, exclude, suffix=""){
  
  script_name <- markdown_path %>%
    basename() %>%
    tools::file_path_sans_ext()
  
  date <- format(Sys.time(), "%Y-%m-%d-%H.%M")
  
  exclude_string <- case_when(!exclude == "none" ~str_glue("Sample{exclude}_excluded"),
                              TRUE ~ "all_samples")
  
  report_name <- stringr::str_glue("{script_name}_{tissue}_{exclude_string}{suffix}_{date}.html")
  
  rmarkdown::render(markdown_path,
                    output_file = report_name,
                    params = list(tissue = tissue,
                                  exclude = exclude,
                                  report_name = report_name),
                    envir = parent.frame())
}

render_tissues_updated <- function (markdown_path, tissue_exclusion_vec, suffix="")
{
  for (i in seq_along(tissue_exclusion_vec)) {
    render_markdown_file_single_tissue_updated(markdown_path, tissue = names(tissue_exclusion_vec[i]),
                                               exclude = tissue_exclusion_vec[i], suffix)
  }
}


render_QCs <- function(markdown_path = "QC.Rmd",
                       tissue_exclusion_vec = c("fat"="none","muscle"="none"),
                       suffix="_BIGTT_SI_in_high_WH") {
  render_tissues_updated(markdown_path, tissue_exclusion_vec, suffix)
}

render_analyses <- function(markdown_path = "analysis.Rmd",
                            tissue_exclusion_vec = c("fat"="none","muscle"="none"),
                            suffix="_BIGTT_SI_in_high_WH") {
  render_tissues_updated(markdown_path, tissue_exclusion_vec, suffix)
}

render_both <- function(){
  render_QCs()
  render_analyses()
}