# Experimental planning (TENX batching and sequencing) -------------------
#' Expects a tibble of libraries to seq with following columns: ATAC/GEX_extra_mil reads, TENX_reaction_ID, ATAC/GEX_index
#'
#' @param type ATAC or GEX
#' @param flow_cell_count numeric
#' @param mil_read_limits vector of read capacity of flow cells in millions corresponding with length equal to flow_cell_count
#' @param libraries_to_seq_tibble_in  tibble of libraries to seq with following columns: ATAC/GEX_extra_mil reads, TENX_reaction_ID, ATAC/GEX_index
#' @param names flow cell names/identifiers
#' @param n_attempts how many attempts to solve the puzzle
#' @param verbose verbosity TRUE/FALSE
#'
#' @return Combined tible with flow cell identifier
#' @export
assign_reactions_to_flow_cells_reaction_centric <- function(type, flow_cell_count, mil_read_limits, libraries_to_seq_tibble_in, names, n_attempts=50, verbose=TRUE) {
  
  attempt_no <- 1
  while (attempt_no <= n_attempts){ # loop over attempts
    # attempt setup
    ## sort list so that largest is assigned first, allowing to more efficient filling of flow cells
    remaining_libraries_to_seq_tibble <- libraries_to_seq_tibble_in %>%
      arrange(desc(.data[[paste0(type,"_extra_mil_reads")]]))
    
    ## clear available flow cells and assigned reactions
    flow_cell_library_list <- map(1:flow_cell_count, ~tibble()) %>% purrr::set_names(names) # list of tibbles where actual libraries are assigned as entries in the tibbles
    available_flow_cells <- 1:flow_cell_count # set all flow cells as availble
    assigned_flow_cell_reads <- rep(0, flow_cell_count) # count the number of reads assigned to each flow cell
    
    if(verbose){print(paste0("Attempt number ", attempt_no, "..."))}
    
    for (i in 1:nrow(remaining_libraries_to_seq_tibble)) { # loop over libraries, adding to random flow cells
      
      library_sample <- slice_head(remaining_libraries_to_seq_tibble, n=1) # draw top library
      library_sample_ID <- library_sample["TENX_reaction_ID"] %>% pull()
      shuffled_available_flow_cells <- sample(available_flow_cells)
      for (flow_cell_choice in shuffled_available_flow_cells) {
        new_flow_cell_reads <- assigned_flow_cell_reads[[flow_cell_choice]] + pull(library_sample, .data[[paste0(type, "_extra_mil_reads")]])
        
        if (nrow(flow_cell_library_list[[flow_cell_choice]])==0) { # do not check stuff since it's first library on flow cell - just add
          flow_cell_library_list[[flow_cell_choice]] <- bind_rows(flow_cell_library_list[[flow_cell_choice]], library_sample) #add to flow cell list
          if(verbose){cat("\t\tAssigned ",library_sample_ID," to flowcell no. ", flow_cell_choice,"...\n")}
          remaining_libraries_to_seq_tibble <- remaining_libraries_to_seq_tibble[-1, ] # remove top library
          assigned_flow_cell_reads[[flow_cell_choice]] <- new_flow_cell_reads
          
          if (nrow(remaining_libraries_to_seq_tibble)==0) { # check if we are done - no more libraries to add
            flow_cell_library_list_combined <- flow_cell_library_list %>%
              bind_rows(.id="flow_cell") %>%
              group_by(flow_cell) %>%
              arrange(flow_cell, `TENX_reaction_ID`) %>%
              mutate(cumm_sum_flow_cell=cumsum(GEX_extra_mil_reads))
            
            print(paste0("Succeeeded after ", attempt_no, " attempts."))
            return(flow_cell_library_list_combined)
          }
          break
        }
        else{ # check conditions
          current_indices <- unique(flow_cell_library_list[[flow_cell_choice]][paste0(type, "_index")])
          index_clash <- pull(library_sample, .data[[paste0(type, "_index")]]) %in% current_indices
          
          if (new_flow_cell_reads<(mil_read_limits[flow_cell_choice]) & (!index_clash))   { # add library if flow-cell capacity not exceeded and no index clash
            flow_cell_library_list[[flow_cell_choice]] <- bind_rows(flow_cell_library_list[[flow_cell_choice]], library_sample) #add to flow cell list
            if(verbose){cat("\t\tAssigned ",library_sample_ID," to flowcell no. ", flow_cell_choice,"...\n")}
            remaining_libraries_to_seq_tibble <- remaining_libraries_to_seq_tibble[-1, ] # remove top library
            assigned_flow_cell_reads[[flow_cell_choice]] <- new_flow_cell_reads
            
            if (nrow(remaining_libraries_to_seq_tibble)==0) { # check if we are done - no more libraries to add
              flow_cell_library_list_combined <- flow_cell_library_list %>%
                bind_rows(.id="flow_cell") %>%
                group_by(flow_cell) %>%
                arrange(flow_cell, `TENX_reaction_ID`) %>%
                mutate(cumm_sum_flow_cell=cumsum(GEX_extra_mil_reads))
              
              print(paste0("Succeeeded after ", attempt_no, " attempts."))
              return(flow_cell_library_list_combined)
            }
            break # library assigned succesfully. Break out of flow cell loop
          }
          else {
            if(verbose){cat("\t\tAssign of ", library_sample_ID, 'to flowcell no. ', flow_cell_choice, " failed (index clash or capacity limit). Trying next flow cell\n.")}
          }
        }
      }
    }
    attempt_no <- attempt_no+1
  }
  return(stop(paste0("Failed after ", n_attempts, " attempts." )))
}

assign_reactions_to_flow_cells_flow_cell_centric <- function(type, flow_cell_count, mil_read_limits, libraries_to_seq_tibble_in, names, n_attempts=50) {
  
  # setup
  flow_cell_library_list <- map(1:flow_cell_count, ~tibble()) %>% purrr::set_names(names)
  
  attempt_no <- 1
  while (attempt_no <= n_attempts){ # loop over attempts
    libraries_to_seq_tibble <- libraries_to_seq_tibble_in
    print(paste0("Attempt number ", attempt_no, "..."))
    
    for (i in 1:flow_cell_count) { # loop over flow cells
      # for (i in 1:nrow(libraries_to_seq_tibble_in)) { # loop over flow cells
      cat("\tFilling up flowcell number ", i, "\n")
      current_reads <- 0
      
      while(TRUE) { # loop to add libraries until break or return condition
        index_clash=FALSE
        library_sample_idx <- sample(1:nrow(libraries_to_seq_tibble), 1) # draw one library at random
        library_sample_ID <- libraries_to_seq_tibble[library_sample_idx, ]["TENX_reaction_ID"] %>% pull()
        cat("\t\tTrying ",library_sample_ID,"...\n")
        new_reads <- current_reads + pull(libraries_to_seq_tibble[library_sample_idx, ][paste0(type, "_extra_mil_reads")])
        
        if (new_reads<(mil_read_limits[i])) {
          if (!nrow(flow_cell_library_list[[i]])==0) { # check for index clash
            current_indices <- unique(flow_cell_library_list[[i]][paste0(type, "_index")])
            index_clash <- libraries_to_seq_tibble[library_sample_idx, ][paste0(type, "_index")] %in% current_indices
          }
          if (!index_clash) { # successfully add library to flow-cell
            flow_cell_library_list[[i]] <- bind_rows(flow_cell_library_list[[i]], libraries_to_seq_tibble[library_sample_idx, ]) #add to flow cell list
            libraries_to_seq_tibble <- libraries_to_seq_tibble[-library_sample_idx, ] # remove from available flow cell list
            current_reads <- new_reads
            
            if (nrow(libraries_to_seq_tibble)==0) { # success - no more libraries to add
              flow_cell_library_list_combined <- flow_cell_library_list %>% bind_rows(.id="flow cell no.")
              print(paste0("Succeeeded after ", attempt_no, " attempts."))
              return(flow_cell_library_list_combined)
            }
          }
          else {
            print("Index clash for library, trying a new sample...")
          }
        }
        else {
          print("Flow-cell went over the limit. Skipped to next one...")
          break
        }
      }
    }
    attempt_no <- attempt_no+1
  }
  return(stop("failed after many tries"))
}

#' Convenience function for demultiplexing by genotype that generates needed
#' .tsv files of each sample in each 10 reaction specified on the larger table.
get_and_write_sample_tsvs_for_multiplexd_10x_reactions <- function(chip_samples_path="/Users/tqb695/GitHub_repositories/SEGMENT/data-raw/samples_present_in_chip.tsv",
                                                                   flinc_translation_excel_sheet="/Users/tqb695/GitHub_repositories/SEGMENT/data-raw/FLINC_all_IDs_samples_31052023_v2_removed_extra_hyphen_1013_8309.xlsx",
                                                                   output_dir="samples", remove_NA=T){
  
  # Load IDs that were not found in the various files
  chip_samples <- read_tsv(file = chip_samples_path, col_names = "ID")
  
  # prepare translation table from FLINC ID database
  translation_table <- readxl::read_xlsx(path = flinc_translation_excel_sheet) %>%
    mutate(genotype_ID=paste0("118x", id)) %>%
    dplyr::filter(genotype_ID %in% chip_samples$ID) %>%
    distinct(record_id, .keep_all = T) %>%
    dplyr::select(subject_ID=record_id, genotype_ID)
  
  # get data base info
  joined_wet_lab_table <- get_joined_wetlab_table() %>%
    dplyr::filter(SCOP_ID=="SCOP_2023_0238") %>%
    dplyr::filter(str_detect(`TENX_reaction_ID`, "multiplex|reaction")) %>%
    left_join(translation_table)
  
  reactions <- joined_wet_lab_table$`TENX_reaction_ID` %>% unique()
  
  map(reactions,
      ~dplyr::filter(joined_wet_lab_table,,
                     `TENX_reaction_ID`==.x) %>%
        distinct(subject_ID, .keep_all = T) %>%
        # distinct(multiplex_batch, .keep_all=T) %>%
        dplyr::select(genotype_ID) %>%
        write_tsv(file=paste0(output_dir,
                              "/",
                              .x,
                              ".tsv"),
                  col_names = F))
}

generate_10x_batches <- function(batch_names = c("D","E","F","G",
                                                 "H","I","J","K"),
                                 phenotype_tibble,
                                 already_sequenced_donor_ids,
                                 random = TRUE,
                                 batch_sizes = c(8,8,8,8,
                                                 8,8,8,8),
                                 extra_batch = TRUE,
                                 strat_col_1_vals =c("Male","Female"),
                                 strat_col_2_vals =c("IR","IS")) {
  
  
  choose_seq_from_df <- phenotype_tibble %>%
    ungroup() %>%
    dplyr::filter(!donor_id %in% already_sequenced_donor_ids) %>%
    arrange(index_Matsuda) %>%
    relocate(group,
             sex_chars,
             index_Matsuda,
             BMI,
             age,
             .after = donor_id) %>%
    mutate(strat=paste(group, sex_chars, sep = "-"), .after=donor_id)
  
  to_seq_per_batch <- list()
  
  for (batch_no in batch_names)  {
    batch_size <- batch_sizes[batch_no]
    batch_name <- batch_names[batch_no]
    batch_tibble <- tibble()
    
    # for each batch we choose:
    for (sex in strat_col_1_vals) { # from each sex
      for (type in strat_col_2_vals) { # top and bottom
        
        try <- 0
        n <- 0
        n_select <- case_when(sex=="Male" ~ 1,
                              sex=="Female" & !(batch_no==9) ~ 3,
                              sex=="Female" & batch_no==9 ~ 1)
        
        while (n<n_select) {  # twice
          
          # get selection IS or IR
          sex_specific_obs <- dplyr::filter(choose_seq_from_df, sex_chars==sex)
          
          if (random) {
            slice_obs  <- dplyr::filter(sex_specific_obs, group==type) %>% slice_sample(n=1)
            
          } else{
            
            n_obs <- nrow(sex_specific_obs)
            
            slice_index <- case_when(type=="IS"~(1+try),
                                     type=="IR"~(n_obs-try))
            
            slice_obs <- slice(sex_specific_obs, slice_index)
          }
          
          
          if (!slice_obs$family_id_fixed %in% batch_tibble$family_id_fixed) {  # IF SUCCESS
            
            batch_tibble <- bind_rows(batch_tibble, slice_obs)
            
            # remove from available options by overwriting
            choose_seq_from_df <- dplyr::filter(choose_seq_from_df, !donor_id %in% batch_tibble$donor_id)
            
            n <- n+1
          }
          
          else {  # IF NOT SUCCESS, TRY AGAIN
            message("Failed.. trying new...")
            try <- try+1
          }
        }
      }
    }
    to_seq_per_batch[[batch_name]] <- batch_tibble
  }
  # the remainder
  to_seq_per_batch[["L"]] <- choose_seq_from_df
  to_seq_per_batch
}

get_seurat_object_sizes <- function(SO) {
  slot_names <- slotNames(SO) %>% purrr::set_names()
  assay_names <- SO@assays %>% names() %>% purrr::set_names()
  
  # General
  purrr::map(slot_names, ~object.size(slot(SO,.x)) %>% format(units='auto')) %>% unlist() %>% print()
  
  # Assays
  purrr::map(assay_names, ~object.size(SO[[.x]]) %>% format(units='auto')) %>% unlist() %>% print()
  
  # ATAC
  purrr::map(slotNames(SO[["ATAC"]]) %>% purrr::set_names(), ~slot(SO[["ATAC"]], .x) %>% object.size() %>% format(units="auto")) %>% unlist() %>% print()
  
}