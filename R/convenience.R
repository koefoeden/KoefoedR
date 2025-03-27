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

#' Scaling function that does not return weird matrix
#' @param vec Vector to scale
#' @return Scaled vector
#' @export
scale_this <- function(vec){
  
  (vec - mean(vec, na.rm=TRUE)) / sd(vec, na.rm=TRUE)
}

#' Inverse-normal transforms a vector with NA's in it
#' @param vec The vector to normal transform
#' @return The inverse-normal transformed vector
#' @export
replace_na_with_INT <- function(vec) {
  
  NAs_idxs <- which(is.na(vec))
  if (length(NAs_idxs >0 )) {
    vec[-NAs_idxs] <- RNOmni::RankNorm(vec[-NAs_idxs])
  }
  else {
    vec <-RNOmni::RankNorm(vec)
  }
  return(vec)
}

#' Write TSV and excel file
#'
#' @param x Object to write to a TSV and excel file
#' @param path Path of the resulting files
#'
#' @return Nothing.
#' @export
write_tsv_and_excel <- function(x, path, ...) {
  write_tsv(x=x, file=paste0(path,".tsv"))
  writexl::write_xlsx(x = x, path = paste0(path,".xlsx"))
}



#'  Switch names and values for named vector
#'
#' @param named_vector The named vector to switch names and values for
#'
#' @return Named vector with switched names and values
#' @export

switch_names_and_values_for_vec <- function(named_vector) {
  my_names <- names(named_vector)
  values <- unname(named_vector)
  
  return(set_names(my_names, values))
}


gzip_and_index_w_sort <- function(frag_path) {
  # 3) Look into parsort from parallel module.
  avail_cores <- Sys.getenv('SLURM_CPUS_PER_TASK')
  out_file <- paste0(frag_path, ".sorted.gz")
  
  "export LC_ALL=C; sort --parallel={avail_cores} -k 1,1 -k2,2n {frag_path} | /maps/projects/scop/apps/bgzip -@ {avail_cores} > {out_file}" %>% str_glue() %>%  system()
  "/maps/projects/scop/apps/tabix -f -p bed {out_file}" %>% str_glue() %>% system()
  
  return(out_file)
}


prefer_std_funcs <- function() {
  conflicts_prefer(dplyr::filter,
                   base::intersect,
                   purrr::set_names,
                   dplyr::select)
}

# Git binds
get_current_commit <- function() {
  
  system('echo "$(git rev-list --count HEAD)-$(git rev-parse --short HEAD)"', intern = TRUE)
}

gcp <- function(message) {
  
  system(str_glue("git commit -a -m '{message}'"))
  system("git push")
}

gp <- function() {
  
  system(str_glue("git pull"))
}
gs <- function() {
  
  system(str_glue("git status"))
}
