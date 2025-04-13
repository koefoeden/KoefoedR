#' Auxillary function that prints a DT::datatable object correctly
#' both when working interactively and in reports
#'
#' @param DT The DT dataframe to plot
#'
#' @return Nothing, invisible object
#' @export
print_DT <- function(DT) {
  if(isTRUE(getOption('knitr.in.progress'))){
    DT %>% knitr::knit_print() %>% cat()
  }
  else {
    print(DT)
  }
  return(invisible())
}

#' Auxillary function that plots a GGplotly object correctly
#' both when working interactively and in reports
#'
#' @param plotly_plot The ggplotly-object to plot
#'
#' @return Nothing, invisible object
#' @export
#'
print_plotly <- function(plotly_plot) {
  plotly_plot %>%  htmltools::tagList() %>% print()
  
  return(invisible())
}

#' Generate appropriate tab-header strings.
#' Useful in a loop to contain each plot in its own tab.
#' @param text The text in the header
#' @param level the level of the header, use one less than parent
#' @export
catHeader <- function(text = "", level = 3) {
  
  cat(paste0("\n\n",
             paste(rep("#", level), collapse = ""),
             " ", text, "\n"))
}


#' Generate appropriate tab-header strings with .tabset attribute
#' Useful in a loop to contain each plot in its own tab.
#'
#' @param text The text in the header
#' @param level the level of the header, use one less than parent
#'
#' @export
catHeader_w_tabset <- function(text = "", level = 3) {
  
  cat(paste0("\n\n",
             paste(rep("#", level), collapse = ""),
             " ", text, " {.tabset}","\n"))
}

#' Required initialization and loading Plotly and DT libaries
#'
#' @return Nothing
#' @export
initialize_plotly_and_DT <- function() {
  
  htmltools::tagList(plotly:::plotly_empty())
  DT::datatable(matrix())
}