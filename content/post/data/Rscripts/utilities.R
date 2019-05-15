## validate ratio values to remove `NA`
validate_rows <-
  function(df) {
    apply(df, 1, function(x) all(!is.na(x)))
  }

## clean raw input
## x refers to the file path
clean_rawdata <- function(x) {
  read_tsv(x) %>%
    as_tibble() %>% 
    ## Remove the contamination and duplicates
    filter(is.na(.$`Only identified by site`) & is.na(.$Reverse) & is.na(.$`Potential contaminant`) & `Razor + unique peptides` >= 1)
}

## clean data and do calculations to obtain consistent comparison between groups 
clean_data <- function(df) {
  cat("Cleaning input dataframe... \n")
  df %>% 
    dplyr::select(order(colnames(.))) %>%
    ## Remove the gene names were not identified.
    ## In the dimethyl labelling, the comparision ratio sometimes were inverted and required inversion.
    mutate_at(vars(contains("invert")), list(~(1 / .))) %>%
    filter(!is.na(Gene_names)) %>%
    ## Some duplicate gene names were removed.
    distinct(Gene_names, .keep_all = TRUE) %>% 
    ## Remove introduced the `invert` in the column names.
    magrittr::set_colnames(sapply(colnames(.), function(x) str_replace(x, "_invert", ""))) %>% 
    dplyr::select(Gene_names, everything())
}

## remove invalid rows and keep the desired columnes
remove_invalid <- function(df) {
  df %>%
    dplyr::filter(rep_tissue_valid == TRUE) %>%
    dplyr::select(-matches("valid$"))
}

validate_df <- function(df, validation_list) {
  map(validation_list,
      function(x) {
        colname <- sym(glue("rep_{x}_valid"))
        grep_string <- ifelse(x == "tissue", "^[A-Za-z]_.+_[0-9]$", glue("_{x}$")) 
        
        df %>% mutate(!!colname := validate_rows(dplyr::select(., matches(grep_string))))
      }) %>% 
    purrr::reduce(left_join)
}