cal_div_table_check <- function(table_list, table_name, col_names) {

  # Check if initial function conditions are met ----
  # table list conditions
  if(!is.list(table_list)) stop("table_list argument must be a LIST with named indices.")
  if(is.null(names(table_list))) stop("table_list argument must be a list with NAMED indices.")
  # table name conditions
  if(!is_scalar_character(table_name)) stop("{{table_name}} argument must be a SCALAR CHARACTER.")
  # Column names conditions
  if(!is.character(col_names)) stop("{{col_names}} argument must be a CHARACTER data type.")

  # Check if table is in table list ----
  # If there is only one table cal div table check fails
  # There is only the parameter metadata in this div
  if(length(table_list) <= 1) return(FALSE)
  # Check if table name is in table list
  if(!pluck_exists(table_list[[table_name]])) return(FALSE)

  # Check if the table has any data ----
  table <- table_list[[table_name]]
  if(nrow(table) == 0) return(FALSE)

  # Check if column names are in the table ----
  table <- table %>%
    pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names)
  col_name_check <- all(col_names %in% names(table))

  # If both table check and column name check pass return TRUE, else return FALSE
  return(col_name_check)

}
