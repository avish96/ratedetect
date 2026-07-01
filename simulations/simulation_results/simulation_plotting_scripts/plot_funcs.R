## Manuscript - Figures
## Plotting script
## Function to select appropriate power results array based on input value

power_select <- function(value) {
  if (value == 'mean') {
    selected_array <- pwr_arr_mean
    test_type = 1
  } else if (value == 'ratio') {
    selected_array <- pwr_arr_ratio
    test_type = 2
  } else if (value == 'tstat') {
    selected_array <- pwr_arr_tstat
    test_type = 3
  } else if (value == 'cv') {
    selected_array <- pwr_arr_cv
    test_type = 4
  } else if (value == 'KST') {
    selected_array <- pwr_arr_KST
    test_type = 5
  } else if (value == 'MWWT') {
    selected_array <- pwr_arr_MWWT
    test_type = 6
  } else if (value == 'DT') {
    selected_array <- pwr_arr_DT
    test_type = 7
  } else {
    stop("Invalid value: Please choose a valid test statistic")
  }
  return(list(selected_array,test_type))
}


