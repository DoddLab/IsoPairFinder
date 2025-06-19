#' @title diff_analysis
#' @author Zhiwei Zhou
#' @param peak_table peak table
#' @param sample_info sample information table
#' @param control_group control_group label. Default: "WT"
#' @param case_group case_group label. Default: "hyuA"
#' @param tracer_condition tracer condition, either '12C' or '13C'. Default: c('12C', '13C')
#' @param p_value_cutoff Default: 0.05
#' @param p_adjust Default: TRUE
#' @param fold_change_cutoff Default: 20
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' diff_analysis(peak_table = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/C12/converted_peak_table.csv',
#'               sample_info = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/sample_info.xlsx',
#'               control_group = "WT",
#'               case_group = 'hyuA',
#'               tracer_condition = '12C',
#'               p_value_cutoff = 0.05,
#'               p_adjust = TRUE,
#'               fold_change_cutoff = 20)
#' }



# peak_table <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/C12/converted_peak_table.csv'
# sample_info <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/sample_info.xlsx'
# tracer_condition <- '12C'

# diff_analysis(peak_table = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/C12/converted_peak_table.csv',
#               sample_info = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/sample_info.xlsx',
#               control_group = "WT",
#               case_group = 'hyuA',
#               tracer_condition = '12C',
#               p_value_cutoff = 0.05,
#               p_adjust = TRUE,
#               fold_change_cutoff = 20)

diff_analysis <- function(peak_table,
                          sample_info,
                          # path = '.',
                          control_group = "WT",
                          case_group = 'hyuA',
                          tracer_condition = c('12C', '13C'),
                          p_value_cutoff = 0.05,
                          p_adjust = TRUE,
                          fold_change_cutoff = 20,
                          abundance_cutoff = 10^4) {
  # browser()
  tracer_condition <- match.arg(tracer_condition)

  # check samples whether included in the sample info
  temp_col_name <- colnames(peak_table)
  sample_id_control <- sample_info %>%
    dplyr::filter(type == 'sample') %>%
    dplyr::filter(tracer_group == tracer_condition) %>%
    dplyr::filter(group == control_group) %>%
    dplyr::pull(sample_id)

  sample_id_case <- sample_info %>%
    dplyr::filter(type == 'sample') %>%
    dplyr::filter(tracer_group == tracer_condition) %>%
    dplyr::filter(group == case_group) %>%
    dplyr::pull(sample_id)

  # extract peak abundance as matrix
  temp_sample_name <- c(sample_id_control, sample_id_case)

  temp_data <- peak_table %>%
    dplyr::select(dplyr::all_of(temp_sample_name))
  idx_control <- match(sample_id_control, colnames(temp_data))
  idx_case <- match(sample_id_case, colnames(temp_data))
  p_values <- apply(temp_data, 1, function(x){
    result <- t.test(x[idx_control], x[idx_case], alternative = 'two.sided', var.equal = TRUE)
    result$p.value
  })
  q_values <- p.adjust(p_values, method = 'BH')

  fold_change <- apply(temp_data, 1, function(x){
    mean(x[idx_case])/mean(x[idx_control])
  })

  # average_abundance = apply(temp_data, 1, mean)
  average_abundance_control <- apply(temp_data, 1, function(x){mean(x[idx_control])})
  average_abundance_case <- apply(temp_data, 1, function(x){mean(x[idx_case])})

  result_table <- peak_table %>%
    dplyr::mutate(p_values = p_values,
                  q_values = q_values,
                  fold_change = fold_change,
                  avg_abundance_control = average_abundance_control,
                  avg_abundance_case = average_abundance_case)


  if (p_adjust) {
    result_table <- result_table %>%
      dplyr::mutate(increase_label = dplyr::case_when(
        q_values <= p_value_cutoff & fold_change >= fold_change_cutoff ~ 'signif_increase',
        !(q_values < p_value_cutoff & fold_change > fold_change_cutoff) ~ 'unsignif')) %>%
      dplyr::mutate(log2_fc = log2(fold_change)) %>%
      dplyr::mutate(abundance_label = dplyr::case_when(
        avg_abundance_control >= abundance_cutoff | avg_abundance_case >= abundance_cutoff ~ 'abundant',
        avg_abundance_control < abundance_cutoff & avg_abundance_case < abundance_cutoff ~ 'low_abundant'))
  } else {
    result_table <- result_table %>%
      dplyr::mutate(increase_label = dplyr::case_when(
        p_values <= p_value_cutoff & fold_change >= fold_change_cutoff ~ 'signif_increase',
        !(p_values < p_value_cutoff & fold_change > fold_change_cutoff) ~ 'unsignif')) %>%
      dplyr::mutate(log2_fc = log2(fold_change)) %>%
      dplyr::mutate(abundance_label = dplyr::case_when(
        avg_abundance_control >= abundance_cutoff | avg_abundance_case >= abundance_cutoff ~ 'abundant',
        avg_abundance_control < abundance_cutoff & avg_abundance_case < abundance_cutoff ~ 'low_abundant'))
  }


  return(result_table)

}
