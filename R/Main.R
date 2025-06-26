################################################################################
# find_intemidates -------------------------------------------------------------
#' @title find_intemidates
#' @author Zhiwei Zhou
#' @param peak_table_unlabel unlabeled peak table
#' @param peak_table_label labeled peak table
#' @param sample_info sample information, e.g. sample_info.xlsx
#' @param path working directory. Default: "."
#' @param polarity ionization polarity, either 'positive' or 'negative'. Default: 'positive'
#' @param control_group label of control groups. e.g. c("WT")
#' @param case_group label of case groups. e.g. c('hyuA')
#' @param raw_data_folder a list of raw data folder names, which should be consist with the group name in sample info, e.g. list('hyuA_12C' = 'hyuA_12C', 'hyuA_13C' = 'hyuA_13C'). Default: NULL. The raw data should be in mzML format.
#' @param mz_tol Default: 10 ppm
#' @param rt_tol Default: 0.05 min
#' @param p_value_cutoff Default: 0.05
#' @param p_adjust whether to adjust p-values using the Benjamini-Hochberg method. Default: TRUE
#' @param fold_change_cutoff Default: 20
#' @param is_recognize_adducts whether to recognize adducts, neutral losses, and in-source fragments. Default: TRUE
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' find_intemidates(peak_table_unlabel = '230830_Csp_YgeX_WT_UA.xlsx',
#'                  peak_table_label = '230830_Csp_YgeX_WT_13CUA.xlsx',
#'                  path = '~/Project/00_Uric_Acid_project/20230831_isotope_tracing/',
#'                  control_group = c("WT_UA", "WT_13CUA"),
#'                  case_group = c('ygeX_UA', 'ygeX_13CUA'),
#'                  mz_tol = 10,
#'                  rt_tol = 0.05,
#'                  p_value_cutoff = 0.05,
#'                  fold_change_cutoff = 20)
#' }




# path <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/'
# peak_table_unlabel <- 'peak_table_C12.csv'
# peak_table_label <- 'peak_table_C13.csv'
# sample_info = 'sample_info.xlsx'
# polarity <- 'positive'
# control_group = c("WT")
# case_group = c('hyuA')
# mz_tol = 10
# rt_tol = 0.05
# p_value_cutoff = 0.05
# fold_change_cutoff = 10
# p_adjust = TRUE
# is_recognize_adducts = TRUE

# raw_data_folder <- list(
#   'hyuA_12C' = 'hyuA_12C',
#   'hyuA_13C' = 'hyuA_13C'
# )

find_intemidates <- function(peak_table_unlabel,
                             peak_table_label,
                             sample_info,
                             path = '.',
                             polarity = c('positive', 'negative'),
                             control_group = c("WT"),
                             case_group = c('hyuA'),
                             raw_data_folder = NULL,
                             mz_tol = 10,
                             rt_tol = 0.05,
                             p_value_cutoff = 0.05,
                             p_adjust = TRUE,
                             fold_change_cutoff = 20,
                             is_recognize_adducts = TRUE) {
  # browser()

  # check parameters
  polarity <- match.arg(polarity)

  if (missing(peak_table_unlabel)) {
    stop('Please provide the peak_table\n')
  }

  if (missing(peak_table_label)) {
    stop('Please provide the peak_table\n')
  }

  if (missing(sample_info)) {
    stop('Please provide the sample_info\n')
  }

  check_parameters(
    peak_table_unlabel = peak_table_unlabel,
    peak_table_label = peak_table_label,
    sample_info = sample_info,
    path = path,
    control_group = control_group,
    case_group = case_group)


  # save the parameters
  list_para <- list(
    'package_verison' = packageVersion('IsoPairFinder'),
    # assign parameter from parameter set
    'peak_table_unlabel' = peak_table_unlabel,
    'peak_table_label' = peak_table_label,
    'sample_info' = sample_info,
    'path' = path,
    'polarity' = polarity,
    'control_group' = control_group,
    'case_group' = case_group,
    'raw_data_folder' = ifelse(is.null(raw_data_folder), '', paste(raw_data_folder, collapse = ';')),
    'mz_tol' = paste(mz_tol, collapse = ';'),
    'rt_tol' = paste(rt_tol, collapse = ';'),
    'p_value_cutoff' = paste(p_value_cutoff, collapse = ';'),
    'p_adjust' = p_adjust,
    'fold_change_cutoff' = fold_change_cutoff,
    'is_recognize_adducts' = is_recognize_adducts
  ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate_all(as.character) %>%
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = 'para')

  path_output <- file.path(path, '00_tracer_result')
  dir.create(file.path(path_output), showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(list_para, path = file.path(path_output, 'para_list.txt'), append = FALSE)

  cat(crayon::bgBlue('Start processing...'), '\n')
  message(crayon::blue('1. Reading the peak tables... \n'))
  raw_data_unlabel <- readr::read_csv(file.path(path, peak_table_unlabel), show_col_types = FALSE)
  raw_data_label <- readr::read_csv(file.path(path, peak_table_label), show_col_types = FALSE)
  if (stringr::str_detect(sample_info, '\\.xlsx$')) {
    sample_info <- readxl::read_xlsx(file.path(path, sample_info))
  } else {
    sample_info <- readr::read_csv(file.path(path, sample_info), show_col_types = FALSE)
  }

  sample_id_control_12C <- sample_info %>%
    dplyr::filter(type == 'sample') %>%
    dplyr::filter(tracer_group == '12C') %>%
    dplyr::filter(group == control_group) %>%
    dplyr::pull(sample_id)

  sample_id_case_12C <- sample_info %>%
    dplyr::filter(type == 'sample') %>%
    dplyr::filter(tracer_group == '12C') %>%
    dplyr::filter(group == case_group) %>%
    dplyr::pull(sample_id)

  sample_id_control_13C <- sample_info %>%
    dplyr::filter(type == 'sample') %>%
    dplyr::filter(tracer_group == '13C') %>%
    dplyr::filter(group == control_group) %>%
    dplyr::pull(sample_id)

  sample_id_case_13C <- sample_info %>%
    dplyr::filter(type == 'sample') %>%
    dplyr::filter(tracer_group == '13C') %>%
    dplyr::filter(group == case_group) %>%
    dplyr::pull(sample_id)

  # change the rt toleracne if the unit is in seconds
  ifelse(max(raw_data_unlabel$rt) > 60, rt_tol <- rt_tol*60, rt_tol <- rt_tol)
  ifelse(max(raw_data_unlabel$rt) > 60, rt_unit <- 'seconds', rt_unit <- 'minutes')


  # 1. Differential analysis
  message(crayon::blue('2. Differential analysis... \n'))
  raw_data_unlabel <- diff_analysis(peak_table = raw_data_unlabel,
                                    sample_info = sample_info,
                                    control_group = control_group,
                                    case_group = case_group,
                                    tracer_condition = '12C',
                                    p_value_cutoff = p_value_cutoff,
                                    p_adjust = p_adjust,
                                    fold_change_cutoff = fold_change_cutoff)

  raw_data_label <- diff_analysis(peak_table = raw_data_label,
                                  sample_info = sample_info,
                                  control_group = control_group,
                                  case_group = case_group,
                                  tracer_condition = '13C',
                                  p_value_cutoff = p_value_cutoff,
                                  p_adjust = p_adjust,
                                  fold_change_cutoff = fold_change_cutoff)

  dir.create(file.path(path, '00_tracer_result', '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
  save(raw_data_unlabel, file = file.path(path, '00_tracer_result', '00_intermediate_data', 'raw_data_unlabel.RData'))
  save(raw_data_label, file = file.path(path, '00_tracer_result', '00_intermediate_data', 'raw_data_label.RData'))

  # Visualize the differential analysis results
  plot_unlabel <- plot_diff_feature(raw_data_unlabel,
                                    p_value_cutoff = p_value_cutoff,
                                    fold_change_cutoff = fold_change_cutoff,
                                    p_adjust = p_adjust,
                                    text_title = paste0(control_group, " vs ", case_group, ' (12C)'))

  ggplot2::ggsave(plot = plot_unlabel,
                  filename = file.path(path, '00_tracer_result', 'vocano_plot_unlabeled.pdf'),
                  width = 8, height = 6)

  plot_label <- plot_diff_feature(raw_data_label,
                                  p_value_cutoff = p_value_cutoff,
                                  fold_change_cutoff = fold_change_cutoff,
                                  p_adjust = p_adjust,
                                  text_title = paste0(control_group, " vs ", case_group, ' (13C)'))

  ggplot2::ggsave(plot = plot_label,
                  filename = file.path(path, '00_tracer_result', 'vocano_plot_labeled.pdf'),
                  width = 8, height = 6)


  # Summary the differential analysis results
  if (p_adjust) {
    feature_sig_unlabel <- raw_data_unlabel %>%
      dplyr::filter(q_values <= p_value_cutoff & fold_change >= fold_change_cutoff)
    feature_sig_label <- raw_data_label %>%
      dplyr::filter(q_values <= p_value_cutoff & fold_change >= fold_change_cutoff)
  } else {
    feature_sig_unlabel <- raw_data_unlabel %>%
      dplyr::filter(p_values <= p_value_cutoff & fold_change >= fold_change_cutoff)
    feature_sig_label <- raw_data_label %>%
      dplyr::filter(p_values <= p_value_cutoff & fold_change >= fold_change_cutoff)
  }

  message(crayon::yellow('Unlabeled group: There are ', nrow(feature_sig_unlabel),
                         'features enriched ( p_value <=', p_value_cutoff,
                         '& fold-change >=', fold_change_cutoff, ')\n'))
  message(crayon::yellow('Labeled group: There are ', nrow(feature_sig_label),
                         'features enriched ( p_value <=', p_value_cutoff,
                         '& fold-change >=', fold_change_cutoff, ')\n'))


  save(feature_sig_unlabel, file = file.path(path, '00_tracer_result', '00_intermediate_data', 'feature_sig_unlabel'))
  save(feature_sig_label, file = file.path(path, '00_tracer_result', '00_intermediate_data', 'feature_sig_label'))

  # print message if no features are significantly enriched
  if (nrow(feature_sig_unlabel) == 0 |
      nrow(feature_sig_label) == 0) {
    message(crayon::underline('No features are significantly enriched in the unlabeled/labeled group!\n'))

    return(NULL)
  }



  # 2. Recognize and merge adducts, neutral losses, in-source fragments
  if (is_recognize_adducts) {
    cat('\n')
    message(crayon::blue('2.1. Recognize adducts, neutral losses, in-source fragements of enriched peaks ... \n'))
    peak_table <- raw_data_unlabel %>%
      dplyr::rename('name' = id) %>%
      dplyr::select(name, mz, rt, dplyr::starts_with(sample_id_control_12C), dplyr::starts_with(sample_id_case_12C))

    peak_table_signif <- raw_data_unlabel %>%
      dplyr::filter(increase_label == 'signif_increase') %>%
      dplyr::rename('name' = id) %>%
      dplyr::select(name, mz, rt, dplyr::starts_with(sample_id_control_12C), dplyr::starts_with(sample_id_case_12C))


    dir.create(file.path(path, '00_tracer_result', '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
    save(peak_table, file = file.path(path, '00_tracer_result', '00_intermediate_data', 'peak_table.RData'))
    save(peak_table_signif, file = file.path(path, '00_tracer_result', '00_intermediate_data', 'peak_table_signif.RData'))

    # assign the raw data folder
    if (length(raw_data_folder) == 0) {
      raw_data_path_12C <- paste0(case_group, '_12C')
      raw_data_path_13C <- paste0(case_group, '_13C')
    } else {
      raw_data_path_12C <- raw_data_folder[[paste0(case_group, '_12C')]]
      raw_data_path_13C <- raw_data_folder[[paste0(case_group, '_13C')]]
    }


    list_peak_group_annotation_merge <- recognize_rela_peak(peak_table = peak_table,
                                                            peak_table_signif = peak_table_signif,
                                                            path_dir = path,
                                                            polarity = polarity,
                                                            raw_data_folder = raw_data_path_12C,
                                                            tol_mz = mz_tol,
                                                            tol_rt = rt_tol,
                                                            cutoff_ssc = -1,
                                                            cutoff_ppc = 0.8)
    rm(peak_table, peak_table_signif);gc()
    result_peak_recognization_unlabel <- generate_recoginize_table(list_peak_group_annotation_merge = list_peak_group_annotation_merge)

    temp_num1 <- nrow(feature_sig_unlabel)
    temp_num2 <- length(unique(result_peak_recognization_unlabel$base_peak))

    cat('\n')
    message(crayon::green(temp_num1 - temp_num2, 'features are merged through recognizing adducts/neutral_loss/in_source_fragments\n'))
    message(crayon::green(temp_num2, 'unlabeled features used to extract pairs','\n'))
    rm(temp_num1, temp_num2);gc()

    feature_sig_unlabel <- feature_sig_unlabel %>% dplyr::filter(id %in% result_peak_recognization_unlabel$base_peak)

    save(feature_sig_unlabel, file = file.path(path, '00_tracer_result', '00_intermediate_data', 'feature_sig_unlabel_merged.RData'))
  }



  # 3. Isotope pair finding
  cat('\n')
  message(crayon::blue('3. Find pairs from these enriched features ... \n'))

  isotope_label_matched <- pbapply::pbmapply(function(x, y, z){
    cat(x, '@', y, '\n')
    find_pair_feature(target_mz = x,
                      target_rt = y,
                      source_feature_id = z,
                      labeled_feature_table = feature_sig_label,
                      mz_tol = 15,
                      rt_tol = rt_tol)
  },
  x = feature_sig_unlabel$mz,
  y = feature_sig_unlabel$rt,
  z = feature_sig_unlabel$id,
  SIMPLIFY = FALSE) %>%
    dplyr::bind_rows()

  pair_table <- feature_sig_unlabel %>%
    dplyr::left_join(isotope_label_matched, by = c('id' = 'source_feature'), multiple = "all") %>%
    dplyr::select(id, matched_feature, mz:rt, actual_mz, actual_rt, mass_shift_label, p_values, q_values, fold_change, avg_abundance_case) %>%
    dplyr::rename('unlabeled_feature_id' = 'id',
                  'labelded_feature_id' = 'matched_feature',
                  'unlabeled_mz' = 'mz',
                  'unlabeled_rt' = 'rt',
                  'labeled_mz' = 'actual_mz',
                  'labeled_rt' = 'actual_rt',
                  'p_value_unlabel' = 'p_values',
                  'q_value_unlabel' = 'q_values',
                  'fold_change_unlabel' = 'fold_change',
                  'average_abundance_unlabel' = 'avg_abundance_case') %>%
    dplyr::filter(!is.na(labelded_feature_id))

  temp_idx <- match(pair_table$labelded_feature_id, raw_data_label$id)
  pair_table <- pair_table %>%
    dplyr::mutate(p_value_label = raw_data_label$p_values[temp_idx],
                  q_value_label = raw_data_label$q_values[temp_idx],
                  fold_change_label = raw_data_label$fold_change[temp_idx],
                  average_abundance_label = raw_data_label$avg_abundance_case[temp_idx])

  save(pair_table,
       file = file.path(path, '00_tracer_result', '00_intermediate_data', 'pair_table.RData'))


  # 4. Export the results
  message(crayon::yellow('Result: There are ', nrow(pair_table), 'labeled pairs are found.\n'))
  print(knitr::kable(pair_table))

  if (is_recognize_adducts) {
    table_export <- list('raw_data_unlabel' = raw_data_unlabel,
                         'raw_data_label' = raw_data_label,
                         'recognized_peaks_unlabel' = result_peak_recognization_unlabel,
                         'paired_table' = pair_table)
    dir.create(file.path(path, '00_tracer_result'), recursive = TRUE, showWarnings = FALSE)
    writexl::write_xlsx(table_export,
                        path = file.path(path, '00_tracer_result', 'tracer_pair_result.xlsx'),
                        format_headers = FALSE)
  } else {
    table_export <- list('raw_data_unlabel' = raw_data_unlabel,
                         'raw_data_label' = raw_data_label,
                         'paired_table' = pair_table)
    dir.create(file.path(path, '00_tracer_result'), recursive = TRUE, showWarnings = FALSE)
    writexl::write_xlsx(table_export,
                        path = file.path(path, '00_tracer_result', 'tracer_pair_result.xlsx'),
                        format_headers = FALSE)
  }


  cat('\n')
  if (nrow(pair_table) == 0) {
    message(crayon::red('No pairs are found! Skip the plots\n'))
    return(NULL)
  }

  message(crayon::blue('4. Plot these pairs ... \n'))
  message(crayon::blue('4.1 Extract EICs for pairs ...\n'))

  if (length(raw_data_folder) == 0) {
    raw_data_path_12C <- paste0(case_group, '_12C')
    raw_data_path_13C <- paste0(case_group, '_13C')
  } else {
    # need to be fixed
    raw_data_path_12C <- raw_data_folder[[paste0(case_group, '_12C')]]
    raw_data_path_13C <- raw_data_folder[[paste0(case_group, '_13C')]]
  }

  if (raw_data_path_12C %in% list.files(path)) {
    temp_files <- list.files(path = file.path(path, raw_data_path_12C), recursive = TRUE)
    eic_data_unlabel <- extract_eic_data(path = path,
                                         files = file.path(raw_data_path_12C, temp_files),
                                         mz_list = pair_table$unlabeled_mz,
                                         mz_tol = mz_tol)
  } else {
    eic_data_unlabel <- extract_eic_data(path = path,
                                         files_pattern = paste0(raw_data_path_12C, '.+\\.mzML|.+\\.mzXML'),
                                         mz_list = pair_table$unlabeled_mz,
                                         mz_tol = mz_tol)
  }

  if (raw_data_path_13C %in% list.files(path)) {
    temp_files <- list.files(path = file.path(path, raw_data_path_13C), recursive = TRUE)
    eic_data_label <- extract_eic_data(path = path,
                                       files = file.path(raw_data_path_13C, temp_files),
                                       mz_list = pair_table$labeled_mz,
                                       mz_tol = mz_tol)
  } else {
    eic_data_label <- extract_eic_data(path = path,
                                       files_pattern = paste0(raw_data_path_13C, '.+\\.mzML|.+\\.mzXML'),
                                       mz_list = pair_table$labeled_mz,
                                       mz_tol = mz_tol)
  }


  cat('\n')
  message(crayon::blue('4.2 Plot pairs\n'))

  # browser()
  temp_plot <- plot_isotope_pairs(signif_feature_unlabeled = feature_sig_unlabel,
                                  signif_feature_labeled = feature_sig_label,
                                  pair_table = pair_table,
                                  eic_data_unlabel = eic_data_unlabel,
                                  eic_data_label = eic_data_label,
                                  mz_tol = mz_tol,
                                  rt_unit = rt_unit)

  dir.create(file.path(path, '00_tracer_result'), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(plot = temp_plot,
                  filename = file.path(path, '00_tracer_result', 'isotope_pair_plot_overview.pdf'),
                  width = 15, height = 9)

  temp_plot_list <- plot_selected_pair(signif_feature_unlabeled = feature_sig_unlabel,
                                       signif_feature_labeled = feature_sig_label,
                                       pair_table = pair_table,
                                       eic_data_unlabel = eic_data_unlabel,
                                       eic_data_label = eic_data_label,
                                       mz_tol = mz_tol,
                                       rt_unit = rt_unit)

  ggplot2::ggsave(
    filename = file.path(path, '00_tracer_result', 'isotope_pair_list.pdf'),
    plot = gridExtra::marrangeGrob(temp_plot_list, nrow=1, ncol=1),
    width = 15, height = 9
  )

  message(crayon::blue('Done!\n'))


}



################################################################################
# find_pair_feature ------------------------------------------------------------

#' @title find_pair_feature
#' @author Zhiwei Zhou
#' @param target_mz target mz value
#' @param target_rt target retention time
#' @param labeled_feature_table labeled feature table
#' @param mz_tol mz tolerance in ppm
#' @param rt_tol retention time tolerance in minutes
#' @return A data frame containing the matched feature information
#' @importFrom dplyr %>%
#' @importFrom rcdk get.formula
#' @importFrom MassToolsMjhelf calcMF
#' @importFrom masstools mz_rt_match
#' @export
#' @examples

# target_mz = 105.0645
# target_rt = 7.629
# mz_tol <- 15
# rt_tol <- 0.05
# labeled_feature_table = feature_sig_label

find_pair_feature <- function(target_mz,
                              target_rt,
                              source_feature_id,
                              labeled_feature_table,
                              mz_tol,
                              rt_tol) {
  # browser()

  # predict target_mz formula
  pred_formula_table <- MassToolsMjhelf::calcMF(mz = target_mz, z = 1, ppm = mz_tol)

  # if the formula can't be predicted, export the NA
  if (length(pred_formula_table) < 1)  {
    output <- data.frame(source_feature = NA,
                         matched_feature = NA,
                         mass_shift_label = NA,
                         theo_mz = NA,
                         actual_mz = NA,
                         mz_error = NA,
                         theo_rt = NA,
                         actual_rt = NA,
                         rt_error = NA)
    return(output)
  }

  carbon_number_max <- pred_formula_table$MF %>%
    lapply(function(x){
      temp_result <- rcdk::get.formula(x)
      temp_result@isotopes %>%
        tibble::as_tibble() %>%
        dplyr::filter(isoto == 'C') %>%
        dplyr::pull(number) %>%
        as.numeric()
    }) %>%
    unlist() %>%
    unique() %>%
    max()

  possible_mz <- target_mz + c(1:carbon_number_max)*1.00335

  possible_table <- data.frame(mz = possible_mz,
                               rt = target_rt,
                               mass_shift_label = paste0('13C*', c(1:carbon_number_max)),
                               stringsAsFactors = FALSE)

  temp_label_table <- labeled_feature_table %>%
    dplyr::select(mz, rt) %>%
    as.data.frame()

  # masstools
  matched_result <- mz_rt_match(
    possible_table,
    temp_label_table,
    mz.tol = mz_tol,
    rt.tol = rt_tol,
    rt.error.type = 'abs'
  )

  if (is.null(matched_result)) {
    output <- data.frame(source_feature = NA,
                         matched_feature = NA,
                         mass_shift_label = NA,
                         theo_mz = NA,
                         actual_mz = NA,
                         mz_error = NA,
                         theo_rt = NA,
                         actual_rt = NA,
                         rt_error = NA)
    return(output)
  }

  # if (nrow(matched_result) == 0) {
  #   return(NA)
  # }

  # index
  idx1 <- matched_result$Index1
  idx2 <- matched_result$Index2

  # export matched pair table
  matched_feature <- labeled_feature_table[idx2,] %>%
    dplyr::pull(id)

  # source_feature <- paste0(target_mz, '@', target_rt)
  source_feature <- source_feature_id

  label_isotope <- possible_table[idx1,] %>% dplyr::pull(mass_shift_label)

  output <- matched_result %>%
    dplyr::mutate(source_feature = source_feature,
                  matched_feature = matched_feature,
                  mass_shift_label = label_isotope) %>%
    dplyr::select(source_feature, matched_feature, mass_shift_label, mz1, mz2, `mz error`, rt1, rt2, `rt error`) %>%
    dplyr::rename('theo_mz' = 'mz1',
                  'actual_mz' = 'mz2',
                  'mz_error' = 'mz error',
                  'theo_rt' = 'rt1',
                  'actual_rt' = 'rt2',
                  'rt_error' = 'rt error')

  return(output)
}








################################################################################
# startup massage --------------------------------------------------------------
.onAttach <- function(libname, pkgname){
  packageStartupMessage("
Version 0.1.3
-------------
Authors: Zhiwei Zhou
Maintainer: Zhiwei Zhou

Updates (20250626)
-------------
- Export parameters
- Update the README and NEWS
- Update pair_isotope_pairs_overview
- Update dependencies
- Add the function module of FunctionOthers to decrease necessary depends
- Check the raw data and ms2 data availabilty
- Support the mfg/mzXML ms2 data


")
}
