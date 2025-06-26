################################################################################
# modify_xcms_table ------------------------------------------------------------

#' @title modify_xcms_table
#' @author Zhiwei Zhou
#' @param table_xcms peak table from xcms
#' @param table_camera isotopes and adducts annotated by the CAMERA
#' @param table_sample_info sample information table
#' @param path working directory. Default: "."
#' @param file_output 'converted_peak_table.csv'
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' table_xcms <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/C12/00_raw_data_processing/Peak-table.csv'
#' table_camera <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/C12/00_raw_data_processing/adduct_result_camera.xlsx'
#' table_sample_info <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/sample_info.xlsx'
#' modify_xcms_table(table_xcms = table_xcms,
#'                   table_camera = table_camera,
#'                   table_sample_info = table_sample_info,
#'                   path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/C12')
#' }


# table_xcms <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/C12/00_raw_data_processing/Peak-table.csv'
# table_camera <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/C12/00_raw_data_processing/adduct_result_camera.xlsx'
# table_sample_info <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/sample_info.xlsx'

# modify_xcms_table(table_xcms = table_xcms,
#                   table_camera = table_camera,
#                   table_sample_info = table_sample_info,
#                   path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/C12')

modify_xcms_table <- function(table_xcms,
                              table_camera,
                              table_sample_info,
                              path = '.',
                              file_output = 'converted_peak_table.csv') {
  # browser()
  cat(crayon::bgBlue('Modifying xcms table...'), '\n')
  table_xcms <- readr::read_csv(table_xcms, show_col_types = FALSE)

  # According to file table to read the table_camera using read_xlsx or read_csv
  if (stringr::str_detect(table_camera, '\\.xlsx')) {
    # table_camera <- readxl::read_xlsx(file.path(path, table_camera))
    table_camera <- readxl::read_xlsx(table_camera)
  } else if (stringr::str_detect(table_camera, '\\.csv')) {
    # table_camera <- readr::read_csv(file.path(path, table_camera), show_col_types = FALSE)
    table_camera <- readr::read_csv(table_camera, show_col_types = FALSE)
  } else {
    stop('The table_camera should be a csv or xlsx file')
  }

  # According to file table to read the table using read_xlsx or read_csv
  if (stringr::str_detect(table_sample_info, '\\.xlsx')) {
    sample_info <- readxl::read_xlsx(table_sample_info)
  } else if (stringr::str_detect(table_sample_info, '\\.csv')) {
    sample_info <- readr::read_csv(table_sample_info, show_col_types = FALSE)
  } else {
    stop('The sample_info should be a csv or xlsx file')
  }

  temp_sample_id <- colnames(table_xcms)[colnames(table_xcms) %in% sample_info$sample_id]
  # select columns name, mzmed, rtmed, and columns after maxint
  table_xcms <- table_xcms %>%
    dplyr::select(name, mzmed, rtmed, any_of(temp_sample_id)) %>%
    dplyr::rename(id = name, mz = mzmed, rt = rtmed)

  # split isotopes columns into multiple columns according to the label "[\\d+]"
  isotope_table <- table_camera %>%
    dplyr::filter(!is.na(isotopes)) %>%
    dplyr::mutate(isotope_group = stringr::str_extract(isotopes, "\\[\\d+\\]"),
                  isotope_label = stringr::str_extract(isotopes, "\\[M.*].+"))

  unique_isotope_group <- unique(isotope_table$isotope_group)

  monoisotope_labels <- sapply(unique_isotope_group, function(x){
    isotope_data <- isotope_table %>%
      dplyr::filter(isotope_group == x) %>%
      dplyr::filter(isotope_label %in% c('[M]+', '[M]2+')) %>%
      dplyr::pull(name)
  })

  isotope_labels <- lapply(unique_isotope_group, function(x){
    isotope_data <- isotope_table %>%
      dplyr::filter(isotope_group == x) %>%
      dplyr::filter(!(isotope_label %in% c('[M]+', '[M]2+'))) %>%
      dplyr::pull(name)
  })

  names(isotope_labels) <- monoisotope_labels

  # remove the isotopes and add a column of isotope label
  table_xcms <- table_xcms %>%
    dplyr::filter(!(id %in% unlist(isotope_labels))) %>%
    dplyr::mutate('ms1_isotopes' = NA) %>%
    dplyr::select(id:rt, ms1_isotopes, everything())

  temp_idx <- match(monoisotope_labels, table_xcms$id)
  temp_label <- sapply(isotope_labels, function(x){paste0(x, collapse = '; ')})

  table_xcms$ms1_isotopes[temp_idx] <- temp_label

  if(!is.null(path)) {
    # write the modified table to a csv file
    readr::write_csv(table_xcms,
                     file.path(path, file_output))
  }


  cat(crayon::green('The xcms table is modified successfully!\n'),
      crayon::blue('The modified table is saved to:'),
      crayon::yellow(file.path(path, file_output)), '\n')
}




################################################################################
# modify_msdial_table ------------------------------------------------------------

#' @title modify_msdial_table
#' @author Zhiwei Zhou
#' @param table_msdial peak table from xcms
#' @param table_sample_info sample information table
#' @param path working directory. Default: "."
#' @param file_output 'converted_peak_table.csv'
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' table_msdial <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/hyuA_UA_48h_area.txt'
#' table_sample_info <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/sample_info.xlsx'
#' modify_msdial_table(table_msdial = table_msdial,
#'                      table_sample_info = table_sample_info,
#'                      path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial',
#'                      file_output = 'converted_peak_table.csv')
#' }


# table_msdial <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/hyuA_UA_48h_area.txt'
# table_sample_info <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/sample_info.xlsx'
#
# modify_msdial_table(table_msdial = table_msdial,
#                      table_sample_info = table_sample_info,
#                      path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial',
#                      file_output = 'peak_table_C12.csv')
#
# table_msdial <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/hyuA_13CUA_48h_area.txt'
# table_sample_info <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/sample_info.xlsx'
#
# modify_msdial_table(table_msdial = table_msdial,
#                     table_sample_info = table_sample_info,
#                     path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial',
#                     file_output = 'peak_table_C13.csv')


modify_msdial_table <- function(
    table_msdial,
    table_sample_info,
    path = '.',
    file_output = 'converted_peak_table.csv') {
  if (stringr::str_detect(table_msdial, '\\.xlsx')) {
    # raw_table <- readxl::read_xlsx(file.path(path, table_msdial))
    raw_table <- readxl::read_xlsx(table_msdial)
  } else {
    # raw_table <- readr::read_tsv(file.path(path, table_msdial))
    raw_table <- readr::read_tsv(table_msdial)
  }

  # According to file table to read the sample_info table using read_xlsx or read_csv
  if (stringr::str_detect(table_sample_info, '\\.xlsx')) {
    sample_info <- readxl::read_xlsx(table_sample_info)
  } else if (stringr::str_detect(table_sample_info, '\\.csv')) {
    sample_info <- readr::read_csv(table_sample_info, show_col_types = FALSE)
  } else {
    stop('The sample_info should be a csv or xlsx file')
  }

  temp_col_name <- raw_table %>% dplyr::slice(4) %>% as.character()
  temp_sample_id <- temp_col_name[temp_col_name %in% sample_info$sample_id]

  # if (!any(stringr::str_detect(temp_col_name, control_group))) {
  #   stop('No control group existed, please check it')
  # }
  # if (!any(stringr::str_detect(temp_col_name, case_group))) {
  #   stop('No case group existed, please check it')
  # }

  raw_table <- raw_table %>%
    dplyr::slice(-c(1:4))
  colnames(raw_table) <- temp_col_name

  raw_table <- raw_table %>%
    dplyr::select(`Average Mz`,
                  `Average Rt(min)`,
                  `MS1 isotopic spectrum`,
                  all_of(temp_sample_id)) %>%
    dplyr::rename('mz' = `Average Mz`,
                  'rt' = `Average Rt(min)`,
                  'ms1_isotopes' = `MS1 isotopic spectrum`) %>%
    dplyr::mutate(dplyr::across(-'ms1_isotopes', as.numeric)) %>%
    dplyr::mutate(mz = round(mz, 4),
                  rt = round(rt, 4)) %>%
    dplyr::mutate(id = paste0(mz, '@', rt)) %>%
    dplyr::select(id, dplyr::everything())


  if(!is.null(path)) {
    # write the modified table to a csv file
    readr::write_csv(raw_table,
                     file.path(path, file_output))
  }
}




################################################################################
# modify_mzmine_table ------------------------------------------------------------

#' @title modify_mzmine_table
#' @author Zhiwei Zhou
#' @param table_mzmine peak table from xcms
#' @param table_sample_info sample information table
#' @param path working directory. Default: "."
#' @param file_output 'converted_peak_table.csv'
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' table_mzmine <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_mzmine/mzmine_12C.csv'
#' table_sample_info <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_mzmine/sample_info.xlsx'
#' modify_mzmine_table(table_mzmine = table_mzmine,
#'                      table_sample_info = table_sample_info,
#'                      path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_mzmine',
#'                      file_output = 'converted_peak_table.csv')
#' }


# table_mzmine <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_mzmine/mzmine_12C.csv'
# table_sample_info <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_mzmine/sample_info.xlsx'
# path <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_mzmine'
# modify_mzmine_table(table_mzmine = table_mzmine,
#                      table_sample_info = table_sample_info,
#                      path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_mzmine',
#                      file_output = 'peak_table_C12.csv')
#
# table_mzmine <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_mzmine/mzmine_13C.csv'
# table_sample_info <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_mzmine/sample_info.xlsx'
# path <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_mzmine'
# modify_mzmine_table(table_mzmine = table_mzmine,
#                      table_sample_info = table_sample_info,
#                      path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_mzmine',
#                      file_output = 'peak_table_C13.csv')


modify_mzmine_table <- function(
    table_mzmine,
    table_sample_info,
    path = '.',
    file_output = 'converted_peak_table.csv') {
  if (stringr::str_detect(table_mzmine, '\\.xlsx')) {
    # raw_table <- readxl::read_xlsx(file.path(path, table_mzmine))
    raw_table <- readxl::read_xlsx(table_mzmine)
  } else if (stringr::str_detect(table_mzmine, '\\.csv')) {
    # raw_table <- readr::read_tsv(file.path(path, table_mzmine))
    raw_table <- readr::read_csv(table_mzmine, show_col_types = FALSE)
  } else {
    raw_table <- readr::read_tsv(table_mzmine, show_col_types = FALSE)
  }

  # According to file table to read the sample_info table using read_xlsx or read_csv
  if (stringr::str_detect(table_sample_info, '\\.xlsx')) {
    sample_info <- readxl::read_xlsx(table_sample_info)
  } else if (stringr::str_detect(table_sample_info, '\\.csv')) {
    sample_info <- readr::read_csv(table_sample_info, show_col_types = FALSE)
  } else {
    stop('The sample_info should be a csv or xlsx file')
  }

  temp_col_name <- raw_table %>% colnames()

  # select
  temp_idx <- stringr::str_detect(temp_col_name, paste(sample_info$sample_id, collapse = '|')) %>%
    which()
  temp_sample_id <- temp_col_name[temp_idx]


  raw_table <- raw_table %>%
    dplyr::select(`row m/z`,
                  `row retention time`,
                  `row identity (all IDs)`,
                  all_of(temp_sample_id)) %>%
    dplyr::rename('mz' = `row m/z`,
                  'rt' = `row retention time`,
                  'ms1_isotopes' = `row identity (all IDs)`) %>%
    dplyr::mutate(dplyr::across(-'ms1_isotopes', as.numeric)) %>%
    dplyr::mutate(mz = round(mz, 4),
                  rt = round(rt, 4)) %>%
    dplyr::mutate(id = paste0(mz, '@', rt)) %>%
    dplyr::select(id, dplyr::everything())

  colnames(raw_table)[-c(1:4)] <- stringr::str_replace(colnames(raw_table)[-c(1:4)], '\\.mzML Peak area', '')


  if(!is.null(path)) {
    # write the modified table to a csv file
    readr::write_csv(raw_table,
                     file.path(path, file_output))
  }
}




################################################################################
# check_parameters -------------------------------------------------------------

#' @title check_parameters
#' @author Zhiwei Zhou
#' @param peak_table_unlabel peak table for unlabeled samples
#' @param peak_table_label peak table for labeled samples
#' @param sample_info sample information table
#' @param path working directory. Default: "."
#' @param control_group control group name. Default: "WT"
#' @param case_group case group name. Default: "hyuA"
#' @param tracer_condition tracer condition, either '12C' or '13C'. Default: '12C'
#' @importFrom magrittr %>%
#' @importFrom crayon bgBlue red yellow green
#' @export
#' @examples
#' \dontrun{
#' check_parameters(peak_table_unlabel = 'peak_table_C12.csv',
#' peak_table_label = 'peak_table_C13.csv',
#' sample_info = 'sample_info.xlsx',
#' path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/',
#' control_group = c("WT"),
#' case_group = c('hyuA'))
#' }


# check_parameters(peak_table_unlabel = 'peak_table_C12.csv',
#                   peak_table_label = 'peak_table_C13.csv',
#                   sample_info = 'sample_info.xlsx',
#                   path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/',
#                   control_group = c("WT"),
#                   case_group = c('hyuA'))

# peak_table_unlabel <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/peak_table_C12.csv'
# peak_table_label <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/peak_table_C13.csv'
# sample_info <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/hyuA/sample_info.xlsx'
# path <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/'

# peak_table_unlabel <- 'peak_table_C12.csv'
# peak_table_label <- 'peak_table_C13.csv'
# sample_info <- 'sample_info.xlsx'
# path <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/'
# tracer_condition <- '12C'

check_parameters <- function(peak_table_unlabel,
                             peak_table_label,
                             sample_info,
                             path = '.',
                             raw_data_folder = NULL,
                             control_group = c("WT"),
                             case_group = c('hyuA')){
  cat(crayon::bgBlue('Checking parameters...'), '\n')

  # check format of peak table
  if (missing(peak_table_unlabel)) {
    stop('Please provide the peak_table\n')
  }

  if (missing(peak_table_label)) {
    stop('Please provide the peak_table\n')
  }

  if (missing(sample_info)) {
    stop('Please provide the sample_info\n')
  }

  # check availability
  file_list <- list.files(path)
  if (!(peak_table_unlabel %in% file_list)) {
    stop(crayon::red('The peak_table_unlabel is not available in the path'))
  }
  if (!(peak_table_label %in% file_list)) {
    stop(crayon::red('The peak_table_label is not available in the path'))
  }
  if (!(sample_info %in% file_list)) {
    stop(crayon::red('The sample_info is not available in the path'))
  }

  # check availability of raw data
  if (length(raw_data_folder) == 0) {
    raw_data_path_12C <- paste0(case_group, '_12C')
    raw_data_path_13C <- paste0(case_group, '_13C')
  } else {
    raw_data_path_12C <- raw_data_folder[[paste0(case_group, '_12C')]]
    raw_data_path_13C <- raw_data_folder[[paste0(case_group, '_13C')]]
  }

  if (!(raw_data_path_12C %in% file_list) | !(raw_data_path_13C %in% file_list)) {
    stop(crayon::red('The raw data folder is not available in the path\n',
                     paste(c(raw_data_path_12C, raw_data_path_13C), collapse = ', ')))
  }
  # check raw data format belong to mzML/mzXML
  raw_data_files_12C <- list.files(file.path(path, raw_data_path_12C))
  raw_data_files_13C <- list.files(file.path(path, raw_data_path_13C))
  if (length(raw_data_files_12C) == 0 | length(raw_data_files_13C) == 0) {
    stop(crayon::red('The raw data folder is empty\n',
                     paste(c(raw_data_path_12C, raw_data_path_13C), collapse = ', ')))
  }
  if (!all(stringr::str_detect(raw_data_files_12C, '\\.mzML$|\\.mzXML$')) |
      !all(stringr::str_detect(raw_data_files_13C, '\\.mzML$|\\.mzXML$'))) {
    stop(crayon::red('The raw data files should be in mzML/mzXML/mgf format\n',
                     paste(c(raw_data_files_12C[!stringr::str_detect(raw_data_files_12C, '\\.mzML$|\\.mzXML$')],
                             raw_data_files_13C[!stringr::str_detect(raw_data_files_13C, '\\.mzML$|\\.mzXML$')]),
                           collapse = ', ')))
  }


  # check availability of ms2 data
  if (!("ms2" %in% file_list)) {
    stop(crayon::red('The ms2 data folder is not available in the path\n',
                     'Please provide the ms2 data folder with name "ms2"'))
  }

  # check ms2 format belong to mzML/mzXML/mgf
  ms2_files <- list.files(file.path(path, 'ms2'))
  if (length(ms2_files) == 0) {
    stop(crayon::red('The ms2 data folder is empty\n',
                     'Please provide the ms2 data folder with name "ms2"'))
  }
  if (!all(stringr::str_detect(ms2_files, '\\.mzML$|\\.mzXML$|\\.mgf$'))) {
    stop(crayon::red('The ms2 data files should be in mzML/mzXML/mgf format\n',
                     paste(ms2_files[!stringr::str_detect(ms2_files, '\\.mzML$|\\.mzXML$|\\.mgf$')], collapse = ', ')))
  }


  # read the tables
  if (stringr::str_detect(peak_table_unlabel, '\\.xlsx')) {
    peak_table_unlabel <- readxl::read_xlsx(file.path(path, peak_table_unlabel))
  } else if (stringr::str_detect(peak_table_unlabel, '\\.csv')) {
    peak_table_unlabel <- readr::read_csv(file.path(path, peak_table_unlabel), show_col_types = FALSE)
  } else {
    stop('The peak_table_unlabel should be a csv or xlsx file')
  }

  if (stringr::str_detect(peak_table_label, '\\.xlsx')) {
    peak_table_label <- readxl::read_xlsx(file.path(path, peak_table_label))
  } else if (stringr::str_detect(peak_table_label, '\\.csv')) {
    peak_table_label <- readr::read_csv(file.path(path, peak_table_label), show_col_types = FALSE)
  } else {
    stop('The peak_table_label should be a csv or xlsx file')
  }

  if (stringr::str_detect(sample_info, '\\.xlsx')) {
    sample_info <- readxl::read_xlsx(file.path(path, sample_info))
  } else if (stringr::str_detect(sample_info, '\\.csv')) {
    sample_info <- readr::read_csv(file.path(path, sample_info), show_col_types = FALSE)
  } else {
    stop('The sample_info should be a csv or xlsx file')
  }

  # check format: the first 4 columns of peak_table should be id, mz, rt, and ms1_isotopes
  if (!all(c('id', 'mz', 'rt', 'ms1_isotopes') %in% colnames(peak_table_unlabel)[1:4])) {
    stop(crayon::red('The peak_table_unlabel should contain id, mz, rt, and ms1_isotopes columns'))
  }
  if (!all(c('id', 'mz', 'rt', 'ms1_isotopes') %in% colnames(peak_table_label)[1:4])) {
    stop(crayon::red('The peak_table_label should contain id, mz, rt, and ms1_isotopes columns'))
  }
  # check format: the first 3 columns of sample_info should be sample_id, group, and tracer_group
  if (!all(c('sample_id', 'group', 'tracer_group') %in% colnames(sample_info)[1:3])) {
    stop(crayon::red('The sample_info should contain sample_id, group, and tracer_group columns'))
  }


  # check control_group and case_group whether included in the sample info
  if (!all(c(control_group, case_group) %in% unique(sample_info$group))) {
    stop(crayon::red('The control_group and case_group do not included in the sample info file\n',
                     paste(c(control_group, case_group), collapse = ', ')))
  }

  # check sample_id in the sample info is consistent with the peak table
  sample_id_unlabel <- sample_info %>%
    dplyr::filter(type == 'sample') %>%
    dplyr::filter(tracer_group == '12C') %>%
    dplyr::pull(sample_id)

  sample_id_label <- sample_info %>%
    dplyr::filter(type == 'sample') %>%
    dplyr::filter(tracer_group == '13C') %>%
    dplyr::pull(sample_id)

  if (!all(sample_id_unlabel %in% colnames(peak_table_unlabel))) {
    temp <- which(!(sample_id_unlabel %in% colnames(peak_table_unlabel))) %>%
      sample_id_unlabel[.]

    stop(length(temp), ' samples are not included in the peak_table_unlabel\n',
         paste(temp, collapse = ', '))
  }

  if (!all(sample_id_label %in% colnames(peak_table_label))) {
    temp <- which(!(sample_id_label %in% colnames(peak_table_label))) %>%
      sample_id_label[.]

    stop(length(temp), ' samples are not included in the peak_table_label\n',
         paste(temp, collapse = ', '))
  }

  # print checked results as a table:
    # peak_table_unlabel and peak_table_label: number of features (rows), number of samples (sample names)
    # sample_info: number of samples (sample names), number of groups, number of tracer groups
    # raw_data_folder: number of raw data files in each folder, file format
    # ms2: number of ms2 files, file format

  # print a split line '------------"

  cat(crayon::silver('----------------------------\n'))
  cat(crayon::green('The parameters are checked successfully!\n'))
  cat(crayon::blue('Peak table (unlabeled):'),
      crayon::yellow(length(colnames(peak_table_unlabel)) - 4), 'samples with',
      crayon::yellow(nrow(peak_table_unlabel)), 'features\n')
  cat(crayon::blue('Peak table (unlabeled):'),
      crayon::yellow(length(colnames(peak_table_label)) - 4), 'samples with',
      crayon::yellow(nrow(peak_table_label)), 'features\n')
  cat(crayon::blue('Sample info:\n'),
      crayon::yellow(length(unique(sample_info$sample_id))), 'samples \n',
      crayon::yellow(length(unique(sample_info$group))), 'groups:',
      paste(unique(sample_info$group), collapse = ', '), '\n',
      crayon::yellow(length(unique(sample_info$tracer_group))), 'tracer groups:',
      paste(unique(sample_info$tracer_group), collapse = ', '), '\n')
  cat(crayon::blue('Raw data folders:\n'),
      crayon::yellow(raw_data_path_12C), ':',
      length(list.files(file.path(path, raw_data_path_12C))), 'files with format:',
      paste(unique(stringr::str_extract(raw_data_files_12C, '\\.mzML$|\\.mzXML$')), collapse = ', '), '\n',
      crayon::yellow(raw_data_path_13C), ':',
      length(list.files(file.path(path, raw_data_path_13C))), 'files with format:',
      paste(unique(stringr::str_extract(raw_data_files_13C, '\\.mzML$|\\.mzXML$')), collapse = ', '), '\n')
  cat(crayon::blue('MS2 data folder:\n'),
      crayon::yellow('ms2'), ':',
      length(ms2_files), 'files with format:',
      paste(unique(stringr::str_extract(ms2_files, '\\.mzML$|\\.mzXML$|\\.mgf$')), collapse = ', '), '\n')
  cat(crayon::silver('----------------------------\n'))



}
