# table_xcms <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/C13/00_raw_data_processing/Peak-table.csv'
# table_camera <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/C13/00_raw_data_processing/adduct_result_camera.xlsx'
# table_sample_info <- '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/sample_info.xlsx'
# modify_xcms_table(table_xcms = table_xcms,
#                   table_camera = table_camera,
#                   table_sample_info = table_sample_info,
#                   path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data/hyuA/C13')
#
#
# modify_xcms_table(table_xcms = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_xcms/00_raw_data_processing_12C/00_raw_data_processing/Peak-table.csv',
#                   table_camera = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_xcms/00_raw_data_processing_12C/00_raw_data_processing/adduct_result_camera.xlsx',
#                   table_sample_info = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_xcms/sample_info.xlsx',
#                   path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_xcms/00_raw_data_processing_12C',
#                   file_output = 'peak_table_C12.csv')
#
# modify_xcms_table(table_xcms = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_xcms/00_raw_data_processing_13C/00_raw_data_processing/Peak-table.csv',
#                   table_camera = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_xcms/00_raw_data_processing_13C/00_raw_data_processing/adduct_result_camera.xlsx',
#                   table_sample_info = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_xcms/sample_info.xlsx',
#                   path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_xcms/00_raw_data_processing_13C',
#                   file_output = 'peak_table_C13.csv')

# find_intemidates(peak_table_unlabel = 'peak_table_C12.csv',
#                  peak_table_label = 'peak_table_C13.csv',
#                  sample_info = 'sample_info.xlsx',
#                  path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_msdial/',
#                  polarity = c('positive', 'negative'),
#                  control_group = c("WT"),
#                  case_group = c('hyuA'),
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  p_adjust = TRUE,
#                  fold_change_cutoff = 20,
#                  is_recognize_adducts = TRUE)

#
# find_intemidates(peak_table_unlabel = 'peak_table_C12.csv',
#                  peak_table_label = 'peak_table_C13.csv',
#                  sample_info = 'sample_info.xlsx',
#                  path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_xcms/',
#                  polarity = 'positive',
#                  control_group = c("WT"),
#                  case_group = c('hyuA'),
#                  raw_data_folder = list(
#                    'hyuA_12C' = 'hyuA_12C',
#                    'hyuA_13C' = 'hyuA_13C'
#                  ),
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  p_adjust = TRUE,
#                  fold_change_cutoff = 20,
#                  is_recognize_adducts = TRUE)

# load('~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_xcms/00_tracer_result/00_intermediate_data/feature_sig_unlabel')
# load('~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_xcms/00_tracer_result/00_intermediate_data/feature_sig_label')
# load('~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_xcms/00_tracer_result/00_intermediate_data/feature_sig_unlabel_merged.RData')
#
# find_intemidates(peak_table_unlabel = 'peak_table_C12.csv',
#                  peak_table_label = 'peak_table_C13.csv',
#                  sample_info = 'sample_info.xlsx',
#                  path = '~/Project/00_Uric_Acid_project/Data/20250606_isopairfind_test/Demo_data_mzmine/',
#                  polarity = 'positive',
#                  control_group = c("WT"),
#                  case_group = c('hyuA'),
#                  raw_data_folder = list(
#                    'hyuA_12C' = 'hyuA_12C',
#                    'hyuA_13C' = 'hyuA_13C'
#                  ),
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  p_adjust = TRUE,
#                  fold_change_cutoff = 20,
#                  is_recognize_adducts = TRUE)


# find_intemidates(peak_table_unlabel = 'peak_table_C12.csv',
#                  peak_table_label = 'peak_table_C13.csv',
#                  sample_info = 'sample_info.csv',
#                  path = '~/Project/00_Uric_Acid_project/Data/20250619_demo_data/test/',
#                  polarity = 'positive',
#                  control_group = c("WT"),
#                  case_group = c('hyuA'),
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  p_adjust = TRUE,
#                  fold_change_cutoff = 20,
#                  is_recognize_adducts = TRUE)
