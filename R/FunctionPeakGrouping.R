#' @title recognize_rela_peak
#' @description recognize related peaks for enriched peaks, e.g. adducts, neutral losses, in-source fragments
#' @author Zhiwei Zhou
#' @param peak_table
#' @param peak_table_signif
#' @param path_dir
#' @param polarity "positive", "negative"
#' @param raw_data_folder
#' @param tol_mz
#' @param tol_rt
#' @param cutoff_ssc
#' @param cutoff_ppc

# raw_data <- readxl::read_xlsx('~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/hyuA/00_tracer_result/tracer_pair_result.xlsx')
# peak_table <- raw_data %>%
#   rename('name' = id) %>%
#   select(name, mz, rt, WT_UA_1:hyuA_UA_3)
#
# peak_table_signif <- raw_data %>%
#   dplyr::filter(increase_label == 'signif_increase') %>%
#   rename('name' = id) %>%
#   select(name, mz, rt, WT_UA_1:hyuA_UA_3)
#
# polarity <- 'positive'
# path_dir <- '~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/hyuA/'
# tol_mz <- 20
# tol_rt <- 0.05
# cutoff_ssc <- 0.8
# cutoff_ppc <- -1

# test <- recognize_rela_peak(peak_table = peak_table,
#                             peak_table_signif = peak_table_signif,
#                             path_dir = '~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/hyuA/',
#                             polarity = 'positive',
#                             raw_data_folder = 'hyuA_UA',
#                             tol_mz = 20,
#                             tol_rt = 0.05,
#                             cutoff_ssc = 0.8,
#                             cutoff_ppc = 0.8)


recognize_rela_peak <- function(peak_table,
                                peak_table_signif,
                                path_dir = '.',
                                polarity = c('positive', 'negative'),
                                raw_data_folder = '.',
                                tol_mz = 20,
                                tol_rt = 0.05,
                                cutoff_ssc = 0.8,
                                cutoff_ppc = 0.8){

  # browser()
  # generate peak group
  temp_adduct <- ifelse(polarity == 'positive', '[M+H]+', '[M-H]-')
  data('lib_adduct_nl', envir = environment())

  cat('Generate peak group ...\n')
  raw_data_files <- list.files(file.path(path_dir, raw_data_folder), recursive = TRUE)
  list_peak_group <- lapply(seq_along(peak_table_signif$name), function(i){
    # cat(i, ' ')
    temp_name <- peak_table_signif$name[i]

    result <- generatePeakGroup(peak_table = peak_table,
                                raw_data_files = raw_data_files,
                                base_peak_name = temp_name,
                                base_peak_adduct = temp_adduct,
                                tol_mz = tol_mz,
                                tol_rt = tol_rt,
                                path = file.path(path_dir, raw_data_folder),
                                cutoff_ssc = cutoff_ssc,
                                cutoff_ppc = -1) # not apply ssc in peak group generation
    return(result)
  })

  names(list_peak_group) <- peak_table_signif$name

  dir.create(file.path(path_dir, '00_tracer_result', "00_intermediate_data"),
             showWarnings = FALSE,
             recursive = TRUE)

  save(list_peak_group,
       file = file.path(path_dir,
                        '00_tracer_result',
                        "00_intermediate_data",
                        'list_peak_group.RData'),
       version = 2)

  cat('\n'); cat('------------------------------------------------------\n')
  cat('Retrieve ms2 ...\n')
  raw_ms2 <- retrieve_ms2_data(peak_table = peak_table,
                               path = path_dir,
                               ms2_type = 'mzML',
                               is_include_precursor = TRUE,
                               int_ms2_min_abs = 50,
                               int_ms2_min_relative = 0.01,
                               mz_tol_combine_ms1_ms2 = 20,
                               rt_tol_combine_ms1_ms2 = 20)

  save(raw_ms2,
       file = file.path(path_dir,
                        '00_tracer_result',
                        "00_intermediate_data",
                        'raw_ms2.RData'),
       version = 2)

  cat('\n'); cat('------------------------------------------------------\n')
  cat('Annotate adduct, neutral loss, in-source fragments ...\n')
  ms2_data <- raw_ms2$spec
  list_peak_group_annotation <- pbapply::pblapply(seq_along(list_peak_group), function(i){
    annotatePeakGroup(peak_group = list_peak_group[[i]],
                      ms2_data = ms2_data,
                      polarity = polarity,
                      tol_mz = tol_mz,
                      is_ms2_check = FALSE,
                      cutoff_ssc = -1)
  })

  names(list_peak_group_annotation) <- names(list_peak_group)

  dir.create(file.path(path_dir,
                       '00_tracer_result',
                       "00_intermediate_data"),
             showWarnings = FALSE,
             recursive = TRUE)

  save(list_peak_group_annotation,
       file = file.path(path_dir,
                        '00_tracer_result',
                        "00_intermediate_data",
                        'list_peak_group_annotation.RData'),
       version = 2, compress = 'gzip')


  cat('\n'); cat('------------------------------------------------------\n')
  cat('Merge adducts, neutral loss, ISF ...\n')
  list_peak_group_annotation_merge <- merge_peak_groups(list_peak_group_annotation = list_peak_group_annotation,
                                                        path_dir = path_dir)

  save(list_peak_group_annotation_merge,
       file = file.path(path_dir,
                        '00_tracer_result',
                        "00_intermediate_data",
                        'list_peak_group_annotation_merge.RData'),
       version = 2, compress = 'gzip')

  cat('\n'); cat('------------------------------------------------------\n')
  cat('Annotate ISF (high-correlation) across peak groups ...\n')

  list_peak_group_annotation_merge <- self_check_isf(list_peak_group_annotation_merge = list_peak_group_annotation_merge,
                                                     ssc_cutoff = cutoff_ssc,
                                                     ppc_cutoff = cutoff_ppc)

  save(list_peak_group_annotation_merge,
       file = file.path(path_dir,
                        '00_tracer_result',
                        "00_intermediate_data",
                        'list_peak_group_annotation_merge.RData'),
       version = 2, compress = 'gzip')

  return(list_peak_group_annotation_merge)

}





################################################################################
# annotatePeakGroup ------------------------------------------------------------

#' @title annotatePeakGroup
#' @author Zhiwei Zhou
#' @description annotate peak group with isotopes, adducts, neutral loss and in-source fragments
#' @param peak_group a object of peak_group
#' @param polarity
#' @param tol_mz Default: 10 ppm
#' @param is_ms2_check whether compare ms2 similarity of neutral loss to base peak ms2. Default: TRUE
#' @param ms2_score_cutoff Default: -1; # -1 represent not filter
#' @param cutoff_ssc Whether apply ssc filter in annotation. Default: 0.3
#' @param cutoff_ssc_int If the peak intensity less than this value, cutoff_ssc_int is invalid.
#' @param is_rule_limitation Whether apply rule limitation in adductAnnotation and nlAnnotation (e.g. '[2M+H]+' needs '[M+H]+'). Default: TRUE
#' @param cutoff_topN Top N fragments in base peak ms2 used for ISF annotation. Default: 5
#' @export
#' @example
#' load(system.file("tempdata", "list_peak_group_200805.RData", package="MetDNA2"))
#' load(system.file("tempdata", "raw_msms_200805.RData", package="MetDNA2"))
#' peak_group_L0076 <- list_peak_group$`M482T929_[M-H]-`
#' test_L0076 <- annotatePeakGroup(peak_group = peak_group_L0076,
#'                                 ms2_data = raw_msms,
#'                                 polarity = 'negative',
#'                                 tol_mz = 10,
#'                                 is_ms2_check = TRUE,
#'                                 ms2_score_cutoff = -1)





# peak_group <- list_peak_group$`173.0669@3.353`
# ms2_data <- raw_ms2$spec
# polarity <- 'positive'
# tol_mz <- 25
# is_ms2_check <- TRUE
# ms2_score_cutoff <- -1
# is_rule_limitation <- TRUE
# cutoff_topN <- 5
# cutoff_ssc = -1
# cutoff_ssc_int = -1
# load('./data/lib_adduct_nl.rda')

setGeneric(name = 'annotatePeakGroup',
           def = function(
    peak_group,
    ms2_data,
    polarity = c('positive', 'negative'),
    tol_mz = tol_mz,
    is_ms2_check = TRUE, # compare ms2 similarity of neutral loss
    ms2_score_cutoff = -1, # -1 represent not filter
    lib_adduct_nl = lib_adduct_nl,
    cutoff_ssc = 0.3,
    cutoff_ssc_int = 3000,
    is_rule_limitation = TRUE,
    cutoff_topN = 5,
    ...
           ){
             # browser()
             # match.arg(polarity)

             if (missing(peak_group)) {
               stop('Please input peak group!\n')
             }

             if (missing(ms2_data)) {
               stop('Please input ms2 data!\n')
             }

             # browser()

             peak_list <- peak_group@peak_list %>% initializePeakList()
             query_mz <- peak_group@base_peak_mz
             query_rt <- peak_group@base_peak_rt
             # query_ccs <- peak_group@base_peak_ccs
             query_adduct <- peak_group@base_peak_adduct
             query_peak_name <- peak_group@base_peak_name
             peak_list_ms2 <- match(peak_list$peak_name, names(ms2_data)) %>%
               ms2_data[.] %>%
               .[!is.na(names(.))]


             # result_isotope <- annotateIsotope(peak_list = peak_list,
             #                                   query_mz = query_mz,
             #                                   query_rt = query_rt,
             #                                   # query_ccs = query_ccs,
             #                                   query_peak_name = query_peak_name,
             #                                   tol_mz = tol_mz,
             #                                   cutoff_ssc = NULL)
             # cat('iso\t')
             # result_isotope <- annotateIsotope(peak_list = peak_list,
             #                                   query_mz = query_mz,
             #                                   query_rt = query_rt,
             #                                   query_ccs = query_ccs,
             #                                   query_peak_name = query_peak_name,
             #                                   tol_mz = tol_mz,
             #                                   cutoff_ssc = NULL,
             #                                   ...)

             # cat('addu\t')
             result_adduct <- annotateAdduct(peak_list = peak_list,
                                             query_mz = query_mz,
                                             query_rt = query_rt,
                                             query_peak_name = query_peak_name,
                                             query_adduct = query_adduct,
                                             polarity = polarity,
                                             tol_mz = tol_mz,
                                             cutoff_ssc = cutoff_ssc,
                                             cutoff_ssc_int = cutoff_ssc_int,
                                             is_rule_limitation = is_rule_limitation)
             # cat('nl\t')
             result_nl <- annotateNeutralLoss(peak_list = peak_list,
                                              peak_list_ms2 = peak_list_ms2,
                                              query_mz = query_mz,
                                              query_rt = query_rt,
                                              query_peak_name = query_peak_name,
                                              query_adduct = query_adduct,
                                              polarity = polarity,
                                              tol_mz = tol_mz,
                                              is_ms2_check = is_ms2_check,
                                              ms2_score_cutoff = ms2_score_cutoff,
                                              cutoff_ssc = cutoff_ssc,
                                              cutoff_ssc_int = cutoff_ssc_int,
                                              is_rule_limitation = is_rule_limitation)


             # cat('isf\t')
             result_isf <- annotateISF(peak_list = peak_list,
                                       peak_list_ms2 = peak_list_ms2,
                                       query_peak_name = query_peak_name,
                                       query_mz = query_mz,
                                       query_rt = query_rt,
                                       tol_mz = tol_mz,
                                       cutoff_ssc = NULL,
                                       cutoff_topN = cutoff_topN)

             peak_list_annotation <- writePeakListAnnotation(pl_annotation_new = result_adduct) %>%
               writePeakListAnnotation(pl_annotation_new = result_nl) %>%
               writePeakListAnnotation(pl_annotation_new = result_isf)

             # # annotate isotopes for adducts, neutral loss, and in-source fragments
             # # cat('PURR\t')
             # temp_new <- extractIsotopeNew(peak_list_annotation)
             #
             # if (length(temp_new) > 0) {
             #   # result_isotope_new <- purrr::map(seq_along(temp_new),
             #   result_isotope_new <- lapply(seq_along(temp_new),
             #                                function(i){
             #                                  temp <- temp_new[i]
             #                                  idx <- which(peak_list$peak_name == temp)
             #                                  # temp_data <- peak_list %>% dplyr::filter(peak_name == temp_new[i])
             #
             #                                  temp_data <- peak_list[idx,]
             #
             #                                  annotateIsotope(peak_list = peak_list,
             #                                                  query_mz = temp_data$mz,
             #                                                  query_rt = temp_data$rt,
             #                                                  query_ccs = temp_data$ccs,
             #                                                  query_peak_name = temp_new[i],
             #                                                  tol_mz = tol_mz,
             #                                                  cutoff_ssc = NULL)
             #
             #                                }) %>% dplyr::bind_rows()
             #
             #   # cat('PURR2\t')
             #   peak_list_annotation <- writePeakListAnnotation(pl_annotation_old = peak_list_annotation,
             #                                                   pl_annotation_new = result_isotope_new)
             #
             # }

             # if the base_peak don't meet the rule,
             #    remove all annotations of this peak group
             is_rule_valid <- query_adduct %in% peak_list_annotation$annotation

             if (!is_rule_valid){
               peak_list_merged <- tibble::tibble(peak_name = character(),
                                                  mz = numeric(),
                                                  rt = numeric(),
                                                  # ccs = numeric(),
                                                  ssc = numeric(),
                                                  ppc = numeric(),
                                                  int = numeric(),
                                                  mz_error = numeric(),
                                                  rt_error = numeric(),
                                                  ccs_error = numeric(),
                                                  int_ratio = numeric(),
                                                  ms2_score = numeric(),
                                                  query_peak = character(),
                                                  isotopeAnnotation = character(),
                                                  adductAnnotation = character(),
                                                  neutralLossAnnotation = character(),
                                                  isfAnnotation = character())

               peak_group@peak_list_annotated <- peak_list_merged

               return(peak_group)

             }

             # cat('merge\t')
             # remove redunancy of annotation
             peak_list_merged <- mergePeakListAnnotation(pl_annotation = peak_list_annotation)

             peak_group@peak_list_annotated <- peak_list_merged
             # cat('\n')
             return(peak_group)
           })

# mergePeakListAnnotation ------------------------------------------------------
# mergePeakListAnnotation(pl_annotation = peak_list_annotation) %>% View()
setGeneric(name = 'mergePeakListAnnotation',
           def = function(
    pl_annotation
           ){

             # browser()
             # when these pair existed: M+_[M+H]+ & M-_[M-H]- & their isotope annotation overlaped
             #   rm isotope annotations from M+ abd M-

             if (all(c('M+', '[M+H]+') %in% pl_annotation$annotation) |
                 all(c('M-', '[M-H]-') %in% pl_annotation$annotation)) {

               temp_peak1 <- which(pl_annotation$annotation == '[M+H]+' | pl_annotation$annotation == '[M-H]-') %>% pl_annotation$peak_name[.]
               temp_peak2 <- which(pl_annotation$annotation == 'M+' | pl_annotation$annotation == 'M-') %>% pl_annotation$peak_name[.]

               temp_peak1_list <- pl_annotation %>% dplyr::filter(query_peak == temp_peak1 & type == 'isotopeAnnotation')
               temp_peak2_list <- pl_annotation %>% dplyr::filter(query_peak == temp_peak2 & type == 'isotopeAnnotation')

               is_overlap <- sum(temp_peak2_list$peak_name %in% temp_peak1_list$peak_name)

               if (is_overlap > 0) {
                 idx_rm <- which((pl_annotation$query_peak == temp_peak2) & (pl_annotation$annotation %in% c('[M+1]', '[M+2]', '[M+3]')))
                 pl_annotation <- pl_annotation[-idx_rm,]
               }

             }

             isotope_list <- pl_annotation %>% dplyr::filter(type == 'isotopeAnnotation')

             # if one peak were annotated as [M] & [M+1]/[M+2], remove reduancy
             #    [M+2] > [M+1] > [M]

             isotope_list <- isotope_list %>%
               dplyr::arrange(peak_name, desc(annotation)) %>%
               dplyr::distinct(peak_name, .keep_all = TRUE)


             other_list <- pl_annotation %>% dplyr::filter(type != 'isotopeAnnotation')

             # if the peak was annotated as [M+1], [M+2], [M+3],
             #    then remove other annotations of this peak
             list_confilict_annotation <- isotope_list %>%
               dplyr::filter(annotation != '[M]') %>%
               dplyr::pull(peak_name)

             # Other conflict was deplicated according to m/z error
             other_list <- other_list %>%
               dplyr::filter(!(peak_name %in% list_confilict_annotation)) %>%
               dplyr::arrange(abs(mz_error)) %>%
               dplyr::distinct(peak_name, .keep_all = TRUE)

             # convert to wider table
             result_list <- isotope_list %>%
               tidyr::pivot_wider(names_from = type, values_from = annotation)

             other_list <- other_list %>%
               tidyr::pivot_wider(names_from = type, values_from = annotation)

             temp_list <- other_list %>%
               dplyr::select(peak_name, dplyr::contains('Annotation'))

             if (nrow(result_list) < 1) {
               result_list <- other_list
             } else {
               result_list <- result_list %>%
                 dplyr::left_join(temp_list, by = 'peak_name') %>%
                 dplyr::arrange(mz)
             }



             # replace mz_error, rt_error int_ratio, ms2_score & query_peak with correct annotation
             idx <- match(temp_list$peak_name, result_list$peak_name)
             result_list$mz_error[idx] <- other_list$mz_error
             result_list$rt_error[idx] <- other_list$rt_error
             # result_list$ccs_error[idx] <- other_list$ccs_error
             result_list$int_ratio[idx] <- other_list$int_ratio
             result_list$ms2_score[idx] <- other_list$ms2_score
             result_list$query_peak[idx] <- other_list$query_peak


             # impute format
             idx <- c('isotopeAnnotation', 'adductAnnotation', 'neutralLossAnnotation', 'isfAnnotation') %>%
               match(., colnames(result_list)) %>%
               is.na() %>%
               which()

             if (length(idx)>0) {
               imputate_list <- matrix(nrow = nrow(result_list), ncol = length(idx))
               colnames(imputate_list) <- c('isotopeAnnotation', 'adductAnnotation', 'neutralLossAnnotation', 'isfAnnotation')[idx]

               options(readr.num_columns = 0)
               imputate_list <- imputate_list %>% tibble::as_tibble()

               result_list <- result_list %>% dplyr::bind_cols(imputate_list)
             }

             result_list <- result_list %>%
               dplyr::select(peak_name:query_peak, isotopeAnnotation, adductAnnotation,
                             neutralLossAnnotation, isfAnnotation)

             return(result_list)
           })



################################################################################
# generate peak group ----------------------------------------------------------
setClass(Class = "PeakGroup",
         slots = list(base_peak_name = "character",
                      base_peak_mz = "numeric",
                      base_peak_rt = "numeric",
                      # base_peak_ccs = 'numeric',
                      base_peak_int = 'numeric',
                      base_peak_adduct = 'character',
                      group_size = 'numeric',
                      peak_list = 'data.frame',
                      peak_list_annotated = 'data.frame')
)


#' @title generatePeakGroup
#' @author Zhiwei Zhou
#' @param peak_table peak table (row--peaks; column--samples). First 3 column should be named as "name", "mz", "rt"
#' @param base_peak_name name of base peak
#' @param base_peak_adduct adduct format of base peak
#' @param base_peak_seed seed of base peak
#' @param tol_rt rt tolerance for peak grouping; Default: 3 (second)
#' @param cutoff_ssc sample-sample correlation; Default: 0
#' @export
#' @examples

# load('~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/hyuA/raw_data_unlabel_231003.RData')

# peak_table <- readr::read_csv('./inst/extdata/peak_table_200STD_neg_200805.csv')
# peak_table_annotated <- readr::read_csv('./inst/extdata/peak_table_annotated_200STD_neg_200805.csv')
# base_peak_name <- 'M482T929'
# base_peak_adduct <- '[M-H]-'
# test <- generatePeakGroup(peak_table = peak_table,
#                           base_peak_name = 'M482T929',
#                           base_peak_adduct = '[M-H]-',
#                           tol_rt = 3,
#                           cutoff_ssc = 0)


# load('~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/hyuA/raw_data_unlabel_231003.RData')
# load('~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ygeY/00_intemidiate_data/raw_data_unlabel.RData')
# peak_table <- raw_data_unlabel %>%
#   rename('name' = id) %>%
#   select(name, mz, rt, WT_UA_1:hyuA_UA_3)
# base_peak_name <- raw_data_unlabel %>% dplyr::filter(increase_label == 'signif_increase') %>% pull(name)
# base_peak_name <- '173.0668@3.511'
# base_peak_adduct <- '[M+H]+'
# path <- '~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/hyuA/hyuA_UA/'
# ms2_data_path <- '~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/hyuA/hyuA_UA/'

# test <- generatePeakGroup(peak_table = peak_table,
#                           base_peak_name = '173.0668@3.511',
#                           base_peak_adduct = '[M+H]+',
#                           tol_rt = 0.05,
#                           cutoff_ssc = 0.8)


setGeneric(name = 'generatePeakGroup',
           def = function(
    peak_table,
    base_peak_name,
    raw_data_files,
    base_peak_adduct = c('[M+H]+', '[M-H]-'),
    tol_mz = 20,
    tol_rt = 0.05,
    cutoff_ssc = 0.8, # cutoff of sample_sample_correlation
    cutoff_ppc = 0.8,
    is_calculate_ppc = TRUE,
    path = '.'
           ){
             # browser()
             base_peak <- peak_table %>%
               dplyr::filter(name == base_peak_name)

             base_rt <- base_peak$rt
             base_mz <- base_peak$mz
             base_int <- apply(base_peak, 1, function(x){
               ssc = mean(as.numeric(x[-c(1:3)]))
             })

             peak_list <- peak_table %>%
               dplyr::filter(rt >= (base_rt-tol_rt) & rt <= (base_rt+tol_rt))

             ssc <- apply(peak_list, 1, function(x){
               ssc = cor(as.numeric(x[-c(1:4)]),
                         as.numeric(base_peak[-c(1:4)]),
                         method = 'pearson')
             })

             int <- apply(peak_list, 1, function(x){
               ssc = mean(as.numeric(x[-c(1:4)]))
             })

             if (is_calculate_ppc) {
               ppc <- calculate_peak_peak_similarity(peak_list = peak_list,
                                                     base_peak_name = base_peak_name,
                                                     base_peak_mz = base_mz,
                                                     base_peak_rt = base_rt,
                                                     path = path,
                                                     files = raw_data_files,
                                                     mz_tol = tol_mz)
             } else {
               ppc <- rep(-1, nrow(peak_list))
             }


             peak_list <- peak_list %>%
               dplyr::mutate(mz = as.numeric(mz),
                             rt = as.numeric(rt)) %>%
               dplyr::mutate(ssc = ssc,
                             ppc = ppc,
                             int = int) %>%
               dplyr::select(name, mz, rt, ssc, ppc, int) %>%
               dplyr::rename(peak_name = name)

             if (cutoff_ssc > 0) {
               peak_list <- peak_list %>%
                 dplyr::filter(ssc >= cutoff_ssc)
             }

             if (cutoff_ppc > 0) {
               peak_list <- peak_list %>%
                 dplyr::filter(ppc >= cutoff_ppc)
             }

             result_peak_group <- new('PeakGroup',
                                      base_peak_name = base_peak_name,
                                      base_peak_mz = as.numeric(base_peak$mz),
                                      base_peak_rt = as.numeric(base_peak$rt),
                                      base_peak_int = as.numeric(base_int),
                                      base_peak_adduct = base_peak_adduct,
                                      group_size = nrow(peak_list),
                                      peak_list = peak_list)

             return(result_peak_group)
           })


setMethod(f = "show",
          signature = "PeakGroup",
          definition = function(object){
            # cat("-----------Meta information------------\n")
            cat("Base peak name:", object@base_peak_name, "\n")
            cat("m/z:", object@base_peak_mz, "\n")
            cat("RT:", object@base_peak_rt, "\n")
            # cat("CCS:", object@base_peak_ccs, "\n")
            cat("Average intensity:", object@base_peak_int, "\n")
            cat("Adduct:", object@base_peak_adduct, "\n")
            # cat("Seed:", object@base_peak_seed, "\n")
            cat("Number of peak in group:", nrow(object@peak_list), "\n")
            cat("Number of annotated peak in group:", nrow(object@peak_list_annotated), "\n")
          }
)


# initializePeakList -----------------------------------------------------------
# setClass(Class = 'PeakGroupAnnotation',
#          slots = list(peak_list_annotated = 'data.frame'),
#          contains = 'PeakGroup')

setGeneric(name = 'initializePeakList',
           def = function(
    peak_list
           ){
             peak_list <- peak_list %>%
               dplyr::mutate(mz_error = NA,
                             rt_error = NA,
                             # ccs_error = NA,
                             int_ratio = NA,
                             ms2_score = NA,
                             type = NA,
                             annotation = NA,
                             query_peak = NA)

             return(peak_list)
           })

# writePeakListAnnotation -----------------------------------------------------------
setGeneric(name = 'writePeakListAnnotation',
           def = function(pl_annotation_old,
                          pl_annotation_new){
             if (missing(pl_annotation_new)) {
               stop('Please input pl_annotation_new')
             }

             if (missing(pl_annotation_old)) {
               pl_annotation_old <- tibble::tibble(peak_name = character(),
                                                   mz = numeric(),
                                                   rt = numeric(),
                                                   # ccs = numeric(),
                                                   ssc = numeric(),
                                                   ppc = numeric(),
                                                   int = numeric(),
                                                   mz_error = numeric(),
                                                   rt_error = numeric(),
                                                   # ccs_error = numeric(),
                                                   int_ratio = numeric(),
                                                   ms2_score = numeric(),
                                                   type = character(),
                                                   annotation = character(),
                                                   query_peak = character())
             }

             result <- pl_annotation_old %>% dplyr::bind_rows(pl_annotation_new)
             return(result)
           })

# annotateIsotope --------------------------------------------------------------
# example L0076 neg in 200STD
# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/list_peak_group_200716.RData')
# peak_group <- list_peak_group$`M482T929_[M-H]-`
# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- peak_group@base_peak_name
# query_mz <- peak_group@base_peak_mz
# query_rt <- peak_group@base_peak_rt
# query_adduct <- peak_group@base_peak_adduct
# annotateIsotope(peak_list = peak_list,
#                 query_mz = query_mz,
#                 query_rt = query_rt,
#                 query_peak_name = query_peak_name,
#                 tol_mz = 10)

# peak_group <- list_peak_group$`M204T364_[M-H]-`
# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- 'M159T364'
# query_mz <- 159.0322
# query_rt <- 364.199
# annotateIsotope(peak_list = peak_list,
#                 query_mz = query_mz,
#                 query_rt = query_rt,
#                 query_peak_name = query_peak_name,
#                 tol_mz = 10,
#                 isotope_int_ratio_check = TRUE,
#                 isotope_int_ratio_cutoff = 500)


# peak_group <- list_peak_group_annotation$M313T547
# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- 'M313T547'
# query_mz <- 313.1431
# query_rt <- 547.098
# annotateIsotope(peak_list = peak_list,
#                 query_mz = query_mz,
#                 query_rt = query_rt,
#                 query_peak_name = query_peak_name,
#                 tol_mz = 10,
#                 isotope_int_ratio_check = TRUE,
#                 isotope_int_ratio_cutoff = 500)

#
# annotateIsotope(peak_list = peak_list,
#                 query_mz = 148.0381,
#                 query_rt = 65,
#                 query_ccs = 149.9,
#                 query_peak_name = "M148T65C150",
#                 tol_mz = 10)


setGeneric(name = 'annotateIsotope',
           def = function(
    peak_list,
    query_mz,
    query_rt,
    # query_ccs,
    query_peak_name,
    # query_adduct,
    tol_mz = 10,
    isotope_delta = 1.003355,
    isotope_max_num = 4,
    isotope_int_ratio_check = TRUE,
    isotope_int_ratio_cutoff = 500,
    monotonic_dec_check = FALSE,
    monotonic_mz_cutoff = 800,
    cutoff_ssc = NULL,
    cutoff_ssc_int = 3000,
    ...
           ){

             # browser()
             mz_isotope_list <- query_mz + seq(0, isotope_max_num-1)*isotope_delta
             isotope_matrix <- getMzRange(mz = mz_isotope_list,
                                          ppm = 25,
                                          mz_ppm_thr = 400)

             # keep peaks with ssc larger than cutoff only
             if (length(cutoff_ssc) > 0) {
               peak_list <- peak_list %>%
                 dplyr::filter((ssc >= cutoff_ssc) | (int <= cutoff_ssc_int))
             }

             result_isotope <- lapply(seq_along(mz_isotope_list), function(i){
               mz_min <- isotope_matrix[i,1]
               mz_max <- isotope_matrix[i,2]

               if (i==1) {
                 label <- "[M]"
               } else {
                 label <- paste0("[M+", i-1, "]")
               }

               result <- peak_list %>%
                 dplyr::filter(mz >= mz_min & mz<= mz_max) %>%
                 dplyr::mutate(mz_error = as.numeric(mz - mz_isotope_list[i]),
                               rt_error = as.numeric(rt - query_rt),
                               # ccs_error = as.numeric(ccs - query_ccs),
                               type = 'isotopeAnnotation',
                               annotation = label,
                               query_peak = query_peak_name)
             })

             result_isotope <- result_isotope %>% dplyr::bind_rows()

             # If multiple isotope has been annotated,
             #    reserve peak with minimum mz error
             result_isotope <- result_isotope %>%
               dplyr::arrange(annotation, abs(mz_error), abs(rt_error)) %>%
               dplyr::distinct(annotation, .keep_all = TRUE) %>%
               dplyr::mutate(int_ratio = int/int[1])

             if (isotope_int_ratio_check) {
               simulate_isotope <- simulateTheoIsotope(mz = query_mz,
                                                       isotope_max_num = isotope_max_num)

               simulate_isotope <- simulate_isotope %>%
                 dplyr::filter(label %in% result_isotope$annotation) %>%
                 dplyr::select(int_theo, label)

               temp_result_isotope <- result_isotope %>%
                 dplyr::left_join(simulate_isotope, by = c('annotation' = 'label'))

               result_isotope <- temp_result_isotope %>%
                 dplyr::rowwise() %>%
                 dplyr::mutate(int_error = abs(int_ratio-int_theo)/int_theo*100) %>%
                 dplyr::filter(int_error <= isotope_int_ratio_cutoff) %>%
                 dplyr::select(-c('int_theo', 'int_error')) %>%
                 tibble::as_tibble()

               # purrr::map2(result_isotope$int_ratio,
               #             simulate_isotope$int_theo,
               #             function(x, y){
               #
               #             })

             }

             # If monotonic_dec_check, the M+1 need to small the M.
             # Note: it is not suitable for molecules with halogen elements!
             if (monotonic_dec_check) {
               if (query_mz <= monotonic_mz_cutoff) {
                 idx <- sapply(seq_along(result_isotope$int_ratio),
                               function(i){
                                 if (i==1) return(TRUE)

                                 if (result_isotope$int_ratio[i] <
                                     result_isotope$int_ratio[i-1]) {
                                   return(TRUE)
                                 } else {
                                   return(FALSE)
                                 }
                               })

                 result_isotope <- result_isotope[idx,]
               }
             }

             return(result_isotope)
           })



setGeneric(name = 'getMzRange',
           def = function(mz,
                          ppm = 10,
                          mz_ppm_thr = 500){
             result <- sapply(mz, function(x) {
               if (x >= mz_ppm_thr) {
                 x * (1 + c(-1, 1) * ppm * 1e-6)
               } else {
                 temp1 <- x + mz_ppm_thr * c(-1, 1) * ppm * 1e-6
               }
             })

             t(result)
           }
)


# simulateTheoIsotope(mz = 481.9783, isotope_max_num = 4)
setGeneric(name = 'simulateTheoIsotope',
           function(
    mz,
    isotope_max_num = 4
           ){
             # calculate simulated carbon number with alkane (CnH2n+2)
             num_carbon <- (mz - 1.0078*2) %/% 14.0156
             simulate_alkane_formula <- paste0('C', num_carbon, 'H', 2*num_carbon+2)

             options(readr.num_columns = 0)
             simulate_isotope <- Rdisop::getMolecule(simulate_alkane_formula, z=1) %>%
               Rdisop::getIsotope(index = seq(isotope_max_num)) %>%
               t() %>%
               tibble::as_tibble() %>%
               dplyr::rename(mz = V1, int_theo = V2)

             if (isotope_max_num > 1) {
               label <- c('[M]', paste0("[M+", seq(isotope_max_num-1), "]"))
             } else {
               label <- '[M]'
             }

             simulate_isotope <- simulate_isotope %>%
               dplyr::mutate(label = label) %>%
               dplyr::mutate(int_theo = int_theo/int_theo[1])

             return(simulate_isotope)

           })


# annotateAdduct ---------------------------------------------------------------

# annotateAdduct(peak_list = peak_list,
#                query_mz = query_mz,
#                query_rt = query_rt,
#                query_peak_name = query_peak_name,
#                query_adduct = query_adduct,
#                polarity = polarity,
#                tol_mz = 25,
#                cutoff_ssc = 0.3,
#                cutoff_ssc_int = 3000)


# load('~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/hyuA/raw_data_unlabel_231003.RData')
# peak_table <- raw_data_unlabel %>%
#   rename('name' = id) %>%
#   select(name, mz, rt, WT_UA_1:hyuA_UA_3)
# base_peak_name <- raw_data_unlabel %>% dplyr::filter(increase_label == 'signif_increase') %>% pull(name)
# base_peak_name <- '173.0668@3.511'
# base_peak_adduct <- '[M+H]+'

# test <- generatePeakGroup(peak_table = peak_table,
#                           base_peak_name = '173.0668@3.511',
#                           base_peak_adduct = '[M+H]+',
#                           tol_rt = 0.05,
#                           cutoff_ssc = 0.8)

# annotateAdduct(peak_list = test@peak_list,
#                query_mz = test@base_peak_mz,
#                query_rt = test@base_peak_rt,
#                query_peak_name = test@base_peak_name,
#                query_adduct = test@base_peak_adduct,
#                polarity = 'positive',
#                tol_mz = 25)

# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- 'M159T364'
# query_mz <- 173.0668
# query_rt <- 3.511
# query_adduct <- '[M+H]+'
# polarity <- 'positive'

setGeneric(name = 'annotateAdduct',
           def = function(
    peak_list,
    query_mz,
    query_rt,
    query_peak_name,
    query_adduct,
    polarity = c('positive', 'negative'),
    tol_mz = 25,
    cutoff_ssc = NULL,
    cutoff_ssc_int = 3000,
    is_rule_limitation = TRUE,
    # lib_adduct_nl = lib_adduct_nl,
    ...
           ){
             # browser()

             adduct_list <- convertMz2Adduct(base_mz = query_mz,
                                             base_adduct = query_adduct,
                                             type = 'adduct',
                                             polarity = polarity)

             mz_adduct_list <- adduct_list$mz
             name_adduct_list <- adduct_list$adduct

             adduct_matrix <- getMzRange(mz = mz_adduct_list,
                                         ppm = tol_mz,
                                         mz_ppm_thr = 400)

             # keep peaks with ssc larger than cutoff only (except int less than cutoff)
             if (length(cutoff_ssc) > 0) {
               peak_list <- peak_list %>%
                 dplyr::filter((ssc >= cutoff_ssc) | (int <= cutoff_ssc_int))
             }

             result_adduct <- lapply(seq_along(mz_adduct_list), function(i){
               mz_min <- adduct_matrix[i,1]
               mz_max <- adduct_matrix[i,2]

               result <- peak_list %>%
                 dplyr::filter(mz >= mz_min & mz<= mz_max) %>%
                 dplyr::mutate(mz_error = as.numeric(mz - mz_adduct_list[i]),
                               rt_error = as.numeric(rt - query_rt),
                               # ccs_error = as.numeric(ccs - query_ccs),
                               type = 'adductAnnotation',
                               annotation = name_adduct_list[i],
                               query_peak = query_peak_name)
             })

             result_adduct <- result_adduct %>% dplyr::bind_rows()

             # If multiple adduct has been annotated,
             #    reserve peak with minimum mz error
             result_adduct <- result_adduct %>%
               dplyr::arrange(annotation, abs(mz_error), abs(rt_error)) %>%
               dplyr::distinct(annotation, .keep_all = TRUE)

             # Check rule requirements for adducts
             if (is_rule_limitation) {
               rule_list <- generateRuleList(polarity = polarity)

               idx <- which(result_adduct$annotation %in% rule_list$rule_adduct)
               if (length(idx) > 0) {
                 idx_eff <- sapply(seq_along(idx), function(i){
                   temp_rule <- rule_list %>%
                     dplyr::filter(rule_adduct == result_adduct$annotation[idx[i]]) %>%
                     dplyr::pull(rule)

                   result <- all(temp_rule %in% result_adduct$annotation)
                 })

                 if (!all(idx_eff)) {
                   adduct_rm <- idx[!idx_eff] %>%
                     result_adduct$annotation[.]

                   result_adduct <- result_adduct %>%
                     dplyr::filter(!(annotation %in% adduct_rm))
                 }
               }

             }

             return(result_adduct)

           })


# lib_adduct_pos <- readr::read_csv('./lib_adduct_pos_200713_manual.csv')
# lib_adduct_neg <- readr::read_csv('./lib_adduct_neg_200713_manual.csv')
# lib_nl_pos <- readr::read_csv('./lib_nl_pos_200713_manual.csv')
# lib_nl_neg <- readr::read_csv('./lib_nl_neg_200713_manual.csv')
#
# lib_pos <- lib_adduct_pos %>% bind_rows(lib_nl_pos)
# lib_neg <- lib_adduct_neg %>% bind_rows(lib_nl_neg)
#
# lib_adduct_nl <- list(positive = lib_pos,
#                       negative = lib_neg)
#
# save(lib_adduct_nl, file = './lib_adduct_nl_200714.RData', version = 2)


# calculateExactMass(formula = "C2H5OH")

setGeneric(name = 'calculateExactMass',
           def = function(
    formula
           ){
             molecule <- Rdisop::getMolecule(formula)
             # getFormula(molecule)
             Rdisop::getMass(molecule)
           })


# exact_mass <- 180.0634
# adduct <- '[M-H]-'
# delta_mz <- -1.0073
# calculateMz(exact_mass = 180.0634,
#             adduct = '[M-H]-',
#             delta_mz = -1.0073)

# calculateMz(exact_mass = 180.0634,
#             adduct = '[2M-H]-',
#             delta_mz = -1.0073)
#
# calculateMz(exact_mass = 180.0634,
#             adduct = '[M-2H]2-',
#             delta_mz = -1.0073)


setGeneric(name = 'calculateMz',
           def = function(
    exact_mass,
    adduct,
    delta_mz,
    nmol = NULL,
    ncharge = NULL
           ){

             if (length(nmol) == 0) {
               if (stringr::str_detect(adduct, pattern = '2M')) {
                 mz <- exact_mass*2 + delta_mz
               } else if (stringr::str_detect(adduct, pattern = '3M')) {
                 mz <- exact_mass*3 + delta_mz
               } else {
                 mz <- exact_mass + delta_mz
               }
             } else {
               mz <- exact_mass*nmol + delta_mz
             }


             if (length(ncharge) == 0) {
               if (stringr::str_detect(adduct, pattern = '\\]2\\-|\\]2\\+')) {
                 mz <- mz/2
               } else if (stringr::str_detect(adduct, pattern = '\\]3\\-|\\]3\\+')) {
                 mz <- mz/3
               } else {
                 mz
               }
             } else {
               mz <- mz/ncharge
             }

             mz
           })


# exact_mass <- 180.0634
# adduct <- '[M-H]-'
# delta_mz <- -1.0073
# transformMz(exact_mass = 180.0634, type = 'adduct', polarity = 'positive')
# transformMz(exact_mass = 180.0634, type = 'nl', polarity = 'positive')
# transformMz(exact_mass = 180.0634, adduct_list = c('[M-H]-', '[M+H]+'))

setGeneric(name = 'transformMz',
           function(
    exact_mass,
    formula = NULL,
    adduct_list = NULL,
    type = c('adduct', 'nl'),
    polarity = c('positive', 'negative'),
    # lib_adduct_nl = lib_adduct_nl,
    ...
           ){

             if (all(is.null(exact_mass), is.null(formula))) {
               stop('Please input exact_mass or formula.')
             }

             if (!is.null(formula)) {
               exact_mass <- calculateExactMass(formula)
             }

             if (is.null(adduct_list)) {
               lib <- switch(polarity,
                             'positive' = {
                               lib_adduct_nl$positive
                             },
                             'negative' = {
                               lib_adduct_nl$negative
                             }
               )

               lib <- switch(type,
                             'adduct' = {
                               lib %>%
                                 dplyr::filter(type == 'Adduct') %>%
                                 dplyr::filter(credential == 'Yes')
                             },
                             'nl' = {
                               lib %>%
                                 dplyr::filter(type == 'NeutralLoss') %>%
                                 dplyr::filter(credential == 'Yes')
                             })
             } else {
               lib <- lib_adduct_nl$positive %>%
                 dplyr::bind_rows(lib_adduct_nl$negative)

               if (!all(adduct_list %in% lib$adduct)) {
                 stop('Sorry, not all adduct included in the adduct list\n')
               }

               lib <- lib %>%
                 dplyr::filter(adduct %in% adduct_list) %>%
                 dplyr::arrange(match(adduct, adduct_list))
             }


             result_mz <- sapply(seq_along(lib$adduct), function(i){
               calculateMz(exact_mass = exact_mass,
                           adduct = lib$adduct[i],
                           delta_mz = lib$delta_mz[i])
             })

             result <- tibble::tibble(exact_mass = exact_mass,
                                      adduct = lib$adduct,
                                      mz = result_mz)


             return(result)

           })




# convertMz2Adduct(base_mz = 181.0707,
#                  base_adduct = '[M+H]+',
#                  adduct_list = NULL,
#                  type = 'adduct',
#                  polarity = 'positive')

# convertMz2Adduct(base_mz = 181.0707,
#                  base_adduct = '[M+H]+',
#                  adduct_list = NULL,
#                  type = 'nl',
#                  polarity = 'negative')

# convertMz2Adduct(base_mz = 181.0707,
#                  base_adduct = '[M+H]+',
#                  adduct_list = c('[M-H2O+H]+', '[M+Na]+'),
#                  type = 'nl',
#                  polarity = 'negative')

# convertMz2Adduct(base_mz = 143.0347,
#                  base_adduct = '[2M-H]-',
#                  adduct_list = c('[M-H2O+H]+', '[M+Na]+'))

# convertMz2Adduct(base_mz = 143.0347,
#                  base_adduct = '[3M-H]-',
#                  adduct_list = c('[M-H2O+H]+', '[M+Na]+'))
# convertMz2Adduct(base_mz = 664.1144, base_adduct = '[M]+', adduct_list = '[M]+')$exact_mass

setGeneric(name = 'convertMz2Adduct',
           def = function(
    base_mz,
    base_adduct,
    # lib_adduct_nl = lib_adduct_nl,
    ...
    # adduct_list = NULL,
    # type = c('adduct', 'nl'),
    # polarity = c('positive', 'negative')
           ){

             lib <- lib_adduct_nl$positive %>%
               dplyr::bind_rows(lib_adduct_nl$negative)

             if (!(base_adduct %in% lib$adduct)) {
               stop('Sorry, base_adduct is not included\n')
             }

             temp_delta_mz <- lib %>%
               dplyr::filter(adduct == base_adduct) %>%
               dplyr::pull(delta_mz)


             if (stringr::str_detect(base_adduct, pattern = '2M')) {
               temp_exact_mass <- (base_mz - temp_delta_mz)/2
             } else if (stringr::str_detect(base_adduct, pattern = '3M')) {
               temp_exact_mass <- (base_mz - temp_delta_mz)/3
             } else {
               temp_exact_mass <- base_mz - temp_delta_mz
             }

             # temp_exact_mass <- base_mz - temp_delta_mz

             result <- transformMz(exact_mass = temp_exact_mass,
                                   ...)

             # if initial seed base_adduct is '[M]+', '[M-2H]-', the adduct would be added for credential
             if ((base_adduct %in% c('[M]+', '[M-2H]-')) & !(base_adduct %in% result$adduct)) {
               temp_result <- transformMz(exact_mass = temp_exact_mass, adduct_list = base_adduct)
               result <- temp_result %>% dplyr::bind_rows(result)
               return(result)
             }

             return(result)

           })


setGeneric(name = 'generateRuleList',
           def = function(
    polarity = c('positive', 'negative')
           ){
             lib <- switch(polarity,
                           'positive' = {
                             lib_adduct_nl$positive
                           },
                           'negative' = {
                             lib_adduct_nl$negative
                           }
             )

             result <- lib %>%
               dplyr::filter(!is.na(rule_limitation)) %>%
               dplyr::select(adduct, rule_limitation) %>%
               tidyr::separate_rows(rule_limitation, sep = ';') %>%
               dplyr::rename(rule_adduct = adduct,
                             rule = rule_limitation)

             return(result)
           })


# annotateNeutralLoss ----------------------------------------------------------

# polarity <- 'positive'
# annotateNeutralLoss(peak_list = test@peak_list,
#                     # peak_list_ms2 = peak_list_ms2,
#                     query_mz = test@base_peak_mz,
#                     query_rt = test@base_peak_rt,
#                     query_peak_name = test@base_peak_name,
#                     query_adduct = test@base_peak_adduct,
#                     polarity = 'positive',
#                     tol_mz = 25,
#                     is_ms2_check = FALSE,
#                     is_rule_limitation = TRUE)

# annotateNeutralLoss(peak_list = test@peak_list,
#                     peak_list_ms2 = peak_list_ms2,
#                     query_mz = query_mz,
#                     query_rt = query_rt,
#                     query_peak_name = query_peak_name,
#                     query_adduct = query_adduct,
#                     polarity = polarity,
#                     tol_mz = 25,
#                     is_ms2_check = TRUE,
#                     is_rule_limitation = TRUE)



setGeneric(name = 'annotateNeutralLoss',
           def = function(
    peak_list,
    peak_list_ms2,
    query_peak_name,
    query_mz,
    query_rt,
    # query_ccs,
    query_adduct,
    polarity = c('positive', 'negative'),
    tol_mz = 25,
    is_ms2_check = FALSE,
    ms2_score_cutoff = -1,
    cutoff_ssc = NULL,
    cutoff_ssc_int = 3000,
    is_rule_limitation = TRUE,
    # lib_adduct_nl = lib_adduct_nl,
    ...
           ){
             nl_list <- convertMz2Adduct(base_mz = query_mz,
                                         base_adduct = query_adduct,
                                         type = 'nl',
                                         polarity = polarity,
                                         lib_adduct_nl = lib_adduct_nl)

             mz_nl_list <- nl_list$mz
             name_nl_list <- nl_list$adduct

             nl_matrix <- getMzRange(mz = mz_nl_list,
                                     ppm = tol_mz,
                                     mz_ppm_thr = 400)

             # keep peaks with ssc larger than cutoff only (except intensity less than cutoff)
             if (length(cutoff_ssc) > 0) {
               peak_list <- peak_list %>%
                 dplyr::filter((ssc >= cutoff_ssc) | (int <= cutoff_ssc_int))
             }

             result_nl <- lapply(seq_along(mz_nl_list), function(i){
               mz_min <- nl_matrix[i,1]
               mz_max <- nl_matrix[i,2]

               result <- peak_list %>%
                 dplyr::filter(mz >= mz_min & mz<= mz_max) %>%
                 dplyr::mutate(mz_error = as.numeric(mz - mz_nl_list[i]),
                               rt_error = as.numeric(rt - query_rt),
                               # ccs_error = as.numeric(ccs - query_ccs),
                               type = 'neutralLossAnnotation',
                               annotation = name_nl_list[i],
                               query_peak = query_peak_name)
             })

             result_nl <- result_nl %>% dplyr::bind_rows()

             # If multiple nl has been annotated,
             #    reserve peak with minimum mz error
             result_nl <- result_nl %>%
               dplyr::arrange(annotation, abs(mz_error), abs(rt_error)) %>%
               dplyr::distinct(annotation, .keep_all = TRUE)

             if (is_ms2_check) {
               if (missing(peak_list_ms2)) {
                 stop('Please input msms check')
               }

               if (nrow(result_nl) == 0) {
                 return(result_nl)
               }

               # if the annotation don't has ms/ms, return result
               temp_idx <- which(names(peak_list_ms2) == query_peak_name)

               if (length(temp_idx) == 0) {
                 return(result_nl)
               } else {
                 base_ms2 <- peak_list_ms2[[temp_idx]]
               }


               result_ms2_score <- lapply(seq_along(result_nl$peak_name),
                                          function(i){
                                            # cat(i, ' ')
                                            nl_mz <- result_nl$mz[i]

                                            # purify ms2 spectra and remove fragments >= NL precursor
                                            base_ms2$spec <- purifyMs2(spec = base_ms2$spec,
                                                                       is_include_precursor = TRUE,
                                                                       mz_range_ms2 = c(0, nl_mz+0.1),
                                                                       is_deisotope = FALSE,
                                                                       int_ms2_min_abs = 30,
                                                                       int_ms2_min_relative = 0.01,
                                                                       mz_precursor = nl_mz,
                                                                       ppm_precursor_filter = 10)

                                            # if no fragment mz less than nl mz, return NA
                                            if (length(base_ms2$spec) == 0) {
                                              return(NA)
                                            }

                                            temp_nl_name <- result_nl$peak_name[i]
                                            temp_nl_ms2 <- try(which(names(peak_list_ms2) == temp_nl_name) %>%
                                                                 peak_list_ms2[[.]], silent = TRUE)

                                            # if no nl ms2, return NA
                                            if (class(temp_nl_ms2) == 'try-error') {
                                              return(NA)
                                            }

                                            # purify ms2 spectra and remove fragments >= NL precursor
                                            temp_nl_ms2$spec <- purifyMs2(spec = temp_nl_ms2$spec,
                                                                          is_include_precursor = TRUE,
                                                                          mz_range_ms2 = c(0, nl_mz+0.1),
                                                                          is_deisotope = FALSE,
                                                                          int_ms2_min_abs = 30,
                                                                          int_ms2_min_relative = 0.01,
                                                                          mz_precursor = nl_mz,
                                                                          ppm_precursor_filter = 10)

                                            # if no fragment mz less than nl mz, return NA
                                            if (length(temp_nl_ms2$spec) == 0) {
                                              return(NA)
                                            }

                                            # modify the msms format for SpectraTools
                                            base_ms2 <- convertSpectraData(ms2_data = base_ms2)
                                            temp_nl_ms2 <- convertSpectraData(ms2_data = temp_nl_ms2)


                                            matchParam <- SpectraTools::MatchParam(ppm = 10,
                                                                                   methodScore = 'dp',
                                                                                   methodMatch = 'direct',
                                                                                   weightIntensity = 1,
                                                                                   weightMZ = 0,
                                                                                   cutoff = 0,
                                                                                   includePrecursor = TRUE,
                                                                                   intensityExpNormed = TRUE,
                                                                                   intensityLibNormed = TRUE,
                                                                                   tuneLibSpectra = FALSE)

                                            result <- try(SpectraTools::MatchSpectra(base_ms2,
                                                                                     temp_nl_ms2,
                                                                                     matchParam),
                                                          silent = TRUE)

                                            if ((class(result) == 'try-error') | (length(result) == 0)) {
                                              return(NA)
                                            }

                                            result <- result@info %>% tibble::as_tibble() %>% dplyr::pull(scoreReverse)

                                          })

               result_nl <- result_nl %>%
                 dplyr::mutate(ms2_score = unlist(result_ms2_score)) %>%
                 dplyr::filter((ms2_score >= ms2_score_cutoff) | is.na(ms2_score))

               # Check rule requirements for adducts
               if (is_rule_limitation) {
                 rule_list <- generateRuleList(polarity = polarity)

                 idx <- which(result_nl$annotation %in% rule_list$rule_adduct)
                 if (length(idx) > 0) {
                   idx_eff <- sapply(seq_along(idx), function(i){
                     temp_rule <- rule_list %>%
                       dplyr::filter(rule_adduct == result_nl$annotation[idx[i]]) %>%
                       dplyr::pull(rule)

                     result <- all(temp_rule %in% result_nl$annotation)
                   })

                   if (!all(idx_eff)) {
                     nl_rm <- idx[!idx_eff] %>%
                       result_nl$annotation[.]

                     result_nl <- result_nl %>%
                       dplyr::filter(annotation != nl_rm)
                   }
                 }

               }


             }

             return(result_nl)

           })


#' @title convertSpectraData
#' @param ms2_data
#' @importClassesFrom SpectraTools 'SpectraData'
#' @export

setGeneric(name = 'convertSpectraData',
           def = function(
    ms2_data
           ){
             options(readr.num_columns = 0)
             temp_info <- ms2_data$info %>%
               dplyr::rename(name = NAME,
                             mz = PRECURSORMZ) %>%
               dplyr::select(name:mz) %>%
               readr::type_convert()

             temp_ms2_data <- ms2_data$spec

             result <- new('SpectraData',
                           info = temp_info,
                           spectra = list(temp_ms2_data))

             return(result)
           })







# annotateISF ------------------------------------------------------------------

# annotateISF via ms2 & correlation --------------------------------------------


# peak_group <- test
# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- peak_group@base_peak_name
# query_mz <- peak_group@base_peak_mz
# query_rt <- peak_group@base_peak_rt
# query_adduct <- peak_group@base_peak_adduct
# raw_msms <- ms2_data_combined$spec
# peak_list_ms2 <- match(peak_list$peak_name, names(raw_msms)) %>%
#   raw_msms[.] %>%
#   .[!is.na(names(.))]
# polarity <- 'positive'
# annotateISF(peak_list = peak_list,
#             peak_list_ms2 = peak_list_ms2,
#             query_peak_name = query_peak_name,
#             query_mz = query_mz,
#             query_rt = query_rt,
#             tol_mz = 10)

setGeneric(name = 'annotateISF',
           function(
    peak_list,
    peak_list_ms2,
    query_peak_name,
    query_mz,
    query_rt,
    # query_ccs,
    tol_mz = 20,
    cutoff_ssc = NULL,
    cutoff_ssc_int = 3000,
    cutoff_topN = 5,
    cutoff_ppc = 0.8,
    ...
           ){
             # annotate ISF using peak shape PPC score

             # annotate ISF using ms2
             # if the annotation don't has ms/ms, return result
             temp_idx <- which(names(peak_list_ms2) == query_peak_name)

             if (length(temp_idx) > 0) {
               base_ms2 <- which(names(peak_list_ms2) == query_peak_name) %>% peak_list_ms2[[.]]

               # limit topN fragments in base peak for ISF annotation
               # mz_isf_list <- base_ms2[,'mz']
               mz_isf_list <- order(base_ms2[,'intensity'], decreasing = TRUE) %>%
                 base_ms2[.,,drop = FALSE] %>%
                 .[,'mz']

               if (length(mz_isf_list) >= cutoff_topN) {
                 mz_isf_list <- mz_isf_list[1:cutoff_topN]
               }

               isf_matrix <- getMzRange(mz = mz_isf_list,
                                        ppm = 25,
                                        mz_ppm_thr = 400)

               # keep peaks with ssc larger than cutoff only (except intensity less than cutoff)
               if (length(cutoff_ssc) > 0) {
                 peak_list <- peak_list %>%
                   dplyr::filter((ssc >= cutoff_ssc) | (int <= cutoff_ssc_int))
               }

               result_isf <- lapply(seq_along(mz_isf_list), function(i){
                 mz_min <- isf_matrix[i,1]
                 mz_max <- isf_matrix[i,2]

                 result <- peak_list %>%
                   dplyr::filter(mz >= mz_min & mz<= mz_max) %>%
                   dplyr::mutate(mz_error = as.numeric(mz - mz_isf_list[i]),
                                 rt_error = as.numeric(rt - query_rt),
                                 # ccs_error = as.numeric(ccs - query_ccs),
                                 type = 'isfAnnotation',
                                 annotation = 'ISF',
                                 query_peak = query_peak_name)
               })

               result_isf <- result_isf %>% dplyr::bind_rows()

             } else {
               return(NULL)
             }


             # If multiple isf has been annotated,
             #    reserve peak with minimum mz error
             result_isf <- result_isf %>%
               dplyr::arrange(annotation, abs(mz_error), abs(rt_error)) %>%
               dplyr::distinct(peak_name, .keep_all = TRUE)

             return(result_isf)
           })





# retrieve_ms2_data ------------------------------------------------------------
# ms2_type <- 'mzML'

# retrieve_ms2_data(peak_table = raw_data_unlabel,
#                   path = '~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/hyuA/',
#                   ms2_type = 'mzML')

retrieve_ms2_data <- function(peak_table,
                              path = '.',
                              ms2_type = c('mzML', 'mzXML'),
                              is_include_precursor = TRUE,
                              int_ms2_min_abs = 50,
                              int_ms2_min_relative = 0.01,
                              mz_tol_combine_ms1_ms2 = 20,
                              rt_tol_combine_ms1_ms2 = 20
) {
  cat('Retrive, purify, and combine MS2 with MS1\n')
  # browser()
  # convert enriched table to a ms1_data

  if (max(peak_table$rt) < 60) {
    variable_data <- peak_table %>%
      dplyr::select(name:rt) %>%
      dplyr::mutate(rt = rt*60) %>%
      as.data.frame()
  } else {
    variable_data <- peak_table %>%
      dplyr::select(name:rt) %>%
      dplyr::mutate(rt = rt) %>%
      as.data.frame()
  }


  expression_profile_data <- peak_table %>%
    dplyr::select(-c('name':'rt')) %>%
    as.data.frame()

  ms1_data <- list(info = variable_data, subject = expression_profile_data)
  rm(variable_data, expression_profile_data);gc()

  ms2_file <- list.files(file.path(path, 'ms2'), pattern = 'mzML', recursive = TRUE)
  ms2_data <- DoddLabMetID::read_ms2(file.path(path, 'ms2', ms2_file), ms2_type = 'mzML')

  ms2_data <- DoddLabMetID::integrate_ms2(ms2_data = ms2_data,
                                          ms2_file = file.path(path, 'ms2', ms2_file),
                                          ms2_type = ms2_type,
                                          is_include_precursor = TRUE,
                                          is_deisotope = FALSE,
                                          int_ms2_min_abs = int_ms2_min_abs,
                                          int_ms2_min_relative = int_ms2_min_relative)


  ms2_data_combined <- DoddLabMetID::combine_ms1_ms2(ms1_data = ms1_data,
                                                     ms2_data = ms2_data,
                                                     ms2_type = ms2_type,
                                                     mz_tol_combine_ms1_ms2 = mz_tol_combine_ms1_ms2,
                                                     rt_tol_combine_ms1_ms2 = rt_tol_combine_ms1_ms2)

  return(ms2_data_combined)
}


# calculate_peak_peak_similarity -----------------------------------------------

# temp_files <- list.files(path = file.path(path, 'hyuA_UA'), recursive = TRUE)
# calculate_peak_peak_similarity(peak_list = peak_list,
#                                base_peak_name = '173.0668@3.511',
#                                base_peak_mz = 173.0668,
#                                base_peak_rt = 3.511,
#                                path = file.path(path, 'hyuA_UA'),
#                                files = file.path(path, 'hyuA_UA', temp_files),
#                                mz_tol = 20)

calculate_peak_peak_similarity <- function(peak_list,
                                           base_peak_name = '173.0668@3.511',
                                           base_peak_mz = 173.0668,
                                           base_peak_rt = 3.511,
                                           path,
                                           files = NULL,
                                           files_pattern = NULL, # 'ygeX_UA.+\\.mzML
                                           mz_tol = 10) {
  # browser()

  # temp_files <- list.files(path = file.path(path, 'hyuA_UA'), recursive = TRUE)
  # eic_data_raw <- extract_eic_data(path = path,
  #                                  files = file.path('hyuA_UA', temp_files),
  #                                  mz_list = peak_list$mz,
  #                                  mz_tol = 20)

  temp_files <- list.files(path = file.path(path), recursive = TRUE)
  eic_data_raw <- extract_eic_data(path = path,
                                   files = temp_files,
                                   mz_list = peak_list$mz,
                                   mz_tol = mz_tol)

  model_file <- eic_data_raw %>% dplyr::count(filename) %>% dplyr::arrange(desc(n)) %>% dplyr::slice(1) %>% dplyr::pull(filename)

  # smooth base peak EIC
  eic_base_data <- eic_data_raw %>%
    dplyr::filter(filename == model_file) %>%
    dplyr::filter(precursor_mz == base_peak_mz) %>%
    # dplyr::filter(rt >= 3.411 & rt <= 3.611) %>%
    dplyr::filter(rt >= base_peak_rt - 0.15 & rt <=  base_peak_rt + 0.15)


  # if the peak is not continuous, add the boundary value for smooth model fitting
  # if the closed point and boundary RT has RT <= 0.01, this peak is continuous; otherwise, it is not continuous

  # if (abs(base_peak_rt - 0.15 - min(eic_base_data$rt)) > 0.01) {
  #   eic_base_data <- eic_base_data %>%
  #     dplyr::add_row(rt = base_peak_rt - 0.15, mz = base_peak_mz, int = 100, filename = model_file, label = model_file, precursor_mz = as.character(base_peak_mz)) %>%
  #     dplyr::arrange(rt)
  # }
  #
  # if (abs(base_peak_rt + 0.15 - max(eic_base_data$rt)) > 0.01) {
  #   eic_base_data <- eic_base_data %>%
  #     dplyr::add_row(rt = base_peak_rt + 0.15, mz = base_peak_mz, int = 100, filename = model_file, label = model_file, precursor_mz = as.character(base_peak_mz)) %>%
  #     dplyr::arrange(rt)
  # }
  #
  # eic_base_smooth <- smooth_peak(eic_base_data)



  list_smoothed_peaks <- lapply(seq_along(peak_list$mz), function(j){
    # cat(j, ' ')
    x <- peak_list$mz[j]
    eic_data <- eic_data_raw %>%
      dplyr::filter(filename == model_file) %>%
      dplyr::filter(precursor_mz == x) %>%
      # dplyr::filter(rt >= 3.411 & rt <= 3.611) %>%
      dplyr::filter(rt >= base_peak_rt - 0.15 & rt <=  base_peak_rt + 0.15)

    if (nrow(eic_data) < 5) {
      return(NULL)
    }

    # if the peak is not continuous, add the boundary value for smooth model fitting
    # if the closed point and boundary RT has RT <= 0.01, this peak is continuous; otherwise, it is not continuous
    if (abs(base_peak_rt - 0.15 - min(eic_data$rt)) > 0.01) {
      eic_data <- eic_data %>%
        dplyr::add_row(rt = base_peak_rt - 0.15, mz = base_peak_mz, int = 0, filename = model_file, label = model_file, precursor_mz = as.character(base_peak_mz)) %>%
        dplyr::arrange(rt)
    }

    if (abs(base_peak_rt + 0.15 - max(eic_data$rt)) > 0.01) {
      eic_data <- eic_data %>%
        dplyr::add_row(rt = base_peak_rt + 0.15, mz = base_peak_mz, int = 0, filename = model_file, label = model_file, precursor_mz = as.character(base_peak_mz)) %>%
        dplyr::arrange(rt)
    }

    eic_smooth <- smooth_peak(eic_data, range = c(base_peak_rt-0.1, base_peak_rt+0.1), step = 0.003, span = 0.25)
  })

  names(list_smoothed_peaks) <- peak_list$name

  eic_smooth_base <- list_smoothed_peaks[[base_peak_name]]

  # if the base peak eic is null,
  #   assign 0 as ppc score for all peaks
  if (!is.null(eic_smooth_base)){
    result_ppc <- sapply(list_smoothed_peaks, function(x){
      if (is.null(x)) {
        return(0)
      } else {
        cor(eic_smooth_base$int, x$int, method = 'pearson')
      }
    })
  } else {
    result_ppc <- rep(0, length(list_smoothed_peaks))
  }

  return(result_ppc)
}

# smooth_peak(eic_data = eic1_data,
#             range = c(3.411, 3.611))

smooth_peak <- function(eic_data,
                        range, # predict RT range
                        step = 0.003, # min
                        span = 0.25,
                        degree = 1,
                        family = 'gaussian') {
  temp_eic <- loess(int~rt, data = eic_data, span = span, degree = degree, family = 'gaussian')
  temp_eic_y <- predict(object = temp_eic, newdata = seq(range[1], range[2], step))

  result <- data.frame(rt = seq(range[1], range[2], step),
                       int = temp_eic_y,
                       stringsAsFactors = FALSE)

  return(result)
}


# eic1 <- eic_data_unlabel %>%
#   dplyr::filter(precursor_mz == 173.0668) %>%
#   # dplyr::filter(rt >= 3.411 & rt <= 3.611) %>%
#   dplyr::filter(rt >= 3.3 & rt <= 3.7) %>%
#   dplyr::filter(label == 'S110.mzML')
#
# eic2 <- eic_data_unlabel %>%
#   dplyr::filter(precursor_mz == 130.061) %>%
#   # dplyr::filter(rt >= 3.411 & rt <= 3.611) %>%
#   dplyr::filter(rt >= 3.3 & rt <= 3.7) %>%
#   dplyr::filter(label == 'S110.mzML')
#
#
# eic_data_unlabel %>%
#   dplyr::filter(precursor_mz == 173.0668) %>%
#   dplyr::filter(rt >= 3.411 & rt <= 3.611) %>%
#   # dplyr::filter(label == 'S110.mzML') %>%
#   ggplot() +
#   geom_line(aes(x=rt, y=int, color=filename))
#
# temp_eic1 <- loess(int~rt, data = eic1, span = 0.25, degree = 1, family = 'gaussian')
# temp_eic1_y <- predict(object = temp_eic1, newdata = seq(3.47, 3.6, 0.003))
#
# plot(seq(3.47, 3.6, 0.003), temp_eic1_y, type = 'l')
#
# temp_eic2 <- loess(int~rt, data = eic2, span = 0.25, family = 'gaussian')
# temp_eic2_y <- predict(object = temp_eic2, newdata = seq(3.47, 3.6, 0.003))
#
# # plot(temp_eic2, type = 'l')
# plot(seq(3.47, 3.6, 0.003), temp_eic2_y, type = 'l')




################################################################################
# merge_peak_groups ------------------------------------------------------------

merge_peak_groups <- function(list_peak_group_annotation,
                              path_dir = '.',
                              is_save_intermidate = TRUE) {
  temp <- statNumListPeakGroup(list_peak_group_annotation)
  num_stat <- temp

  temp <- statTypeListPeakGroup(list_peak_group_annotation)
  type_stat <- temp

  # 01 remove peak groups not meet the rule requirement
  cat('Remove peak_group filtered by rules\n')
  idx_rm <- sapply(list_peak_group_annotation, function(x){
    nrow(x@peak_list_annotated)
  })
  idx_rm <- which(idx_rm == 0)
  if (length(idx_rm) > 0) {
    list_peak_group_annotation_concised <- list_peak_group_annotation[-idx_rm]
    record_rule_filter_peak_group <- names(list_peak_group_annotation)[idx_rm]

    if (is_save_intermidate) {
      save(record_rule_filter_peak_group,
           file = file.path(path_dir,
                            '00_tracer_result',
                            "00_intermediate_data",
                            'record_rule_filter_peak_group.RData'))
    }
  } else {
    list_peak_group_annotation_concised <- list_peak_group_annotation
    record_rule_filter_peak_group <- NULL

    if (is_save_intermidate) {
      save(record_rule_filter_peak_group,
           file = file.path(path_dir,
                            '00_tracer_result',
                            "00_intermediate_data",
                            'record_rule_filter_peak_group.RData'))
    }
  }

  # 02 remove overlap peak group
  cat('Remove overlap peak group\n')
  list_peak_group_annotation_concised <- concisePeakGroup2PeakGroup(
    list_peak_group_annotation = list_peak_group_annotation_concised,
    path_dir = path_dir,
    is_save_intermidate = is_save_intermidate)

  temp <- statNumListPeakGroup(list_peak_group_annotation_concised)
  num_stat <- num_stat %>% dplyr::bind_rows(temp)

  temp <- statTypeListPeakGroup(list_peak_group_annotation_concised)
  type_stat <- type_stat %>% dplyr::bind_rows(temp)

  rownames(num_stat) <- c('initial annotation',
                          'overlap removal')

  rownames(type_stat) <- c('initial annotation',
                           'overlap removal')

  cat('General statistics:')
  print(knitr::kable(num_stat)); cat('\n\n')

  cat('Annotation type summary:')
  print(knitr::kable(type_stat)); cat('\n\n')

  if (is_save_intermidate) {
    dir.create(file.path(path_dir,
                         '00_tracer_result',
                         "00_intermediate_data"),
               showWarnings = FALSE,
               recursive = TRUE)

    save(list_peak_group_annotation_concised,
         file = file.path(path_dir,
                          '00_tracer_result',
                          "00_intermediate_data",
                          'list_peak_group_annotation_concised.RData'))
  }

  return(list_peak_group_annotation_concised)

}






################################################################################
# concisePeakGroup2PeakGroup ---------------------------------------------------

#' @title concisePeakGroup2PeakGroup
#' @author Zhiwei Zhou
#' @description concise overlap peak groups (Peak group - peak group, e.g. M175T752_[M-H]-, ISF of M175T752 (M132T752_[M+NH4-2H]-, M132T752_[M-H]-), Isotope of M175T752 (M177T752_[M-H]-))
#' @param list_peak_group_annotation
#' @export
#' @examples
#' load(system.file("tempdata", "list_peak_group_annotation_concised_p2pg_200805.RData", package="MetDNA2"))
#' test <- concisePeak2PeakGroup(list_peak_group_annotation = list_peak_group_annotation_concised)

# load('./inst/tempdata/list_peak_group_annotation_concised_p2pg_200805.RData')
# test <- concisePeakGroup2PeakGroup(list_peak_group_annotation = list_peak_group_annotation_concised)

setGeneric(name = 'concisePeakGroup2PeakGroup',
           def = function(
    list_peak_group_annotation,
    path_dir,
    is_save_intermidate = TRUE,
    ...
           ){

             stat_peak_groups <- purrr::map(seq_along(list_peak_group_annotation), function(i){
               tibble::tibble(name = list_peak_group_annotation[[i]]@base_peak_name,
                              mz = list_peak_group_annotation[[i]]@base_peak_mz,
                              rt = list_peak_group_annotation[[i]]@base_peak_rt,
                              int = list_peak_group_annotation[[i]]@base_peak_int,
                              # ccs = list_peak_group_annotation[[i]]@base_peak_ccs,
                              adduct = list_peak_group_annotation[[i]]@base_peak_adduct,
                              # seed = list_peak_group_annotation[[i]]@base_peak_seed,
                              num_peak_group = nrow(list_peak_group_annotation[[i]]@peak_list),
                              num_peak_group_anno = nrow(list_peak_group_annotation[[i]]@peak_list_annotated))

             }) %>%
               dplyr::bind_rows() %>%
               dplyr::arrange(dplyr::desc(num_peak_group_anno), dplyr::desc(int))


             list_unconcised_peak_group <- stat_peak_groups$name
             list_concised_peak_group <- vector('list',
                                                length = length(list_unconcised_peak_group))
             record_overlap_peak_group <- vector('list',
                                                 length = length(list_unconcised_peak_group))
             i <- 1
             while (length(list_unconcised_peak_group) > 0) {
               # cat(i, ' ')

               list_unconcised_peak_name <- list_unconcised_peak_group %>%
                 sapply(., function(x){
                   stringr::str_split(x, pattern = '_\\[|_\\M')[[1]][1]
                 }) %>% unname()

               # select largest peak group
               base_peak_group_name <- list_unconcised_peak_group[1]
               base_peak_name <- (base_peak_group_name %>% stringr::str_split(pattern = '_\\[|_\\M'))[[1]][1]
               base_peak_group_peak_list <- list_peak_group_annotation[[base_peak_group_name]]@peak_list_annotated
               base_peak_adduct <- list_peak_group_annotation[[base_peak_group_name]]@base_peak_adduct

               # if dimer dectected and moner existed,
               #   skip the base_peak_adduct
               if (stringr::str_detect(base_peak_adduct, '2M|3M')) {
                 temp <- stat_peak_groups %>%
                   dplyr::filter(paste(name, adduct, sep = '_') %in% list_unconcised_peak_group) %>%
                   # dplyr::filter(name %in% list_unconcised_peak_name) %>%
                   dplyr::filter(name %in% base_peak_group_peak_list$peak_name) %>%
                   dplyr::filter(name != base_peak_name)

                 if (nrow(temp) > 0) {
                   if (sum(!stringr::str_detect(temp$adduct, '2M|3M')) > 0) {
                     list_unconcised_peak_group <- c(list_unconcised_peak_group[-1], base_peak_group_name)
                     next()
                   }
                 }

               }

               overlap_peak_group_name <- stat_peak_groups %>%
                 dplyr::filter(name %in% list_unconcised_peak_group) %>%
                 # dplyr::filter(name %in% list_unconcised_peak_name) %>%
                 dplyr::filter(name %in% base_peak_group_peak_list$peak_name) %>%
                 dplyr::filter(name != base_peak_name) %>%
                 dplyr::pull(name)

               # if overlap peak group existed,
               #    remove overlap peak groups & add to list_concised_peak_group
               #  else
               #    add to list_concised_peak_group

               if (length(overlap_peak_group_name) > 0) {
                 # remove overlaped peak groups
                 idx <- match(overlap_peak_group_name, names(list_peak_group_annotation))
                 list_peak_group_annotation <- list_peak_group_annotation[-idx]

                 # modify concised and unconcised peak groups
                 list_concised_peak_group[[i]] <- c(overlap_peak_group_name, base_peak_group_name)
                 record_overlap_peak_group[[i]] <- paste(overlap_peak_group_name, base_peak_group_name, sep = '@')
               } else {
                 list_concised_peak_group[[i]] <- c(base_peak_group_name)
                 record_overlap_peak_group[[i]] <- NULL
               }

               idx <- list_concised_peak_group[[i]] %>%
                 match(., list_unconcised_peak_group)
               list_unconcised_peak_group <- list_unconcised_peak_group[-idx]

               i <- i + 1
             }

             if (is_save_intermidate) {
               save(record_overlap_peak_group,
                    file = file.path(path_dir,
                                     '00_tracer_result',
                                     "00_intermediate_data",
                                     'record_overlap_peak_group.RData'))
             }

             return(list_peak_group_annotation)
           })



# self_check_isf ---------------------------------------------------------------
# load('~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/hyuA/00_tracer_result/00_intermediate_data/list_peak_group_annotation_merge.RData')

self_check_isf <- function(list_peak_group_annotation_merge,
                           ssc_cutoff = 0.8,
                           ppc_cutoff = 0.8) {

  stat_peak_groups <- purrr::map(seq_along(list_peak_group_annotation_merge), function(i){
    result <- tibble::tibble(name = list_peak_group_annotation_merge[[i]]@base_peak_name,
                             mz = list_peak_group_annotation_merge[[i]]@base_peak_mz,
                             rt = list_peak_group_annotation_merge[[i]]@base_peak_rt,
                             # ccs = list_peak_group_annotation[[i]]@base_peak_ccs,
                             adduct = list_peak_group_annotation_merge[[i]]@base_peak_adduct,
                             # seed = list_peak_group_annotation[[i]]@base_peak_seed,
                             num_peak_group = nrow(list_peak_group_annotation_merge[[i]]@peak_list),
                             num_peak_group_anno = nrow(list_peak_group_annotation_merge[[i]]@peak_list_annotated))

    base_peak_int <- list_peak_group_annotation_merge[[i]]@peak_list %>%
      dplyr::filter(peak_name == list_peak_group_annotation_merge[[i]]@base_peak_name) %>%
      dplyr::pull(int)

    result <- result %>% dplyr::mutate(int = base_peak_int)

  }) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(dplyr::desc(num_peak_group_anno), dplyr::desc(int))




  checked_peak_groups <- vector()
  unchecked_peak_groups <- stat_peak_groups$name
  # list_peak_group_annotation_merge$`173.0668@3.511`@peak_list

  while (length(unchecked_peak_groups) > 0) {
    temp_peak_group_name <- unchecked_peak_groups[1]
    cat(temp_peak_group_name, ' ')
    temp_base_mz <- list_peak_group_annotation_merge[[temp_peak_group_name]]@base_peak_mz
    temp_peak_group <- list_peak_group_annotation_merge[[temp_peak_group_name]]

    peak_high_cor_isf <- temp_peak_group@peak_list %>%
      dplyr::filter(ssc >= ssc_cutoff & ppc >= ppc_cutoff) %>%
      dplyr::filter(mz < temp_base_mz) %>%
      dplyr::pull(peak_name)

    # peak_high_cor_isf condition
    # 1. it is not a existed annotation
    # 2. it need to be in the peak group list
    # 3. it can't merge peak_group that have been groupped
    peak_high_cor_isf <- peak_high_cor_isf[!(peak_high_cor_isf %in% temp_peak_group@peak_list_annotated$peak_name)]
    peak_high_cor_isf <- peak_high_cor_isf[!(peak_high_cor_isf %in% checked_peak_groups)]
    peak_high_cor_isf <- peak_high_cor_isf[peak_high_cor_isf %in% names(list_peak_group_annotation_merge)]

    if (length(peak_high_cor_isf) > 0) {
      # add ISF(HighCor) annotation
      temp_result <- temp_peak_group@peak_list %>%
        initializePeakList() %>%
        dplyr::filter(peak_name %in% peak_high_cor_isf) %>%
        dplyr::select(peak_name:ms2_score, query_peak) %>%
        dplyr::mutate(query_peak = temp_peak_group_name,
                      isotopeAnnotation = NA,
                      adductAnnotation = NA,
                      neutralLossAnnotation = NA,
                      isfAnnotation = 'ISF(HighCor)')


      list_peak_group_annotation_merge[[temp_peak_group_name]]@peak_list_annotated <-
        list_peak_group_annotation_merge[[temp_peak_group_name]]@peak_list_annotated %>%
        dplyr::bind_rows(temp_result)

      # remove the ISF(HighCor) peak group
      temp_idx <- which(names(list_peak_group_annotation_merge) == peak_high_cor_isf)
      list_peak_group_annotation_merge <- list_peak_group_annotation_merge[-temp_idx]
    }

    checked_peak_groups <- c(checked_peak_groups, temp_peak_group_name)
    unchecked_peak_groups <- unchecked_peak_groups[!(unchecked_peak_groups %in% c(checked_peak_groups, peak_high_cor_isf))]
  }


  return(list_peak_group_annotation_merge)
  # get_summary_peak_group(peak_group_list = list_peak_group_annotation_merge)

}




################################################################################
# statistic peak group ---------------------------------------------------------

#' @title statNumListPeakGroup
#' @author Zhiwei Zhou
#' @description statistics number of peak groups
#' @param list_pak_group
#' @export
#' @examples
#' load(system.file("tempdata", "list_peak_group_annotation_200805.RData", package="MetDNA2"))
#' statNumListPeakGroup(list_peak_group = list_peak_group_annotation)


# load('./inst/tempdata/list_peak_group_200805.RData')
# statNumListPeakGroup(list_peak_group = list_peak_group)

# load('./inst/tempdata/list_peak_group_annotation_200805.RData')
# statNumListPeakGroup(list_peak_group = list_peak_group_annotation)


setGeneric(name = 'statNumListPeakGroup',
           def = function(
    list_peak_group
           ){
             num_peak <- sapply(list_peak_group, function(x){
               x@base_peak_name
             }) %>%
               unname() %>%
               unique() %>%
               length()

             num_peak_group <- sapply(list_peak_group, function(x){
               paste0(x@base_peak_name, '_', x@base_peak_adduct)
             }) %>%
               unname() %>%
               unique() %>%
               length()

             num_annotated_peak <- sapply(list_peak_group, function(x){
               nrow(x@peak_list_annotated)
             }) %>%
               unname() %>%
               sum()

             result <- tibble::tibble('No. of peak' = num_peak,
                                      'No. of peak group' = num_peak_group,
                                      'No. of annotated peak' = num_annotated_peak)

             return(result)

           })


#' @title statTypeListPeakGroup
#' @author Zhiwei Zhou
#' @description
#' @param list_peak_group
#' @export
#' @examples

# load('./inst/tempdata/list_peak_group_annotation_200805.RData')
# statTypeListPeakGroup(list_peak_group = list_peak_group_annotation)

setGeneric(name = 'statTypeListPeakGroup',
           function(
    list_peak_group
           ){
             table_annotated <- purrr::map(list_peak_group, function(x){
               x@peak_list_annotated %>%
                 dplyr::mutate(peak_group = paste0(x@base_peak_name,
                                                   '_',
                                                   x@base_peak_adduct))
             }) %>% dplyr::bind_rows()

             if (nrow(table_annotated) == 0) {
               result <- tibble::tibble('No. of unique peak' = 0,
                                        'No. of peak' = 0,
                                        'No. of isotope' = 0,
                                        'No. of adduct' = 0,
                                        'No. of NL' = 0,
                                        'No. of ISF' = 0)

               return(result)
             }

             table_annotated <- table_annotated %>%
               tidyr::pivot_longer(cols = c('adductAnnotation',
                                            'neutralLossAnnotation',
                                            'isfAnnotation'),
                                   names_to = 'adduct_nl_isf',
                                   values_to = 'label') %>%
               dplyr::filter(!is.na(label) | isotopeAnnotation %in% c('[M+1]',
                                                                      '[M+2]',
                                                                      '[M+3]',
                                                                      '[M+4]'))

             temp_isotope <- table_annotated %>%
               dplyr::filter(isotopeAnnotation %in% c('[M+1]',
                                                      '[M+2]',
                                                      '[M+3]',
                                                      '[M+4]')) %>%
               dplyr::distinct(peak_name,
                               isotopeAnnotation,
                               peak_group,
                               .keep_all = TRUE)

             temp_adduct_nl_isf <- table_annotated %>%
               dplyr::filter(!is.na(label))


             num_unique_peak <- table_annotated %>%
               dplyr::distinct(peak_name) %>%
               dplyr::count() %>%
               dplyr::pull(n)

             num_peak <- nrow(temp_isotope) + nrow(temp_adduct_nl_isf)
             num_isotope <- nrow(temp_isotope)
             num_adduct <- temp_adduct_nl_isf %>%
               dplyr::filter(adduct_nl_isf == 'adductAnnotation') %>% nrow()
             num_nl <- temp_adduct_nl_isf %>%
               dplyr::filter(adduct_nl_isf == 'neutralLossAnnotation') %>% nrow()
             num_isf <- temp_adduct_nl_isf %>%
               dplyr::filter(adduct_nl_isf == 'isfAnnotation') %>% nrow()

             result <- tibble::tibble('No. of unique peak' = num_unique_peak,
                                      'No. of peak' = num_peak,
                                      'No. of isotope' = num_isotope,
                                      'No. of adduct' = num_adduct,
                                      'No. of NL' = num_nl,
                                      'No. of ISF' = num_isf)

             return(result)
           })


#' @title convertSemiTargetedResult2LongTable
#' @author Zhiwei Zhou
#' @description convert wide semi targeted clustering result to long table
#' @param result_semi_table
#' @export
#' @examples
#' result_semi_table <- readr::read_csv('I:/00_projects/03_MetDNA2/00_data/20200805_targeted_annotation_evaluation_200STD/01_targeted_annotation/200STD_neg_reverse/annot_credential/semi_targeted_annotation_result.csv')
#' test <- convertSemiTargetedResult2LongTable(result_semi_table)

# result_semi_table <- readr::read_csv('I:/00_projects/03_MetDNA2/00_data/20200805_targeted_annotation_evaluation_200STD/01_targeted_annotation/200STD_neg_reverse/annot_credential/semi_targeted_annotation_result.csv')
#
# test <- convertSemiTargetedResult2LongTable(result_semi_table)

setGeneric(name = 'convertSemiTargetedResult2LongTable',
           def = function(
    result_semi_table
           ){
             result <- result_semi_table %>%
               tidyr::separate_rows(ssc,
                                    int,
                                    mz_error,
                                    rt_error,
                                    # ccs_error,
                                    int_ratio,
                                    ms2_score,
                                    isotope,
                                    type_annotation,
                                    annotation,
                                    peak_group,
                                    peak_group_scale,
                                    sep = ';') %>%
               readr::type_convert()

             return(result)
           })


setGeneric(name = 'get_summary_peak_group',
           def = function(peak_group_list){
             stat_peak_groups <- purrr::map(seq_along(peak_group_list), function(i){
               result <- tibble::tibble(name = peak_group_list[[i]]@base_peak_name,
                                        mz = peak_group_list[[i]]@base_peak_mz,
                                        rt = peak_group_list[[i]]@base_peak_rt,
                                        # ccs = list_peak_group_annotation[[i]]@base_peak_ccs,
                                        adduct = peak_group_list[[i]]@base_peak_adduct,
                                        # seed = list_peak_group_annotation[[i]]@base_peak_seed,
                                        num_peak_group = nrow(peak_group_list[[i]]@peak_list),
                                        num_peak_group_anno = nrow(peak_group_list[[i]]@peak_list_annotated))

               base_peak_int <- peak_group_list[[i]]@peak_list %>%
                 dplyr::filter(peak_name == peak_group_list[[i]]@base_peak_name) %>%
                 dplyr::pull(int)

               result <- result %>% dplyr::mutate(int = base_peak_int)

             }) %>%
               dplyr::bind_rows() %>%
               dplyr::arrange(dplyr::desc(num_peak_group_anno), dplyr::desc(int))


             return(stat_peak_groups)
           })

################################################################################

generate_recoginize_table <- function(list_peak_group_annotation_merge) {
  peak_group_summary <- get_summary_peak_group(list_peak_group_annotation_merge) %>%
    dplyr::mutate(group_order = paste0('[', seq(dplyr::n()), ']'))

  temp_result_table <- lapply(list_peak_group_annotation_merge, function(x){
    x@peak_list_annotated %>%
      tidyr::pivot_longer(-c('peak_name':'isotopeAnnotation'), names_to = 'type', values_to = 'relationship') %>%
      dplyr::select(-c('isotopeAnnotation', 'type', 'int_ratio', 'ms2_score')) %>%
      dplyr::filter(!is.na(relationship)) %>%
      dplyr::rename('base_peak' = query_peak)
  }) %>% dplyr::bind_rows()

  idx <- match(temp_result_table$base_peak, peak_group_summary$name)
  result <- temp_result_table %>%
    dplyr::mutate(num_annotated_peak = peak_group_summary$num_peak_group_anno[idx],
                  group_order = peak_group_summary$group_order[idx]) %>%
    dplyr::arrange(match(base_peak, peak_group_summary$name), mz)

  return(result)
}


################################################################################
# Plots in AnnotationCredential ------------------------------------------------
#   plotPseudoMs1Spec ----------------------------------------------------------
#' @title plotPseudoMs1Spec
#' @description generate Pseudo MS1 spectrum for annotated peak group
#' @author Zhiwei Zhou
#' @param peak_group
#' @export
#' @examples
#' load(system.file("tempdata", "list_peak_group_annotation_200805.RData", package="MetDNA2"))
#' list_peak_group_annotation[[35]] %>% plotPseudoMs1Spec()

# load('./inst/tempdata/list_peak_group_annotation_200805.RData')
# list_peak_group_annotation[[35]] %>% plotPseudoMs1Spec()

setGeneric(name = 'plotPseudoMs1Spec',
           def = function(
    peak_group
           ){
             if (nrow(peak_group@peak_list_annotated) == 0) {
               stop('No annotated MS1 peaks\n')
             }
             temp_data <- peak_group@peak_list_annotated %>%
               tidyr::pivot_longer(cols = c('isotopeAnnotation', 'adductAnnotation',
                                            'neutralLossAnnotation', 'isfAnnotation'),
                                   names_to = 'type',
                                   values_to = 'label') %>%
               dplyr::filter(!is.na(label)) %>%
               dplyr::arrange(peak_name,
                              match(type, c('adductAnnotation',
                                            'neutralLossAnnotation',
                                            'isfAnnotation',
                                            'isotopeAnnotation'))) %>%
               dplyr::distinct(peak_name, .keep_all = TRUE) %>%
               dplyr::mutate(int = int/max(int))


             base_peak_mz <- peak_group@base_peak_mz
             base_peak_int <- temp_data %>%
               dplyr::filter(peak_name == peak_group@base_peak_name) %>%
               dplyr::pull(int)
             base_peak_adduct <- peak_group@base_peak_adduct


             temp_plot <- ggplot2::ggplot(temp_data) +
               ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz,
                                                  y = 0, yend = int,
                                                  colour = type)) +
               ggplot2::geom_point(ggplot2::aes(x = base_peak_mz, y = base_peak_int))+
               ggplot2::scale_colour_manual(values = c(
                 'adductAnnotation' = 'tomato',
                 'neutralLossAnnotation' = 'dodgerblue',
                 'isfAnnotation' = 'orange',
                 'isotopeAnnotation' = 'limegreen'
               )) +
               ggplot2::geom_hline(yintercept = 0) +
               ggplot2::xlim(0.95*min(temp_data$mz),
                             1.05*max(temp_data$mz)) +
               ggplot2::ylim(0, 1.05) +
               ggplot2::xlab('m/z') +
               ggplot2::ylab('relative intensity') +
               ZZWTheme() +
               ggplot2::theme(legend.position = c(0.8, 0.8))

             # add labels in the spectrum
             temp_adduct <- temp_data %>%
               dplyr::filter(type %in% c('adductAnnotation'))

             if (nrow(temp_adduct)>0) {
               temp_plot <- temp_plot +
                 ggplot2::annotate('text',
                                   x = temp_adduct$mz,
                                   y = temp_adduct$int + 0.04,
                                   label = temp_adduct$label)
             }

             temp_nl <- temp_data %>%
               dplyr::filter(type %in% c('neutralLossAnnotation'))

             if (nrow(temp_nl)>0) {

               # if protonation/deprotpnation base peak
               #    plot loss formula
               #   else only add label
               if (base_peak_adduct %in% c('[M+H]+', '[M-H]-')) {
                 temp_label <- temp_nl$label %>%
                   stringr::str_extract(pattern = 'M\\-(\\w+)') %>%
                   stringr::str_replace(pattern = 'M', replacement = '')

                 temp_nl <- temp_nl %>%
                   dplyr::select(peak_name:int, type, label) %>%
                   dplyr::mutate(xend = base_peak_mz)

                 temp_plot <- temp_plot +
                   ggplot2::geom_segment(ggplot2::aes(x = mz,
                                                      xend = xend,
                                                      y = 0.8*int,
                                                      yend = 0.8*int),
                                         linetype = 'dashed',
                                         data = temp_nl) +
                   ggplot2::annotate('text',
                                     x = (base_peak_mz - temp_nl$mz)/2 + temp_nl$mz,
                                     y = temp_nl$int*0.8 + 0.04,
                                     label = temp_label)

                 # temp_plot <- temp_plot +
                 #   ggplot2::geom_segment(aes(x = temp_nl$mz,
                 #                             xend = rep(base_peak_mz, nrow(temp_nl)),
                 #                             y = 0.8*temp_nl$int,
                 #                             yend = 0.8*temp_nl$int),
                 #                         linetype = 'dashed') +
                 #   ggplot2::annotate('text',
                 #                     x = (base_peak_mz - temp_nl$mz)/2 + temp_nl$mz,
                 #                     y = temp_nl$int*0.8 + 0.04,
                 #                     label = temp_label)
               } else {
                 temp_plot <- temp_plot +
                   ggplot2::annotate('text',
                                     x = temp_nl$mz,
                                     y = temp_nl$int + 0.04,
                                     label = temp_nl$label)
               }

             }

             temp_isf <- temp_data %>%
               dplyr::filter(type %in% c('isfAnnotation'))

             if (nrow(temp_adduct)>0) {
               temp_plot <- temp_plot +
                 ggplot2::annotate('text',
                                   x = temp_isf$mz,
                                   y = temp_isf$int + 0.04,
                                   label = temp_isf$label)

             }


             return(temp_plot)

           })






################################################################################
# ZZWTheme ---------------------------------------------------------------------
#' @title ZZWTheme
#' @author Zhiwei Zhou
#' \email{zhouzw@@stanford.edu}
#' @description ggplot theme
#' @param type: 'common', 'classic'
#' @export
#' @examples
#' mydata <- data.frame(
#'  Lebal  = c("linerange1","linerange2","linerange3","linerange4","linerange5"),
#'  xstart = c(3.5,7,12,16,20),
#'  ymin   = c(2.5,6.5,3,4.5,3.8),
#'  ymax   = c(7.5,9.5,9,13.5,4.2),
#'  class  = c("A","A","A","C","C")
#' )
#' test <- ggplot(mydata) +
#'  geom_linerange(aes(x = xstart, ymin = ymin , ymax = ymax , colour = class) , size = 1.5)
#' test
#' mytheme <- ZZWTheme()
#' test + mytheme
#' mytheme2 <- ZZWTheme('classic')
#' test +
#'   mytheme2 +
#'   theme(axis.text.y = ggplot2::element_text(size = 10, angle = 0,
#'                                             colour = 'black'))



# mydata <- data.frame(
#   Lebal  = c("linerange1","linerange2","linerange3","linerange4","linerange5"),
#   xstart = c(3.5,7,12,16,20),
#   ymin   = c(2.5,6.5,3,4.5,3.8),   ymax   = c(7.5,9.5,9,13.5,4.2),
#   class  = c("A","A","A","C","C")
# )
# test <- ggplot(mydata) +
#   geom_linerange(aes(x = xstart, ymin = ymin , ymax = ymax , colour = class) , size = 1.5)
# test
# mytheme <- ZZWTheme()
# test + mytheme
#
# mydata <- data.frame(
#   Lebal  = c("linerange1","linerange2","linerange3","linerange4","linerange5"),
#   xstart = c(3.5,7,12,16,20),
#   ymin   = c(2.5,6.5,3,4.5,3.8),   ymax   = c(7.5,9.5,9,13.5,4.2),
#   class  = c("A","A","A","C","C")
# )
# test <- ggplot(mydata) +
#   geom_linerange(aes(x = xstart, ymin = ymin , ymax = ymax , colour = class) , size = 1.5)
# test
# mytheme1 <- ZZWTheme()
# test + mytheme1
# mytheme2 <- ZZWTheme('classic')
# test +
#   mytheme2 +
#   theme(axis.text.y = ggplot2::element_text(size = 10, angle = 0,
#                                             colour = 'black'))

setGeneric(name = 'ZZWTheme',
           def = function(
    type = c('common', 'classic')
           ){

             type <- match.arg(type)
             switch (type,
                     'common' = {
                       result <- ggplot2::theme_bw() +
                         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.4,
                                                                           face = "bold",
                                                                           size = 14),
                                        panel.grid.major = ggplot2::element_blank(),
                                        panel.grid.minor = ggplot2::element_blank(),
                                        panel.background = ggplot2::element_blank(),
                                        panel.border = ggplot2::element_rect(fill = NA, colour = 'black'),
                                        legend.background = ggplot2::element_blank(),
                                        axis.text.x = ggplot2::element_text(size = 10,
                                                                            colour = 'black',
                                                                            hjust = 0.5,
                                                                            vjust = 0.5),
                                        axis.text.y = ggplot2::element_text(size = 10,
                                                                            angle = 90,
                                                                            colour = 'black',
                                                                            hjust = 0.5,
                                                                            vjust = 0.5),
                                        axis.title.x = ggplot2::element_text(size = 12, colour = 'black'),
                                        axis.title.y = ggplot2::element_text(size = 12, colour = 'black'),
                                        axis.ticks = ggplot2::element_line(colour = 'black'),
                                        axis.ticks.length = ggplot2::unit(2, "mm"))
                     },
                     'classic' = {
                       result <- ggplot2::theme_classic() +
                         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.4,
                                                                           face = "bold",
                                                                           size = 14),
                                        panel.grid.major = ggplot2::element_blank(),
                                        panel.grid.minor = ggplot2::element_blank(),
                                        panel.background = ggplot2::element_rect(fill = NA, colour = NA),
                                        # panel.border = ggplot2::element_rect(fill = NA, colour = 'black'),
                                        legend.background = ggplot2::element_blank(),
                                        axis.line = ggplot2::element_line(colour = 'black'),
                                        axis.text.x = ggplot2::element_text(size = 10,
                                                                            colour = 'black',
                                                                            hjust = 0.5,
                                                                            vjust = 0.5),
                                        axis.text.y = ggplot2::element_text(size = 10,
                                                                            angle = 90,
                                                                            colour = 'black',
                                                                            hjust = 0.5,
                                                                            vjust = 0.5),
                                        axis.title.x = ggplot2::element_text(size = 12, colour = 'black'),
                                        axis.title.y = ggplot2::element_text(size = 12, colour = 'black'),
                                        axis.ticks = ggplot2::element_line(colour = 'black'),
                                        axis.ticks.length = ggplot2::unit(2, "mm"))
                     }
             )



             return(result)
           }
)

