################################################################################
# plot_diff_feature ------------------------------------------------------------

plot_diff_feature <- function(raw_table,
                              p_value_cutoff,
                              fold_change_cutoff,
                              p_adjust = TRUE,
                              text_title) {
  # browser()
  if (p_adjust) {
    temp_plot1 <- ggplot2::ggplot(raw_table) +
      ggplot2::geom_point(ggplot2::aes(x = log2_fc, y = -log10(q_values),
                                       colour = increase_label)) +
      ggrepel::geom_text_repel(
        data = subset(raw_table, increase_label == 'signif_increase'),
        ggplot2::aes(x = log2_fc, y = -log10(q_values), label = id),
        # position = position_jitter(seed = 1),
        box.padding = 1
      ) +
      ggplot2::geom_vline(xintercept = log2(fold_change_cutoff), linetype = "dashed") +
      ggplot2::geom_hline(yintercept = -log10(p_value_cutoff), linetype = "dashed") +
      ZZWTheme()
  } else {
      temp_plot1 <- ggplot2::ggplot(raw_table) +
        ggplot2::geom_point(ggplot2::aes(x = log2_fc, y = -log10(p_values),
                                         colour = increase_label)) +
        ggrepel::geom_text_repel(
          data = subset(raw_table, increase_label == 'signif_increase'),
          ggplot2::aes(x = log2_fc, y = -log10(p_values), label = id),
          # position = position_jitter(seed = 1),
          box.padding = 1
        ) +
        ggplot2::geom_vline(xintercept = log2(fold_change_cutoff), linetype = "dashed") +
        ggplot2::geom_hline(yintercept = -log10(p_value_cutoff), linetype = "dashed") +
        ZZWTheme()
  }

  if (!missing(text_title)) {
    temp_plot1 <- temp_plot1 +
      ggplot2::ggtitle(text_title)
  }

  return(temp_plot1)

}





################################################################################
# plot_isotope_pairs -----------------------------------------------------------

plot_isotope_pairs <- function(signif_feature_unlabeled,
                               signif_feature_labeled,
                               pair_table,
                               eic_data_unlabel,
                               eic_data_label,
                               col_abundance = 'avg_abundance_case',
                               rt_unit = c('minutes', 'seconds'),
                               mz_tol = 10) {

  match.arg(rt_unit)
  # browser()

  # point table
  # temp_labeled <- signif_feature_unlabeled %>%
  #   dplyr::select(id, mz, rt, get(average_abundance), q_value, fold_change) %>%
  #   dplyr::mutate(group_label = 'Unlabeled',
  #                 delta_mz = 0.1)
  #
  # temp_unlabeled <- signif_feature_labeled %>%
  #   dplyr::select(id, mz, rt, average_abundance, q_value, fold_change) %>%
  #   dplyr::mutate(group_label = 'Labeled',
  #                 delta_mz = -0.1)

  temp_labeled <- signif_feature_unlabeled %>%
    dplyr::select(id, mz, rt, col_abundance, q_values, fold_change) %>%
    # dplyr::rename('average_abundance' = col_abundance) %>%
    dplyr::mutate(group_label = 'Unlabeled',
                  delta_mz = 0.1)

  temp_unlabeled <- signif_feature_labeled %>%
    dplyr::select(id, mz, rt, col_abundance, q_values, fold_change) %>%
    dplyr::mutate(group_label = 'Labeled',
                  delta_mz = -0.1)

  point_table <- temp_unlabeled %>%
    dplyr::bind_rows(temp_labeled)

  # pair table
  pair_table <- pair_table %>%
    dplyr::mutate(avg_mz = (unlabeled_mz + labeled_mz)/2) %>%
    dplyr::mutate(mz_adj1 = unlabeled_mz - avg_mz,
                  mz_adj2 = labeled_mz - avg_mz)

  # replace delta_mz in point table
  temp_idx1 <- match(pair_table$unlabeled_feature_id, point_table$id)
  temp_idx2 <- match(pair_table$labelded_feature_id, point_table$id)

  point_table$delta_mz[temp_idx1] <- pair_table$mz_adj1
  point_table$delta_mz[temp_idx2] <- pair_table$mz_adj2

  # segment table1
  segment_table <- pair_table %>%
    dplyr::mutate(rt_mean = (unlabeled_rt + labeled_rt)/2) %>%
    dplyr::mutate(rt_adj = rt_mean - 0.1) %>%
    dplyr::select(mz_adj1, mz_adj2, rt_mean, rt_adj, mass_shift_label)

  # rescale EICs
  temp_scale <- round(pair_table$labeled_mz - pair_table$unlabeled_mz) %>% max()
  eic_data_unlabel <- eic_data_unlabel %>%
    dplyr::mutate(int_adj = int/max(int)*(-temp_scale/2) - 0.1)
  eic_data_label <- eic_data_label %>%
    dplyr::mutate(int_adj = int/max(int)*(temp_scale/2) + 0.1)

  if (rt_unit == 'seconds') {
    eic_data_unlabel$rt <- eic_data_unlabel$rt*60
    eic_data_label$rt <- eic_data_label$rt*60
    xlab_text <- 'RT (seconds)'
  } else {
    eic_data_unlabel$rt <- eic_data_unlabel$rt
    eic_data_label$rt <- eic_data_label$rt
    xlab_text <- 'RT (minutes)'
  }


  temp_plot <- ggplot2::ggplot(point_table) +
    ggplot2:: geom_line(data = eic_data_label,
                        ggplot2::aes(x = rt, y = int_adj, color = precursor_mz),
                        alpha = 0.6) +
    ggplot2::geom_line(data = eic_data_unlabel,
                       ggplot2::aes(x = rt, y = int_adj, color = precursor_mz), alpha = 0.6) +
    ggplot2::geom_segment(data = segment_table,
                          ggplot2::aes(x = rt_adj,
                                       y = mz_adj1,
                                       xend = rt_adj,
                                       yend = mz_adj2),
                          color = 'gray') +
    ggplot2::geom_segment(data = segment_table,
                          ggplot2::aes(x = rt_adj,
                                       y = mz_adj1,
                                       xend = rt_mean,
                                       yend = mz_adj1),
                          color = 'gray') +
    ggplot2::geom_segment(data = segment_table,
                          ggplot2::aes(x = rt_adj,
                                       y = mz_adj2,
                                       xend = rt_mean,
                                       yend = mz_adj2),
                          color = 'gray') +
    ggplot2::geom_point(ggplot2::aes(x = rt,
                                     y = delta_mz,
                                     fill = group_label,
                                     size = fold_change),
                        alpha = 0.5, shape = 21) +
    ggrepel::geom_text_repel(
      data = subset(point_table, abs(delta_mz) > 0.5 ),
      ggplot2::aes(x = rt, y = delta_mz, label = id),
      position = position_jitter(seed = 1)
    ) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::scale_fill_manual(values = c('Labeled' = 'red',
                                          'Unlabeled' = 'black'),
                               label = c('Labeled' = '13C',
                                         'Unlabeled' = '12C'),
                               name = 'Group') +
    ggplot2::scale_colour_discrete(name = 'Precursor m/z') +
    ggplot2::scale_size_continuous(range = c(1, 10),
                                   name = paste0('Fold_change (Mutant/WT)')) +
    ggplot2::scale_y_continuous(breaks = c((seq(-floor(temp_scale/2), -1, by = 1)) - 0.1,
                                           0,
                                           (seq(1, floor(temp_scale/2), by = 1)) + 0.1),
                                labels = seq(-floor(temp_scale/2), floor(temp_scale/2), by = 1)) +
    ggplot2::ylab('Relative Intensity') +
    ggplot2::xlab(xlab_text) +
    ZZWTheme() +
    ggplot2::theme(legend.position = 'right')

  return(temp_plot)

}


################################################################################
# extract_eic_data -------------------------------------------------------------


# path <- '.'
# files_pattern <- 'ygeX_UA.+\\.mzML'
# mz_list <- pair_table$unlabeled_mz
# mz_tol = 10

extract_eic_data <- function(path,
                             files = NULL,
                             files_pattern = NULL, # 'ygeX_UA.+\\.mzML'
                             mz_list,
                             mz_tol = 10) {
  # browser()
  require(data.table)
  require(RaMS)

  if (is.null(files)) {
    if (is.null(files_pattern)) {
      stop('Please provide files or files_pattern\n')
    }

    files <- list.files(path = path, pattern = files_pattern, recursive = TRUE)

    if (length(files) == 0) {
      stop('No effective mzML data found\n')
    }
  }

  msdata <- RaMS::grabMSdata(file.path(path, files), grab_what = 'MS1')

  # unlabeled EIC
  mz_list <- unique(round(mz_list, 4))

  temp_mz_tol <- lapply(mz_list, function(x){
    RaMS::pmppm(x, ppm = mz_tol)
  })

  idx_list <- lapply(temp_mz_tol, function(x){
    which(msdata$MS1$mz >= x[1] & msdata$MS1$mz <= x[2])
  })

  eic_data_export <- mapply(function(x, y){
    msdata$MS1[x,] %>%
      dplyr::mutate(label = stringr::str_replace(filename, '\\_\\d+\\.mzML', ''),
                    precursor_mz = as.character(y))
  },
  x = idx_list,
  y = mz_list,
  SIMPLIFY = FALSE) %>%
    dplyr::bind_rows()

  return(eic_data_export)
}


