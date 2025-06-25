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




################################################################################
# plot_selected_pair -----------------------------------------------------------

plot_selected_pair <- function(signif_feature_unlabeled,
                               signif_feature_labeled,
                               pair_table,
                               eic_data_unlabel,
                               eic_data_label,
                               col_abundance = 'avg_abundance_case',
                               rt_unit = c('minutes', 'seconds'),
                               mz_tol = 10) {

  match.arg(rt_unit)

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
  max_int_unlabel <- eic_data_unlabel %>%
    dplyr::pull(int) %>%
    max()
  max_int_label <- eic_data_label %>%
    dplyr::pull(int) %>%
    max()
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


  # return the list pair table
  list_plots <- lapply(seq_along(pair_table$unlabeled_feature_id), function(i){
    temp_pair_table <- pair_table %>% dplyr::slice(i)
    temp_feature_unlabel <- pair_table$unlabeled_feature_id[i]
    temp_feature_label <- pair_table$labelded_feature_id[i]

    temp_point_table_unlabel <- point_table %>%
      dplyr::filter(id == temp_feature_unlabel & group_label == 'Unlabeled')

    temp_point_table_label <- point_table %>%
      dplyr::filter(id == temp_feature_label & group_label == 'Labeled')

    temp_point_table <- temp_point_table_unlabel %>%
      dplyr::bind_rows(temp_point_table_label)

    temp_eic_data_unlabel <- eic_data_unlabel %>%
      dplyr::filter(precursor_mz == pair_table$unlabeled_mz[i])

    temp_eic_data_label <- eic_data_label %>%
      dplyr::filter(precursor_mz == pair_table$labeled_mz[i])

    temp_segment_table <- segment_table %>%
      dplyr::slice(i)


    temp_max_int_unlabel <- temp_eic_data_unlabel %>%
      dplyr::filter(abs(rt - temp_pair_table$unlabeled_rt) <= 0.1) %>%
      dplyr::pull(int_adj) %>%
      min()

    temp_max_int_label <- temp_eic_data_label %>%
      dplyr::filter(abs(rt - temp_pair_table$labeled_rt) <= 0.1) %>%
      dplyr::pull(int_adj) %>%
      max()

    temp_point_table <- temp_point_table %>%
      dplyr::mutate(rela_int = c(temp_max_int_unlabel, temp_max_int_label))

    text_unlabel <- paste0('Unlabeled feature: ', temp_feature_unlabel, '\n',
                          'm/z: ', temp_pair_table$unlabeled_mz, '\n',
                          'RT: ', temp_pair_table$unlabeled_rt, '\n',
                          "Abundance: ", format(signif(temp_pair_table$average_abundance_unlabel, 3), scientific = TRUE), '\n',
                          'Fold change: ', round(temp_pair_table$fold_change_unlabel), '\n',
                          'q-value: ', format(signif(temp_pair_table$q_value_unlabel, 2), scientific = TRUE))

    text_label <- paste0('Labeled feature: ', temp_feature_label, '\n',
                         'm/z: ', temp_pair_table$labeled_mz, '\n',
                         'RT: ', temp_pair_table$labeled_rt, '\n',
                         "Abundance: ", format(signif(temp_pair_table$average_abundance_label, 3), scientific = TRUE), '\n',
                         'Fold change: ', round(temp_pair_table$fold_change_label), '\n',
                         'q-value: ', format(signif(temp_pair_table$q_value_label, 2), scientific = TRUE))

    text_title <- paste0('Isotope pair: ', temp_feature_unlabel, ' - ', temp_feature_label, '\n',
                         '#C ', temp_pair_table$mass_shift_label[i], '\n')

    temp_plot <- ggplot2::ggplot(temp_point_table) +
      ggplot2::geom_line(data = temp_eic_data_label,
                          ggplot2::aes(x = rt, y = int_adj),
                          size = 1,
                          alpha = 0.6,
                          color = '#fdb562') +
      ggplot2::geom_line(data = temp_eic_data_unlabel,
                         ggplot2::aes(x = rt, y = int_adj, color = precursor_mz),
                         alpha = 0.6,
                         size = 1,
                         color = "#8cd4c8") +
      ggplot2::geom_point(ggplot2::aes(x = rt,
                                       y = rela_int,
                                       fill = group_label),
                          size = 4,
                          alpha = 0.5,
                          shape = 21) +
      ZZW_annotate_text2(label = text_label, x = 0.9, y = 0.9) +
      ZZW_annotate_text2(label = text_unlabel, x = 0.9, y = 0.1) +
      ggplot2::scale_y_continuous(breaks = c((seq(-(temp_scale/2), -1, by = 1)) - 0.1,
                                             0,
                                             (seq(1, (temp_scale/2), by = 1)) + 0.1),
                                  labels = c((seq(-(temp_scale/2), -1, by = 1)) - 0.1,
                                             0,
                                             (seq(1, (temp_scale/2), by = 1)) + 0.1),
                                  limits = c(-(temp_scale/2) - 0.1,
                                             (temp_scale/2) + 0.1)) +
      ggplot2::ylab('Relative Intensity') +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::xlab(xlab_text) +
      ggtitle(text_title) +
      ZZWTheme() +
      ggplot2::theme(legend.position = 'right')


    return(temp_plot)
  })

  return(list_plots)

}


################################################################################
#' @title ZZW_annotate_text2
#' @author Zhiwei Zhou
#' @param label
#' @param x
#' @param y
#' @export
#' @examples
#' library(ggplot2)
#' ggplot(diamonds, aes(carat, price, color = cut)) +
#'   geom_point() +
#'   annotate(geom = 'text', x = 4, y = 5000,
#'            label = "a\nb\nc")
#' qplot(1:10,1:10) + ZZW_annotate_text2('some long text\nx = 1', x=1, y=0, hjust=1)

# library(ggplot2)
#
# ggplot2::ggplot(ggplot2::diamonds, ggplot2::aes(carat, price, color = cut)) +
#   ggplot2::geom_point() +
#   ggplot2::annotate(geom = 'text', x = 4, y = 5000,
#            label = "atop(a,b,c)",
#            parse = TRUE, color = 'blue') +
# ZZW_annotate_text2('some long text\nx = 1', x=1, y=0, hjust=1)
#
# ggplot2::ggplot(ggplot2::diamonds, ggplot2::aes(carat, price, color = cut)) +
#   ggplot2::geom_point() +
#   ggplot2::annotate(geom = 'text', x = 4, y = 5000,
#                     label = "atop(a,b,c)",
#                     parse = TRUE, color = 'blue') +
#   ZZW_annotate_text2('some long text\nx = 1', x=0, y=1, hjust=1)
#
#
# ggplot(diamonds, aes(carat, price, color = cut)) +
#   geom_point() +
#   annotate(geom = 'text', x = 4, y = 5000,
#            label = "a\nb\nc")
#
#
# library("ggplot2")
ZZW_annotate_text2 <- function(label, x, y, facets=NULL, hjust=0, vjust=0, color='black', alpha=NA,
                               family=thm$text$family, size=thm$text$size, fontface=1, lineheight=1.0,
                               box_just=ifelse(c(x,y)<0.5,0,1), margin=grid::unit(size/2, 'pt'), thm=ggplot2::theme_get()) {
  x <- scales::squish_infinite(x)
  y <- scales::squish_infinite(y)
  data <- if (is.null(facets)) data.frame(x=NA) else data.frame(x=NA, facets)

  tg <- grid::textGrob(
    label, x=0, y=0, hjust=hjust, vjust=vjust,
    gp=grid::gpar(col=ggplot2::alpha(color, alpha), fontsize=size, fontfamily=family, fontface=fontface, lineheight=lineheight)
  )
  ts <- grid::unit.c(grid::grobWidth(tg), grid::grobHeight(tg))
  vp <- grid::viewport(x=x, y=y, width=ts[1], height=ts[2], just=box_just)
  tg <- grid::editGrob(tg, x=ts[1]*hjust, y=ts[2]*vjust, vp=vp)
  inner <- grid::grobTree(tg, vp=grid::viewport(width=grid::unit(1, 'npc')-margin*2, height=grid::unit(1, 'npc')-margin*2))

  ggplot2::layer(
    data = NULL,
    stat = ggplot2::StatIdentity,
    position = ggplot2::PositionIdentity,
    geom = ggplot2::GeomCustomAnn,
    inherit.aes = TRUE,
    params = list(
      grob=grid::grobTree(inner),
      xmin=-Inf,
      xmax=Inf,
      ymin=-Inf,
      ymax=Inf
    )
  )
}


