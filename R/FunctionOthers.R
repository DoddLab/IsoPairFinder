################################################################################
# some functions developed by other packages/tools -----------------------------
# Please refer to the original packages/tools for more details
# Please cite the original packages/tools in your work

# DoddLabMetID:
# Link: https://github.com/DoddLab/DoddLabMetID
#   read_ms2
#   readMGF, ListMGF
#   readMSP
#   readCEF, ListCEF
#   readMzXML, readMzML, readMzML2
#   integrate_ms2
#   purifyMs2
#   combine_ms1_ms2
#   split_ms2
#   add_has_ms2_to_SpecAnnotationClass
#   add_ms2_purifty_to_SpecAnnotationClass


# masstools:
# Link: https://github.com/tidymass/masstools
#   mz_rt_match
#   keep_one


################################################################################
# read_ms2 ----------------------------------------------------------------------

#' @title read_ms2
#' @author Zhiwei Zhou
#' @param ms2_file MS/MS file names
#' @param ms2_type 'mgf', 'cef', 'msp', 'mzXML', 'mzML'
#' @export

#
# path <- '~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData_DoddLabMetabolomics/ms2/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mzML') %>% dir(path)[.]
# ms2_data <- read_ms2(ms2_file = file.path(path, temp_files),
#                      ms2_type = 'mzML')
#
# path <- '~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData_mgf/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# ms2_data <- read_ms2(ms2_file = file.path(path, temp_files),
#                      ms2_type = 'mgf')
#
setGeneric(name = 'read_ms2',
           def = function(
    ms2_file,
    ms2_type = c('mgf', 'cef', 'msp', 'mzXML', 'mzML'),
    ...
           ){
             ms2_type <- match.arg(ms2_type)

             cat('Read MS2...\n\n')

             switch(ms2_type,
                    'mgf' = {
                      pbapply::pboptions(type='timer', char='+')
                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readMGF(x, ...)
                      })
                    },

                    'cef' = {
                      pbapply::pboptions(type='timer', char='+')
                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readCEF(x, ...)
                      })
                    },

                    'msp' = {
                      pbapply::pboptions(type='timer', char='+')
                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readMSP(x, ...)
                      })
                    },

                    'mzXML' = {
                      pbapply::pboptions(type='timer', char='+')
                      # ms2_data <- pbapply::pblapply(ms2_file, function(x){
                      #   readMzXML(x, ...)
                      # })
                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readMzML2(x, ...)
                      })

                    },

                    'mzML' = {
                      pbapply::pboptions(type='timer', char='+')
                      # ms2_data <- pbapply::pblapply(ms2_file, function(x){
                      #   readMzML(x, ...)
                      # })

                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readMzML2(x, ...)
                      })

                    }
             )

             if (!(ms2_type %in% c('mzML', 'mzXML'))) {
               file_name <- basename(ms2_file)
               ms2_data <- lapply(seq_along(ms2_data), function(i){
                 temp_ms2 <- ms2_data[[i]]
                 result <- lapply(temp_ms2, function(x){
                   x$info <- c(x$info, 'file_name_idx' = i)
                   return(x)
                 })
                 return(result)
               })
             }


             ms2_data <- do.call(c, ms2_data)
           }
)



#   readMGF ----------------------------------------------------------------------

#' @title readMGF
#' @author Zhiwei Zhou
#' @param file the name of mgf file
#' @export

# path <- 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_201007_metdna2'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# test <- readMGF(file.path(path, temp_files[1]))

# path <- '~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData_mgf/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# test <- readMGF(file = file.path(path, temp_files[1]))

setGeneric('readMGF',
           def = function(file,
                          ...){

             mgf.data <- ListMGF(file)
             whether_has_ccs <- stringr::str_detect(mgf.data[[1]], pattern = '^CCS') %>% any()

             if (!whether_has_ccs) {
               # mgf.data <- ListMGF(file)
               # nl.spec <- grep('^\\d', mgf.data)
               nl.spec <- lapply(mgf.data, function(x) grep('^\\d', x)) # select spec index
               info.mz <- lapply(mgf.data, function(x) grep('^PEPMASS', x, value = T))
               info.rt <- lapply(mgf.data, function(x) grep('^RTINSECONDS', x, value = T))

               info.mz <- sapply(info.mz, function(x){
                 temp <- gsub(pattern = "\\w+=", "", x)
                 if (stringr::str_detect(temp, pattern = "\ ")) {
                   temp <- as.numeric(substr(temp, start = 1, stop = regexpr("\ ", temp)))
                 }

                 temp <- as.numeric(temp)
                 # temp <- round(as.numeric(temp), digits = 4)
               })

               info.rt <- unlist(info.rt)
               # info.rt <- round(as.numeric(gsub(pattern = "\\w+=", "", info.rt)), digits = 0)
               info.rt <- as.numeric(gsub(pattern = "\\w+=", "", info.rt))

               spec <- mapply(function(x, y){
                 list(do.call(rbind, strsplit(x[y], split = " ")))
               },
               x = mgf.data,
               y = nl.spec)

               spec <- lapply(seq(length(info.mz)), function(i){
                 # data.frame(mz=as.numeric(spec[1:10, i]), intensity=as.numeric(spec[11:20, i]), stringsAsFactors = F)
                 data.frame(mz=as.numeric(spec[[i]][, 1]), intensity=as.numeric(spec[[i]][, 2]), stringsAsFactors = F)
               })

               ms2 <- lapply(seq(length(info.mz)), function(i){
                 temp.info <- c(info.mz[i], info.rt[i])
                 names(temp.info) <- c("mz", "rt")
                 temp.spec <- spec[[i]]
                 temp.result <- list(info=temp.info, spec=temp.spec)
               })

               return(ms2)
             } else {
               # mgf.data <- ListMGF(file)
               # nl.spec <- grep('^\\d', mgf.data)
               nl.spec <- lapply(mgf.data, function(x) grep('^\\d', x)) # select spec index
               info.mz <- lapply(mgf.data, function(x) grep('^PEPMASS', x, value = T))
               info.rt <- lapply(mgf.data, function(x) grep('^RTINSECONDS', x, value = T))
               info.ccs <- lapply(mgf.data, function(x) grep('^CCS', x, value = T))

               info.mz <- sapply(info.mz, function(x){
                 temp <- gsub(pattern = "\\w+=", "", x)
                 if (stringr::str_detect(temp, pattern = "\ ")) {
                   temp <- as.numeric(substr(temp, start = 1, stop = regexpr("\ ", temp)))
                 }

                 temp <- as.numeric(temp)
                 # temp <- round(as.numeric(temp), digits = 4)
               })

               info.rt <- unlist(info.rt)
               # info.rt <- round(as.numeric(gsub(pattern = "\\w+=", "", info.rt)), digits = 0)
               info.rt <- as.numeric(gsub(pattern = "\\w+=", "", info.rt))

               info.ccs <- unlist(info.ccs)
               info.ccs <- as.numeric(gsub(pattern = "\\w+=", "", info.ccs))


               spec <- mapply(function(x, y){
                 list(do.call(rbind, strsplit(x[y], split = " ")))
               },
               x = mgf.data,
               y = nl.spec)

               spec <- lapply(seq(length(info.mz)), function(i){
                 # data.frame(mz=as.numeric(spec[1:10, i]), intensity=as.numeric(spec[11:20, i]), stringsAsFactors = F)
                 data.frame(mz=as.numeric(spec[[i]][, 1]), intensity=as.numeric(spec[[i]][, 2]), stringsAsFactors = F)
               })

               ms2 <- lapply(seq(length(info.mz)), function(i){
                 temp.info <- c(info.mz[i], info.rt[i], info.ccs[i])
                 names(temp.info) <- c("mz", "rt", 'ccs')
                 temp.spec <- spec[[i]]
                 temp.result <- list(info=temp.info, spec=temp.spec)
               })

               return(ms2)
             }

           })


#' @title ListMGF
#' @author Zhiwei Zhou
#' @param file the name of mgf file

setGeneric('ListMGF',
           def = function(file){
             mgf.data <- readLines(file)
             nl.rec.new <- 1
             idx.rec <- 1
             rec.list <- list()
             for(nl in 1:length(mgf.data))
             {
               if(mgf.data[nl]=="END IONS")
               {
                 rec.list[idx.rec] <- list(Compound = mgf.data[nl.rec.new : nl])
                 nl.rec.new <- nl + 1
                 idx.rec <- idx.rec + 1
               }
             }
             rec.list
           })



#   readMSP ----------------------------------------------------------------------

#' @title readMSP
#' @description read MSP spectra files
#' @author Zhiwei Zhou
#' \email{zhouzw@@stanford.edu}
#' @param file the file name
#' @param mode standard: extract name, mz and RT from the file, which fit for MSP data exported from QI software; all: extract all information in the file
#' @return
#' A list object. info: MS1 information, including ms1, rt etc.; spec: the MS/MS spectrum
#' @examples
#' test <- readMSP(file = 'F:/01 MetIMMS/00 data processing/190515 external validation msms data extraction/zhumetlib_validation_pos_20v_190520.msp', mode = 'all')

setGeneric('readMSP', function(file,
                               mode = c('all', 'standard'),
                               source = c('MetAnalyzer', 'MSDIAL', 'Other')) {
  # devtools::use_package('dplyr')

  mode <- match.arg(mode)
  source <- match.arg(source)
  msp.data.list <- ListDB(file)
  nr.num.pk <- grep('Num Peaks', stringr::str_to_title(msp.data.list[[1]]))
  info.spec <- lapply(msp.data.list, function(msp.data) {
    info.list <- strsplit(msp.data[1:nr.num.pk], split = ': ')
    info <- do.call(cbind, info.list)
    colnames(info) <- info[1, ]

    if (mode=='standard') {
      mz <- round(as.numeric(info[2,3]), digits = 4)
      rt <- round(as.numeric(strsplit(info[2, "Comment"], split = "_")[[1]][1])*60, digits = 0)
      name <- info[2, "Comment"]
      # info <- matrix(c(mz, rt), ncol = 2)
      info <- data.frame(mz=mz, rt=rt)
      rownames(info) <- name
      # colnames(info) <- c("mz", "rt")
    } else {
      info <- as.data.frame(tibble::as.tibble(info))
      info <- info[-1,,drop=F]

      if ('comment' %in% tolower(colnames(info))) {
        info$Comment <- stringr::str_replace(info$Comment, pattern = '\t', replacement = '')
      }

      if (info$NAME == 'Unknown') {
        source <- 'MSDIAL'
      }

      # change colnames for MSP file from MetAnalyzer
      switch (source,
              'MetAnalyzer' = {
                info <- info %>%
                  dplyr::rename(name = NAME,
                                mz = PRECURSORMZ,
                                npeaks = `Num Peaks`) %>%
                  dplyr::mutate(mz = as.numeric(mz))
              },
              'MSDIAL' = {
                info <- info %>%
                  dplyr::rename(name = NAME,
                                mz = PRECURSORMZ,
                                npeaks = `Num Peaks`) %>%
                  dplyr::mutate(mz = as.numeric(mz),
                                name = Comment)
              }
      )

      rownames(info) <- NULL
    }

    # if NULL spectra exported, return NULL (changed for MSDIAL)
    if (length(msp.data) <= nr.num.pk) {
      return(NULL)
    }

    spec.list <- strsplit(msp.data[(nr.num.pk+1):length(msp.data)], ' |\\t')
    spec.list <- lapply(spec.list, as.numeric)
    spec <- do.call(rbind, spec.list)[, c(1, 2), drop = FALSE]
    colnames(spec) <- c('mz', 'intensity')
    # spec <- list('spec' = spec)

    # list('info' = info[-1, , drop = FALSE],
    #      'spec' = spec)

    list('info' = info,
         'spec' = spec)

  })

  # remove NULL spectra
  idx_rm <- sapply(info.spec, function(x){length(x) == 0}) %>% which()
  if (length(idx_rm) > 0) {
    info.spec <- info.spec[-idx_rm]
  }

  return(info.spec)
})



setGeneric('ListDB', function(file) {
  msp.data <- readLines(file)
  nl.db.new <- 1
  idx.db <- 1
  db.list <- list()
  len.data <- length(msp.data)
  for(nl in 1:len.data)
  {
    if(msp.data[nl]=="") {
      db.list[idx.db] <- list(Compound=msp.data[nl.db.new:(nl-1)])
      nl.db.new <- nl + 1
      idx.db <- idx.db + 1
    } else if (nl == len.data){
      db.list[idx.db] <- list(Compound=msp.data[nl.db.new:nl])
    }
  }
  db.list
})


#   readCEF ----------------------------------------------------------------------
#' @title readCEF
#' @description  Read a CEF file and return a list of spectra for all feature
#'  with feature information
#' @param file path of the CEF file
#' @export

# file <- 'H:/00_projects/03_MetDNA2/00_data/20201012_metdna2_developmemt/CEF/mice_liver_pos/QC_01_MA-d3-c3_SR.cef'
# test <- readCEF(file)

setGeneric('readCEF', function(file) {
  cef.data.list <- ListCEF(file)
  # nr.num.pk: the number of
  # nr.num.pk <- grep('Num Peaks', cef.data.list[[1]])
  info.spec <- lapply(cef.data.list, function(cef.data) {

    # Extract basic information
    #  line 13, extract mz
    #  line 2, extract string rt, ccs "rt=\"\\d+.\\d+", example "rt=\"5.934"
    #    rt in seconds

    mz <- stringr::str_extract(string = cef.data[13], "\\d+.\\d+")
    mz <- round(as.numeric(mz), 4)

    rt <- stringr::str_extract(string = cef.data[2], "rt=\"\\d+.\\d+")
    ccs <- stringr::str_extract(string = cef.data[2], "ccs=\"\\d+.\\d+")

    rt <- round(as.numeric(stringr::str_extract(string = rt, pattern = "\\d+.\\d+"))*60, 0)
    ccs <- round(as.numeric(stringr::str_extract(string = ccs, pattern = "\\d+.\\d+")), 1)

    # info <- data.frame(mz=mz, rt=rt, ccs=ccs)
    info <- c(mz, rt, ccs)
    names(info) <- c('mz', 'rt', 'ccs')

    # extract MS/MS spectrum -------------------
    idx.spec.start <- stringr::str_locate(string = cef.data, pattern = "\\s+<MSPeaks>")
    idx.spec.end <- stringr::str_locate(string = cef.data, pattern = "\\s+</MSPeaks>")

    # search the label "<MSPeaks>"
    # 1st is the product ion, 2nd is the isopote

    idx.spec.start <- which(idx.spec.start[,1]>0)[1]+1
    idx.spec.end <- which(idx.spec.end[,1]>0)[1]-1

    spec.list.mz<- stringr::str_extract(cef.data[idx.spec.start:idx.spec.end],
                                        pattern = "x=\"\\d+.\\d+")
    spec.list.int <- stringr::str_extract(cef.data[idx.spec.start:idx.spec.end],
                                          pattern = "y=\"\\d+.\\d+")

    spec.list.mz <- round(as.numeric(stringr::str_extract(spec.list.mz, "\\d+.\\d+")),4)
    spec.list.int <- round(as.numeric(stringr::str_extract(spec.list.int, "\\d+.\\d+")),0)


    spec <- data.frame(mz=spec.list.mz, intensity=spec.list.int)
    spec <- as.matrix(spec)

    list('info' = info,
         'spec' = spec)

  })

  return(info.spec)
})



#' @title ListCEF
#' @author Zhiwei Zhou
#' @param file the name of CEF file

setGeneric('ListCEF',
           def = function(file){
             # cef.data <- readLines(file)
             cef.data <- readr::read_lines(file)

             idx.start <- stringr::str_locate(string = cef.data, "<Compound mppid=\"")
             idx.end <- stringr::str_locate(string = cef.data, "</Compound>")

             idx.start <- which(as.numeric(idx.start[,1])>0)
             idx.end <- which(as.numeric(idx.end[,1])>0)

             rec.list <- mapply(function(x, y){
               result <- cef.data[x:y]
               return(result)
             },
             x = idx.start,
             y = idx.end)

             return(rec.list)
           })


#   readMzXML ----------------------------------------------------------------------

#' @title readMzXML
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param file
#' @export
# file <- 'H:/00_projects/03_MetDNA2/00_data/20201012_metdna2_developmemt/mzXML/POS/Sample2-W30 POS06-POS-W30 POS06.mzXML'
# test <- readMzXML(file)

setGeneric('readMzXML',
           def = function(file){
             mzxml.data <- mzR::openMSfile(file)
             mzxml.info <- mzR::header(mzxml.data)
             mzxml.peak <- mzR::peaks(mzxml.data)

             ms2.idx <- which(mzxml.info$msLevel == 2)
             ms2.info <- mzxml.info[ms2.idx, c("precursorMZ", "retentionTime")]
             ms2.info <- apply(ms2.info, 1, list)
             ms2.spec <- mzxml.peak[ms2.idx]

             result.ms2 <- mapply(function(x, y){
               # temp_info <- data.frame(mz = x[[1]][1],
               #                         rt = x[[1]][2],
               #                         stringsAsFactors = FALSE)

               temp_info <- as.numeric(c(x[[1]][1], x[[1]][2]))
               names(temp_info) <- c("mz", "rt")

               temp_spec <- y
               colnames(temp_spec) <- c("mz", "intensity")

               result <- list('info' = temp_info,
                              'spec' = temp_spec)

               return(result)
             },
             x = ms2.info,
             y = ms2.spec,
             SIMPLIFY = FALSE)

             return(result.ms2)

           })


#   readMzML -------------------------------------------------------------------

#' @title readMzML
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param file
#' @export
# test <- readMzXML(file)

# file <- '~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData_DoddLabMetabolomics/ms2/B011_PoolQC-DDA-r001.mzML'
# test <- readMzXML(file)

setGeneric('readMzML',
           def = function(file){
             mzxml.data <- mzR::openMSfile(file)
             mzxml.info <- mzR::header(mzxml.data)
             mzxml.peak <- mzR::peaks(mzxml.data)

             ms2.idx <- which(mzxml.info$msLevel == 2)
             ms2.info <- mzxml.info[ms2.idx, c("precursorMZ", "retentionTime")]
             ms2.info <- apply(ms2.info, 1, list)
             ms2.spec <- mzxml.peak[ms2.idx]

             result.ms2 <- mapply(function(x, y){
               # temp_info <- data.frame(mz = x[[1]][1],
               #                         rt = x[[1]][2],
               #                         stringsAsFactors = FALSE)

               temp_info <- as.numeric(c(x[[1]][1], x[[1]][2]))
               names(temp_info) <- c("mz", "rt")

               temp_spec <- y
               colnames(temp_spec) <- c("mz", "intensity")

               result <- list('info' = temp_info,
                              'spec' = temp_spec)

               return(result)
             },
             x = ms2.info,
             y = ms2.spec,
             SIMPLIFY = FALSE)

             return(result.ms2)

           })




#   readMzML2 --------------------------------------------------------------------
#' @title readMzML2
#' @author Zhiwei Zhou
#' @param file mzml or mzxml file. ms2 files
#' @param isolation_window 'narrow', 'medium', 'wide'
#' @param is_purification whether do the ms2 purification
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' file <- '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_mzXML/Stanford_HILIC_pos/SU-A01_1.mzXML'
#' file <- '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzml_hilic_pos/StanfordPool1-2.mzML'
#'
#' test <- retrieve_raw_ms2(file = 'SU-A01_1.mzXML',
#'                          isolation_window = 'narrow',
#'                          dir_path = '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_mzXML/Stanford_HILIC_pos')
#' }


# file <- '~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData_DoddLabMetabolomics/ms2/B011_PoolQC-DDA-r001.mzML'
# test3 <- readMzML2(file = file)
# file <- '~/Project/05_coworker/Aedan/230727_untargeted_data_processing/HILIC_pos/group1/S1.mzML'
# test3 <- readMzML2(file = file)

readMzML2 <- function(file,
                      isolation_window = c('narrow', 'medium', 'wide'),
                      # dir_path = '.',
                      is_purification = FALSE) {
  isolation_window <- match.arg(isolation_window)
  switch (isolation_window,
          'narrow' = {
            isolation_range <- 0.65
          },
          'medium' = {
            isolation_range <- 4
          },
          'wide' = {
            isolation_range <- 9
          }
  )

  mzml.data <- mzR::openMSfile(file)
  mzml.info <- mzR::header(mzml.data)
  mzml.peak <- mzR::peaks(mzml.data)

  mzml_info_ms1 <- mzml.info %>%
    dplyr::filter(msLevel == 1)
  mzml_info_ms2 <- mzml.info %>%
    dplyr::filter(msLevel == 2)


  # In targeted ms/ms, no precursorScanNum provided. Assign the closed ms1 scan as precursor scan
  precursor_scan_num <- mzml_info_ms2$precursorScanNum
  if (all(is.na(precursor_scan_num))) {
    temp_position <- outer(mzml_info_ms1$seqNum, mzml_info_ms2$seqNum, function(x, y){
      x - y
    })
    # the ms1 scan number must less than ms2 scan. assign 10000 if the ms1 scan larger than ms2 scan
    precursor_idx <- apply(temp_position, 2, function(x){
      x[x>0] <- 10000
      which.min(abs(x))
    })

    precursor_seq_num <- mzml_info_ms1$seqNum[precursor_idx]
    precursor_rt <- mzml_info_ms1$retentionTime[precursor_idx]
  } else {
    precursor_seq_num <- match(precursor_scan_num, mzml_info_ms1$acquisitionNum) %>% mzml_info_ms1$seqNum[.]
    precursor_rt <- match(precursor_scan_num, mzml_info_ms1$acquisitionNum) %>% mzml_info_ms1$retentionTime[.]
  }


  precursor_info <- mapply(function(seq_num, mz_precursor){
    temp_spec <- mzml.peak[[seq_num]]

    temp <- which(abs(temp_spec[,'mz'] - mz_precursor) <= 0.0015)
    if (length(temp) < 1) {
      precursor_info <- tibble::tibble(precursor_mz = mz_precursor,
                                       intensity = 0,
                                       purity = -1,
                                       ms1_scan = seq_num,
                                       mz_error = -1)

      return(precursor_info)
    }

    # select the correct precursor
    idx <- which.min(abs(temp_spec[,'mz'] - mz_precursor))
    int_precursor <- temp_spec[,'intensity'][idx]
    mz_exp <- temp_spec[,'mz'][idx]
    mz_error <- abs(temp_spec[,'mz'] - mz_precursor)[idx]

    # calculate the precursor purity
    temp_idx <- which((temp_spec[,'mz'] >= (mz_precursor - isolation_range)) &
                        (temp_spec[,'mz'] <= (mz_precursor + isolation_range)))
    purity_precursor <- int_precursor/sum(temp_spec[,'intensity'][temp_idx])

    # generate precursor info
    precursor_info <- tibble::tibble(precursor_mz = mz_precursor,
                                     intensity = int_precursor,
                                     purity = purity_precursor,
                                     ms1_scan = seq_num,
                                     mz_error = mz_error)

    return(precursor_info)
  },
  seq_num = precursor_seq_num,
  mz_precursor = mzml_info_ms2$precursorMZ,
  SIMPLIFY = FALSE)

  precursor_info <- precursor_info %>% dplyr::bind_rows()

  ms2_info <- mzml_info_ms2 %>%
    dplyr::bind_cols(precursor_info) %>%
    dplyr::mutate(precursor_rt = precursor_rt) %>%
    dplyr::mutate(precursor_mz = precursor_mz,
                  ms2_rt = retentionTime,
                  ce = collisionEnergy,
                  ms2_scan = seqNum,
                  precuror_int = intensity,
                  ms2_file = file) %>%
    dplyr::select(precursor_mz,
                  precursor_rt,
                  precuror_int,
                  ms2_rt,
                  ce,
                  purity,
                  ms1_scan,
                  ms2_scan,
                  ms2_file) %>%
    dplyr::rename(mz = precursor_mz,
                  rt = precursor_rt)

  ms2_spec <- mzml.peak[ms2_info$ms2_scan]

  if (is_purification) {
    purified_spec <- mapply(function(x, y){
      DoddLabMetID::purifyMs2(x,
                              mz_precursor = y,
                              is_include_precursor = TRUE,
                              is_remove_ring_effect = TRUE,
                              ppm_precursor_filter = 10,
                              int_ms2_min_abs = 50,
                              int_ms2_min_relative = 0.01)
    },
    x = ms2_spec,
    y = ms2_info$precursor_mz,
    SIMPLIFY = FALSE)

    idx_eff <- which(sapply(purified_spec, length) > 0)
    if (length(idx_eff) > 0) {
      purified_spec <- purified_spec[idx_eff]
      ms2_info <- ms2_info[idx_eff,]
      result <- list('info' = ms2_info,
                     'spec' = purified_spec)
      return(result)
    } else {
      result <- list(spec_info = NULL,
                     best_spec = NULL)
      return(result)
    }
  }

  result <- lapply(seq(nrow(ms2_info)), function(i){
    temp_info <- ms2_info %>% dplyr::slice(i)
    temp_spec <- ms2_spec[[i]]
    result <- list('info' = temp_info,
                   'spec' = temp_spec)
  })

  return(result)
}



################################################################################
# integrate_ms2 ----------------------------------------------------------------

#' @title integrate_ms2
#' @description Purify and integrate multiple MS/MS files
#' @author Zhiwei Zhou
#' \email{zhouzw@@stanford.edu}
#' @param ms2_data ms2 data
#' @param is_include_precursor remove precursor peaks in raw spectra or not. Default: FALSE
#' @param is_deisotope remove isotope peaks in raw spectra or not. Default: TRUE
#' @param int_ms2_min_abs the minmium intensity of ms2 fragment. Default: 10
#' @param int_ms2_min_relative the minmium relative of ms2 fragment. Default: 0.01 ((1%))
#' @param ppm_precursor_filter the m/z tolerance to rexmove precursor ion in spectra. Default: 20
#' @param extract_ms2_file extract ms2 file. Default: NULL
#' @param mz_range_ms2 the range of ms2 data.
#' @export

# path <- '~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData_mgf/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# ms2_data <- read_ms2(ms2_file = file.path(path, temp_files),
#                      ms2_type = 'mgf')
# test <- integrate_ms2(ms2_data = ms2_data,
#                       ms2_file = file.path(path, temp_files),
#                       ms2_type = 'mgf',
#                       is_include_precursor = TRUE,
#                       is_deisotope = FALSE,
#                       int_ms2_min_abs = 0,
#                       int_ms2_min_relative = 0.01,
#                       ppm_precursor_filter = 25,
#                       mz_range_ms2 = NULL)


# path <- '~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData_DoddLabMetabolomics/ms2/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mzML') %>% dir(path)[.]
# ms2_data <- read_ms2(ms2_file = file.path(path, temp_files),
#                      ms2_type = 'mzML')
# ms2_data <- integrate_ms2(ms2_data = ms2_data,
#                       ms2_file = file.path(path, temp_files),
#                       ms2_type = 'mzML',
#                       is_include_precursor = TRUE,
#                       is_deisotope = FALSE,
#                       int_ms2_min_abs = 0,
#                       int_ms2_min_relative = 0.01,
#                       ppm_precursor_filter = 25,
#                       mz_range_ms2 = NULL)


setGeneric(name = 'integrate_ms2',
           def = function(
    ms2_data,
    ms2_file,
    ms2_type = c("mgf", "msp", "cef", 'mzXML', 'mzML'),
    is_include_precursor = TRUE,
    is_deisotope = FALSE,
    int_ms2_min_abs = 0,
    int_ms2_min_relative = 0.01,
    ppm_precursor_filter = 20,
    mz_range_ms2=NULL,
    ...
           ){
             match.arg(ms2_type)
             # browser()

             cat('Purify and integrate MS/MS spectra\n')
             tune_MS2_spec <- pbapply::pblapply(seq_along(ms2_data), function(i){
               # cat(i, ' ')
               result <- purifyMs2(spec = ms2_data[[i]]$spec,
                                   mz_precursor = as.numeric(ms2_data[[i]]$info['mz']),
                                   is_include_precursor = is_include_precursor,
                                   is_deisotope = is_deisotope,
                                   int_ms2_min_abs = int_ms2_min_abs,
                                   int_ms2_min_relative = int_ms2_min_relative,
                                   ppm_precursor_filter = ppm_precursor_filter,
                                   mz_range_ms2 = mz_range_ms2)


               return(result)
             })

             # add file names
             info <- lapply(seq(length(ms2_data)), function(i){
               ms2_data[[i]]$info
             })

             info <- do.call(rbind, info)

             if (ms2_type == "mgf") {
               info <- data.frame(mz=as.numeric(info[,1]),
                                  rt=as.numeric(info[,2]),
                                  filename = as.numeric(info[,3]),
                                  stringsAsFactors = F)

               file_name_list <- basename(ms2_file)
               info <- info %>%
                 dplyr::mutate(filename = file_name_list[info$filename])
             }

             if (ms2_type == "msp") {
               temp <- rownames(info)
               info <- data.frame(name=as.character(info[,1]),
                                  mz=as.numeric(info[,2]),
                                  filename = as.numeric(info[,3]),
                                  stringsAsFactors = F)
               rownames(info) <- temp
             }

             if (ms2_type == 'mzXML') {
               # temp <- rownames(info)
               # info <- data.frame(mz=as.numeric(info[,1]),
               #                    rt=as.numeric(info[,2]),
               #                    filename = as.numeric(info[,3]),
               #                    stringsAsFactors = F)

               info <- info %>%
                 dplyr::rename(filename = ms2_file) %>%
                 dplyr::mutate(filename = basename(filename))
             }

             if (ms2_type == 'mzML') {
               # temp <- rownames(info)
               # info <- data.frame(mz=as.numeric(info[,1]),
               #                    rt=as.numeric(info[,2]),
               #                    filename = as.numeric(info[,3]),
               #                    stringsAsFactors = F)

               info <- info %>%
                 dplyr::rename(filename = ms2_file) %>%
                 dplyr::mutate(filename = basename(filename))
             }

             idx_null <- which(sapply(tune_MS2_spec, is.null))

             if (length(idx_null) > 0) {
               info <- info[-idx_null,,drop=F]
               tune_MS2_spec <- tune_MS2_spec[-idx_null]
             }

             result <- list(info=info, spec=tune_MS2_spec)
             return(result)

           }
)


# purifyMs2 -------------------------------------------------------------------


#' @title purifyMs2
#' @author Zhiwei Zhou
#' @param spec the spectra matrix. column1: mz; column2: intensity
#' @param mz_precursor the precursor m/z. Numeric.
#' @param is_include_precursor remove precursor peak in raw spectra or not_ Default: False
#' @param mz_range_ms2 the range of m/z. Default: NULL
#' @param int_ms2_min_abs the minmium intensity of ms2 fragment_ Default: 10
#' @param int_ms2_min_relative the minmium relative of ms2 fragment_ Default: 0.03 ((3%))
#' @param ppm_precursor_filter the m/z tolerance to remove precursor ion in spectra_ Default: 20
#' @export

setGeneric(name = 'purifyMs2',
           def = function(
    spec,
    mz_precursor,
    is_include_precursor = FALSE,
    is_deisotope = FALSE,
    is_remove_ring_effect = TRUE,
    mz_range_ms2 = NULL,
    int_ms2_min_abs = 10,
    int_ms2_min_relative = 0.03,
    ppm_precursor_filter = 20
           ){
             # browser()
             spec <- spec[order(spec[, 'mz']), , drop = FALSE]

             # considering precursor ion
             if (missing(mz_precursor)) {
               mz_precursor <- max(spec[, 'mz'])
               mz_precursor <- as.numeric(mz_precursor)
             }

             mz_precursor_range <- getMzRange(mz_precursor, ppm_precursor_filter, mz_ppm_thr = 400)
             idx_mz_precursor_range <- ifelse(is_include_precursor, 2, 1)

             #change mz range depend precusor include or not
             mz_cutoff <- mz_precursor_range[idx_mz_precursor_range]
             spec <- spec[spec[,'mz'] < mz_cutoff, , drop = FALSE]

             if (!is.null(mz_range_ms2)) {
               nr_keep <- which(spec[, 'mz'] >= mz_range_ms2[1] &
                                  spec[, 'mz'] <= mz_range_ms2[2])
               if (length(nr_keep) > 0) {
                 spec <- spec[nr_keep, , drop = FALSE]
               }
               else {
                 return()
               }
             }

             # discarding low intensity spec (1% highest int and int.ms2.min.abs)
             int_cutoff <- max(max(spec[, 'intensity']) *
                                 int_ms2_min_relative,
                               int_ms2_min_abs)
             spec <- spec[spec[, 'intensity'] >= int_cutoff, , drop = FALSE]
             if (nrow(spec) == 0) {
               return()
             }

             if (is_deisotope) {
               spec <- removeIsotopes(spec)
             }

             # discarding ring effects
             if (is_remove_ring_effect) {
               spec <- removeRingEffect(spec)
             }

             return(spec)

           }
)



#' @title removeIsotopes
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param spec the matrix spectra. 1st column: "mz", 2nd column: "intensity"

setGeneric('removeIsotopes',
           def = function(
    spec
           ){
             # transform spec to data.frame
             temp_spec <- data.frame(spec, annotation='M', stringsAsFactors = F)

             is_filter_isotope <- TRUE

             while (is_filter_isotope) {
               idx_max <- which.max(temp_spec[,'intensity'])
               # cat(idx.max); cat(' ')

               mono_mz <- temp_spec[idx_max, 'mz']
               mono_int <- temp_spec[idx_max, 'intensity']

               temp_spec$intensity[idx_max] <- 0

               # generate potential isotope list
               isotopes_list <- mono_mz + c(1.0003, 2.0044)

               isotopes_list <- lapply(isotopes_list, function(x){
                 c(x, x + c(-1, 1) * 0.003)
               })

               isotopes_list <- do.call(rbind, isotopes_list)
               colnames(isotopes_list) <- c('mz', 'min', 'max')

               # m/z match
               # if multiple spectra were

               temp_idx_isotope <- lapply(seq(nrow(isotopes_list)), function(i){
                 idx_isotope <- which(temp_spec[,'mz'] >= isotopes_list[i,'min'] & temp_spec[,'mz'] <= isotopes_list[i,'max'])

                 if (length(idx_isotope) > 1) {
                   temp_mz_error <- abs(temp_spec[idx_isotope, 'mz'] - isotopes_list[i,'mz'])
                   idx_isotope <- idx_isotope[which.min(temp_mz_error)]
                 }

                 idx_isotope
               })


               # judge the intensity relationship
               # if matched M+1 and M+2, intensity relationship: M > M+1 > M+2
               # if matched M+1, intensity relationship: M > M+1

               if (length(temp_idx_isotope[[1]]) > 0 & length(temp_idx_isotope[[2]]) > 0) {
                 temp_m_1_int <- temp_spec$intensity[temp_idx_isotope[[1]]]
                 temp_m_2_int <- temp_spec$intensity[temp_idx_isotope[[2]]]

                 # intensity decrease
                 if (mono_int > temp_m_1_int & temp_m_1_int > temp_m_2_int) {
                   idx_isotope <- c(temp_idx_isotope[[1]], temp_idx_isotope[[2]])
                   temp_spec$annotation[idx_isotope] <- c('M+1', 'M+2')
                   temp_spec$intensity[idx_isotope] <- 0

                   # names(idx_isotope) <- c('M+1', 'M+2')
                 }

                 if (mono_int > temp_m_1_int & temp_m_1_int <= temp_m_2_int) {
                   idx_isotope <- c(temp_idx_isotope[[1]])
                   temp_spec$annotation[idx_isotope] <- c('M+1')
                   temp_spec$intensity[idx_isotope] <- 0
                   # names(idx_isotope) <- c('M+1', 'M+2')
                 }

               }

               if (length(temp_idx_isotope[[1]]) > 0 & length(temp_idx_isotope[[2]]) == 0) {
                 temp_m_1_int <- temp_spec$intensity[temp_idx_isotope[[1]]]

                 if (mono_int > temp_m_1_int) {
                   idx_isotope <- c(temp_idx_isotope[[1]])
                   temp_spec$annotation[idx_isotope] <- c('M+1')
                   temp_spec$intensity[idx_isotope] <- 0
                 }

               }


               if (all(temp_spec$intensity==0)) {
                 is_filter_isotope <- FALSE
               }

             }

             idx_remove <- which(temp_spec$annotation != 'M')
             spec <- spec[-idx_remove, , drop=F]
             spec
           }
)




#' @title removeRingEffect
#' @author Yandong Yin; Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param spec 2-column table. column 1: mz; column 2: intensity
#' @param mz_diff_thr the threold of mz difference. Default: 0.3
#' @param int_rel_thr the threold of relative intensity. Default: 0.2
#' @export

setGeneric(name = 'removeRingEffect',
           def = function(spec,
                          mz_diff_thr = 0.3,
                          int_rel_thr = 0.2){
             nr_ring <- nrow(spec) + 1
             mz <- spec[, 'mz']

             mz_diff <- diff(mz)
             idx_mzdiff <- which(mz_diff <= mz_diff_thr)
             if (length(idx_mzdiff) == 0) {
               return(spec)
             }

             nr_ring_possible <- unique(c(idx_mzdiff, idx_mzdiff + 1))

             # remove ringeffect loop
             while (TRUE) {

               idx_int_max <- which.max(spec[nr_ring_possible, 2])
               nr_int_max <- nr_ring_possible[idx_int_max] # the index of possible Ringeffect ions with maxium intensity
               int_thr <- spec[nr_int_max, 2] * int_rel_thr # the threshold = 0_2*max_int (possible ring)

               mz_diff <- abs(mz[nr_ring_possible[-idx_int_max]] - mz[nr_int_max])
               int <- spec[nr_ring_possible[-idx_int_max], 2]
               nr_ring <- append(nr_ring, nr_ring_possible[-idx_int_max][which(mz_diff <= mz_diff_thr & int <= int_thr)])
               nr_ring_possible <- nr_ring_possible[!nr_ring_possible %in% c(nr_ring, nr_int_max)]
               if (length(nr_ring_possible) == 0) {
                 break # break loop untill satisfy the nr_ring_possible==0
               }
             }

             return(spec[-nr_ring, , drop = FALSE])
           }
)


################################################################################
# combine_ms1_ms2 --------------------------------------------------------------
#' @title combine_ms1_ms2
#' @author Zhiwei Zhou
#' \email{zhouzw@@stanford.edu}
#' @param ms1_data The name of ms1 peak table. Column 1 is "name", Column 2 is "mz" and column 3 is "rt".
#' @param ms2_data The vector of names of ms2 files. MS2 file must be mzXML.
#' @param mz_tol_combine_ms1_ms2 mz tol for ms1 and ms2 data matching.
#' @param rt_tol_combine_ms1_ms2 RT tol for ms1 and ms2 data matching.
#' @param ms2_type The type of MS2 file, default is mzXML.
#' @return Return ms1 and ms2 data.
#' @export


# load('~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData/01_metabolite_annotation_dodd_mz_rt_ms2/00_intermediate_data/ms1_data')

# path <- '~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData_DoddLabMetabolomics/ms2/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mzML') %>% dir(path)[.]
# ms2_data <- read_ms2(ms2_file = file.path(path, temp_files),
#                      ms2_type = 'mzML')
# ms2_data <- integrate_ms2(ms2_data = ms2_data,
#                       ms2_file = file.path(path, temp_files),
#                       ms2_type = 'mzML',
#                       is_include_precursor = TRUE,
#                       is_deisotope = FALSE,
#                       int_ms2_min_abs = 0,
#                       int_ms2_min_relative = 0.01,
#                       ppm_precursor_filter = 25,
#                       mz_range_ms2 = NULL)
# ms2_data_combined <- combine_ms1_ms2(ms1_data = ms1_data,
#                                      ms2_data = ms2_data,
#                                      ms2_type = 'mzML',
#                                      mz_tol_combine_ms1_ms2 = 15,
#                                      rt_tol_combine_ms1_ms2 = 15)

setGeneric(name = "combine_ms1_ms2",
           function(ms1_data,
                    ms2_data,
                    mz_tol_combine_ms1_ms2 = 15,
                    rt_tol_combine_ms1_ms2 = 15,
                    ms2_type = c("mgf", 'msp', "mzXML", 'mzML')){
             # browser()

             ms2_info <- ms2_data$info
             ms2_spec <- ms2_data$spec

             if (ms2_type %in% c("mzXML", 'mzML', 'mgf')) {
               cat("\n"); cat("Match MS1 and MS2 data.(mz tolerance ",
                              mz_tol_combine_ms1_ms2," ppm, RT tolerance ",
                              rt_tol_combine_ms1_ms2, " second).\n",
                              sep = "")

               temp_ms1_info <- ms1_data$info %>% dplyr::select(mz, rt)
               # match_result <- SXTMTmatch(data1 = ms2_info,
               #                            data2 = temp_ms1_info,
               #                            mz.tol = mz_tol_combine_ms1_ms2,
               #                            rt.tol = rt_tol_combine_ms1_ms2,
               #                            rt.error.type = "abs")

               match_result <- mz_rt_match(data1 = ms2_info,
                                           data2 = temp_ms1_info,
                                           mz.tol = mz_tol_combine_ms1_ms2,
                                           rt.tol = rt_tol_combine_ms1_ms2,
                                           rt.error.type = "abs")

               match_result <- match_result %>%
                 tibble::as_tibble() %>%
                 dplyr::mutate(peak_name = ms1_data$info$name[Index2])

               # remove duplicated MS2 spectrum, if one peak has more than 1 ms2 spectrum,
               #    select the biggest of the sum intensity of top 10 fragments.

               unique_name <- match_result$peak_name %>% unique()
               cat("\n")
               cat("Select the most intense MS2 spectrum for one peak.\n")
               remain_idx <- pbapply::pblapply(seq_along(unique_name), function(i){
                 name <- unique_name[i]
                 temp_idx <- match_result %>% dplyr::filter(peak_name == name) %>% dplyr::pull(Index1)
                 # temp_idx <- which(ms2_name == name)
                 if(length(temp_idx) == 1) return(temp_idx)
                 temp_ms2 <- ms2_spec[temp_idx]
                 # temp_ms2 <- lapply(temp_ms2, function(x) x[[2]])
                 temp_int <- lapply(temp_ms2, function(x) {
                   x <- as.data.frame(x)
                   x <- x[order(x[,2], decreasing = TRUE),]
                   if(nrow(x) >= 10) {sum(x[1:10,2])}else{sum(x[,2])}
                 })
                 temp_int <- unlist(temp_int)
                 temp_idx <- temp_idx[which.max(temp_int)]
                 temp_idx
               })

               names(remain_idx) <- unique_name
               remain_idx <- unlist(remain_idx)

               ms2_name <- names(remain_idx)
               ms2_info <- ms2_info[remain_idx,]
               ms2_spec <- ms2_spec[remain_idx]
               rownames(ms2_info) <- ms2_name
               names(ms2_spec) <- ms2_name

               ms2_list <- list(info = ms2_info, spec = ms2_spec)
               return(ms2_list)
             }

             if(ms2_type == "msp"){
               cat("\n")
               cat("Match MS1 and MS2 data according to the name.\n")

               # remove spectra not in peak table
               ms1_info <- ms1_data$info
               idx <- match(ms2_info$name, ms1_info$name) %>% is.na() %>% which()

               if (length(idx) > 0) {
                 ms2_info <- ms2_info[-idx,,drop = FALSE]
                 ms2_spec <- ms2_spec[-idx]
               }

               idx <- match(ms2_info$name, ms1_info$name)
               ms2_name <- ms2_info$name
               ms2_mz <- ms1_info$mz[idx]
               ms2_rt <- ms1_info$rt[idx]
               ms2_info <- data.frame(mz = ms2_mz,
                                      rt = ms2_rt,
                                      filename = ms2_name,
                                      stringsAsFactors = FALSE)
               names(ms2_spec) <- ms2_name

               ms2_list <- list(info = ms2_info, spec = ms2_spec)
               return(ms2_list)

               # # # change ms2 info to MetDNA1 default format
               # ms2 <- lapply(seq_along(ms2_name), function(i){
               #   temp_info <- data.frame('NAME' = ms2_name[i],
               #                           'PRECURSORMZ' = ms2_mz[i],
               #                           'PRECURSORRT' = ms2_rt[i],
               #                           stringsAsFactors = FALSE)
               #   temp_info <- t(temp_info)
               #   rownames(temp_info) <- c('NAME', 'PRECURSORMZ', 'PRECURSORRT')
               #   result <- list(info = temp_info,
               #                  spec = ms2_spec[[i]])
               #
               #   return(result)
               # })
               #
               # names(ms2) <- ms2_name
               # # save(ms2, file = file.path(path, "ms2"), compress = "xz", version = 2)
               #
               # return(ms2)
             }

           })






# split_ms2 --------------------------------------------------------------------

setGeneric(name = 'split_ms2',
           def = function(ms2_data_combined){
             ms2_name <- ms2_data_combined$spec %>% names()
             ms2 <- pbapply::pblapply(seq_along(ms2_name), function(i){
               NAME <- ms2_name[i]
               PRECURSORMZ <- ms2_data_combined$info$mz[i]
               PRECURSORRT <- ms2_data_combined$info$rt[i]
               result_info <- data.frame(NAME, PRECURSORMZ,
                                         PRECURSORRT, stringsAsFactors = FALSE)
               result_info <- t(result_info)
               rownames(result_info) <- c('NAME', 'PRECURSORMZ', 'PRECURSORRT')

               result <- list(info = result_info,
                              spec = ms2_data_combined$spec[[i]])

               return(result)
             })

             names(ms2) <- ms2_name
             return(ms2)
           })



# add_has_ms2_to_SpecAnnotationClass -------------------------------------------
#' @title add_has_ms2_to_SpecAnnotationClass
#' @author Zhiwei Zhou
#' @param obj_class a pecAnnotationClass object
#' @param ms2
#' @export

setGeneric(name = "add_has_ms2_to_SpecAnnotationClass",
           def = function(obj_class, ms2){
             peaks_with_ms2 <- names(ms2)
             result <- lapply(obj_class, function(x){
               if (x@peak_info$name %in% peaks_with_ms2) {
                 x@peak_info$has_ms2 <- 1
               } else {
                 x@peak_info$has_ms2 <- -1
               }

               return(x)
             })

             names(result) <- names(obj_class)
             return(result)
           })


# add_ms2_purifty_to_SpecAnnotationClass ---------------------------------------
#' @title add_ms2_purifty_to_SpecAnnotationClass
#' @author Zhiwei Zhou
#' @param obj_class a pecAnnotationClass object
#' @param ms2
#' @export
setGeneric(name = "add_ms2_purifty_to_SpecAnnotationClass",
           def = function(obj_class, ms2_data_combined){

             # rownames(ms2_data_combined$info) <- names(ms2_data_combined$spec)

             ms2_info <- ms2_data_combined$info
             peaks_with_ms2 <- rownames(ms2_info)
             result <- lapply(obj_class, function(x){
               temp_feature_name <- x@peak_info$name
               if (temp_feature_name %in% peaks_with_ms2) {
                 temp_idx <- which(peaks_with_ms2 == temp_feature_name)
                 x@peak_info$ms2_purity <- ms2_info$purity[temp_idx]
               } else {
                 x@peak_info$ms2_purity <- -1
               }

               return(x)
             })

             names(result) <- names(obj_class)
             return(result)
           })


################################################################################
# mz_rt_match ------------------------------------------------------------------

#' @title mz_rt_match
#' @description Match peaks according to m/z and RT.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data1 First data for matching, first column must be mz
#' and seconod column must be rt.
#' @param data2 Second data for matching, first column must be mz
#' and seconod column must be rt.
#' @param mz.tol mz tol for ms1 and ms2 data matching.
#' @param rt.tol RT tol for ms1 and ms2 data matching.
#' @param rt.error.type RT error is calculated with relative or absolute.
#' @return Return a result which give the matching
#' result of data1 and database.
# #' @export
#' @examples
#' data1 <- data.frame(mz = 1:10, rt = 1:10)
#' data2 <- data.frame(mz = 1:10, rt = 1:10)
#' mz_rt_match(data1, data2, mz.tol = 10)

mz_rt_match <-
  function(data1,
           data2,
           mz.tol,
           rt.tol = 30,
           rt.error.type = c("relative", "abs")) {
    UseMethod("mz_rt_match")
  }

#' @export
mz_rt_match.data.frame <-
  function(data1,
           data2,
           mz.tol,
           rt.tol = 30,
           rt.error.type = c("relative", "abs")) {
    rt.error.type <- match.arg(rt.error.type)
    mz_rt_match_default(data1 = data1,
                        data2 = data2,
                        mz.tol = mz.tol,
                        rt.tol = rt.tol,
                        rt.error.type = rt.error.type)
  }

#' @export
mz_rt_match.matrix <-
  function(data1,
           data2,
           mz.tol,
           rt.tol = 30,
           rt.error.type = c("relative", "abs")) {
    rt.error.type <- match.arg(rt.error.type)
    mz_rt_match_default(data1 = data1,
                        data2 = data2,
                        mz.tol = mz.tol,
                        rt.tol = rt.tol,
                        rt.error.type = rt.error.type)
  }


mz_rt_match_default <-
  function(data1,
           data2,
           mz.tol,
           rt.tol = 30,
           rt.error.type = c("relative", "abs")) {
    rt.error.type <- match.arg(rt.error.type)
    #
    if (nrow(data1) == 0 | nrow(data2) == 0) {
      result <- NULL
      return(result)
    }
    # mz1 <- as.numeric(data1[, 1])
    # rt1 <- as.numeric(data1[, 2])
    info1 <- data1[, c(1, 2)]
    info1 <- apply(info1, 1, list)

    mz2 <- as.numeric(data2[, 1])
    rt2 <- as.numeric(data2[, 2])

    result <- pbapply::pblapply(info1, function(x) {
      temp.mz1 <- x[[1]][[1]]
      temp.rt1 <- x[[1]][[2]]
      mz.error <- abs(temp.mz1 - mz2) * 10 ^ 6 / temp.mz1
      if (rt.error.type == "relative") {
        rt.error <- abs(temp.rt1 - rt2) * 100 / temp.rt1
      } else {
        rt.error <- abs(temp.rt1 - rt2)
      }

      j <- which(mz.error <= mz.tol & rt.error <= rt.tol)
      if (length(j) == 0) {
        matrix(NA, ncol = 7)
      } else {
        cbind(j, temp.mz1, mz2[j],
              mz.error[j], temp.rt1, rt2[j], rt.error[j])
      }
    })

    if (length(result) == 1) {
      result <- cbind(1, result[[1]])
    } else {
      result <- mapply(function(x, y) {
        list(cbind(x, y))
      },
      x <- seq_along(info1),
      y = result)
      result <- do.call(rbind, result)
    }

    result <-
      matrix(result[which(!apply(result, 1, function(x) {
        any(is.na(x))
      })), ], ncol = 8)
    if (nrow(result) == 0) {
      return(NULL)
    }
    colnames(result) <-
      c("Index1",
        "Index2",
        "mz1",
        "mz2",
        "mz error",
        "rt1",
        "rt2",
        "rt error")
    result <- as.data.frame(result)
    result
  }





# keep_one ---------------------------------------------------------------------

#' @title keep_one
#' @description Remove multiple vs. one in result from mz_rt_match function.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param result result from mz_rt_match function.
#' @param according.to According to mz error or rt error?
#' @return Return a result without multiple vs. one.
# #' @export
#' @examples
#' data1 <- data.frame(mz = 1:10, rt = 1:10)
#' data2 <- data.frame(mz = 1:10, rt = 1:10)
#' result <- mz_rt_match(data1, data2, mz.tol = 10)
#'
#' keep_one(result)
keep_one <- function(result,
                     according.to = c("mz.error", "rt.error")) {
  according.to <- match.arg(according.to)
  if (is.null(result)) {
    return(result)
  }

  if (!is.matrix(result) & !is.data.frame(result)) {
    stop("result must be matrix or data.frame.")
  }

  if (ncol(result) != 8) {
    stop("result must from mz_rt_match.")
  }
  if (paste(colnames(result), collapse = ";") !=
      "Index1;Index2;mz1;mz2;mz error;rt1;rt2;rt error") {
    stop("result must from mz_rt_match.")
  }

  duplicated.idx <-
    unique(result$index1[duplicated(result$index1)])
  if (length(duplicated.idx) == 0) {
    return(result)
  }

  for (i in seq_along(duplicated.idx)) {
    temp.idx <- which(result$index1 == duplicated.idx[i])
    temp.result <- result[temp.idx, ]
    if (according.to == "mz.error") {
      temp.idx1 <- temp.idx[which.min(temp.result$mz.error)]
      temp.idx2 <- setdiff(temp.idx, temp.idx1)
      result <- result[-temp.idx2, ]
    }

    if (according.to == "rt.error") {
      temp.idx1 <- temp.idx[which.min(temp.result$rt.error)]
      temp.idx2 <- setdiff(temp.idx, temp.idx1)
      result <- result[-temp.idx2, ]
    }
  }
  result
}
