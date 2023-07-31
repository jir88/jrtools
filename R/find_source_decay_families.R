#' Find putative source decay families in mass spec alignments
#'
#' Some compounds fragment when ionized into the gas phase, forming a source
#' decay family. Mass spec features in these families are essentially duplicates
#' of the parent compound and should be identified and treated specially. This
#' function identifies such families based on very similar retention times,
#' correlated peak areas across all samples, and identical peak shapes. This
#' approach also finds instances where the alignment software has missed adducts,
#' isotopologues, or multiply charged versions of a single compound. It also
#' corrects for occasional alignment software errors.
#'
#' @param msa An ms_alignment object containing the data to be analyzed.
#' @param raw_data_folder A folder where the raw MS data files are stored in mzML format.
#' @param ref_align_file The base file name (no path, no extension) of the MS data file to use as reference for retention time alignment.
#' @param drt_max Maximum retention time difference (in minutes) to consider.
#' @param area_cor_min Minimum peak area correlation to consider.
#' @param xic_tol Mass tolerance for extracted ion chromatogram calculation (XIC) in parts per million (ppm).
#' @param xic_error_max Maximum XIC error to consider (range 0-1).
#' @param bp_param Optional BiocParallelParam object describing parallel processing to be used. Defaults to serial processing if not supplied.
#'
#' @details
#' When comparing peak shapes, this algorithm uses XICs from the file where
#' either compound is highest. This process only uses files marked as 'Sample' or
#' 'Quality Control' in the original alignment so that the peak shapes represent
#' behavior in real samples.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{matched_compound_pairs} \tab Table of all matching pairs of compounds. \cr
#'    \code{sdf_graph} \tab Network object constructed from compound pairs. \cr
#'    \code{sdf_cluster_members} \tab Table of cluster member IDs and associated cluster IDs. \cr
#'    \code{sdf_cluster_centers} \tab Compound IDs of the central member of each compound cluster. \cr
#'    \tab \cr
#'}
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
find_source_decay_families <- function(msa, raw_data_folder, ref_align_file,
                                       drt_max = 0.02, area_cor_min = 0.997,
                                       xic_tol = 5, xic_error_max = 0.3,
                                       bp_param = BiocParallel::SerialParam()) {
  # TODO: Adjust method to use exported CD data instead of ripping it straight
  # out of the alignments? Is there info I can only get from direct access?
  # Peak RT windows might be hardest to get.

  # Pull feature info from alignment ----

  aln_compounds <- msa$unknown_compound_items %>%
    # drop compounds that were excluded by the workflow
    # in this case, these were always dropped by the 'group compounds' node
    # probably because they had poor peak quality or only showed up once or twice
    dplyr::filter(.data$ExcludedBy == -1) %>%
    # drop features called as background by CD
    dplyr::filter(.data$BackgroundStatus == 0) %>%
    # add generic compound names
    # normal CD output truncates these, so we do too
    dplyr::mutate(RoundMW = round(.data$MolecularWeight, digits = 5),
           RoundRT = round(.data$RetentionTime, digits = 3),
           Name = paste0(.data$RoundMW, "@", .data$RoundRT))

  all_areas <- jrtools::extract_peak_areas(msa = msa, ids = aln_compounds$ID)
  # all of these have peak areas, so we can ignore the flags
  all_areas <- all_areas$areas

  # Find features with similar retention times ----
  message("Finding features with similar retention times...")

  delta_rt <- abs(outer(aln_compounds$RetentionTime, aln_compounds$RetentionTime,
                        FUN = "-"))
  delta_rt[lower.tri(delta_rt, diag = TRUE)] <- NA
  rt_matches <- which(delta_rt < drt_max, arr.ind = TRUE)

  # Find features with correlated peak areas ----
  message("Finding features with correlated peak areas...")

  # grab all features that were part of a match and extract peak areas
  unique_idx <- unique(c(unique(rt_matches[, "row"]), unique(rt_matches[, "col"])))
  unique_cids <- aln_compounds$ID[unique_idx]
  unique_areas <- jrtools::extract_peak_areas(msa = msa,
                                              ids = aln_compounds$ID[unique_idx])$areas
  match_idx1 <- match(aln_compounds$ID[rt_matches[, "row"]], colnames(unique_areas))
  match_idx2 <- match(aln_compounds$ID[rt_matches[, "col"]], colnames(unique_areas))

  # stupid gibberish to push the unique areas over to the worker threads
  area_cor_factory <- function() {
    # convert to log10 area because areas are very skewed
    # this also pulls ua into the local environment, so it'll get passed to workers
    ua <- dplyr::mutate(unique_areas, dplyr::across(-.data$StudyFileID, log10))
    f <- function(idx1, idx2) {
      return(cor(x = ua[, idx1], y = ua[, idx2]))
    }
    return(f)
  }
  suppressWarnings({
  df2 <- BiocParallel::bpmapply(FUN = area_cor_factory(),
                  idx1 = match_idx1, idx2 = match_idx2,
                  BPPARAM = bp_param)
  })

  area_cor_matches <- tibble::tibble(CID1 = aln_compounds$ID[rt_matches[, "row"]],
                             CID2 = aln_compounds$ID[rt_matches[, "col"]],
                             Correlation = df2) %>%
    dplyr::filter(.data$Correlation > area_cor_min)

  # Find features with matching peak shapes ----
  message("Finding features with matching peak shapes...")

  ## Pull retention time windows for each match ----

  # get list of files in run time order
  ms_file_order <- msa$input_files %>%
    # drop extra files, only want ones used for integration
    dplyr::filter(.data$SampleType %in% c("Sample", "Quality Control")) %>%
    dplyr::mutate(CreationDate = strptime(.data$CreationDate, "%m/%d/%Y %H:%M:%S")) %>%
    dplyr::arrange(.data$CreationDate)


  # TODO: set this up to run for old data!
  # NOTE: for older alignments, need to use heuristic to choose reference ions
  # match_ions <- jrtools::get_compound_ions(msa = msa, ids = unique(area_cor_matches$CID1)) %>%
  #   # in case of ties, use the more abundant ion
  #   group_by(ConsolidatedUnknownCompoundItemsID, IonDescription) %>%
  #   summarise(IonCount = n(),
  #             IonTotal = sum(Area.ion),
  #             .groups = "drop_last") %>%
  #   # sort more abundant first
  #   arrange(desc(IonTotal), .by_group = TRUE) %>%
  #   # take most used/most abundant if tied
  #   slice_max(order_by = IonCount, n = 1, with_ties = FALSE)


  # get ion info so we can select only reference ions for consideration.
  # for efficiency, we grab all of it at once
  match_ions <- jrtools::get_compound_ions(msa = msa, ids = unique(c(area_cor_matches$CID1,
                                                                     area_cor_matches$CID2))) %>%
    # keep only ID columns and the ion description and charge
    dplyr::select(ConsolidatedUnknownCompoundItemsID, UnknownCompoundInstanceItemsWorkflowID,
           UnknownCompoundInstanceItemsID, UnknownCompoundIonInstanceItemsWorkflowID,
           UnknownCompoundIonInstanceItemsID, IonDescription, Charge)

  # get most intense peaks from the set of files we're willing to use
  match_pks1 <- jrtools::get_chromatogram_peaks(msa = msa,
                                                ids = unique(area_cor_matches$CID1)) %>%
    # keep only reference isotopologues
    dplyr::filter(.data$IsRefPeak == 1) %>%
    dplyr::left_join(match_ions, by = c("ConsolidatedUnknownCompoundItemsID",
                                        "UnknownCompoundInstanceItemsWorkflowID",
                                        "UnknownCompoundInstanceItemsID",
                                        "UnknownCompoundIonInstanceItemsWorkflowID",
                                        "UnknownCompoundIonInstanceItemsID")) %>%
    # only consider peaks using the reference ion for that feature
    dplyr::filter(.data$ReferenceIon == .data$IonDescription) %>%
    dplyr::filter(.data$StudyFileID.peak %in% ms_file_order$StudyFileID) %>%
    dplyr::group_by(.data$ConsolidatedUnknownCompoundItemsID) %>%
    dplyr::slice_max(order_by = .data$Area.peak, n = 1) %>%
    dplyr::mutate(Polarity = sign(.data$Charge))
  match_pks2 <- jrtools::get_chromatogram_peaks(msa = msa,
                                                ids = unique(area_cor_matches$CID2)) %>%
    # keep only reference isotopologues
    dplyr::filter(.data$IsRefPeak == 1) %>%
    dplyr::left_join(match_ions, by = c("ConsolidatedUnknownCompoundItemsID",
                                        "UnknownCompoundInstanceItemsWorkflowID",
                                        "UnknownCompoundInstanceItemsID",
                                        "UnknownCompoundIonInstanceItemsWorkflowID",
                                        "UnknownCompoundIonInstanceItemsID")) %>%
    # only consider peaks using the reference ion for that feature
    dplyr::filter(.data$ReferenceIon == .data$IonDescription) %>%
    dplyr::filter(.data$StudyFileID.peak %in% ms_file_order$StudyFileID) %>%
    dplyr::group_by(.data$ConsolidatedUnknownCompoundItemsID) %>%
    dplyr::slice_max(order_by = .data$Area.peak, n = 1) %>%
    dplyr::mutate(Polarity = sign(.data$Charge))

  # # get most intense peaks from the set of files we're willing to use
  # match_pks1 <- jrtools::get_chromatogram_peaks(msa = msa,
  #                                               ids = unique(area_cor_matches$CID1)) %>%
  #   dplyr::filter(.data$IsRefPeak == 1) %>%
  #   # only consider peaks using the reference ion for that feature.
  #   # reference ion is the one CD used to calculate peak areas.
  #   dplyr::filter(.data$Area.peak == .data$AreaReferenceIon) %>%
  #   dplyr::filter(.data$StudyFileID.peak %in% ms_file_order$StudyFileID) %>%
  #   dplyr::group_by(.data$ConsolidatedUnknownCompoundItemsID) %>%
  #   dplyr::slice_max(order_by = .data$Area.peak, n = 1) %>%
  #   dplyr::mutate(Polarity = sign(as.numeric(stringr::str_sub(.data$ReferenceIon, start = -2))))
  # match_pks2 <- jrtools::get_chromatogram_peaks(msa = msa,
  #                                               ids = unique(area_cor_matches$CID2)) %>%
  #   dplyr::filter(.data$IsRefPeak == 1) %>%
  #   # only consider peaks using the reference ion for that feature.
  #   # reference ion is the one CD used to calculate peak areas.
  #   dplyr::filter(.data$Area.peak == .data$AreaReferenceIon) %>%
  #   dplyr::filter(.data$StudyFileID.peak %in% ms_file_order$StudyFileID) %>%
  #   dplyr::group_by(.data$ConsolidatedUnknownCompoundItemsID) %>%
  #   dplyr::slice_max(order_by = .data$Area.peak, n = 1) %>%
  #   dplyr::mutate(Polarity = sign(as.numeric(stringr::str_sub(.data$ReferenceIon, start = -2))))

  match_pks_bounds <- area_cor_matches %>%
    dplyr::left_join(dplyr::select(match_pks1, "ConsolidatedUnknownCompoundItemsID",
                                   "Polarity", "Mass", "LeftRT", "RightRT",
                                   StudyFileID = "StudyFileID.peak",
                                   MaxArea = "Area.peak"),
                     by = c("CID1" = "ConsolidatedUnknownCompoundItemsID")) %>%
    dplyr::left_join(dplyr::select(match_pks2, "ConsolidatedUnknownCompoundItemsID",
                                   "Polarity", "Mass", "LeftRT", "RightRT",
                                   StudyFileID = "StudyFileID.peak",
                                   MaxArea = "Area.peak"),
                     by = c("CID2" = "ConsolidatedUnknownCompoundItemsID"),
                     suffix = c(".1", ".2")) %>%
    # set bounds to cover both peaks completely
    dplyr::rowwise() %>%
    dplyr::mutate(LeftRT = min(.data$LeftRT.1, .data$LeftRT.2),
           RightRT = max(.data$RightRT.1, .data$RightRT.2)) %>%
    dplyr::ungroup() %>%
    # only consider the file where the more intense feature is most intense
    dplyr::mutate(StudyFileID = dplyr::if_else(.data$MaxArea.1 > .data$MaxArea.2,
                                               .data$StudyFileID.1, .data$StudyFileID.2))

  ## Load raw mass spec data ----
  message("Loading raw mass spec data...")

  base_name <- regexec("\\\\([^\\]+)\\.", ms_file_order$PhysicalFileName)
  base_name <- regmatches(ms_file_order$PhysicalFileName, m = base_name)
  # second part of each match is the actual file name
  base_name <- purrr::map_chr(base_name, 2)
  base_name <- paste0(base_name, ".mzML")

  ms_files <- file.path(raw_data_folder, base_name)

  # load all MS1 data, accessing scans from disk as needed
  ms_data_all <- MSnbase::readMSData(files = ms_files, msLevel. = 1, mode = "onDisk")

  # cwp <- xcms::CentWaveParam(ppm = 25, peakwidth = c(5, 50), noise = 1000)
  # mfp <- xcms::MatchedFilterParam(binSize = 1, fwhm = 5)

  # pull positive and negative data and centroid the spectra
  xdata <- ms_data_all %>%
    MSnbase::filterPolarity(polarity = 1) %>%
    MSnbase::smooth(method = "SavitzkyGolay", halfWindowSize = 4) %>%
    MSnbase::pickPeaks(refineMz = "kNeighbors", k = 5, halfWindowSize = 3) #%>%
    # findChromPeaks(param = mfp, BPPARAM = bp_param)

  xdata_neg <- ms_data_all %>%
    MSnbase::filterPolarity(polarity = 0) %>%
    MSnbase::smooth(method = "SavitzkyGolay", halfWindowSize = 4) %>%
    MSnbase::pickPeaks(refineMz = "kNeighbors", k = 5, halfWindowSize = 3) #%>%
    # findChromPeaks(param = cwp, BPPARAM = par_params)

  # alignment in XCMS is slow, and we've already done it on the instrument software anyway
  # we can grab the preexisting alignment data and shove it into the XCMS objects
  # also saves time because we don't have to align pos/neg mode separately

  # convert MSnbase experiments to XCMSnExp objects so adjustment works right
  xdata_aln <- methods::as(xdata, "XCMSnExp")
  xdata_aln_neg <- methods::as(xdata_neg, "XCMSnExp")

  # grab preexisting alignment data
  rt_adj <- jrtools::get_file_rt_corrections(msa, file_id = ms_file_order$FileID) %>%
    # convert time to seconds
    dplyr::mutate(OriginalRT = 60*.data$OriginalRT,
                  CorrectedRT = 60*.data$CorrectedRT)

  # grab XCMS scan times -- these don't quite align with raw data
  xrt <- xcms::rtime(xdata_aln, adjusted = FALSE, bySample = TRUE)
  xrt_neg <- xcms::rtime(xdata_aln_neg, adjusted = FALSE, bySample = TRUE)
  names(xrt) <- ms_file_order$StudyFileID
  names(xrt_neg) <- ms_file_order$StudyFileID

  # interpolate existing alignment to match XCMS scan times
  cd_rt_adj <- lapply(X = ms_file_order$StudyFileID, function(fid) {
    # get compound discoverer adjustments for this file
    file_rt_adj <- dplyr::filter(rt_adj, .data$StudyFileID == fid)
    # get XCMS scan retention times
    file_xrt <- xrt[[fid]]
    # the new alignment method uses a reference file, so it doesn't get adjusted
    if(nrow(file_rt_adj) == 0) {
      # return the original retention times unaltered
      return(file_xrt)
    } else {
      # interpolate within CD adjusted RTs
      file_xrt_adj <- stats::approx(x = file_rt_adj$OriginalRT, y = file_rt_adj$CorrectedRT,
                             xout = file_xrt, yleft = 0, rule = 1:2)$y
      # apply scan names
      names(file_xrt_adj) <- names(file_xrt)
      return(file_xrt_adj)
    }
  })
  # apply adjusted times
  xcms::adjustedRtime(xdata_aln) <- cd_rt_adj

  # repeat for negative mode
  # interpolate existing alignment to match XCMS scan times
  cd_rt_neg_adj <- lapply(X = ms_file_order$StudyFileID, function(fid) {
    # get compound discoverer adjustments for this file
    file_rt_adj <- dplyr::filter(rt_adj, .data$StudyFileID == fid)
    # get XCMS scan retention times
    file_xrt <- xrt_neg[[fid]]
    # the new alignment method uses a reference file, so it doesn't get adjusted
    if(nrow(file_rt_adj) == 0) {
      # return the original retention times unaltered
      return(file_xrt)
    } else {
      # interpolate within CD adjusted RTs
      file_xrt_adj <- stats::approx(x = file_rt_adj$OriginalRT, y = file_rt_adj$CorrectedRT,
                                    xout = file_xrt, yleft = 0, rule = 1:2)$y
      # apply scan names
      names(file_xrt_adj) <- names(file_xrt)
      return(file_xrt_adj)
    }
  })
  # apply adjusted times
  xcms::adjustedRtime(xdata_aln_neg) <- cd_rt_neg_adj

  # # align retention times across samples
  # # we won't use blanks for RT alignment
  # sub_idx <- which(ms_file_order$SampleType != "Blank")
  # # reference sample index is within the subset indices
  # ref_idx <- match(paste0(ref_align_file, ".mzML"), base_name[sub_idx])
  # obw_par <- xcms::ObiwarpParam(binSize = 0.4,
  #                         subset = sub_idx,
  #                         centerSample = ref_idx)
  # # make sure these run serially -- parallel processing doesn't work for some
  # # unknown reason
  # # get the registry state
  # bpps <- BiocParallel::registered()
  # # force serial operatioin
  # BiocParallel::register(BPPARAM = BiocParallel::SerialParam(), default = TRUE)
  #
  # # have to force alignments into XCMSnExp objects so adjustment works right
  # suppressWarnings({
  # xdata_aln <- xcms::adjustRtime(methods::as(xdata, "XCMSnExp"), param = obw_par)
  # xdata_aln_neg <- xcms::adjustRtime(methods::as(xdata_neg, "XCMSnExp"), param = obw_par)
  # })
  #
  # # now fix the parallel processing registry
  # for (param in rev(bpps)) {
  #   BiocParallel::register(param)
  # }

  ## Calculate XICs ----
  message("Calculating extracted ion chromatograms...")

  # for each compound found in 1+ matches, get its RT and mz windows
  match_mz1 <- match_pks_bounds %>%
    select(CID = "CID1", Polarity = "Polarity.1", Mass = "Mass.1", "LeftRT", "RightRT", "StudyFileID")
  match_mz2 <- match_pks_bounds %>%
    select(CID = "CID2", Polarity = "Polarity.2", Mass = "Mass.2", "LeftRT", "RightRT", "StudyFileID")

  match_mz <- bind_rows(match_mz1, match_mz2) %>%
    dplyr::group_by(.data$CID) %>%
    dplyr::summarise(Polarity = dplyr::first(.data$Polarity),
                     Mass = dplyr::first(.data$Mass),
                     LeftRT = min(.data$LeftRT)*60,
                     RightRT = max(.data$RightRT)*60,
                     # use file where more intense feature is highest
                     StudyFileID = dplyr::first(.data$StudyFileID)) %>%
    dplyr::mutate(mzmin = .data$Mass - xic_tol/1e6*.data$Mass,
                  mzmax = .data$Mass + xic_tol/1e6*.data$Mass)

  match_mz_pos <- match_mz %>%
    dplyr::filter(.data$Polarity == 1)
  match_mz_neg <- match_mz %>%
    dplyr::filter(.data$Polarity == -1)

  # for each compound in 1+ matches, get XICs from all alignment files
  # these throw a zillion warnings about zero values, so we suppress those
  suppressMessages({
    match_xics_pos <- xcms::chromatogram(xdata_aln,
                                          mz = as.matrix(dplyr::select(match_mz_pos, "mzmin", "mzmax")),
                                          rt = as.matrix(dplyr::select(match_mz_pos, "LeftRT", "RightRT")),
                                          msLevel = 1, BPPARAM = bp_param,
                                          missing = 0)
    match_xics_neg <- xcms::chromatogram(xdata_aln_neg,
                                          mz = as.matrix(dplyr::select(match_mz_neg, "mzmin", "mzmax")),
                                          rt = as.matrix(dplyr::select(match_mz_neg, "LeftRT", "RightRT")),
                                          msLevel = 1, BPPARAM = bp_param,
                                          missing = 0)
  })

  # Now pull XICs for each match pair from the file where the more abundant
  # feature is most abundant. This will give us the best peaks for comparing profiles.

  # make 2-column matrix of (compound, file) indices
  xic_file_order <- colnames(match_xics_pos) %>%
    stringr::str_sub(end = -6)
  study_file_key <- msa$input_files %>%
    # drop extra files, only want ones used for integration
    dplyr::filter(.data$SampleType %in% c("Sample", "Quality Control")) %>%
    dplyr::select("FileName", "StudyFileID") %>%
    dplyr::mutate(File = stringr::str_match(.data$FileName, "[\\\\/]([^\\\\/]*)\\.raw")[, 2],
                  XIC_Idx = match(.data$File, xic_file_order))

  # pull XICs from each file of interest for every CID
  high_xics <- bind_rows(match_mz1, match_mz2) %>%
    select(CID, Polarity, StudyFileID) %>%
    distinct()
  # high_xics$XIC <- mapply(cid = high_xics$CID, pol = high_xics$Polarity, SIMPLIFY = FALSE, FUN = function(cid, pol) {
  high_xics$XIC <- lapply(X = 1:nrow(high_xics), FUN = function(r) {
    cid <- high_xics$CID[[r]]
    if(high_xics$Polarity[[r]] == 1) {
      # index of compound (row in XIC object)
      cid_idx <- match(cid, match_mz_pos$CID)
      # index of file (column in XIC object)
      file_idx <- match(high_xics$StudyFileID[[r]], ms_file_order$StudyFileID)
      xic <- match_xics_pos[cid_idx, file_idx]
    } else {
      # index of compound (row in XIC object)
      cid_idx <- match(cid, match_mz_neg$CID)
      # index of file (column in XIC object)
      file_idx <- match(high_xics$StudyFileID[[r]], ms_file_order$StudyFileID)
      xic <- match_xics_neg[cid_idx, file_idx]
    }
    return(xic)
  })

  ## Calculate Euclidean error norms on scaled XICs ----
  message("Comparing peak shapes...")

  # p-norm calculation assumes that both vectors are non-negative!
  scaled_Lp_error_norm <- function(v1, v2, p) {
    # NA values should actually be 0
    v1[is.na(v1)] <- 0
    v2[is.na(v2)] <- 0
    # scale both vectors to max == 1
    v1 <- v1/max(v1)
    v2 <- v2/max(v2)
    vd <- sum(abs(v2 - v1)^p)^(1/p)/length(v1)^(1/p)
    return(vd)
  }

  xic_cors <- lapply(X = 1:nrow(match_pks_bounds), FUN = function(r) {
    # get the relevant XICs
    xic1 <- filter(high_xics, CID == match_pks_bounds$CID1[r],
                   StudyFileID == match_pks_bounds$StudyFileID[r])$XIC[[1]]
    xic2 <- filter(high_xics, CID == match_pks_bounds$CID2[r],
                   StudyFileID == match_pks_bounds$StudyFileID[r])$XIC[[1]]
    # filter to only the relevant peak
    xic1 <- MSnbase::filterRt(xic1, rt = c(match_pks_bounds$LeftRT[r]*60, match_pks_bounds$RightRT[r]*60))
    xic2 <- MSnbase::filterRt(xic2, rt = c(match_pks_bounds$LeftRT[r]*60, match_pks_bounds$RightRT[r]*60))
    # if(sum(rtime(xic2) != rtime(xic1)) > 0) {
    #   errorCondition("Non matching XICs!")
    # }
    xc <- MSnbase::compareChromatograms(x = xic1, y = xic2, FUN = scaled_Lp_error_norm,
                                        FUNARGS = list(p = 2))
  })


  # # for each compound found in 1+ matches, get its RT and mz windows
  # match_mz1 <- match_pks_bounds %>%
  #   dplyr::group_by(.data$CID1) %>%
  #   dplyr::summarise(Polarity.1 = dplyr::first(.data$Polarity.1),
  #             Mass.1 = dplyr::first(.data$Mass.1),
  #             LeftRT = min(.data$LeftRT)*60,
  #             RightRT = max(.data$RightRT)*60,
  #             # use file where more intense feature is highest
  #             StudyFileID = dplyr::first(.data$StudyFileID)) %>%
  #   dplyr::mutate(mzmin.1 = .data$Mass.1 - xic_tol/1e6*.data$Mass.1,
  #          mzmax.1 = .data$Mass.1 + xic_tol/1e6*.data$Mass.1)
  # match_mz1_pos <- match_mz1 %>%
  #   dplyr::filter(.data$Polarity.1 == 1)
  # match_mz1_neg <- match_mz1 %>%
  #   dplyr::filter(.data$Polarity.1 == -1)
  #
  # match_mz2 <- match_pks_bounds %>%
  #   dplyr::group_by(.data$CID2) %>%
  #   dplyr::summarise(Polarity.2 = dplyr::first(.data$Polarity.2),
  #                    Mass.2 = dplyr::first(.data$Mass.2),
  #                    LeftRT = min(.data$LeftRT)*60,
  #                    RightRT = max(.data$RightRT)*60,
  #                    # use file where more intense feature is highest
  #                    StudyFileID = dplyr::first(.data$StudyFileID)) %>%
  #   dplyr::mutate(mzmin.2 = .data$Mass.2 - xic_tol/1e6*.data$Mass.2,
  #                 mzmax.2 = .data$Mass.2 + xic_tol/1e6*.data$Mass.2)
  # match_mz2_pos <- match_mz2 %>%
  #   dplyr::filter(.data$Polarity.2 == 1)
  # match_mz2_neg <- match_mz2 %>%
  #   dplyr::filter(.data$Polarity.2 == -1)
  #
  # # for each compound in 1+ matches, get XICs from all alignment files
  # # these throw a zillion warnings about zero values, so we suppress those
  # suppressMessages({
  # match_xics1_pos <- xcms::chromatogram(xdata_aln,
  #                                 mz = as.matrix(dplyr::select(match_mz1_pos, "mzmin.1", "mzmax.1")),
  #                                 rt = as.matrix(dplyr::select(match_mz1_pos, "LeftRT", "RightRT")),
  #                                 msLevel = 1, BPPARAM = bp_param,
  #                                 missing = 0)
  # match_xics1_neg <- xcms::chromatogram(xdata_aln_neg,
  #                                 mz = as.matrix(dplyr::select(match_mz1_neg, "mzmin.1", "mzmax.1")),
  #                                 rt = as.matrix(dplyr::select(match_mz1_neg, "LeftRT", "RightRT")),
  #                                 msLevel = 1, BPPARAM = bp_param,
  #                                 missing = 0)
  # match_xics2_pos <- xcms::chromatogram(xdata_aln,
  #                                 mz = as.matrix(dplyr::select(match_mz2_pos, "mzmin.2", "mzmax.2")),
  #                                 rt = as.matrix(dplyr::select(match_mz2_pos, "LeftRT", "RightRT")),
  #                                 msLevel = 1, BPPARAM = bp_param,
  #                                 missing = 0)
  # match_xics2_neg <- xcms::chromatogram(xdata_aln_neg,
  #                                 mz = as.matrix(dplyr::select(match_mz2_neg, "mzmin.2", "mzmax.2")),
  #                                 rt = as.matrix(dplyr::select(match_mz2_neg, "LeftRT", "RightRT")),
  #                                 msLevel = 1, BPPARAM = bp_param,
  #                                 missing = 0)
  # })
  #
  # # Now pull XICs for each match pair from the file where the more abundant
  # # feature is most abundant. This will give us the best peaks for comparing profiles.
  #
  # # make 2-column matrix of (compound, file) indices
  # xic_file_order <- colnames(match_xics1_pos) %>%
  #   stringr::str_sub(end = -6)
  # study_file_key <- msa$input_files %>%
  #   dplyr::select("FileName", "StudyFileID") %>%
  #   dplyr::mutate(File = stringr::str_match(.data$FileName, "[\\\\/]([^\\\\/]*)\\.raw")[, 2],
  #          XIC_Idx = match(.data$File, xic_file_order))
  #
  # # for each compound, get the index of the file where it is most abundant
  # match_file_idx1_pos <- match_mz1_pos %>%
  #   dplyr::left_join(study_file_key, by = "StudyFileID") %>%
  #   dplyr::select("XIC_Idx") %>%
  #   # same row order in match_mz1_pos and match_xics1_pos
  #   dplyr::mutate(Row = dplyr::row_number(), .before = 1) %>%
  #   as.matrix()
  # match_file_idx1_neg <- match_mz1_neg %>%
  #   dplyr::left_join(study_file_key, by = "StudyFileID") %>%
  #   dplyr::select("XIC_Idx") %>%
  #   dplyr::mutate(Row = dplyr::row_number(), .before = 1) %>%
  #   as.matrix()
  # match_file_idx2_pos <- match_mz2_pos %>%
  #   dplyr::left_join(study_file_key, by = "StudyFileID") %>%
  #   dplyr::select("XIC_Idx") %>%
  #   dplyr::mutate(Row = dplyr::row_number(), .before = 1) %>%
  #   as.matrix()
  # match_file_idx2_neg <- match_mz2_neg %>%
  #   dplyr::left_join(study_file_key, by = "StudyFileID") %>%
  #   dplyr::select("XIC_Idx") %>%
  #   dplyr::mutate(Row = dplyr::row_number(), .before = 1) %>%
  #   as.matrix()
  #
  # # for each compound, get the XIC from the file where it is most abundant
  # match_xic_sub1_pos <- lapply(X = 1:nrow(match_file_idx1_pos), FUN = function(r) {
  #   return(match_xics1_pos[match_file_idx1_pos[r, 1], match_file_idx1_pos[r, 2]])
  # })
  # match_xic_sub1_neg <- lapply(X = 1:nrow(match_file_idx1_neg), FUN = function(r) {
  #   return(match_xics1_neg[match_file_idx1_neg[r, 1], match_file_idx1_neg[r, 2]])
  # })
  # match_xic_sub2_pos <- lapply(X = 1:nrow(match_file_idx2_pos), FUN = function(r) {
  #   return(match_xics2_pos[match_file_idx2_pos[r, 1], match_file_idx2_pos[r, 2]])
  # })
  # match_xic_sub2_neg <- lapply(X = 1:nrow(match_file_idx2_neg), FUN = function(r) {
  #   return(match_xics2_neg[match_file_idx2_neg[r, 1], match_file_idx2_neg[r, 2]])
  # })
  #
  # # compare chromatograms for each putative pair
  # cid_xic1_pos <- tibble::tibble(CID = match_mz1_pos$CID1,
  #                        XIC = match_xic_sub1_pos)
  # cid_xic2_pos <- tibble::tibble(CID = match_mz2_pos$CID2,
  #                        XIC = match_xic_sub2_pos)
  # cid_xic1_neg <- tibble::tibble(CID = match_mz1_neg$CID1,
  #                        XIC = match_xic_sub1_neg)
  # cid_xic2_neg <- tibble::tibble(CID = match_mz2_neg$CID2,
  #                        XIC = match_xic_sub2_neg)
  # # bind all XICs for each comparison side into single tibbles
  # cid_xic1 <- dplyr::bind_rows(cid_xic1_neg, cid_xic1_pos)
  # cid_xic2 <- dplyr::bind_rows(cid_xic2_neg, cid_xic2_pos)
  #
  # ## Calculate Euclidean error norms on scaled XICs ----
  # message("Comparing peak shapes...")
  #
  # # p-norm calculation assumes that both vectors are non-negative!
  # scaled_Lp_error_norm <- function(v1, v2, p) {
  #   # NA values should actually be 0
  #   v1[is.na(v1)] <- 0
  #   v2[is.na(v2)] <- 0
  #   # scale both vectors to max == 1
  #   v1 <- v1/max(v1)
  #   v2 <- v2/max(v2)
  #   vd <- sum(abs(v2 - v1)^p)^(1/p)/length(v1)^(1/p)
  #   return(vd)
  # }
  #
  # xic_cors <- lapply(X = 1:nrow(match_pks_bounds), FUN = function(r) {
  #   # get the relevant XICs
  #   xic1 <- cid_xic1$XIC[[match(match_pks_bounds$CID1[r], cid_xic1$CID)]]
  #   xic2 <- cid_xic2$XIC[[match(match_pks_bounds$CID2[r], cid_xic2$CID)]]
  #   # filter to only the relevant peak
  #   xic1 <- MSnbase::filterRt(xic1, rt = c(match_pks_bounds$LeftRT[r]*60, match_pks_bounds$RightRT[r]*60))
  #   xic2 <- MSnbase::filterRt(xic2, rt = c(match_pks_bounds$LeftRT[r]*60, match_pks_bounds$RightRT[r]*60))
  #   # if(sum(rtime(xic2) != rtime(xic1)) > 0) {
  #   #   errorCondition("Non matching XICs!")
  #   # }
  #   xc <- MSnbase::compareChromatograms(x = xic1, y = xic2, FUN = scaled_Lp_error_norm,
  #                              FUNARGS = list(p = 2))
  # })

  df <- aln_compounds %>%
    dplyr::select(CID = "ID", RT = "RetentionTime")

  match_xic_cors <- match_pks_bounds %>%
    dplyr::mutate(XIC_Cor = as.numeric(xic_cors)) %>%
    dplyr::filter(.data$XIC_Cor < xic_error_max) %>%
    dplyr::left_join(df, by = c("CID1" = "CID")) %>%
    dplyr::left_join(df, by = c("CID2" = "CID"), suffix = c(".1", ".2")) %>%
    dplyr::mutate(DeltaRT = abs(.data$RT.1 - .data$RT.2))

  # Construct related compound networks ----
  message("Constructing feature networks...")

  el <- match_xic_cors %>%
    dplyr::select("CID1", "CID2") %>%
    dplyr::mutate(CID1 = as.character(.data$CID1),
           CID2 = as.character(.data$CID2)) %>%
    as.matrix()

  sdf_graph <- igraph::graph_from_edgelist(el = el, directed = FALSE)
  # edges are in the order we fed them in, so we can apply the attributes directly
  # for error norms, we subtract from 1 to make small errors large weights
  igraph::E(sdf_graph)$weight <- 1 - match_xic_cors$XIC_Cor
  igraph::E(sdf_graph)$peak_shape_cor <- 1 - match_xic_cors$XIC_Cor
  igraph::E(sdf_graph)$area_cor <- match_xic_cors$Correlation
  igraph::E(sdf_graph)$rt_diff <- match_xic_cors$DeltaRT

  # vertices are not in an obvious order, but we can replicate it with pivot_longer
  df <- match_xic_cors %>%
    dplyr::mutate(Name.1 = paste(round(.data$Mass.1, digits = 4), round(.data$RT.1, digits = 4), sep = "@"),
           Name.2 = paste(round(.data$Mass.2, digits = 4), round(.data$RT.2, digits = 4), sep = "@")) %>%
    dplyr::select(-dplyr::starts_with("StudyFile"), -dplyr::starts_with("Left"), -dplyr::starts_with("Right"),
           -"Correlation", -"XIC_Cor", -"DeltaRT") %>%
    dplyr::rename(CID.1 = "CID1", CID.2 = "CID2") %>%
    tidyr::pivot_longer(cols = dplyr::everything(), names_sep = "\\.", names_to = c(".value", NA)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(CID = as.character(.data$CID))

  # use sensible labels instead of CIDs
  igraph::V(sdf_graph)$label <- df$Name
  igraph::V(sdf_graph)$cid <- df$CID
  igraph::V(sdf_graph)$RT <- df$RT
  igraph::V(sdf_graph)$coreness <- igraph::coreness(sdf_graph)
  igraph::V(sdf_graph)$strength <- igraph::strength(sdf_graph, weights = igraph::E(sdf_graph)$peak_shape_cor)

  df <- igraph::incident_edges(sdf_graph, v = igraph::V(sdf_graph))
  df <- sapply(X = df, FUN = function(es) {
    return(mean(es$peak_shape_cor))
  })
  igraph::V(sdf_graph)$avg_peak_shape_cor <- df

  df <- igraph::incident_edges(sdf_graph, v = igraph::V(sdf_graph))
  df <- sapply(X = df, FUN = function(es) {
    return(stats::median(es$peak_shape_cor))
  })
  igraph::V(sdf_graph)$med_peak_shape_cor <- df

  ## Calculate central member of each family ----

  # which compounds belong to which clusters?
  sdf_clusters <- tibble::enframe(igraph::components(sdf_graph)$membership,
                          name = "CID", value = "ClusterID") %>%
    dplyr::mutate(CID = as.numeric(.data$CID))

  # calculate central members
  sdf_comps <- igraph::components(sdf_graph)

  sdf_eigen_centers <- sapply(X = seq.int(from = 1, to = sdf_comps$no, by = 1), FUN = function(comp) {
    # get members
    comp_idx <- which(sdf_comps$membership == comp)
    # get component
    comp_graph <- igraph::induced_subgraph(sdf_graph, vids = comp_idx)
    # calculate eigenvector centrality (unweighted, b/c already ~1, other thresholds stringent too)
    comp_eigen_centrality <- igraph::eigen_centrality(comp_graph)
    # pick the most central node
    return(igraph::V(comp_graph)$cid[which.max(comp_eigen_centrality$vector)])
  })

  sdf_eigen_centers <- dplyr::tibble(ClusterID = seq.int(from = 1, to = sdf_comps$no, by = 1),
                                     CenterCID = as.numeric(sdf_eigen_centers))

  message("Done!")

  return(list(matched_compound_pairs = match_xic_cors, sdf_graph = sdf_graph,
              sdf_cluster_members = sdf_clusters,
              sdf_cluster_centers = sdf_eigen_centers))
}
