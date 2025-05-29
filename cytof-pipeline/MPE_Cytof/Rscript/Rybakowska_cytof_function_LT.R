#' extract quantile (0.01, 0.25, 0.5, 0.75, and 0.99) of  
#' marker specified by markers_to_plot or all (markers_to_plot=NULL
# in_dir = getwd(),
extract_marker_quantiles_4SDL <- function(select_fcs_files, 
                                  in_dir = NULL,                                  
                                  batch_pattern,
                                  arcsine_transform = TRUE, 
                                  markers_to_plot = NULL){ 
                                  
  if (is.null(in_dir)){
     fcs_files <- select_fcs_files 
  } else {
     fcs_files <- file.path(in_dir,select_fcs_files)  #*****
  }
   
  if (!file.exists(fcs_files[1])){
    stop("incorrect file path, the fcs file does not exist")
  }
  
  if(check(norm_markers) == 0){
    
    o <- capture.output(ff_tmp <- read.FCS(file.path(fcs_files[1])))
    
    if (!is.null(markers_to_plot)){
      
      if(!is.character(markers_to_plot)){
        stop ("markers are not a character vector")
      }
      
      matches <- paste(markers_to_plot, collapse="|")
      
      norm_markers <- grep(matches,
                FlowSOM::GetMarkers(ff_tmp, find_mass_ch(ff_tmp,value = TRUE)), 
                value = TRUE, ignore.case = F)
    } else {
      norm_markers <- find_mass_ch(ff_tmp, value = TRUE)
      norm_markers <- FlowSOM::GetMarkers(ff_tmp, norm_markers)
    }
  }
  
  quantile_values <-  c(0.01, 0.25, 0.5, 0.75, 0.99)
  quantiles <- expand.grid(File = fcs_files,
                           Marker = norm_markers,
                           Quantile = quantile_values,
                           Value = NA)
  quantiles <- cbind(quantiles, "Batch" = stringr::str_match(
    basename(as.character(quantiles$File)), batch_pattern)[,2])

  for (file in fcs_files) {
    print(file)
    
    o <- capture.output(ff <- read.FCS(file))  
    
    if(arcsine_transform == TRUE){
      ff <- flowCore::transform(ff, transformList(grep("Di", colnames(ff), value = TRUE),
                                        arcsinhTransform(a = 0, b = 1/5, c = 0)))
    }
         
    for (marker in names(norm_markers)) {
      quantiles_res <-stats::quantile(exprs(ff)[, marker],
                                quantile_values)
      for (i in seq_along(quantiles_res)) {
        quantiles <- quantiles %>%
          dplyr::mutate(Value = replace(Value, 
                                        File == file & 
                                          Marker == norm_markers[marker] &
                                          Quantile == quantile_values[i],
                                        quantiles_res[i]))
        
      }
    }
  }

  quantiles$Sample <- gsub(pattern = "Norm_", replacement = "", ignore.case = TRUE,
                             x = gsub(pattern = "_CC_gated.fcs|_gated.fcs|_beadNorm.fcs|.FCS|.fcs",
                                  replacement = "", ignore.case = TRUE,
                                  x = basename(as.character(quantiles$File))))

  
  quantiles$Sample <- gsub(pattern = "_CC_gated.fcs|_gated.fcs|_beadNorm.fcs|.FCS|.fcs",
                                  replacement = "", ignore.case = TRUE,
                                  x = basename(as.character(quantiles$File)))
 
  Ftab_out <- reshape(quantiles, idvar = c("File","Marker","Batch","Sample"),
                      timevar = "Quantile", direction = "wide") 

  return(Ftab_out)
}

## add nCellsPerSample
file_quality_check_4SDL <- function(fcs_files, 
                               file_batch_id = NULL, 
                               out_dir = getwd(), 
                               nCellsPerSample = 1000,
                               phenotyping_markers = NULL, 
                               arcsine_transform = TRUE, 
                               sd = 3, 
                               nClus = 10){
                               
  
  if(!dir.exists(out_dir)) dir.create(out_dir)
  
  if (!is.null(file_batch_id)) {
    scores <- list()
    for (batch in unique(file_batch_id)){
      print(batch)
      
      tmp_files <- fcs_files[file_batch_id == batch]
      nCells <- nCellsPerSample*length(tmp_files)
      fsom <- fsom_aof_4SDL(fcs_files = tmp_files, 
                       phenotyping_markers = phenotyping_markers, 
                       nCells = nCells,
                       out_dir = out_dir, 
                       arcsine_transform = arcsine_transform,
                       nClus = nClus,
                       batch = batch)
    
      scores[[batch]] <- aof_scoring_4SDL(fcs_files = tmp_files, 
                                     phenotyping_markers = phenotyping_markers,
                                     fsom = fsom, out_dir = out_dir, batch = batch)
    }
    
  } else {
    tmp_files <- fcs_files
    fsom <- fsom_aof_4SDL(fcs_files = tmp_files, phenotyping_markers = phenotyping_markers, 
                     out_dir = out_dir, arcsine_transform = arcsine_transform, 
                     nClus = nClus,
                     batch = NULL)
    
    scores <- aof_scoring_4SDL(fcs_files = tmp_files, 
                          phenotyping_markers = phenotyping_markers,
                          fsom = fsom, out_dir = out_dir, batch = NULL)
  }
  
  final_score <- file_outlier_detecion(scores = scores, out_dir = out_dir, 
                                       sd = sd)
}


fsom_aof_4SDL <- function(fcs_files, 
                     phenotyping_markers,
                     nCells = length(fcs_files)*10000,
                     xdim = 10,
                     ydim = 10,
                     nClus = 10,
                     out_dir, 
                     batch = NULL,
                     arcsine_transform = TRUE, 
                     seed = 1){
  
  
  if(check(phenotyping_channels) == 0){
    o <- capture.output(ff_tmp <- flowCore::read.FCS(file.path(fcs_files[1])))  #****error in original code
    markers <- FlowSOM::GetMarkers(ff_tmp, colnames(ff_tmp))
    phenotyping_channels <- grep(paste(phenotyping_markers, 
                                       collapse = ("|")), markers, value = TRUE)
  }
  
  if(arcsine_transform  == TRUE){
    trans <- transformList(names(phenotyping_channels), CytoNorm::cytofTransform)
  } else {
    trans <- NULL
  }
  
  fsom <- CytoNorm::prepareFlowSOM(file = fcs_files,
                                   colsToUse = names(phenotyping_channels),
                                   seed = seed, 
                                   nCells = nCells,
                                   transformList = trans,
                                   FlowSOM.params = list(xdim = xdim, 
                                                         ydim = ydim, 
                                                         nClus = nClus, 
                                                         scale = FALSE))
  
  myCol <- c("tomato", "violet", "grey50", "slateblue1", "yellow","turquoise2",
             "yellowgreen", "skyblue", "wheat2","steelblue", "blue2", "navyblue",
             "orange", "violetred", "red4", "springgreen2",  "palegreen4",
             "tan", "tan2", "tan3", "brown", "grey70", "grey30")
  
  if(max(as.numeric(fsom$metaclustering)) > length(myCol)){
    backgroundColors <- NULL
  } else {
    backgroundColors <- myCol
  }
  
  if(!is.null(batch)){
    filename <- paste0(batch, "_FlowSOM_clustering.pdf")
  } else {
    filename <- "FlowSOM_clustering.pdf"
  }
  
  # pdf(file.path(out_dir, filename), width = 14, height = 10)
  fsomPlot <- FlowSOM::PlotStars(fsom = fsom,
                                 title = "FlowSOM clustering",
                                 backgroundValues = fsom$metaclustering,
                                 maxNodeSize = 3,
                                 backgroundColors = backgroundColors)

  fsomTsne <- FlowSOM::PlotDimRed(fsom = fsom, plotFile = NULL, seed = seed, cTotal = 20000,  
                                  title = "tSNE visualization of FlowSOM metaclusters")
  
  figure <- ggarrange(fsomPlot, fsomTsne,
                      # labels = c("FlowSOM clustering", "tsne"),
                      ncol = 2, nrow = 1)
  
  ggplot2::ggsave(filename = filename, plot = figure, device = "pdf", path = out_dir,
         width =24, height = 10)
  # dev.off()
  
  return(fsom)
}


aof_scoring_4SDL <- function(fcs_files,
                        phenotyping_markers,
                        fsom,
                        out_dir,
                        batch = NULL){

  if(check(phenotyping_channels) == 0){

    o <- capture.output(ff_tmp <- read.FCS(file.path(fcs_files[1])))    #***error in original code
    markers <- FlowSOM::GetMarkers(ff_tmp, colnames(ff_tmp))
    phenotyping_channels <- grep(paste(phenotyping_markers,
                                       collapse = ("|")), markers, value = TRUE)

    if(length(grep("Ir", phenotyping_channels)) > 1){
      phenotyping_channels <- phenotyping_channels[-(grep("Ir",
                                                          phenotyping_channels)[2])]
    }
  }

  aof_scores <- matrix(NA,
                       nrow = length(fcs_files),
                       ncol = length(phenotyping_channels),
                       dimnames = list(fcs_files,
                                       names(phenotyping_channels)))

  for(file in fcs_files){
    print(paste("calculating AOF", file))
    File_ID <- which(fcs_files == file)
    idx <- which(fsom$data[,"File"] == File_ID)
    fcs_data <- fsom$data[idx,]
    MC <- fsom$metaclustering[fsom$map$mapping[idx, 1]]

    aof_tmp <- cytutils::greedyCytometryAof(fcs_data = fcs_data,
                                            y = MC,
                                            channel_names = names(phenotyping_channels),
                                            width = 0.05,
                                            cofactor = 5,
                                            verbose = TRUE)
    aof_scores[file, ] <- aof_tmp$Aof
  }

  scaled_aof_score(aof_scores = aof_scores,
                   out_dir = out_dir,
                   aof_channels = phenotyping_channels,
                   batch = batch)
}


change_fcs_FIL <- function(fcs_files,
                     out_dir = NULL){

  nFiles <- length(fcs_files)

  input.dir <- dirname(fcs_files[1])
  if (is.null(out_dir)){ 
     out_dir <- file.path(input.dir,"updateFileName")
  }

  if (!dir.exists(out_dir)) { 
     dir.create(out_dir)
  }

  for (i in seq_len(nFiles)) {
    fileName <- basename(fcs_files[i])

    f <- flowCore::read.FCS(fcs_files[i])
    f@description$`$FIL` <- fileName
    flowCore::write.FCS(x = f, 
         filename = file.path(out_dir, fileName), 
         endian = "big")
    message(i)
   }
}
   
#-----------------------------------------------
# using lower thres of Ir191 and 193    
# before 4/29/2025
gate_intact_cells_4SDL_old <- function(flow_frame, 
                              file_name = NULL,
                              hard_cutoff = 4,
                              tinypeak_removal1 = 0.8,
                              tinypeak_removal2 = 0.8,
                              alpha1 = 0.05,
                              alpha2 = 0.1, 
                              arcsine_transform = TRUE, ...){
  
  ff <- flow_frame
  
  if (is.null(file_name)){
    file_name <- ff@description$FIL   #****ff@description$GUID.original
  } else {
    file_name 
  }
  
  if(arcsine_transform == TRUE){
    
    ff_t <- flowCore::transform(ff, 
                                flowCore::transformList(colnames(ff)[grep("Di", colnames(ff))], 
                                                        CytoNorm::cytofTransform))
  } else {
    ff_t <- ff
  }
  
  selection <- matrix(TRUE,
                      nrow = nrow(ff),
                      ncol = 1,
                      dimnames = list(NULL,
                                      c("intact")))
  
  tr <- list()
  for(m in c("Ir193Di", "Ir191Di")){
    
    tr[[m]] <- c(flowDensity::deGate(ff_t, m,
                                     tinypeak.removal = tinypeak_removal1, 
                                     upper = FALSE, use.upper = TRUE,
                                     alpha = alpha1, verbose = F, count.lim = 3), 
                 flowDensity::deGate(ff_t, m,
                                     tinypeak.removal = tinypeak_removal2, 
                                     upper = TRUE, use.upper = TRUE,
                                     alpha = alpha2, verbose = F, count.lim = 3)) 
  }
  tr_thres <- tr 

  # ******only using high_thres if low_thres<2 by LT 3/2/25
  ##for(m in c("Ir193Di", "Ir191Di")){
  ##  if (tr[[m]][1]<0.1 | is.na(tr[[m]][1])){  ###*******if low_thres is so low
  ##       selection[ff_t@exprs[,m] < tr[[m]][2], "intact"] <- FALSE
  ##       tr_thres[[m]][3] <- tr[[m]][2]
  ##  }else{
  ##       selection[ff_t@exprs[,m] < tr[[m]][1], "intact"] <- FALSE
  ##       tr_thres[[m]][3] <- tr[[m]][1]
  ##  }
  ##}
  
  for(m in c("Ir193Di", "Ir191Di")){
     if (tr[[m]][1]<0.1 | is.na(tr[[m]][1])){
        if (tr[[m]][2]>hard_cutoff){
           tr_thres[[m]][3]  <- hard_cutoff
           selection[ff_t@exprs[,m] < tr_thres[[m]][3], "intact"] <- FALSE
        } else {
           tr_thres[[m]][3] <- tr[[m]][2]
           selection[ff_t@exprs[,m] < tr[[m]][2], "intact"] <- FALSE
        }
     }else{
         tr_thres[[m]][3] <- tr[[m]][1]
         selection[ff_t@exprs[,m] < tr[[m]][1], "intact"] <- FALSE
      }
   }

 
  percentage <- (sum(selection)/length(selection))*100
  flowDensity::plotDens(ff_t, c("Ir193Di", "Ir191Di"), 
                        main = paste0(basename(file_name)," ( ", format(round(percentage, 2), 
                                                                        nsmall = 2), "% )"))
  
  abline(h = c(tr_thres[["Ir191Di"]]))
  abline(v = c(tr_thres[["Ir193Di"]]))
  points(ff_t@exprs[!selection[,"intact"], c("Ir193Di", "Ir191Di")], pch = ".")

  tr_thres <- as.data.frame(tr_thres)
  tr_thres$info <- c("low","high","selected")
  fout4thres <- gsub(".fcs","_intact_thres.txt",file_name)
  write.table(tr_thres,fout4thres,sep="\t",quote=F,row.names=F)
 
  ff <- ff[selection[,"intact"], ]
  
  return(ff)
}


gate_intact_cells_4SDL <- function(flow_frame, 
                              file_name = NULL,
                              accepted_thres = c(1.5, 4),
                              adjusted_thres = c(2, 3),
                              tinypeak_removal1 = 0.8,
                              tinypeak_removal2 = 0.8,
                              alpha1 = 0.05,
                              alpha2 = 0.1, 
                              arcsine_transform = TRUE, ...){
  
  ff <- flow_frame
  
  if (is.null(file_name)){
    file_name <- ff@description$FIL   #****ff@description$GUID.original
  } else {
    file_name 
  }
  
  if(arcsine_transform == TRUE){
    
    ff_t <- flowCore::transform(ff, 
                                flowCore::transformList(colnames(ff)[grep("Di", colnames(ff))], 
                                                        CytoNorm::cytofTransform))
  } else {
    ff_t <- ff
  }
  
  selection <- matrix(TRUE,
                      nrow = nrow(ff),
                      ncol = 1,
                      dimnames = list(NULL,
                                      c("intact")))
  
  tr <- list()
  
  for(m in c("Ir193Di", "Ir191Di")){
    
    tr[[m]] <- c(flowDensity::deGate(ff_t, m,
                                     tinypeak.removal = tinypeak_removal1, 
                                     upper = FALSE, use.upper = TRUE,
                                     alpha = alpha1, verbose = F, count.lim = 3), 
                 flowDensity::deGate(ff_t, m,
                                     tinypeak.removal = tinypeak_removal2, 
                                     upper = TRUE, use.upper = TRUE,
                                     alpha = alpha2, verbose = F, count.lim = 3)) 
  }
  tr_thres <- tr 
  
  for (m in c("Ir193Di", "Ir191Di")) {

    if (tr[[m]][1] < accepted_thres[1] | is.na(tr[[m]][1])) { # lower end < accepted_thres[1]

      if (tr[[m]][2] < accepted_thres[1]) {                   # upper < accepted_thres[1], completely lower accepted interval

        tr_thres[[m]][3] <- adjusted_thres[1]
        selection[ff_t@exprs[, m] < tr_thres[[m]][3], "intact"] <- FALSE

      } else if (tr[[m]][2] > accepted_thres[2]) {            # upper > accepted_thres[2], i.e. (lower,end) interval >  accepted_thres interval

        tr_thres[[m]][3] <- adjusted_thres[2]
        selection[ff_t@exprs[, m] < tr_thres[[m]][3], "intact"] <- FALSE

      } else {                                                #  upper is inside the accepted_thres interval

        tr_thres[[m]][3] <- tr[[m]][2]
        selection[ff_t@exprs[, m] < tr_thres[[m]][3], "intact"] <- FALSE

      }
    } else if (tr[[m]][1] < accepted_thres[2]) {              # lower is inside the accepted_thres interval
      tr_thres[[m]][3] <- tr[[m]][1]
      selection[ff_t@exprs[, m] < tr[[m]][1], "intact"] <- FALSE
      
    } else {                                                  # lower > accepted_thres[2], i.e. completely outside accepted interval
      tr_thres[[m]][3] <- tr[[m]][1]
      selection[ff_t@exprs[, m] < tr[[m]][1], "intact"] <- FALSE
    }
  }
 
  percentage <- (sum(selection)/length(selection))*100
  flowDensity::plotDens(ff_t, c("Ir193Di", "Ir191Di"), 
                        main = paste0(basename(file_name)," ( ", format(round(percentage, 2), 
                                                                        nsmall = 2), "% )"))
  
  abline(h = c(tr_thres[["Ir191Di"]]))
  abline(v = c(tr_thres[["Ir193Di"]]))
  points(ff_t@exprs[!selection[,"intact"], c("Ir193Di", "Ir191Di")], pch = ".")

  tr_thres <- as.data.frame(tr_thres)
  tr_thres$info <- c("low","high","selected")

  # add file name column
  tr_thres$file_name <- gsub("\\.fcs$", "_gated.fcs", basename(file_name), ignore.case = TRUE)
 
  # filter for rows that have "intact" == TRUE
  ff <- ff[selection[,"intact"], ]

  # output the dataframe for concatenation
  return(list(flowFrame = ff, info_df = tr_thres))
}


# -------------------------------------------------------------------------------
# no Ir191-gating
gate_live_cells_4SDL <- function(flow_frame, 
                            file_name = NULL,
                            viability_channel,
                            tinypeak_removal_viability = 0.8,
                            alpha_viability = 0.1,
                            tinypeak_removal_Iridium = 0.8,
                            alpha_Iridium = 0.05,
                            arcsine_transform = TRUE,
                            accepted_thres = c(3, 6),
                            adjusted_thres = c(4, 6),
                            ... ){
  
  ff <- flow_frame
  
  if (is.null(file_name)){
    file_name <- ff@description$FIL
  } else {
    file_name 
  }
  
  if(arcsine_transform == TRUE){
    
    ff_t <- flowCore::transform(ff, 
                                transformList(colnames(ff)[grep("Di", colnames(ff))], 
                                              CytoNorm::cytofTransform))
  } else {
    ff_t <- ff
  }
  
  selection <- matrix(TRUE,
                      nrow = nrow(ff),
                      ncol = 1,
                      dimnames = list(NULL,
                                      c("live")))
  
  
  v_ch <- grep(viability_channel, colnames(ff), value = T)
  
  tr <- list()
  ###for(m in c("Ir191Di", v_ch)){
  for(m in v_ch){
     upper = TRUE
     alpha = alpha_viability
     tr[[m]] <- flowDensity::deGate(ff_t, m,
                      tinypeak.removal = tinypeak_removal_viability, 
                            upper = upper, use.upper = TRUE,
                            alpha = alpha, verbose = F, count.lim = 3)
  }

  tr_thres <- tr

  # check & adjust viability channel threshold
  for (m in v_ch) {

    if (tr[[m]][1] < accepted_thres[1]) {

      tr_thres[[m]][1] <- adjusted_thres[1]
      selection[ff_t@exprs[, m] > tr_thres[[m]][1], "live"] <- FALSE

    } else if (tr[[m]][1] > accepted_thres[2]) {

      tr_thres[[m]][1] <- adjusted_thres[2]
      selection[ff_t@exprs[, m] > tr_thres[[m]][1], "live"] <- FALSE

    } else {
      tr_thres[[m]][1] <- tr[[m]][1]
      selection[ff_t@exprs[, m] > tr_thres[[m]][1], "live"] <- FALSE
    }
  }

  percentage <- (sum(selection)/length(selection))*100
  flowDensity::plotDens(ff_t, c(v_ch, "Ir191Di"), 
                        main = paste0(file_name," ( ", format(round(percentage, 2), nsmall = 2), "% )"),
                        xlim = c(0, 8), ylim = c(0, 8), ...)
  
  abline(h = tr[["Ir191Di"]])
  # abline(v = tr[[v_ch]])
  abline(v = tr_thres[[v_ch]])
  
  points(ff_t@exprs[!selection[,"live"], c(v_ch, "Ir191Di")], pch = ".") 
  
  ff <- ff[selection[,"live"], ]

  # make dataframe to output
  tr_thres <- as.data.frame(tr_thres)
  tr_thres$file_name <- gsub("\\.fcs$", "_gated.fcs", basename(file_name), ignore.case = TRUE)

  
  return(list(flowFrame = ff, info_df = tr_thres))
  
}



remove_mad_outliers_4SDL <- function(flow_frame, 
                                channels = "Event_length", 
                                n_mad = 2,
                                mad_f = mad,
                                plot = TRUE,
                                center = "center",
                                main = "",
                                event_length_cutoff = 50,
                                ...){
  boundaries <- matrix(NA,
                       nrow = 5,
                       ncol = length(channels),
                       dimnames = list(c("median", "center", "mad", "l_lim", "u_lim"),
                                       channels))
  for (channel in channels) {
    x <- flow_frame@exprs[, channel]
    boundaries["median", channel] <- median(x)
    boundaries["center", channel] <- density(x)$x[which.max(density(x)$y)]
    boundaries["mad", channel] <- mad_f(x,
                                        center = boundaries[center, channel] )
    
    # boundaries["l_lim", channel] <- boundaries[center, channel] - n_mad * boundaries["mad", channel]
    # boundaries["u_lim", channel] <- boundaries[center, channel] + n_mad * boundaries["mad", channel]

    # determine event length boundaries
    l_lim <- boundaries[center, channel] - n_mad * boundaries["mad", channel]
    u_lim <- boundaries[center, channel] + n_mad * boundaries["mad", channel]

    # adjust u_lim to < 50 if original is > 50
    if (u_lim > event_length_cutoff) {
      u_lim <- event_length_cutoff
    }

    boundaries["l_lim", channel] <- l_lim
    boundaries["u_lim", channel] <- u_lim
  }
  
  selection <- rep(TRUE, nrow(flow_frame))
  for (channel in channels) {
    selection <- selection & (flow_frame@exprs[, channel] > boundaries["l_lim", channel])
    selection <- selection & (flow_frame@exprs[, channel] < boundaries["u_lim", channel])
  }
  percentage <- (sum(selection)/length(selection))*100
  if (plot) {
    flowDensity::plotDens(flow_frame, 
                          c(channels, "Ir191Di"), 
                          main = paste0(main, " ( ", format(round(percentage, 2), 
                                                            nsmall = 2), "% )"),
                          ...)
    if(length(channels) == 2) {
      points(flow_frame@exprs[!selection, channels], col = "red", pch = ".")
      abline(v = boundaries[c("l_lim", "u_lim"), channels[1]], col = "grey")
      abline(h = boundaries[c("l_lim", "u_lim"), channels[2]], col = "grey")
    } else if(length(channels) == 1) {
      points(flow_frame@exprs[!selection, c(channels, "Ir191Di")], pch = ".")
      abline(v = boundaries[c("l_lim", "u_lim"), channels[1]], col = "grey")
    }
  }

  # generate dataframe
  df_boundaries <- as.data.frame(t(boundaries))
  colnames(df_boundaries) <- paste0(rep(colnames(boundaries), each = nrow(boundaries)),
                                    "_",
                                    rep(rownames(boundaries), times = ncol(boundaries)))
  # return(selection)
  return(list(selection = selection, info_df = df_boundaries))
}



gate_singlet_cells_4SDL <- function(flow_frame, 
                               channels = "Event_length", 
                               arcsine_transform = TRUE,
                               file_name = NULL,
                               n_mad = 2,
                               event_length_cutoff = 50,
                               ...){
  
  if (is.null(file_name)){
    file_name <- flow_frame@description$FIL
  } else {
    file_name 
  }
  
  if(arcsine_transform == TRUE){
    
    flow_frame_t <- flowCore::transform(flow_frame, 
                                        flowCore::transformList(colnames(flow_frame)[grep("Di", colnames(flow_frame))], 
                                                                CytoNorm::cytofTransform))
  } else {
    flow_frame_t <- flow_frame
  }
  
  selection <- matrix(TRUE,
                      nrow = nrow(flow_frame),
                      ncol = 1,
                      dimnames = list(NULL,
                                      c("singlets")))  
  
  results <- remove_mad_outliers_4SDL(flow_frame = flow_frame_t, 
                                      channels = channels,
                                      main = paste("Singlets", file_name),
                                      n_mad = n_mad,
                                      xlim = c(0, 100),
                                      ylim = c(0, 8),
                                      event_length_cutoff = event_length_cutoff,
                                      ...)

  selection[, "singlets"] <- results$selection

  # extract gating info dataframe & add file name column
  gate_info_df <- results$info_df
  gate_info_df$file_name <- gsub("\\.fcs$", "_gated.fcs", basename(file_name), ignore.case = TRUE)
  
  flow_frame <- flow_frame[selection[,"singlets"], ]
  
  return(list(flowFrame = flow_frame, info_df = gate_info_df))
  
}






# extract sample info from fcs files
getCoreID <- function(x) {
  xlist <- unlist(strsplit(x, split = "_"))[1:4]
  xcore <- paste0(xlist, collapse = "_")
}






#' Title
#'
#' @param iterations 
#'
#' @returns
#' @export
#'
#' @examples
reloadProgressBar <- function(iterations) {
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = iterations, clear = FALSE, width = 60
  )
  return(pb)
}






#' Title
#'
#' @param sce 
#' @param k 
#' @param assay 
#' @param fun 
#' @param scale 
#' @param q 
#' @param pal 
#'
#' @returns
#' @export
#'
#' @examples
dotplot <- function(sce,
                    k,
                    assay = "exprs",
                    fun = "median",
                    scale = TRUE,
                    q = 0.01,
                    pal = hcl.colors(11, "viridis"),
                    output_dir = NULL) {

  sce$cluster_id <- cluster_ids(sce, k)
  es <- assay(sce, assay)   # expression matrix
  th <- rowMedians(es)

  cs <- seq_len(ncol(sce))
  cs <- split(cs, sce$cluster_id)
  fq <- sapply(cs, function(i)
    rowMeans(es[, i, drop = FALSE] > th))

  # compute median expression by cluster
  lab <- paste(fun, assay)
  ms <- CATALYST:::.agg(sce, by = "cluster_id", assay = assay, fun = fun)
  if (scale) {
    lab <- paste("scaled", lab)
    ms <- CATALYST:::.scale_exprs(ms, q = q)
  }

  # do hierarchical clustering on rows & columns
  cluster_order <- function(x) order.dendrogram(as.dendrogram(hclust(dist(x))))
  ro <- colnames(ms)[cluster_order(t(ms))]
  co <- rownames(ms)[cluster_order(ms)]

  ms_ordered <- ms[co, ro]
  fq_ordered <- fq[co, ro]
  df <- cbind(melt(ms_ordered), fq = melt(fq_ordered)$value)

  df_wide <- dcast(df, Var1 ~ Var2, value.var = "value")

  # write.table(df, file = paste0("dotplot_", fun, "_expression_matrix.txt"), sep = "\t", quote = FALSE, col.names = NA)
  # write.table(df_wide, file = paste0("dotplot_", fun, "_expression_matrix_wide.txt"), sep = "\t", quote = FALSE, col.names = NA)

  # Write the *raw* (unclustered) matrix
  df_raw <- cbind(melt(ms), fq = melt(fq)$value)
  df_raw_wide <- dcast(df_raw, Var1 ~ Var2, value.var = "value")

  write.table(df_raw, file = file.path(output_dir, paste0(k, "dotplot_", fun, "_expression_matrix.txt")), sep = "\t", quote = FALSE, col.names = NA)

  write.table(df_raw_wide, file = file.path(output_dir, paste0(k, "dotplot_", fun, "_expression_matrix_wide.txt")), sep = "\t", quote = FALSE, col.names = NA)


  dot_PLOT <- ggplot(df, aes(Var1, Var2, col = value, size = fq, label = sprintf("%.2f", value))) +
    geom_point() +
    # geom_text(color = "black", size = 2.5, vjust = 0.5) +  # show mean/median value; comment out line if desired
    scale_x_discrete("marker", limits = co, expand = c(0, 0.5)) +
    scale_y_discrete("cluster_id", limits = ro, expand = c(0, 0.5)) +
    scale_color_gradientn(lab, breaks = seq(0, 1, 0.5), colors = pal) +
    scale_size_continuous(
      range = c(0, 5),
      labels = formatC(seq(0, 1, 0.25), 2, format = "f")
    ) +
    guides(
      color = guide_colorbar(order = 1),
      size = guide_legend("% cells with expr. above\n global marker median")
    ) +
    coord_equal() +
    theme_linedraw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  return(dot_PLOT)
}




get_cluster_prop_within_sample <- function(sce, k) {
  # get cluster ids & sample ids
  clusters <- cluster_ids(sce, k)
  sample_ids <- colData(sce)$sample_id
  
  # set column names
  cluster_col_name <- paste0(k, "_cluster")
  
  # build dataframe
  df <- as.data.frame(table(sample_id = sample_ids,
                            cluster = clusters))
  
  # rename cluster column
  colnames(df)[2] <- cluster_col_name
  
  # calculate cluster proportion per sample
  df <- df %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(proportion = Freq / sum(Freq)) %>%
    dplyr::ungroup()
  
  return(df)
}



get_cluster_mapping <- function(sce, meta_vec = c("meta8", "meta10", "meta20")) {
  # Validate inputs
  if (!"cluster_id" %in% names(colData(sce))) {
    stop("The SCE object must contain a 'cluster_id' column in colData.")
  }
  if (is.null(metadata(sce)$cluster_codes)) {
    stop("The SCE object does not contain cluster_codes in metadata.")
  }
  
  cluster_codes <- metadata(sce)$cluster_codes
  
  # Ensure all meta columns exist
  missing_metas <- setdiff(meta_vec, colnames(cluster_codes))
  if (length(missing_metas) > 0) {
    stop(paste("The following metaclusters do not exist:", paste(missing_metas, collapse = ", ")))
  }
  
  # Extract columns
  mapping <- cluster_codes[, c("som100", meta_vec), drop = FALSE]
  colnames(mapping)[1] <- "cluster_id"
  
  # Sort by cluster_id
  mapping$cluster_id <- as.integer(as.character(mapping$cluster_id))
  mapping <- mapping[order(mapping$cluster_id), ]
  
  return(as.data.frame(mapping))
}


























