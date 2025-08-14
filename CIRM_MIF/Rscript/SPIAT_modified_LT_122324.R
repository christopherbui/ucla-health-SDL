## ******************************************************************************************************
## Convert halo output files to single cell object
## adapt from SPIAT (https://www.biorxiv.org/content/10.1101/2020.05.28.122614v1)
## orginial function is format_image_to_sce
## by LT 09/13/21 
## input files are:
## 1. image_fpath: full path to the halo ouput csv file
## 2. channelInfo_fpath: full path to the csv file: 
##       column 1 = column names of intensity (i.e. xxx Cell Intensity) --> continuous expression
##       column 2 = column names of binary value based on intensity (i.e. xxxx Positive Classification) --> 1,0
##       column 3 = markers, including DAPI 
##       First row is header
##       Oder in each row must be identical to halo ouput, 
##                  e.g. Opal 480, Opal 520, etc. 
##       --> copied from halo header file
## 3. HaloDefinedPhenotype_fpath: full path to the csv file listing column names in halo ouput
##       recording phenotypes in column1 and modified names in Col2--> NO HEADER
##       It is an option.
## 4. Prefix: sample/slide ID
## Differeing from original code:
## 1. no removal of cells with DAPI intensity=0
## 2. add prefix to cell/object ID
## 3. additional information based on HaloPhenotype
## 4. add id for ROI
## 5. add cell area in (micron unit)
## *****************************************************************************************************************
format_halo_to_sce_4DL <- function(image_fpath, channelInfo_fpath, haloPhenotype_fpath=NULL, prefixID="Cell") {

    # process the data based on data format
    # read in the image
    image <- read.csv(image_fpath,stringsAsFactors=FALSE,header=T)

    # read channel/marker information
    channelInfo <- read.csv(channelInfo_fpath,stringsAsFactors=FALSE,header=T)

    # import column names recording phenotypes defined by halo - NO HEADER
    if (!is.null(haloPhenotype_fpath)){
       haloPhenotype <- as.matrix(read.csv(haloPhenotype_fpath,stringsAsFactors=FALSE,header=F))      
    }

    # replace the spaces and non-alphanumeric characters as a '.' for column selection
    intensity_columns_interest <- gsub("[^[:alnum:]]", ".", channelInfo[,1])  ##**continuous variables
    dye_columns_interest <- gsub("[^[:alnum:]]", ".", channelInfo[,2])        ##**binary variables  
        
    # CHECK - if image contains all the columns specified and vectors of same length
    image_colnames <- colnames(image)
    if (!all(intensity_columns_interest %in% image_colnames)) {
        stop("One or more Intensity_columns_interest not found in image")
    }
    if (!all(dye_columns_interest %in% image_colnames)) {
        stop("One or more dye_columns_interest not found in image")
    }

    # extract ROIs 
    roi_column_interest <- c("Analysis.Region")
    if (is.element(roi_column_interest,image_colnames)){
        roi_info <- image[,roi_column_interest]
    } else {
        stop("roi_column_interest not found in image")
    }

    # extract Area
    area_column_interest <- c("Cell.Area")
    area_column_idx <- grepl(area_column_interest,image_colnames)
    if (sum(area_column_idx)==1){
        area_data <- image[,area_column_idx==TRUE]
    } else {
        stop("area_column_interest not found in image")
    }

    
    markers <- channelInfo[,3]
    markers <- gsub(" ","",markers)  ##remove white space

    # CHECK - if DAPI is in the dyes
    idx <- which(markers == "DAPI")
    if (length(idx) == 0) {
        stop("Please include DAPI in the markers")
    }       
        
    # extract intensities
    intensity_of_markers <- image[,intensity_columns_interest]
    colnames(intensity_of_markers) <- markers
    intensity_of_markers[intensity_of_markers == "#N/A"] <- NA
    ## skip below code since import data with stringsAsFactors=FALSE
    ##intensity_of_markers <- apply(intensity_of_markers, 2, function(x){   
    ##        as.numeric(as.character(x))
    ##})

    # get the expression status columns
    expression_status_cols <- image[,dye_columns_interest]
    colnames(expression_status_cols) <- markers

    # start reading in the Phenotypes of every cell based on binary status 
    # ****Phenotype as concatenate of positive markers******
    expression_status_cols$Phenotype <- ""
    for (marker in markers) {
        if (marker == "DAPI") {
            phenotype <- "OTHER,"
        } else {
            phenotype <- paste(marker, ",", sep = "")
        }

        # get the row idx of the cells that express the specific marker, and paste the phenotype
        rows_true_exp <- which(expression_status_cols[,marker] != 0)
        if (length(rows_true_exp) != 0) {
            expression_status_cols[rows_true_exp,]$Phenotype <- paste(expression_status_cols[rows_true_exp,]$Phenotype, 
                    phenotype, sep="")
        }
    }

    # now clean the phenotype column - remove DAPI ("OTHER") marker in the concatenated list
    if (nrow(expression_status_cols[expression_status_cols$Phenotype == "", ]) != 0) {
        expression_status_cols[expression_status_cols$Phenotype == "", ]$Phenotype <- "Negative_all"
    }
    if (nrow(expression_status_cols[expression_status_cols$Phenotype == "OTHER,", ]) != 0) {
        expression_status_cols[expression_status_cols$Phenotype == "OTHER,", ]$Phenotype <- "OTHER"
    }
    expression_status_cols$Phenotype <- gsub("OTHER,", "", expression_status_cols$Phenotype)
    expression_status_cols$Phenotype <- gsub(",OTHER", "", expression_status_cols$Phenotype)
    expression_status_cols$Phenotype <- gsub(",$", "", expression_status_cols$Phenotype)

    # start reading in the HALO Phenotypes of every cell based on specified column
    # *********Define by HALO users*************************************
    if (!is.null(haloPhenotype_fpath)){
        ctype_columns_interest <- gsub("[^[:alnum:]]", ".", haloPhenotype[,1])

        if (!all(ctype_columns_interest %in% image_colnames)) {
            stop("One or more ctype_columns_interest not found in image")
        }
        
        # get halo defined Phenotype
        ctype_status_cols <- image[,ctype_columns_interest]
        if (ncol(haloPhenotype)==2){
            colnames(ctype_status_cols) <- haloPhenotype[,2]  ###sapply(haloPhenotype,function(x) gsub(" ","",x)) 
        }else{
            colnames(ctype_status_cols) <- apply(haloPhenotype[,1],function(x) gsub(" ","",x))
        }
        cellTypes <- colnames(ctype_status_cols)
        tmp <- apply(ctype_status_cols[,1:length(cellTypes)],1,
                function(x) {paste0(cellTypes[x==1],collapse="_",sep="")})
        ctype_status_cols$HaloPhenotype <- as.character(factor(tmp))
        ctype_status_cols[ctype_status_cols$HaloPhenotype == "", ]$HaloPhenotype <- "NotSpecified"
    }

    # grab relavant columns
    image <- image[,c("Object.Id", "XMin", "XMax", "YMin", "YMax")]

    # rename Object.ID to Cell.ID
    colnames(image)[1] <- "Cell.ID"

    # add roi, "Cell_" in front of Cell.ID
    image$Cell.ID <- paste(roi_info,"_",prefixID,image$Cell.ID, sep="")
    image$Loc.ID <- apply(image[,c("XMin", "XMax", "YMin", "YMax")],1,
                       function(x) paste0(x,collapse="_"))

    # add averaged X and Y position
    image$Cell.X.Position <- (image$XMin + image$XMax)/2
    image$Cell.Y.Position <- (image$YMin + image$YMax)/2

    # grab the phenotype column and cbind to image
    phenotype_column <- data.frame(expression_status_cols$Phenotype)
    colnames(phenotype_column) <- "Phenotype"
    image <- cbind(image, phenotype_column)
    image$Phenotype <- as.character(image$Phenotype)
    image$ROI <- as.character(roi_info)
    image$CellArea <- area_data
    image <- image[,c("Cell.ID", "Loc.ID", "ROI", "Phenotype", "Cell.X.Position", "Cell.Y.Position","CellArea")]
    if (!is.null(haloPhenotype_fpath)){
        ##Halophenotype_column <- ctype_status_cols$HaloPhenotype   ####data.frame(ctype_status_cols$HaloPhenotype)
        ##colnames(Halophenotype_column) <- "HaloPhenotype" 
        image$HaloPhenotype <- ctype_status_cols$HaloPhenotype   ###Halophenotype_column       
    }
        
    # create the formatted_data with intensity levels
    formatted_data <- cbind(image, intensity_of_markers)

    ## *******DO NOT REMOVE CELLS WTH DAPI INTENSITY=0 YET   
    ##idx <- which(markers == "DAPI")
    ##DAPI_col_name <- dye_columns_interest[idx]
    ##DAPI_non_zero_rows <- which(image[,DAPI_col_name] != 0)
    ##formatted_data <- formatted_dat
    # now create the SCE object...
    # grab the expression level, markers and cell IDs
    assay_data <- formatted_data[,markers]
    assay_rownames <- markers
    assay_colnames <- formatted_data[,"Cell.ID"]

    # transpose the matrix so every column is a cell and every row is a marker
    assay_data_matrix <- as.matrix(assay_data)
    colnames(assay_data_matrix) <- NULL
    rownames(assay_data_matrix) <- NULL
    assay_data_matrix_t <- t(assay_data_matrix)

    sce <- SingleCellExperiment(assays = list(counts = assay_data_matrix_t))

    rownames(sce) <- assay_rownames
    colnames(sce) <- assay_colnames

    # Assign the phenotype, X and Y positions as the colData
    coldata_phenotype <- formatted_data[,"Phenotype"]    
    coldata_Xpos <- formatted_data[,"Cell.X.Position"]
    coldata_Ypos <- formatted_data[,"Cell.Y.Position"]
    colData(sce)$Cell.X.Position <- coldata_Xpos
    colData(sce)$Cell.Y.Position <- coldata_Ypos
    colData(sce)$Phenotype <- coldata_phenotype
    if (!is.null(haloPhenotype_fpath)){
       coldata_halophenotype <- formatted_data[,"HaloPhenotype"]
       colData(sce)$HaloPhenotype<- coldata_halophenotype
    }
    colData(sce)$ROI <- formatted_data[,"ROI"]
    colData(sce)$CellArea <- formatted_data[,"CellArea"]
    colData(sce)$Loc.ID <- formatted_data[,"Loc.ID"]
    return(sce)
}

# ***********************************************************************
# filter DAPI negative cells
format_halo_to_sce_4DL_DAPIfilter <- function(image_fpath, channelInfo_fpath, haloPhenotype_fpath=NULL, prefixID="Cell") {

    # process the data based on data format
    # read in the image
    image <- read.csv(image_fpath,stringsAsFactors=FALSE,header=T)

    # read channel/marker information
    channelInfo <- read.csv(channelInfo_fpath,stringsAsFactors=FALSE,header=T)

    # import column names recording phenotypes defined by halo - NO HEADER
    if (!is.null(haloPhenotype_fpath)){
       haloPhenotype <- as.matrix(read.csv(haloPhenotype_fpath,stringsAsFactors=FALSE,header=F))      
    }

    # replace the spaces and non-alphanumeric characters as a '.' for column selection
    intensity_columns_interest <- gsub("[^[:alnum:]]", ".", channelInfo[,1])  ##**continuous variables
    dye_columns_interest <- gsub("[^[:alnum:]]", ".", channelInfo[,2])        ##**binary variables  
        
    # CHECK - if image contains all the columns specified and vectors of same length
    image_colnames <- colnames(image)
    if (!all(intensity_columns_interest %in% image_colnames)) {
        stop("One or more Intensity_columns_interest not found in image")
    }
    if (!all(dye_columns_interest %in% image_colnames)) {
        stop("One or more dye_columns_interest not found in image")
    }

    # extract ROIs 
    roi_column_interest <- c("Analysis.Region")
    if (is.element(roi_column_interest,image_colnames)){
        roi_info <- image[,roi_column_interest]
    } else {
        stop("roi_column_interest not found in image")
    }

    # extract Area
    area_column_interest <- c("Cell.Area")
    area_column_idx <- grepl(area_column_interest,image_colnames)
    if (sum(area_column_idx)==1){
        area_data <- image[,area_column_idx==TRUE]
    } else {
        stop("area_column_interest not found in image")
    }

    
    markers <- channelInfo[,3]
    markers <- gsub(" ","",markers)  ##remove white space

    # CHECK - if DAPI is in the dyes
    idx <- which(markers == "DAPI")
    if (length(idx) == 0) {
        stop("Please include DAPI in the markers")
    }       

    # ***add in DAPI filter file
    idx <- which(markers == "DAPI")
    DAPI_col_name <- dye_columns_interest[idx]
    DAPI_non_zero_rows <- which(image[,DAPI_col_name] != 0)
    # ---------
    
    # extract intensities
    intensity_of_markers <- image[,intensity_columns_interest]
    colnames(intensity_of_markers) <- markers
    intensity_of_markers[intensity_of_markers == "#N/A"] <- NA
    ## skip below code since import data with stringsAsFactors=FALSE
    ##intensity_of_markers <- apply(intensity_of_markers, 2, function(x){   
    ##        as.numeric(as.character(x))
    ##})

    # get the expression status columns
    expression_status_cols <- image[,dye_columns_interest]
    colnames(expression_status_cols) <- markers

    # start reading in the Phenotypes of every cell based on binary status 
    # ****Phenotype as concatenate of positive markers******
    expression_status_cols$Phenotype <- ""
    for (marker in markers) {
        if (marker == "DAPI") {
            phenotype <- "OTHER,"
        } else {
            phenotype <- paste(marker, ",", sep = "")
        }

        # get the row idx of the cells that express the specific marker, and paste the phenotype
        rows_true_exp <- which(expression_status_cols[,marker] != 0)
        if (length(rows_true_exp) != 0) {
            expression_status_cols[rows_true_exp,]$Phenotype <- paste(expression_status_cols[rows_true_exp,]$Phenotype, 
                    phenotype, sep="")
        }
    }

    # now clean the phenotype column - remove DAPI ("OTHER") marker in the concatenated list
    if (nrow(expression_status_cols[expression_status_cols$Phenotype == "", ]) != 0) {
        expression_status_cols[expression_status_cols$Phenotype == "", ]$Phenotype <- "Negative_all"
    }
    if (nrow(expression_status_cols[expression_status_cols$Phenotype == "OTHER,", ]) != 0) {
        expression_status_cols[expression_status_cols$Phenotype == "OTHER,", ]$Phenotype <- "OTHER"
    }
    expression_status_cols$Phenotype <- gsub("OTHER,", "", expression_status_cols$Phenotype)
    expression_status_cols$Phenotype <- gsub(",OTHER", "", expression_status_cols$Phenotype)
    expression_status_cols$Phenotype <- gsub(",$", "", expression_status_cols$Phenotype)

    # start reading in the HALO Phenotypes of every cell based on specified column
    # *********Define by HALO users*************************************
    if (!is.null(haloPhenotype_fpath)){
        ctype_columns_interest <- gsub("[^[:alnum:]]", ".", haloPhenotype[,1])

        if (!all(ctype_columns_interest %in% image_colnames)) {
            stop("One or more ctype_columns_interest not found in image")
        }
        
        # get halo defined Phenotype
        ctype_status_cols <- image[,ctype_columns_interest]
        if (ncol(haloPhenotype)==2){
            colnames(ctype_status_cols) <- haloPhenotype[,2]  ###sapply(haloPhenotype,function(x) gsub(" ","",x)) 
        }else{
            colnames(ctype_status_cols) <- apply(haloPhenotype[,1],function(x) gsub(" ","",x))
        }
        cellTypes <- colnames(ctype_status_cols)
        tmp <- apply(ctype_status_cols[,1:length(cellTypes)],1,
                function(x) {paste0(cellTypes[x==1],collapse="_",sep="")})
        ctype_status_cols$HaloPhenotype <- as.character(factor(tmp))
        ctype_status_cols[ctype_status_cols$HaloPhenotype == "", ]$HaloPhenotype <- "NotSpecified"
    }

    # grab relavant columns
    image <- image[,c("Object.Id", "XMin", "XMax", "YMin", "YMax")]

    # rename Object.ID to Cell.ID
    colnames(image)[1] <- "Cell.ID"

    # add roi, "Cell_" in front of Cell.ID
    image$Cell.ID <- paste(roi_info,"_",prefixID,image$Cell.ID, sep="")
    image$Loc.ID <- apply(image[,c("XMin", "XMax", "YMin", "YMax")],1,
                       function(x) paste0(x,collapse="_"))

    # add averaged X and Y position
    image$Cell.X.Position <- (image$XMin + image$XMax)/2
    image$Cell.Y.Position <- (image$YMin + image$YMax)/2

    # grab the phenotype column and cbind to image
    phenotype_column <- data.frame(expression_status_cols$Phenotype)
    colnames(phenotype_column) <- "Phenotype"
    image <- cbind(image, phenotype_column)
    image$Phenotype <- as.character(image$Phenotype)
    image$ROI <- as.character(roi_info)
    image$CellArea <- area_data
    image <- image[,c("Cell.ID", "Loc.ID", "ROI", "Phenotype", "Cell.X.Position", "Cell.Y.Position","CellArea")]
    if (!is.null(haloPhenotype_fpath)){
        ##Halophenotype_column <- ctype_status_cols$HaloPhenotype   ####data.frame(ctype_status_cols$HaloPhenotype)
        ##colnames(Halophenotype_column) <- "HaloPhenotype" 
        image$HaloPhenotype <- ctype_status_cols$HaloPhenotype   ###Halophenotype_column       
    }
        
    # create the formatted_data with intensity levels
    formatted_data <- cbind(image, intensity_of_markers)

    ## *******REMOVE CELLS WTH DAPI INTENSITY=0 YET   
    formatted_data <- formatted_data[DAPI_non_zero_rows,]

    # now create the SCE object...
    # grab the expression level, markers and cell IDs
    assay_data <- formatted_data[,markers]
    assay_rownames <- markers
    assay_colnames <- formatted_data[,"Cell.ID"]

    # transpose the matrix so every column is a cell and every row is a marker
    assay_data_matrix <- as.matrix(assay_data)
    colnames(assay_data_matrix) <- NULL
    rownames(assay_data_matrix) <- NULL
    assay_data_matrix_t <- t(assay_data_matrix)

    sce <- SingleCellExperiment(assays = list(counts = assay_data_matrix_t))

    rownames(sce) <- assay_rownames
    colnames(sce) <- assay_colnames

    # Assign the phenotype, X and Y positions as the colData
    coldata_phenotype <- formatted_data[,"Phenotype"]    
    coldata_Xpos <- formatted_data[,"Cell.X.Position"]
    coldata_Ypos <- formatted_data[,"Cell.Y.Position"]
    colData(sce)$Cell.X.Position <- coldata_Xpos
    colData(sce)$Cell.Y.Position <- coldata_Ypos
    colData(sce)$Phenotype <- coldata_phenotype
    if (!is.null(haloPhenotype_fpath)){
       coldata_halophenotype <- formatted_data[,"HaloPhenotype"]
       colData(sce)$HaloPhenotype<- coldata_halophenotype
    }
    colData(sce)$ROI <- formatted_data[,"ROI"]
    colData(sce)$CellArea <- formatted_data[,"CellArea"]
    colData(sce)$Loc.ID <- formatted_data[,"Loc.ID"]
    return(sce)
}

# ***********************************************************************
# filter DAPI negative cells and include tissue clasification information if exist
format_haloTissueClass_to_sce_4DL_DAPIfilter <- function(image_fpath, channelInfo_fpath, haloPhenotype_fpath=NULL, prefixID="Cell") {

    # process the data based on data format
    # read in the image
    image <- read.csv(image_fpath,stringsAsFactors=FALSE,header=T)

    # read channel/marker information
    channelInfo <- read.csv(channelInfo_fpath,stringsAsFactors=FALSE,header=T)

    # import column names recording phenotypes defined by halo - NO HEADER
    if (!is.null(haloPhenotype_fpath)){
       haloPhenotype <- as.matrix(read.csv(haloPhenotype_fpath,stringsAsFactors=FALSE,header=F))      
    }

    # replace the spaces and non-alphanumeric characters as a '.' for column selection
    intensity_columns_interest <- gsub("[^[:alnum:]]", ".", channelInfo[,1])  ##**continuous variables
    dye_columns_interest <- gsub("[^[:alnum:]]", ".", channelInfo[,2])        ##**binary variables  
        
    # CHECK - if image contains all the columns specified and vectors of same length
    image_colnames <- colnames(image)
    if (!all(intensity_columns_interest %in% image_colnames)) {
        stop("One or more Intensity_columns_interest not found in image")
    }
    if (!all(dye_columns_interest %in% image_colnames)) {
        stop("One or more dye_columns_interest not found in image")
    }

    # extract ROIs 
    roi_column_interest <- c("Analysis.Region")
    if (is.element(roi_column_interest,image_colnames)){
        roi_info <- image[,roi_column_interest]
    } else {
        stop("roi_column_interest not found in image")
    }

    # extract Area
    area_column_interest <- c("Cell.Area")
    area_column_idx <- grepl(area_column_interest,image_colnames)
    if (sum(area_column_idx)==1){
        area_data <- image[,area_column_idx==TRUE]
    } else {
        stop("area_column_interest not found in image")
    }

    # extract Classifier.Label (tissue vs. glass vs. stromal), added on 5/1/2024
    classLabel_column_interest <- c("Classifier.Label")
    classLabel_column_idx <- grepl(classLabel_column_interest,image_colnames)
    if (sum(area_column_idx)==1){
        classLabel_data <- image[,classLabel_column_idx==TRUE]
    } else {
        stop("Classifier.Label column not found in image")
    }
    
    markers <- channelInfo[,3]
    markers <- gsub(" ","",markers)  ##remove white space

    # CHECK - if DAPI is in the dyes
    idx <- which(markers == "DAPI")
    if (length(idx) == 0) {
        stop("Please include DAPI in the markers")
    }       

    # ***add in DAPI filter file
    idx <- which(markers == "DAPI")
    DAPI_col_name <- dye_columns_interest[idx]
    DAPI_non_zero_rows <- which(image[,DAPI_col_name] != 0)
    # ---------
    
    # extract intensities
    intensity_of_markers <- image[,intensity_columns_interest]
    colnames(intensity_of_markers) <- markers
    intensity_of_markers[intensity_of_markers == "#N/A"] <- NA
    ## skip below code since import data with stringsAsFactors=FALSE
    ##intensity_of_markers <- apply(intensity_of_markers, 2, function(x){   
    ##        as.numeric(as.character(x))
    ##})

    # get the expression status columns
    expression_status_cols <- image[,dye_columns_interest]
    colnames(expression_status_cols) <- markers

    # start reading in the Phenotypes of every cell based on binary status 
    # ****Phenotype as concatenate of positive markers******
    expression_status_cols$Phenotype <- ""
    for (marker in markers) {
        if (marker == "DAPI") {
            phenotype <- "OTHER,"
        } else {
            phenotype <- paste(marker, ",", sep = "")
        }

        # get the row idx of the cells that express the specific marker, and paste the phenotype
        rows_true_exp <- which(expression_status_cols[,marker] != 0)
        if (length(rows_true_exp) != 0) {
            expression_status_cols[rows_true_exp,]$Phenotype <- paste(expression_status_cols[rows_true_exp,]$Phenotype, 
                    phenotype, sep="")
        }
    }

    # now clean the phenotype column - remove DAPI ("OTHER") marker in the concatenated list
    if (nrow(expression_status_cols[expression_status_cols$Phenotype == "", ]) != 0) {
        expression_status_cols[expression_status_cols$Phenotype == "", ]$Phenotype <- "Negative_all"
    }
    if (nrow(expression_status_cols[expression_status_cols$Phenotype == "OTHER,", ]) != 0) {
        expression_status_cols[expression_status_cols$Phenotype == "OTHER,", ]$Phenotype <- "OTHER"
    }
    expression_status_cols$Phenotype <- gsub("OTHER,", "", expression_status_cols$Phenotype)
    expression_status_cols$Phenotype <- gsub(",OTHER", "", expression_status_cols$Phenotype)
    expression_status_cols$Phenotype <- gsub(",$", "", expression_status_cols$Phenotype)

    # start reading in the HALO Phenotypes of every cell based on specified column
    # *********Define by HALO users*************************************
    if (!is.null(haloPhenotype_fpath)){
        ctype_columns_interest <- gsub("[^[:alnum:]]", ".", haloPhenotype[,1])

        if (!all(ctype_columns_interest %in% image_colnames)) {
            stop("One or more ctype_columns_interest not found in image")
        }
        
        # get halo defined Phenotype
        ctype_status_cols <- image[,ctype_columns_interest]
        if (ncol(haloPhenotype)==2){
            colnames(ctype_status_cols) <- haloPhenotype[,2]  ###sapply(haloPhenotype,function(x) gsub(" ","",x)) 
        }else{
            colnames(ctype_status_cols) <- apply(haloPhenotype[,1],function(x) gsub(" ","",x))
        }
        cellTypes <- colnames(ctype_status_cols)
        tmp <- apply(ctype_status_cols[,1:length(cellTypes)],1,
                function(x) {paste0(cellTypes[x==1],collapse="_",sep="")})
        ctype_status_cols$HaloPhenotype <- as.character(factor(tmp))
        ctype_status_cols[ctype_status_cols$HaloPhenotype == "", ]$HaloPhenotype <- "NotSpecified"
    }

    # grab relavant columns
    image <- image[,c("Object.Id", "XMin", "XMax", "YMin", "YMax")]

    # rename Object.ID to Cell.ID
    colnames(image)[1] <- "Cell.ID"

    # add roi, "Cell_" in front of Cell.ID
    image$Cell.ID <- paste(roi_info,"_",prefixID,image$Cell.ID, sep="")
    image$Loc.ID <- apply(image[,c("XMin", "XMax", "YMin", "YMax")],1,
                       function(x) paste0(x,collapse="_"))

    # add averaged X and Y position
    image$Cell.X.Position <- (image$XMin + image$XMax)/2
    image$Cell.Y.Position <- (image$YMin + image$YMax)/2

    # grab the phenotype column and cbind to image
    phenotype_column <- data.frame(expression_status_cols$Phenotype)
    colnames(phenotype_column) <- "Phenotype"
    image <- cbind(image, phenotype_column)
    image$Phenotype <- as.character(image$Phenotype)
    image$ROI <- as.character(roi_info)
    image$CellArea <- area_data
    image$ClassifierLabel <- classLabel_data 
    image <- image[,c("Cell.ID", "Loc.ID", "ROI", "Phenotype", "Cell.X.Position", "Cell.Y.Position","CellArea","ClassifierLabel")]  ##mod on 5/1/24
    if (!is.null(haloPhenotype_fpath)){
        ##Halophenotype_column <- ctype_status_cols$HaloPhenotype   ####data.frame(ctype_status_cols$HaloPhenotype)
        ##colnames(Halophenotype_column) <- "HaloPhenotype" 
        image$HaloPhenotype <- ctype_status_cols$HaloPhenotype   ###Halophenotype_column       
    }
        
    # create the formatted_data with intensity levels
    formatted_data <- cbind(image, intensity_of_markers)

    ## *******REMOVE CELLS WTH DAPI INTENSITY=0 YET   
    formatted_data <- formatted_data[DAPI_non_zero_rows,]

    # now create the SCE object...
    # grab the expression level, markers and cell IDs
    assay_data <- formatted_data[,markers]
    assay_rownames <- markers
    assay_colnames <- formatted_data[,"Cell.ID"]

    # transpose the matrix so every column is a cell and every row is a marker
    assay_data_matrix <- as.matrix(assay_data)
    colnames(assay_data_matrix) <- NULL
    rownames(assay_data_matrix) <- NULL
    assay_data_matrix_t <- t(assay_data_matrix)

    sce <- SingleCellExperiment(assays = list(counts = assay_data_matrix_t))

    rownames(sce) <- assay_rownames
    colnames(sce) <- assay_colnames

    # Assign the phenotype, X and Y positions as the colData
    coldata_phenotype <- formatted_data[,"Phenotype"]    
    coldata_Xpos <- formatted_data[,"Cell.X.Position"]
    coldata_Ypos <- formatted_data[,"Cell.Y.Position"]
    colData(sce)$Cell.X.Position <- coldata_Xpos
    colData(sce)$Cell.Y.Position <- coldata_Ypos
    colData(sce)$Phenotype <- coldata_phenotype
    if (!is.null(haloPhenotype_fpath)){
       coldata_halophenotype <- formatted_data[,"HaloPhenotype"]
       colData(sce)$HaloPhenotype<- coldata_halophenotype
    }
    colData(sce)$ROI <- formatted_data[,"ROI"]
    colData(sce)$CellArea <- formatted_data[,"CellArea"]
    colData(sce)$Loc.ID <- formatted_data[,"Loc.ID"]
    colData(sce)$ClassifierLabel <-  formatted_data[,"ClassifierLabel"]  #added on 5/1/24
    return(sce)
}

# ******************************************************************************************
# Re-order list of markers in case the marker-fluorophore relationship was altered
reorder_markerList <- function(old_list){
   new_list <- old_list
   for (i in c(1:length(old_list))){
      tmp_ele <- sort(unlist(strsplit(old_list[i],split=",")))
      new_list[i] <- paste0(tmp_ele,collapse=",")
   }
   tmp_out <- cbind(old_list,new_list)
   colnames(tmp_out) <- c("old_ored","new_order")
   return(tmp_out)
}

# ****************************************************************************************************
# Import Halo output of TIL analysis (bands/ring from tumor border)
# Iputs:
# 1. Halo output of immune filtration analysis
# 2. colnames associated with cell type ID: 
#        col 1: column names in Halo file (copy and paste)
#        col 2: shortnames of cell types (i.e. CD8, tumor), must be IDENTICAL to those in 
#               Halo defined phenotype in single-cell data
# 3. PrefixID: must be IDENTICAL to those speficified when importing single-cell data
# 9/13/21 by LT
# ****************************************************************************************************
import_bandInfo <- function(tilArea_fpath,tilArea_var,prefixID="Cell"){

    # import data
    data_in <- read.csv(tilArea_fpath,stringsAsFactors=FALSE,header=T)
    #cellTypeInfo <- as.matrix(read.csv(tilArea_var,stringsAsFactors=FALSE,header=T))

    ##til_colnames <- colnames(data_in)
 
    ##ctype_columns_interest <- gsub("[^[:alnum:]]", ".", cellTypeInfo[,1])
    ##cellTypes <- gsub(" ","",cellTypeInfo[,2])	

    ##if (!all(ctype_columns_interest %in% til_colnames)) {
    ##     stop("One or more ctype_columns_interest not found in til")
    ##}

    # extract cell IDs and types (specified in single-cell data)
    ##ctype_mat <- data_in[,ctype_columns_interest]
    ##colnames(ctype_mat) <- cellTypes
    ##cellID <- rep(NA,nrow(ctype_mat))  #NULL
    ##cellPhenotype <- rep(NA,nrow(ctype_mat))
    ##for (i in c(1:length(cellTypes))){
    ##    ctype <- cellTypes[i]
    ##    tmpidx <- which(!is.na(ctype_mat[,ctype]))
    ##    tmp_cellID <- paste(prefixID,ctype_mat[tmpidx,ctype],sep="_")
    ##    cellID[tmpidx] <- tmp_cellID
    ##    cellPhenotype[tmpidx] <- ctype 
    ##}  
    ##cellID_info <- data.frame(cell.ID=cellID,HaloPhenotype=cellPhenotype)
    ##rownames(cellID_info) <- cellID
 
    # extract location
    radius <- data_in[,6]
    cell_position <- data_in[,c("XMin", "XMax", "YMin", "YMax")]
    cell_position$Loc.ID <- apply(data_in[,c("XMin", "XMax", "YMin", "YMax")],1,
                       function(x) paste0(x,collapse="_"))
    cell_position$Cell.X.Position <- (cell_position$XMin + cell_position$XMax)/2
    cell_position$Cell.Y.Position <- (cell_position$YMin + cell_position$YMax)/2
    cell_position$radius <- radius

    # remove duplicates
    keptR <- !duplicated(cell_position$Loc.ID) 
    cellID_info <- cell_position[keptR,]
    #cellID_info <- cbind(cellID_info,cell_position)

    return(cellID_info)
}
# *************************************************************************************************
# There are bugs in the marker_expression_boxplot --> corrected by LT 9/15/21
# **************************************************************************************************

marker_expression_boxplot_4DL <- function(sce_object, marker){

    formatted_data <- data.frame(colData(sce_object))

    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    expression_matrix <- assay(sce_object)

    cell_ids <- colnames(expression_matrix)
    markers  <- rownames(expression_matrix)   #***add by LT on 9/15/21

    rownames(expression_matrix) <- NULL
    colnames(expression_matrix) <- NULL
    expression_matrix_t <- t(expression_matrix)
    expression_df <- data.frame(expression_matrix_t)
    colnames(expression_df) <- markers

    formatted_data <- cbind(formatted_data, expression_df)
    formatted_data <- formatted_data[complete.cases(formatted_data),]

    #selecting cells that do express the marker
    expression_true <- formatted_data[grepl(marker, formatted_data$Phenotype), ]

    #selecting cells that do not contain the marker
    #for multiple entries that does not contain marker
    expression_false <- formatted_data[!grepl(marker, formatted_data$Phenotype), ] 

    #select the specific marker and add a boolean of expression
    if (nrow(expression_true) != 0) {
        expression_true$expression <- "P"
    } else{
        stop(paste("There are no cells positive for ", marker, sep=""))
    }

    if (nrow(expression_false) != 0){
        expression_false$expression <- "N"
    } else{
        stop(paste("There are no cells negative for ", marker, sep=""))
    }

    #bind the 2 dataframes together
    expression_by_marker <- rbind(expression_true, expression_false)
    expression_by_marker$expression <- factor(expression_by_marker$expression, levels=c("P", "N"))
    
    #plot boxplot
    title <- paste("Level of ", marker, sep="")

    #code from ggplot2 to show the number in each group
    give.n <- function(x){
        return(c(y = max(x)+1, label = length(x)))
    }

    p <- ggplot(expression_by_marker, aes(x = expression, y = expression_by_marker[, marker]))
    p <- p + geom_boxplot()
    p <- p + labs(title = title, x = "Marker status (Positive/Negative cell)", y = "Marker level")
    p <- p + stat_summary(fun.data = give.n, geom = "text")
    print(p)
}

# ************************************************************************************************
# select cells by specified cell types
# 9/13/21 by LT
# ***********************************************************************************************
select_cellsByVar <- function(sce_object, selVariableName, phenotypes){
  if (!all(selVariableName %in% colnames(colData(sce_object)))) {
        stop("selVariableName not found in meta-data")
    }
  tmpsel <- is.element(as.matrix(colData(sce_object)[,selVariableName]),phenotypes)
  sce_object2 <- sce_object[,tmpsel]
  return(sce_object2)
}

# ******************************************************************************************************
# Modify function to export 1 plot/marker on png, instead of pdf, on 09/15/21 by LT
#' plot_cell_marker_levels
#'
#' @description Produces a scatter plot of the level of every marker in each cell.
#' Cells that were not phenotyped as being positive for the particular marker are excluded.
#'
#' @param sce_object Singlecellexperiment object in the form of the output of format_image_to_sce
#' @param print TRUE for plots to be printed, FALSE otherwise
#' @param filename Path and name of output pdf file if to be printing to a file, if required
#' @param return_data TRUE if the function should return the formatted data for plotting
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @import qpdf
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#' @export
# %>% operator is in package 'magrittr' but imported by dplyr
# colData() is in package 'SummarizedExperiment' but imported by SingleCellExperiment
# *****************************************************************************************************
plot_cell_marker_levels_4DL <- function(sce_object, print = TRUE, filenamePrefix=c("tmp_"), return_data=TRUE) {

    formatted_data <- data.frame(colData(sce_object))

    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    expression_matrix <- assay(sce_object)

    markers <- rownames(expression_matrix)
    cell_ids <- colnames(expression_matrix)

    rownames(expression_matrix) <- NULL
    colnames(expression_matrix) <- NULL
    expression_matrix_t <- t(expression_matrix)
    expression_df <- data.frame(expression_matrix_t)
    colnames(expression_df) <- markers

    formatted_data <- cbind(formatted_data, expression_df)
    formatted_data <- formatted_data[complete.cases(formatted_data),]

    for (marker in markers) {

        #every cell is stained by DAPI, so no need to remove intensity
        if (marker == "DAPI"){
            next
        }

        #selecting cells that do not contain the marker
        rows <- formatted_data[formatted_data$Phenotype != marker, ] #for one entry that is not marker
        rows <- rows[!grepl(marker, rows$Phenotype), ] #for multiple entries that does not contain marker

        #for those cell without the marker, set marker intensity to 0
        #and merge the formatted_data
        rows[, marker] <- 0
        formatted_data[match(rows$Cell.ID,formatted_data$Cell.ID),]<-rows

    }

    if(print){
   
        for (marker in markers){
            # create filename
            filename <- paste(filenamePrefix,marker,".png",sep="")
            png(filename,width = 5, height = 5, units = "in",res=72)       

            #selecting the cells that have intensity for a specific marker
            column <- which(colnames(formatted_data) == marker)
            rows_non_zero <- which(formatted_data[,column] != 0)
            intensity_by_marker <- formatted_data[rows_non_zero,]

            if (nrow(intensity_by_marker) == 0) {
                print(paste("There are no true expression for: ", marker, sep=""))
            }

            #log the intensity to improve contrast
            intensity_by_marker[,marker] <- log10(intensity_by_marker[,marker])
            #print(intensity_by_marker)

            p <- ggplot(intensity_by_marker, aes(x = Cell.X.Position, 
                           y = Cell.Y.Position, colour = eval(parse(text = marker)))) +
                 geom_point(aes(colour=eval(parse(text = marker))),size = 0.1) +
                 ggtitle(marker)+
                 guides(alpha = F) + scale_colour_viridis_c(direction = -1) +
                 labs(colour = paste("log10","(", as.character(marker)," Intensity", ")", sep="")) +
                 theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.background = element_rect(fill = "white"),
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(), legend.key.height = unit(2.5, "cm"))

            print(p)
            dev.off()
        }
    }

    if(isTRUE(return_data)){
         return(formatted_data)
    }
}

# ************************************************************************************************
# modified plot_cell_marker_levels_4DL to plot a specific marker --> standard output
# 1/11/22
# ************************************************************************************************
plot_cell_specific_marker_levels_4DL <- function(sce_object,selMarker,isMakerPositiveCells) {

    formatted_data <- data.frame(colData(sce_object))

    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    expression_marker <- assay(sce_object)[selMarker,]

    formatted_data <- cbind(formatted_data, expression_marker)
    colnames(formatted_data)[ncol(formatted_data)] <- selMarker
    formatted_data <- formatted_data[complete.cases(formatted_data),]

    #selecting the cells that have intensity for a specific marker
    intensity_by_marker <- formatted_data[isMakerPositiveCells,]
    
    marker <- selMarker
    p <- ggplot(intensity_by_marker, aes(x = Cell.X.Position, 
                           y = Cell.Y.Position, colour = eval(parse(text = marker)))) +
         geom_point(aes(colour=eval(parse(text = marker))),size = 0.1) +
         ggtitle(marker)+
         guides(alpha = "none") + scale_colour_viridis_c(direction = -1) +
         labs(colour = paste("log10","(", as.character(marker)," Intensity", ")", sep="")) +
         theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_rect(fill = "white"),
               axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(), legend.key.height = unit(2.5, "cm"))
    return(p)
}



# ************************************************************************************
# Modify function to plot selected cell types
#' plot_cell_proportions
#'
#' @description Plots cells proportions as barplots
#' @param cell_proportions Output from calculate_cell_proportions
#' @param tumour_marker Tumour marker to exclude if needed
#' @import ggplot2
#' @export
# ************************************************************************************

plot_cell_percentages_selectCT4DL <- function(cell_proportions, selectCT=NULL){

  cell_proportions$Percentage_label <- round(cell_proportions$Percentage, digits=1)
  cell_proportions <- cell_proportions[cell_proportions$Cell_type != "OTHER",]

  if (is.null(selectCT)){
      cell_proportions_selectCT <- cell_proportions
  } else {
      cell_proportions_selectCT <- cell_proportions[cell_proportions$Cell_type %in% selectCT,]
      cell_proportions_selectCT$Proportion <- NULL
      cell_proportions_selectCT$Percentage <- (cell_proportions_selectCT$Number_of_cells/sum(cell_proportions_selectCT$Number_of_cells))*100
      cell_proportions_selectCT$Percentage_label <- round(cell_proportions_selectCT$Percentage, digits=1)
  }
  cell_percentages_selectCT_plot <-
       ggplot(cell_proportions_selectCT,aes(x = Cell_type,y = Percentage,
                                     fill = Cell_type)) +
       geom_bar(stat = 'identity') +
       ggtitle("Excluding tumor cells")+
       theme_classic() +
       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none") +
       xlab("Cell Type") + ylab("Proportion of total cells") +
        geom_text(aes(label = Percentage_label), position=position_dodge(width=0.9), vjust=-0.25, size = 2)
 
  print(cell_percentages_selectCT_plot)
}

# **************************************************************************************************************
# Convert data frame based on phenotypes (output of composition_of_cluster_and_community)
# into matrix communities x phenotypes whose element is % phenotype/communities (rowSum=1)
# Reuse the code from SPIAT plot_composition_heatmap
# *************************************************************************************************************
PercentPerCluster_table <- function(composition, column_to_consider) {

  if(column_to_consider == "Community"){
    cluster_size <- composition[!duplicated(composition$Community),
                                c("Community","Total_number_of_cells")]
    rownames(cluster_size) <- cluster_size$Community
    cluster_size$Community <- NULL

    composition2 <- composition[,c("Phenotype", "Community", "Percentage")]
    composition2 <- dcast(composition2,  Phenotype ~ Community, value.var="Percentage")

  }else if(column_to_consider == "Cluster"){
    cluster_size <- composition[!duplicated(composition$Cluster),
                                c("Cluster","Total_number_of_cells")]
    rownames(cluster_size) <- cluster_size$Cluster
    cluster_size$Cluster <- NULL

    composition2 <- composition[,c("Phenotype", "Cluster", "Percentage")]
    composition2 <- dcast(composition2,  Phenotype ~ Cluster, value.var="Percentage")
  }else{
    stop("Only Community and Cluster are accepted as valid column names")
  }
  
  rownames(composition2) <- composition2$Phenotype
  composition2$Phenotype <- NULL
  composition2[is.na(composition2)] <- 0  #***original = 0 
  composition2 <- as.matrix(composition2)/100

  colnames(cluster_size) <- c("Total_cells")
  list(cluster_size=cluster_size,cluster_composition=t(composition2))

}

# ****************************************************************************************************
# Determine meta-clusters/communities by clustering cell proportions of communitities
# The input is $cluster_composition (communities x phenotypes) of the output of PercentPerCluster_table 
# requires pheatmap
# 09/221/21
# ****************************************************************************************************
identify_metaclusters_4DL <- function(compositionMatrix,number_clusters,method=c("hclust")){
   
  no_phenotypes  <- ncol(compositionMatrix)
  if (method==c("hclust")){
     hm4_metaClusters <- pheatmap(as.matrix(compositionMatrix), scale = "none",
           cluster_cols=TRUE, cluster_rows=TRUE,
           show_rownames= FALSE,silent=TRUE)   #no plotting
     meta_clusters <- cutree(hm4_metaClusters$tree_row,number_clusters)
  } else if (method=="kmeans"){
     km.res <- kmeans(as.matrix(compositionMatrix), number_clusters, nstart = 25)
     meta_clusters <- km.res$cluster
  } else {
      stop("Current function supports hlcuts and kmeans")
  }
  tmp <- aggregate(compositionMatrix,by=list(meta_clusters),mean)
  metaID <- paste("metaC_",tmp$Group.1,sep="")
  rownames(tmp) <- metaID
  tmp$Group.1 <- NULL

  meta_cluster_mems <- data.frame(Community=names(meta_clusters),
                                  metaID= paste("metaC_",meta_clusters,sep=""))
  meta_cluster_info <- list(meta_members = meta_cluster_mems,
                            meta_profile = tmp)
  return(meta_cluster_info)
}

# ************************************************************************************************
# Transfer metaCluster to single-cell data
# Inputs: (1) colData, which is the output of identify_cell_communities
#         (2) metaCluster info, which is the $meta_members of identify_metaclusters_4DL
#         (3) specify column id for cluster/communities (shared between (1) and (2))
# (1) --> composition_of_clusters_and_communities --> PercentPerCluster_table -->
# identify_metaclusters_4DL --> (2)
# 09/22/21 by LT
# ************************************************************************************************
transfer_metaCluster2sc_4DL <- function(sc_Data,metaCluster_mems,column_to_consider){

  # Check input data
  if (!is.element(column_to_consider,colnames(sc_Data))){
      stop("column_to_consider not found in single-cell data")
  }
  if (!is.element(column_to_consider,colnames(metaCluster_mems))){
      stop("column_to_consider not found in meta-cluster data")
  }

  sc_metaCluster <- rep(NA,nrow(sc_Data))
  metaCluster_id <- names(table(metaCluster_mems$metaID))
  no_meta <- length(metaCluster_id)
  for (metaId in metaCluster_id){
     idx <- metaCluster_mems$metaID==metaId
     commid <- metaCluster_mems[idx,column_to_consider] 
     is_ele <- sc_Data[,column_to_consider] %in% commid
     sc_metaCluster[is_ele] <- metaId
  }
  names(sc_metaCluster) <- sc_Data$Cell.ID
  return(sc_metaCluster)
}

# ***********************************************************************************
# Add new phenotype to sce object
# input: original sce
#        newPhenotype: data frame of new phenotype
# output: sce with new slot "newPhenoID", 
#         "OTHER" indicates that newPheno is not available
# ************************************************************************************
add_newPheno_2sce_4DL <- function(sce_object,newPheno){
  
  original_sc <- data.frame(Cell.ID=colnames(sce_object))
  newPheno_df <- data.frame(Cell.ID=names(newPheno),
                          newPhenoID=newPheno)
  newPheno_sc <- left_join(original_sc,newPheno_df,
                   by = c("Cell.ID"="Cell.ID"))
  if (sum(is.na(newPheno_sc$newPhenoID))>0){
     newPheno_sc$newPhenoID[is.na(newPheno_sc$newPhenoID)] <- c("OTHER")
  }
  colData(sce_object)$newPhenoID <- newPheno_sc$newPhenoID
  return(sce_object)
}



# *********************************************************************************************
# plot_cell_categories has a bug in case no "OTHER" phenotype is specified
# 9/27/21
# ********************************************************************************************
# Modified spiat plot_cell_categories
# ********************************************************************************************
plot_cell_categories_4DL <- function(sce_object, phenotypes_of_interest, colour_vector) {
  
  formatted_data <- data.frame(colData(sce_object))
  
  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column
  
  formatted_data <- formatted_data[complete.cases(formatted_data),]
  
  #CHECK
  if (length(phenotypes_of_interest) != length(colour_vector)) {
    stop("The colour vector is not the same length as the phenotypes of interest")
  }
  for (phenotype in phenotypes_of_interest) {
    if (!(phenotype %in% unique(formatted_data$Phenotype))) {
      stop(paste(phenotype, "cells were not found"), sep="")
    }
  }
    
  #set all phenotypes of those that aren't in phenotypes_of_interest to be "OTHER"
  if (sum(!formatted_data$Phenotype %in% phenotypes_of_interest)>0) {
      formatted_data[!formatted_data$Phenotype %in% phenotypes_of_interest,]$Phenotype <- "OTHER"
  }
  
  #Assign the colour to corresponding phenotypes in df
  formatted_data$color <- ""
  for (phenotype in phenotypes_of_interest) {
    idx <- which(phenotypes_of_interest == phenotype)
    formatted_data[formatted_data$Phenotype == phenotype, ]$color <- colour_vector[idx]
  }
  if (sum(formatted_data$Phenotype == "OTHER")>0){
    formatted_data[formatted_data$Phenotype == "OTHER", ]$color <- "lightgrey"
  }
  
  
  all_phenotypes <- c(phenotypes_of_interest, "OTHER")
  all_colours <- c(colour_vector, "lightgrey")
  
  p <- ggplot(formatted_data, aes(x = Cell.X.Position, y = Cell.Y.Position, colour = Phenotype)) +
    geom_point(aes(colour = Phenotype), size = 1) +
    guides(alpha = F) +
    labs(colour = "Phenotypes") + 
    scale_color_manual(breaks = all_phenotypes, values=all_colours) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  print(p)
}



# ***********************************************************************************************
# Modified plot_cell_categories with input as data frame instead of sce object.
# Input
# sc_Data: data frame of meta-data, e.g. colData(sce)
# column_to_consider: name of phenotype column
# phenotypes_of_interest: subset of phenotype groups
# colour_vector: color code
# require RColorBrewer
# 9/21/21
# ***********************************************************************************************
plot_cell_categories_byDF <- function(sc_Data, column_to_consider,phenotypes_of_interest=NULL, colour_vector=NULL) {
  
  formatted_data <- sc_Data[,c("Cell.X.Position","Cell.Y.Position")]
 
  # CHECK 
  if (!is.element(c("Cell.ID"),colnames(sc_Data))) {
      formatted_data$Cell.ID <- rownames(sc_Data)
  }

  if (!is.element(column_to_consider,colnames(sc_Data)) ){
      stop("column_to_consider not found in single-cell data")
  }
    
  if (is.null(phenotypes_of_interest)){
      phenotypes_of_interest <- names(table(sc_Data[,column_to_consider]))
  } 

  if (is.null(colour_vector)){
     colour_vector <- colorRampPalette(brewer.pal(8, "Accent"))(length(phenotypes_of_interest))
  }

  formatted_data$Phenotype <- sc_Data[,column_to_consider]

  #set all phenotypes of those that aren't in phenotypes_of_interest to be "OTHER"
  if (sum(!formatted_data$Phenotype %in% phenotypes_of_interest)>0){
     formatted_data[!formatted_data$Phenotype %in% phenotypes_of_interest,]$Phenotype <- "OTHER"
  }
  
  #Assign the colour to corresponding phenotypes in df
  formatted_data$color <- ""
  for (phenotype in phenotypes_of_interest) {
    idx <- which(phenotypes_of_interest == phenotype)
    formatted_data[formatted_data$Phenotype == phenotype, ]$color <- colour_vector[idx]
  }
  if (sum(formatted_data$Phenotype=="OTHER")>0){
    formatted_data[formatted_data$Phenotype == "OTHER", ]$color <- "lightgrey"
  }
  
  all_phenotypes <- c(phenotypes_of_interest, "OTHER")
  all_colours <- c(colour_vector, "lightgrey")
  
  p <- ggplot(formatted_data, aes(x = Cell.X.Position, y = Cell.Y.Position, colour = Phenotype)) +
    geom_point(aes(colour = Phenotype), size = 1) +
    guides(alpha = F) +
    labs(colour = "Phenotypes") + 
    scale_color_manual(breaks = all_phenotypes, values=all_colours) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  print(p)

}

# ********************************************************************************
# calculate nearest distance among cells without contraint of cell type
# output: vector of no_cells
# modify from SPIAT average_minimal_distance to explore distribution instead of mean
# using RANN nn2 function
# 09/23/21 by LT
# ********************************************************************************
calculate_minimum_distance4DL <- function(sce_object) {

    formatted_data <- data.frame(colData(sce_object))
    #formatted_data <- formatted_data[complete.cases(formatted_data),]
    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    #extract the cell coordinates
    all_cell_cords <- formatted_data[,c("Cell.X.Position", "Cell.Y.Position")]
    
    #CHECK
    if (nrow(all_cell_cords) == 0) {
        stop("No cells found in average minimum distance calculation")
    }
    
    #calculate the closest 2 neighbours, 1st being itself
    all_closest <- nn2(data = all_cell_cords, k = 2)

    #grab the distances and find the average
    all_closest_dist <- all_closest$nn.dists[,2]
    
    return(all_closest_dist)
}

# ********************************************************************************
# calculate nearest distance from target to reference cell type
# output: vector of no_cells of reference cells
# modify from SPIAT calculate_summary_distances_between_phenotypes to get distribution
# instead of mean
# using RANN nn2 function
# 09/23/21 by LT
# ********************************************************************************
calculate_min_dist_bet2CTs_4DL <- function(sce_object,reference_phenotype, target_phenotype) {

    formatted_data <- data.frame(colData(sce_object))
    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    formatted_data <- formatted_data[,c("Cell.ID","Phenotype", "Cell.X.Position", "Cell.Y.Position")]
    formatted_data <- formatted_data[formatted_data$Phenotype != "",]

    # Add a new column Cell_type which duplicates the phenotype
    formatted_data$Cell_type <- as.character(formatted_data$Phenotype)

    # Get coordinates
    all_celltype1_cord <- formatted_data[formatted_data$Phenotype == reference_phenotype, 
                                  c("Cell.X.Position", "Cell.Y.Position")]
    all_celltype2_cord <- formatted_data[formatted_data$Phenotype == target_phenotype,
                                  c("Cell.X.Position", "Cell.Y.Position")]
   
    # Find nearest distance of target to each reference cells: k=1 means nearest
    all_closest_2to1 <- nn2(data = all_celltype2_cord, query = all_celltype1_cord, k = 1)
    min_dist_2to1 <- all_closest_2to1$nn.dists

    return(min_dist_2to1)
}  

# **********************************************************************************************************
# calculate nearest distance from target to reference cell type
# similar to calculate_min_dist_bet2CTs_4DL, but input is matrix
# reference ct is the same as target ct, calculate distance to the nearest neighbor
# 5/19/22
# **********************************************************************************************************
calc_min_dist_bet2CTs_4DL_matInput <- function(formatted_data,reference_phenotype, target_phenotype){
    
    formatted_data <- formatted_data[,c("Cell.ID","Phenotype", "Cell.X.Position", "Cell.Y.Position")]

    # keep only specified cell types
    ##formatted_data <- formatted_data[formatted_data$Phenotype %in% cell_phenotypes_of_interest,]

    # Add a new column Cell_type which duplicates the phenotype
    formatted_data$Cell_type <- as.character(formatted_data$Phenotype)

    # Get coordinates
    all_celltype1_cord <- formatted_data[formatted_data$Phenotype == reference_phenotype, 
                                  c("Cell.X.Position", "Cell.Y.Position")]
    all_celltype2_cord <- formatted_data[formatted_data$Phenotype == target_phenotype,
                                  c("Cell.X.Position", "Cell.Y.Position")]
   
    # Find nearest distance of target to each reference cells: k=1 means nearest
    if (target_phenotype==reference_phenotype){
        all_closest_2to1 <- nn2(data = all_celltype1_cord,k = 2)
        min_dist_2to1 <- all_closest_2to1$nn.dists[,2]
    } else {
        if (is.null(all_celltype2_cord)){
            min_dist_2to1 <- rep(NA,nrow(all_celltype1_cord))
        }else{
            all_closest_2to1 <- nn2(data = all_celltype2_cord, query = all_celltype1_cord, k = 1)
            min_dist_2to1 <- all_closest_2to1$nn.dists
        }
    }  
    return(min_dist_2to1)
}

# *****************************************************************************
swapper_calc_nndist_bet2CTs_4DL_matInput <- function(formatted_data,centroidCT,selTargetCT){
   vec4sampleID <- names(table(formatted_data$ROI))
   nsamples <- length(vec4sampleID)

   distNN_2selCT <- NULL
   distNN_2selCT.info <- NULL
   for (i in c(1:nsamples)){
      tmpin <- subset(formatted_data,ROI==vec4sampleID[i])
      tmpref_idx <- which(tmpin$Phenotype==centroidCT)
      if (length(tmpref_idx)>1){
         tmpout_id <- tmpin[tmpref_idx,c("Cell.ID","Phenotype","ROI")]
         tmpout <- matrix(NA,length(tmpref_idx),length(selTargetCT))         
         for (j in c(1:length(selTargetCT))){
            istarget_exist <- sum(tmpin$Phenotype==selTargetCT[j])
            if (istarget_exist>0){
               tmpout[,j] <- calc_min_dist_bet2CTs_4DL_matInput(tmpin,
                               reference_phenotype=centroidCT, 
                               target_phenotype=selTargetCT[j])    
            }   
         }
         colnames(tmpout) <- selTargetCT
         tmpout <- bind_cols(tmpout_id,as.data.frame(tmpout))
         distNN_2selCT <- bind_rows(distNN_2selCT,tmpout)
      }      
   } 
   rownames(distNN_2selCT) <- NULL
   return(distNN_2selCT)    
}


# *********************************************************************************************************************
# calculate nearest distance between each cell to other cell types --> similar to akoya compute_all_nearest_distance
# inputs
# 1) sce object
# 2) list of phenotypes, e.g. names(table(colData(sce_object)$Phenotype))
# requires RANN::nn2 function
# modified from SPIAT calculate_summary_distance_between_phenotypes
# 09/23/21 LT
# **********************************************************************************************************************
calc_all_nndist_4DL <- function(sce_object,cell_phenotypes_of_interest){
    
    formatted_data <- data.frame(colData(sce_object))
    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    formatted_data <- formatted_data[,c("Cell.ID","Phenotype", "Cell.X.Position", "Cell.Y.Position")]

    # keep only specified cell types
    formatted_data <- formatted_data[formatted_data$Phenotype %in% cell_phenotypes_of_interest,]

    # Add a new column Cell_type which duplicates the phenotype
    formatted_data$Cell_type <- as.character(formatted_data$Phenotype)

    unique_cells <- names(table(formatted_data$Cell_type))

    # Create matrix distance between cell types
    result_dist  <- matrix(0,nrow=nrow(formatted_data),ncol=length(unique_cells))
    result_id  <- matrix(NA,nrow=nrow(formatted_data),ncol=length(unique_cells))
    nct <- length(unique_cells)
    for (i in c(1:nct)){
        ref_cell_idx <- which(formatted_data$Cell_type == unique_cells[i]) 
        if (length(ref_cell_idx)>0){
            all_celltype1_cord <- formatted_data[ref_cell_idx, 
                                  c("Cell.X.Position", "Cell.Y.Position")]
            for (j in c(1:nct)){
                if(j==i){    #with same cell types for target and reference
                    if (nrow(all_celltype1_cord)>1){
                        dist_2to1 <- nn2(data = all_celltype1_cord, k = 2)
                        result_dist[ref_cell_idx,j] <- dist_2to1$nn.dists[,2]
                        result_id[ref_cell_idx,j]  <- formatted_data[ref_cell_idx[dist_2to1$nn.idx[,2]],c("Cell.ID")]
                    }else{
                        result_id[ref_cell_idx,j] <- NA
                    }
                } else {
                    target_cell_idx <- which(formatted_data$Cell_type == unique_cells[j]) 
                    all_celltype2_cord <- formatted_data[target_cell_idx, 
                                  c("Cell.X.Position", "Cell.Y.Position")]
                    dist_2to1 <- nn2(data = all_celltype2_cord, query = all_celltype1_cord, k = 1)
                    result_dist[ref_cell_idx,j] <- dist_2to1$nn.dists[,1]
                    result_id[ref_cell_idx,j]  <- formatted_data[target_cell_idx[dist_2to1$nn.idx[,1]],c("Cell.ID")]
                }
             }
         }
     }
     colnames(result_dist) <- unique_cells
     rownames(result_dist) <- formatted_data[,c("Cell.ID")]

     colnames(result_id) <- unique_cells
     rownames(result_id) <- formatted_data[,c("Cell.ID")]
    
     ref_pheno <- formatted_data$Cell_type

     results <- list(min_dist = result_dist, nn_id = result_id,  ref_pheno = ref_pheno)
     results
}

# *********************************************************************************************************************
# calculate nearest distance between each cell to other cell types --> similar to akoya compute_all_nearest_distance
# inputs
# 1) matrix/data frame (colData(sce)
# 2) list of phenotypes, e.g. names(table(colData(sce_object)$Phenotype))
# requires RANN::nn2 function
# modified from SPIAT calculate_summary_distance_between_phenotypes
# 01/06/22 LT
# **********************************************************************************************************************
calc_all_nndist_matInput <- function(formatted_data,cell_phenotypes_of_interest){
    
    formatted_data <- formatted_data[,c("Cell.ID","Phenotype", "Cell.X.Position", "Cell.Y.Position")]

    # keep only specified cell types
    formatted_data <- formatted_data[formatted_data$Phenotype %in% cell_phenotypes_of_interest,]

    # Add a new column Cell_type which duplicates the phenotype
    formatted_data$Cell_type <- as.character(formatted_data$Phenotype)

    unique_cells <- names(table(formatted_data$Cell_type))

    # Create matrix distance between cell types
    result_dist  <- matrix(0,nrow=nrow(formatted_data),ncol=length(unique_cells))
    result_id  <- matrix(NA,nrow=nrow(formatted_data),ncol=length(unique_cells))
    nct <- length(unique_cells)
    for (i in c(1:nct)){
        ref_cell_idx <- which(formatted_data$Cell_type == unique_cells[i]) 
        if (length(ref_cell_idx)>0){
            all_celltype1_cord <- formatted_data[ref_cell_idx, 
                                  c("Cell.X.Position", "Cell.Y.Position")]
            for (j in c(1:nct)){
                if(j==i){    #with same cell types for target and reference
                    if (nrow(all_celltype1_cord)>1){
                        dist_2to1 <- nn2(data = all_celltype1_cord, k = 2)
                        result_dist[ref_cell_idx,j] <- dist_2to1$nn.dists[,2]
                        result_id[ref_cell_idx,j]  <- formatted_data[ref_cell_idx[dist_2to1$nn.idx[,2]],c("Cell.ID")]
                    }else{
                        result_id[ref_cell_idx,j] <- NA
                    }
                } else {
                    target_cell_idx <- which(formatted_data$Cell_type == unique_cells[j]) 
                    all_celltype2_cord <- formatted_data[target_cell_idx, 
                                  c("Cell.X.Position", "Cell.Y.Position")]
                    dist_2to1 <- nn2(data = all_celltype2_cord, query = all_celltype1_cord, k = 1)
                    result_dist[ref_cell_idx,j] <- dist_2to1$nn.dists[,1]
                    result_id[ref_cell_idx,j]  <- formatted_data[target_cell_idx[dist_2to1$nn.idx[,1]],c("Cell.ID")]
                }
             }
         }
     }
     colnames(result_dist) <- unique_cells
     rownames(result_dist) <- formatted_data[,c("Cell.ID")]

     colnames(result_id) <- unique_cells
     rownames(result_id) <- formatted_data[,c("Cell.ID")]
    
     ref_pheno <- formatted_data$Cell_type

     results <- list(min_dist = result_dist, nn_id = result_id,  ref_pheno = ref_pheno)
     results
}


# *****************************************************************************************************
# Statistical summary of min distance between 2 cells types.
# row is reference cells of interest (COI), column is other cells to COI
# input is the output of calc_all_nndist_4DL
# 9/23/21 by LT
# *****************************************************************************************************
stat_summary_nndist_4DL <- function(all_nndist){
   
    unique_cells <- names(table(all_nndist$ref_pheno))
    nct <- length(unique_cells)

    avg_dist <- matrix(0,nct,nct)
    median_dist <- matrix(0,nct,nct)
    std_dist <- matrix(0,nct,nct)
    for (j in c(1:nct)){
        tmp_avg <- aggregate(all_nndist$min_dist[,j],by = list(all_nndist$ref_pheno),mean,na.rm=T)
        avg_dist[,j] <- tmp_avg[,2]
        tmp_std <- aggregate(all_nndist$min_dist[,j],by = list(all_nndist$ref_pheno),sd,na.rm=T)
        std_dist[,j] <- tmp_std[,2]
        tmp_median <- aggregate(all_nndist$min_dist[,j],by = list(all_nndist$ref_pheno),median,na.rm=T)
        median_dist[,j] <- tmp_median[,2]
    }
    rownames(avg_dist) <- unique_cells
    colnames(avg_dist) <- colnames(all_nndist$min_dist)

    rownames(std_dist) <- unique_cells
    colnames(std_dist) <- colnames(all_nndist$min_dist)

    rownames(median_dist) <- unique_cells
    colnames(median_dist) <- colnames(all_nndist$min_dist)

    results <- list(avg_dist=avg_dist,std_dist=std_dist,median_dist=median_dist)
    results
}

# *********************************************************************************************************************
# calculate components of reference all + k=10 nearest neighbors of each cell --> Nolan's lab approach
# inputs
# 1) sce object
# 2) knn: set number of nearest neighbors
# output is the list of
# ref_knnCT: matrix in which each row is the component of knn of each cell in ROI
# ref_pheno: phenotype of reference cells (i.e. colData(sce_object)$Phenotype
# requires RANN::nn2 function
# modified from SPIAT calculate_summary_distance_between_phenotypes
# 09/23/21 LT
# **********************************************************************************************************************
calc_component_knn_4DL <- function(sce_object,knn = 10){
    
    formatted_data <- data.frame(colData(sce_object))
    formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

    formatted_data <- formatted_data[,c("Cell.ID","Phenotype", "Cell.X.Position", "Cell.Y.Position")]

    # keep only specified cell types
    ##if (!is.null(cell_phenotypes_of_interest)){
    ##    formatted_data <- formatted_data[formatted_data$Phenotype %in% cell_phenotypes_of_interest,]
    ##}

    # Add a new column Cell_type which duplicates the phenotype
    formatted_data$Cell_type <- as.character(formatted_data$Phenotype)

    unique_cells <- names(table(formatted_data$Cell_type))

    # Create matrix components: ncells x cellTypes
    result_component  <- matrix(0,nrow=nrow(formatted_data),ncol=length(unique_cells))
    nct <- length(unique_cells)

    all_celltype_cord <- formatted_data[, 
                                  c("Cell.X.Position", "Cell.Y.Position")]
    knn_info <- nn2(data = all_celltype_cord, k = (knn+1))
    knn_id <- knn_info$nn.idx
    knn_ct <- t(apply(knn_id,1,function(x) formatted_data$Cell_type[x]))
    knn_component <- sapply(unique_cells,function(x) rowSums(knn_ct==x))
    colnames(knn_component) <- unique_cells
    rownames(knn_component) <- formatted_data$Cell.ID
    #knn_component <- knn_component/(1+knn)

    # create sce object with assay as cell components
    sce_out <- SingleCellExperiment(assays = list(counts = t(knn_component)))
    rownames(sce_out) <- colnames(knn_component) 
    colnames(sce_out) <- rownames(knn_component) 
    colData(sce_out)$Cell.X.Position  <- formatted_data$Cell.X.Position
    colData(sce_out)$Cell.Y.Position  <- formatted_data$Cell.Y.Position 
    colData(sce_out)$Phenotype  <- formatted_data$Phenotype
 
    knn_out <- list(ref_info=formatted_data[,c("Cell.ID","Phenotype")],
                    ref_knnCT=knn_component,sce=sce_out)
    return(knn_out)
}

# ****************************************************************************
# modify  spiat percentage_of_cells_within_radius to calculate % of other cell
# type (i.e. others != ref_type) within setting radius. 
# spiat percentage_of_cells_within_radius: pairwise relationship between target
#     and reference cell type (how namy target within r from reference cell type
# this is analogous to calc_component_knn_4DL, which is based on k nearest eneighbors
#     --> using RANN:nn function and output as NUMBER instead of %
# require dbscan frNN
# inputs
# 1) sce object
# 2) radius: set distance within each cell
# output is the list of
# otherCT_component: matrix in which each row is the component of all CT in ROI
# ref_pheno: phenotype of reference cells (i.e. colData(sce_object)$Phenotype
# sce: component data in sce object format
# requires RANN::nn2 function
# 9/29/21 LT
# *****************************************************************************
calc_component_radius_4DL <- function(sce_object, radius = 100){

  formatted_data <- data.frame(colData(sce_object))
  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

  # Add a new column Cell_type which duplicates the phenotype
  formatted_data$Cell_type <- as.character(formatted_data$Phenotype)

  unique_cells <- names(table(formatted_data$Cell_type))

  # Create matrix components: ncells x cellTypes
  result_component  <- matrix(0,nrow=nrow(formatted_data),ncol=length(unique_cells))
  nct <- length(unique_cells)

  #fixed radius nearest neighbors for each point. Note: self-matches are not returned for frNN (different from RANN::nn
  search_result <- frNN(formatted_data[,c("Cell.X.Position", "Cell.Y.Position")],
                        eps = radius, sort = FALSE)
  search_result_no <- sapply(search_result$id,length)
  if (min(search_result_no)<2){
      stop("***Some cells have < 2 neighbors based on set up radius, calculation aborted***")
  } else {
      search_result_id <- lapply(search_result$id,function(x) formatted_data$Cell_type[x])
      search_result_component <- matrix(0,nrow=nrow(formatted_data),ncol=length(unique_cells))
      for (i in c(1:length(unique_cells))){
          search_result_component[,i] <- unlist(lapply(search_result_id,
                                                  function(x) sum(x==unique_cells[i])))
      }
      #search_result_component <- sapply(unique_cells,function(x) rowSums(search_result_id==x))
      colnames(search_result_component) <- unique_cells
      rownames(search_result_component) <- formatted_data$Cell.ID
  }

  # create sce object with assay as cell components
  sce_out <- SingleCellExperiment(assays = list(counts = t(search_result_component)))
  rownames(sce_out) <- colnames(search_result_component) 
  colnames(sce_out) <- rownames(search_result_component) 
  colData(sce_out)$Cell.X.Position  <- formatted_data$Cell.X.Position
  colData(sce_out)$Cell.Y.Position  <- formatted_data$Cell.Y.Position 
  colData(sce_out)$Phenotype  <- formatted_data$Phenotype
 
  search_result_out  <- list(ref_info=formatted_data[,c("Cell.ID","Phenotype")],
                    otherCT_component=search_result_component,sce=sce_out)

  return(search_result_out)
}

# ******************************************************************************

calc_component_radiusSubnet_4DL <- function(sce_object, radius = 100){

  formatted_data <- data.frame(colData(sce_object))
  formatted_data <- formatted_data %>% rownames_to_column("Cell.ID") #convert rowname to column

  # Add a new column Cell_type which duplicates the phenotype
  formatted_data$Cell_type <- as.character(formatted_data$Phenotype)

  unique_cells <- names(table(formatted_data$Cell_type))

  # Create matrix components: ncells x cellTypes
  result_component  <- matrix(0,nrow=nrow(formatted_data),ncol=length(unique_cells))
  nct <- length(unique_cells)

  #fixed radius nearest neighbors for each point. Note: self-matches are not returned for frNN (different from RANN::nn
  search_result <- frNN(formatted_data[,c("Cell.X.Position", "Cell.Y.Position")],
                        eps = radius, sort = FALSE)
  search_result_no <- sapply(search_result$id,length)
  keptR <- which(search_result_no>2)
  kept_result_id <- lapply(search_result$id[keptR],function(x) formatted_data$Cell_type[x])
  kept_result_component <- matrix(0,nrow=length(keptR),ncol=length(unique_cells))
  for (i in c(1:length(unique_cells))){
      kept_result_component[,i] <- unlist(lapply(kept_result_id,
                               function(x) sum(x==unique_cells[i])))
  }
  colnames(kept_result_component) <- unique_cells
  rownames(kept_result_component) <- formatted_data$Cell.ID[keptR]
 
  # create sce object with assay as cell components
  sce_out <- SingleCellExperiment(assays = list(counts = t(kept_result_component)))
  rownames(sce_out) <- colnames(kept_result_component) 
  colnames(sce_out) <- rownames(kept_result_component) 
  colData(sce_out)$Cell.X.Position  <- formatted_data$Cell.X.Position[keptR]
  colData(sce_out)$Cell.Y.Position  <- formatted_data$Cell.Y.Position[keptR] 
  colData(sce_out)$Phenotype  <- formatted_data$Phenotype[keptR]
 
  search_result_out  <- list(ref_info=formatted_data[keptR,c("Cell.ID","Phenotype")],
                    otherCT_component=kept_result_component,sce=sce_out)

  return(search_result_out)
}

# ****************************************************************************
# Calculate component of neighbors with setting radius
# modify  from 
# 1. spiat percentage_of_cells_within_radius to calculate % of other cell
# type (i.e. others != ref_type) within setting radius. 
# 2. calc_component_radius_4DL
# 3. calc_component_radiusSubnet_4DL
# New version: input as matrix of cell location and phenotype
#              radius = 40 (90% epi cells with diameter = 20 in human data)
# 1/26/22
# ****************************************************************************
calc_component_radiusSubnet_matInput <- function(formatted_data, radius = 40,cell_phenotypes_of_interest=NULL){

  formatted_data <- formatted_data[,c("Cell.ID","Phenotype", "Cell.X.Position", "Cell.Y.Position")]

  # Add a new column Cell_type which duplicates the phenotype
  formatted_data$Cell_type <- as.character(formatted_data$Phenotype)
  
  if (is.null(cell_phenotypes_of_interest)){
      unique_cells <- names(table(formatted_data$Cell_type))
  } else {
      unique_cells <- cell_phenotypes_of_interest
  }

  # Create matrix components: ncells x cellTypes
  nct <- length(unique_cells)
    
  #fixed radius nearest neighbors for each point. Note: self-matches are not returned for frNN (different from RANN::nn
  search_result <- frNN(formatted_data[,c("Cell.X.Position", "Cell.Y.Position")],
                        eps = radius, sort = FALSE)
  search_result_no <- sapply(search_result$id,length)
  keptR <- which(search_result_no>0)
  kept_result_id <- lapply(search_result$id[keptR],function(x) formatted_data$Cell_type[x])
  kept_result_component <- matrix(0,nrow=length(keptR),ncol=nct)
  for (i in c(1:nct)){
      kept_result_component[,i] <- unlist(lapply(kept_result_id,
                               function(x) sum(x==unique_cells[i])))
  }
  colnames(kept_result_component) <- unique_cells
  rownames(kept_result_component) <- formatted_data$Cell.ID[keptR]
 
  # create output list
  out_tab <- matrix(0,nrow=nrow(formatted_data),ncol=nct)
  out_tab[keptR,] <- kept_result_component
  colnames(out_tab) <- unique_cells
  rownames(out_tab) <- formatted_data$Cell.ID
  out_tab <- cbind(search_result_no,out_tab)
  ##search_result_out  <- list(ref_info=formatted_data[keptR,],
  ##                  otherCT_component=kept_result_component)

  return(out_tab)
}


# **************************************************************************
## color code from CATALYST
gg_color_catalyst <- function(number_of_colours){
   #kpal <- CATALYST:::.cluster_cols
   kpal <- c("#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72",
             "#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",
             "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D",
             "#E6AB02","#7570B3","#BEAED4","#666666","#999999",
             "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000",
             "#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00")
   if (number_of_colours > length(kpal)){
      kpal <- colorRampPalette(kpal)(number_of_colours)
   }
   kpal[seq_len(number_of_colours)]
}

# ************************************************************************
#Title: gg_color_hue.R
#Author: tinyheero
#Date: 2nd December, 2019
#Availability: https://github.com/tinyheero/tinyutils/blob/master/R/gg_color_hue.R
gg_color_hue <- function(number_of_colours) {
  hues = seq(15, 375, length = number_of_colours + 1)
  hcl(h = hues, l = 65, c = 100)[1:number_of_colours]
}

# ************************************************************************************
dist2ref <- function(x,ref_pt){
   sqrt(sum((x-ref_pt)^2))
}