# WCS Effectiveness Chapter

# Objectives and Hypotheses
## Objectives
  ## 1. Determine how community composition of WCS changes through time
  ## 2. Identify which factors are most important for predicting WCS beta diversity
  ## 3. Develop regression model for predicting WCS use
## Hypotheses
  ## 1. WCS beta diversity will change over time
  ## 2. Structural characteristics will be the most influential characteristics in predicting WCS use 
  ## 3. A regression model containing structural, anthropogenic, and environmental characteristics of 
      ##WCS will effectively model the beta diversity of a site

# Response Variables
## 1. Total number of detections of each species at the WCS
## 2. Number of successful interactions at a WCS
## 3. Number of unsuccessful interactions at a WCS

devtools::install_github("tomyamashita/cameraTrapping")

rm(list = ls()); gc()

#######################################################################################################################
######################################## Part 0: Functions used in this script ########################################
#######################################################################################################################

# Used in Part 1: Prepping raw data ####
## Process the interactions data
processInts <- function(x, species){
  #x <- ints1
  #species <- spec_analysis
  
  ## Subset only those species that will be analysed (wild mammals)
  ints2 <- x[x$Species %in% species,]
  #sort(unique(ints2$Species))
  
  ## Separate date_time information to calculate actual date-times
  ints3 <- tidyr::separate(ints2, Date_Time, sep = " ", into = c("year", "month", "day", "hour", "minute", "second", "serial"), fill = "right")
  ints3$dt <- with(ints3, lubridate::ymd_hms(paste(year, month, day, hour, minute, second, sep = " ")))
  
  ## Assign the appropriate identifiers for each to identify experimental units
  ints3$WCS_ID <- with(ints3, paste(Highway, Site, sep = "_"))
  ints3$TimeID <- with(ints3, paste(year, formatC(month, width = 2, flag = "0"), sep = ""))
  ints3$WCS_Time_ID <- with(ints3, paste(WCS_ID, TimeID, sep = "_"))
  #ints3[1:5,]
  
  ## Summarize the data by experimental unit, species, and interaction type (Class) and calculate the number of instances of each class and total individuals for each class
  ints4 <- tidyr::pivot_wider(dplyr::summarise(dplyr::group_by(ints3, WCS_Time_ID, WCS_ID, TimeID, Species, Class), ints_sum = sum(Individuals)), 
                              names_from = Class, values_from = ints_sum, names_sort = TRUE, values_fill = 0)
  ints4$Sort <- paste(ints4$WCS_Time_ID, ints4$Species, sep = "__")
  #ints4[1:5,]
  
  ## For some of the experimental units, we did not assess number of interactions so we need to address this issue
  noints <- c("FM106_WCS1_201907", "FM106_WCS1_201908", "FM106_WCS1_201909", "FM106_WCS1_201910", "FM106_WCS1_201911", 
              "FM106_WCS2_201907", "FM106_WCS2_201908", "FM106_WCS2_201909", "FM106_WCS2_201910", "FM106_WCS2_201911", 
              "FM106_WCS3_201907", "FM106_WCS3_201908", "FM106_WCS3_201909", "FM106_WCS3_201910", "FM106_WCS3_201911", 
              "FM106_WCS4_201907", "FM106_WCS4_201908", "FM106_WCS4_201909", "FM106_WCS4_201910", "FM106_WCS4_201911", 
              "FM106_WCS5_201907", "FM106_WCS5_201908", "FM106_WCS5_201909", "FM106_WCS5_201910", "FM106_WCS5_201911", 
              "FM106_WCS6_201907", "FM106_WCS6_201908", "FM106_WCS6_201909", "FM106_WCS6_201910", "FM106_WCS6_201911", 
              "FM106_WCS7_201907", "FM106_WCS7_201908", "FM106_WCS7_201909", "FM106_WCS7_201910", "FM106_WCS7_201911", 
              "FM106_WCS8_201907", "FM106_WCS8_201908", "FM106_WCS8_201909", "FM106_WCS8_201910", "FM106_WCS8_201911")
  ints4$noints <- ints4$WCS_Time_ID %in% noints
  #ints4[1:5,]
  ints4$A <- ifelse(ints4$noints, NaN, ints4$A)
  ints4$B <- ifelse(ints4$noints, NaN, ints4$B)
  ints4$C <- ifelse(ints4$noints, NaN, ints4$C)
  ints4$D <- ifelse(ints4$noints, NaN, ints4$D)
  ints4$E <- ifelse(ints4$noints, NaN, ints4$E)
  ints4$`NA` <- ifelse(!ints4$noints & ints4$`NA` == 0, NaN, ints4$`NA`)
  #ints4[1:5,]
  
  ## Now we can calculate our response variables
  ints4$Total <- apply(ints4[,c("A", "B", "C", "D", "E", "NA")], 1, sum, na.rm = TRUE)
  ints4$Success <- with(ints4, A)
  ints4$Failure <- with(ints4, B + C + E)
  ints4$A_total <- with(ints4, Success / Total)
  ints4$A_ints <- with(ints4, Success / (Success + Failure))
  ints4$DaysInMonth <- lubridate::days_in_month(lubridate::ym(ints4$TimeID))  # Assumes that at least one camera was active at each site for the duration of the study
  
  ## Calculate the number of detections, successes, and failures per day as a standardization technique for monthly data
  ints4$TotalPD <- with(ints4, Total/DaysInMonth)
  ints4$SuccessPD <- with(ints4, Success/DaysInMonth)
  ints4$FailurePD <- with(ints4, Failure/DaysInMonth)
  #ints4[1:5,]
  
  ## View a simplified output
  #ints4[1:5,c("WCS_Time_ID", "Species", "TotalPD", "SuccessPD", "FailurePD", "A_total", "A_ints")]
  
  ## Create site x species table for each response variable
  fields <- c("TotalPD", "SuccessPD", "FailurePD")
  ints5 <- lapply(fields, function(x){
    x1 <- ints4[,c("WCS_Time_ID", "Species", x)]
    tidyr::pivot_wider(x1, names_from = "Species", names_sort = TRUE, values_from = all_of(x), values_fill = 0)
  })
  names(ints5) <- fields
  #ints5$TotalPD[1:5,]
  #ints5$SuccessPD[1:20,]
  
  return(ints5)
  rm(ints2, ints3, ints4, ints5, noints, fields)
  rm(x, species)
}

## Calculation of a distance matrix  (function modified from Paul and Anderson 2013)
distance <- function(Y, measure = "BC", trans = "none", adj = 0){
  #Y <- Yme  # Data matrix
  #measure <- "adjBC"  # Dissimilarity type
  #trans <- "none"  # Transformation type (before calculation of dissimilarity)
  
  # Determine number of dimensions of the input data matrix
  n <- dim(Y)[1]; p <- dim(Y)[2]
  Y <- switch(trans, 
              none = Y, 
              sqrt = sqrt(Y), 
              fourthroot = Y^0.25, 
              pa = (Y>0),
              rowpropns = Y/apply(Y, 1, sum),
              stop("\n", "Unrecognized transformation",
                   "\n", "Must be one of: none, sqrt, fourthroot, pa, rowpropns", "\n"))
  
  # Check that there are at least 2 rows (sites)
  if(n < 2){stop("Must be at least 2 sites")}
  
  # Calculate the dissimilarities
  ## Create an empty data matrix
  D <- matrix(0, nrow = n, ncol = n)
  ## Fill in the upper triangle of the matrix
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      D[i,j] <- switch(measure,
                       # Bray-curtis and zero-adjusted bray-curtis (adjusted bray-curtis added by me)
                       BC = (sum(abs(Y[i,] - Y[j,]))/(adj * 2 + sum(Y[i,] + Y[j,]))), 
                       # I don't know what the next 3 are
                       sqrtBC = (sum(abs(Y[i,]-Y[j,]))/sum(Y[i,]+Y[j,]))^0.5,
                       Can = sum(ifelse(Y[i,]+Y[j,]>0,abs(Y[i,]-Y[j,])/(Y[i,]+Y[j,]),0)),
                       sqrtCan = sum(ifelse(Y[i,]+Y[j,]>0,abs(Y[i,]-Y[j,])/(Y[i,]+Y[j,]),0))^0.5,
                       # Euclidean distance
                       Eucl = sqrt(sum((Y[i,]-Y[j,])^2)),
                       stop("\n", "Unrecognized distance measure", 
                            "\n", "Must be one of ", "BC, sqrtBC, Can, sqrtCan, or Eucl", 
                            "\n")
      )
    }
  }
  
  ## Fill in the lower triangle of the matrix with the appropriate data
  D <- D + t(D)
  return(D)
}

## Conduct distance-based Moran's Eigenvector map processing
dbmemFun <- function(x, type){
  #x <- cov.all.list$temporal
  #type <- "spatial"
  
  # Define site information as a matrix
  x1 <- as.matrix(x[,-1])
  dimnames(x1)[[1]] <- x[,1]
  
  # Compute a Euclidean distance matrix on the cartesian coordinates
  d.x <- distance(Y = x1, measure = "Eucl", trans = "none", adj = 0)
  
  ## NOTE: 
    ## We could run an AEM on the temporal data by replacing all values that are greater than 0 in the distance matrix with 1's 
    ## See Legendre's paper on AEMs to see how this makes sense
    ## Then you would just run a PCO on this matrix to compute the AEM
  
  # Examine the minimum spanning tree to determine threshold distance for dbMEM
  mst.x <- vegan::spantree(d = d.x)
  print(summary(mst.x$dist))
  
  dbmem <- vegan::pcnm(dis = d.x, threshold = max(mst.x$dist), dist.ret = TRUE)
  
  x2 <- data.frame(WCS_Time_ID = row.names(x1), dbmem$vectors)
  if(type == "spatial"){
    colnames(x2) <- c("WCS_Time_ID", paste("s_", colnames(x2)[-1], sep = ""))
  }else if(type == "temporal"){
    colnames(x2) <- c("WCS_Time_ID", paste("t_", colnames(x2)[-1], sep = ""))
  }else{
    message("You did not correctly specify type. No adjustment to the axis names are being made. Use type = c('spatial', 'temporal').")
  }
  
  return(list(dbmem = dbmem, df = x2))
  rm(x1, d.x, mst.x, dbmem, x2)
  rm(x, type)
}


# Used in Part 2: Plotting PCO and dbRDA ####
## Function for prepping the scores data
ordPlotPrep <- function(files, plot.type, score.type, axes){
  #files <- files_plot
  #plot.type <- "dbRDA"
  #score.type <- "pred"
  #axes <- 4
  
  if(plot.type == "PCO"){
    fs1 <- files[grep("PCO", files)]
  }else if(plot.type == "dbRDA"){
    fs1 <- files[grep("dbRDA", files)]
  }else{
    stop("You must choose a valid plot type. Choose one of c('PCO', 'dbRDA').")
  }
  
  if(score.type == "site"){
    fs2 <- fs1[grep("site", fs1)]
    scores1 <- pbapply::pblapply(1:length(fs2), function(i){
      name1 <- sub(paste(plot.type, "_", sep = ""), "", sub(paste("_", score.type, sep = ""), "", basename(fs::path_ext_remove(fs2[i]))))
      x <- openxlsx::read.xlsx(fs2[i])
      x1 <- data.frame(scores = "site", data = name1, label = x[,1], 
                       x[,c("WCS_ID", "TimeID", "Name", "Road", "ConPeriod", "ConMonths", "CatWalks", "Substrate", "Fencing")], 
                       x[,which(colnames(x) %in% paste(plot.type, 1:axes, sep = ""))])
      x1$Response <- factor(x1$data, levels = c("total", "succ", "fail"), labels = c("Total", "Successful", "Unsuccessful"))
      x1$WCS_ID <- factor(x1$WCS_ID, levels = c("FM106_WCS1", "FM106_WCS2", "FM106_WCS3", "FM106_WCS4", "FM106_WCS5", "FM106_WCS6", "FM106_WCS7", "FM106_WCS8", 
                                                "FM1847_WCS1", "FM1847_WCS2", "FM1847_WCS3", "FM1847_WCS4", "FM1847_WCS5", 
                                                "SH100_WCS1", "SH100_WCS2", "SH100_WCS3", "SH100_WCS3A", "SH100_WCS4"), 
                          labels = c("FM106 WCS1", "FM106 WCS2", "FM106 WCS3", "FM106 WCS4", "FM106 WCS5", "FM106 WCS6", "FM106 WCS7", "FM106 WCS8", 
                                     "FM1847 WCS1", "FM1847 WCS2", "FM1847 WCS3", "FM1847 WCS4", "FM1847 WCS5", 
                                     "SH100 WCS1", "SH100 WCS2", "SH100 WCS3", "SH100 WCS3A", "SH100 WCS4"))
      x1$RoadCon <- factor(paste(x1$Road, x1$ConPeriod, sep = "_"), 
                           levels = c("FM106_Con", "FM106_Postcon", "FM1847_Con", "FM1847_Postcon", "SH100_Con", "SH100_Postcon"), 
                           labels = c("FM106 Con", "FM106 Post-con", "FM1847 Con", "FM1847 Post-con", "SH100 Con", "SH100 Post-con"))
      return(x1)
    })
  }else if(score.type == "spec"){
    fs2 <- fs1[grep("spec", fs1)]
    scores1 <- pbapply::pblapply(1:length(fs2), function(i){
      name1 <- sub(paste(plot.type, "_", sep = ""), "", sub(paste("_", score.type, sep = ""), "", basename(fs::path_ext_remove(fs2[i]))))
      x <- openxlsx::read.xlsx(fs2[i])
      x1 <- data.frame(t(x[,-1]))
      colnames(x1) <- c(x[,1])
      x2 <- data.frame(scores = "spec", data = name1, label = rownames(x1),
                       apply(x1[,which(colnames(x1) %in% paste(plot.type, 1:axes, sep = ""))], 2, as.numeric))
      x2$Response <- factor(x2$data, levels = c("total", "succ", "fail"), labels = c("Total", "Successful", "Unsuccessful"))
      x2$Species <- factor(x2$label, 
                           levels = c("armadillo", "beaver", "bobcat", "cottontail", "coyote", "feral_hog", "fox_squirrel", "gray_squirrel", 
                                      "grey_fox", "jackrabbit", "javelina", "mexican_ground_squirrel", "nilgai", "nutria", "ocelot", "opossum", 
                                      "raccoon", "rodent", "striped_skunk", "unk_mammal", "weasel", "white_tailed_deer"), 
                           labels = c("Armadillo", "Beaver", "Bobcat", "Cottontail", "Coyote", "Feral hog", "Fox squirrel", "Gray squirrel", 
                                      "Grey fox", "Jackrabbit", "Javelina", "Mexican ground squirrel", "Nilgai", "Nutria", "Ocelot", "Opossum", 
                                      "Raccoon", "Rodent", "Striped skunk", "Unk mammal", "Weasel", "White-tailed deer"))
      x2$SpecNum <- factor(x2$label, 
                           levels = c("armadillo", "beaver", "bobcat", "cottontail", "coyote", "feral_hog", "fox_squirrel", "gray_squirrel", 
                                      "grey_fox", "jackrabbit", "javelina", "mexican_ground_squirrel", "nilgai", "nutria", "ocelot", "opossum", 
                                      "raccoon", "rodent", "striped_skunk", "unk_mammal", "weasel", "white_tailed_deer"), 
                           labels = seq(1:length(unique(x2$label))))
      return(x2)
    })
  }else if(score.type == "pred"){
    if(plot.type == "PCO"){
      stop("There should not be any covariate scores for plot.type == 'PCO'.")
    }else{
      fs2 <- fs1[grep("pred", fs1)]
    }
    scores1 <- pbapply::pblapply(1:length(fs2), function(i){
      name1 <- sub(paste(plot.type, "_", sep = ""), "", sub(paste("_", score.type, sep = ""), "", basename(fs::path_ext_remove(fs2[i]))))
      x <- openxlsx::read.xlsx(fs2[i])
      x1 <- data.frame(t(x[,-1]))
      colnames(x1) <- c(x[,1])
      x2 <- data.frame(scores = "pred", data = name1, label = rownames(x1), 
                       apply(x1[,which(colnames(x1) %in% paste(plot.type, 1:axes, sep = ""))], 2, as.numeric))
      x2$Response <- factor(x2$data, levels = c("total", "succ", "fail"), labels = c("Total", "Successful", "Unsuccessful"))
      
      str <- c("CatWalks.No", "CatWalks.Yes", "log(Openness_m+1)", "Substrate.Concrete", "Substrate.Dirt", "Substrate.Water")
      env <- c("PropNatural", "PropWoody", "PropWater", "Precip_total_mm")
      ant <- c("log(HumanActivity+1)", "log(DomesticActivity+1)", "log(LivestockActivity+1)", "log(PropBuilding+1)", "SpeedLimit", "VehicleTraffic")
      spat <- paste("s_PCNM", seq(1:11), sep = "")
      temp <- paste("t_PCNM", seq(1:22), sep = "")
      
      x2$CovType <- factor(ifelse(x2$label %in% str, "structural", ifelse(x2$label %in% env, "environmental", ifelse(x2$label %in% ant, "anthropogenic", ifelse(x2$label %in% spat, "spatial", "temporal")))), 
                           levels = c("structural", "environmental", "anthropogenic", "spatial", "temporal"))
      
      x2$Covariate <- factor(x2$label, 
                             levels = c("log(Openness_m+1)", "CatWalks.No", "CatWalks.Yes", "Substrate.Water", "Substrate.Concrete", "Substrate.Dirt", 
                                        "PropNatural", "PropWater", "PropWoody", "Precip_total_mm", 
                                        "SpeedLimit", "VehicleTraffic", "log(HumanActivity+1)", "log(DomesticActivity+1)", "log(LivestockActivity+1)", "log(PropBuilding+1)", 
                                        x2$label[x2$CovType == "spatial"], 
                                        x2$label[x2$CovType == "temporal"]), 
                             labels = c("Openness", "No CatWalk", "Catwalk", "Water Sub.", "Concrete Sub.", "Dirt Sub.", 
                                        "Prop Natural", "Prop Water", "Prop Woody", "Precipitation", 
                                        "Speed Limit", "Vehicle Traffic", "Human Activity", "Domestic Activity", "Livestock Activity", "Prop Buildings", 
                                        x2$label[x2$CovType == "spatial"], 
                                        x2$label[x2$CovType == "temporal"]))
      x2$CovNum <- factor(x2$label, 
                          levels = c("log(Openness_m+1)", "CatWalks.No", "CatWalks.Yes", "Substrate.Water", "Substrate.Concrete", "Substrate.Dirt", 
                                     "PropNatural", "PropWater", "PropWoody", "Precip_total_mm", 
                                     "SpeedLimit", "VehicleTraffic", "log(HumanActivity+1)", "log(DomesticActivity+1)", "log(LivestockActivity+1)", "log(PropBuilding+1)", 
                                     x2$label[x2$CovType == "spatial"], 
                                     x2$label[x2$CovType == "temporal"]), 
                          labels = c(1:length(str), 
                                     1:length(env), 
                                     1:length(ant), 
                                     sub("s_PCNM", "", x2$label[x2$CovType == "spatial"]), 
                                     sub("t_PCNM", "", x2$label[x2$CovType == "temporal"])))
      
      return(x2)
    })
  }else{
    stop("You must choose a valid scores type. Choose one of c('site', 'spec', 'pred').")
  }
  scores2 <- do.call(rbind, scores1)
  return(scores2)
  rm(fs1, fs2, scores1, scores2)
  #rm(files, plot.type, score.type, axes)
}

## Function for making plots
plotOrdination <- function(ds, plot.type, response, X, Y, spec.lab, cov.lab, include.lines = TRUE, average = FALSE, size.text, size.lab, size.point){
  #ds <- dbrda.full
  #plot.type <- "dbRDA"
  #response <- "Road_WCS"
  #X <- "PCO1"
  #Y <- "PCO2"
  #X <- "dbRDA1"
  #Y <- "dbRDA2"
  #spec.lab <- "SpecNum"
  #cov.lab <- "CovNum"
  #include.lines <- TRUE
  #average <- FALSE
  #size.text <- 3
  #size.lab <- 9
  #size.point <- 3
  
  # Extract the site and species scores
  sites <- ds$sitescores
  specs <- ds$specscores
  
  spec_X <- sub("av_", "", X)
  spec_Y <- sub("av_", "", Y)
  
  # Load the ggplot2 package
  require(ggplot2)
  
  # The theme to make the graph look how we want it to
  theme.plot <- theme(text = element_text(family = "serif")) + 
    theme(plot.title = element_text(hjust = 0.5, size = size.lab*1.50, margin = margin(b = 0.5, unit = "inch")), 
          plot.background = element_blank(), 
          plot.margin = unit(c(.1,.1,.1,.1), "inch")) +
    theme(axis.ticks = element_line(color = NA, linewidth = 1, linetype = "solid"), 
          axis.line = element_line(color = NA, linewidth = .1, linetype = "solid"), 
          axis.title=element_text(size=size.lab*1.15, margin = margin(t = 0.25, unit="inch")),  
          axis.title.x = element_text(vjust = -1.0, hjust = 0.50), 
          axis.title.y = element_text(angle = 90, hjust = 0.50, vjust = 2.0), 
          axis.text = element_text(size = size.lab, color = "black"), 
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.0), 
          axis.text.y = element_text(angle = 0, hjust = 0)) + 
    theme(panel.border = element_rect(fill = NA, color = "black"), 
          panel.background = element_rect(fill = NA, color = NA), 
          panel.grid.major = element_line(color = NA), 
          panel.grid.major.x = element_line(color = NA), 
          panel.grid.major.y = element_line(color = NA), 
          panel.grid.minor = element_line(color = NA), 
          panel.grid.minor.x = element_line(color = NA), 
          panel.grid.minor.y = element_line(color = NA)) +  
    theme(legend.margin=margin(c(0.15,0.15,0.15,0.15), unit = "inch"), 
          legend.background = element_rect(fill = NA, color = NA), 
          legend.text=element_text(size = size.lab), 
          legend.title=element_text(size = size.lab), 
          legend.position = "right", 
          #legend.key = element_rect(color = "black", fill = NA), 
          legend.key.height = unit(0.25,"inch"), 
          legend.key.width = unit(0.25, "inch")) + 
    theme(strip.background = element_rect(fill = "gray85", color = "black"), 
          strip.text = element_text(size = size.lab), 
          strip.text.x = element_text(margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")), 
          strip.text.y = element_text(angle = -90, margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")))
  
  # Adjust fill colors and shapes
  if(response == "WCS_ID"){
    COL <- "WCS_ID"
    SHAPE <- NULL
    scale_color <- NULL
    #scale_color <- scale_color_manual("Site", values = c())
    scale_shape <- NULL
  }else if(response == "Road"){
    COL <- "Road"
    SHAPE <- NULL
    scale_color <- scale_color_manual("Road", values = c("red", "yellow", "blue"))
    scale_shape <- NULL
  }else if(response == "WCS_Con"){
    COL <- "WCS_ID"
    SHAPE <- "ConPeriod"
    scale_color <- NULL
    #scale_color <- scale_color_manual("Site", values = c())
    scale_shape <- scale_shape_manual("Con Period", values = c(15, 16))
  }else if(response == "Road_WCS"){
    COL <- "Name"
    SHAPE <- "Road"
    scale_color <- scale_color_manual("WCS", values = c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255", "#EE9D28"))
    scale_shape <- scale_shape_manual("Road", values = c(15, 16, 17))
  }else if(response == "Road_Con"){
    COL <- "Name"
    SHAPE <- "RoadCon"
    scale_color <- scale_color_manual("WCS", values = RColorBrewer::brewer.pal(9, "Set1"))
    scale_shape <- scale_shape_manual("Road Con", values = c(0,15, 1,16, 2,17))
  }else{
    stop("You must choose an appropriate combination of plots. Choose one of c('WCS_ID', 'Road', 'WCS_Con', 'Road_WCS')")
  }
  
  # If averaging the values
  if(isTRUE(average)){
    if(is.null(SHAPE)){
      sites2 <- dplyr::summarise(dplyr::group_by(sites, Response,  !!sym(COL)), X = mean(!!sym(X)), Y = mean(!!sym(Y)))
    }else{
      sites2 <- dplyr::summarise(dplyr::group_by(sites, Response, !!sym(COL), !!sym(SHAPE)), X = mean(!!sym(X)), Y = mean(!!sym(Y)))
    }
  }else{
    sites2 <- sites
    colnames(sites2)[colnames(sites2) %in% c(X, Y)] <- c("X", "Y")
  }
  
  # Set the limits of the plots based on the data
  limits <- list(xlims = c(min = plyr::round_any(min(sites2[,"X"]), accuracy = 5, f = floor), 
                           max = plyr::round_any(max(sites2[,"X"]), accuracy = 5, f = ceiling)), 
                 ylims = c(min = plyr::round_any(min(sites2[,"Y"]), accuracy = 5, f = floor), 
                           max = plyr::round_any(max(sites2[,"Y"]), accuracy = 5, f = ceiling)))
  limits <- list(xlims = c(min = -55, max = 75), ylims = c(min = -50, max = 60))
  
  # Standardize the correlations to the maximum and minimum values
  ## Site scores
  site_split <- split(sites2, f = sites2$Response)
  site_split2 <- site_split[which(sapply(site_split, nrow) > 0)]
  ## Species scores
  spec_split <- split(specs, f = specs$Response)[names(site_split)]
  spec_split1 <- spec_split[which(sapply(spec_split, nrow) > 0)]
  spec_split2 <- lapply(1:length(spec_split1), function(i){
    out <- data.frame(spec_split1[[i]][,c(spec.lab, "Response", spec_X, spec_Y)], 
                      ifelse(spec_split1[[i]][,spec_X] >= 0, spec_split1[[i]][,spec_X] * limits$xlims["max"], spec_split1[[i]][,spec_X] * -limits$xlims["min"]), 
                      ifelse(spec_split1[[i]][,spec_Y] >= 0, spec_split1[[i]][,spec_Y] * limits$ylims["max"], spec_split1[[i]][,spec_Y] * -limits$ylims["min"]))
    colnames(out) <- c("Species", "Response", paste(c(spec_X,spec_Y), "orig", sep = "_"), "X", "Y")
    return(out)
  })
  names(spec_split2) <- names(spec_split1)
  ## Covariate scores
  if(plot.type == "dbRDA"){
    covs <- ds$covscores
    cov_split1 <- split(covs, f = ~Response)
    cov_split2 <- lapply(cov_split1, split, f = ~CovType)
    cov_split3 <- lapply(1:length(cov_split2), function(i){
      x <- cov_split2[[i]]
      
      x1 <- lapply(1:length(x), function(j){
        y <- x[[j]]
        
        y1 <- data.frame(y[,c(cov.lab, "Response", spec_X, spec_Y)], 
                         ifelse(y[,spec_X] >= 0, y[,spec_X] * limits$xlims["max"], y[,spec_X] * -limits$xlims["min"]), 
                         ifelse(y[,spec_Y] >= 0, y[,spec_Y] * limits$ylims["max"], y[,spec_Y] * -limits$ylims["min"]))
        colnames(y1) <- c("Covariate", "Response", paste(c(spec_X, spec_Y), "orig", sep = "_"), "X", "Y")
        return(y1)
      })
      names(x1) <- names(x)
      return(x1)
    })
    names(cov_split3) <- names(cov_split2)
    
    ds_in <- lapply(names(site_split2), function(x){
      c(site = list(site_split2[[x]]), spec = list(spec_split2[[x]]), cov_split3[[x]])
    })
    names(ds_in) <- names(site_split2)
  }else if(plot.type == "PCO"){
    ds_in <- lapply(names(site_split2), function(x){
      c(site = list(site_split2[[x]]), spec = list(spec_split2))
    })
    names(ds_in) <- names(site_split2)
  }else{
    stop("You must choose an appropriate plot type. Choose one of c('PCO','dbRDA').")
  }
  
  plot.out <- lapply(1:length(ds_in), function(i){
    x <- ds_in[[i]]
    
    spec.col <- "grey30"
    cov.col <- "grey30"
    
    # Labels and lines for species scores
    mult <- 1.60
    maxover <- 20
    gtext.spec <- ggrepel::geom_text_repel(data = x$spec, mapping = aes(x = X * mult, y = Y * mult, label = Species), 
                                           size = size.text, color = spec.col, max.overlaps = maxover)
    
    # Labels and lines for covariate scores
    if(plot.type == "dbRDA"){
      # Structural
      gtext.str <- ggrepel::geom_text_repel(data = x$structural, mapping = aes(x = X * mult, y = Y * mult, label = Covariate), size = size.text, color = cov.col, max.overlaps = maxover)
      gline.str <- geom_segment(data = x$structural, mapping = aes(x = 0, y = 0, xend = X * mult, yend = Y * mult), linetype = 1, color = cov.col, arrow = arrow(length = unit(0.10, "inches"), type = "closed"))
      # Environmental
      gtext.env <- ggrepel::geom_text_repel(data = x$environmental, mapping = aes(x = X * mult, y = Y * mult, label = Covariate), size = size.text, color = cov.col, max.overlaps = maxover)
      gline.env <- geom_segment(data = x$environmental, mapping = aes(x = 0, y = 0, xend = X * mult, yend = Y * mult), linetype = 1, color = cov.col, arrow = arrow(length = unit(0.10, "inches"), type = "closed"))
      # Anthropogenic
      gtext.ant <- ggrepel::geom_text_repel(data = x$anthropogenic, mapping = aes(x = X * mult, y = Y * mult, label = Covariate), size = size.text, color = cov.col, max.overlaps = maxover)
      gline.ant <- geom_segment(data = x$anthropogenic, mapping = aes(x = 0, y = 0, xend = X * mult, yend = Y * mult), linetype = 1, color = cov.col, arrow = arrow(length = unit(0.10, "inches"), type = "closed"))
      # Spatial 
      gtext.spat <- ggrepel::geom_text_repel(data = x$spatial, mapping = aes(x = X * mult, y = Y * mult, label = Covariate), size = size.text, color = cov.col, max.overlaps = maxover)
      gline.spat <- geom_segment(data = x$spatial, mapping = aes(x = 0, y = 0, xend = X * mult, yend = Y * mult), linetype = 1, color = cov.col, arrow = arrow(length = unit(0.10, "inches"), type = "closed"))
      # Temporal
      gtext.temp <- ggrepel::geom_text_repel(data = x$temporal, mapping = aes(x = X * mult, y = Y * mult, label = Covariate), size = size.text, color = cov.col, max.overlaps = maxover)
      gline.temp <- geom_segment(data = x$temporal, mapping = aes(x = 0, y = 0, xend = X * mult, yend = Y * mult), linetype = 1, color = cov.col, arrow = arrow(length = unit(0.10, "inches"), type = "closed"))
    }else{
      gtext.str <- gline.str <- gtext.env <- gline.env <- gtext.ant <- gline.ant <- gtext.spat <- gline.spat <- gtext.temp <- gline.temp <- NULL
    }
    
    # Should lines for the correlations between species scores/covariate scores and site scores be included
    if(isTRUE(include.lines)){
      gline.spec <- geom_segment(data = x$spec, mapping = aes(x = 0, y = 0, xend = X * mult, yend = Y * mult), linetype = 2, color = spec.col, arrow = arrow(length = unit(0.10, "inches"), type = "closed"))
    }else{
      gline.spec <- gline.str <- gline.env <- gline.ant <- gline.spat <- gline.temp <- NULL
    }
    
    # If shape is being included
    if(is.null(SHAPE)){
      gpoint <- geom_point(data = x$site, mapping = aes(x = X, y = Y, col = !!sym(COL)), size = size.point)
    }else{
      gpoint <- geom_point(data = x$site, mapping = aes(x = X, y = Y, col = !!sym(COL), shape = !!sym(SHAPE)), size = size.point)
    }
    
    # Define the break distance
    break.dist <- 15
    
    hline <- geom_hline(yintercept = 0, color = "black")
    vline <- geom_vline(xintercept = 0, color = "black")
    scale_x <- scale_x_continuous(name = spec_X, breaks = seq(-120,120,break.dist), limits = c(-120,120))
    scale_y <- scale_y_continuous(name = spec_Y, breaks = seq(-120,120,break.dist), limits = c(-120,120))
    
    out1 <- list(
         spec = ggplot(data = x$site) + 
           hline + vline + gline.spec + gtext.spec + 
           scale_x + scale_y + 
           scale_color + scale_shape + coord_cartesian(xlim = limits$xlims, ylim = limits$ylim) + theme.plot, 
         str = ggplot(data = x$site) + 
           hline + vline + gline.str + gtext.str + 
           scale_x + scale_y + 
           scale_color + scale_shape + coord_cartesian(xlim = limits$xlims, ylim = limits$ylim) + theme.plot, 
         env = ggplot(data = x$site) + 
           hline + vline + gline.env + gtext.env + 
           scale_x + scale_y + 
           scale_color + scale_shape + coord_cartesian(xlim = limits$xlims, ylim = limits$ylim) + theme.plot, 
         ant = ggplot(data = x$site) + 
           hline + vline + gline.ant + gtext.ant + 
           scale_x + scale_y + 
           scale_color + scale_shape + coord_cartesian(xlim = limits$xlims, ylim = limits$ylim) + theme.plot, 
         spat = ggplot(data = x$site) + 
           hline + vline + gline.spat + gtext.spat + 
           scale_x + scale_y + 
           scale_color + scale_shape + coord_cartesian(xlim = limits$xlims, ylim = limits$ylim) + theme.plot, 
         temp = ggplot(data = x$site) + 
           hline + vline + gline.temp + gtext.temp + 
           scale_x + scale_y + 
           scale_color + scale_shape + coord_cartesian(xlim = limits$xlims, ylim = limits$ylim) + theme.plot, 
         points = ggplot(data = x$site) + 
           hline + vline + gpoint + 
           scale_x + scale_y + 
           scale_color + scale_shape + coord_cartesian(xlim = limits$xlims, ylim = limits$ylim) + theme.plot
         )
    out2 <- ggpubr::ggarrange(plotlist = out1, ncol = 3, nrow = 3, labels = "AUTO", font.label = list(family = "serif"), common.legend = TRUE) + ggpubr::bgcolor("white")
    out2
    return(out2)
    #rm(x, spec.col, cov.col, gtext.spec, gline.spec, gtext.str, gline.str, gtext.env, gline.env, gtext.ant, gline.ant, gtext.spat, gline.spat, gtext.temp, gline.temp, gpoint, out1, out2)
    #rm(i)
  })
  names(plot.out) <- names(ds_in)
  plot.out[[1]]
  
  if(length(plot.out) == 1){
    plot.out <- plot.out[[1]]
  }else{
    message("The output will be a list because multiple plots are created")
  }
  
  return(plot.out)
  rm(sites, specs, spec_X, spec_Y, 
     theme.plot, COL, SHAPE, scale_color, scale_shape, 
     sites2, limits, site_split, site_split2, spec_split, spec_split1, spec_split2, covs, cov_split1, cov_split2, cov_split3, 
     ds_in, plot.out)
  #rm(ds, plot.type, response, X, Y, spec.lab, cov.lab, include.lines, average, size.text, size.lab, size.point)
}


# Used in Part 3: Variance partitioning ####
## Function for calculation of fractional variation
varpart5 <- function(parts, Xs, type){
  #parts <- vp2[[3]]
  Xs <- c("str", "env", "ant", "spat", "temp")
  #type <- "AdjustedR2"
  
  # Define the sets
  names(Xs) <- paste("X", 1:length(Xs), sep = "")
  Xs2 <- data.frame(t(Xs))
  
  # The original output RDA values
  ps1 <- as.list(parts[,type])
  names(ps1) <- parts$RDA
  
  # Intermediary functions
  Is <- with(ps1, list(A29 - A25, 
                       A28 - A24, 
                       A27 - A23, 
                       A26 - A22, 
                       A28 - A21, 
                       A27 - A20, 
                       A26 - A19, 
                       A27 - A18, 
                       A27 - A17, 
                       A26 - A16, 
                       A21 - A15, 
                       A20 - A14, 
                       A19 - A13, 
                       A18 - A12, 
                       A17 - A11, 
                       A16 - A10, 
                       A18 - A9, 
                       A17 - A8, 
                       A16 - A7, 
                       A16 - A6, 
                       A9 - A5, 
                       A8 - A4, 
                       A7 - A3, 
                       A6 - A2, 
                       A6 - A1))
  names(Is) <- paste("I", seq(1,length(Is)), sep = "")
  
  # Computing the proportion of variance explained
  ## Single set and two set partitions
  ps2 <- with(ps1, with(Is, list(A31 - A30, 
                                 A31 - A29, 
                                 A31 - A28, 
                                 A31 - A27, 
                                 A31 - A26, 
                                 I1 - A1, 
                                 I2 - A1, 
                                 I3 - A1, 
                                 I4 - A1, 
                                 I5 - A2, 
                                 I6 - A2, 
                                 I7 - A2, 
                                 I8 - A3, 
                                 I9 - A3, 
                                 I10 - A4)))
  names(ps2) <- paste("B", seq(1,length(ps2)), sep = "")
  ## 3 set partitions
  ps2b <- with(ps1, with(Is, with(ps2, list(I11 - I1 - B7, 
                                            I12 - I1 - I8, 
                                            I13 - I1 - B9, 
                                            I14 - I2 - B8, 
                                            I15 - I2 - B9, 
                                            I16 - I3 - B9, 
                                            I17 - I5 - B11, 
                                            I18 - I5 - B12, 
                                            I19 - I6 - B12, 
                                            I20 - I8 - B14))))
  names(ps2b) <- paste("B", seq(length(ps2)+1,length(ps2) + length(ps2b)), sep = "")
  ps2 <- c(ps2, ps2b)
  ## 4 set partitions
  ps2c <- with(ps1, with(Is, with(ps2, list(I21 - I11 - B8 - B17 - B19, 
                                            I22 - I11 - B9 - B18 - B20, 
                                            I23 - I12 - B9 - B20 - B21, 
                                            I24 - I17 - B12 - B23 - B24, 
                                            I25 - I17 - B12 - B23 - B24))))
  names(ps2c) <- paste("B", seq(length(ps2) + 1, length(ps2) + length(ps2c)), sep = "")
  ps2 <- c(ps2, ps2c)
  ## 5 set partitions
  ps2d <- with(ps1, with(Is, with(ps2, list(A1 - I11 - B8 - B9 - B17 - B18 - B19 - B20 - B21 - B26 - B27 - B28 - B29, 
                                            1 - A31))))
  names(ps2d) <- paste("B", seq(length(ps2) + 1, length(ps2) + length(ps2d)), sep = "")
  ps2 <- c(ps2, ps2d)
  
  # The groups that each set belongs to
  Bs <- with(Xs2, list(X1, X2, X3, X4, X5, 
                       c(X1, X2), c(X1, X3), c(X1, X4), c(X1, X5), 
                       c(X2, X3), c(X2, X4), c(X2, X5), 
                       c(X3, X4), c(X3, X5), 
                       c(X4, X5), 
                       c(X1, X2, X3), c(X1, X2, X4), c(X1, X2, X5), c(X1, X3, X4), c(X1, X3, X5), c(X1, X4, X5), 
                       c(X2, X3, X4), c(X2, X3, X5), c(X2, X4, X5), c(X3, X4, X5), 
                       c(X1, X2, X3, X4), c(X1, X2, X3, X5), c(X1, X2, X4, X5), c(X1, X3, X4, X5), c(X2, X3, X4, X5), 
                       c(X1, X2, X3, X4, X5), 
                       "residual"
  )
  )
  
  # Output as a data.frame
  out <- data.frame(name = names(ps2), 
                    no.sets = c(sapply(Bs[-length(Bs)], length), 0), 
                    partition = sapply(Bs, function(x){paste(x, collapse = " ")}), 
                    type = do.call(c, ps2))
  colnames(out) <- c("name", "no.sets", "partition", type)
  row.names(out) <- NULL
  
  return(out)
}

## Function for plotting
varpartPlot <- function(data, size.text, size.lab){
  #data <- vp4$Total
  #size.text <- 3
  #size.lab <- 8
  
  library(ggplot2)
  
  ## Modify the theme for the plot
  theme.plot <- theme(text = element_text(family = "serif")) + 
    theme(plot.title = element_text(hjust = 0.5, size = size.lab*1.50, margin = margin(b = 0.5, unit = "inch")), 
          plot.background = element_blank(), 
          plot.margin = unit(c(.1,.1,.1,.1), "inch")) +
    theme(axis.ticks = element_blank(), 
          axis.line = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank()) + 
    theme(panel.border = element_rect(fill = NA, color = "black"), 
          panel.background = element_rect(fill = NA, color = NA), 
          panel.grid = element_blank()) + 
    theme(legend.margin=margin(c(0.15,0.15,0.15,0.15), unit = "inch"), 
          legend.background = element_rect(fill = NA, color = NA), 
          legend.text=element_text(size = size.lab), 
          legend.title=element_text(size = size.lab), 
          legend.position = "right", 
          #legend.key = element_rect(color = "black", fill = NA), 
          legend.key.height = unit(0.25,"inch"), 
          legend.key.width = unit(0.25, "inch")) + 
    theme(strip.background = element_rect(fill = "gray85", color = "black"), 
          strip.text = element_text(size = size.lab), 
          strip.text.x = element_text(margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")), 
          strip.text.y = element_text(angle = -90, margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")))
  
  ## Create a circle for plotting in ggplot
  circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  
  # Create circles
  cp <- matrix(c(0.05, 0, 0.95, 0, 0.5, -sqrt(3/4)), nrow = 3, ncol = 2, byrow = TRUE)
  dimnames(cp) <- list(c("c1", "c2", "c3"), c("x", "y"))
  cp
  
  diam <- 0.725 * 2
  
  circ1 <- circleFun(center = cp[1,], diameter = diam, npoints = 100)
  circ2 <- circleFun(center = cp[2,], diameter = diam, npoints = 100)
  circ3 <- circleFun(center = cp[3,], diameter = diam, npoints = 100)
  
  # Create labels
  l.name <- data.frame(x = c(cp[,1]), 
                       y = c(cp[1:2,2]+diam/2, cp[3,2]-diam/2), 
                       label = c("Structural", "Environmental", "Anthropogenic"))
  
  l.value <- data.frame(rbind(cp, c("x" = 0.5, "y" = 0.0), c("x" = 0.81, "y" = -0.4930127), c("x" = 0.19, "y" = -0.4930127), c("x" = 0.5, "y" = -0.2886751), c(1.3, -1.65)), 
                        value = data[,4], 
                        label = c(paste(round(data[-8,4] * 100, digits = 2), "%", sep = ""), paste("Residual = ", round(data[8,4] * 100, digits = 2), "%", sep = "")))
  
  #l.value <- data.frame(rbind(cp, colMeans(cp[1:2,]), colMeans(cp[2:3,]), colMeans(cp[c(1,3),]), colMeans(cp), c(1.5, -1.65)), 
  #                      value = data[,4], 
  #                      label = c(paste(round(data[-8,4] * 100, digits = 2), "%", sep = ""), paste("Residual = ", round(data[8,4] * 100, digits = 2), "%", sep = "")))
  
  # Create the plot
  plot.out <- ggplot(mapping = aes(x, y)) + 
    geom_path(data = circ1) + 
    geom_path(data = circ2) + 
    geom_path(data = circ3) + 
    geom_text(data = l.value, mapping = aes(x = x, y = y, label = label), size = size.text, size.unit = "mm") + 
    geom_label(data = l.name, mapping = aes(x = x, y = y, label = label), size = size.text, size.unit = "mm") + 
    coord_fixed() + 
    theme.plot
  plot.out
  
  return(plot.out)
  rm(theme.plot, circleFun, cp, diam, circ1, circ2, circ3, l.name, l.value, plot.out)
  #rm(data, size.lab, size.text)
}


# Used in Part 4: Regression Analyses ####
## Calculation of a distance matrix  (function modified from Paul and Anderson 2013)
## distance()   # ABOVE

## PCO analysis (function modified from Paul and Anderson 2013)
### NOTE: PCO must be done on dissimilarity matrices so similarity matrices must be converted before calculation
### 1 - (similarity matrix) = (dissimilarity matrix)
pco <- function(D, varplot=T){
  #D <- Dme  # Distance matrix
  #varplot <- TRUE # Plot the proportion of explained variance
  
  # Number of dimensions of the distance matrix
  N <- dim(D)[1]
  
  # Computations for the centered Gower's Matrix (G)
  ## Do the initial conversion
  A <- -0.5 * D^2 
  
  # Center the A matrix
  centre.matrix <- function(A) {
    n <- dim(A)[1]
    return(A - apply(A, 1, mean) %o% rep(1, n) - rep(1, n) %o% apply(A, 2, mean) + mean(A)) 
  }
  
  # Center the matrix
  G <- centre.matrix(A) 
  # Eigenvalue decomposition analysis on the centered Gower's matrix
  eigenG <- eigen(G, symmetric=T)
  # Plot the proportion of explained variance
  if(varplot==T){
    propnvar <- c(0, cumsum(eigenG$val^2))/sum(eigenG$val^2)
    plot(0:length(eigenG$val), propnvar,type="l",
         ylim = c(0,1), xlab="# pco's", ylab="Propn of variation explained")
    abline(h = 1, v = seq(0, length(propnvar), 5), lty = 3)
  }
  
  Q <- eigenG$vectors
  L <- eigenG$values
  pos.neg <- L >= 0.0
  PCO <- Q[,1:N] %*% sqrt(diag(abs(L[1:N])))
  
  return(list(vectors = PCO, values = L))
}

## Function for data loading for cleanliness purposes
prepData <- function(files, type, mo = "all", site = "all", drop.pred = NULL, drop_site = FALSE, drop_mo = FALSE, autocorrelation = FALSE){
  #files <- files_reg
  #type <- "total"
  #mo <- "all"
  #site <- "all"
  #drop.pred <- c("VehicleTraffic")
  #autocorrelation <- FALSE
  #drop_site <- TRUE  # Should the function iteratively drop a single site from the initial processing
  #drop_mo <- FALSE    # Should the function iteratively drop a single month from the initial processing
  
  # Get only the Response variables
  Ys <- files[!grepl("cov_", files)]
  # Get only the predictor variables
  Xs <- files[grepl("cov_", files)]
  
  # Read in the files
  if(type == "total"){
    Y.all <- openxlsx::read.xlsx(Ys[grep(type, Ys, ignore.case = TRUE)], rowNames = TRUE)
    X.all <- openxlsx::read.xlsx(Xs[grep(type, Xs, ignore.case = TRUE)], rowNames = TRUE, rows = 2:349)
  }else if(type == "succ"){
    Y.all <- openxlsx::read.xlsx(Ys[grep(type, Ys, ignore.case = TRUE)], rowNames = TRUE)
    X.all <- openxlsx::read.xlsx(Xs[grep(type, Xs, ignore.case = TRUE)], rowNames = TRUE, rows = 2:309)
  }else if(type == "fail"){
    Y.all <- openxlsx::read.xlsx(Ys[grep(type, Ys, ignore.case = TRUE)], rowNames = TRUE)
    X.all <- openxlsx::read.xlsx(Xs[grep(type, Xs, ignore.case = TRUE)], rowNames = TRUE, rows = 2:309)
  }else{
    stop("You must choose an appropriate type. Choose one of c('total','succ','fail').")
  }
  
  # Select specific months or sites
  if(mo == "all"){
    mos <- sort(unique(Y.all$ConMonths))
  }else{
    mos <- mo
  }
  if(site == "all"){
    sites <- sort(unique(Y.all$WCS_ID))
  }else{
    sites <- site
  }
  
  # Determine which variables to include
  if(isTRUE(autocorrelation)){
    cols <- c("log(Openness_m+1)", "CatWalks.Yes", "Substrate.Concrete", "Substrate.Dirt", "log(Fencing_m+1)", 
              "PropNatural", "PropWater", "PropWoody", "Precip_total_mm", 
              "SpeedLimit", "VehicleTraffic", "log(HumanActivity+1)", "log(DomesticActivity+1)", "log(LivestockActivity+1)", "log(PropBuilding+1)", 
              paste("s_PCNM", seq(1:11), sep = ""), paste("t_PCNM", seq(1:22), sep = ""))
  }else if(isFALSE(autocorrelation)){
    cols <- c("log(Openness_m+1)", "CatWalks.Yes", "Substrate.Concrete", "Substrate.Dirt", "log(Fencing_m+1)", 
              "PropNatural", "PropWater", "PropWoody", "Precip_total_mm", 
              "SpeedLimit", "VehicleTraffic", "log(HumanActivity+1)", "log(DomesticActivity+1)", "log(LivestockActivity+1)", "log(PropBuilding+1)", 
              "ConMonths", "UTM_X", "UTM_Y")
  }else{
    stop("autocorrelation must be logical.")
  }
  
  if(!is.null(drop.pred)){
    message(paste(length(drop.pred), " predictor(s) is being dropped from the model.\nIt(They) is(are): ", paste(drop.pred, collapse = ", "), sep = ""))
    cols <- cols[!(cols %in% drop.pred)]
  }
  
  Y.full <- Y.all[Y.all$ConMonths %in% mos & Y.all$WCS_ID %in% sites,]
  X.full <- X.all[X.all$ConMonths %in% mos & X.all$WCS_ID %in% sites,]
  
  Y <- as.matrix(Y.full[,c(1:22)])
  X <- as.matrix(X.full[,colnames(X.all) %in% cols])
  D <- distance(Y = Y, measure = "BC", trans = "none", adj = 0.03)
  N <- dim(D)[1]
  
  full <- list(Y = Y, X = X, D = D, N = N)
  
  if(isTRUE(drop_site)){
    sdrop <- pbapply::pblapply(sites, function(s){
      Y <- as.matrix(Y.full[Y.full$WCS_ID != s, c(1:22)])
      X <- as.matrix(X.full[X.full$WCS_ID != s, colnames(X.all) %in% cols])
      Y.drop <- as.matrix(Y.full[Y.full$WCS_ID == s, c(1:22)])
      X.drop <- as.matrix(X.full[X.full$WCS_ID == s, colnames(X.all) %in% cols])
      D <- distance(Y = Y, measure = "BC", trans = "none", adj = 0.03)
      N <- dim(D)[1]
      list(Y = Y, X = X, D = D, N = N, Y.drop = Y.drop, X.drop = X.drop)
    })
    names(sdrop) <- paste("drop_", sites, sep = "")
  }else{
    sdrop <- NULL
  }
  
  if(isTRUE(drop_mo)){
    tdrop <- pbapply::pblapply(mos, function(m){
      Y <- as.matrix(Y.full[Y.full$ConMonths != m, c(1:22)])
      X <- as.matrix(X.full[X.full$ConMonths != m, colnames(X.all) %in% cols])
      Y.drop <- as.matrix(Y.full[Y.full$ConMonths == m, c(1:22)])
      X.drop <- as.matrix(X.full[X.full$ConMonths == m, colnames(X.all) %in% cols])
      D <- distance(Y = Y, measure = "BC", trans = "none", adj = 0.03)
      N <- dim(D)[1]
      list(Y = Y, X = X, D = D, N = N, Y.drop = Y.drop, X.drop = X.drop)
    })
    names(tdrop) <- paste("drop_", mos, sep = "")
  }else{
    tdrop <- NULL
  }
  
  out <- c(list(full = full), sdrop, tdrop)
  
  return(out)
  rm(Ys, Xs, Y.all, X.all, Y.full, X.full, Y, X, D, N, full, mos, sites, cols, sdrop, tdrop, out)
  #rm(files, type, mo, site, drop_mo, drop_site, autocorrelation)
}

## Function for running PCO analyses
pcoOut <- function(input, out.diag = FALSE){
  #input <- total
  
  ### Conduct PCO Analysis
  PCO.result <- pco(D = input$D, varplot = TRUE)
  
  ### Extract eigenvectors
  PCO <- PCO.result$vectors
  colnames(PCO) <- paste("PCO", 1:ncol(PCO), sep = "")
  row.names(PCO) <- row.names(input$Y)
  
  ### Convert the PCO matrix to a data.frame
  PCO.df <- data.frame(site = row.names(input$Y), PCO)
  pos.neg <- PCO.result$values >= 0.0
  
  if(isTRUE(out.diag)){
    print("Percent of variability explained by each PCO axis:")
    print(round(100 * PCO.result$values/sum(PCO.result$values), digits = 2))
  }
  
  out <- list(values = PCO.result$values, vectors = PCO, df = PCO.df, pos.neg = pos.neg)
  return(out)
  rm(PCO.result, PCO, PCO.df, pos.neg)
}

## Function for running regression analyses
runRegression <- function(input, PCO, method = "hand"){
  #input <- in.total.cor$full
  #PCO <- pco.total.cor$full
  #method <- "hand"
  #method <- "lm"
  
  # Calcluate the means of the input X's
  X.mean <- with(input, apply(X, 2, mean))
  
  # Center the X's around their mean
  Xc <- with(input, X - rep(1, N) %o% X.mean)
  
  if(method == "hand"){
    # Compute the regression coefficients
    ## First compute the hat matrix
    hat.matrix <- Xc %*% solve(t(Xc) %*% Xc) %*% t(Xc)
    hat.diag <- diag(hat.matrix)
    
    ### NOTE: These are the beta coefficients for each predictor in each experimental unit given a particular "response" variable from PCO
    ## Compute the beta coefficients
    pco.beta <- solve(t(Xc) %*% Xc) %*% t(Xc) %*% PCO$vectors
    
    # Compute the fitted PCO values
    ## NOTE: These are the fitted values of this regression for all PCO axes and is calculated by multiplying and summing the values of X with their appropriate beta coefficient
    ## One can estimate the fitted values from the beta coefficients or the hat matrix or neither
    pco.fitted <- Xc %*% pco.beta
    #pco.fitted <- hat.matrix %*% PCO$vectors  # Using the hat matrix
    #pco.fitted <- Xc %*% solve(t(Xc) %*% Xc) %*% t(Xc) %*% PCO$vectors  # Using no shorthand 
    
    # Compute the residuals of the PCO regression
    pco.resid <- PCO$vectors - pco.fitted
  }else if(method == "lm"){
    ## Create a linear model using R's linear modelling framework
    model1 <- lm(PCO$vectors ~ Xc)
    
    ## Compute the hat matrix
    message("hatvalues from a linear model are only the diagnonal of the hat matrix...which is unfortunate")
    hat.matrix <- stats::hatvalues(model1)
    hat.diag <- hat.matrix
    
    ## Compute the beta coefficients
    pco.beta <- stats::coef(model1)[-1,]
    row.names(pco.beta) <- colnames(Xc)
    
    # Compute the fitted PCO values
    pco.fitted <- stats::fitted(model1)
    
    # Compute the residuals of the PCO regression
    pco.resid <- stats::residuals(model1)
  }else{
    stop("You must specify a method. Choose either 'hand' or 'lm'. ")
  }
  
  # Degrees of freedom
  ## Number of samples
  n <- input$N
  ## Number of parameters in the model
  k <- dim(Xc)[2]
  ## Total df
  df.total <- n - 1
  ## Error df
  df.resid <- n - k - 1
  
  # Calculate Sums of Squares (Equations from Kutner et al. pgs. 204-205)
  J <- matrix(1, nrow = n, ncol = n)
  ## Total Sums of Squares
  SST <- diag(t(PCO$vectors) %*% PCO$vectors - 1/n * t(PCO$vectors) %*% J %*% PCO$vectors)
  ## Error Sums of Squares
  SSE <- diag(t(pco.resid) %*% pco.resid)
  ## Regression Sums of Squares
  if(any(is.na(pco.beta))){
    # SSR can also be calculated as the difference between SST and SSE
    message("SSR was calculated by substraction instead matrix multiplication due to NA values in one of the predictors")
    SSR <- SST - SSE
  }else{
    SSR <- diag(t(pco.beta) %*% t(Xc) %*% PCO$vectors - 1/n * t(PCO$vectors) %*% J %*% PCO$vectors)
  }
  
  # Compile the model diagnostics into a single data frame
  ## Calculate R^2, adjusted R^2, AIC, AICc, and the PRESS statistic
  diagn <- data.frame(SST = SST, SSR = SSR, SSE = SSE, 
                      n = input$N, k = dim(Xc)[2], df.total = n - 1, df.resid = n - k - 1, 
                      R2 = 1 - SSE/SST, adjustedR2 = 1 - (SSE/SST)*(df.total/df.resid), AIC = n * log(SSE/n) + 2*k, AICc = (n * log(SSE/n) + 2*k) + (2*k*(k+1))/(n-k-1), 
                      PRESS = apply(pco.resid, 2, function(x){sum((x/(1 - hat.diag))^2)}), 
                      pos.eigen = PCO$pos.neg)
  
  # Output the result
  out <- list(data = PCO$vectors, X = input$X, X.mean = X.mean, Xc = Xc, hat = hat.matrix, beta = pco.beta, fitted = pco.fitted, residuals = pco.resid, validation = diagn)
  return(out)
  rm(X.mean, Xc, hat.matrix, pco.beta, pco.fitted, pco.resid, J, SST, SSE, SSR, n, k, df.total, df.resid, diagn, out)
}

## Create a function for computing the diagnostics
pcoDiagnostics <- function(Y, D, iterations){
  #Y <- Y.total
  #D <- D.total
  
  N <- dim(D)[1]
  
  ## Run PCO
  PCO.result <- pco(D, varplot = FALSE)
  PCO <- PCO.result$vectors
  L <- PCO.result$values
  pos.neg <- L >= 0.0
  
  ## Scree plot with expectations for eigenvalues under the broken stick model
  message("Running the broken stick model method...")
  broken.stick <- sapply(1:N, function(k){ifelse(L[k] >= 0, sum(1/(k:N)), 0)})
  
  ## Bootstrap eigenvalue method
  message("Broken stick model complete. Running the bootstrap eigenvalue method...")
  nboot <- iterations
  
  L.boot <- do.call(rbind, pbapply::pblapply(1:nboot, function(iboot){
    index <- sample(1:N, replace = TRUE)
    D.boot <- D[index, index]
    pco(D.boot, varplot = FALSE)$values
  }))
  L.perc.boot <- t(apply(L.boot, 1, function(x){100 * x/sum(x)}))
  
  # Get the empirical 95% confidence interval on the bootstrap eigenvalues
  low <- function(x) quantile(x,probs=0.025)
  high <- function(x) quantile(x,probs=0.975)
  
  ## Permutation method (e.g., McCune et al. 2002; Clarke et al. 2008) with
  ## variation explained by each axis as either:
  ##    1. a fraction of the total (holistically), or
  ##    2. a fraction of the variation remaining given prior axes (conditionally)
  message("Bootstrap eigenvalue method complete. Running the permutation method...")
  nperm <- iterations
  nvars <- dim(Y)[2]   # no. of variables in original Y matrix
  
  L.perm <- do.call(rbind, pbapply::pblapply(1:nperm, function(iperm){
    # Reorganize the original data matrix
    Y.perm <- sapply(1:nvars, function(i){index <- sample(1:N, replace = FALSE); Y[index,i]})
    # Calculate the dissimilarity
    D.perm <- distance(Y = Y.perm, measure = "BC", trans = "none", adj = 0.03)
    # Run PCO
    lambda.perm <- pco(D.perm, varplot = F)$values
    return(lambda.perm)
  }))
  
  L.perc.perm1 <- t(apply(L.perm, 1, function(x){100 * x/sum(x)}))
  L.perc.perm2 <- t(apply(L.perm, 1, function(x){sapply(1:N, function(k){100 * x[k]/sum(x[k:N])})}))
  
  # Get the empirical 95% confidence interval on the permutation eigenvalues, in each case
  message("All models complete. Compiling results...")
  df <- data.frame(PC = 1:N, 
                   PC.name = paste("PCO", 1:N, sep = ""), 
                   Pos.Eigen = pos.neg, 
                   real = 100*L/sum(L), 
                   real.cum = cumsum(100*L/sum(L)), 
                   conditional = sapply(1:N, function(k){100*L[k]/sum(L[k:N])}), 
                   broken.stick = 100 * broken.stick/sum(broken.stick), 
                   lower.boot = apply(L.perc.boot, MARGIN = 2, low), 
                   upper.boot = apply(L.perc.boot, MARGIN = 2, high), 
                   lower.p1 = apply(L.perc.perm1, MARGIN = 2, low), 
                   upper.p1 = apply(L.perc.perm1, MARGIN = 2, high),
                   middle.p1 = apply(L.perc.perm1, MARGIN = 2, median),
                   lower.p2 = apply(L.perc.perm2, MARGIN = 2, low),
                   upper.p2 = apply(L.perc.perm2, MARGIN = 2, high),
                   middle.p2 = apply(L.perc.perm2, MARGIN =2, median))
  return(df)
  rm(N, PCO.result, Q, L, pos.neg, PCO, broken.stick, nboot, L.boot, L.perc.boot, low, high, nperm, nvars, L.perm, L.perc.perm1, L.perc.perm2, df)
  #rm(Y, D)
}
pcoDiagPlots <- function(pco.diag, max.axes = NULL, size.point, size.text, size.lab){
  #pco.diag <- pcodiagn.total
  
  library(ggplot2)
  
  ## Modify the theme for the plot
  theme.plot <- theme(text = element_text(family = "serif")) + 
    theme(plot.title = element_text(hjust = 0.5, size = size.lab*1.25, margin = margin(b = 0.05, unit = "inch")), 
          #plot.background = element_blank(), 
          plot.margin = unit(c(.1,.1,.1,.1), "inch")) +
    theme(axis.ticks = element_line(color = NA, linewidth = 1, linetype = "solid"), 
          axis.line = element_line(color = NA, linewidth = .1, linetype = "solid"), 
          axis.title=element_text(size=size.lab, margin = margin(t = 0.25, unit="inch")),  
          axis.title.x = element_text(vjust = 0), 
          axis.title.y = element_text(angle = 90, vjust = 1.5), 
          axis.text = element_text(size = size.lab*0.75, color = "black"), 
          axis.text.x = element_text(angle = 0, hjust = 0.5), 
          axis.text.y = element_text(angle = 0, hjust = 0)) + 
    theme(panel.border = element_rect(fill = NA, color = "black"), 
          panel.background = element_rect(fill = NA, color = NA), 
          panel.grid.major = element_line(color = "grey50"), 
          panel.grid.minor = element_line(color = NA), 
          panel.spacing.x = unit(0.25, units = "inch")) + 
    theme(legend.margin=margin(c(0.15,0.15,0.15,0.15), unit = "inch"), 
          legend.background = element_rect(fill = NA, color = NA), 
          legend.text=element_text(size = size.lab*0.75), 
          legend.title=element_text(size=size.lab*0.75), 
          #legend.position = "top", 
          legend.key = element_rect(color = NA, fill = NA), 
          legend.key.height = unit(0.25,"inch"), 
          legend.key.width = unit(0.25, "inch")) + 
    theme(strip.background = element_rect(fill = "gray90", color = "black"), 
          strip.text = element_text(size = size.lab*0.75), 
          strip.text.x = element_text(margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")), 
          strip.text.y = element_text(angle = -90, margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")))
  
  # Limit the number of PCO axes to show
  if(is.null(max.axes)){
    diag.plot <- pco.diag
  }else{
    diag.plot <- pco.diag[1:max.axes,]
  }
  
  
  # Broken stick method
  plot.bs <- ggplot(data = diag.plot, mapping = aes(x = PC, y = real)) + 
    geom_line() + 
    geom_point(size = size.point) + 
    geom_point(mapping = aes(y = broken.stick), size = size.point, shape = 1) + 
    coord_cartesian(ylim = c(-10,100)) + 
    xlab("PCO Axis Number") + ylab("Percent variation explained") + 
    ggtitle("Broken Stick Model") + 
    theme.plot
  
  # Bootstrapping method
  plot.boot <- ggplot(data = diag.plot, mapping = aes(x = PC, y = real)) + 
    geom_line() + 
    geom_point(size = size.point) + 
    geom_errorbar(mapping = aes(ymin = lower.boot, ymax = upper.boot)) + 
    coord_cartesian(ylim = c(-10,100)) + 
    xlab("PCO Axis Number") + ylab("Percent variation explained") + 
    ggtitle("Bootstrapped Eigenvalues Model") + 
    theme.plot
  
  # Holistic permutation model
  plot.perm1 <- ggplot(data = diag.plot, mapping = aes(x = PC, y = real)) + 
    geom_line() + 
    geom_point(size = size.point) + 
    geom_point(mapping = aes(y = middle.p1), size = size.point, shape = 1) + 
    geom_errorbar(mapping = aes(ymin = lower.p1, ymax = upper.p1)) + 
    coord_cartesian(ylim = c(-10,100)) + 
    xlab("PCO Axis Number") + ylab("Percent variation explained") + 
    ggtitle("Permutation Model") + 
    theme.plot
  
  # Conditional permutation model
  plot.perm2 <- ggplot(data = diag.plot, mapping = aes(x = PC, y = conditional)) + 
    geom_line() + 
    geom_point(size = size.point) + 
    geom_point(mapping = aes(y = middle.p2), size = size.point, shape = 1) + 
    geom_errorbar(mapping = aes(ymin = lower.p2, ymax = upper.p2)) + 
    coord_cartesian(ylim = c(-10,100)) + 
    xlab("PCO Axis Number") + ylab("Percent variation explained") + 
    ggtitle("Conditional Permutation Model") + 
    theme.plot
  
  plot.out <- ggpubr::ggarrange(plot.bs, plot.boot, plot.perm1, plot.perm2)
  plot.out
  
  return(plot.out)
  rm(theme.plot, plot.bs, plot.boot, plot.perm1, plot.perm2, plot.out)
  #rm(pco.diag)
}
pcoPlots <- function(ds, PCO, PCO.fit = NULL, X, Y, add.fitted, size.text, size.lab){
  #ds <- PCO.diagn.total
  #PCO <- PCO.df.total
  #X <- "PCO1"
  #Y <- "PCO2"
  
  library(ggplot2)
  
  ## Modify the theme for the plot
  theme.plot <- theme(text = element_text(family = "serif")) + 
    theme(plot.title = element_text(hjust = 0.5, size = size.lab*1.50, margin = margin(b = 0.5, unit = "inch")), 
          plot.background = element_blank(), 
          plot.margin = unit(c(.1,.1,.1,.1), "inch")) +
    theme(axis.ticks = element_line(color = NA, linewidth = 1, linetype = "solid"), 
          axis.line = element_line(color = NA, linewidth = .1, linetype = "solid"), 
          axis.title=element_text(size=size.lab*1.15, margin = margin(t = 0.25, unit="inch")),  
          axis.title.x = element_text(vjust = -1.0, hjust = 0.50), 
          axis.title.y = element_text(angle = 90, hjust = 0.50, vjust = 2.0), 
          axis.text = element_text(size = size.lab, color = "black"), 
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.0), 
          axis.text.y = element_text(angle = 0, hjust = 0)) + 
    theme(panel.border = element_rect(fill = NA, color = "black"), 
          panel.background = element_rect(fill = NA, color = NA), 
          panel.grid.major = element_line(color = NA), 
          panel.grid.major.x = element_line(color = NA), 
          panel.grid.major.y = element_line(color = NA), 
          panel.grid.minor = element_line(color = NA), 
          panel.grid.minor.x = element_line(color = NA), 
          panel.grid.minor.y = element_line(color = NA)) +  
    theme(legend.margin=margin(c(0.15,0.15,0.15,0.15), unit = "inch"), 
          legend.background = element_rect(fill = NA, color = NA), 
          legend.text=element_text(size = size.lab), 
          legend.title=element_text(size = size.lab), 
          legend.position = "right", 
          #legend.key = element_rect(color = "black", fill = NA), 
          legend.key.height = unit(0.25,"inch"), 
          legend.key.width = unit(0.25, "inch")) + 
    theme(strip.background = element_rect(fill = "gray85", color = "black"), 
          strip.text = element_text(size = size.lab), 
          strip.text.x = element_text(margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")), 
          strip.text.y = element_text(angle = -90, margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")))
  
  PCO.fit <- as.data.frame(PCO.fit)
  
  ## Create the plot
  if(isFALSE(add.fitted)){
    plot.out <- ggplot(data = PCO, mapping = aes(x = !!sym(X), y = !!sym(Y))) + 
      geom_hline(yintercept = 0, linewidth = 1.5, color = "grey80") + geom_vline(xintercept = 0, linewidth = 1.5, color = "grey80") + 
      geom_point(size = 2) + 
      xlab(paste(X, " (", round(ds[ds$PC.name == X, "real"], digits = 2), "% of total variation)", sep = "")) + 
      ylab(paste(Y, " (", round(ds[ds$PC.name == Y, "real"], digits = 2), "% of total variation)", sep = "")) + 
      theme.plot
  }else if(isTRUE(add.fitted)){
    plot.out <- ggplot(data = PCO, mapping = aes(x = !!sym(X), y = !!sym(Y))) + 
      geom_hline(yintercept = 0, linewidth = 1.5, color = "grey80") + geom_vline(xintercept = 0, linewidth = 1.5, color = "grey80") + 
      geom_point(size = 2) + 
      geom_point(data = PCO.fit, shape = 1, size = 2) + 
      xlab(paste(X, " (", round(ds[ds$PC.name == X, "real"], digits = 2), "% of total variation)", sep = "")) + 
      ylab(paste(Y, " (", round(ds[ds$PC.name == Y, "real"], digits = 2), "% of total variation)", sep = "")) + 
      theme.plot
  }
  
  plot.out
  
  return(plot.out)
  rm(theme.plot, plot.out)
  #rm(ds, PCO, X, Y, size.text, size.lab)
}

## Create a function for predicting a dropped site or time period
predictReg <- function(reg, axes, describe = TRUE){
  #reg <- reg.total
  #axes <- 1:3   ## Which axes should be used for calculation
  
  reg.drop <- reg[grepl("drop", names(reg))]
  reg.full <- do.call(c, reg[grepl("full", names(reg))])
  
  preds <- lapply(reg.drop, function(x){
    #x <- reg.drop$drop_FM106_WCS1
    
    # Combine the full data with the dropped data
    x1 <- c(reg.full, x)
    
    ## Extract only the dropped X and Y data
    dropped.y <- matrix(x1$full.data[!(row.names(x1$full.data) %in% row.names(x1$data)),], 
                        nrow = nrow(x1$full.data) - nrow(x1$data), 
                        dimnames = list(row.names(x1$full.data)[!(row.names(x1$full.data) %in% row.names(x1$data))], colnames(x1$full.data)))
    dropped.x <- matrix(x1$full.X[!(row.names(x1$full.X) %in% row.names(x1$X)),], 
                        nrow = nrow(x1$full.X) - nrow(x1$X), 
                        dimnames = list(row.names(x1$full.X)[!(row.names(x1$full.X) %in% row.names(x1$X))], colnames(x1$full.X)))
    
    full.fit.pred <- matrix(x1$full.fitted[!(row.names(x1$full.fitted) %in% row.names(x1$data)),], 
                            nrow = nrow(x1$full.fitted) - nrow(x1$data), 
                            dimnames = list(row.names(x1$full.fitted)[!(row.names(x1$full.fitted) %in% row.names(x1$data))], colnames(x1$full.fitted)))
    
    ## Center the dropped X values around the mean used in the regression
    Xc.pred <- dropped.x - rep(1, nrow(dropped.x)) %o% x1$X.mean
    
    ## Calculated fitted values based on the centered data
    Yhat.pred <- Xc.pred %*% x1$beta
    
    x2 <- c(x, list(data.pred = dropped.y, orig.fit = full.fit.pred, X.pred = dropped.x, Xc.pred = Xc.pred, Yhat.pred = Yhat.pred))
    
    return(x2)
    rm(x1, dropped.y, dropped.x, Xc.pred, Yhat.pred, x2)
    #rm(x)
  })
  
  
  # Calculate the mantel correlations used for model validation
  mantels1 <- lapply(preds, function(x){
    #x <- preds$drop_FM106_WCS1
    #x <- preds$`drop_-20`
    
    ## Define the objects used for Mantel correlation calculations
    full.data <- reg.full$full.data
    full.fit <- reg.full$full.fitted
    drop.data <- x$data
    drop.fit <- x$fitted
    full.fit.pred <- x$orig.fit
    fitted.pred <- rbind(drop.fit, x$Yhat.pred)[order(row.names(full.data)),]
    dropped.y <- x$data.pred
    yhat.pred <- x$Yhat.pred
    
    ## Correlation between the full original data and the fitted data from the full model
    m1 <- cor(dist(as.matrix(full.data[,axes])), dist(as.matrix(full.fit[,axes])))
    ## Correlation between the full drop.one model and the fitted data from that model
    m2 <- cor(dist(as.matrix(drop.data[,axes])), dist(as.matrix(drop.fit[,axes])))
    ## Correlation between the original data from the full model and the predicted value from the drop.one model
    m3 <- cor(dist(as.matrix(dropped.y[,axes])), dist(as.matrix(yhat.pred[,axes])))
    ## Correlation between the fitted data from the full model and predicted value from the drop.one model
    m4 <- cor(dist(as.matrix(full.fit.pred[,axes])), dist(as.matrix(yhat.pred[,axes])))
    ## Correlation between the full original data and fitted data from the drop.one model
    m5 <- cor(dist(as.matrix(full.data[,axes])), dist(as.matrix(fitted.pred[,axes])))
    
    
    out.m <- data.frame(m1, m2, m3, m4, m5)
    return(out.m)
    rm(full.data, full.fit, full.fit.pred, drop.data, drop.fit, fitted.pred, dropped.y, yhat.pred, m1, m2, m3, m4, m5)
  })
  
  if(isTRUE(describe)){
    message("Explanation of the Mantel Correlations: 
        m1 is the correlation between the full original data and the fitted data from the full model. This does not change across models.
        m2 is the correlation between the full drop.one model and the fitted data from that model.
        m3 is the correlation between the original data from the full model and predicted value from the drop.one model.
        m4 is the correlation between the fitted data from the full model and predicted value from the drop.one model.
        m5 is the correlation between the full original data and the fitted data from the drop.one model.")
  }
  
  mantels2 <- do.call(rbind, mantels1)
  
  preds2 <- c(full = list(reg.full), preds)
  
  out <- list(pred = preds2, mantel = mantels2)
  return(out)
  
  rm(reg.drop, reg.full, preds, mantels1, mantels2, out)
  #rm(reg, axes, describe)
}
predictPlots <- function(ds, plot.type, X, Y, fill, facet.plot = FALSE, average = FALSE, add.cor, cor.file, size.text, size.lab, size.point){
  #ds <- plot.site
  #plot.type <- "site"
  #X <- "PCO1"
  #Y <- "PCO2"
  #fill <- "drop"
  #facet.plot <- FALSE
  #average <- FALSE
  #add.cor <- TRUE
  #cor.file <- plot.spec
  #size.text <- 4
  #size.lab <- 9
  #size.point <- 3
  
  # Load the ggplot2 package
  require(ggplot2)
  
  # The theme to make the graph look how we want it to
  theme.plot <- theme(text = element_text(family = "serif")) + 
    theme(plot.title = element_text(hjust = 0.5, size = size.lab*1.50, margin = margin(b = 0.5, unit = "inch")), 
          #plot.background = element_blank(), 
          plot.background = element_rect(fill = "white"), 
          plot.margin = unit(c(.1,.1,.1,.1), "inch")) +
    theme(axis.ticks = element_line(color = NA, linewidth = 1, linetype = "solid"), 
          axis.line = element_line(color = NA, linewidth = .1, linetype = "solid"), 
          axis.title=element_text(size=size.lab*1.15, margin = margin(t = 0.25, unit="inch")),  
          axis.title.x = element_text(vjust = -1.0, hjust = 0.50), 
          axis.title.y = element_text(angle = 90, hjust = 0.50, vjust = 2.0), 
          axis.text = element_text(size = size.lab, color = "black"), 
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.0), 
          axis.text.y = element_text(angle = 0, hjust = 0)) + 
    theme(panel.border = element_rect(fill = NA, color = "black"), 
          panel.background = element_rect(fill = NA, color = NA), 
          panel.grid.major = element_line(color = NA), 
          panel.grid.major.x = element_line(color = NA), 
          panel.grid.major.y = element_line(color = NA), 
          panel.grid.minor = element_line(color = NA), 
          panel.grid.minor.x = element_line(color = NA), 
          panel.grid.minor.y = element_line(color = NA)) +  
    theme(legend.margin=margin(c(0.15,0.15,0.15,0.15), unit = "inch"), 
          legend.background = element_rect(fill = NA, color = NA), 
          legend.text=element_text(size = size.lab), 
          legend.title=element_text(size = size.lab), 
          legend.position = "right", 
          #legend.key = element_rect(color = "black", fill = NA), 
          legend.key.height = unit(0.25,"inch"), 
          legend.key.width = unit(0.25, "inch")) + 
    theme(strip.background = element_rect(fill = "gray85", color = "black"), 
          strip.text = element_text(size = size.lab), 
          strip.text.x = element_text(margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")), 
          strip.text.y = element_text(angle = -90, margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")))
  
  # Define the scheme for colors and shapes
  if(plot.type == "site"){
    scale_color <- scale_color_manual("Site", values = c("#332288", "#117733", "#44AA99", "#DDCC77", "#CC6677", "#AA4499", "#882255", "#EE9D28", 
                                                         "#332288", "#117733", "#44AA99", "#DDCC77", "#CC6677", 
                                                         "#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77"))
    scale_shape <- scale_shape_manual("Site", values = c(rep(15, 8), rep(16, 5), rep(17, 5)))
  }else if(plot.type == "time"){
    scale_color <- scale_color_manual("Month", values = c("red", "yellow", "blue", "purple", 
                                                          "red", "yellow", "blue", "purple", 
                                                          "red", "yellow", "blue", "purple"))
    scale_shape <- scale_shape_manual("Month", values = rep(c(15:18), each = 3))
  }else{
    stop("You must choose an appropriate plot.type. Choose one of c('site','time').")
  }
  
  if(isTRUE(average)){
    sites2 <- dplyr::summarise(dplyr::group_by(ds, response, ds, drop), X = mean(!!sym(X)), Y = mean(!!sym(Y)))
  }else{
    sites2 <- ds
    colnames(sites2)[colnames(sites2) %in% c(X, Y)] <- c("X", "Y")
  }
  
  # Should the plot be a facet plot or multi-plot
  if(isTRUE(facet.plot)){
    ds2 <- list(sites2)
    facets <- facet_grid(rows = vars(ds), cols = vars(response))
  }else{
    ds2 <- split(sites2, f = ~ ds + response)
    facets <- NULL
  }
  
  if(isTRUE(add.cor)){
    if(isTRUE(facet.plot)){
      spec1 <- list(plot.spec)
    }else{
      spec1 <- rep(split(plot.spec, f = ~ response), each = 3)
    }
  }
  
  plot.out <- lapply(1:length(ds2), function(i){
    ds.plot <- ds2[[i]]
    
    if(X == "PCO1" & Y == "PCO2"){
      if(unique(ds.plot$response) == "Total detections"){
        coords.x <- c(-0.8, 0.8)
        coords.y <- c(-0.5, 0.6)
      }else if(unique(ds.plot$response) == "Successful crossings"){
        coords.x <- c(-0.7, 0.8)
        coords.y <- c(-0.9, 0.6)
      }else if(unique(ds.plot$response) == "Unsuccessful crossings"){
        coords.x <- c(-0.7, 0.6)
        coords.y <- c(-0.6, 0.7)
      }
    }else if(X == "PCO1" & Y == "PCO3"){
      if(unique(ds.plot$response) == "Total detections"){
        coords.x <- c(-0.8, 0.8)
        coords.y <- c(-0.6, 0.7)
      }else if(unique(ds.plot$response) == "Successful crossings"){
        coords.x <- c(-0.6, 1.6)
        coords.y <- c(-0.7, 0.7)
      }else if(unique(ds.plot$response) == "Unsuccessful crossings"){
        coords.x <- NULL
        coords.y <- NULL
      }
    }else if(X == "PCO2" & Y == "PCO3"){
      if(unique(ds.plot$response) == "Total detections"){
        coords.x <- c(-0.6, 0.6)
        coords.y <- c(-0.6, 0.7)
      }else if(unique(ds.plot$response) == "Successful crossings"){
        coords.x <- c(-0.9, 0.5)
        coords.y <- c(-0.6, 0.6)
      }else if(unique(ds.plot$response) == "Unsuccessful crossings"){
        coords.x <- NULL
        coords.y <- NULL
      }
    }else{
      coords.x <- NULL
      coords.y <- NULL
    }
    
    if(isTRUE(add.cor)){
      spec.plot <- spec1[[i]]
      mult <- 1.25
      spec.line <- geom_segment(data = spec.plot, mapping = aes(x = 0, y = 0, xend = !!sym(X)*mult, yend = !!sym(Y)*mult), color = "grey80", arrow = arrow(length = unit(0.15, "inches"), type = "closed"))
      spec.label <- geom_text(data = spec.plot, mapping = aes(x = !!sym(X)*mult, y = !!sym(Y)*mult, label = species), color = "grey30", size = size.text, family = "serif")
    }else{
      spec.line <- NULL
      spec.label <- NULL
    }
    
    ggplot(data = ds2[[i]], mapping = aes(x = X, y = Y)) + 
      geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
      spec.line + 
      geom_point(mapping = aes(color = !!sym(fill), shape = !!sym(fill)), size = size.point) + 
      spec.label + 
      scale_x_continuous(X, breaks = seq(-2.0, 2.0, 0.2), limits = c(-5.0, 5.0)) + 
      scale_y_continuous(Y, breaks = seq(-5.0, 5.0, 0.2), limits = c(-5.0, 5.0)) + 
      scale_color + scale_shape + 
      coord_fixed(xlim = coords.x, ylim = coords.y) + 
      #coord_fixed() + 
      facets + theme.plot
    #rm(ds.plot, spec.line, spec.label, spec.plot)
  })
  #plot.out[[1]]
  #plot.out[[2]]
  #plot.out[[3]]
  #plot.out[[4]]
  #plot.out[[5]]
  #plot.out[[6]]
  #plot.out[[7]]
  #plot.out[[8]]
  #plot.out[[9]]
  
  if(isFALSE(facet.plot)){
    names(plot.out) <- names(ds2)
  }
  
  if(length(plot.out) == 1){
    plot.out <- plot.out[[1]]
  }else{
    message("The output will be a list because multiple plots are created")
  }
  
  return(plot.out)
  rm(theme.plot, scale_color, scale_shape, ds2, facets, plot.out, spec1)
  #rm(ds, plot.type, X, Y, fill, facet.plot, size.text, size.lab, size.point)
}


#######################################################################################################################
###################################### Part 1: Prep Raw Data for Primer Analyses ###################################### 
#######################################################################################################################

rm(list = ls()); gc()

# Functions used in this section: ####
  ## processInts()
  ## distance()
  ## dbmemFun()  


# Step 0: What are all the files? ####
files_data <- fs::dir_ls(path = file.path(getwd(), "Data"), type = )
files_data


# Step 1: Load the data ####
## Interactions data
files_int <- files_data[grep("Interactions", files_data)]
ints1 <- do.call(rbind, lapply(openxlsx::getSheetNames(files_int)[1:3], function(x){openxlsx::read.xlsx(files_int, sheet = x)}))
ints1[1:5,]

## DataOrganize data
files_DO <- files_data[grep("DataOrganize", files_data)]
DO1 <- pbapply::pblapply(files_DO, read.table, header = TRUE)
names(DO1) <- sub("_all_20240405.txt", "", sub("DataOrganize_", "", basename(names(DO1))))
DO1[[1]][1:5,]

## CT table
files_CT <- files_data[grep("CTtable", files_data)]
CT1 <- lapply(openxlsx::getSheetNames(files_CT), function(x){openxlsx::read.xlsx(files_CT, sheet = x, detectDates = TRUE)})
names(CT1) <- openxlsx::getSheetNames(files_CT)
CT1[[1]][1:5,]

## Environmental data
### Camera Level
files_env_cam <- files_data[grep("EnvData_cam", files_data)]
envdata_cam <- lapply(openxlsx::getSheetNames(files_env_cam), function(x){openxlsx::read.xlsx(files_env_cam, sheet = x)})
names(envdata_cam) <- openxlsx::getSheetNames(files_env_cam)
### Site Level
files_env_site <- files_data[grep("EnvData_site", files_data)]
envdata_site <- lapply(openxlsx::getSheetNames(files_env_site)[1:3], function(x){openxlsx::read.xlsx(files_env_site, sheet = x)})
names(envdata_site) <- openxlsx::getSheetNames(files_env_site)[1:3]

## Species list
files_spec <- files_data[grep("Species_List", files_data)]
spec_list <- openxlsx::read.xlsx(files_spec, sheet = 1)


# Step 2: Process the interactions data into monthly interactions by species ####
## Select only those species that will be analysed (wild mammals)
spec.analyse <- spec_list$Timelapse[spec_list$Analysis == 1]

## Process the interactions data so it is usable
ints2 <- processInts(x = ints1, species = spec.analyse)


# Step 3: Process environmental data so it can be combined with interactions ####
## Create the full envdata file
env.site.all <- merge.data.frame(merge.data.frame(envdata_site$Spatiotemporal, envdata_site$Spatial, by = "WCS_ID", all.x = TRUE), envdata_site$Temporal, by = "TimeID", all.x = TRUE)
env.site.all[1:5,]

## For combining with the Species data
env.shell <- env.site.all[order(env.site.all$WCS_Time_ID), c("WCS_Time_ID", "WCS_ID", "TimeID")]
env.shell[1:5,]

## Define which fields are the factors
factors <- env.site.all[order(env.site.all$WCS_Time_ID), 
                         c("WCS_Time_ID", "WCS_ID", "TimeID", 
                           "Name", "Road", "Year", "Month", 
                           "ConPeriod", "ConMonths", "CatWalks", "Substrate", "Fencing", 
                           "IntsAnalysis")]
factors[1:5,]

## Define which fields are covariates
cov.all.list <- list(structural = env.site.all[order(env.site.all$WCS_Time_ID), 
                                            c("WCS_Time_ID", 
                                              "Openness_m", "CatWalks.No", "CatWalks.Yes", "Substrate.Water", "Substrate.Concrete", "Substrate.Dirt", "Fencing_m")], 
                 environmental = env.site.all[order(env.site.all$WCS_Time_ID), 
                                               c("WCS_Time_ID", 
                                                 "PropNatural", "PropWater", "PropWoody", "Precip_total_mm")], 
                 anthropogenic = env.site.all[order(env.site.all$WCS_Time_ID), 
                                               c("WCS_Time_ID", 
                                                 "SpeedLimit", "VehicleTraffic", "HumanActivity", "DomesticActivity", "LivestockActivity", "PropBuilding")], 
                 spatial = env.site.all[order(env.site.all$WCS_Time_ID),
                                         c("WCS_Time_ID", 
                                           "UTM_14N_X", "UTM_14N_Y")], 
                 temporal = env.site.all[order(env.site.all$WCS_Time_ID),
                                          c("WCS_Time_ID", 
                                            "ConMonths_cov")]
)

## Combine into single file of covariates
### Include spatial and temporal for completeness
cov.all <- do.call(cbind, lapply(1:length(cov.all.list), function(i){
  if(i > 1){
    c1 <- data.frame(cov.all.list[[i]][,-1])
    colnames(c1) <- colnames(cov.all.list[[i]])[-1]
  }else{
    c1 <- data.frame(cov.all.list[[i]])
    colnames(c1) <- colnames(cov.all.list[[i]])
  }
  return(c1)
}))
## Define indicators for use in Primer
### Including spatial and temporal
cov.ind.all <- with(cov.all.list, data.frame(cov = c(colnames(spatial)[-1], colnames(temporal)[-1], colnames(structural)[-1], colnames(environmental)[-1], colnames(anthropogenic)[-1]), 
                           group = c(rep("spatial", ncol(spatial)-1), 
                                     rep("temporal", ncol(temporal)-1), 
                                     rep("structural", ncol(structural)-1), 
                                     rep("environmental", ncol(environmental)-1), 
                                     rep("anthropogenic", ncol(anthropogenic)-1)), 
                           single = c("X", "Y", 
                                      "ConMonths", 
                                      "openness", "catwalk", "catwalk", "substrate", "substrate", "substrate", "fencing", 
                                      "Natural", "Water", "Woody", "Precip", 
                                      "SpeedLimit", "VehicleTraffic", "HumanAct", "DomesticAct", "LivestockAct", "Buildings")))

## Examine variance inflation factor among sets of covariates
with(cov.all.list, diag(solve(cor(as.matrix(structural[,c(2,4,6,7,8)])))))
with(cov.all.list, diag(solve(cor(as.matrix(environmental[,c(2:5)])))))
with(cov.all.list, diag(solve(cor(as.matrix(anthropogenic[,c(2:7)])))))


# Step 4: Process spatial and temporal variables to properly model spatial and temporal autocorrelation ####
## Spatial Autocorrelation
### Total Detections
spat.out.tot <- dbmemFun(x = cov.all.list$spatial, type = "spatial")
spat.pcnm.tot <- with(spat.out.tot, do.call(rbind, lapply(2:ncol(df), function(i){data.frame(cov.all.list$spatial, axis = factor(colnames(df)[i]), value = df[,i])})))
### Interactions
spat.out.int <- dbmemFun(x = cov.all.list$spatial[factors$IntsAnalysis == 1,], type = "spatial")
spat.pcnm.int <- with(spat.out.int, do.call(rbind, lapply(2:ncol(df), function(i){data.frame(cov.all.list$spatial[factors$IntsAnalysis == 1,], axis = factor(colnames(df)[i]), value = df[,i])})))

## Temporal Autocorrelation
### Total Detections
temp.out.tot <- dbmemFun(x = cov.all.list$temporal, type = "temporal")
temp.pcnm.tot <- with(temp.out.tot, do.call(rbind, lapply(2:ncol(df), function(i){data.frame(cov.all.list$temporal, axis = factor(colnames(df)[i]), value = df[,i])})))
### Interactions
temp.out.int <- dbmemFun(x = cov.all.list$temporal[factors$IntsAnalysis == 1,], type = "temporal")
temp.pcnm.int <- with(temp.out.int, do.call(rbind, lapply(2:ncol(df), function(i){data.frame(cov.all.list$temporal[factors$IntsAnalysis == 1,], axis = factor(colnames(df)[i]), value = df[,i])})))

## Create plots plots of the spatial and temporal autocorrelation
### Spatial Autocorrelation
#### Total detections
#### Interactions analyses

### Create the plots
library(ggplot2)
plot.pcnm.spat.tot <- ggplot(data = spat.pcnm.tot, mapping = aes(x = UTM_14N_X, y = UTM_14N_Y, shape = factor(value/abs(value)), size = abs(value))) + 
  geom_point() + 
  scale_shape_manual("pos.neg", values = c(1,16), labels = c("negative", "positive")) + 
  scale_size_binned("Value", n.breaks = 5) + 
  xlab("X") + ylab("Y") + 
  facet_wrap(facets = vars(axis), nrow = 4, ncol = 3) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90))
#### Interactions detections
plot.pcnm.spat.int <- ggplot(data = spat.pcnm.int, mapping = aes(x = UTM_14N_X, y = UTM_14N_Y, shape = factor(value/abs(value)), size = abs(value))) + 
  geom_point() + 
  scale_shape_manual("pos.neg", values = c(1,16), labels = c("negative", "positive")) + 
  scale_size_binned("Value", n.breaks = 5) + 
  xlab("X") + ylab("Y") + 
  facet_wrap(facets = vars(axis), nrow = 4, ncol = 3) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90))
#### Combined plots
plot.pcnm.spat <- ggpubr::ggarrange(plot.pcnm.spat.tot, plot.pcnm.spat.int, ncol = 2, common.legend = TRUE, labels = "AUTO", font.label = list(family = "serif")) + 
  ggpubr::bgcolor("white") + ggpubr::border(color = NA)
plot.pcnm.spat
### Temporal Autocorrelation
#### Total detections
plot.pcnm.temp.tot <- ggplot(data = temp.pcnm.tot, mapping = aes(x = ConMonths_cov, y = value)) + 
  geom_line() + geom_point(mapping = aes(shape = factor(value/abs(value)), size = abs(value))) + 
  scale_shape_manual("pos.neg", values = c(1,16), labels = c("negative", "positive")) + 
  scale_size_binned("Value", n.breaks = 5) + 
  xlab("Month since construxction") + ylab("dbMEM eigenvalue") + 
  facet_wrap(facets = vars(axis), ncol = 3, nrow = 8) + 
  theme_bw()
#### Interactions detections
plot.pcnm.temp.int <- ggplot(data = temp.pcnm.int, mapping = aes(x = ConMonths_cov, y = value)) + 
  geom_line() + geom_point(mapping = aes(shape = factor(value/abs(value)), size = abs(value))) + 
  scale_shape_manual("pos.neg", values = c(1,16), labels = c("negative", "positive")) + 
  scale_size_binned("Value", n.breaks = 5) + 
  xlab("Month since construxction") + ylab("dbMEM eigenvalue") + 
  facet_wrap(facets = vars(axis), ncol = 3, nrow = 8) + 
  theme_bw()
#### Combined plots
plot.pcnm.temp <- ggpubr::ggarrange(plot.pcnm.temp.tot, plot.pcnm.temp.int, ncol = 2, common.legend = TRUE, labels = "AUTO", font.label = list(family = "serif")) + 
  ggpubr::bgcolor("white") + ggpubr::border(color = NA)
plot.pcnm.temp

## Save the plots
#w.pcnm <- 14
#h.pcnm <- 6
ggsave(filename = file.path("Results", paste("PCNM_spat_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = plot.pcnm.spat, device = "tiff", width = 10, height = 8, dpi = 600, compression = "lzw")
ggsave(filename = file.path("Results", paste("PCNM_temp_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = plot.pcnm.temp, device = "tiff", width = 10, height = 12, dpi = 600, compression = "lzw")


# Step 5: Combine interactions data with environmental data ####
## Merge Interactions data with full shell
ints3 <- lapply(ints2, function(x){
  x1 <- merge.data.frame(env.shell, x, by = "WCS_Time_ID", all.x = TRUE)
  x1[is.na(x1)] <- 0
  x2 <- x1[order(x1$WCS_Time_ID),]
  d1 <- x2[,colnames(x)]
  return(d1)
})

## Create a list object with all factor and covariate information
envs1 <- list(factors = factors, 
              spatial_all = spat.pcnm.tot, spatial_ints = spat.pcnm.int, 
              temporal_all = temp.pcnm.tot, temporal_ints = temp.pcnm.int, 
              covariates = cov.all, indicators = cov.ind.all)

## Save the output for use in Primer
#openxlsx::write.xlsx(ints3, file = file.path("Output_to_Others", paste("Ints_PRIMER_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = "")))
#openxlsx::write.xlsx(envs1, file = file.path("Output_to_Others", paste("Covariates_PRIMER_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = "")))


# Step 6: Diagnostic information for the camera data ####
## Calculate the number of active trap nights from DataOrganize data
DO2 <- lapply(1:length(DO1), function(i){
  x1 <- cameraTrapping::calculateEvents(do = DO1[[i]], envdata = envdata_cam[[i]], sort.col = "Site", start_date = "2010-01-01", interval = "30 min")
})
names(DO2) <- names(DO1)

camop <- list(stationCol = "Site", setupCol = "Setup_date", retrievalCol = "Retrieval_date", hasProblems = T, cameraCol = "Camera", byCamera = F, allCamsOn = F, camerasIndependent = F, sessionCol = "Session")
DO3 <- lapply(1:length(DO2), function(i){
  name <- names(DO2)[i]
  x1 <- cameraTrapping::summarizeEvents(DO2[[i]], ct = CT1[[i]], unit = "1 month", include = c("bobcat"), camOP = camop, out_form = "long", out_data = "AB", out_correction = "raw")
  
  if(name == "SH100"){
    WCS1 <- x1[x1$Site == "WCS1" & x1$Interval >= "2017-10-01" & x1$Interval <= "2019-05-01",]
    WCS2 <- x1[x1$Site == "WCS2" & x1$Interval >= "2017-10-01" & x1$Interval <= "2019-05-01",]
    WCS3 <- x1[x1$Site == "WCS3" & x1$Interval >= "2017-12-01" & x1$Interval <= "2019-05-01",]
    WCS3a <- x1[x1$Site == "WCS3A" & x1$Interval >= "2016-09-01" & x1$Interval <= "2019-05-01",]
    WCS4 <- x1[x1$Site == "WCS4" & x1$Interval >= "2017-12-01" & x1$Interval <= "2019-05-01",]
    x2 <- rbind(WCS1, WCS2, WCS3, WCS3a, WCS4)
  }else if(name == "FM106"){
    WCS1 <- x1[x1$Site == "WCS1" & x1$Interval >= "2019-07-01" & x1$Interval <= "2020-12-01",]
    WCS2 <- x1[x1$Site == "WCS2" & x1$Interval >= "2019-07-01" & x1$Interval <= "2020-12-01",]
    WCS3 <- x1[x1$Site == "WCS3" & x1$Interval >= "2019-07-01" & x1$Interval <= "2020-12-01",]
    WCS4 <- x1[x1$Site == "WCS4" & x1$Interval >= "2019-07-01" & x1$Interval <= "2020-12-01",]
    WCS5 <- x1[x1$Site == "WCS5" & x1$Interval >= "2019-07-01" & x1$Interval <= "2020-12-01",]
    WCS6 <- x1[x1$Site == "WCS6" & x1$Interval >= "2019-07-01" & x1$Interval <= "2020-12-01",]
    WCS7 <- x1[x1$Site == "WCS7" & x1$Interval >= "2019-07-01" & x1$Interval <= "2020-12-01",]
    WCS8 <- x1[x1$Site == "WCS8" & x1$Interval >= "2019-07-01" & x1$Interval <= "2020-12-01",]
    x2 <- rbind(WCS1, WCS2, WCS3, WCS4, WCS5, WCS6, WCS7, WCS8)
  }else{
    x2 <- x1
  }
  return(x2)
  rm(name, x1, WCS1, WCS2, WCS3, WCS3a, WCS4, WCS5, WCS6, WCS7, WCS8, x2)
})
names(DO3) <- names(DO2)

DO3$FM106[1:15,]
DO3$FM1847[1:15,]
DO3$SH100

DO4 <- do.call(rbind, lapply(1:length(DO3), function(i){
  name <- names(DO3)[i]
  x <- data.frame(road = name, DO3[[i]])
  x1 <- dplyr::summarise(dplyr::group_by(x, road, Site), total.nights = sum(totald))
  return(x1)
}))
DO4


## Calculate the number of cameras per site
CamsPerSite <- lapply(1:length(CT1), function(i){
  name <- names(CT1)[i]
  x <- CT1[[i]]
  
  if(name == "FM106"){
    C1 <- x[x$ConPeriod == "ConLate" & x$Retrieval_date >= "2019-12-01",]
    C2 <- x[x$ConPeriod == "PostCon" & x$Setup_date <= "2020-12-01",]
    x2 <- rbind(C1, C2)
  }else if(name == "SH100"){
    C1 <- x[x$Setup_date <= "2019-05-01",]
    x2 <- C1
  }else{
    x2 <- x
  }
  
  dplyr::summarise(dplyr::group_by(x2, ConPeriod, Site), n = dplyr::n())
})
CamsPerSite

#######################################################################################################################
#################################### Part 2: Create Plots of PCO and dbRDA Results #################################### 
#######################################################################################################################

#rm(list = ls()); gc()

# Step 1: Load the data from Primer ####
## View all files
files_primer <- fs::dir_ls(path = file.path("Data_from_Primer"), type = "file", recurse = TRUE)
files_plot <- files_primer[grep("dbRDA_detrended", files_primer)]


# Step 2: Prepare the data for plotting ####
## PCO Plots
### Prepping the site scores for plotting
#pco.sitescores <- ordPlotPrep(files = files_plot, plot.type = "PCO", score.type = "site", axes = 4)
### Prep species scores for plotting
#pco.specscores <- ordPlotPrep(files = files_plot, plot.type = "PCO", score.type = "spec", axes = 4)
### Combine site and species scores
#pco.full <- list(sitescores = pco.sitescores, specscores = pco.specscores)

## dbRDA Plots
### Prepping the site scores for plotting
dbrda.sitescores <- ordPlotPrep(files = files_plot, plot.type = "dbRDA", score.type = "site", axes = 4)
### Prepping the species scores for plotting
dbrda.specscores <- ordPlotPrep(files = files_plot, plot.type = "dbRDA", score.type = "spec", axes = 4)
dbrda.specscores$SpecShort <- factor(dbrda.specscores$Species, 
                                     levels = c("Armadillo", "Beaver", "Bobcat", "Cottontail", "Coyote", 
                                                "Feral hog", "Fox squirrel", "Gray squirrel", "Grey fox", "Jackrabbit", 
                                                "Javelina", "Mexican ground squirrel", "Nilgai", "Nutria", "Ocelot", 
                                                "Opossum", "Raccoon", "Rodent", "Striped skunk", "Unk mammal", 
                                                "Weasel", "White-tailed deer"), 
                                     labels = c("ARMA", "BEAV", "BOBC", "ECOT", "COYO", 
                                                "FHOG", "FSQR", "GSQR", "GFOX", "JACK", 
                                                "JAVE", "MSQR", "NILG", "NUTR", "OCEL", 
                                                "OPOS", "RACC", "UROD", "SSKU", "UMAM", 
                                                "WEAS", "DEER"))
### Prepping the covariate scores for plotting
dbrda.covscores <- ordPlotPrep(files = files_plot, plot.type = "dbRDA", score.type = "pred", axes = 4)
### Combine site, species, and covariate scores
dbrda.full <- list(sitescores = dbrda.sitescores, specscores = dbrda.specscores, covscores = dbrda.covscores)
#dbrda.covs <- list(sitescores = dbrda.sitescores, specscores = dbrda.specscores, covscores = dbrda.cov.covs)
#dbrda.spat <- list(sitescores = dbrda.sitescores, specscores = dbrda.specscores, covscores = dbrda.cov.spat)
#dbrda.temp <- list(sitescores = dbrda.sitescores, specscores = dbrda.specscores, covscores = dbrda.cov.temp)


# Step 3: Make the plots ####

## PCO Plots
#pco.plot.12 <- plotOrdination(ds = pco.full, plot.type = "PCO", response = Response, X = "PCO1", Y = "PCO2", spec.lab = "SpecNum", cov.lab = NULL, include.lines = FALSE, facet.plot = TRUE, size.text = sizeText, size.lab = sizeLab, size.point = sizePoint)
#pco.plot.13 <- plotOrdination(ds = pco.full, plot.type = "PCO", response = Response, X = "PCO1", Y = "PCO3", spec.lab = "SpecNum", cov.lab = NULL, include.lines = FALSE, facet.plot = TRUE, size.text = sizeText, size.lab = sizeLab, size.point = sizePoint)
#pco.plot.23 <- plotOrdination(ds = pco.full, plot.type = "PCO", response = Response, X = "PCO2", Y = "PCO3", spec.lab = "SpecNum", cov.lab = NULL, include.lines = FALSE, facet.plot = TRUE, size.text = sizeText, size.lab = sizeLab, size.point = sizePoint)
#pco.plot.12
#pco.plot.13
#pco.plot.23

## dbRDA Plots
### Shared arguments
ARGS <- list(ds = dbrda.full, plot.type = "dbRDA", response = "Road_WCS", 
             spec.lab = "SpecShort", cov.lab = "Covariate", include.lines = TRUE, average = FALSE, 
             size.text = 4, size.lab = 12, size.point = 3)
### Vary the axes
dbrda.plot.12 <- do.call(plotOrdination, c(X = "dbRDA1", Y = "dbRDA2", ARGS))
dbrda.plot.13 <- do.call(plotOrdination, c(X = "dbRDA1", Y = "dbRDA3", ARGS))
dbrda.plot.23 <- do.call(plotOrdination, c(X = "dbRDA2", Y = "dbRDA3", ARGS))
dbrda.plot.12$Total
dbrda.plot.13$Total
dbrda.plot.23$Total

## Save the plots
w.dbrda <- 12
h.dbrda <- 12.75
pbapply::pblapply(1:length(dbrda.plot.12), function(i){
  x <- dbrda.plot.12[[i]]
  name <- names(dbrda.plot.12)[i]
  ggsave(filename = file.path("Results", paste("dbRDA_", name, "_12_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
         plot = x, device = "tiff", width = w.dbrda, height = h.dbrda, dpi = 600, compression = "lzw")
})
pbapply::pblapply(1:length(dbrda.plot.13), function(i){
  x <- dbrda.plot.13[[i]]
  name <- names(dbrda.plot.13)[i]
  ggsave(filename = file.path("Results", paste("dbRDA_", name, "_13_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
         plot = x, device = "tiff", width = w.dbrda, height = h.dbrda, dpi = 600, compression = "lzw")
})
pbapply::pblapply(1:length(dbrda.plot.23), function(i){
  x <- dbrda.plot.23[[i]]
  name <- names(dbrda.plot.23)[i]
  ggsave(filename = file.path("Results", paste("dbRDA_", name, "_23_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
         plot = x, device = "tiff", width = w.dbrda, height = h.dbrda, dpi = 600, compression = "lzw")
})

#######################################################################################################################
############################################# Part 3: Variance Partitioning ########################################### 
#######################################################################################################################

#rm(list = ls()); gc()

# Step 1: Load the required data ####
## Load the variance partitioning results
vp1 <- openxlsx::read.xlsx(file.path("Data_from_Primer", "varpart_20240514.xlsx"), sheet = 2)
vp2 <- split(vp1, f = vp1$Response)
### Check the order of models
openxlsx::read.xlsx(file.path("Data_from_Primer", "varpart_20240514.xlsx"), sheet = 4)

# Step 2: Calculate the fraction of variation explained by each set of variables ####
## Calculate the variance explained
vp3 <- lapply(vp2, varpart5, Xs = c("str", "env", "ant", "spat", "temp"), type = "AdjustedR2")

## View the results
vp3$Total
vp3$Success
vp3$Failure

## Select only those responses needed for plotting
vp4 <- lapply(vp3, function(x){x[x$name %in% c("B1", "B2", "B3", "B6", "B7", "B10", "B16", "B32"),]})

## Combine variance partitioning results into a file that can be outputted
vp5 <- do.call(rbind, lapply(1:length(vp3), function(i){data.frame(dataset = names(vp3)[i], vp3[[i]])}))

## Save the variance partitioning results
#openxlsx::write.xlsx(x = vp5, file = file.path("Results", paste("varpart_result_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = "")))


# Step 3: Plot the amount of variance explained by each set of variables of interest ####
## Create the plots
vp.plot1 <- lapply(vp4, function(x){varpartPlot(data = x, size.text = 2.75, size.lab = 6)})
vp.plot1
vp.plot2 <- vp.plot1[c("Total", "Success", "Failure")]

vp.plot3 <- ggpubr::ggarrange(plotlist = vp.plot2, nrow = 2, ncol = 2, labels = "AUTO", font.label = list(family = "serif", size = 12)) + 
  ggpubr::bgcolor("white") + ggpubr::border(color = NA)
vp.plot3

## Save the plots
### Full plots
ggsave(filename = file.path("Results", paste("varpart_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), plot = vp.plot3, 
       device = "tiff", width = 7.0, height = 6.5, dpi = 600, compression = "lzw")
### Just total detections (for presentation purposes)
ggsave(filename = file.path("Results", paste("varpart_total_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = vp.plot1$Total, device = "tiff", width = 10, height = 10, dpi = 600, compression = "lzw")


#######################################################################################################################
############################################ Part 4: Regression on PCO axes ########################################### 
#######################################################################################################################

#rm(list = ls()); gc()

# Functions used in this section: ####
  ## distance()
  ## pco()
  ## pcoDiagnostics()
  ## pcoDiagPlots()
  ## pcoPlots()


# Step 1: Load the required data ####
## View all files
files_primer <- fs::dir_ls(path = file.path("Data_from_Primer"), type = "file", recurse = TRUE)
files_reg <- files_primer[grep("regression", files_primer)]

## Load in all the data and prep it for analyses
in.total.all <- prepData(files = files_reg, type = "total", mo = "all", site = "all", drop.pred = NULL, drop_mo = FALSE, drop_site = FALSE, autocorrelation = FALSE)
in.succ.all <- prepData(files = files_reg, type = "succ", mo = "all", site = "all", drop.pred = NULL, drop_mo = FALSE, drop_site = FALSE, autocorrelation = FALSE)
in.fail.all <- prepData(files = files_reg, type = "fail", mo = "all", site = "all", drop.pred = NULL, drop_mo = FALSE, drop_site = FALSE, autocorrelation = FALSE)

## Assess multicollinearity of the full model
### Visual assessment
with(in.total.all$full, pairs(X))
with(in.succ.all$full, pairs(X))
with(in.fail.all$full, pairs(X))
### Initially remove dummy variables in calculation of VIF
with(in.total.all$full, diag(solve(cor(X[,-c(2:4)]))))
with(in.succ.all$full, diag(solve(cor(X[,-c(2:4)]))))
with(in.fail.all$full, diag(solve(cor(X[,-c(2:4)]))))
### Remove vehicle traffic due to high VIF
with(in.total.all$full, diag(solve(cor(X[,-c(2:4,11)]))))
with(in.succ.all$full, diag(solve(cor(X[,-c(2:4,11)]))))
with(in.fail.all$full, diag(solve(cor(X[,-c(2:4,11)]))))
### Examine effect of dummy variables on VIF scores
#### Add cat walks
with(in.total.all$full, diag(solve(cor(X[,-c(3:4,11)]))))
with(in.succ.all$full, diag(solve(cor(X[,-c(3:4,11)]))))
with(in.fail.all$full, diag(solve(cor(X[,-c(3:4,11)]))))
#### Add substrate
with(in.total.all$full, diag(solve(cor(X[,-c(2,11)]))))
with(in.succ.all$full, diag(solve(cor(X[,-c(2,11)]))))
with(in.fail.all$full, diag(solve(cor(X[,-c(2,11)]))))
#### Add both dummy variables
with(in.total.all$full, diag(solve(cor(X[,-c(11)]))))
with(in.succ.all$full, diag(solve(cor(X[,-c(11)]))))
with(in.fail.all$full, diag(solve(cor(X[,-c(11)]))))
### Suggestions
#### Remove the variable "Fencing" from all models due to high multicollinearity

## Load in all the data and prep it for analyses
drop.preds1 <- c("log(Fencing_m+1)", "log(DomesticActivity+1)", "log(LivestockActivity+1)")
drop.preds2 <- c("log(Fencing_m+1)")
in.total <- prepData(files = files_reg, type = "total", mo = "all", site = "all", drop.pred = drop.preds1, drop_mo = TRUE, drop_site = TRUE, autocorrelation = FALSE)
in.succ <- prepData(files = files_reg, type = "succ", mo = "all", site = "all", drop.pred = drop.preds2, drop_mo = TRUE, drop_site = TRUE, autocorrelation = FALSE)
in.fail <- prepData(files = files_reg, type = "fail", mo = "all", site = "all", drop.pred = drop.preds1, drop_mo = TRUE, drop_site = TRUE, autocorrelation = FALSE)

## Load in the detrended data and prep it for analyses
in.total.cor <- prepData(files = files_reg, type = "total", mo = "all", site = "all", drop.pred = drop.preds1, drop_mo = FALSE, drop_site = FALSE, autocorrelation = TRUE)
in.succ.cor <- prepData(files = files_reg, type = "succ", mo = "all", site = "all", drop.pred = drop.preds2, drop_mo = FALSE, drop_site = FALSE, autocorrelation = TRUE)
in.fail.cor <- prepData(files = files_reg, type = "fail", mo = "all", site = "all", drop.pred = drop.preds1, drop_mo = FALSE, drop_site = FALSE, autocorrelation = TRUE)


# Step 2: Conduct PCO analysis ####
## Run the PCOs on the full in data
pco.total <- pbapply::pblapply(in.total, pcoOut, out.diag = FALSE)
pco.succ <- pbapply::pblapply(in.succ, pcoOut, out.diag = FALSE)
pco.fail <- pbapply::pblapply(in.fail, pcoOut, out.diag = FALSE)

## Run the PCOs on the detrended datasets
pco.total.cor <- pbapply::pblapply(in.total.cor, pcoOut, out.diag = FALSE)
pco.succ.cor <- pbapply::pblapply(in.succ.cor, pcoOut, out.diag = FALSE)
pco.fail.cor <- pbapply::pblapply(in.fail.cor, pcoOut, out.diag = FALSE)


# Step 3: Determine appropriate number of PCO axes ####
## Set the size of the labels and text in the plot
sizePoint <- 2
sizeText <- 4
sizeLab <- 12

## Run the diagnostics
iter <- 999
#PCO.diagn.total <- with(in.total$full, pcoDiagnostics(Y = Y, D = D, iterations = iter))
PCO.diagn.total <- readRDS(file.path("Results", "PCO_diagnostics_total_20240508.RDS"))
#PCO.diagn.succ <- with(in.succ$full, pcoDiagnostics(Y = Y, D = D, iterations = iter))
PCO.diagn.succ <- readRDS(file.path("Results", "PCO_diagnostics_succ_20240508.RDS"))
#PCO.diagn.fail <- with(in.fail$full, pcoDiagnostics(Y = Y, D = D, iterations = iter))
PCO.diagn.fail <- readRDS(file.path("Results", "PCO_diagnostics_fail_20240508.RDS"))
### Save these things becasue they take several minutes to run
#saveRDS(object = PCO.diagn.total, file = file.path("Results", paste("PCO_diagnostics_total_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = "")))
#saveRDS(object = PCO.diagn.succ, file = file.path("Results", paste("PCO_diagnostics_succ_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = "")))
#saveRDS(object = PCO.diagn.fail, file = file.path("Results", paste("PCO_diagnostics_fail_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = "")))

## Number of nontrivial axes
### Broken-stick method
sum(with(PCO.diagn.total, real > broken.stick))
sum(with(PCO.diagn.succ, real > broken.stick))
sum(with(PCO.diagn.fail, real > broken.stick))
### Bootstrap Eigenvalue method
with(PCO.diagn.total, lower.boot[1:346] > upper.boot[2:347])[1:20]
with(PCO.diagn.succ, lower.boot[1:306] > upper.boot[2:307])[1:20]
with(PCO.diagn.fail, lower.boot[1:306] > upper.boot[2:307])[1:20]
### Permutation method based on the full set of PCO axes
with(PCO.diagn.total, real > upper.p1)[1:20]
with(PCO.diagn.succ, real > upper.p1)[1:20]
with(PCO.diagn.fail, real > upper.p1)[1:20]
### Conditional Permutation method
with(PCO.diagn.total, conditional > upper.p2)[1:20]
with(PCO.diagn.succ, conditional > upper.p2)[1:20]
with(PCO.diagn.fail, conditional > upper.p2)[1:20]

## Plot the diagnostics
### Make the diagnostic plots
plot.diagn.total <- pcoDiagPlots(pco.diag = PCO.diagn.total, max.axes = 20, size.point = sizePoint, size.text = sizeText, size.lab = sizeLab)
plot.diagn.succ <- pcoDiagPlots(pco.diag = PCO.diagn.succ, max.axes = 20, size.point = sizePoint, size.text = sizeText, size.lab = sizeLab)
plot.diagn.fail <- pcoDiagPlots(pco.diag = PCO.diagn.fail, max.axes = 20, size.point = sizePoint, size.text = sizeText, size.lab = sizeLab)
### View the diagnostic plots
plot.diagn.total
plot.diagn.succ
plot.diagn.fail
### Save the plots
w.diagn <- 6.5
h.diagn <- 6.5
ggsave(filename = file.path("Results", paste("diagnostics_total_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = plot.diagn.total, device = "tiff", width = w.diagn, height = h.diagn, dpi = 600, compression = "lzw")
ggsave(filename = file.path("Results", paste("diagnostics_succ_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = plot.diagn.succ, device = "tiff", width = w.diagn, height = h.diagn, dpi = 600, compression = "lzw")
ggsave(filename = file.path("Results", paste("diagnostics_fail_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = plot.diagn.fail, device = "tiff", width = w.diagn, height = h.diagn, dpi = 600, compression = "lzw")

## Plot the raw PCOs
### Make the PCO plots
plot.PCO.total.12 <- pcoPlots(ds = PCO.diagn.total, PCO = pco.total$full$df, X = "PCO1", Y = "PCO2", add.fitted = FALSE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.total.13 <- pcoPlots(ds = PCO.diagn.total, PCO = pco.total$full$df, X = "PCO1", Y = "PCO3", add.fitted = FALSE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.total.23 <- pcoPlots(ds = PCO.diagn.total, PCO = pco.total$full$df, X = "PCO2", Y = "PCO3", add.fitted = FALSE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.succ.12 <- pcoPlots(ds = PCO.diagn.succ, PCO = pco.succ$full$df, X = "PCO1", Y = "PCO2", add.fitted = FALSE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.succ.13 <- pcoPlots(ds = PCO.diagn.succ, PCO = pco.succ$full$df, X = "PCO1", Y = "PCO3", add.fitted = FALSE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.succ.23 <- pcoPlots(ds = PCO.diagn.succ, PCO = pco.succ$full$df, X = "PCO2", Y = "PCO3", add.fitted = FALSE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.fail.12 <- pcoPlots(ds = PCO.diagn.fail, PCO = pco.fail$full$df, X = "PCO1", Y = "PCO2", add.fitted = FALSE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.fail.13 <- pcoPlots(ds = PCO.diagn.fail, PCO = pco.fail$full$df, X = "PCO1", Y = "PCO3", add.fitted = FALSE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.fail.23 <- pcoPlots(ds = PCO.diagn.fail, PCO = pco.fail$full$df, X = "PCO2", Y = "PCO3", add.fitted = FALSE, size.text = sizeText, size.lab = sizeLab)
### View the PCO plots
plot.PCO.total.12
plot.PCO.total.13
plot.PCO.total.23
plot.PCO.succ.12
plot.PCO.succ.13
plot.PCO.succ.23
plot.PCO.fail.12
plot.PCO.fail.13
plot.PCO.fail.23

## Make a Determination on number of PCO axes
### For the Total detections model, 3 axes are likely enough to explain the variation in beta diversity for the full dataset
### For the Successful detections model, 3 axes are likely enough to explain the variation in beta diversity for the full dataset
### For the Failed detections model, 2 axes are likely enough to explain the variation in beta diversity for the full dataset


# Step 4: Conduct regression analysis on PCO axes ####
## Run the regressions on the full data
reg.total <- pbapply::pblapply(1:length(in.total), function(i){runRegression(input = in.total[[i]], PCO = pco.total[[i]], method = "hand")})
names(reg.total) <- names(in.total)
reg.succ <- pbapply::pblapply(1:length(in.succ), function(i){runRegression(input = in.succ[[i]], PCO = pco.succ[[i]], method = "hand")})
names(reg.succ) <- names(in.succ)
reg.fail <- pbapply::pblapply(1:length(in.fail), function(i){runRegression(input = in.fail[[i]], PCO = pco.fail[[i]], method = "hand")})
names(reg.fail) <- names(in.fail)

## Examine model diagnostics
reg.total$full$validation[1:3,]
reg.succ$full$validation[1:3,]
reg.fail$full$validation[1:2,]

## Run the regressions on the detrended data
reg.total.cor <- lapply(1:length(in.total.cor), function(i){runRegression(input = in.total.cor[[i]], PCO = pco.total.cor[[i]], method = "lm")})
names(reg.total.cor) <- names(in.total.cor)
reg.succ.cor <- lapply(1:length(in.succ.cor), function(i){runRegression(input = in.succ.cor[[i]], PCO = pco.succ.cor[[i]], method = "lm")})
names(reg.succ.cor) <- names(in.succ.cor)
reg.fail.cor <- lapply(1:length(in.fail.cor), function(i){runRegression(input = in.fail.cor[[i]], PCO = pco.fail.cor[[i]], method = "hand")})
names(reg.fail.cor) <- names(in.fail.cor)

## Examine the model diagnostics
reg.total.cor$full$validation[1:3,]
reg.succ.cor$full$validation[1:3,]
reg.fail.cor$full$validation[1:2,]

# Step 5: Diagnostics of the regression analyses from the non-detrended data ####
## Plot the original PCOs and fitted PCOs
### Make the plots
plot.PCO.total.12.fit <- pcoPlots(ds = PCO.diagn.total, PCO.fit = reg.total$full$fitted, PCO = pco.total$full$df, X = "PCO1", Y = "PCO2", add.fitted = TRUE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.total.13.fit <- pcoPlots(ds = PCO.diagn.total, PCO.fit = reg.total$full$fitted, PCO = pco.total$full$df, X = "PCO1", Y = "PCO3", add.fitted = TRUE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.total.23.fit <- pcoPlots(ds = PCO.diagn.total, PCO.fit = reg.total$full$fitted, PCO = pco.total$full$df, X = "PCO2", Y = "PCO3", add.fitted = TRUE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.succ.12.fit <- pcoPlots(ds = PCO.diagn.succ, PCO.fit = reg.succ$full$fitted, PCO = pco.succ$full$df, X = "PCO1", Y = "PCO2", add.fitted = TRUE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.succ.13.fit <- pcoPlots(ds = PCO.diagn.succ, PCO.fit = reg.succ$full$fitted, PCO = pco.succ$full$df, X = "PCO1", Y = "PCO3", add.fitted = TRUE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.succ.23.fit <- pcoPlots(ds = PCO.diagn.succ, PCO.fit = reg.succ$full$fitted, PCO = pco.succ$full$df, X = "PCO2", Y = "PCO3", add.fitted = TRUE, size.text = sizeText, size.lab = sizeLab)
plot.PCO.fail.12.fit <- pcoPlots(ds = PCO.diagn.fail, PCO.fit = reg.fail$full$fitted, PCO = pco.fail$full$df, X = "PCO1", Y = "PCO2", add.fitted = TRUE, size.text = sizeText, size.lab = sizeLab)
#plot.PCO.fail.13.fit <- pcoPlots(ds = PCO.diagn.fail, PCO.fit = reg.fail$fitted, PCO = pco.fail$df, X = "PCO1", Y = "PCO3", add.fitted = TRUE, size.text = sizeText, size.lab = sizeLab)
#plot.PCO.fail.23.fit <- pcoPlots(ds = PCO.diagn.fail, PCO.fit = reg.fail$fitted, PCO = pco.fail$df, X = "PCO2", Y = "PCO3", add.fitted = TRUE, size.text = sizeText, size.lab = sizeLab)
### View the plots
plot.PCO.total.12.fit
plot.PCO.total.13.fit
plot.PCO.total.23.fit
plot.PCO.succ.12.fit
plot.PCO.succ.13.fit
plot.PCO.succ.23.fit
plot.PCO.fail.12.fit
#plot.PCO.fail.13.fit
#plot.PCO.fail.23.fit

## Compute the strength of the association in 2D PCO space
### Strength of the association in 2d PCO space can be examined by calculating the Mantel correlation between the two 2d Euclidean distance matrices.
### This is just the Mantel correlation between the original PCO points and their fitted values
### Total
with(reg.total$full, cor(dist(as.matrix(data[,1:3])), dist(as.matrix(fitted[,1:3]))))
### Successful
with(reg.succ$full, cor(dist(as.matrix(data[,1:3])), dist(as.matrix(fitted[,1:3]))))
### Failure
with(reg.fail$full, cor(dist(as.matrix(data[,1:3])), dist(as.matrix(fitted[,1:3]))))

## Test for identically distributed errors for assumption of the permutation test
### Autocorrelation functions
#### Calculate autocorrelation functions of the residuals
acf.total <- with(reg.total$full, list(PC1 = stats::acf(residuals[,1], plot = FALSE), 
                                       PC2 = stats::acf(residuals[,2], plot = FALSE), 
                                       PC3 = stats::acf(residuals[,3], plot = FALSE)))
acf.succ <- with(reg.succ$full, list(PC1 = stats::acf(residuals[,1], plot = FALSE), 
                                     PC2 = stats::acf(residuals[,2], plot = FALSE), 
                                     PC3 = stats::acf(residuals[,3], plot = FALSE)))
acf.fail <- with(reg.fail$full, list(PC1 = stats::acf(residuals[,1], plot = FALSE), 
                                     PC2 = stats::acf(residuals[,2], plot = FALSE)))
#### Combine the acf functions for plotting purposes
acf.all <- rbind(with(acf.total$PC1, data.frame(lag = lag, acf = acf, dataset = "total", PC = "PC1")), 
                 with(acf.total$PC2, data.frame(lag = lag, acf = acf, dataset = "total", PC = "PC2")), 
                 with(acf.total$PC3, data.frame(lag = lag, acf = acf, dataset = "total", PC = "PC3")), 
                 with(acf.succ$PC1, data.frame(lag = lag, acf = acf, dataset = "succ", PC = "PC1")), 
                 with(acf.succ$PC2, data.frame(lag = lag, acf = acf, dataset = "succ", PC = "PC2")), 
                 with(acf.succ$PC3, data.frame(lag = lag, acf = acf, dataset = "succ", PC = "PC3")), 
                 with(acf.fail$PC1, data.frame(lag = lag, acf = acf, dataset = "fail", PC = "PC1")), 
                 with(acf.fail$PC2, data.frame(lag = lag, acf = acf, dataset = "fail", PC = "PC2")))
acf.all$dataset <- factor(acf.all$dataset, levels = c("total", "succ", "fail"), labels = c("Total", "Successful", "Unsuccessful"))
acf.all$PC <- factor(acf.all$PC, levels = c("PC1", "PC2", "PC3"), labels = c("PC1", "PC2", "PC3"))
### Calculate Mantel Correlogram of the residuals
#### Do initial calculations to calculate the mantel correlograms
mantel.in <- list(total = with(reg.total$full, list(resid.d = sqrt(abs(dist(residuals[,validation$pos.eigen], method = "euclidean")^2 - 
                                                                       dist(residuals[,!validation$pos.eigen], method = "euclidean")^2)), 
                                                    times.d = dist(seq(1, validation$n[1]), method = "euclidean"))), 
                  succ = with(reg.succ$full, list(resid.d = sqrt(abs(dist(residuals[,validation$pos.eigen], method = "euclidean")^2 - 
                                                                     dist(residuals[,!validation$pos.eigen], method = "euclidean")^2)), 
                                                  times.d = dist(seq(1, validation$n[1]), method = "euclidean"))), 
                  fail = with(reg.fail$full, list(resid.d = sqrt(abs(dist(residuals[,validation$pos.eigen], method = "euclidean")^2 - 
                                                                     dist(residuals[,!validation$pos.eigen], method = "euclidean")^2)), 
                                                  times.d = dist(seq(1, validation$n[1]), method = "euclidean"))))
#### Compute the mantel correlograms
mantel.r <- pbapply::pblapply(mantel.in, function(x){do.call(what = vegan::mantel.correlog, 
                                                             args = list(D.eco = x$resid.d, D.geo = x$times.d, 
                                                                         nperm = 999, mult = "bonferroni", progressive = TRUE))})
#### Prepare the outputs for plotting
mantel.r.all <- rbind(with(mantel.r$total, data.frame(class.index = mantel.res[!is.na(mantel.res[,3]),1], Mantel.cor = mantel.res[!is.na(mantel.res[,3]),3], p.value = mantel.res[!is.na(mantel.res[,3]),4], dataset = "total")), 
                      with(mantel.r$succ, data.frame(class.index = mantel.res[!is.na(mantel.res[,3]),1], Mantel.cor = mantel.res[!is.na(mantel.res[,3]),3], p.value = mantel.res[!is.na(mantel.res[,3]),4], dataset = "succ")), 
                      with(mantel.r$fail, data.frame(class.index = mantel.res[!is.na(mantel.res[,3]),1], Mantel.cor = mantel.res[!is.na(mantel.res[,3]),3], p.value = mantel.res[!is.na(mantel.res[,3]),4], dataset = "fail")))
mantel.r.all$dataset <- factor(mantel.r.all$dataset, levels = c("total", "succ", "fail"), labels = c("Total", "Successful", "Unsuccessful"))
mantel.r.all$sig <- factor(ifelse(mantel.r.all$p.value < 0.05, "sig", "not"), levels = c("not", "sig"), labels = c("Not Significant", "Significant"))
### Plot the relationship between the fitted values and residuals
#### Combine the fitted and residual values for plotting
fit.resid <- rbind(with(reg.total$full, 
                        rbind(data.frame(fitted = fitted[,1], residuals = residuals[,1], dataset = "total", PC = "PC1"), 
                              data.frame(fitted = fitted[,2], residuals = residuals[,2], dataset = "total", PC = "PC2"), 
                              data.frame(fitted = fitted[,3], residuals = residuals[,3], dataset = "total", PC = "PC3"))), 
                   with(reg.succ$full, 
                        rbind(data.frame(fitted = fitted[,1], residuals = residuals[,1], dataset = "succ", PC = "PC1"), 
                              data.frame(fitted = fitted[,2], residuals = residuals[,2], dataset = "succ", PC = "PC2"), 
                              data.frame(fitted = fitted[,3], residuals = residuals[,3], dataset = "succ", PC = "PC3"))), 
                   with(reg.fail$full, 
                        rbind(data.frame(fitted = fitted[,1], residuals = residuals[,1], dataset = "fail", PC = "PC1"), 
                              data.frame(fitted = fitted[,2], residuals = residuals[,2], dataset = "fail", PC = "PC2"))))
fit.resid$dataset <- factor(fit.resid$dataset, levels = c("total", "succ", "fail"), labels = c("Total", "Successful", "Unsuccessful"))
fit.resid$PC <- factor(fit.resid$PC, levels = c("PC1", "PC2", "PC3"), labels = c("PC1", "PC2", "PC3"))
### Make plots of the residuals
#### Autocorrelation functions
plot.reg.acf <- ggplot(acf.all, mapping = aes(x = lag, y = acf, xend = lag, yend = 0)) + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(2/sqrt(in.total$full$N), -2/sqrt(in.total$full$N)), color = "blue", linetype = 2) + 
  geom_hline(yintercept = c(2/sqrt(in.succ$full$N), -2/sqrt(in.succ$full$N)), color = "green", linetype = 2) + 
  geom_segment() + 
  xlab("Lag") + ylab("ACF") + 
  facet_grid(rows = vars(PC), cols = vars(dataset)) + 
  theme_bw()
plot.reg.acf
#### Mantel correlograms
plot.reg.mantel <- ggplot(data = mantel.r.all, mapping = aes(x = class.index, y = Mantel.cor)) + 
  geom_hline(yintercept = 0, color = "red", linewidth = 1.0) + 
  geom_point(mapping = aes(shape = sig), size = 3) + 
  geom_line() + 
  scale_shape_manual("Significance", values = c(1,16)) + 
  xlab("Distance class index") + ylab("Mantel correlation") + 
  facet_grid(rows = vars(dataset)) + 
  theme_bw()
plot.reg.mantel
#### Fitted values vs. residuals
plot.reg.fitted <- ggplot(data = fit.resid, mapping = aes(x = fitted, y = residuals)) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_point() + 
  facet_grid(rows = vars(PC), cols = vars(dataset)) + 
  theme_bw()
plot.reg.fitted
### Save the plots
w.reg <- 6.5
h.reg <- 6.5
ggsave(filename = file.path("Results", paste("reg_diag_acf_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = plot.reg.acf, device = "tiff", width = w.reg, height = h.reg, dpi = 600, compression = "lzw")
ggsave(filename = file.path("Results", paste("reg_diag_mantel_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = plot.reg.mantel, device = "tiff", width = w.reg, height = h.reg, dpi = 600, compression = "lzw")
ggsave(filename = file.path("Results", paste("reg_diag_fitted_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = plot.reg.fitted, device = "tiff", width = w.reg, height = h.reg, dpi = 600, compression = "lzw")


# Step 6: Diagnostics of the regression analyses from the detrended data (For comparison purposes) ####
## Compute Mantel correlation between between the original PCO points and their fitted values
### Total
with(reg.total.cor$full, cor(dist(as.matrix(data[,1:3])), dist(as.matrix(fitted[,1:3]))))
### Successful
with(reg.succ.cor$full, cor(dist(as.matrix(data[,1:3])), dist(as.matrix(fitted[,1:3]))))
### Failure
with(reg.fail.cor$full, cor(dist(as.matrix(data[,1:3])), dist(as.matrix(fitted[,1:3]))))

## Test for independent and identically distributed samples
### Independence test via autocorrelation functions
acf.total.cor <- with(reg.total.cor$full, list(PC1 = stats::acf(residuals[,1], plot = FALSE), 
                                               PC2 = stats::acf(residuals[,2], plot = FALSE), 
                                               PC3 = stats::acf(residuals[,3], plot = FALSE)))
acf.succ.cor <- with(reg.succ.cor$full, list(PC1 = stats::acf(residuals[,1], plot = FALSE), 
                                             PC2 = stats::acf(residuals[,2], plot = FALSE), 
                                             PC3 = stats::acf(residuals[,3], plot = FALSE)))
acf.fail.cor <- with(reg.fail.cor$full, list(PC1 = stats::acf(residuals[,1], plot = FALSE), 
                                             PC2 = stats::acf(residuals[,2], plot = FALSE)))
acf.all.cor <- rbind(with(acf.total.cor$PC1, data.frame(lag = lag, acf = acf, dataset = "total", PC = "PC1")), 
                     with(acf.total.cor$PC2, data.frame(lag = lag, acf = acf, dataset = "total", PC = "PC2")), 
                     with(acf.total.cor$PC3, data.frame(lag = lag, acf = acf, dataset = "total", PC = "PC3")), 
                     with(acf.succ.cor$PC1, data.frame(lag = lag, acf = acf, dataset = "succ", PC = "PC1")), 
                     with(acf.succ.cor$PC2, data.frame(lag = lag, acf = acf, dataset = "succ", PC = "PC2")), 
                     with(acf.succ.cor$PC3, data.frame(lag = lag, acf = acf, dataset = "succ", PC = "PC3")), 
                     with(acf.fail.cor$PC1, data.frame(lag = lag, acf = acf, dataset = "fail", PC = "PC1")), 
                     with(acf.fail.cor$PC2, data.frame(lag = lag, acf = acf, dataset = "fail", PC = "PC2")))
acf.all.cor$dataset <- factor(acf.all.cor$dataset, levels = c("total", "succ", "fail"), labels = c("Total", "Successful", "Unsuccessful"))
acf.all.cor$PC <- factor(acf.all.cor$PC, levels = c("PC1", "PC2", "PC3"), labels = c("PC1", "PC2", "PC3"))
### Independence test via the Mantel correlogram
#### Do initial calculations to calculate the mantel correlograms
mantel.in.cor <- list(total = with(reg.total.cor$full, list(resid.d = sqrt(abs(dist(residuals[,validation$pos.eigen], method = "euclidean")^2 - 
                                                                               dist(residuals[,!validation$pos.eigen], method = "euclidean")^2)), 
                                                            times.d = dist(seq(1, validation$n[1]), method = "euclidean"))), 
                      succ = with(reg.succ.cor$full, list(resid.d = sqrt(abs(dist(residuals[,validation$pos.eigen], method = "euclidean")^2 - 
                                                                             dist(residuals[,!validation$pos.eigen], method = "euclidean")^2)), 
                                                          times.d = dist(seq(1, validation$n[1]), method = "euclidean"))), 
                      fail = with(reg.fail.cor$full, list(resid.d = sqrt(abs(dist(residuals[,validation$pos.eigen], method = "euclidean")^2 - 
                                                                             dist(residuals[,!validation$pos.eigen], method = "euclidean")^2)), 
                                                          times.d = dist(seq(1, validation$n[1]), method = "euclidean"))))
#### Compute the mantel correlograms
mantel.r.cor <- pbapply::pblapply(mantel.in.cor, function(x){do.call(what = vegan::mantel.correlog, 
                                                                     args = list(D.eco = x$resid.d, D.geo = x$times.d, 
                                                                                 nperm = 999, mult = "bonferroni", progressive = TRUE))})
#### Prepare the outputs for plotting
mantel.r.all.cor <- rbind(with(mantel.r.cor$total, data.frame(class.index = mantel.res[!is.na(mantel.res[,3]),1], Mantel.cor = mantel.res[!is.na(mantel.res[,3]),3], p.value = mantel.res[!is.na(mantel.res[,3]),4], dataset = "total")), 
                          with(mantel.r.cor$succ, data.frame(class.index = mantel.res[!is.na(mantel.res[,3]),1], Mantel.cor = mantel.res[!is.na(mantel.res[,3]),3], p.value = mantel.res[!is.na(mantel.res[,3]),4], dataset = "succ")), 
                          with(mantel.r.cor$fail, data.frame(class.index = mantel.res[!is.na(mantel.res[,3]),1], Mantel.cor = mantel.res[!is.na(mantel.res[,3]),3], p.value = mantel.res[!is.na(mantel.res[,3]),4], dataset = "fail")))
mantel.r.all.cor$dataset <- factor(mantel.r.all.cor$dataset, levels = c("total", "succ", "fail"), labels = c("Total", "Successful", "Unsuccessful"))
mantel.r.all.cor$sig <- factor(ifelse(mantel.r.all.cor$p.value < 0.05, "sig", "not"), levels = c("not", "sig"), labels = c("Not Significant", "Significant"))
### Identically distributed errors via plots of fitted vs. residuals
fit.resid.cor <- rbind(with(reg.total.cor$full, 
                        rbind(data.frame(fitted = fitted[,1], residuals = residuals[,1], dataset = "total", PC = "PC1"), 
                              data.frame(fitted = fitted[,2], residuals = residuals[,2], dataset = "total", PC = "PC2"), 
                              data.frame(fitted = fitted[,3], residuals = residuals[,3], dataset = "total", PC = "PC3"))), 
                   with(reg.succ.cor$full, 
                        rbind(data.frame(fitted = fitted[,1], residuals = residuals[,1], dataset = "succ", PC = "PC1"), 
                              data.frame(fitted = fitted[,2], residuals = residuals[,2], dataset = "succ", PC = "PC2"), 
                              data.frame(fitted = fitted[,3], residuals = residuals[,3], dataset = "succ", PC = "PC3"))), 
                   with(reg.fail.cor$full, 
                        rbind(data.frame(fitted = fitted[,1], residuals = residuals[,1], dataset = "fail", PC = "PC1"), 
                              data.frame(fitted = fitted[,2], residuals = residuals[,2], dataset = "fail", PC = "PC2"))))
fit.resid.cor$dataset <- factor(fit.resid.cor$dataset, levels = c("total", "succ", "fail"), labels = c("Total", "Successful", "Unsuccessful"))
fit.resid.cor$PC <- factor(fit.resid.cor$PC, levels = c("PC1", "PC2", "PC3"), labels = c("PC1", "PC2", "PC3"))
### Make plots of the residuals
#### Autocorrelation functions
ggplot(acf.all.cor, mapping = aes(x = lag, y = acf, xend = lag, yend = 0)) + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(2/sqrt(in.total$full$N), -2/sqrt(in.total$full$N)), color = "blue", linetype = 2) + 
  geom_hline(yintercept = c(2/sqrt(in.succ$full$N), -2/sqrt(in.succ$full$N)), color = "green", linetype = 2) + 
  geom_segment() + 
  xlab("Lag") + ylab("ACF") + 
  facet_grid(rows = vars(PC), cols = vars(dataset)) + 
  theme_bw()
#### Mantel correlograms
ggplot(data = mantel.r.all.cor, mapping = aes(x = class.index, y = Mantel.cor)) + 
  geom_hline(yintercept = 0, color = "red", linewidth = 1.0) + 
  geom_point(mapping = aes(shape = sig), size = 3) + 
  geom_line() + 
  scale_shape_manual("Significance", values = c(1,16)) + 
  xlab("Distance class index") + ylab("Mantel correlation") + 
  facet_grid(rows = vars(dataset)) + 
  theme_bw()
#### Fitted values vs. residuals
ggplot(data = fit.resid.cor, mapping = aes(x = fitted, y = residuals)) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_point() + 
  facet_grid(rows = vars(PC), cols = vars(dataset)) + 
  theme_bw()


# Step 7: Make a prediction from the regression equation ####
## Run the prediction function
pred.total <- predictReg(reg = reg.total, axes = 1:3, describe = FALSE)
pred.succ <- predictReg(reg = reg.succ, axes = 1:3, describe = FALSE)
pred.fail <- predictReg(reg = reg.fail, axes = 1:2, describe = TRUE)

## Examine the Mantel Correlations
pred.total$mantel[-c(19:51),]
pred.succ$mantel[-c(19:51),]
pred.fail$mantel[-c(19:51),]
### Summary statistics
apply(pred.total$mantel[-c(19:51),], 2, function(x){c(mean = mean(x), sd = sd(x), min = min(x), max = max(x))})
apply(pred.succ$mantel[-c(19:51),], 2, function(x){c(mean = mean(x), sd = sd(x), min = min(x), max = max(x))})
apply(pred.fail$mantel[-c(19:51),], 2, function(x){c(mean = mean(x), sd = sd(x), min = min(x), max = max(x))})
### Save the Mantel correlations 
openxlsx::write.xlsx(x = list(total = pred.total$mantel[-c(19:51),], succ = pred.succ$mantel[-c(19:51),], fail = pred.fail$mantel[-c(19:51),]), 
                     file = file.path("Results", paste("mantel_correlations_all_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = "")), rowNames = TRUE)

# Step 8: Plot the predicted values ####
## Prepare data for creation of plots
### Combine the prediction outputs
pred.all <- c(total = pred.total, succ = pred.succ, fail = pred.fail)
### Create a plotting object for site
plot.site <- do.call(dplyr::bind_rows, lapply(grep("pred", names(pred.all)), function(i){
  x <- pred.all[[i]]
  name1 <- sub(".pred", "", names(pred.all)[i])
  x1 <- lapply(2:19, function(j){
    y <- x[[j]]
    name2 <- sub("drop_", "", names(x)[j])
    
    data.pred <- data.frame(ds = "Raw PCO", drop = name2, label = 1:nrow(y$data.pred), y$data.pred)
    yhat.pred <- data.frame(ds = "Predicted", drop = name2, label = 1:nrow(y$Yhat.pred), y$Yhat.pred)
    orig.fit <- data.frame(ds = "Fitted", drop = name2, label = 1:nrow(y$orig.fit), y$orig.fit)
    
    y1 <- dplyr::bind_rows(data.pred, yhat.pred, orig.fit)
    y1$ds <- factor(y1$ds, levels = c("Raw PCO", "Fitted", "Predicted"), labels = c("Raw PCO Scores", "Fitted from Full Model", "Predicted from Drop Models"))
    row.names(y1) <- NULL
    return(y1)
  })
  x2 <- data.frame(response = factor(name1, levels = c("total", "succ", "fail"), labels = c("Total detections", "Successful crossings", "Unsuccessful crossings")), do.call(rbind, x1))
}))
### Craete a plotting object for month
plot.time <- do.call(dplyr::bind_rows, lapply(grep("pred", names(pred.all)), function(i){
  x <- pred.all[[i]]
  name1 <- sub(".pred", "", names(pred.all)[i])
  x1 <- lapply(40:51, function(j){
    y <- x[[j]]
    name2 <- sub("drop_", "", names(x)[j])
    
    data.pred <- data.frame(ds = "Raw PCO", drop = name2, label = 1:nrow(y$data.pred), y$data.pred)
    yhat.pred <- data.frame(ds = "Predicted", drop = name2, label = 1:nrow(y$Yhat.pred), y$Yhat.pred)
    orig.fit <- data.frame(ds = "Fitted", drop = name2, label = 1:nrow(y$orig.fit), y$orig.fit)
    
    y1 <- dplyr::bind_rows(data.pred, yhat.pred, orig.fit)
    y1$ds <- factor(y1$ds, levels = c("Raw PCO", "Fitted", "Predicted"), labels = c("Raw PCO Scores", "Fitted from Full Model", "Predicted from Drop Models"))
    y1$drop <- factor(y1$drop, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8" , "9", "10", "11"))
    row.names(y1) <- NULL
    return(y1)
  })
  x2 <- data.frame(response = factor(name1, levels = c("total", "succ", "fail"), labels = c("Total detections", "Successful crossings", "Unsuccessful crossings")), do.call(rbind, x1))
}))

## Calculate species scores from the full models
### First combine the full in files
in.full <- list(total = in.total$full, succ = in.succ$full, fail = in.fail$full)
### Calculate the species correlations with the fitted dbRDA values
plot.spec <- do.call(dplyr::bind_rows, lapply(grep("pred", names(pred.all)), function(i){
  name <- names(pred.all)[i]
  
  # Get the fitted PCO scores
  x <- pred.all[[i]]$full$full.fitted
  
  # Extract the species abundance information
  Y1 <- in.full[[grep(sub(".pred", "", name), names(in.full))]]$Y
  Y2 <- Y1[,which(colSums(Y1) > 0)]
  
  # Calculate the Pearson correlation
  cors <- cor(Y2, x)
  
  out <- data.frame(response = factor(sub(".pred", "", name), levels = c("total", "succ", "fail"), labels = c("Total detections", "Successful crossings", "Unsuccessful crossings")), 
                    spec.full = factor(row.names(cors)), 
                    species = NA, 
                    cors, row.names = NULL)
  
  return(out)
  #rm(name, x, Y1, Y2, cors, out)
}))
plot.spec$species <- factor(plot.spec$spec.full, 
                            levels = c("armadillo", "beaver", "bobcat", "cottontail", "coyote", "feral_hog", "fox_squirrel", 
                                       "gray_squirrel", "grey_fox", "jackrabbit", "javelina", "mexican_ground_squirrel", "nilgai", 
                                       "nutria", "ocelot", "opossum", "raccoon", "rodent", "striped_skunk", "unk_mammal", "weasel", 
                                       "white_tailed_deer"), 
                            labels = c("ARMA", "BEAV", "BOBC", "ECOT", "COYO", "FHOG", "FSQR", 
                                       "GSQR", "GFOX", "JACK", "JAVE", "MSQR", "NILG", 
                                       "NUTR", "OCEL", "OPOS", "RACC", "UROD", "SSKU", "UMAM", "WEAS", 
                                       "DEER"))
## Four letter codes based on scientific name
c("DANO", "CACA", "LYRU", "SYFL", "CALA", "SUSC", "SCNI", 
  "SCCA", "URCI", "LECA", "PETA", "SPME", "BOTR", 
  "MYCO", "LEPA", "DIVI", "PRLO", "UROD", "MEME", "UMAM", "MUFR", 
  "ODVI")


## Create plots of raw PCO scores, fitted values from the full model, and predicted values from the drop model
### Define plotting parameters
sizeText <- 3.5
sizeLab <- 12
sizePoint <- 3
### Create the plots
#### Arguments
ARGS.no <- list(ds = plot.site, plot.type = "site", fill = "drop", facet.plot = FALSE, average = FALSE, add.cor = FALSE, size.text = sizeText, size.lab = sizeLab, size.point = sizePoint)
ARGS.cor <- list(ds = plot.site, plot.type = "site", fill = "drop", facet.plot = FALSE, average = FALSE, add.cor = TRUE, cor.file = plot.spec, size.text = sizeText, size.lab = sizeLab, size.point = sizePoint)
#### The plots
plot.site.no.12 <- do.call(predictPlots, args = c(list(X = "PCO1", Y = "PCO2"), ARGS.no))
plot.site.no.13 <- do.call(predictPlots, args = c(list(X = "PCO1", Y = "PCO3"), ARGS.no))
plot.site.no.23 <- do.call(predictPlots, args = c(list(X = "PCO2", Y = "PCO3"), ARGS.no))
plot.site.cor.12 <- do.call(predictPlots, args = c(list(X = "PCO1", Y = "PCO2"), ARGS.cor))
plot.site.cor.13 <- do.call(predictPlots, args = c(list(X = "PCO1", Y = "PCO3"), ARGS.cor))
plot.site.cor.23 <- do.call(predictPlots, args = c(list(X = "PCO2", Y = "PCO3"), ARGS.cor))
### Create the multi-plots
ARGS.arrange <- list(ncol = 2, nrow = 2, labels = "AUTO", font.label = list(size = 12, family = "serif"), common.legend = TRUE, legend = "right")
plot.site.total.12 <- do.call(ggpubr::ggarrange, args = c(list(plotlist = c(plot.site.no.12[1:3], plot.site.cor.12[3])), ARGS.arrange)) + 
  ggpubr::bgcolor("white") + ggpubr::border(color = NA)
plot.site.succ.12 <- do.call(ggpubr::ggarrange, args = c(list(plotlist = c(plot.site.no.12[4:6], plot.site.cor.12[6])), ARGS.arrange)) + 
  ggpubr::bgcolor("white") + ggpubr::border(color = NA)
plot.site.fail.12 <- do.call(ggpubr::ggarrange, args = c(list(plotlist = c(plot.site.no.12[7:9], plot.site.cor.12[9])), ARGS.arrange)) + 
  ggpubr::bgcolor("white") + ggpubr::border(color = NA)
plot.site.total.13 <- do.call(ggpubr::ggarrange, args = c(list(plotlist = c(plot.site.no.13[1:3], plot.site.cor.13[3])), ARGS.arrange)) + 
  ggpubr::bgcolor("white") + ggpubr::border(color = NA)
plot.site.succ.13 <- do.call(ggpubr::ggarrange, args = c(list(plotlist = c(plot.site.no.13[4:6], plot.site.cor.13[6])), ARGS.arrange)) + 
  ggpubr::bgcolor("white") + ggpubr::border(color = NA)
plot.site.total.23 <- do.call(ggpubr::ggarrange, args = c(list(plotlist = c(plot.site.no.23[1:3], plot.site.cor.23[3])), ARGS.arrange)) + 
  ggpubr::bgcolor("white") + ggpubr::border(color = NA)
plot.site.succ.23 <- do.call(ggpubr::ggarrange, args = c(list(plotlist = c(plot.site.no.23[4:6], plot.site.cor.23[6])), ARGS.arrange)) + 
  ggpubr::bgcolor("white") + ggpubr::border(color = NA)
### View the plots
plot.site.total.12
plot.site.succ.12
plot.site.fail.12
plot.site.total.13
plot.site.succ.13
plot.site.total.23
plot.site.succ.23
### Save the plots
ARGS.save <- list(device = "tiff", width = 9, height = 6, dpi = 600, compression = "lzw")
#w.reg <- 15
#h.reg <- 9
do.call(ggsave, args = c(list(filename = file.path("Results", paste("pred_total_12_c_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
                              plot = plot.site.total.12), ARGS.save))
do.call(ggsave, args = c(list(filename = file.path("Results", paste("pred_succ_12_c_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
                              plot = plot.site.succ.12), ARGS.save))
do.call(ggsave, args = c(list(filename = file.path("Results", paste("pred_fail_12_c_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
                              plot = plot.site.fail.12), ARGS.save))
do.call(ggsave, args = c(list(filename = file.path("Results", paste("pred_total_13_c_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
                              plot = plot.site.total.13), ARGS.save))
do.call(ggsave, args = c(list(filename = file.path("Results", paste("pred_succ_13_c_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
                              plot = plot.site.succ.13), ARGS.save))
do.call(ggsave, args = c(list(filename = file.path("Results", paste("pred_total_23_c_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
                              plot = plot.site.total.23), ARGS.save))
do.call(ggsave, args = c(list(filename = file.path("Results", paste("pred_succ_23_c_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
                              plot = plot.site.succ.23), ARGS.save))


## Create plots of fitted vs observed data from the predicted models
### Combine the prediction models
plot.pred.site <- do.call(dplyr::bind_rows, lapply(grep("pred", names(pred.all)), function(i){
  x <- pred.all[[i]]
  name1 <- sub(".pred", "", names(pred.all)[i])
  
  x1 <- lapply(2:19, function(j){
    y <- x[[j]]
    name2 <- sub("drop_", "", names(x)[j])
    
    y1 <- data.frame(drop = name2, label = 1:nrow(y$data.pred), axis = "PCO1", raw = y$data.pred[,1], pred = y$orig.fit[,1], fit = y$Yhat.pred[,1])
    y2 <- data.frame(drop = name2, label = 1:nrow(y$data.pred), axis = "PCO2", raw = y$data.pred[,2], pred = y$orig.fit[,2], fit = y$Yhat.pred[,2])
    y3 <- data.frame(drop = name2, label = 1:nrow(y$data.pred), axis = "PCO3", raw = y$data.pred[,3], pred = y$orig.fit[,3], fit = y$Yhat.pred[,3])
    
    y4 <- rbind(y1, y2, y3)
    y5 <- tidyr::separate(data = y4, col = drop, into = c("Road", "WCS"))
    return(y5)
    #rm(y, name2, y1)
  })
  x2 <- data.frame(response = factor(name1, levels = c("total", "succ", "fail"), labels = c("Total detections", "Successful crossings", "Unsuccessful crossings")), do.call(rbind, x1))
  rownames(x2) <- NULL
  return(x2)
  #rm(x, name1, x1, x2)
}))

ggplot(data = plot.pred.site, mapping = aes(x = raw, y = pred)) + 
  geom_segment(mapping = aes(x = -Inf, xend = Inf, y = -Inf, yend = Inf), color = "grey50", linetype = 2, linewidth = 1.0) + 
  geom_point(mapping = aes(color = WCS, shape = Road)) + 
  coord_fixed() + 
  scale_color_manual("WCS", values = c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255", "#EE9D28")) + 
  scale_shape_manual("Road", values = c(15, 16, 17)) +
  xlab("Observed PCO Scores") + ylab("Predicted PCO Scores") + 
  facet_grid(cols = vars(response), rows = vars(axis)) + 
  theme_bw()
ggsave(filename = file.path("Results", paste("Obs_vs_Pred_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       device = "tiff", width = 6.5, height = 4.5, dpi = 600, compression = "lzw")

## Create plots for presentation purposes
### Plots that do not include species scores
ARGS_pres <- list(ds = plot.site, plot.type = "site", fill = "drop", facet.plot = FALSE, average = FALSE, add.cor = FALSE, cor.file = NULL, size.text = 8, size.lab = 16, size.point = 3)
plot.site.pres <- do.call(predictPlots, args = c(list(X = "PCO1", Y = "PCO2"), ARGS_pres))
plot.site.pres[grep("Total", names(plot.site.pres))]

lapply(grep("Total", names(plot.site.pres)), function(i){
  x1 <- plot.site.pres[[i]]
  name <- sub(".Total detections", "", names(plot.site.pres)[i])
  
  ggsave(filename = file.path("Results", paste("prediction_pres_", name, "_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
         plot = x1, device = "tiff", width = 10, height = 7, dpi = 900, compression = "lzw")
})

ARGS_pres_spec <- list(ds = plot.site, plot.type = "site", fill = "drop", facet.plot = FALSE, average = FALSE, add.cor = TRUE, cor.file = plot.spec, size.text = 8, size.lab = 16, size.point = 3)
plot.site.pres.spec <- do.call(predictPlots, args = c(list(X = "PCO1", Y = "PCO2"), ARGS_pres_spec))
plot.site.pres.spec[grep("Total", names(plot.site.pres.spec))]

lapply(grep("Total", names(plot.site.pres.spec)), function(i){
  x1 <- plot.site.pres.spec[[i]]
  name <- sub(".Total detections", "", names(plot.site.pres.spec)[i])
  
  ggsave(filename = file.path("Results", paste("prediction_pres_spec_", name, "_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
         plot = x1, device = "tiff", width = 10, height = 7, dpi = 900, compression = "lzw")
})


#######################################################################################################################
################################################ OLD AND UNNEEDED CODE ################################################ 
#######################################################################################################################


# Creating unnecessary data in the initial data processing...####
## Need to create a shell with all possible combinations of WCS_Time_ID and Species
shell <- data.frame(WCS_Time_ID = rep(sort(unique(ints4$WCS_Time_ID)), each = length(unique(ints4$Species))), 
                    Species = rep(sort(unique(ints4$Species)), times = length(unique(ints4$WCS_Time_ID))))
shell <- tidyr::separate(shell, col = WCS_Time_ID, sep = "_", into = c("Highway", "Name", "TimeID"), remove = FALSE)
shell$WCS_ID <- paste(shell$Highway, shell$Name, sep = "_")
shell$Sort <- paste(shell$WCS_Time_ID, shell$Species, sep = "__")
shell <- shell[,-c(2,3)]
shell <- shell[,c(1,4,2,3,5)]
shell[1:5,]

## Now merge the interactions data with the shell to get an "all possible values" table
ints5 <- merge.data.frame(shell, ints4[,!(colnames(ints4) %in% c("WCS_Time_ID", "WCS_ID", "TimeID", "Species"))], by = "Sort", all.x = TRUE)
ints5[is.na(ints5)] <- 0
ints5[1:5,]

## For some of the experimental units, we did not assess number of interactions so we need to address this issue
noints <- c("FM106_WCS1_201907", "FM106_WCS1_201908", "FM106_WCS1_201909", "FM106_WCS1_201910", "FM106_WCS1_201911", 
            "FM106_WCS2_201907", "FM106_WCS2_201908", "FM106_WCS2_201909", "FM106_WCS2_201910", "FM106_WCS2_201911", 
            "FM106_WCS3_201907", "FM106_WCS3_201908", "FM106_WCS3_201909", "FM106_WCS3_201910", "FM106_WCS3_201911", 
            "FM106_WCS4_201907", "FM106_WCS4_201908", "FM106_WCS4_201909", "FM106_WCS4_201910", "FM106_WCS4_201911", 
            "FM106_WCS5_201907", "FM106_WCS5_201908", "FM106_WCS5_201909", "FM106_WCS5_201910", "FM106_WCS5_201911", 
            "FM106_WCS6_201907", "FM106_WCS6_201908", "FM106_WCS6_201909", "FM106_WCS6_201910", "FM106_WCS6_201911", 
            "FM106_WCS7_201907", "FM106_WCS7_201908", "FM106_WCS7_201909", "FM106_WCS7_201910", "FM106_WCS7_201911", 
            "FM106_WCS8_201907", "FM106_WCS8_201908", "FM106_WCS8_201909", "FM106_WCS8_201910", "FM106_WCS8_201911")
ints5$noints <- ints5$WCS_Time_ID %in% noints
ints5[1:5,]
ints5$A <- ifelse(ints5$noints, NaN, ints5$A)
ints5$B <- ifelse(ints5$noints, NaN, ints5$B)
ints5$C <- ifelse(ints5$noints, NaN, ints5$C)
ints5$D <- ifelse(ints5$noints, NaN, ints5$D)
ints5$E <- ifelse(ints5$noints, NaN, ints5$E)
ints5$`NA` <- ifelse(!ints5$noints & ints5$`NA` == 0, NaN, ints5$`NA`)
ints5[1:5,]

## Now we can calculate our response variables
ints5$Total <- apply(ints5[,c("A", "B", "C", "D", "E", "NA")], 1, sum, na.rm = TRUE)
ints5$Success <- with(ints5, A)
ints5$Failure <- with(ints5, B + C + E)
ints5$A_total <- with(ints5, Success / Total)
ints5$A_ints <- with(ints5, Success / (Success + Failure))
ints5$DaysInMonth <- lubridate::days_in_month(lubridate::ym(ints5$TimeID))  # Assumes that at least one camera was active at each site for the duration of the study

## Calculate the number of detections, successes, and failures per day as a standardization technique for monthly data
ints5$TotalPD <- with(ints5, Total/DaysInMonth)
ints5$SuccessPD <- with(ints5, Success/DaysInMonth)
ints5$FailurePD <- with(ints5, Failure/DaysInMonth)
ints5[1:5,]

## View a simplified output
ints5[1:5,c("WCS_Time_ID", "Species", "TotalPD", "SuccessPD", "FailurePD", "A_total", "A_ints")]

## Save output for use in Primer
lapply(1:length(ints_env), function(i){openxlsx::write.xlsx(ints_env[[i]], file = file.path("Output_to_Others", paste("Ints_PERMANOVA_", names(ints_env)[i], "_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = "")))})

### Indicator information for covariates
cov_indicators <- openxlsx::read.xlsx(files_env_site, sheet = "Metadata", rowNames = TRUE)

### Exclude spatial and temporal variables (because they need additional processing)
cov.wcs <- do.call(cbind, lapply(c("structural", "environmental", "anthropogenic"), function(i){
  if(i != "structural"){
    c1 <- data.frame(cov.all.list[[i]][,-1])
    colnames(c1) <- colnames(cov.all.list[[i]])[-1]
  }else{
    c1 <- data.frame(cov.all.list[[i]])
    colnames(c1) <- colnames(cov.all.list[[i]])
  }
  return(c1)
}))
### Excluding spatial and temporal
cov.ind.wcs <- with(cov.all.list, data.frame(cov = c(colnames(structural)[-1], colnames(environmental)[-1], colnames(anthropogenic)[-1]), 
                                             group = c(rep("structural", ncol(structural)-1), 
                                                       rep("environmental", ncol(environmental)-1), 
                                                       rep("anthropogenic", ncol(anthropogenic)-1)), 
                                             single = c("openness", "catwalk", "catwalk", "substrate", "substrate", "substrate", "fencing", 
                                                        "Natural", "Water", "Woody", "Precip", 
                                                        "SpeedLimit", "VehicleTraffic", "HumanAct", "DomesticAct", "LivestockAct", "Buildings")))



# Modeling the temporal variability in the data ####
## Using the AEM approach does not reasonably model the variability in time because there is only 1 time
## I'm probably doing it wrong but dbMEM seems to reasonably model the potential temporal structure in the data


# Brute force approaches to various tasks ####
### Structural variables
structural <- envdata_site2[order(envdata_site2$WCS_Time_ID), 
                            c("WCS_Time_ID", 
                              "Openness_m", "CatWalks.No", "CatWalks.Yes", "Substrate.Water", "Substrate.Concrete", "Substrate.Dirt", "Fencing_m")]
### Environmental variables
environ <- envdata_site2[order(envdata_site2$WCS_Time_ID), 
                         c("WCS_Time_ID", 
                           "PropNatural", "PropWater", "PropWoody", "Precip_total_mm")]
### Anthropogenic variables
anthro <- envdata_site2[order(envdata_site2$WCS_Time_ID), 
                        c("WCS_Time_ID", 
                          "SpeedLimit", "VehicleTraffic", "HumanActivity", "DomesticActivity", "LivestockActivity", "PropBuilding")]

## Define which fields are spatial data
spatial <- envdata_site2[order(envdata_site2$WCS_Time_ID),
                         c("WCS_Time_ID", 
                           "UTM_14N_X", "UTM_14N_Y")]

## Define which fields are temporal data
temporal <- envdata_site2[order(envdata_site2$WCS_Time_ID),
                          c("WCS_Time_ID", 
                            "ConMonths_cov")]



# Spatial and temporal autocorrelation functions based on dbMEM
## Spatial variables for the Total analysis
### Define site information as a matrix
xy1 <- as.matrix(spatial[,-1])
dimnames(xy1)[[1]] <- spatial[,1]
xy1[1:5,]
### Compute a Euclidean distance matrix on the cartesian coordinates
d.xy1 <- distance(Y = xy1, measure = "Eucl", trans = "none", adj = 0)
### Examine the minimum spanning tree to determine threshold distance for dbMEM
mst.xy1 <- vegan::spantree(d = d.xy1)
summary(mst.xy1$dist)
### Calculate the distance-based Moran's Eigenvector Maps (dbMEM)
dbmem.xy1 <- vegan::pcnm(dis = d.xy1, threshold = max(mst.xy1$dist), dist.ret = FALSE)
### Visualize the data to determine the spatial scale each dbMEM axis operates on
#### Axes 1-4 are operating on broad scales
#### Axes 5-11 are operating on fine scales
### Extract the relevant data
spatial.out1 <- data.frame(WCS_Time_ID = row.names(xy1), dbmem.xy1$vectors)
colnames(spatial.out1) <- c("WCS_Time_ID", paste("s_", colnames(spatial.out1)[-1], sep = ""))
spatial.out1[1:5,]

## Spatial variables for the Interactions analyses
xy2 <- as.matrix(spatial[factors$IntsAnalysis == 1,-1])
dimnames(xy2)[[1]] <- spatial[factors$IntsAnalysis == 1,1]
xy2[1:5,]
### Compute a Euclidean distance matrix on the cartesian coordinates
d.xy2 <- distance(Y = xy2, measure = "Eucl", trans = "none", adj = 0)
### Examine the minimum spanning tree to determine threshold distance for dbMEM
mst.xy2 <- vegan::spantree(d = d.xy2)
summary(mst.xy2$dist)
### Calculate the distance-based Moran's Eigenvector Maps (dbMEM)
dbmem.xy2 <- vegan::pcnm(dis = d.xy2, threshold = max(mst.xy2$dist), dist.ret = FALSE)
### Visualize the data to determine the spatial scale each dbMEM axis operates on
#### Axes 1-4 are operating on broad scales
#### Axes 5-11 are operating on fine scales
### Extract the relevant data
spatial.out2 <- data.frame(WCS_Time_ID = row.names(xy2), dbmem.xy2$vectors)
colnames(spatial.out2) <- c("WCS_Time_ID", paste("s_", colnames(spatial.out2)[-1], sep = ""))
spatial.out2[1:5,]

## Temporal variables for the Total analysis
###Define the temporal variables as a matrix
mo1 <- as.matrix(temporal[,-1])
dimnames(mo1) <- list(temporal[,1], colnames(temporal)[-1])
mo1[1:5,]
### Compute a distance matrix based on Euclidean distance
d.mo1 <- distance(Y = mo1, measure = "Eucl", trans = "none", adj = 0)
### Calculate the distance-based Moran's Eigenvector Maps (dbMEM)
dbmem.mo1 <- vegan::pcnm(dis = d.mo1, dist.ret = FALSE)
dbmem.mo1$vectors[1:5,]
### Visualize the data to determine the temporal scale each axis operates on
#### Axes 1-12 are operating on the majority of sites and time periods
#### Axes 13-22 are operating solely on WCS3A in the early time periods (< -10)
### Extract the relevant data
temporal.out1 <- data.frame(WCS_Time_ID = row.names(mo1), dbmem.mo1$vectors)
colnames(temporal.out1) <- c("WCS_Time_ID", paste("t_", colnames(temporal.out1)[-1], sep = ""))
temporal.out1[1:5,]

## Temporal variables for the Interactions analyses
###Define the temporal variables as a matrix
mo2 <- as.matrix(temporal[factors$IntsAnalysis == 1,-1])
dimnames(mo2) <- list(temporal[factors$IntsAnalysis == 1,1], colnames(temporal)[-1])
mo2[1:5,]
### Compute a distance matrix based on Euclidean distance
d.mo2 <- distance(Y = mo2, measure = "Eucl", trans = "none", adj = 0)
### Calculate the distance-based Moran's Eigenvector Maps (dbMEM)
dbmem.mo2 <- vegan::pcnm(dis = d.mo2, dist.ret = FALSE)
dbmem.mo2$vectors[1:5,]
### Visualize the data to determine the temporal scale each axis operates on
#### Axes 1-12 are operating on the majority of sites and time periods
#### Axes 13-22 are operating solely on WCS3A in the early time periods (< -10)
### Extract the relevant data
temporal.out2 <- data.frame(WCS_Time_ID = row.names(mo2), dbmem.mo2$vectors)
colnames(temporal.out2) <- c("WCS_Time_ID", paste("t_", colnames(temporal.out2)[-1], sep = ""))
temporal.out2[1:5,]

# Step 1: Load the required data
## Species data
### Total (data is already square root transformed)
Y.total.all <- openxlsx::read.xlsx(files_reg[grep("total", files_reg)], rowNames = TRUE)
Y.total <- as.matrix(Y.total.all[,c(1:22)])
#Y.total <- as.matrix(Y.total.all[Y.total.all$ConMonths == 11,c(1,3:8,10:13,15:17,19:22)])
which(colSums(Y.total) == 0)
### Successful (data is already square root transformed)
Y.succ.all <- openxlsx::read.xlsx(files_reg[grep("succ", files_reg)], rowNames = TRUE)
Y.succ <- as.matrix(Y.succ.all[,c(1:22)])
#Y.succ <- as.matrix(Y.succ.all[Y.succ.all$ConMonths == 11,c(1,3:6,10:11,15:17,19,21:22)])
which(colSums(Y.succ) == 0)
### Failure (data is already square root transformed)
Y.fail.all <- openxlsx::read.xlsx(files_reg[grep("fail", files_reg)], rowNames = TRUE)
Y.fail <- as.matrix(Y.fail.all[,c(1:22)])
#Y.fail <- as.matrix(Y.fail.all[,c(1,3:8,10:17,19:22)])  # Excluding species with no data does not affect the PCO computation
#Y.fail <- as.matrix(Y.fail.all[Y.fail.all$ConMonths == 11,c(1,3:5,8,11,16:17,19,21:22)])
which(colSums(Y.fail) == 0)

## Covariates
### Total
X.total.all <- openxlsx::read.xlsx(files_reg[grep("Total", files_reg)], rowNames = TRUE, rows = 2:349)
X.total <- as.matrix(X.total.all[,c(1,3,5:15,36)])
#X.total <- as.matrix(X.total.all[X.total.all$ConMonths == 11,c(1,3,5:15,36)])
### Successful
X.succ.all <- openxlsx::read.xlsx(files_reg[grep("Success", files_reg)], rowNames = TRUE, rows = 2:309)
X.succ <- as.matrix(X.succ.all[,c(1,3,5:17,39)])
#X.succ <- as.matrix(X.succ.all[X.succ.all$ConMonths == 11,c(1,3,5:17,39)])
### Failure
X.fail.all <- openxlsx::read.xlsx(files_reg[grep("Fail", files_reg)], rowNames = TRUE, rows = 2:309)
X.fail <- as.matrix(X.fail.all[,c(1,3,5:15,34)])
#X.fail <- as.matrix(X.fail.all[X.fail.all$ConMonths == 11,c(1,3,5:15,34)])

# Step 2: Calculate an appropriate dissimilarity matrix
## Total
D.total <- distance(Y.total, measure = "BC", trans = "none", adj = 0.03)
N.total <- dim(D.total)[1]

## Successful
D.succ <- distance(Y.succ, measure = "BC", trans = "none", adj = 0.03)
N.succ <- dim(D.succ)[1]

## Failure
D.fail <- distance(Y.fail, measure = "BC", trans = "none", adj = 0.03)
N.fail <- dim(D.fail)[1]


# Step 3: Conduct PCO analysis
## Total
### Conduct PCO Analysis
PCO.result.total <- pco(D = D.total, varplot = TRUE)
### Extract eigenvectors
PCO.total <- PCO.result.total$vectors
colnames(PCO.total) <- paste("PCO", 1:ncol(PCO.total), sep = "")
### Convert the PCO matrix to a data.frame
PCO.df.total <- data.frame(site = row.names(Y.total), PCO.total)
PCO.df.total[1:5,1:5]
### Identify positive and negative PCO axes
pos.neg.total <- PCO.result.total$values >= 0.0 
### Percent of variability explained
round(100 * PCO.result.total$values/sum(PCO.result.total$values), digits = 2)

## Successful
### Conduct PCO Analysis
PCO.result.succ <- pco(D = D.succ, varplot = TRUE)
### Standardized PCO scores to var = lambda
PCO.succ <- PCO.result.succ$vectors
colnames(PCO.succ) <- paste("PCO", 1:ncol(PCO.succ), sep = "")
### Convert the PCO matrix to a data.frame
PCO.df.succ <- data.frame(site = row.names(Y.succ), PCO.succ)
PCO.df.succ[1:5,1:5]
### Identify positive and negative PCO axes
pos.neg.succ <- PCO.result.succ$values >= 0.0 
### Percent of variability explained
round(100 * PCO.result.succ$values/sum(PCO.result.succ$values), digits = 2)

## Failure
### Conduct PCO Analysis
PCO.result.fail <- pco(D = D.fail, varplot = TRUE)
### Standardized PCO scores to var = lambda
PCO.fail <- PCO.result.fail$vectors
colnames(PCO.fail) <- paste("PCO", 1:ncol(PCO.fail), sep = "")
### Convert the PCO matrix to a data.frame
PCO.df.fail <- data.frame(site = row.names(Y.fail), PCO.fail)
PCO.df.fail[1:5,1:5]
### Identify positive and negative PCO axes
pos.neg.fail <- PCO.result.fail$values >= 0.0 
### Percent of variability explained
round(100 * PCO.result.fail$values/sum(PCO.result.fail$values), digits = 2)


# Step 5: Conduct regression analysis on PCO axes
## Center the predictor variables around their mean
### Total
Xc.total <- X.total - rep(1, dim(X.total)[1]) %o% apply(X.total, 2, mean)
### Successful
Xc.succ <- X.succ - rep(1, dim(X.succ)[1]) %o% apply(X.succ, 2, mean)
### Failure
Xc.fail <- X.fail - rep(1, dim(X.fail)[1]) %o% apply(X.fail, 2, mean)

## Compute the regression coefficients
### NOTE: These are the beta coefficients for each predictor in each experimental unit given a particular "response" variable from PCO
### Total
pco.beta.total <- solve(t(Xc.total) %*% Xc.total) %*% t(Xc.total) %*% PCO.total
### Successful
pco.beta.succ <- solve(t(Xc.succ) %*% Xc.succ) %*% t(Xc.succ) %*% PCO.succ
### Failure
pco.beta.fail <- solve(t(Xc.fail) %*% Xc.fail) %*% t(Xc.fail) %*% PCO.fail

## Compute the fitted PCO values
### NOTE: These are the fitted values of this regression for all PCO axes and is calculated by multiplying and summing the values of X with their appropriate beta coefficient
### Total
pco.fit.total <- Xc.total %*% pco.beta.total
### Successful
pco.fit.succ <- Xc.succ %*% pco.beta.succ
### Failure
pco.fit.fail <- Xc.fail %*% pco.beta.fail

## Compute the residuals of the PCO regression
### Total
pco.res.total <- PCO.total - pco.fit.total
### Successful
pco.res.succ <- PCO.succ - pco.fit.succ
### Failure
pco.res.fail <- PCO.fail - pco.fit.fail


# Old version of the function for data loading for cleanliness purposes ####
prepData <- function(files, type, mo = "all", site = "all", drop_site = FALSE, drop_mo = FALSE, autocorrelation = FALSE){
  #files <- files_reg
  #type <- "total"
  #mo <- "all"
  #site <- "all"
  #autocorrelation <- FALSE
  #drop_site <- TRUE  # Should the function iteratively drop a single site from the initial processing
  #drop_mo <- FALSE    # Should the function iteratively drop a single month from the initial processing
  
  # Get only the Response variables
  Ys <- files[!grepl("cov_", files)]
  # Get only the predictor variables
  Xs <- files[grepl("cov_", files)]
  
  # Read in the files
  if(type == "total"){
    Y.all <- openxlsx::read.xlsx(Ys[grep(type, Ys, ignore.case = TRUE)], rowNames = TRUE)
    X.all <- openxlsx::read.xlsx(Xs[grep(type, Xs, ignore.case = TRUE)], rowNames = TRUE, rows = 2:349)
  }else if(type == "succ"){
    Y.all <- openxlsx::read.xlsx(Ys[grep(type, Ys, ignore.case = TRUE)], rowNames = TRUE)
    X.all <- openxlsx::read.xlsx(Xs[grep(type, Xs, ignore.case = TRUE)], rowNames = TRUE, rows = 2:309)
  }else if(type == "fail"){
    Y.all <- openxlsx::read.xlsx(Ys[grep(type, Ys, ignore.case = TRUE)], rowNames = TRUE)
    X.all <- openxlsx::read.xlsx(Xs[grep(type, Xs, ignore.case = TRUE)], rowNames = TRUE, rows = 2:309)
  }else{
    stop("You must choose an appropriate type. Choose one of c('total','succ','fail').")
  }
  
  # Select specific months or sites
  if(mo == "all"){
    mos <- sort(unique(Y.all$ConMonths))
    if(isTRUE(drop_mo)){
      mos2 <- lapply(1:length(mos), function(i){mos[-i]})
      names(mos2) <- paste("drop", mos, sep = "_")
    }else{
      mos2 <- list(mos)
    }
  }else{
    mos <- mo
    mos2 <- list(mo)
  }
  if(site == "all"){
    sites <- sort(unique(Y.all$WCS_ID))
    if(isTRUE(drop_site)){
      sites2 <- lapply(1:length(sites), function(i){sites[-i]})
      names(sites2) <- paste("drop", sites, sep = "_")
    }else{
      sites2 <- list(sites)
    }
  }else{
    sites <- site
    sites2 <- list(site)
  }
  
  if(isTRUE(autocorrelation)){
    cols <- c("log(Openness_m+1)", "CatWalks.Yes", "Substrate.Concrete", "Substrate.Dirt", "log(Fencing_m_+1)", 
              "PropNatural", "PropWater", "PropWoody", "Precip_total_mm", 
              "SpeedLimit", "VehicleTraffic", "log(HumanActivity+1)", "log(DomesticActivity+1)", "log(LivestockActivity+1)", "log(PropBuilding+1)", 
              paste("s_PCNM", seq(1:11), sep = ""), paste("t_PCNM", seq(1:22), sep = ""))
  }else if(isFALSE(autocorrelation)){
    cols <- c("log(Openness_m+1)", "CatWalks.Yes", "Substrate.Concrete", "Substrate.Dirt", "log(Fencing_m_+1)", 
              "PropNatural", "PropWater", "PropWoody", "Precip_total_mm", 
              "SpeedLimit", "VehicleTraffic", "log(HumanActivity+1)", "log(DomesticActivity+1)", "log(LivestockActivity+1)", "log(PropBuilding+1)", 
              "ConMonths", "UTM_X", "UTM_Y")
  }else{
    stop("autocorrelation must be logical.")
  }
  
  
  out <- pbapply::pblapply(mos2, function(m){
    M1 <- lapply(sites2, function(s){
      Y <- as.matrix(Y.all[Y.all$ConMonths %in% m & Y.all$WCS_ID %in% s, c(1:22)])
      X <- as.matrix(X.all[X.all$ConMonths %in% m & X.all$WCS_ID %in% s, colnames(X.all) %in% cols])
      D <- distance(Y = Y, measure = "BC", trans = "none", adj = 0.03)
      N <- dim(D)[1]
      list(Y = Y, X = X, D = D, N = N)
    })
    if(length(M1) == 1){
      M1 <- M1[[1]]
    }
    return(M1)
    #rm(M1)
    #rm(m)
  })
  if(length(out) == 1){
    out <- out[[1]]
  }
  
  return(out)
  rm(Ys, Xs, Y.all, X.all, Y, X, D, N, mos, mos2, sites, sites2, cols, out)
  #rm(files, type, mo, site, drop_mo, drop_site, autocorrelation)
}


# Development of the drop.one.site and drop.one.month methods for validation ####
## Load in only 1 month of data to avoid temporal autocorrelation issues
in.total.11 <- prepData(type = "total", mo = 11)
in.succ.11 <- prepData(type = "succ", mo = 11)
in.fail.11 <- prepData(type = "fail", mo = 11)

## Run the PCOs on the full in data
pco.total <- pcoOut(input = in.total$full)
pco.succ <- pcoOut(input = in.succ$full)
pco.fail <- pcoOut(input = in.fail$full)

## Run the PCOs on the subset of data
pco.total.11 <- pcoOut(input = in.total.11)
pco.succ.11 <- pcoOut(input = in.succ.11)
pco.fail.11 <- pcoOut(input = in.fail.11)

## Run the regressions on the full data
reg.total <- runRegression(input = in.total.all, PCO = pco.total)
reg.succ <- runRegression(input = in.succ.all, PCO = pco.succ)
reg.fail <- runRegression(input = in.fail.all, PCO = pco.fail)

## Run the regressions on the subsetted data
reg.total.11 <- runRegression(input = in.total.11, PCO = pco.total.11) ## THIS DOESN'T WORK

## Load in data for site drop test
in.total.site <- prepData(files = files_reg, type = "total", mo = "all", site = "all", drop_mo = FALSE, drop_site = TRUE, autocorrelation = FALSE)
in.succ.site <- prepData(files = files_reg, type = "succ", mo = "all", site = "all", drop_mo = FALSE, drop_site = TRUE, autocorrelation = FALSE)
in.fail.site <- prepData(files = files_reg, type = "fail", mo = "all", site = "all", drop_mo = FALSE, drop_site = TRUE, autocorrelation = FALSE)

## Run the PCOs on the dropped site data
pco.total.site <- lapply(in.total.site, function(x){pcoOut(input = x)})
pco.succ.site <- lapply(in.succ.site, function(x){pcoOut(input = x)})
pco.fail.site <- lapply(in.fail.site, function(x){pcoOut(input = x)})

reg.total.site <- lapply(1:length(in.total.site), function(i){runRegression(input = in.total.site[[i]], PCO = pco.total.site[[i]])})
names(reg.total.site) <- names(pco.total.site)
reg.succ.site <- lapply(1:length(in.succ.site), function(i){runRegression(input = in.succ.site[[i]], PCO = pco.succ.site[[i]])})
names(reg.succ.site) <- names(pco.succ.site)
reg.fail.site <- lapply(1:length(in.fail.site), function(i){runRegression(input = in.fail.site[[i]], PCO = pco.fail.site[[i]])})
names(reg.fail.site) <- names(pco.fail.site)

## Calculate the predicted values of the holdout samples
### Center the X's in the holdout dataset around the means of the test dataset
Xc.2 <- in.total.all$X[1:18,] - rep(1,18) %o% reg.total.site$drop_FM106_WCS1$X.mean
### Compute fitted values of the holdout dataset
YhatH <- Xc.2 %*% reg.total.site$drop_FM106_WCS1$beta

### There is relatively high correlation between the original PCO axes and the test-PCO axes for all relevant axes
cor(pco.total$vectors[-c(1:18),1:4], pco.total.site$drop_FM106_WCS1$vectors[,1:4])

### There is generally high correlation between the test-PCO axes and the fitted PCO axes of the test dataset
cor(pco.total$vectors[-c(1:18),1:4], reg.total.site$drop_FM106_WCS1$fitted[,1:4])
cor(pco.total.site$drop_FM106_WCS1$vectors[,1:4], reg.total.site$drop_FM106_WCS1$fitted[,1:4])

### There is low correlation between the true PCO axes and the fitted axes of the holdout dataset
cor(pco.total$vectors[1:18,1:4], YhatH[,1:4])

### Calculate the residuals of the holdout samples
residH <- pco.total$vectors[1:18,1:4] - YhatH[,1:4]

### Calculate the correlation between the overall fitted values and the holdout fitted values
#### Mantel correlation between the full model fitted values and predicted model fitted values with the drop-one-site approach
cor(dist(as.matrix(reg.total$fitted[1:18,1:3])), dist(as.matrix(YhatH[,1:3])))
#### Combine the predicted values with the original fitted values and compare to the raw values
cor(dist(as.matrix(pco.total$vectors[,1:3])), dist(as.matrix(rbind(YhatH[,1:3], reg.total.site$drop_FM106_WCS1$fitted[,1:3]))))

### Calculate mantel correlation between real values and predicted values of the holdout data
cor(dist(as.matrix(pco.total$vectors[1:18,1:3])), dist(as.matrix(YhatH[,1:3])))

## Calculate the mantel correlation between real values and fitted values of the training dataset
cor(dist(as.matrix(pco.total.site$drop_FM106_WCS1$vectors[,1:3])), dist(as.matrix(reg.total.site$drop_FM106_WCS1$fitted[,1:3])))


cor(pco.total.site$drop_FM106_WCS1$vectors[,1:4], reg.total.site$drop_FM106_WCS1$fitted[,1:4])

cor(pco.total.site$drop_FM106_WCS1$vectors[,1:4], reg.total.site$drop_FM106_WCS1$fitted[,1:4])^2 - cor(pco.total$vectors[1:18,1:4], YhatH[,1:4])^2

pco.total$vectors[-c(1:18),1:4]
pco.total.site$drop_FM106_WCS1$vectors[,1:4]

cor(pco.total$vectors[-c(1:18),1:4], pco.total.site$drop_FM106_WCS1$vectors[,1:4])
cor(pco.total$vectors[-c(1:18),1:4], reg.total.site$drop_FM106_WCS1$fitted[,1:4])

ggplot(data = pred.total$pred$drop_FM106_WCS1$data.pred, mapping = aes(x = PCO1, y = PCO2, label = 1:18)) + 
  geom_text() + 
  geom_text(data = pred.total$pred$drop_FM106_WCS1$Yhat.pred, color = "blue") + 
  geom_text(data = pred.total$pred$drop_FM106_WCS1$orig.fit, color = "red") + 
  theme_bw()


ggplot(data = pco.total$vectors[1:18,], mapping = aes(x = PCO1, y = PCO2, label = 1:18)) + 
  geom_text() + 
  geom_text(data = YhatH, color = "blue") + 
  geom_text(data = reg.total$fitted[1:18,], color = "red")

apply(pco.total$vectors[1:18,1:4] - (reg.total$Xc[1:18,] %*% reg.total.site$drop_FM106_WCS1$beta)[,1:4], 2, function(x){sum(x^2)})


# Preliminary prediction plots ####
ggplot(data = plot.site, mapping = aes(x = PCO1, y = PCO2, color = drop, shape = drop)) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  geom_point() + 
  scale_color_manual("Site", values = c("red", "orange", "yellow", "green", "blue", "purple", "violet", "black", 
                                        "red", "orange", "yellow", "green", "blue",
                                        "red", "orange", "yellow", "green", "blue")) + 
  scale_shape_manual("Site", values = c(rep(15, 8), rep(16, 5), rep(17, 5))) + 
  coord_fixed() + 
  facet_grid(rows = vars(ds), cols = vars(response)) + 
  theme_bw()


ggplot(data = plot.time, mapping = aes(x = PCO1, y = PCO2, color = drop, shape = drop)) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  geom_point() + 
  scale_color_manual("Month", values = c("red", "yellow", "blue", "purple", 
                                         "red", "yellow", "blue", "purple", 
                                         "red", "yellow", "blue", "purple")) + 
  scale_shape_manual("Month", values = rep(c(15:18), each = 3)) + 
  coord_fixed() + 
  facet_grid(rows = vars(ds), cols = vars(response)) + 
  theme_bw()

ggsave(filename = paste("Preliminary_", "site_prediction_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = ""), device = "tiff", width = 6.5, height = 11, dpi = 600, compression = "lzw")

