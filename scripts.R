suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(mosaic))
suppressPackageStartupMessages(library(rlang))

## Given a Data Frame of Categories and Counts, Create the Extended data version
countsToCases <- function(x, countcol = "Freq") {
  idx <- rep.int(seq_len(nrow(x)), x[[countcol]])
  x[[countcol]] <- NULL
  x[idx, ]
}

## Given an aggregate table in csv, convert the table into long form
table_to_df <- function(data = table){
  temp <- data |>
    pivot_longer(cols = 2:5, values_to = "Freq") |>
    separate(
      col = name,
      into = c("sex", "status")
    ) |> mutate(
      group = factor(group),
      sex = factor(sex),
      status = factor(toupper(status))
    )
  return(countsToCases(temp))
}

#### Checkerboard Copulas Code

## Joint PMF of Data
getPMF_joint <- function(data = df, type = "table"){
  ## Usual table version
  if(type == "table"){
    return(prop.table(table(data)))
  }
  ## Get frequency version of table for Joint CDF calculation
  else if(type == "dataframe"){
    return(data |>
             group_by_all() |>
             dplyr::count() |>
             ungroup() |>
             mutate(prop = n / sum(n)) |>
             select(-n)
           )
  }
  else{
    print(paste("ERROR INVALIDE TYPE:", type))
  }
}

findPMF_joint <- function(range_vec = c(0.5,0.5), jointPMF = pmf, 
                          all_cdfs = cdfs){
  #jointPMF <- getPMF_joint(data = data)
  
  #all_cdfs <- lapply(1:ncol(data), getCDF_margin, data = data)
  
  ## Find category they belong in
  idx_vec <- mapply(FUN = function(cdf, vec){
    findInterval(vec, cdf, left.open = TRUE, rightmost.closed = TRUE)
  }, cdf = all_cdfs, vec = range_vec)
  
  return(jointPMF[idx_vec[1], idx_vec[2]])
}

## Returns marginal PMF of response variable
getPMF_margin <- function(var_index = 1, data = df){
  if ("data.frame" %in% class(data)) {
    return(data |>
             count(data[, var_index]) |>
             mutate(prop = n / sum(n)) |>
             select(prop) |>
             unlist())
  } else {
    if (var_index == 1) {
      return(rowSums(data))
    } else {
      return(colSums(data))
    }
  }
}

## Returns CDF
getCDF_margin <- function(var_index = 1, data = df){
  return(
    ## append 0 at the beginning
    append(cumsum(getPMF_margin(var_index = var_index, data = data)), value = 0, 
           after = 0)
  )
}

getCDF_joint <- function(range_vec = c(0.5, 0.5), all_cdfs = cdfs, data = data){
  ## Marginal CDFs
  #all_cdfs <- lapply(1:ncol(data), getCDF_margin, data = data)
  
  ## Find category they belong in
  idx_vec <- mapply(FUN = function(cdf, vec){
    findInterval(vec, cdf, left.open = TRUE, rightmost.closed = TRUE)
  }, cdf = all_cdfs, vec = range_vec)
  ## Find combinations of all categories
  combinations <- expand.grid(lapply(data, levels))
  
  ## ensure group_by is successful
  colnames(combinations) <- colnames(data)
  
  ## Calculate the joint PMF
  jointPMF <- getPMF_joint(data = data, type = "dataframe")
  combined <- suppressMessages(
    left_join(combinations, jointPMF, by = colnames(combinations))
  )
  ## If any values are NA, make it 0 for calculations
  combined[is.na(combined)] <- 0
  combined <- combined |>
    mutate(across(.cols = everything(), as.numeric))
  
  ## Filter for desired values and then just add it all
  for(i in 1:length(idx_vec)){
    combined <- combined[combined[, i] <= idx_vec[i], ]
  }
  return(sum(combined$prop))
}

getCDF_joint_idx <- function(idx_vec = c(1, 2), all_cdfs = cdfs, data = data){
  ## Marginal CDFs
  #all_cdfs <- lapply(1:ncol(data), getCDF_margin, data = data)
  
  ## Find combinations of all categories
  combinations <- expand.grid(lapply(data, levels))
  
  ## ensure group_by is successful
  colnames(combinations) <- colnames(data)
  
  ## Calculate the joint PMF
  jointPMF <- getPMF_joint(data = data, type = "dataframe")
  combined <- suppressMessages(
    left_join(combinations, jointPMF, by = colnames(combinations))
  )
  ## If any values are NA, make it 0 for calculations
  combined[is.na(combined)] <- 0
  combined <- combined |>
    mutate(across(.cols = everything(), as.numeric))
  
  ## Filter for desired values and then just add it all
  for(i in 1:length(idx_vec)){
    combined <- combined[combined[, i] <= idx_vec[i], ]
  }
  return(sum(combined$prop))
}

getCheckerboardCopula <- function(range_vec = c(0.5,0.5),
                                         data = df){
  ## Get all CDFs of each variable
  all_cdfs <- lapply(1:ncol(data), getCDF_margin, data = data)
  
  ## infimums for each variable
  inf_values <- mapply(FUN = function(cdf, vec){max(cdf[cdf<=vec])},
                       cdf = all_cdfs, vec = range_vec)
  ## supremums
  sup_values <- mapply(FUN = function(cdf, vec){max(cdf[cdf>=vec])},
                       cdf = all_cdfs, vec = range_vec)
  lambdas <- mapply(FUN = function(sup, inf, vec){
                      ifelse(inf < sup, (vec - inf) / (sup - inf), 1)
                    }, sup = sup_values, inf = inf_values, vec = range_vec)
  
  ## Now find the categories that match with the range_vec inputted
  sup_idx <- mapply(FUN = function(cdf, vec){
    findInterval(vec, cdf, left.open = TRUE, rightmost.closed = TRUE)
  }, cdf = all_cdfs, vec = sup_values)
  
  inf_idx <- mapply(FUN = function(cdf, vec){
    findInterval(vec, cdf, left.open = TRUE, rightmost.closed = TRUE)
  }, cdf = all_cdfs, vec = inf_values)
  
  ## Create all combinations of values up to number of variables in the data
  ## To line up with calculations
  subsets <- sets::set_power(1:ncol(data))
  
  ## calculate the checkerboard Copula Density
  sum <- 0
  for(i in subsets){
    prod <- 1
    idx_vector <- c()
    for(j in c(1:ncol(data))){
      if(j %in% i) {
        prod <- prod * lambdas[j]
        idx_vector <- append(idx_vector, sup_idx[j])
      } 
      else{
        prod <- prod * (1- lambdas[j])
        idx_vector <- append(idx_vector, inf_idx[j])
      }
    }
    sum <- sum + (prod * getCDF_joint_idx(idx_vec = idx_vector,
                                          all_cdfs = all_cdfs,
                                          data = data))
  }
  return(sum)
}

getCCDensityFunction <- function(data = df){
  ## unique row combinations
  data_int <- unique(data) |>
    mutate(across(.cols = everything(), as.numeric))
  
  ## get all marginal PMFS
  all_pmfs <- lapply(1:ncol(data), getPMF_margin, data = data)
  
  ##Code to calculate the the density values for each combination seen 
  cc_prod <- data.frame()
  for(i in 1:nrow(data_int)){
    c <- c()
    for(j in 1:length(all_pmfs)){
      c <- append(c, all_pmfs[[j]][pull(data_int[i, j])])
    }
    cc_prod <- rbind(cc_prod, c)
  }
  cc_prod <- cc_prod |>
    mutate(prod = Reduce('*', cc_prod))
  
  ## Return the actual value
  return_data <- getPMF_joint(data = data, type = "dataframe") |>
    mutate(cc_value = prop / cc_prod$prod)
  
  return(return_data)
}

getCCDensity <- function(range_vec = c(0.5, 0.5), data = df, density = density){
  ## Get all marginal CDFs
  all_cdfs <- lapply(1:ncol(data), getCDF_margin, data = data)
  
  ## find the idx which categories belong in
  values <- mapply(FUN = function(x, y) {
    findInterval(y, x, left.open = TRUE, rightmost.closed = TRUE)
  }, x = all_cdfs, y = range_vec)
  
  # get density
  density2 <- density
  
  ## return density value
  for (i in 1:length(range_vec)) {
    density2 <- density2[density2[, i] == values[i], ]
  }
  if(nrow(density2) == 0){
    return(0)
  }
  else{
    return(density2$cc_value)
  }
}

## Takes in the CDF, returns a vector of scores
getScores <- function(data = Xcdf){
  ## Ridit Calculation
  zoo::rollmean(data, 2)
}

## Very similar to getPMF_margins() function, except now we are grouping by
## explanatory variables at the end
## Returns Conditional PMF of each response category based off of explanatory vars
getPMF_conditional <- function(var_index = 1, data = df) {
  grouping_cols <- c(1:ncol(data))[-var_index]
  cond <- data |>
    group_by_all() |>
    count() |>
    ungroup() |>
    mutate(prop = n / sum(n)) |>
    ## These lines of code below differentiates this function with getPMF_margins() 
    group_by(across(all_of(grouping_cols))) |>
    ## calculate conditional proportion of each response
    mutate(prop_cond = n / sum(n)) |>
    ungroup()
  return(cond)
}


## Returns the BECCR data frame given data and response var
getRegression <- function(var_index = 1, data = df) {
  cdf <- getCDF_margin(var_index = var_index, data = data)
  scores <- getScores(data = cdf)
  
  # getting conditional pmfs
  conditional_pmf <- getPMF_conditional(
    var_index = var_index,
    data = data
  ) |>
    ## Using as.numeric because it makes calculations easier
    mutate(across(.cols = everything(), as.numeric))
  
  ## Explanatory Vars
  grouping_cols <- c(1:ncol(data))[-var_index]
  
  ## Given the scores vector, find what the Checkerboard Copula Score is for each 
  ## response based on what category it is. Then combine it with the conditional_pmf
  ## data frame to get a final data frame
  
  scores_vec <- lapply(conditional_pmf[, var_index],
                       FUN = function(x) scores[x])
  
  combined <- cbind(conditional_pmf, scores_vec)
  
  ## Rename the last column (which is the scores_vec) to score
  colnames(combined)[ncol(combined)] <- "score"
  
  ## Calculate what the regression should be, which is the 
  ## sum of conditional proportion * score
  ## across all explanatory variables (hence the group_by)
  regression_data <- combined |> 
    mutate(regression = prop_cond * score) |>
    group_by(across(all_of(grouping_cols))) |>
    summarize(reg = sum(regression), .groups = "drop")
  return(regression_data)
}

## Takes in regression and cdf of the data, as well as the response categories
## that are in order that they are present in the regression df columns
## and returns the predicted category of the response
predictCategory <- function(idx_vec = vecs, regression = BECCR, cdf = Xcdf) {
  ##
  for (i in 1:length(idx_vec)) {
    # print(i)
    regression <- regression[regression[, i] == idx_vec[i], ]
  }
  
  pred_val <- findInterval(regression$reg, cdf, left.open = TRUE, rightmost.closed = TRUE)
  if (!length(pred_val)) {
    #return("no prediction")
    return(NA)
  }
  return(pred_val)
}

## Returns BECCR but takes in the CDF as an input as well
## Will be used in the Bootstrapping
getBootstrapRegression <- function(var_index = 1, cdf = Xcdf, data = df){
  scores <- getScores(data = cdf)
  # getting conditional pmfs
  conditional_pmf <- getPMF_conditional(
    var_index = var_index,
    data = data
  ) |>
    mutate(across(.cols = everything(), as.numeric))
  
  grouping_cols <- c(1:ncol(data))[-var_index]
  
  scores_vec <- lapply(conditional_pmf[, var_index],
                       FUN = function(x) scores[x]
  )
  combined <- cbind(conditional_pmf, scores_vec)
  colnames(combined)[ncol(combined)] <- "score"
  regression_data <- combined |> 
    mutate(regression = prop_cond * score) |>
    group_by(across(all_of(grouping_cols))) |>
    summarize(reg = sum(regression), .groups = "drop")
  return(regression_data)
}

runBootstrapRegression <- function(var_index = var_index, data = df,
                                   seed = 713, runs = 100){
  set.seed(seed)
  if(typeof(data[,var_index]) != "list"){
    print("Please enter a tibble data frame")
    return(NA)
  }
  ## Making the factor columns all numeric, as it will cause errors otherwise
  temp_df <- data |>
    mutate(across(.cols = everything(), as.numeric))
  
  ## Explainatory variables
  grouping_cols <- c(1:ncol(data))[-var_index]
  
  ## Generic way to find all combinations of the categories of predictors
  list_cols <- list()
  for(col in grouping_cols){
    list_cols <- append(list_cols, unique(temp_df[,col])) 
  }
  prediction_data <- expand.grid(list_cols)
  
  ## Iterating for each run
  for(n in 1:runs){
    #print(n)
    ################# 
    ##IMPORTANT: THIS IS HOW I AM IMPLEMENTING KEEPING TRACK OF THE RESAMPLED CDF
    ## Find the CDF of the actual data set first so the length is set
    resample_cdf <- getCDF_margin(var_index = var_index, data = temp_df)
    
    ## Resample with replacement of the data
    resample_data <- resample(temp_df, nrow(temp_df), replace = TRUE)[,1:length(temp_df)]
    
    ## If resampling somehow resamples s.t. a category is missing,
    ## will cause errors so first will see what all the observed categories are
    ## Sorting is important as we will see later
    observed_categories <- sort(unique(resample_data[,var_index])|> pull())
    
    ## temporary cdf of the observed categories (may potentially not line up with
    ## actual CDF)
    temp_cdf <- getCDF_margin(var_index = var_index, data = resample_data)
    
    ## temporary index, as this will always be the cdf value of the "first" category observed
    temp_i = 2
    
    ## Recall that observed_categories are all the categories of the response variable
    ## and that each category is a number starting from 1 to k
    ## for all responses detected, their index in the actual "cdf" vector is 1 shifted
    ## to the right because we start from 0.
    for(i in observed_categories){
      
      ## Shifting idx to the right 1 to match up with the actual CDF placement
      idx = i+1
      ## update the resample cdf for the actual values calculated by the temporary cdf
      resample_cdf[idx] = as.numeric(temp_cdf[temp_i])
      
      ## now keep on shifting the temp_index because next category
      temp_i <- temp_i + 1
      
      ## it is important to note that this works because it is SORTED from least to greatest
    }
    
    ## If a category of the response doesn't show up, will affect the generated cdf
    ## So will manually adjust
    
    ## vector containing all the categories that don't show up in observed_categories
    to_fix <- setdiff(1:(length(resample_cdf)-1), observed_categories)
    
    for(i in to_fix){
      ## shifting idx to the right for cdf again
      idx = i + 1
      
      ## these are the idx of the categories to the "right" and "left" of our missing category
      right_end_point = findInterval(i, observed_categories) + 1
      left_end_point = findInterval(i, observed_categories)
      
      ## this is the index of the category to the right of i in the resampled_cdf
      right_idx <- observed_categories[right_end_point] + 1
      
      ## idx of the category to the left of i in the resampled_cdf
      ## if the left_end_point is 0, then make left_idx 1 because start of the first vec
      if(left_end_point == 0){
        left_idx = 1
      } else{
        left_idx = observed_categories[left_end_point] + 1
      }
      
      ## if right_idx DNE due to observed_categories not being big enough, then 
      ## right_idx is simply the right most value
      if(is.na(right_idx)){
        right_idx = length(resample_cdf)
      }
      
      ## If a category is missing, basically draw a straight line from last two
      ## known categories and "extrapolate" the cdf at that point
      
      ## difference between right_idx and left_idx in case more than one response is missing
      ## inbetween values
      delta = right_idx - left_idx
      
      ## so manually adjusting the resample_cdf based off our values
      resample_cdf[idx] <- (idx-left_idx)*(resample_cdf[right_idx] - resample_cdf[left_idx])/delta +
        resample_cdf[left_idx]
      
    }
    ## End of my method to manually adjust CDF
    #################
    
    ## using the manually adjusted cdf as well as resampled data, calculate regression
    resample_regression <- getBootstrapRegression(var_index = var_index, 
                                                  cdf = resample_cdf, 
                                                  data = resample_data)
    
    ## Generic method to find to perform predictions given n-1 predictors
    ## Creates a list where each "index" of the list is a vector of the categories for that desired
    ## predictor
    listz = list()
    for(col in 1:(length(resample_regression) - 1)){
      listz <- append(listz, resample_regression[,col])
    }
    
    ## This essentially "pairs" or groups the elements in the same index of each vector 
    ## in the list
    ## i.e. for every vector in listz, the first element is grouped together in a vector, 
    ## then the every second, etc.
    
    vec_list <- list()
    for(row in 1:nrow(resample_regression)){
      vec_list[[row]] <- sapply(listz, "[[", row)
    }
    
    
    ## as we now have a list of vectors (which are predictions)m just lapply our 
    ## predictCategory() function
    estimates <- lapply(vec_list, predictCategory,
                        regression = resample_regression,
                        cdf = resample_cdf)
    
    ## group by all unique combinations of explanatory variables that show up
    resample_prediction <- resample_data |>
      group_by(across(all_of(grouping_cols))) |>
      summarise(.groups = "drop")
    
    ## creating new column name known as index_1 to index_n
    name = paste0("index_", n)
    
    ## create new column in our prediction
    resample_prediction[[name]] = estimates |> unlist()
    
    ## suppressing because I know column names will all match up
    ## left join because prediction_data has all possible explanatory combinations
    ## and resample_prediction may be missing some
    prediction_data <- suppressMessages(left_join(prediction_data, resample_prediction))
  }
  return(prediction_data)
}

getBubblePlotData <- function(regression_data = predictions, data = df, 
                              var_index = var_index){
  
  ## Explanatory Vars
  grouping_cols <- c(1:ncol(data))[-var_index]
  
  ## want to generically find the ending column based off parameters
  end_num <- length(regression_data) - (length(data) - 1)
  end_col <- paste0("index_",end_num)
  
  ## so can refer to the response variable
  response_var <- colnames(data)[var_index]
  
  ## Pivot all columns from index_1 to index_1000
  long_data <- regression_data |> 
    pivot_longer(cols = "index_1":!!(end_col), 
                 names_to = "sample",
                 values_to = response_var) |>
    select(-sample)
  
  ## Find the predicted counts given resampling
  counts <- long_data |>
    dplyr::count(across(.cols = everything())) 
  counts[response_var] = unlist(counts[response_var]) |> as.numeric()
  
  ## Generic way to to find combinations of all explanatory categories
  list_cols <- list()
  for(col in 1:length(data)){
    list_cols <- append(list_cols, unique(data[,col])) 
  }
  combinations <- expand.grid(list_cols) |>
    mutate(across(.cols = everything(), as.numeric))
  
  ## Combine by left joining, and replacing any NA values with 0
  full_data <- suppressMessages(left_join(combinations,counts)) |>
    replace_na(list(n = 0)) 
  return(full_data)
}


createBubblePlotUngroup <- function(var_index = 1, 
                                    regression_data = predictions, data = df){
  full_data <- getBubblePlotData(var_index = var_index, regression_data = regression_data,
                                 data = data) |>
    mutate(Y_name = factor(paste0("y", Y), 
                           levels = c("y5", "y4", "y3", "y2", "y1")),
           X_name = paste0("x", X))
  
  if(var_index == 1 || var_index == 2){
    p <- ggplot(full_data, aes(x = X_name, y = Y_name)) +
      ## Drawwing the "bubbles" 
      geom_point(aes(
        size = ifelse(n == 0, NA, n * 100), color = X_name),
        shape = 21, colour = "black",
        fill = "white", stroke = .5
      ) +
      ## making the scaling of the "bubbles"
      scale_size_continuous(range = c(1, 20)) +
      ## drawing the "dots"
      geom_point(data = full_data, aes(
        x = X_name,
        y = Y_name,
        size =
          ifelse(n == 0 | n <= 500,
                 NA, .1
          )
      )) +
      geom_text(aes(label = ifelse(n == 0, "", n)),
                vjust = -1.2, size = 2
      ) +
      theme(
        legend.position = "none",
        text = element_text(size = 13),
        axis.text.y = element_text(
          vjust = 0.5,
          hjust = 1, size = 8
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
      ) +
      labs(
        x = "X",
        y = "Y",
        #title = "Aggregate Bubble Plot"
      )
  } else {
    p <- ggplot(full_data, aes(x = Y_name, y = X_name)) +
      ## Drawwing the "bubbles"
      geom_point(aes(
        size = ifelse(n == 0, NA, n * 100), color = Y_name),
        shape = 21, colour = "black",
        fill = "white", stroke = .5
      ) +
      ## making the scaling of the "bubbles"
      scale_size_continuous(range = c(1, 20)) +
      ## drawing the "dots"
      geom_point(data = full_data, aes(
        x = Y_name,
        y = X_name,
        size =
          ifelse(n == 0 | n <= 500,
                 NA, .1
          )
      )) +
      geom_text(aes(label = ifelse(n == 0, "", n)),
                vjust = -1.2, size = 2
      ) +
      theme(
        legend.position = "none",
        text = element_text(size = 13),
        axis.text.y = element_text(
          vjust = 0.5,
          hjust = 1, size = 8
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
      ) +
      labs(
        x = "X",
        y = "Y",
        #title = "Aggregate Bubble Plot"
      )
  }
  
  return(p)
}

createBubblePlot <- function(data = df, x = x_combo, y = y_name, B = 1000){
    ggplot(data = data, aes(x = !! sym(x), y = !! sym(y))) +
    ## Drawing the "bubbles"
    geom_point(aes(
      size = ifelse(n == 0, NA, n * 100), color = sex),
      shape = 21, colour = "black",
      fill = "white", stroke = .5
    ) +
    ## making the scaling of the "bubbles"
    scale_size_continuous(range = c(1, 20)) +
    ## drawing the "dots"
    geom_point(data = data, mapping = aes(x = !! sym(x), y = !!sym(y),
      size =
        ifelse(n == 0 | n <= B/2,
               NA, .1
        )
    )) +
    geom_text(aes(label = ifelse(n == 0, "", n)),
              vjust = -1.2, size = 2
    ) +
    theme(
      legend.position = "none",
      text = element_text(size = 13),
      axis.text.y = element_text(
        vjust = 0.5,
        hjust = 1, size = 8
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    )
}

## For specific example
createBubblePlot2 <- function(var_index = 1, regression_data = predictions, data = df){
  full_data <- getBubblePlotData(var_index = var_index, regression_data = regression_data,
                                 data = data) |>
    mutate(Y_name = factor(paste0("y", Y), 
                           levels = c("y5", "y4", "y3", "y2", "y1")),
           X_name = paste0("x", X),
           G_name = paste0("g", G))
  
  if(var_index == 1){
    full_data$combination_xg <- paste0(
      "(",
      full_data$G_name,
      ",",
      full_data$X_name,
      ")"
    )
    p <- ggplot(full_data, aes(x = combination_xg, y = Y_name)) +
      ## Drawwing the "bubbles"
      geom_point(aes(
        size = ifelse(n == 0, NA, n * 100), color = G),
        shape = 21, colour = "black",
        fill = "white", stroke = .5
      ) +
      ## making the scaling of the "bubbles"
      scale_size_continuous(range = c(1, 20)) +
      ## drawing the "dots"
      geom_point(data = full_data, aes(
        x = combination_xg,
        y = Y_name,
        size =
          ifelse(n == 0 | n <= 500,
                 NA, .1
          )
      )) +
      geom_text(aes(label = ifelse(n == 0, "", n)),
                vjust = -1.2, size = 2
      ) +
      theme(
        legend.position = "none",
        text = element_text(size = 13),
        axis.text.y = element_text(
          vjust = 0.5,
          hjust = 1, size = 8
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
      ) +
      labs(
        x = "Combination of Group and X",
        y = "Y",
        #title = "Bubble Plot"
      )
  } else{
    
    full_data$combination_yg <- paste0(
      "(",
      full_data$G_name,
      ",",
      full_data$Y_name,
      ")"
    )
    p <- ggplot(full_data, aes(x = combination_yg, y = X_name)) +
      ## Drawwing the "bubbles"
      geom_point(aes(
        size = ifelse(n == 0, NA, n * 100), color = G),
        shape = 21, colour = "black",
        fill = "white", stroke = .5
      ) +
      ## making the scaling of the "bubbles"
      scale_size_continuous(range = c(1, 20)) +
      ## drawing the "dots"
      geom_point(data = full_data, aes(
        x = combination_yg,
        y = X_name,
        size =
          ifelse(n == 0 | n <= 500,
                 NA, .1
          )
      )) +
      geom_text(aes(label = ifelse(n == 0, "", n)),
                vjust = -1.2, size = 2
      ) +
      theme(
        legend.position = "none",
        text = element_text(size = 13),
        axis.text.y = element_text(
          vjust = 0.5,
          hjust = 1, size = 8
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
      ) +
      labs(
        x = "Combination of Group and Y",
        y = "X",
        #title = "Bubble Plot"
      )
  }
  
  return(p)
}

createMultiBKPlot <- function(data = df){
  data_bk <- data |>
    mutate(new_Y1 = factor(ifelse(Y == "y1", "y1", "y_other")),
           new_Y2 = factor(ifelse(Y == "y2", "y2", "y_other")),
           new_Y3 = factor(ifelse(Y == "y3", "y3", "y_other")),
           new_Y4 = factor(ifelse(Y == "y4", "y4", "y_other")),
           new_Y5 = factor(ifelse(Y == "y5", "y5", "y_other")),
    )
  
  data_bk_temp <- data_bk |>
    select(new_Y1, X, G)
  
  p1 <- createBKPlot(data = data_bk_temp, var_index = 1, by_index = 2)
  
  data_bk_temp <- data_bk |>
    select(new_Y2, X, G) |>
    arrange(desc(new_Y2))
  
  p2 <- createBKPlot(data = data_bk_temp, var_index = 1, by_index = 2)
  
  data_bk_temp <- data_bk |>
    select(new_Y3, X, G) |>
    arrange(desc(new_Y3))
  
  p3 <- createBKPlot(data = data_bk_temp, var_index = 1, by_index = 2)
  
  data_bk_temp <- data_bk |>
    select(new_Y4, X, G) |>
    arrange(desc(new_Y4))
  
  p4 <- createBKPlot(data = data_bk_temp, var_index = 1, by_index = 2)
  
  data_bk_temp <- data_bk |>
    select(new_Y5, X, G) |>
    arrange(desc(new_Y5))
  
  p5 <- createBKPlot(data = data_bk_temp, var_index = 1, by_index = 2)
  return(invisible(list(p1,p2,p3,p4,p5)))
}

