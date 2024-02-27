require("fixest")
require("plm")
require("tidyverse")

get_weights <- function(df, j, time_index, lim){
  
  df[,paste0("group_h", j)] <- NA
  df[,paste0("group_h", j)] <- ifelse(lim, df[,time_index], df[,paste0("group_h", j)])
  
  frmla <- as.formula(paste0("exp_var ~ 1 | ", time_index))
  tmp <- feols(frmla, data = df[lim,])
  
  df[,paste0("num_weights_", j)] <- NA
  df[lim,paste0("num_weights_", j)] <- tmp$residuals
  df[is.na(df$exp_var) | df$exp_var == 0,paste0("num_weights_", j)] <- NA
  den_weights = sum(df[,paste0("num_weights_", j)], na.rm = T)
  df[,paste0("weight_", j)] <- df[,paste0("num_weights_", j)] / den_weights
  suppressWarnings(df[,paste0("gweight_", j)] <- ave(df[,paste0("weight_", j)], df[,paste0("group_h", j)], FUN = function(x) max(x, na.rm = TRUE)))
  
  lim <- is.na(df[,paste0("weight_", j)])
  df[lim,paste0("weight_", j)] <- df[lim,paste0("gweight_", j)]
  
  tmp <- 1 / df[,paste0("weight_", j)]
  df <- df[,colnames(df)[!grepl("group_|weights_|weight_[^0]|.weight_", colnames(df))]]
  return(tmp)
}

lpdid <- function(df, window = c(NA, NA), y,
                  unit_index, time_index,
                  exp_var = "",
                  controls = NULL, outcome_lags = 0,
                  reweight = FALSE,
                  nonabsorbing = TRUE, nonabsorbing_lag
                  ){
  
  pre_window <- -1*window[1]; post_window <- window[2]
  
  #create difference of exp_var
  exp_var_sym <- rlang::sym(exp_var)  # Convert the string to a symbol
  df<-df %>%
    group_by(countrycode) %>%
    mutate(exp_var := {{exp_var_sym}} - dplyr::lag({{exp_var_sym}}),
           exp_var=case_when(is.na(exp_var)~0,
                             TRUE~exp_var))%>%ungroup()

  # Convert df to pdata.frame
  if(!inherits(df, "pdata.frame")) df <- pdata.frame(df, index=c(unit_index,time_index), drop.index=FALSE, row.names=FALSE)
  df[,unit_index] <- as.character(df[,unit_index]); df[,time_index] <- as.numeric(df[,time_index])

  # create "rel_time" variable for nonabsorbing
  rel_time <- "rel_time"
  df[,rel_time] <- NA
  for(i in unique(df[,unit_index])){
    
    id_logic <- which(df[,unit_index] == i)
    #change from ==1 to >0
    valz <- df[id_logic[df[id_logic, exp_var] > 0], time_index]
    if(!is.na(valz[1])){
      
      mat <- matrix(NA, ncol = length(valz), nrow = length(id_logic))
      
      for(j in valz){
        
        mat[,which(j == valz)] <- df[id_logic, time_index] - j
      }
      
      df[id_logic, rel_time] <- apply(mat, 1, function(x) x[which(abs(x) == min(abs(x)))[1]])
    } else {
      
      df[id_logic, rel_time] <- -1000
    }
  }


df$treat <- ifelse(df[,rel_time] >= 0, 1, 0)


  # Calculate lags of the outcome
  if(outcome_lags>0){

    for(outcome_lag in 1:outcome_lags){

      df[,paste0("y_diff_lag", outcome_lag)] <- lag(df[,y] - lag(df[,y], 1), outcome_lag)
    }
    controls <- c(controls, colnames(df)[grepl("y_diff_lag.", colnames(df))])
  }

  lpdid_betaz <- rep(0, length(-pre_window:post_window))
  lpdid_sez <- rep(0, length(-pre_window:post_window))
  lpdid_nz <- rep(0, length(-pre_window:post_window))

  loop_bound <- max(post_window, pre_window)
for(j in 0:loop_bound){

    if(!reweight) df$reweight_0 <- df$reweight_use <- 1

    # Post
    if(j <= post_window){

      df$Dy <- lead(df[,y], j) - lag(df[,y], 1)

      # Create Formula
      if(!is.null(controls)) controls <- paste0(" + ", paste(controls, collapse = " + "))
      frmla <- as.formula(paste0("Dy ~ exp_var", controls, " | ", time_index))
      df$cluster_var <- df[,unit_index]

      # Create "Limit"

        ## Non-Absorbing Limit
        lim_ctrl <- TRUE; lim_treat <- TRUE
        for(i in -nonabsorbing_lag:j){

          lim_ctrl <- lim_ctrl & if(i >= 0) lead(df$exp_var, i) == 0 else lag(df$exp_var, -i) == 0
          lim_treat <- lim_treat & if(i >= 0) lead(df$treat, i) == 1 else lag(df$treat, -i) == 0
        }
        lim <- lim_ctrl | lim_treat

      lim <- !is.na(lim) & lim

      # Calculate weights for j
      if(reweight){

        df$reweight_use <- get_weights(df = df, j = j, time_index = time_index, lim = lim)
        if(j == 0) df$reweight_0 <- df$reweight_use
      }

      lim <- lim & df$reweight_use > 0


        # Estimate and Save
        tmp <- suppressMessages(feols(frmla, data = df[lim,], cluster = ~cluster_var, weights = ~reweight_use))
        lpdid_betaz[match(j, -pre_window:post_window)] <- tmp$coeftable[1,1]
        lpdid_sez[match(j, -pre_window:post_window)] <- tmp$coeftable[1,2]
        lpdid_nz[match(j, -pre_window:post_window)] <- nobs(tmp)

    }

    # Pre
    if(j>1 & j<=pre_window){

      df$Dy <- lag(df[,y], j) - lag(df[,y], 1)

      if(!is.null(controls)) controls <- paste0(" + ", paste(controls, collapse = " + "))
      frmla <- as.formula(paste0("Dy ~ exp_var", controls, " | ", time_index))
      df$cluster_var <- df[,unit_index]

      # Create Limits
        lim_ctrl <- TRUE; lim_treat <- TRUE
        for(i in -nonabsorbing_lag:j){

          lim_ctrl <- lim_ctrl & if(i >= 0) lead(df$exp_var, i) == 0 else lag(df$exp_var, -i) == 0
          lim_treat <- lim_treat & if(i >= 0) lead(df$treat, i) == 1 else lag(df$treat, -i) == 0
        }
        lim <- lim_ctrl | lim_treat


      lim <- !is.na(lim) & lim & df$reweight_0 > 0


    # Estimate
        suppressMessages(tmp <- feols(frmla, data = df[lim,],vcov = DK~rel_time, weights = ~reweight_0, panel.id = ~countrycode+rel_time))
        lpdid_betaz[match(-j, -pre_window:post_window)] <- tmp$coeftable[1,1]
        lpdid_sez[match(-j, -pre_window:post_window)] <- tmp$coeftable[1,2]
        lpdid_nz[match(-j, -pre_window:post_window)] <- nobs(tmp)

    }
   }

  coeftable <- data.frame(Estimate = lpdid_betaz,
                          "Std. Error" = lpdid_sez,
                          "t value" = lpdid_betaz/lpdid_sez,
                          "Pr(>|t|)" = NA,
                          check.names = FALSE)
  coeftable[,4] <- pnorm(coeftable$`t value`, lower.tail = F)
  return(list(coeftable = coeftable[!is.na(coeftable$Estimate),],
              # df = df,
              window = c(-pre_window:post_window)[!is.na(coeftable$Estimate)],
              nobs = data.frame(window = c(-pre_window:post_window)[!is.na(coeftable$Estimate)],
                                nobs = lpdid_nz)))
}
