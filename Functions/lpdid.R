require("fixest")
require("plm")
require("tidyverse")

lpdid <- function(df, window = c(-5,7), y,
                  unit_index, time_index,
                  exp_var = "",
                  controls = NULL, outcome_lags = 0,
                  reweight = FALSE,
                  nonabsorbing = TRUE, nonabsorbing_lag
                  ){
  
  pre_window <- -1*window[1]; post_window <- window[2]
  
  #create difference of exp_var
  exp_var_sym <- rlang::sym(exp_var) 
  df<-df %>%
    group_by(countrycode) %>%
    mutate(exp_var := {{exp_var_sym}} - dplyr::lag({{exp_var_sym}})#Delta FR
           )%>%
    mutate(across(matches(controls),~.-dplyr::lag(.)))%>% #Delta Controls
    ungroup()

  # Convert df to pdata.frame
  if(!inherits(df, "pdata.frame")) df <- pdata.frame(df, index=c(unit_index,time_index), drop.index=FALSE, row.names=FALSE)


  # Lagged dep var
  if(outcome_lags>0){

    for(outcome_lag in 1:outcome_lags){

      df[,paste0("y_diff_lag", outcome_lag)] <- lag(df[,y] - lag(df[,y], 1), outcome_lag)
    }
    controls <- c(controls, colnames(df)[grepl("y_diff_lag.", colnames(df))])
  }

#saving  
  lpdid_beta <- rep(0, length(-pre_window:post_window))
  lpdid_se <- rep(0, length(-pre_window:post_window))
  lpdid_n <- rep(0, length(-pre_window:post_window))

  loop_bound <- max(post_window, pre_window)
for(j in 0:loop_bound){

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
          lim_treat <- lim_treat & if (i==0) abs(df$exp_var)>0 else {if(i > 0) lead(df$exp_var, i) == 0 else lag(df$exp_var, -i) == 0}
        }
        lim <- lim_ctrl | lim_treat

        # Estimate and Save
        tmp <- suppressMessages(feols(frmla, data = df[lim,], vcov = DK~countrycode))
        lpdid_beta[match(j, -pre_window:post_window)] <- tmp$coeftable[1,1]
        lpdid_se[match(j, -pre_window:post_window)] <- tmp$coeftable[1,2]
        lpdid_n[match(j, -pre_window:post_window)] <- nobs(tmp)
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
          lim_treat <- lim_treat & if (i==0) abs(df$exp_var)>0 else {if(i > 0) lead(df$exp_var, i) == 0 else lag(df$exp_var, -i) == 0}
        }
        lim <- lim_ctrl | lim_treat

    # Estimate
        suppressMessages(tmp <- feols(frmla, data = df[lim,],vcov = DK~countrycode))
        lpdid_beta[match(-j, -pre_window:post_window)] <- tmp$coeftable[1,1]
        lpdid_se[match(-j, -pre_window:post_window)] <- tmp$coeftable[1,2]
        lpdid_n[match(-j, -pre_window:post_window)] <- nobs(tmp)
    }
   }

  coeftable <- data.frame(Estimate = lpdid_beta,
                          "Std. Error" = lpdid_se,
                          "t value" = lpdid_beta/lpdid_se,
                          "Pr(>|t|)" = NA,
                          check.names = FALSE)
  coeftable[,4] <- pnorm(coeftable$`t value`, lower.tail = F)
  return(list(coeftable = coeftable[!is.na(coeftable$Estimate),],
              # df = df,
              window = c(-pre_window:post_window)[!is.na(coeftable$Estimate)],
              nobs = data.frame(window = c(-pre_window:post_window)[!is.na(coeftable$Estimate)],
                                nobs = lpdid_n)))
}
