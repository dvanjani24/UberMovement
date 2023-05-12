###### Set-up
libraries <- c("tidyverse", "lubridate", "ggplot2", "ggfortify", "fpp",
               "stringr", "forecast", "scales", "zoo", "imputeTS", "vars",
               "MuMIn", "keras", "TSLSTM")
load <- lapply(libraries, require, character.only = TRUE)


###### 1. Load Data
get_file_names <- function(data_path){
  file_names <- list.files(path = data_path, pattern = "*.csv")
  file_times <- file.mtime(file.path(data_path, file_names))
  file_names_sorted <- file_names[order(file_times, decreasing = FALSE)]
}

load_data <- function(data_path){
  destination_dirs <- list.files(path = data_path)
  df_full <- data.frame()
  for (destination in destination_dirs){
    dest_path = paste0(data_path,"/",destination)
    file_list <- get_file_names(dest_path)
    for (file in file_list){
      tryCatch({
        temp <- read.csv(paste0(dest_path,"/",file))
        re <- "(?<=\\.\\s)\\D+(?=\\s*\\()|(?<=\\.\\s)\\D+$"
        temp$dest_name <- rep(trimws(str_extract_all(destination, re)),nrow(temp))
        df_full <- rbind(df_full, temp)
      }, error = function(e) {
        cat(paste("Warning:", "No data present in this file, continuing..."), "\n")
      })
    }
  }
  df <- df_full[,c(1,2,4,6,ncol(df_full))]
  names(df) <- c("date", "origin", "destination", "mean_time", "neighborhood")
  df$date <- as.Date(df$date, format = "%m/%d/%Y")
  df_list <- split(df, df$neighborhood)
  df_list_sorted <- lapply(df_list, function(x) x[order(x$date),])
  return(df_list_sorted)
}

###### 2. Preliminary Data Cleaning
# Apply weekly averaging
get_weekly_avg <- function(route){
  route$week <- floor_date(route$date, "week")
  weekly_df <- route %>%
    group_by(week) %>%
    summarize(mean_travel_time = mean(mean_time, na.rm = TRUE))
  return(weekly_df)
}

# Filter routes that have a high number of missing values
missing_breach <- function(weekly_df, threshold = 0.15){
  prop_missing = sum(is.na(weekly_df$mean_travel_time))/nrow(weekly_df)
  cat(paste0("Missing Values (%): ", round(prop_missing,4)*100, "%", "\n"))
  if (prop_missing > threshold){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Apply some basic cleaning
clean_df <- function(weekly_df){
  # 1. Subset series to end before March 2020 to avoid covid impact
  end <- "2020-02-23"
  weekly_df <- weekly_df[weekly_df$week <= end, ]
  
  # 2. Convert to data to minutes from seconds
  weekly_df$mean_travel_time <- weekly_df$mean_travel_time/60
  
  # # 3. Impute missing values (if necessary) using last value carried forward
  # cat(paste0("Missing Values (%): ", round(sum(is.na(weekly_df$mean_travel_time))/nrow(weekly_df),4)*100,"%","\n"))
  data_missing <- sum(is.na(weekly_df$mean_travel_time)) > 0
  weekly_df$mean_travel_time <- if(data_missing){na_locf(weekly_df$mean_travel_time)} else {weekly_df$mean_travel_time}
  cat("Imputing missing values using LOCF, if necessary...\n")
  cat(paste("Confirm No Missing Values:", sum(is.na(weekly_df$mean_travel_time))==0), "\n")
  return(weekly_df)
}

# Plot series
plot_series <- function(weekly_df){
  p <- ggplot(weekly_df) + 
    geom_line(aes(x = week, y = mean_travel_time)) +
    scale_x_date(labels = function(x) format(x, "%b-%Y")) + 
    ggtitle(paste0("Weekly Average Travel Times from LAX (734) to ",
                   unique(route))) +
    xlab("") + 
    ylab("Travel Time (Minutes)") +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 9))
  print(p)
}

###### 3. Diagnosing Stationarity
# Check for stationarity
check_stationarity <- function(weekly_df, adf_thresh = 0.05){
  
  is_stationary <- (adf.test(weekly_df$mean_travel_time)$p.value < adf_thresh)
  return(is_stationary)
}

model <- function(df, features, stationary){
  ts_data <- ts(df$mean_travel_time,
                start = c(year(df$week[1]), as.numeric(format(as.Date(df$week[1]), "%U"))),
                frequency = 52)
  
  # ARIMA
  cat("Running ARIMA...\n")
  arima_fit <- auto.arima(ts_data)
  arima_order <- arimaorder(arima_fit)
  arima_spec <- paste0("ARIMA(",arima_order[1], ",", arima_order[2], ",", arima_order[3],")",
                 "(", arima_order[4], ",", arima_order[5], ",", arima_order[6], ")",
                 "[", arima_order[7], "]")
  arima_df <- data.frame(`Model Type` = "ARIMA",
                         "Specification" = arima_spec,
                         AICc = arima_fit$aicc,
                         AIC = arima_fit$aic, 
                         BIC = arima_fit$bic,
                         RMSE = accuracy(arima_fit)[,'RMSE'])
  
  # Exponential Smoothing
  cat("Running Exponential Smoothing...\n")
  ets_fit <- ets(df$mean_travel_time)
  ets_spec <- ets_fit$method
  ets_df <- data.frame(`Model Type` = "Exponential Smoothing",
                       "Specification" = ets_spec,
                       AICc = ets_fit$aicc,
                       AIC = ets_fit$aic, 
                       BIC = ets_fit$bic, 
                       RMSE = accuracy(ets_fit)[,'RMSE'])

  # Linear Regression
  cat("Running Linear Regression...\n")
  regression_df <- left_join(df, features, by = 'week')
  full_lm <- lm(mean_travel_time~.-week, data = regression_df)
  best_lm <- step(full_lm, direction = 'both', trace = FALSE)
  chosen_vars <- all.vars(best_lm$call)[3:(length(all.vars(best_lm$call))-1)]
  chosen_vars_read <- paste("X =", paste(chosen_vars, collapse = ", "))
  lm_df <- data.frame(`Model Type` = "Linear Regression",
                      "Specification" = chosen_vars_read,
                      AICc = AICc(best_lm),
                      AIC = AIC(best_lm),
                      BIC = BIC(best_lm), 
                      RMSE = accuracy(best_lm)[,'RMSE'])
  
  # Linear Regression with Arima Errors
  regression_ts <- ts(regression_df,
                      start = c(year(regression_df$week[1]),
                                as.numeric(format(as.Date(regression_df$week[1]), "%U"))),
                      frequency = 52)
  cat("Running Linear Regression with Arima Errors...\n")
  arimax_fit <- auto.arima(regression_ts[,'mean_travel_time'],
                           xreg = regression_ts[,chosen_vars])
  arimax_order <- arimaorder(arimax_fit)
  arimax_spec <- paste0("ARIMA(",arimax_order[1], ",", arimax_order[2], ",", arimax_order[3],")",
                        "(", arimax_order[4], ",", arimax_order[5], ",", arimax_order[6], ")",
                        "[", arimax_order[7], "]")
  arimax_df <- data.frame(`Model Type` = "Linear Regression with Arima Errors",
                          "Specification" = paste0(chosen_vars_read, " (Errors: ",arimax_spec, ")"),
                          AICc = AICc(arimax_fit),
                          AIC = AIC(arimax_fit),
                          BIC = BIC(arimax_fit),
                          RMSE = accuracy(arimax_fit)[,'RMSE'])
  
  # Spectral Analysis
  tbats_fit <- tbats(ts_data)
  tbats_spec <- paste0("TBATS(",
                       ifelse(is.null(tbats_fit$lambda), 1, tbats_fit$lambda),
                       ",{", tbats_fit$parameters$control$p, ",",
                       tbats_fit$parameters$control$q, "},",
                       ifelse(is.null(tbats_fit$damping.parameter), "-", round(tbats_fit$damping.parameter,2)), ",",
                       ifelse(is.null(tbats_fit$seasonal.periods), "NULL", tbats_fit$seasonal.periods), ")")
  tbats_df <- data.frame(`Model Type` = "TBATS",
                         "Specification" = tbats_spec,
                         AICc = NA,
                         AIC = tbats_fit$AIC,
                         BIC = NA,
                         RMSE = accuracy(tbats_fit)[,'RMSE'])
  
  model_summary_df <- rbind(arima_df, ets_df, lm_df, arimax_df, tbats_df)
  return(model_summary_df)
}

###### Execute Pipeline
data_path <- paste0(getwd(),"/Data")
df_list_all <- load_data(data_path)

features <- read.csv("features.csv")
features$week <- as.Date(features$week, format = "%m/%d/%y")

# Subset for testing
df_list <- df_list_all[c(1:10)]
# df_list <- df_list_all

neighborhoods <- names(df_list)

models_df <- rbind()
for (route in neighborhoods){
  cat("***************************************************\n")
  cat(paste0("Building models for LAX to ", route, "..."), "\n")
  weekly_df <- get_weekly_avg(df_list[[route]])
  if(missing_breach(weekly_df)){
    cat(paste0("Removing ", route, " from analysis, as too many missing values\n"))
    next
  }
  else {
    weekly_df <- clean_df(weekly_df)
    plot_series(weekly_df)
    cat(paste("Series is stationary:",check_stationarity(weekly_df)), "\n")
    temp <- cbind(Destination = route,
                  Stationary = check_stationarity(weekly_df),
                  model(weekly_df, features = features))
    models_df <- rbind(models_df, temp)
  }
}

write.csv(models_df, "models_df.csv")
write_excel_csv(models_df, "models_df.xslx")


#### LSTM

df <- read.csv(paste0(getwd(),"/Sample data to model/Adams-Normandie.csv"))
(is_stationary <- adf.test(df$mean_travel_time)$p.value < 0.05)
(is_stationary <- adf.test(diff(df$mean_travel_time))$p.value < 0.05)

series <- diff(df$mean_travel_time)
lstm_model <- ts.lstm(diff(df$mean_travel_time), tsLag = 12, LSTMUnits = 10, Epochs = 50)

## Running CV
arima_cv <- tsCV(ts_data, function(x, h){forecast(arima_fit, h=h)}, h = h)
rmse <- sqrt(colMeans(arima_cv^2, na.rm = TRUE))
