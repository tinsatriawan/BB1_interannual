


dir <- "G:/.shortcut-targets-by-id/1txCh9lZ7VGCujXGvBaJCMVnuxT-65q4K/Micromet Lab/Projects/2014-BB1 Burns\ Bog/Flux-tower"

output_dir <- paste0(dir, "/flux_data")


# variable we need for LE & H gap-filling
input <- read.table(paste0(dir, "/flux_data/REddyProc_input/for_gap_filling_partitioning.txt"), header = T)

# Define predictors
predictors <- c("Rg", "PAR", "Ustar", "Tair", "rH", "VPD","DoY", "WTH") 

variable <- "LE"


gapfill_RF <- function(input,
                       n_iteration, 
                       predictors, 
                       variable, 
                       variable_name) {# Delete first row that contains units
  input <- input[-1, ]
  input <- data.frame(lapply(input, function(x) as.numeric(as.character(x))))
  
  input$HH <- floor(input$Hour)
  input$MM <- (input$Hour-input$HH)*60
  
  # Create time stamp
  input$TIMESTAMP_END <- make_datetime(input$Year, 1, input$DoY, input$HH, input$MM)
  
  # select predictor
  ML.df <- input %>% select(predictors)
  
  # Replace all -9999 with NA
  ML.df[ML.df == -9999] <- NA
  
  
  # Add sine and cosine functions
  ML.df$s <- sin((ML.df$DoY-1)/365*2*pi)
  ML.df$c <- cos((ML.df$DoY-1)/365*2*pi)
  
  # period when dependent variable is not missing
  wm_only <- ML.df[!is.na(ML.df[ , variable]), ]
  
  
  ############### Random forest run 20x ###############
  
  #Add parallel processing for the fast processing if you want to
  # library(parallel)
  # library(doParallel)
  
  combined_result <- list()
  
      for(i in 1:n_iteration){
        start_time <- Sys.time()
        # setting seed
        set.seed(i)
        train_rows <- sample(1:nrow(wm_only), 0.75*nrow(wm_only))
        # select the training set
        train_set <- wm_only %>% slice(train_rows)
        # select the validation set
        test_set <- anti_join(wm_only, train_set)
        #### option 1. random forest model with mtry tuning
        # cluster <- makeCluster(6)
        # cluster <- parallel::makeCluster(5, setup_timeout = 0.5)
        # registerDoParallel(cluster)
        
        formula <- as.formula(paste0(variable, "~ ."))
        RF_result <- train(formula, data = train_set[, predictors], 
                             method = "rf",
                             preProcess = c("medianImpute"),                #impute missing met data with median
                             trControl=trainControl(method = "cv",   #three-fold cross-validation for model parameters 3 times
                                                    number = 3),
                             na.action = na.pass,
                             allowParallel=FALSE, # This requires parallel packages. Otherwise you can choose FALSE.
                             ntree=400, # can generate more trees
                             importance = TRUE)
        # eval(parse(text = paste0('RF_variable <- train(', variable ,' ~ ., data = train_set[,predictors],
        #                                               method = "rf",
        #                                               preProcess = c("medianImpute"),           
        #                                               trControl=trainControl(method = "cv",  
        #                                                                      number = 3),
        #                                               na.action = na.pass,
        #                                               allowParallel=FALSE, 
        #                                               ntree=400, 
        #                                               importance = TRUE)')))
        # 
        RF_result$bestTune
        RF_result$results
        
        
        # save result
        # saveRDS(RF_result, paste0(output_dir, "/RF_results/RF_NEE_",i))
        eval(parse(text = paste0('saveRDS(RF_result, paste0(output_dir, "/RF_results/RF_', variable ,'_",i))')))
        
        ############### Results
        # variable importance
        
        eval(parse(text = paste0('png(paste0(output_dir, "/RF_results/', variable,'_variable_importance_",i,".png"))')))
        print(plot(varImp(RF_result, scale = FALSE), main="variable importance"))
        dev.off()
        
        #generate NEE_rf predictions for testset
        test_set <- test_set %>% 
          mutate(gf = predict(RF_result, test_set, na.action = na.pass))
        regrRF <- lm(as.formula(paste0("gf ~", variable)), data = test_set); 
        print(summary(regrRF))
        test_set %>% 
          ggplot(aes(x={{variable_name}}, y=test_set$gf)) + 
          geom_abline(slope = 1, intercept = 0)+
          geom_point() + 
          geom_smooth(method = "lm") + 
          ggtitle("testset")
        
        # whole dataset
        RF_result_whole <- data.frame(ngf = ML.df[, variable]) # you can add datetime column here if you want to.
        RF_result_whole$model <- predict(RF_result, ML.df, na.action = na.pass) # NEE RF model
        RF_result_whole$filled <- ifelse(is.na(RF_result_whole$ngf),RF_result_whole$model, RF_result_whole$ngf) # gap-filled column (true value when it is, gap-filled value when missing)
        RF_result_whole$residual <- ifelse(is.na(RF_result_whole$ngf),NA,RF_result_whole$model - RF_result_whole$ngf) # residual (model - obs). can be used for random uncertainty analysis
        
        # time series
        RF_result_whole$DateTime <- Input$TIMESTAMP_END
        
        RF_result_whole %>% 
          ggplot(aes(DateTime,ngf)) + 
          geom_point() + 
          theme_bw() + 
          ylab(expression(paste(variable,"(umol ", m^-2,s^-1,")"))) %>% 
          print()
        
        RF_result_whole %>% 
          ggplot(aes(DateTime,filled)) + 
          geom_point(color="red",alpha=0.5) +
          geom_point(aes(DateTime,ngf),color="black") +
          theme_bw() + 
          ylab(expression(paste(variable, "(umol ", m^-2,s^-1,")"))) %>% 
          print()
        
          # whole data comparison
          print(ggplot(RF_result_whole, 
                       aes(x = ngf, y = model)) + 
                  geom_abline(slope = 1, intercept = 0) +
                  geom_point() + 
                  geom_smooth(method = "lm") + 
                  ggtitle("whole dataset"))
          
          
          regrRF_whole <- lm(RF_result_whole$model ~ RF_result_whole$ngf);
          print(summary(regrRF_whole))
          
          
          RF_result_whole$iteration <- i
          combined_result[[i]] <- RF_result_whole
          end_time <- Sys.time()
          print(end_time - start_time)
        }
    
  
  rf_result_df <- data.table::rbindlist(combined_result)
  eval(parse(text = paste0('saveRDS(rf_result_df, paste0(output_dir, "/RF_results/BB_', variable, '_rf_result"))')))
  eval(parse(text = paste0('write_csv(rf_result_df, paste0(output_dir, "/RF_results/BB_', variable, '_rf_result.csv"))')))
  }


gapfill_RF(input = read.table(paste0(dir, "/flux_data/REddyProc_input/for_gap_filling_partitioning.txt"), header = T),
           n_iteration = 1,
           predictors = c("Rg", "PAR", "Ustar", "Tair", "rH", "VPD","DoY", "WTH", "LE"), 
           variable = "LE", 
           variable_name = LE)

gapfill_RF(input = read.table(paste0(dir, "/flux_data/REddyProc_input/for_gap_filling_partitioning.txt"), header = T),
           n_iteration = 1,
           predictors = c("Rg", "PAR", "Ustar", "Tair", "rH", "VPD","DoY", "WTH", "H"), 
           variable = "H", 
           variable_name = H)

