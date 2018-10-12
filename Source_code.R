
### Source code for:
  ### Schneider PP, van Gool CJAW, Spreeuwenberg P, Hooiveld M, Donker GA, Barnett DJ, Paget J. 
  ###
  ### Using digital epidemiology methods to monitor influenza-like illness 
  ### in the Netherlands in real time: the 2017-2018 season. 
  ###
  ### Until publication, please see the PREPRINT version 
  ### available at: https://doi.org/10.1101/440867 
  ### 
  ### Source code license: CC BY
  ### Contact: project.flutrend@gmail.com
  ### 


# SETUP 
  cat(" \n Initiate Dutch Flu Trend \n")
  # Keep track of computation time
    time.start = Sys.time()

  # Convinient functions
    LoadPackages <- function(package){
      # install (if neccessary) and load all required packages
    for(i in 1:length(package)){
      if(eval(parse(text=paste("require(",package[i],")")))==0) {
        install.packages(package)}}
    return (eval(parse(text=paste("require(",package,")"))))
  }
    RescaleData = function(x,scaled.y = y.full){
      # re-transform data that were scaled and centered
      x.rescaled = x * attr(scaled.y, 'scaled:scale') + attr(scaled.y, 'scaled:center')
      return(x.rescaled)
    }
    PlotRelationships = function(outcome.df = ili, x.df = x.single){
      # Plot individual relationships between predictors and ili incidence
      # Not used in the paper! #
      x.y.df = cbind(outcome.df,x.df)
      plot.single.cor.list = list()
      plot.co.cor.list = list()
      
      for(p in 3:length(x.y.df)){
        temp.df = x.y.df[,c(2,p)]
        names(temp.df) = c("ili","x")
        temp = ggplot(data=temp.df,aes(x=x,y=ili))+
          geom_point()+ 
          ylab("") +
          geom_smooth() +
          xlab(paste(names(x.y.df)[p]))
        
        plot.single.cor.list[[p-2]] = temp
        
      }
      
      relationship.plot = plot_grid(plotlist = plot.single.cor.list,ncol=5)
      return(relationship.plot)
    }
    CustomMetrics = function (data, lev = NULL, model = NULL) {
      # Function to pass to caret::trainControl to retrive values from cross validation holdouts
      # Neccessary, since we decided to use the maximum, insetad of the commonly used mean absolute error
      data.frame(median=median(abs(data$pred- data$obs)),
                 MAE=mean(abs(data$pred- data$obs)),
                 min=min(abs(data$pred- data$obs)),
                 max=max(abs(data$pred- data$obs)),
                 p1 = data$pred[1],
                 p2 = data$pred[2],
                 p3 = data$pred[3],
                 p4 = data$pred[4],
                 o1 = data$obs[1],
                 o2 = data$obs[2],
                 o3 = data$obs[3],
                 o4 = data$obs[4])
      
    }
    AssessCriterionPlot = function(fit = ith.model.fit,criterion.index = "MAE"){
      # Plots predicted values from cross-validation against observed values over the entire tuning horizon 
      # and provides the mean absolute prediction error
      # 'fit' must be a list of class 'train' (from caret)
      # 'criterion.index' is used to specify which metric should be optimized
      crit.pos = which(names(fit$results) == criterion.index)
      lam = fit$results$lambda[crit.pos]
      alp = fit$results$alpha[crit.pos]
      best.index = which(fit$results[,crit.pos] == min(fit$results[,crit.pos]))[1]
      resample.best = fit$resample[fit$resample$lambda == fit$results$lambda[best.index] &
                                     fit$resample$alpha ==fit$results$alpha[best.index], ]
      resample.best = resample.best[,names(resample.best) %in% c(paste("p",1:4,sep=""),paste("o",1:4,sep=""))]
      resample.best = RescaleData(resample.best)
      len = 1:dim(resample.best)[1]
      train.length = length(fit$trainingData[,1])
      init.horizon = fit$control$initialWindow
      
      model = glmnet(y=y.train[1:init.horizon],x=as.matrix(x.train[1:init.horizon,]),lambda = lam, alpha = alp)
      train.preds = predict(model, newx= as.matrix(x.train[1:init.horizon,]))
      train.preds = RescaleData(train.preds)
      
      y.train.re = RescaleData(y.train)
      
      plot.it = 
        ggplot() +
        geom_line(aes(x=date.train,y=y.train.re,col="observed")) +
        geom_line(aes(x=date.train[c(1:init.horizon)],y=train.preds,col="Training")) +
        geom_line(aes(x=date.train[c((init.horizon+1):(train.length-3))],y=resample.best$p1,col="predicted p1")) +
        geom_line(aes(x=date.train[c((init.horizon+2):(train.length-2))],y=resample.best$p2,col="predicted p2")) +
        geom_line(aes(x=date.train[c((init.horizon+3):(train.length-1))],y=resample.best$p3,col="predicted p3")) +
        geom_line(aes(x=date.train[c((init.horizon+4):(train.length))],y=resample.best$p4,col="predicted p4")) +
        scale_color_manual(values = c(1:6)) +
        ggtitle(paste("lambda:",
                      round(fit$results$lambda[criterion.index],2),
                      ",alpha:",
                      round(fit$results$alpha[criterion.index],2),
                      sep=""))
      
      error.df = data.frame(
        W1.error = mean(abs(resample.best$p1-resample.best$o1)), # week 1 prediction error 
        W2.error = mean(abs(resample.best$p2-resample.best$o2)), # week 2...
        W3.error = mean(abs(resample.best$p3-resample.best$o3)),
        W4.error = mean(abs(resample.best$p4-resample.best$o4)))
      
      list(plot.it = plot.it,
           error.df = error.df)
    }
    take.plot.results.snapshot.hot = function(item = results.list$model_2018_01_06, set.xlim = NULL,set.ylim = NULL,sec.col = "white",tit.i=""){
      # Function to plot the real-time prediction results from an individual model
      # specific to the data structure used in this analysis
      plot.results.snapshot.hot.df = item
      date.test = plot.results.snapshot.hot.df$date.test = as.Date(plot.results.snapshot.hot.df$date.test)
      date.train = plot.results.snapshot.hot.df$date.train = as.Date(plot.results.snapshot.hot.df$date.train)
      y.train = plot.results.snapshot.pred.train = predict(plot.results.snapshot.hot.df$ith.model.fit,newdata=as.matrix(plot.results.snapshot.hot.df$x.train))
      x.train = plot.results.snapshot.hot.df$x.train
      plot.results.snapshot.pred.test = predict(plot.results.snapshot.hot.df$ith.model.fit,newdata=as.matrix(plot.results.snapshot.hot.df$x.test))
      plot.results.snapshot.combi.line = c(plot.results.snapshot.pred.train[length(plot.results.snapshot.pred.train)],plot.results.snapshot.pred.test)
      plot.results.snapshot.combi.date = as.Date(c(plot.results.snapshot.hot.df$date.train[length(plot.results.snapshot.hot.df$date.train)],plot.results.snapshot.hot.df$date.test))
      y.combi = c(plot.results.snapshot.hot.df$y.train[length(plot.results.snapshot.hot.df$y.train)],plot.results.snapshot.hot.df$y.test)
      
      plot.results.snapshot.eval = AssessCriterionPlot(plot.results.snapshot.hot.df$ith.model.fit,criterion.index = "MAE")
      title.nw = format(plot.results.snapshot.hot.df$date.test[length(plot.results.snapshot.hot.df$date.test)], format = "Week %V - %Y")
      
      if(is.null(set.xlim)){
           set.xlim = as.Date(c(plot.results.snapshot.hot.df$date.train[length(plot.results.snapshot.hot.df$date.train)-20],plot.results.snapshot.hot.df$date.test[length(plot.results.snapshot.hot.df$date.test)]))
      } else {
        set.xlim = set.xlim
        }
       
      
      date.breakings = as.Date(split.list)
      date.breakings = date.breakings[date.breakings>=set.xlim[1] & date.breakings <= set.xlim[2]]
      current.index = date.breakings == plot.results.snapshot.hot.df$date.test[4]
      current.week = format(plot.results.snapshot.hot.df$date.test[4],format="W%V - %m/%Y") 
      retain.label = rep_len(c(T,F,F,F), length.out=length(date.breakings))
      date.labels = as.character(format(date.breakings ,format="W%V - %m/%Y") )
      date.labels[current.index] = "Current week"
      keep.index = current.index | retain.label
      date.labels = date.labels[keep.index]
      date.breakings.maj = date.breakings[keep.index]
      
   
      plot.results.snapshot=
        ggplot() +
        geom_vline(aes(xintercept = ( plot.results.snapshot.hot.df$date.train[length(plot.results.snapshot.hot.df$date.train)]), 
                       linetype="Last update"),size = 0.5, col ="cyan") +
        geom_line(aes(x=plot.results.snapshot.combi.date,
                      y=RescaleData(y.combi),
                      col = "Observed (validation)"), size =1.5, alpha = 1) + 
        geom_line(aes(x=plot.results.snapshot.hot.df$date.train,
                      y=RescaleData(plot.results.snapshot.hot.df$y.train),
                      col = "Observed (training)"), size =1.5) + 
        geom_line(aes(x=  plot.results.snapshot.combi.date,
                      y = RescaleData(plot.results.snapshot.combi.line),
                      col = "Predicted (validation)"),size =1,alpha = 1) +
        geom_point(aes(x=  plot.results.snapshot.hot.df$date.test[4],
                       y = RescaleData(plot.results.snapshot.pred.test[4])),col=c("red"),alpha = 0.8,
                       size =2) +
      geom_line(aes(x=  plot.results.snapshot.hot.df$date.train,
                    y = RescaleData(plot.results.snapshot.pred.train),
                    col = "Predicted (Training)"),size =1) +
          scale_linetype(name = "") +
          xlab("Week (W) month/year") +
          ylab("Ili incidence per 10,000") +
          theme_minimal()  +
          theme(axis.text.x = element_text(angle = 25, hjust = 1))+
        
          scale_color_manual(values=c("black","gray","orange","red"),
                             name="") +
        scale_x_date(labels = date.labels, #date.labels,
                     minor_breaks = date.breakings,
                     breaks = date.breakings.maj #, limits = set.xlim
                     ) +
          guides(colour = guide_legend(order = 1), 
                 linetype = guide_legend(order = 2)) +
        annotate("text",x=set.xlim[1]+28,y=17,label=paste(tit.i, title.nw,sep=". "),col=c("black"),size=5) 
      
      plot.results.snapshot= plot.results.snapshot+
        annotate("rect",xmin=plot.results.snapshot.hot.df$date.test[4],xmax = plot.results.snapshot.hot.df$date.test[4],
                 ymin= -1,ymax=1,col=c("red")) +
        coord_cartesian(ylim=set.ylim,xlim=set.xlim) +
        theme(axis.text.x=
                element_text(colour=ifelse(date.labels=="Current week","red",sec.col)))
      
      return(plot.results.snapshot)
    }
    
    
  # Load packages
    
    cat("Install and load required packages \n")
    cat("CAVE: Please Check - gtrendsR developer version required! \n")
    required.packages<-c("devtools",
                         "gtrendsR", 
                         "ggplot2",
                         "cowplot",
                         "gridExtra",
                         "corrplot",
                         "scales",
                         "shiny",
                         "lattice",
                         "ISOweek",
                         "reshape2",
                         "Matrix",
                         "caret",
                         "glmnet",
                         "doParallel",
                         "iterators",
                         "shiny")
    LoadPackages(required.packages)
    # Developer version of gtrendsR is required
    # Can be installed from github:
    # install_github("PMassicotte/gtrendsR") 
  
    
# Load ili data 
    
    cat("Load ili data \n")
    ili = read.csv("input/ili_ecdc.csv")
    ili = ili[!is.na(ili$date),] # Week 53 data are removed
    ili$date = as.Date(as.character(ili$date))
    # head(ili)
  
    # Trim ili data (Public Google Trends API only works for periods >= 5 years)
    if((max(ili$date,na.rm=T) - min(ili$date,na.rm=T))/365>5){
      max.span = max(ili$date) - 5*365
      ili = ili[ili$date >= max.span,]
      time.span = paste(min(ili$date),max(ili$date),sep=" ")
    } else {
      time.span = paste(min(ili$date),max(ili$date),sep=" ")
    }
  
    
# Retrieve data from Google Trends 
    
    cat("Load Google Trends data -  \n")
    # Initial keyword to retrieve linked pages 
    init.keyword = "Griep" # (= Dutch word for 'flu')
    cat("Retrieve pages related to '",init.keyword,"' \n",sep="")
    initial_dat = gtrends(keyword = init.keyword,
                          geo="NL",
                          time=time.span)
    related.keywords = initial_dat$related_queries$value[initial_dat$related_queries$related_queries=="top"]
    related.keywords = c(init.keyword,related.keywords)
    # Remove keywords with a year in it for further analysis
    related = related.keywords[!grepl("\\d",related.keywords)] 
    # print(related)
    
    cat("Retrieve interest over time \n")
    raw.gt.data.df = data.frame(date=initial_dat$interest_over_time$date)
    for(key in related){
      cat("\n Retrieve", key)
      tryCatch({
        temp = gtrends(keyword = key,    # Keyword
                       geo="NL",         # ISO code for the Netherlabds
                       time= time.span)  # Time span 
        
        temp_hits = as.numeric(gsub("<","",temp$interest_over_time$hits))
        date.temp = temp$interest_over_time$date
        temp.df = data.frame(date.temp,temp_hits)
        names(temp.df) = c("date",temp$interest_over_time$keyword[1]) 
        raw.gt.data.df = merge(raw.gt.data.df,temp.df, by="date")
      }, error = function(e)cat("\n Couldn't retrieve information on",key," - Continue"))
    }
    
    # Save retrieved data as csv file
    write.csv(raw.gt.data.df,file="output/Google_Trends_data.csv",row.names = F)

    
# Data processing
    
    cat("Processing data \n")
    # single variables
      x.single = raw.gt.data.df[,-1]
      # Plotting individual relationships between predictors and outcome (not used in the paper)
      # x.single.relationship.plot = PlotRelationships(outcome.df = ili, x.df = x.single)
      
    # Assess correlations
      cor.M = cor(raw.gt.data.df[,-1])
      par(mar=c(0,0,0,0),mfrow=c(1,1))
      multi.correlation = cor.M
      for(i in 1:length(multi.correlation[,1])){
        multi.correlation[i,i:length(multi.correlation[1,])] = NA
      }
      multi.correlation.info = t(data.frame(mean = mean(multi.correlation,na.rm=T), 
                                            median = median(multi.correlation,na.rm=T),
                                            min = min(multi.correlation,na.rm=T), 
                                            max= max(multi.correlation,na.rm=T)))
      # plot correlations (not shown in the paper)
      # corrplot(cor.M, method="pie",mar=c(3,0,1,0),tl.cex=0.5,title="Multicolinearity",type="lower")
      
    # creating interactions + quadratic polynomials
      mat.ext = matrix(data=0,
                       ncol= sum(length(x.single):1),
                       nrow=length(raw.gt.data.df[,1]))
      colnames(mat.ext) = rep(1:sum(length(x.single):1))
      cc = 0
      for(i in  1:(length(x.single))){
        for(i2 in i:length(x.single)){
          cc = cc + 1
          mat.ext[,cc] = x.single[,i]*x.single[,i2]
          colnames(mat.ext)[cc] = paste(names(x.single)[i],"XX",names(x.single)[i2],sep="")
        }}
      mat.ext = data.frame(mat.ext)
      
    # combine dataframes
      date.full = raw.gt.data.df$date
      x.full = as.matrix(cbind(x.single,mat.ext))
      y.full = scale(ili$y)
      
    # save processed data
    write.csv(x.full,file="output/Predictor_data.csv",row.names = F)
    
    
## MODEL BUILDING ##
    
# Setup of time series cross-validation
    cat("Modelling setup \n")
    controlObject <- trainControl(method = "timeslice",
                                  allowParallel = TRUE,   
                                  savePredictions = "final", 
                                  returnResamp = "all",   # Save all information from hold out sets
                                  summaryFunction = CustomMetrics, # using custom metric: 'Minimize maximum abs error'
                                  skip = 0,         # Skip for faster computation
                                  horizon = 4, # CV-holdout optimized for a 4 week reporting lag
                                  initialWindow = 52,  # start with 52 weeks of training data
                                  fixedWindow = FALSE) # increasing training data window
    
    # Define tuning grid for tuning paramter lambda
    lasso_grid <- expand.grid(.lambda = c(10 ^ seq(0.25, -8, length= 100) ),.alpha = 1)
    # Evaluation loop index 
    split.list = date.full[(length(raw.gt.data.df$date)-52):length(raw.gt.data.df$date)]
    
# Running outer model loop
    cat("Running outer modelling loop \n")
    # list & index
    results.list = list()
    index = 0
    for(split.i in split.list)
    {
      tryCatch({
        #--> try this one:
        cat("\n Running model for split:",as.character(as.POSIXct(split.i, origin=as.Date("1970-01-01"))))
        index = index + 1
        cat("\n Status:",round(index/length(split.list),4)*100,"% \n")
        split = which(date.full<split.i) 
        
        date.train = raw.gt.data.df$date[split]
        date.test = raw.gt.data.df$date[-split]
        if(length(date.test)>4){
          date.test = date.test[1:4]
        }
        
        x.train = x.full[split,] # Predictor training data set
        y.train = y.full[split] # Outcome for training data set
        
        x.test  = x.full[-split,] # Predictors for testing/evaluation data set
        if(length(date.test)>=4){
          x.test = x.test[1:4,]  
        }
        y.test = y.full[-split] # Outcome for testing data set
        if(length(date.test)>=4){
          y.test = y.test[1:4]
        }
        y.test = y.test[!is.na(y.test)]
        
        # Scaling, centering, transformation and imputation of remaining NAs by K-nearest neighbours
        preprocess.df.train = preProcess(x.train, method=c("scale","center","nzv")) 
        x.train = predict(preprocess.df.train, newdata = x.train)
        x.test = predict(preprocess.df.train,newdata = x.test)
        
        # Model Cross-validation
        ith.model.fit  = train(y= y.train ,
                       x = x.train,
                       method = "glmnet",
                       tuneGrid = lasso_grid,
                       metric = "max", 
                       maximize = FALSE,
                       trControl = controlObject)
        # lambda - CV hold-out results
        cv.error.assessment = plot(ith.model.fit)
        # best lambda
        Best.Tune = ith.model.fit$bestTune 
        # extract coeffficients from best lasso model
        coefs.ith.model.fit = data.frame(var=rownames(coef(ith.model.fit$finalModel,s=ith.model.fit$bestTune$lambda))[as.numeric(coef(ith.model.fit$finalModel,s=ith.model.fit$bestTune$lambda)) !=0],
                                 value=coef(ith.model.fit$finalModel,s=ith.model.fit$bestTune$lambda)[coef(ith.model.fit$finalModel,s=ith.model.fit$bestTune$lambda)!=0])
        coefs.ith.model.fit = coefs.ith.model.fit[order(coefs.ith.model.fit$value,decreasing = T),]
        # plot predicted versus 
        ith.model.fit.eval = AssessCriterionPlot(ith.model.fit,"max")
        # fit final model on training data
        glmnet.fit = glmnet(y = y.train,
                            x = x.train,
                            alpha = ith.model.fit$bestTune$alpha, 
                            lambda = ith.model.fit$bestTune$lambda)
        
        # test model's predictions on validation data set
        preds.test = predict(glmnet.fit,newx = x.test)
        preds.test = data.frame(Week.without = 1:length(preds.test),
                                date=date.test,
                                preds.test = preds.test, 
                                observed = y.test)
        
        # gather results for the ith prediction model
        results.list[[index]] = list(
          ith.model.fit = ith.model.fit,
          cv.error.assessment = cv.error.assessment,
          Best.Tune = Best.Tune,
          coefs.ith.model.fit = coefs.ith.model.fit,
          ith.model.fit.eval = ith.model.fit.eval,
          preds.test = preds.test,
          glmnet.fit = glmnet.fit,
          x.train = x.train,
          y.train = y.train,
          date.train = date.train,
          x.test = x.test,
          y.test = y.test,
          date.test = date.test
        )
        
        # Assign index name
        name.i = as.Date(as.POSIXct(split.i,origin = as.Date("1970-01-01")))
        name.i = paste("model_",gsub("-","_",as.character(name.i)),sep="")
        names(results.list)[index] = name.i
        
      },error = function(e){cat("Something went wrong at",split.i,"\n")}
      )
      
    }

    # save modelling results
    save(results.list,file="output/analysis_loop_result.rdata")
    
  
# Analysis of model results
  results.validation.prediction.df = data.frame(Week.without = NA, date = NA, s0 = NA,observed=NA)
  cv.error.list = results.list[[1]]$ith.model.fit.eval$error.df[1,]
  for(i in 1:length(results.list)){
    results.validation.prediction.df = rbind(results.validation.prediction.df,results.list[[i]]$preds.test)
    cv.error.list = rbind(cv.error.list,results.list[[i]]$ith.model.fit.eval$error.df)
    rownames(cv.error.list)[i+1] = paste("cv",i,sep="")
  }
  
  # Assess prediction errors
    # Error in validation set 
    results.validation.prediction.df = results.validation.prediction.df[-1,]
    results.validation.prediction.df$Week.without = as.factor(results.validation.prediction.df$Week.without)
    results.validation.prediction.df$date = as.POSIXct(results.validation.prediction.df$date,origin = as.Date("1970-01-01"))
    results.validation.prediction.df$predicted = RescaleData(results.validation.prediction.df$s0)
    results.validation.prediction.df$observed.rescaled = RescaleData(results.validation.prediction.df$observed)
    results.validation.prediction.df$error.rescaled = results.validation.prediction.df$predicted - results.validation.prediction.df$observed.rescaled
    results.validation.prediction.df$abs.error.rescaled = abs(results.validation.prediction.df$predicted - results.validation.prediction.df$observed.rescaled)
    # Mean absolute error in validation set
    mae.validation = as.numeric(by(data=c(results.validation.prediction.df$predicted - results.validation.prediction.df$observed.rescaled),INDICES = results.validation.prediction.df$Week.without,FUN = function(x) mean(abs(x))))
    mae.validation = data.frame(week = 1:4,abs.error.rescaled= round(mae.validation,2))
    # Maximum absolute error in validation set
    maxae.validation = results.validation.prediction.df[results.validation.prediction.df$abs.error.rescaled ==  max(results.validation.prediction.df$abs.error.rescaled),]
    maxae.validation = data.frame(value = maxae.validation$abs.error.rescaled,date = format(maxae.validation$date,format="%V/%Y"),week = maxae.validation$Week.without)
    # Error in CV hold outs
    cv.error.list = cv.error.list[-1,]
    cv.error.list$cv.loop = 1:length(cv.error.list$W1.error)
    cv.error.list = melt(cv.error.list,id.vars=c("cv.loop"))
    cv.error.list$value.rescaled =  RescaleData(cv.error.list$value)
    cv.error.stats = aggregate(value ~ variable, cv.error.list, mean)
    cv.error.stats$value = round(cv.error.stats$value,3)
    cv.error.max = aggregate(value ~ variable, cv.error.list, max)
    cv.error.max$value = round(cv.error.max$value,2)
    # Training error in first iteration model
    mae.training = mean(abs(RescaleData(results.list[[1]]$y.train) - RescaleData(predict(results.list[[1]]$glmnet.fit, newx = results.list[[1]]$x.train))))
    maxae.training = max(abs(RescaleData(results.list[[1]]$y.train) - RescaleData(predict(results.list[[1]]$glmnet.fit, newx = results.list[[1]]$x.train))))
    # Week 4 predictions as ili forecast (MAE in validation set)
    results.future.prediction.df = results.validation.prediction.df[results.validation.prediction.df$Week.without==4,]
    l1 = length(results.future.prediction.df$observed.rescaled)
    mae.future = mean(abs(results.future.prediction.df$observed.rescaled[2:l1] - results.future.prediction.df$predicted[1:(l1-1)]))
    # Correlations between observed and predicted values in the validation set
    cor.validation.p = numeric(); cor.validation.sm = numeric()
    for(w in levels(results.validation.prediction.df$Week.without)){
      temp = cor(results.validation.prediction.df$predicted[results.validation.prediction.df$Week.without==w],
                 results.validation.prediction.df$observed.rescaled[results.validation.prediction.df$Week.without==w],
                 method = "pearson")
      cor.validation.p = c(cor.validation.p,temp)
      temp = cor(results.validation.prediction.df$predicted[results.validation.prediction.df$Week.without==w],
                 results.validation.prediction.df$observed.rescaled[results.validation.prediction.df$Week.without==w],
                 method = "spearman") 
      cor.validation.sm = c(cor.validation.sm, temp) 
    }
    cor.validation = data.frame(
      cor.validation.p  = formatC(cor.validation.p,digits=2,format="f"),
      cor.validation.sm = formatC(cor.validation.sm,digits=2,format="f")
    )
  
# First loop training predictions
    first.loop.pred = predict(results.list[[1]]$glmnet.fit, newx = results.list[[1]]$x.train)
    first.loop.pred.rescaled = RescaleData(first.loop.pred)
  
# Create Figure 1
  # Figure 1 - Overview
  plot.results.overview = 
    ggplot() +
    geom_line(aes(x=date.full,y=RescaleData( y.full),col="Observed")) +
    geom_line(aes(x=date.full[1:length(first.loop.pred)],y=first.loop.pred.rescaled,col="Training")) +
    geom_line(aes(x=results.validation.prediction.df$date,y=results.validation.prediction.df$predicted,col=results.validation.prediction.df$Week.without)) +
    geom_vline(xintercept = min(split.list),col="cyan") +
    xlab("Year") +
    ylab("Ili incidence per 10,000") +
    scale_color_manual(values=c("darkgreen","blue","purple","red",1,"orange"),
                       labels=c("Week 1 prediction","Week 2 prediction",
                                "Week 3 prediction","Week 4 prediction",
                                "Observed ili incidence","Training prediction"),
                       name="Legend") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),plot.margin = margin(1,1,10,30)) +
    scale_x_datetime(labels = date_format("%Y")) +
    scale_y_continuous(breaks=c(0,10,20))
  
  # Figure 1 - Zoom in 2017/18 season
  plot.results.zoom = 
    plot.results.overview + 
    theme(axis.text.x = element_text(angle = 35, hjust = 0.7, vjust = 0.8),plot.margin = margin(1,1,40,30)) +
    scale_y_continuous(breaks=c(0,5,10,15,20)) +
    scale_x_datetime(labels = date_format("W%V - %m/%Y"),
    minor_breaks = c(split.list),
    limits = c(min(split.list),max(split.list))) +
    xlab("Week (W) - month/year") 
  
  # set combined limits
  f.p.z = plot.results.zoom + theme(legend.position = "right") + ggtitle("Predicting the 2017/18 ili epidemic in real time")
  f.p.o = plot.results.overview  + ggtitle("Overview: Training and validation period") + theme_minimal() + theme(axis.text.x = element_text(angle = 30, hjust = 1),plot.margin = margin(1,30,10,35), legend.position = "none") 
  get.legend2 = get_legend(plot.results.zoom + theme(legend.position = "top")) # retrieve legend from plot
  g1=ggplot_gtable(ggplot_build(f.p.z))
  g2=ggplot_gtable(ggplot_build(f.p.o))
  max.width2 = grid::unit.pmax(g1$widths[2:3], g2$widths[2:3])
  g1$widths[2:3] <- as.list(max.width2)
  g2$widths[2:3] <- as.list(max.width2)
  
  # combining plots
  combined.plots = 
    grid.arrange(
      arrangeGrob(g1,g2,
                  #get.legend2,
                  nrow=2,heights=c(4,2))
    )
  
  
# Figure 2: Temporal dynamics - model plot.results.snapshot.hots
  # create 'plot.results.snapshot.hots'
  plot.results.snapshot.1 = take.plot.results.snapshot.hot(tit.i = "A",results.list$model_2017_12_02,set.ylim=c(0,20), set.xlim = c(as.Date("2017-11-05"),as.Date("2018-05-14"))) + xlab("") + ylab("") + theme(legend.position = "none", axis.text.x = element_blank()) 
  plot.results.snapshot.2 = take.plot.results.snapshot.hot(tit.i = "B",results.list$model_2018_01_13,set.ylim=c(0,20), set.xlim = c(as.Date("2017-11-05"),as.Date("2018-05-14"))) + xlab("") + ylab("") + theme(legend.position = "none", axis.text.x = element_blank()) 
  plot.results.snapshot.3 = take.plot.results.snapshot.hot(tit.i = "C",results.list$model_2018_02_24,set.ylim=c(0,20), set.xlim = c(as.Date("2017-11-05"),as.Date("2018-05-14"))) + xlab("") + ylab("") + theme(legend.position = "none", axis.text.x = element_blank())
  plot.results.snapshot.4 = take.plot.results.snapshot.hot(tit.i = "D",results.list$model_2018_03_24,set.ylim=c(0,20), set.xlim = c(as.Date("2017-11-05"),as.Date("2018-05-14"))) + xlab("") + ylab("") + theme(legend.position = "none", axis.text.x = element_blank()) 
  plot.results.snapshot.5 = take.plot.results.snapshot.hot(tit.i = "E",results.list$model_2018_04_07,set.ylim=c(0,20), set.xlim = c(as.Date("2017-11-05"),as.Date("2018-05-14")),sec.col = "black") + xlab("") + ylab("") + theme(legend.position = "none")
  get.legend = get_legend(plot.results.snapshot.5 + theme(legend.position = "top") + guides(colour = guide_legend(nrow = 2)))
  
  # Set combined limits 
  gA=ggplot_gtable(ggplot_build(plot.results.snapshot.1))
  gB=ggplot_gtable(ggplot_build(plot.results.snapshot.2))
  gC=ggplot_gtable(ggplot_build(plot.results.snapshot.3))
  gD=ggplot_gtable(ggplot_build(plot.results.snapshot.4))
  gE=ggplot_gtable(ggplot_build(plot.results.snapshot.5))
  max.width = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3],gC$widths[2:3],gD$widths[2:3],gE$widths[2:3])
  gA$widths[2:3] <- as.list(max.width)
  gB$widths[2:3] <- as.list(max.width)
  gC$widths[2:3] <- as.list(max.width)
  gD$widths[2:3] <- as.list(max.width)
  gE$widths[2:3] <- as.list(max.width)
  
  # combine plots
  # grid.newpage()
  # pdf('graphs/test.pdf',width=8,height=6)
  plot.results.snapshot.hots = 
    grid.arrange(
      arrangeGrob(gA,gB,gC,gD,gE,get.legend,nrow=6,heights=c(1,1,1,1,1.2,0.5))
    )
  
  # Temporal dynamics
  # begin
  epidemic.start.date = as.Date(results.validation.prediction.df$date[ min(which(results.validation.prediction.df$observed.rescaled>=5))])
  epidemic.start.week = format(epidemic.start.date, format = "%V/%Y")
  predicted.start.date = as.Date(results.validation.prediction.df$date[ min(which(results.validation.prediction.df$predicted>=5))])
  predicted.start.week = format(predicted.start.date, format = "%V/%Y")
  # peak
  epidemic.peak.df = results.validation.prediction.df[!duplicated(results.validation.prediction.df$date),]
  epidemic.peak.df = epidemic.peak.df[order(epidemic.peak.df$observed.rescaled, decreasing = T),]
  epidemic.peak.date = as.Date(epidemic.peak.df$date[1])
  epidemic.peak.value = round(epidemic.peak.df$observed.rescaled[1],2)
  epidemic.peak.week = format(epidemic.peak.date, format = "%V/%Y")
  # 2nd peak
  epidemic.2nd.peak.date = as.Date(epidemic.peak.df$date[2])
  epidemic.2nd.peak.value = round(epidemic.peak.df$observed.rescaled[2],2)
  epidemic.2nd.peak.week = format(epidemic.2nd.peak.date, format = "%V/%Y")
  predicted.peak.date = min(as.Date(results.validation.prediction.df$date[ which(results.validation.prediction.df$predicted == max(results.validation.prediction.df$predicted))]))
  predicted.peak.week = format(predicted.peak.date, format = "%V/%Y")
  # end
  epidemic.end.date = min(as.Date(results.validation.prediction.df$date)[as.Date(results.validation.prediction.df$date) > epidemic.peak.date & results.validation.prediction.df$observed.rescaled <=5])
  epidemic.end.week = format(epidemic.end.date, format = "%V/%Y")
  predicted.end.date = min(as.Date(results.validation.prediction.df$date)[as.Date(results.validation.prediction.df$date) > predicted.peak.date & results.validation.prediction.df$predicted <=5])
  predicted.end.week = format(predicted.end.date, format = "%V/%Y")
  
  
# Figure 3: Coefficients in final models
  # lambda ad coefs over hold outs
    lambda.list = numeric()
    coef.list = matrix(coef(results.list[[1]]$glmnet.fit))
    for(i in 1:length(results.list)){
      lambda.list[i] = results.list[[i]]$Best.Tune$lambda
      coef.list = cbind(coef.list,coef(results.list[[i]]$glmnet.fit))
    }
    coef.list = coef.list[,-1]
    coef.list = coef.list[as.logical(apply(coef.list,1,function(x)sum((x!=0))>0)),]
    coef.list[coef.list == 0] = 0
    coef.list = as.matrix(coef.list)
    coef.list = data.frame(coef.list)
    coef.list$var = rownames(coef.list)
    names(coef.list ) = c(1:52,"Variable")
    coef.list = melt(coef.list)
    names(coef.list)= c("Variable","Hold.out","value")
    coef.list$Hold.out = as.numeric(coef.list$Hold.out)
    coef.list = coef.list[-which(coef.list$value==0),]
    n.coef.list = aggregate(Variable ~ as.factor(Hold.out), data = coef.list, FUN= function(x) sum(x>0.000001) )
    coef.list$Variable = ifelse(grepl("XX",coef.list$Variable),paste("Interaction:",gsub("XX","/",coef.list$Variable)),coef.list$Variable)
    coef.list$Variable[coef.list$Variable=="(Intercept)"] = "Intercept"
    
    # Create figure 3 - coefficients retained
    coefficient.importance = ggplot(coef.list) +
      geom_line(aes(x=Hold.out,y=value,col=Variable))+
      geom_point(aes(x=Hold.out,y=value,col=Variable)) +
      ylab("Coefficient value") +
      xlab("Iteration/week in validation period") +
      theme_minimal() 
  
    rle.lambda = rle(lambda.list) # number of changes in lambda over time
    rle.lambda = length(rle.lambda$lengths)-1

  
# Time information objects for the manuscript
    time.span_total = c(min(date.train),max(split.list))
    time.span_total = format(as.Date(time.span_total), format="%V/%Y")
    time.span_training_only = as.Date(date.train[!date.train %in% split.list])
    time.span_training_only = c(min(time.span_training_only),max(  time.span_training_only))
    time.span_training_only = format(as.Date(time.span_training_only), format="%V/%Y")
    time.span_loop = c(min(split.list),max(split.list))
    time.span_loop = format(as.Date(time.span_loop), format="%V/%Y")

#### Finishing
  ## SAVE the entire working enviroment as an R-object
  cat("Saving everything \n")
  save(list=ls(),file=paste("output/Savegame.",Sys.Date(),".rdata",sep=""))
  
  # finish up 
  time.end = Sys.time() 
  comp.duration = time.end - time.start
  comp.duration = paste(round(as.numeric(comp.duration),2),attributes(comp.duration)$units)
  cat("\n \n Analysis results ready! \n   Time elapsed",comp.duration," \n ")

#### Produce the manuscript from RMarkdown
  cat("\n Producing the manuscript \n ")
  
  # render PDF
  rmarkdown::render("manuscript/manuscript_source_code.Rmd", 
                    output_format = c("pdf_document"),
                    output_file = c("manuscript.pdf"),
                    output_dir = c("manuscript"),
                    clean = T)
  
  cat("\n Finished.") 
  
##
## Please send your comments and feedback to: project.flutrend@gmail.com
##  
  
