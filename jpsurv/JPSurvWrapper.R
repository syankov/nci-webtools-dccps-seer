library('rjson')
library('JPSurv')
VERBOSE=TRUE

getDictionary <- function (inputFile, path, tokenId) {
  fqFileName = file.path(path, inputFile)
  outputFileName = paste("form-", tokenId, ".json", sep="")
  fqOutputFileName = file.path(path, outputFileName)
  seerFormData = dictionary.overview(fqFileName)
  cat(toJSON(seerFormData), file = fqOutputFileName)
  return(tokenId)
}

#Creates the subset expression for fitted result
getSubsetStr <- function (yearOfDiagnosisVarName, yearOfDiagnosisRange, cohortVars, cohortValues) {

  yearOfDiagnosisVarName=paste0("`",getCorrectFormat(yearOfDiagnosisVarName), "`")
  startYearStr=paste(yearOfDiagnosisVarName, ">=", yearOfDiagnosisRange[1])
  endYearStr=paste(yearOfDiagnosisVarName, "<=", yearOfDiagnosisRange[2])
  yearStr=paste(startYearStr, endYearStr, sep='&')
  cohortVars=paste0("`",getCorrectFormat(cohortVars), "`")
  
  subsetStr=paste(paste(cohortVars, cohortValues, sep="=="), collapse='&')
  subsetStr=paste(subsetStr, yearStr, sep='&')


  return (subsetStr)

}

#Creates the model.form expression for fitted result
getFactorStr <- function (covariateVars) {
  factorStr=""
  if (nchar(covariateVars)!=0) {
    covariateVars=paste0("`", getCorrectFormat(covariateVars), "`")
    factorStr=paste("~-1+", paste(gsub("$", ")", gsub("^", "factor(", covariateVars)), collapse="+"), sep='')
  }

  return (factorStr)
}

#replace empty space with _, strip anything other than alphanumeric _ /
getCorrectFormat <-function(variable) {
  variable=gsub("[^[:alnum:]_/]", "", gsub(" ", "_", variable))
  return (variable)
}


jpsurvData = list()

#Parses the JSON string and sends to getFittedResult to create the SEER Data and the Fitted Results
getFittedResultWrapper <- function (filePath, jpsurvDataString) {

  jpsurvData <<- fromJSON(jpsurvDataString)

  seerFilePrefix = jpsurvData$calculate$static$seerFilePrefix
  yearOfDiagnosisVarName = jpsurvData$calculate$static$yearOfDiagnosisVarName
  yearOfDiagnosisRange = jpsurvData$calculate$form$yearOfDiagnosisRange
  allVars=jpsurvData$calculate$static$allVars
  cohortVars=jpsurvData$calculate$form$cohortVars
  cohortValues=jpsurvData$calculate$form$cohortValues
  covariateVars=jpsurvData$calculate$form$covariateVars
  numJP=jpsurvData$calculate$form$maxjoinPoints
  numbetwn=jpsurvData$calculate$form$numbetwn 
  numfromstart=jpsurvData$calculate$form$numfromstart 
  numtoend=jpsurvData$calculate$form$numtoend
  adanced_options=list(numbetwn,numfromstart,numtoend)
  
  fileName = paste('output', jpsurvData$tokenId, sep="-" )
  fileName = paste(fileName, "rds", sep="." )
  outputFileName = fileName

  fileName = paste('output', jpsurvData$tokenId, sep="-" )
  fileName = paste(fileName, "rds", sep="." )

  outputFileName =paste(filePath, fileName, sep="/" )
  
  return (getFittedResult(filePath, seerFilePrefix, yearOfDiagnosisVarName, yearOfDiagnosisRange, allVars, cohortVars, cohortValues, covariateVars, numJP,adanced_options, outputFileName))
 # getFullDataDownloadWrapper(filePath,jpsurvDataString)
  
}





#Creates the SEER Data and Fitted Result
getFittedResult <- function (filePath, seerFilePrefix, yearOfDiagnosisVarName, yearOfDiagnosisRange, allVars, cohortVars, cohortValues, covariateVars, numJP, adanced_options,outputFileName) {

  file=paste(filePath, seerFilePrefix, sep="/" )
  
  varLabels=getCorrectFormat(allVars)
  seerdata = joinpoint.seerdata(seerfilename=file,
                                newvarnames=varLabels,
                                NoFit=T,
                                UseVarLabelsInData=varLabels)

  subsetStr=getSubsetStr(yearOfDiagnosisVarName, yearOfDiagnosisRange, cohortVars, cohortValues)
  #assign subsetStr in the global in order for eval(parse(text=)) to work
  assign("subsetStr", subsetStr, envir = .GlobalEnv)

  factorStr=getFactorStr(covariateVars)
  assign("factorStr", factorStr, envir= .GlobalEnv)
  fittedResult=joinpoint(seerdata,
                       subset = eval(parse(text=subsetStr)),
                       year=getCorrectFormat(yearOfDiagnosisVarName),
                       observedrelsurv="Relative_Survival_Cum",
                       model.form = eval(parse(text=factorStr)),
                       op=adanced_options,
                       maxnum.jp=numJP);

  #save seerdata and fit.result as RData
  cat("***outputFileName")
  cat(outputFileName)
  #cat("\n")
  print(outputFileName)
  outputData=list("seerdata"=seerdata, "fittedResult"=fittedResult)
  
  saveRDS(outputData, outputFileName)
  downloadFile=paste("link-", jpsurvData$tokenId,".csv", sep="")
 # getFullDataDownload(outputData, subsetStr, downloadFile, 0) 
  return (outputData)

}
getFullDataDownload <- function(outputData, subsetStr, downloadFile, jpInd) {
  downloadOutput = output.overview(outputData$seerdata, outputData$fittedResult$FitList[[jpInd+1]], subsetStr);
  write.csv(downloadOutput, downloadFile)
}


getRelativeSurvivalByYearWrapper <- function (filePath,filename,jpsurvDataString,jpInd) {
  jpsurvData = fromJSON(jpsurvDataString)
  file=paste(filePath, paste("output-", jpsurvData$tokenId,".rds", sep=""), sep="/")
  outputData=readRDS(file)
  intervals=c(5,10)
#  intervals = jpsurvData$plot$form$intervals
  covariateValues = c("Localized", "Distant")
#  covariateValues = jpsurvData$plot$form$covariateVars

  #take the nth from FitList
  fit.result=outputData$FitList[jpInd+1]
#  png(filename = paste(filePath, paste("plot_Year-", jpsurvData$tokenId, "-",jpsurvData$plot$static[[1]], ".png", sep=""), sep="/"))
#  downloadFile = paste(filePath, paste("data_Year-", jpsurvData$tokenId, "-",jpsurvData$plot$static[[1]], ".csv", sep=""), sep="/")
#  survData=plot.relsurv.year(outputData,intervals,c(NA, NA, NA),covariateValues)
  dev.off()
  return(graphData)
  
}
#Graphs the Survival vs Time graph and saves a csv file of the data
getRelativeSurvivalByIntWrapper <- function (filePath,jpsurvDataString) {
  

  jpsurvDataPath= paste(filePath, jpsurvDataString, sep="/" )
  jpsurvData=fromJSON(jpsurvDataString)
  jpInd=0
  # jpind=jpsurvData$calculate$form$jpInd #<-----new
  
  fileName=paste("output-", jpsurvData$tokenId,".rds", sep="")
  file=paste(filePath, fileName, sep="/" )
  outputData=readRDS(file)
  yearOfDiagnosisVarName = jpsurvData$calculate$static$yearOfDiagnosisVarName
  yearOfDiagnosis = jpsurvData$calculate$form$yearOfDiagnosisRange[[1]]
  
  downloadFile = paste(filePath, paste("data_Int-", jpsurvData$tokenId, "-",jpsurvData$plot$static[[1]], ".csv", sep=""), sep="/") #CSV file to download
  graphFile=(filename = paste(filePath, paste("plot_Int-", jpsurvData$tokenId, "-",jpsurvData$plot$static[[1]], ".png", sep=""), sep="/")) #png file for graph
  
  survData=plot.relsurv.int(outputData$fittedResult$FitList[[jpInd+1]], yearOfDiagnosisVarName, yearOfDiagnosis);
#  write.csv(survData, downloadFile) #<----need to fix this
  
  dev.off()
  survDataJSON=paste(toJSON(survData))
  jsonl =c(survDataJSON,graphFile,downloadFile) #returns 
  
return (jsonl)
}





# getFullDataDownloadWrapper <- function (filePath, jpsurvDataString) {
#   #print("R: getDownloadOutputWrapper")
#   jpsurvData = fromJSON(jpsurvDataString)
#   #print(jpsurvData)
# 
#   #seerFilePrefix = jpsurvData$calculate$static$seerFilePrefix
#   yearOfDiagnosisVarName = jpsurvData$calculate$static$yearOfDiagnosisVarName
#   yearOfDiagnosisRange = jpsurvData$calculate$form$yearOfDiagnosisRange
#   #allVars=jpsurvData$calculate$static$allVars
#   cohortVars=jpsurvData$calculate$form$cohortVars
#   cohortValues=jpsurvData$calculate$form$cohortValues
#   #covariateVars=jpsurvData$calculate$form$covariateVars
#   #maxJP=jpsurvData$calculate$form$joinPoints
#   subsetStr=getSubsetStr(yearOfDiagnosisVarName, yearOfDiagnosisRange, cohortVars, cohortValues)
#   assign("subsetStr", subsetStr, envir = .GlobalEnv)
# 
#   fileName = paste('link', jpsurvData$tokenId, sep="-" )
#   fileName = paste(fileName, "rds", sep="." )
#  
#   #outputFileName = fileName
# 
#   fileName = paste('link', jpsurvData$tokenId, sep="-" )
#   fileName = paste(fileName, "rds", sep="." )
#   outputFileName =paste(filePath, fileName, sep="/" )
# 
#   outFile=paste(filePath, paste("output-", jpsurvData$tokenId,".rds", sep=""), sep="/")
#   downloadFile=paste("link-", jpsurvData$tokenId,".csv", sep="")
# 
#   downloadFile=paste(filePath, paste(filePath, downloadFile, sep="/" ), sep="/")
# 
#   return (downloadFile)
# }



getModelSelectionInfo<-function(filepath,jpsurvDataString) {
  outputData=readRDS(outFile)
  fitList=outputData$fit.result$FitList
  
}
#Gets the coefficients table in the Model Estimates tab
getcoefficientsWrapper <- function (filePath,jpsurvDataString) {
  jpsurvDataPath= paste(filePath, jpsurvDataString, sep="/" )
  jpsurvData=fromJSON(jpsurvDataString)
  fileName=paste("output-", jpsurvData$tokenId,".rds", sep="")
  jpInd=0
  # jpind=jpsurvData$calculate$form$jpInd #<-----new
  file=paste(filePath, fileName, sep="/" )
  outputData=readRDS(file)
  print (file)
  coefficientsJson=paste(toJSON(outputData$fittedResult$coefficients))
  return (coefficientsJson)
}

#gets all the model selection info for all joint points
geALLtModelWrapper <- function (filePath,jpsurvDataString) {
  jpsurvDataPath= paste(filePath, jpsurvDataString, sep="/" )
  jpsurvData=fromJSON(jpsurvDataString)
  fileName=paste("output-", jpsurvData$tokenId,".rds", sep="")
  jpInd=0
  # jpind=jpsurvData$calculate$form$jpInd #<-----new
  file=paste(filePath, fileName, sep="/" )
  outputData=readRDS(file)
  jsonl=c()
  saved=outputData$fittedResult$FitList
  
  for(i in 1:length(saved)) 
  {
    aicJson=paste(toJSON(saved[[i]]$aic))
    bicJson=paste(toJSON(saved[[i]]$bic))
    llJson=paste(toJSON(saved[[i]]$ll))
    convergedJson=paste(toJSON(saved[[i]]$converged))
    jsonl =c(jsonl, aicJson, bicJson, llJson, convergedJson)
    
  }
  return (jsonl)
}

#gets all the model selection info for all joint points
getJointtModelWrapper <- function (filePath,jpsurvDataString) {
  jpsurvDataPath= paste(filePath, jpsurvDataString, sep="/" )
  jpsurvData=fromJSON(jpsurvDataString)
  fileName=paste("output-", jpsurvData$tokenId,".rds", sep="")
  jpInd=0
 # jpind=jpsurvData$calculate$form$jpInd #<-----new
  file=paste(filePath, fileName, sep="/" )
  outputData=readRDS(file)
  jsonl=c()
  saved=outputData$fittedResult$FitList
  aicJson=paste(toJSON(saved[[jpInd+1]]$aic))
  bicJson=paste(toJSON(saved[[jpInd+1]]$bic))
  llJson=paste(toJSON(saved[[jpInd+1]]$ll))
  convergedJson=paste(toJSON(saved[[jpInd+1]]$converged))
  jsonl =c(aicJson, bicJson, llJson, convergedJson)
    
  return (jsonl)
}


getTrendWrapper<- function (filePath,jpsurvDataString,trend_type) {
  
  jpsurvDataPath= paste(filePath, jpsurvDataString, sep="/" )
  jpsurvData=fromJSON(jpsurvDataString)
  fileName=paste("output-", jpsurvData$tokenId,".rds", sep="")
  jpInd=0
  # jpind=jpsurvData$calculate$form$jpInd #<-----new
  file=paste(filePath, fileName, sep="/" )
  outputData=readRDS(file)
  jpInd=0
  trend_type="CS_AAPC"
  # jpind=jpsurvData$calculate$form$jpInd #<-----new
  # trend_type=CS_AAPC, HAZ_APC#<-----new
  
  file=paste(filePath, fileName, sep="/" )
  outputData=readRDS(file)
  trend=toJSON(aapc(outputData$fittedResult$FitList[[jpInd+1]],type=trend_type))
  return(trend)
  
}
