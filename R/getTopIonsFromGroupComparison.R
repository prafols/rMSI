#' getTopIonsFromGroupComparison
#'
#' @param groupsComparison
#' @param massVector 
#' @param thresh_p 
#' @param thresh_FC 
#' @param maxNumberOfTopIons 
#' @param pAdjust 
#' @param pAdjustMethod 
#'
#' @return
#' @export
#'
#' @examples
getTopIonsFromGroupComparison <- function(groupsComparison, massVector, thresh_p = 0.05, thresh_FC = 1, maxNumberOfTopIons = 10, pAdjust  = F, pAdjustMethod  = "fdr" )
{
  results <- list()
  for(i in 1:length(groupsComparison))
  {
    downIons <- which((groupsComparison[[i]]$pV >= -log10(thresh_p)) & (groupsComparison[[i]]$FC <= -thresh_FC))
    upIons <- which((groupsComparison[[i]]$pV >= -log10(thresh_p)) & (groupsComparison[[i]]$FC >= thresh_FC))
    
    if(length(upIons) > maxNumberOfTopIons)
    {
      maxpVs <- upIons[order(groupsComparison[[i]]$pV[upIons], decreasing = T)[1:maxNumberOfTopIons]]
      maxFCs <- upIons[order(groupsComparison[[i]]$FC[upIons], decreasing = T)[1:maxNumberOfTopIons]]
      maxVolcanoIons <- unique(c(maxpVs, maxFCs))
      if(length(maxVolcanoIons) > maxNumberOfTopIons)
      {
        upIons <- (maxVolcanoIons[order(sqrt(((groupsComparison[[i]]$pV[maxVolcanoIons]/max(groupsComparison[[i]]$pV[maxVolcanoIons]))^2) +
                                               ((groupsComparison[[i]]$FC[maxVolcanoIons]/max(groupsComparison[[i]]$FC[maxVolcanoIons]))^2)), decreasing = T)])[1:maxNumberOfTopIons]
      }
      else 
      {
        upIons <- maxVolcanoIons
      }
    }
    
    if(length(downIons) > maxNumberOfTopIons)
    {    
      maxpVs <- downIons[order(groupsComparison[[i]]$pV[downIons], decreasing = T)[1:maxNumberOfTopIons]]
      maxFCs <- downIons[order(abs(groupsComparison[[i]]$FC[downIons]), decreasing = T)[1:maxNumberOfTopIons]]
      maxVolcanoIons <- unique(maxpVs, maxFCs)
      if(length(maxVolcanoIons) > maxNumberOfTopIons)
      {
        downIons <- (maxVolcanoIons[order(sqrt(((groupsComparison[[i]]$pV[maxVolcanoIons]/max(groupsComparison[[i]]$pV[maxVolcanoIons]))^2) +
                                                 (abs(groupsComparison[[i]]$FC[maxVolcanoIons]/max(groupsComparison[[i]]$FC[maxVolcanoIons]))^2)), decreasing = T)])[1:maxNumberOfTopIons]
      }
      else 
      {
        downIons <- maxVolcanoIons
      }
    }
    
    upIons <- sort(upIons)
    downIons <- sort(downIons)
    
    results[[names(groupsComparison)[i]]]$Up <- data.frame(ion = upIons,
                                                 mz = massVector[upIons],
                                                 pV = groupsComparison[[i]]$pV[upIons],
                                                 FC = groupsComparison[[i]]$FC[upIons])
    
    
    results[[names(groupsComparison)[i]]]$Down <- data.frame(ion = downIons,
                                                   mz = massVector[downIons],
                                                   pV = groupsComparison[[i]]$pV[downIons],
                                                   FC = groupsComparison[[i]]$FC[downIons])
    
  }
  
  if(pAdjust)
  {
    for(i in 1:length(groupsComparison))
    {
      if((length(results[[names(groupsComparison)[i]]]$Up$pV) > 0) & (length(results[[names(groupsComparison)[i]]]$Down$pV) > 0))
      {
        adjusted.pValues <- p.adjust((10^(-c(results[[names(groupsComparison)[i]]]$Up$pV, results[[names(groupsComparison)[i]]]$Down$pV))), method = pAdjustMethod)
        results[[names(groupsComparison)[i]]]$Up$pV <- -log10(adjusted.pValues[1:length(results[[names(groupsComparison)[i]]]$Up$pV)])
        results[[names(groupsComparison)[i]]]$Down$pV <- -log10(adjusted.pValues[length(results[[names(groupsComparison)[i]]]$Up$pV)+(1:length(results[[names(groupsComparison)[i]]]$Down$pV))])
      }
      else if((length(results[[names(groupsComparison)[i]]]$Up$pV) > 0))
      {
        adjusted.pValues <- p.adjust(10^(-results[[names(groupsComparison)[i]]]$Up$pV), method = pAdjustMethod)
        results[[names(groupsComparison)[i]]]$Up$pV <- -log10(adjusted.pValues)
      }
      else if((length(results[[names(groupsComparison)[i]]]$Down$pV) > 0))
      {
        adjusted.pValues <- p.adjust(10^(-results[[names(groupsComparison)[i]]]$Down$pV), method = pAdjustMethod)
        results[[names(groupsComparison)[i]]]$Down$pV <- -log10(adjusted.pValues)
      }
    }
  }
  
  for(i in 1:length(groupsComparison))
  {
    colnames(results[[names(groupsComparison)[i]]]$Up) <- c("index", "m/z", "-log10(p.value)", "log2(FC)")
    colnames(results[[names(groupsComparison)[i]]]$Down) <- c("index", "m/z", "-log10(p.value)", "log2(FC)")
  }
  
  return(results)
}