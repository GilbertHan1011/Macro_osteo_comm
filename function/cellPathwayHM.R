#== signal pathway by ident-------------------

#' cellPathwayHM - Visualization function with cellchat object
#'
#' visualize pathway sent or recieved by a specific ident
#' 
#' @encoding UTF-8
#' @param object cellchat object
#' 
#' @param ident ident in cellchat
#' 
#' @param role one of 'send' or 'recieve'
#' 
#' @param row_scale Boolean, whether to scale row
#' 
#' @param color_fun color function in ComplexHeatmap
#' 
#' @return a plot by ComplexHeatmap
#'

cellPathwayHM <- function(object,ident,role="send",row_scale=TRUE,
                          color_fun=colorRamp2(c(-1, 0, 3),c("#4DBBD5FF", "white", "#E64B35FF"))){
  require(ComplexHeatmap)
  netProbility <- object@netP$prob
  identList <- rownames(netProbility[1,,])
  identIndex <- which(identList==ident)
  if (length(identIndex) == 0) {
    stop("Error: please provide a true ident in cellchat's ident")
  }
  if (!(role %in% c("send","recieve"))) {
    stop("Error: role must be one of 'send' or 'revieve'")
  }
  # subset pathway sended or recieved by ident 
  message(paste0("Find pathway ",role, " by ",ident))
  if (role=="send"){
    sigMatrix <- apply(netProbility,3,function(x) x[identIndex,])
  }
  if (role=="recieve"){
    sigMatrix <- apply(netProbility,3,function(x) x[identIndex,])
  }
  print(head(sigMatrix))
  
  # calculate column sums
  col_sums <- colSums(sigMatrix)
  
  # subset data frame
  sigMatrix <- sigMatrix[, col_sums != 0]
  
  # scale
  if (row_scale=TRUE){
    sigMatrix = t(scale(t(sigMatrix)))
  }
  
  # plot heatmap
  hm <- Heatmap(scaled_mat,col = color_fun,
          name = "scaled prob",
          column_title  = paste0("Signaling Pathway",role," by ",ident))
  draw(hm)
}
