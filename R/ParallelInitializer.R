#' @title Creates and Uses Physio Spaces as a dimension reduction mapping
#'
#' @description This package uses 'Big Data' to make robust 'Physiological Vectors' in N dimensions spaces, with which you will map new data to extract infromation from a big high dimensional confusing new data.
#'
#' @param CLUSTER,machineAddresses,rscript,NumbrOfCores
#'
#' @return NULL
#'
#' @examples parallelInitializer()
#'
#' @export parallelInitializer

parallelInitializer <- function(NumbrOfCores){
  # suppressPackageStartupMessages(require(doParallel)) #for now this is move to'calculatePhysioMap.R'
  # if(CLUSTER) {
  #   cl <- parallelInitializer.Cluster(machineAddresses, rscript)
  # } else {
    cl <- parallelInitializer.SingleMachine(NumbrOfCores)
  registerDoParallel(cl)
  return(cl)
}

# #' @export parallelInitializer.Cluster
#
# parallelInitializer.Cluster <- function(machineAddresses, rscript){
#   ##My ip is dynamic I have to get it everytime!:
#   primary <- system("dig +short myip.opendns.com @resolver1.opendns.com", intern = T)
#
#   spec <- lapply(machineAddresses,
#                  function(machine) {
#                    rep(list(list(host=machine$host,
#                                  user=machine$user)),
#                        machine$ncore)
#                  })
#   spec <- unlist(spec,recursive=FALSE)
#   cl <- makeCluster(type='PSOCK',
#                               master=primary,
#                               spec=spec,rscript=rscript, outfile="")
#   return(cl)
# }

#' @export parallelInitializer.SingleMachine
#
parallelInitializer.SingleMachine <- function(NumbrOfCores){
  if(is.na(NumbrOfCores)) NumbrOfCores=detectCores()-1 # not to overload your computer
  cl <- makeCluster(NumbrOfCores, outfile="") # outfile is for having progressbar working in foreach
  return(cl)
}

