#' @title Initializing and Preparing for Parallel Processing
#'
#' @description parallelInitializer is an internal function used by
#' calculatePhysioMap. It is called when user wants to
#' run calculatePhysioMap in parallel.
#'
#' @param NumbrOfCores Number of cores (Threads) to use in parallel.
# #' @param CLUSTER logical value, setting if the calculation should be
# #' done on a Cluster. Not implemented in this version of PhysioSpaceMethods.
# #' @param machineAddresses Address of Cluster nodes, only relevant when
# #' CLUSTER==TRUE. Not implemented in this version of PhysioSpaceMethods.
# #' @param rscript Location of rscript on each node of Cluster, only
# #' relevant when CLUSTER==TRUE. Not implemented in this version of
# #' PhysioSpaceMethods.
#'
#' @return An object of class c("SOCKcluster", "cluster"), made by
#' makeCluster of package parallel.
#'
#' @examples parallelInitializer(NumbrOfCores = 2)
#'
#' @export
parallelInitializer <- function(NumbrOfCores){
    # suppressPackageStartupMessages(require(doParallel))
    # if(CLUSTER) {
    #   cl <- parallelInitializer.Cluster(machineAddresses, rscript)
    # } else {
    cl <- parallelInitializerSolo(NumbrOfCores)
    registerDoParallel(cl)
    return(cl)
}

# #' @export parallelInitializer.Cluster
#
# parallelInitializer.Cluster <- function(machineAddresses, rscript){
#   ##My ip is dynamic I have to get it everytime!:
#   primary <- system("dig +short myip.opendns.com @resolver1.opendns.com",
#                                                                intern = T)
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

#' @title Initializing and Preparing for Parallel Processing on a Single Machine
#'
#' @description parallelInitializerSolo is an internal
#' function of parallelInitializer.
#'
#' @param NumbrOfCores Number of cores (Threads) to use in parallel.
#'
#' @return An object of class c("SOCKcluster", "cluster"), made by
#' makeCluster of package parallel.
#'
#' @examples parallelInitializerSolo(NumbrOfCores = 2)
#'
#' @export
parallelInitializerSolo <- function(NumbrOfCores){
    if(is.na(NumbrOfCores)) {
        NumbrOfCores=detectCores()-1
        #not to overload your computer
        if(NumbrOfCores==0) NumbrOfCores <- 1
        #In case the computer is single core
    }
    if(NumbrOfCores > detectCores())
        stop("'NumbrOfCores' can not be higher than",
                "the available number of cores, which is",
                as.character(detectCores())
                )
    cl <- makeCluster(NumbrOfCores, outfile="")
    #outfile is for having progressbar working in foreach
    return(cl)
}

