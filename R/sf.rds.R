#' Save variable as .rds file
#'
#' This function saves a variable as an .rds file in a specified directory.
#'
#' @param variable The R object to be saved.
#' @param file Character string specifying the directory path where the file will be saved.
#'
#' @return A message indicating the file name and directory where the .rds file has been saved.
#' @export
#'
#'
sf.rds <- function(variable,file = "viarablename"){
  dir.create(file)
  namefile <- gsub(" ","",paste0("CWAS_result",file,Sys.time(),".rds"))
  namefile <- gsub(":","",namefile)
  #print(variable)
  saveRDS(variable,file = paste0(file,"/",namefile),compress = F)
  output <- paste0("CWAS_result ",namefile,"Has been saved in",file,"(Folder)")
  return(output)
}


