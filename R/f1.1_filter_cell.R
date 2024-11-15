filter_cell <- function(sparkle.data){

  celltypeall <- unique(sparkle.data$Celltype)

  filter_cellname <- list()
  for(i in 1:length(celltypeall)){

   df <-   sparkle.data%>% dplyr::select_all() %>% dplyr::filter(Celltype==celltypeall[i])
   num <- length(unique(df$Phenotype))

   if(num<2){filter_cellname[[i]] <- celltypeall[i]}

  }

  filter_cell <- unlist(filter_cellname)

  print(paste(filter_cell,"is filter due to the incomplete phenotype information"))

  sparkle.data.new <- sparkle.data %>% dplyr::select_all() %>% dplyr::filter(!Celltype%in%filter_cell)

  return(sparkle.data.new)
}
