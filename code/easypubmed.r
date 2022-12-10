filenames_1 <- list.files('./deg/AB transform/',pattern = ".xls",full.name=TRUE)
filenames_1
ab1 <- lapply(filenames_1, function(fl) 
  read.delim(fl))
