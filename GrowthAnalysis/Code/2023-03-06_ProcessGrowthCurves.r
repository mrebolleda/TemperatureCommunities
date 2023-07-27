## Jackie Folmar and Maria Rebolleda-Gomez
## April 6th 2023
## Functions to process growth curves


######RESHAPE FUNCTIONS#####
shape_gc_batch <- function (file, n.skip, cols, string_to_replace = " ", r = NULL) {
  y <- fread(file, skip = n.skip, blank.lines.skip = TRUE, fill = TRUE) # read in data
  
  if (!is.null(r)){
     y <- y[r,] %>% select_if(~ !any(is.na(.))) # remove all NAs
  } else {
    y <- y %>% select_if(~ !any(is.na(.)))
  }

  # Change to be more generic in size if necessary
  y[,2:ncol(y)] <- lapply(y[,2:ncol(y)],as.numeric) # make class numeric 
  y <- y %>% melt(measure.vars = cols, variable.name = "well", value.name = "abs")
 
  #Change to find the plate number that way can work with any code
  y[,plate:=substring(file,nchar(file)-5,nchar(file)-4)] # CHANGE TO VARY WITH FILE finds plate number, changed from last time
  vec <- gsub(string_to_replace, "", y$plate) # string_to_replace = " " by default
  y$plate <- as.numeric(vec) # make class numeric
  y
}


shape_gc <- function(file, nskip, c, r, t = "Time [s]"){
  y <- fread(file, skip = nskip, blank.lines.skip = TRUE, fill = TRUE) %>%
  select(-c)

  # select rows without artifact information
  y <- y[r,] 

  # since artifact information was str, create as.numeric list
  y[,2:ncol(y)] <- lapply(y[,2:ncol(y)],as.numeric) # make class numeric 
  y <- melt(data.table(y),id.vars = t, variable.name='well', value.name='abs')
  y
}
