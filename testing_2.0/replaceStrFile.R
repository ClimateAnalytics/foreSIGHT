#####################################

replaceStrFile = function(fname,oldStr,newStr){

  x <- readLines(fname)
  x <- gsub(oldStr,newStr,x)
  cat(x, file=fname, sep="\n")
}

#####################################

