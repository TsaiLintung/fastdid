silence = function(expr) {
  #from : https://stackoverflow.com/questions/6177629/how-to-silence-the-output-from-this-r-package
  
  #temp file
  f = file()
  
  #write output to that file
  sink(file = f)
  
  #evaluate expr in original environment
  y = eval(expr, envir = parent.frame())
  
  #close sink
  sink()
  
  #get rid of file
  close(f)
  
  y
}