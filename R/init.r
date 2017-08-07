

.onLoad <-function(libname,pkgname) {
  #print( "Calling .onLoad in genepathway");
  library.dynam("SNPath", pkgname, libname)
}