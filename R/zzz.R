".First.lib" <- function(lib, pkg)
{
  library.dynam("Bayesthresh", package = pkg, lib.loc = lib)
  return(invisible(0))
}

