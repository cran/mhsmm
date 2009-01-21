.First.lib <- function(lib,pkg)
{
	#	require(mvtnorm)
   library.dynam("mhsmm",pkg,lib)
   cat("hsmm shared library loaded\n")
}

