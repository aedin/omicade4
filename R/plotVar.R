plotVar <-
function(x, var=NA, axes=1:2, var.col="red",
         var.lab=FALSE, bg.var.col="gray", 
         nlab=0, 
         sepID.data=NULL, # NULL, 1, 2, 1:2
         sepID.sep="_", ...) {
  UseMethod("plotVar")
}
