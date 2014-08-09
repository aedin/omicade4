.checklim <-
function(axis.lim) {
  if (!is.numeric(axis.lim)) 
    stop("only numeric allowed for defining plotted region")
  if (length(axis.lim) != 2)
    stop("only two numbers needed in axis.lim, first for the bottom boundary, 
           second for upper boundary ")
  if (axis.lim[1] >= axis.lim[2])
    stop("axis.lim[1] should be lower than axis.lim[2]")
}
