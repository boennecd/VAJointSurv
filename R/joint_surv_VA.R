joint_ms_ptr <- function(markers){
  if(inherits(markers, "marker_term"))
    markers <- list(markers)
  else
    stopifnot(all(sapply(markers, inherits, "marker_term")))

  # compute starting values
  stop("todo")

  .joint_ms_ptr(markers)
}
