# Set of axuxiliary functions for the TP inference.
# Author: Costas Bouyioukos

matrix_threshold <- function(m, t, both = TRUE, ...) {
# A function to apply a threshold value on a matrix m
  diag(m) <- 0;
  m[m >= -t & m <= t] <-0;
  return(m)
}


export_cytoscape <- function(m, filename, ...) {
# Function to export a correlation matrix as a cytoscape network file.
  #if () {
  #  stop("Pleaase specify a new file name there is not much sense in appending a network in an existing filename.")
  #}
  fn <- filename;
  write("from\tto\tedgeWeight\tsign", file = fn);
  for (i in 1:nrow(m)) {
    for (j in i:ncol(m)) {
      if (m[i,j] != 0) {
        write(sprintf("%s\t%s\t%f\t%s", colnames(m)[i], rownames(m)[j], abs(m[i,j]), sign(m[i,j])), file = fn, append = TRUE);
      }
    }
  }
}
