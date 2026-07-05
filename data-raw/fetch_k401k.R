# Fetch and clean the 401(k) dataset from the 'wooldridge' package
if (!requireNamespace("wooldridge", quietly = TRUE)) {
  install.packages("wooldridge", repos = "https://cloud.r-project.org")
}

# Load the dataset
data("k401k", package = "wooldridge")
fracreg_k401k <- k401k

# The participation rate 'prate' is a percentage (0-100). 
# We need to scale it to [0, 1] to use it as a fractional response variable.
fracreg_k401k$prate <- fracreg_k401k$prate / 100

# Calculate squared terms that are often used in empirical analysis
fracreg_k401k$age_sq <- fracreg_k401k$age^2
fracreg_k401k$mrate_sq <- fracreg_k401k$mrate^2

# Save the data to the 'data/' directory
save(fracreg_k401k, file = "data/fracreg_k401k.rda", compress = "xz")
