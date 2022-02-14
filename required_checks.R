## Required checks


# 0. User input parameters correct
# correct algorithm names (which are supported by us)
# maybe warning for very large iterations settings
# segments >= 1

# 1. Check if NAs are present
# I actually think we should allow inputs without NAs - users also want to compare methods on complete data
if (!anyNA(x)) {
  return(x)
}

# 2. Check if input is numeric
if (!is.numeric(x)) {
  stop("Input x is not numeric")
}

# 3. Check if input is univariate
# what to do with multivariate data

# 4. Check if input is zoo / ts object

# 5. Check if a whole segment is all NA

# Minimum length segment (100 Segments with 1 observations each wouldn't make sense)

# 6. Check if imputation algortihm worked correctly
# Try-catch block around imputation call - since we can't check/mitigate all conditions a package/algorithms
# might not provide complete imputation
# starting imputation has to work though - maybe fallback to default
