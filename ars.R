## STAT 243 Final Project: adaptive rejection sampling ##
## Jinze Gu, Haavard Kvamme, Reid Stevens and Yang Wang ##

## This is the primary function that performs the sampling from a log-concave probability density function using the adaptive rejection sampling algorithm developed by Gilks and Wild (Appl. Statist. 41, 337, 1992) ##

#f is the target density function, n is the number of required sample, test indicates whether to perform the test, k is the initial number of abscissas

ars <- function(f, n, left_bound = -Inf, right_bound = Inf, x_init, ...) {
  #Generating the initial abscissas x
  x <- x_init
  x <- sort(x)
  hx <- log(f(x, ...))
  hpx <- diag(attributes(numericDeriv(quote(log(f(x, ...))), "x"))$gradient)
  if (((hpx[1] < 0) & (left_bound==-Inf)) | ((hpx[length(hpx)]>0) & (right_bound==Inf))) stop("The derivative at the first/last initial point must be positive/negative")
  sample <- rep(NA, n)
  count <- 0
  z <- make_z(x, hx, hpx, left_bound, right_bound) 

  while (count < n) {
    #Make the lower_bound, upper bound functions
    lower_bound <- make_lower_bound(x, hx)
    upper_bound <- make_upper_bound(x, hx, hpx, z)
    
    #Draw samples from the upper bound function
    cand <- sample_upper_bound(n - count, x, hx, hpx, z)
    
    #"update" records the first point in the sample that can not be accepted based on lowerbound
    #"accepted" records the candidates drawn from upper_bound function that are accepted by only using the lowerbound until the "update" point
    cand_filtered <- filter(cand, lower_bound, upper_bound, n-count)
    accepted <- cand_filtered$accepted
    update <- cand_filtered$update
    
    if (is.na(update[1])) {
      sample[(count+1):n] <- cand
      return(sample)
    }
    
    #Check the log-concaveness by checking if the log(f(update)) falls in between lower and upper bound function 
    if ((log(f(update$cand, ...)) < lower_bound(update$cand)) | (log(f(update$cand, ...)) > upper_bound(update$cand))) stop("The sample function is not log-concaved!")
    
    #Update the sample using cand_filtered
    index <- update_sample(cand, accepted, update, f, upper_bound)
    if (index != 0) {
      sample[(count+1):(count+index)] <- cand[1:index]
      count <- count + index
    }
    #Update the abscissas x
    update_absci <- update_x(f, x, hx, hpx, update$cand, ...)
    x <- update_absci$x
    hx <- update_absci$hx
    hpx <- update_absci$hpx
    z <- make_z(x, hx, hpx, left_bound, right_bound)
  }
  return (sample)
}

# -----------------------------------------------------------------------
# The function filter() used in adaptive_rejection_sampling_main.R
# It performs a squeezing test between the lower and upper bound, on the candidates in cand.

# Returns a list with
# $accepted: integer: position of last accepted candidate
# $update:   list: $cand is fist value not accepted
# 		   $w is corresponding draw form uniform dist
#                  Or if all candidates are accepted returns NA
# If all the cand are accepted, the update value is NA
# If no cand is accepted, the accepted value is NA

filter <- function(cand, lower_bound, upper_bound, m){
  w <- runif(m) 
  squeeze <- lower_bound(cand) - upper_bound(cand)

  # Gives NA if all are accepted
  i <- which(log(w) > squeeze)[1]

  if (is.na(i)) {
    return(list(accepted=m, update=NA))
  } 
  else{
    if (i==1) {
      return(list(accepted=NA, update=list(cand=cand[i], w=w[i])))
    }
    else {
      return(list(accepted=(i-1), update=list(cand=cand[i], w=w[i])))
    }
  }
}

#----------------------------------------------------------------------------
# update_sample 
# Tests if the first candidate rejected by filer should be accepted after all
# Returns the index of the last accepted candidate.
update_sample <- function(cand, accepted, update, f, upper_bound, ...) {
  flag <- (log(update$w) < log(f(update$cand, ...))/upper_bound(update$cand))
  if (is.na(accepted)) {
    if (flag){
      return(1)
    }
  }
  else {
    if (flag) {
      return(accepted + 1)
    }
    else {
      return(accepted)
    }
  }
  return(0)
}


# --------------------------------------------------------------
#make_lower_bound function returns the lower_bound function
#lower_bound function is set to be -Inf when argument falls outside the abscissas range.

make_lower_bound <- function(x, hx) {
  min <- min(x)
  max <- max(x)
  left <- min - abs(min*.Machine$double.eps)
  right <- max + abs(max*.Machine$double.eps)
  x <- c(left, x, right)
  hx <- c(-Inf, hx, -Inf)
  lower_bound <- approxfun(x, hx, rule = c(2,2))
  return(lower_bound)
}

# ------------------------------------------------------------------
make_upper_bound <- function(x, hx, hpx, z) {
  upper_bound_eval_z <- function(z) {
    index <- rowSums(outer(z, z, function(x1,x2) x1>x2))
    index[index==0] <- 1
    index[index==length(z)] <- length(z) - 1
    ux <- hx[index] + (z-x[index])*hpx[index]
    return (ux)
  }
  # Because of numerical issues
  # Need to find a better way to do to this. Is x sorted? 
  m <- length(z)
  if (z[1] == -Inf) z[1] <- -10*max(abs(x))
  if (z[m] == Inf) z[m] <- 10*max(abs(x))
  y <- upper_bound_eval_z(z)
  upper_bound <- approxfun(z,y)
  return(upper_bound)
}


# ---------------------------------------------------------------------------
# Calculating the z based on abscissae x. z[1]=left_bound, z[k+1]=right_bound
# Returns values for z
make_z <- function(x, fx, fpx, left_bound, right_bound) {
  k <- length(x)
  xm1 <- c(0, x)
  fxm1 <- c(0, fx)
  fpxm1 <- c(0, fpx)
  x <- c(x, 0)
  fx <- c(fx, 0)
  fpx <- c(fpx, 0)
  
  z <- (fx - fxm1 - x*fpx + xm1*fpxm1) / (fpxm1 - fpx)
  z[1] <- left_bound
  z[k+1] <- right_bound
  
  return (z)
}


# ----------------------------------------------------------------------------
# sample_upper_bound function samples m points from the upper bound funtion using inverse cdf method. The inverse cdf is calculated analytically and implemented in the auxillary function inversecdf.
sample_upper_bound <- function(m, x, hx, hpx, z) {
  k <- length(x)
  
  #Calcualte the normalized integration of upper_bound funtion at each linear interval
  factor1 <- exp(hx - x*hpx) 
  factor2 <- (exp(c(0, hpx) * z)[2:(k+1)] - exp(hpx*z[1:k])) / hpx
  I <- factor1 * factor2
  Inormalize <- sum(I)
  
  intervalcdf <- cumsum(I)
  sample <- runif(m, min=0, max=Inormalize)
  #calculate which interval the sample falls into
  #sample_interval <- sapply(sample, function(t) sum(intervalcdf<t)) + 1
  sample_interval <- rowSums(outer(sample, intervalcdf, function(x,y) x>y)) + 1
  t <- sample - c(0, intervalcdf)[sample_interval]
  
  sample_x <- inversecdf(t, sample_interval, factor1, hpx, z)
  return (sample_x)
}


inversecdf <- function(t, j, factor1, hpx, z) {
  return ( log(t/factor1[j]*hpx[j] + exp(hpx[j]*z[j])) / hpx[j] )
}

# ----------------------------------------------------------------------------
#update_x function adds the update point into the abscissas vector, and also checks if the derivative at the update point is between its neighbors(i.e., checking the concaveness)

update_x <- function(f, x, hx, hpx, update, ...) {
  if (is.na(update)){
    return (list(x=x, hx=hx, hpx=hpx))
  }
  else {
    hx_update <- log(f(update, ...))
    xx <- update
    hpx_update <- attributes(numericDeriv(quote(log(f(xx, ...))), "xx"))$gradient[1, 1]
    rank_x <- rank(c(x, update))
    new_x <- rep(NA, length(x)+1)
    new_hx <- rep(NA, length(x)+1)
    new_hpx <- rep(NA, length(x)+1)
    
    new_x[rank_x] <- c(x, update)
    new_hx[rank_x] <- c(hx, hx_update)
    new_hpx[rank_x] <- c(hpx, hpx_update)
    
    #Check concaveness by testing if h'(x) is decreasing
    if (is.unsorted(-new_hpx)) stop("The sample function is not log-concaved!")
    
    return(list(x=new_x, hx=new_hx, hpx=new_hpx))
  }
}

#############################################################################
# ARS function test
set.seed(0)
hist(ars(dnorm, 10000, x_init = c(-1,1)), main = "ARS sampling from a truncated normal distribution", freq = 0)

hist(ars(dnorm, 10000, x_init = c(-1,1)), main = "ARS sampling from a truncated normal distribution", freq = 0)

cat("The following does not satisfy our input requirements. \n")
hist(ars(dnorm, 2000, x_init = runif(2)), main = "ARS sampling from a normal distribution")

# qqplot
ars_sample <- ars(dnorm, 10000, x_init=c(-1,1))
qqnorm(ars_sample)
qqline(ars_sample)

#Test for extra arguments
ars_sample <- ars(dnorm, 10000, x_init=c(1,3), mean = 2, sd = 5)
qqnorm(ars_sample)
qqline(ars_sample)

#Test for bounded distribution
ars_sample <- ars(dnorm, 10000, left_bound=-2, right_bound=2,x_init=c(-1,1))
hist(ars_sample)

#To check if it can catch bad initial points
ars_sample <- ars(dnorm, 10000, x_init=c(1,2))
ars_sample <- ars(dnorm, 10000, left_bound=0.1, x_init=c(1,2))
hist(ars_sample)

#############################################################################
# Make_z function test

set.seed(0)
x <- sort(rnorm(100))
fx <- dnorm(x, log = 1)
fpx <- diag(attributes(numericDeriv(quote(dnorm(x, log = 1)), "x"))$gradient)
plot(x, fx, type = "p", pch = 16)
points(make_z(x, fx, fpx, -12, 12), dnorm(make_z(x, fx, fpx, -12, 12),log = 1), pch = 17, col = "red")
legend("topright", c("log-normal plot", "make_z generation plot"), pch = c(16, 17))


cat("Since z is just the point from the density function(in this case, lognormal), so we can just plot two set of points together to see if they fall into the same distribution trend. \n")

#############################################################################
# make_upper/lower bound tests
# Test 1
set.seed(0)
x <- runif(100)
hx <- dnorm(x, log = 1)
hpx <- diag(attributes(numericDeriv(quote(dnorm(x, log = 1)), "x"))$gradient)
z <- make_z(x,hx,hpx,-12.0, 12.0)
lower <- make_lower_bound(x, hx)
upper <- make_upper_bound(x, hx, hpx, z)
xx <- seq(-5, 5, by=0.1)
cat("This shows that all lower bound should less than upper bound by our generation. \n")
all((lower(xx) <= upper(xx)) == TRUE)

#Test 2
cat("Now we use -x^2 as a function to plot upper and lower bound. \n")
x <- seq(-10, 10, by=4)
cat("Notice here we need to define a function and input the name of this function \n")
f<-function(x){-x^2}
fpx <- diag(attributes(numericDeriv(quote(f(x)), "x"))$gradient)
z <- make_z(x,f(x),fpx,-Inf, Inf)
lower <- make_lower_bound(x, f(x))
upper <- make_upper_bound(x, f(x), fpx, z)
xx <- seq(-10, 10, by=0.2)
plot(xx, -xx^2, type='l')
lines(xx, lower(xx), lty = 2)
lines(xx, upper(xx), lty = 6)
legend("topright", c("upper_bound", "lower_bound"), lty = c(6,2))


cat("In test 1, we basically just check if the corresponding lower bound is less than upper bound, and the function passes the test as long as the output is TRUE. Then we tried to plot our upper/lower bound using previous function we defined: $make_z$.We use a logconcave function $f(x)=x^2$, and it turns out the plot is what we had expected. \n")

#############################################################################
# sample_upper_bound test

set.seed(0)
x <- sort(rnorm(10000))
cat("Notice here we need to define a function and input the name of this function \n")
f<-function(x){-x^2}
fpx <- diag(attributes(numericDeriv(quote(f(x)), "x"))$gradient)
z <- make_z(x,f(x),fpx,-Inf, Inf)
upper <- make_upper_bound(x, f(x), fpx, z)
xx <- seq(-10, 10, by=0.2)
hist(sample_upper_bound(2000, x, f(x), fpx, z), freq = 0)

cat("Since the histogram of sample_upper_bound function matchs the trend of $f(x) = -x^2$, so it is obvious that our sample_upper_bound function works as we expected. \n")

#############################################################################
# filter function test
set.seed(0)
f <- dnorm
n <- 10
x <- sort(rnorm(7))
left_bound   <- -10
right_bound  <- 10
hx <- log(f(x))
hpx <- diag(attributes(numericDeriv(quote(log(f(x))), "x"))$gradient)
sample <- rep(NA, n)
count <- 0
z <- make_z(x, hx, hpx, left_bound, right_bound) 
lower_bound <- make_lower_bound(x, hx)
upper_bound <- make_upper_bound(x, hx, hpx, z)
cand <- sample_upper_bound(n - count, x, hx, hpx, z)
cand_filtered <- filter(cand, lower_bound, upper_bound, n - count)

cat("Basically, we have our candidates and cand_filtered, and the function does what it is suppose to do. \n")
cand
cand_filtered

#############################################################################
# update_sample test

set.seed(0)
f<-dnorm
x <- sort(rnorm(7))
hx <- log(f(x))
hpx <- diag(attributes(numericDeriv(quote(log(f(x))), "x"))$gradient)
sample <- c(1,2,3, rep(NA, 10))
cat("When accepted point is not 'NA' \n")
accepted <- cand_filtered$accepted
update <- cand_filtered$update
count <- length(na.omit(sample))
upper_bound <- make_upper_bound(x, hx, hpx, z)
update_sample(cand = cand, accepted=accepted, update = update, f = f, upper_bound = upper_bound)
cat("When accepted point is 'NA' \n")
accepted <- NA
update_sample(cand, accepted, update,  f, upper_bound)

cat("In this case, the update_sample_function returns the index of the last accepted point, and it returns 1 if the point satisfies the test flag and we accepted no point. It returns 0 if the point does not satisfy flag and we accepted no point. According to the test, it works correctly with specific output. \n")

#############################################################################
# update_x test

set.seed(0)
x<-seq(0.0001,1, by = 0.2)
f<-dnorm
hx <- log(f(x))
hpx <- diag(attributes(numericDeriv(quote(log(f(x))), "x"))$gradient)
update <- cand_filtered$update
upper_bound <- make_upper_bound(x, hx, hpx, z)
update_x(f = f, x = x, hx = hx, hpx = hpx, update = update$cand)

cat("The update_x outputs the updated points, their evaluations under function h(x) and their corresponding derivatives \n")


#############################################################################
