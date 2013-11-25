#Approximate derivatives:
x = c(0.3,0.4,0.5)
h<-10^(-(1:15))
digamma(x)
test<-function(x){
  return(c(lgamma(x),digamma(x),lbeta(x,x+1)))
}
(lgamma(x+h) - lgamma(x))/h
fpx<-numericDeriv(quote(lgamma(x)),"x")
diag(attr(fpx,"gradient"))

#The input of fx should be a expression, say 'dnorm','lgamma'..
# x_star should come from the S function, initial is the sample used to create the lower_bound and upper_bound 
initial<-seq(0.1,0.5,by = 0.004)
x_star<-seq(0.07,0.37, by = 0.003)
#Check length(x_star) = length(initial), becasue I have output h_prime_x_star as part of my data.frame
#Otherwise, it is ok to have different length which is general in most of cases
#h_f represents h(x), h_prime represents 1st deriv of h(x), T_k and initial all represent the
#points that have already been evaluated as the assumption in paper.
make_z<-function(initial,h_f = dlnorm, left_bound = -Inf,right_bound = Inf){
  T_k = sort(initial)
  h_prime<-diag(attr(numericDeriv(quote(h_f(T_k)),"T_k"),"gradient"))
  if(all(diff(h_prime) <= 0)){# Add sth to test differentiability)
    k<-length(T_k)
    #Need some test here for the bounds
    z_mid<-(h_f(T_k[2:k])-h_f(T_k[1:(k-1)])-T_k[2:k]*h_prime[2:k]+T_k[1:(k-1)]*h_prime[1:(k-1)]) / (h_prime[1:(k-1)]-h_prime[2:k])
    z<-c(left_bound, z_mid, right_bound)
    return(z)
  }
  else {print("The log-concativity is not good in this interval")}
}
make_z(initial, h_f = dnorm)

make_upper_bound<-function(x, initial = initial,h_f = dlnorm,left_bound = -Inf,right_bound = Inf){
  x_star<-sort(x); length_x<-length(x_star); length_initial<-length(initial)
  T_k = sort(initial)
  z<-make_z(initial, h_f, left_bound = left_bound, right_bound = right_bound); length_z<-length(z)
  h_prime_x_star<-diag(attr(numericDeriv(quote(h_f(x_star)),"x_star"),"gradient"))
  h_prime_initial<-diag(attr(numericDeriv(quote(h_f(T_k)),"T_k"),"gradient"))
  u_z <- h_f(initial[1:(length_initial-1)]) + (z[2:(length_z -1)] - initial[1:length_initial-1])*h_prime_initial[1:length_initial - 1]
  u_z <- c(u_z,NA)
  upper_matrix<-data.frame(lower_limit = z[1:(length_z - 1)], upper_limit = z[2:length_z], h_prime_x_star = h_prime_x_star, 
                           h_f_initial = h_f(T_k), h_prime_initial = h_prime_initial,upper_hull_zvalue = u_z)
  # Add this part after sunday discussion to evaluate x_star in upper bound
  y_upper <- vector("numeric",length_x)
  j_store = 1
  for(i in 1:length_x){
    for(j in j_store:(length_initial-1)){
      if(x_star[i]>=upper_matrix[j,1] & x_star[i]<=upper_matrix[j,2]){
        #T_k[j]s are the x[j]s in the paper, which are the initial points.
        y_upper[i] <- h_f(T_k)[j] + (x_star[i] - T_k[j])*h_prime_initial[j]
        j_store <- j# Since we have sorted x_star, it saves some computation
      }
    }
  }
  return(data.frame(upper_matrix,upper_value = y_upper))
}
a<-make_upper_bound(x_star, initial, left_bound = 0.001, right_bound = 0.5)
make_lower_bound<-function(x, initial, h_f = dlnorm){
  x_star<-sort(x);length_x<-length(x_star)
  T_k = sort(initial);length_initial<-length(initial)
  lower_matrix<-data.frame(lower_limit = T_k[1:(length_initial-1)], upper_limit = T_k[2:length_initial])
  y_lower <- vector("numeric",length_x)
  #The following for-loop is not efficient enough and should be changed
  #Note: x_star is increasing, so we can at least reduce some computation
  j_store = 1
  for(i in 1:length_x){
    for(j in j_store:(length_initial-1)){
      if(x_star[i]>=lower_matrix[j,1] & x_star[i]<=lower_matrix[j,2]){
        #T_k[j]s are the x[j]s in the paper, which are the initial points.
        y_lower[i] <- ((T_k[j+1]-x_star[i])*h_f(T_k[j]) + (x_star[i]-T_k[j])*h_f(T_k[j+1]))/(T_k[j+1]-T_k[j])
        j_store <- j
      }
      if(x_star[i]<lower_matrix[1,1]||x_star[i]>lower_matrix[(length_initial-1),2]){
        y_lower[i] <- -Inf
      }
    }
  }
  return(data.frame(x_star = x_star,lower_value = y_lower))
}
b<-make_lower_bound(x_star,initial)

#To compare the difference between lower and upper bound
a$upper_value - b$lower_value >= 0
# I want to draw the upper_bound function, do it later
seg<-make_z(seq, h_f = dnorm)$z[1:40]
segments(seg[1:39],dnorm(seg[1:39]),seg[2:40],dnorm(seg[2:40]))
plot(seg[1:39],dnorm(seg[1:39]),type = "b")
points(seg[2:40],dnorm(seg[2:40]))

x<-rpois(10,3)
a<-c(1:10)
b<-c(2:11)
f<-function(x) {x[x>=a[1:10]&x<=b[1:10]] +1} #beautiful function

