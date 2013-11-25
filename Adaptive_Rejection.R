y <- vector("numeric",length_x)
#The following for-loop is not efficient enough and should be changed
#Note: x_star is increasing, so we can at least reduce some computation
j_store = 1
for(i in 1:length_x){
  for(j in j_store:(length_initial-1)){
    if(x_star[i]>=lower_matrix[j,1] & x_star[i]<=lower_matrix[j,2]){
      #T_k[j]s are the x[j]s in the paper, which are the initial points.
      y[i] <- ((T_k[j+1]-x_star[i])*h_f(T_k[j]) + (x_star[i]-T_k[j])*h_f(T_k[j+1]))/(T_k[j+1]-T_k[j])
      j_store <- j
    }
    if(x_star[i]<lower_matrix[1,1]||x_star[i]>lower_matrix[(length_initial-1),2]){
      y[i] <- -Inf
    }
  }
}