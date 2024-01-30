hotelling_onesample<-function(x,mu=rep( 0,times= ncol(x))){

  n<-nrow(x)
  p<-ncol(x)
  x_bar<-apply(X = x,MARGIN = 2,FUN = mean)
  sigma_inv<-solve(var(x))
  t2<-n*t(x_bar-mu) %*% sigma_inv %*% (x_bar-mu)
  f<- (n-p)/(p*(n-1))*t2
  p_val<-pf(q = f,df1 = p,df2 = n-p,lower.tail = F)

  return(list(t2=t2,F=f,p_value=p_val))

}


hotelling_twosample <- function(data1, data2, delta) {
  data1 <- as.matrix(data1)
  data2 <- as.matrix(data2)
  p <- dim(data1)[2]
  n1 <- dim(data1)[1]
  n2 <- dim(data2)[1]
  x1_bar = as.matrix(apply(X = data1, MARGIN = 2, FUN = mean))
  x2_bar = as.matrix(apply(X = data2, MARGIN = 2, FUN = mean))
  s1 = as.matrix(cov(data1))
  s2 = as.matrix(cov(data2))

  s_pool = as.matrix((s1 * (n1 - 1) + s2 * (n2 - 1)) / (n1 + n2 - 2))
  m <- x1_bar - x2_bar - as.matrix(delta)

  t2 <- t(m) %*% solve(s_pool * ((1 / n1) + (1 / n2))) %*% m
  f = (n1 + n2 - p - 1) / (p * (n1 + n2 - 2)) * t2
  p_val = pf(
    q = f,
    df1 = p,
    df2 = n1 + n2 - p - 1,
    lower.tail = F
  )

  return(c(
    T2 = t2,
    F_stat = f,
    p_value = p_val
  ))
}
