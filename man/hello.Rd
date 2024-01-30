\name{HT}
\alias{HT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hotelling’s T-Square
%%  ~~function to do ... ~~
}
\description{ A more preferable test statistic is Hotelling’s $$T^2$$ and we will focus on this test.

To motivate Hotelling's T^2, consider the square of the t-statistic for testing a hypothesis regarding a univariate mean. Recall that under the null hypothesis, t has a distribution with n-1 degrees of freedom. Now consider squaring this test statistic as shown below:


When you square a t-distributed random variable with n-1 degrees of freedom, the result is an F-distributed random variable with 1 and n-1 degrees of freedom. We reject H_0 at level $$\alpha$$ if t^2 is greater than the critical value from the F-table with 1 and n-1 degrees of freedom, evaluated at level
.
 \deqn {F = \frac{(n_1 + n_2 − p − 1) T^2}{(n_1 + n_2 − 2) p}∼ F_{p,n_1+n_2−p−1}}{%F = {(n_1 + n_2 − p − 1) × T^2}/{(n_1 + n_2 − 2) × p}∼ F_{p,n_1+n_2−p−1}}


%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hote_one(x)
hote_Two(data1 = data1,data2 = data2,delta = delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data}
  \item{n}{Number of row}
  \item{p}{Number of coloum}
  \item{x_bar}{Mean}
  \item{Sigma_inv}{inverce of variance matrix }
  \item{t2}{Value of t square Test}
  \item{f}{F Test}
  \item{p_val}{p value}

____________________________________________
 \item{data1}{First data}
 \item{data2}{Second data}
 \item{p}{Number of column}
 \item{n1}{First data Number of row}
 \item{n2}{Second data Number of row}
 \item{x1_bar}{First data mean}
 \item{x2_bar}{Second data mean}
 \item{s1}{First data variance}
 \item{s2}{Second data variance}
 \item{s_pool}{Poll variance}
 \item{t2}{Value of t square Test}
 \item{f}{F Test}
 \item{p_val}{p value}

%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
https://online.stat.psu.edu/stat505/lesson/7/7.1/7.1.3

https://online.stat.psu.edu/stat505/lesson/7/7.1/7.1.12
%% ~put references to the literature/web site here ~
}
\author{
Abu Sufiun Rhyme

Md Monowar Hossain(Apon)

%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


## The function is currently defined as
#### Q.1
A company wants to assess whether there is a significant difference in the purchasing behavior of its customers across multiple product categories. They collect data from a random sample of 80 customers and record their average spending in three different product categories: electronics (X1), clothing (X2), and home appliances (X3). The data is contained in the expenditure.csv file. The company wants to determine if the average spending is significantly different from the hypothesized population mean vector [300, 100, 250 250]. Perform the appropriate test to answer this.


hote_one<-function(x,mu=rep( 0,times= ncol(x))){

n<-nrow(x)
p<-ncol(x)
x_bar<-apply(X = x,MARGIN = 2,FUN = mean)
sigma_inv<-solve(var(x))
t2<-n*t(x_bar-mu) %*% sigma_inv %*% (x_bar-mu)
f<- (n-p)/(p*(n-1))*t2
p_val<-pf(q = f,df1 = p,df2 = n-p,lower.tail = F)

return(list(t2=t2,F=f,p_value=p_val))

  }


library(mvtnorm)
n <- 100 # observations
p <- 18 # parameter dimension
mu <- rep(0,p) # no signal: mu=0
x <- rmvnorm(n = n, mean = mu)
dim(x)

hote_one(x)

#install.packages("DEoptimR")
#install.packages("rrcov")
library(rrcov)
library(DEoptimR)
#DescTools::HotellingsT2Test(x = x)
rrcov::T2.test(x)

# Ques 1.
library(readxl)
expenditure <- read_excel("D:/3rd Year/AST-332 batch-27/expenditure.xlsx")
#View(expenditure)

hote_one(expenditure[,-1], mu=c(300,100,250))
rrcov::T2.test(expenditure[,-1], mu=c(300,100,250))

### Two sample Hotelling T square test
https://online.stat.psu.edu/stat505/lesson/7/7.1/7.1.12

#### Q.2
Test whether the mean differs between the following two groups:

||Group 01|| Group 02|
|:----:|:----:|:----:|:----:|
|x11| x12| x21| x22|
|6 |27 |25 |15 |
|6 |23 |28 |13 |
|18 |64 |36 |22 |
|8 |44 |35 |29 |
|11 |30 |15 |31 |
|34 |75 |44 |64 |
|28 |26 |42 |30 |
|71 |124 |54 |64 |
|43 |54 |34 |56 |
|33 |30 |29 |20 |
|20 |14 |39 |21 |
Help: Consider the following hypotheses are to be tested:
$$\begin{array}{c}
H_0:\left[\begin{array}{l}
\mu_{11} \\
\mu_{12}
\end{array}\right]-\left[\begin{array}{l}
\mu_{21} \\
\mu_{22}
\end{array}\right]=\left[\begin{array}{l}
\delta_1 \\
\delta_2
\end{array}\right]=\left[\begin{array}{l}
0 \\
0
\end{array}\right] \\
H_1:\left[\begin{array}{l}
\mu_{11} \\
\mu_{12}
\end{array}\right]-\left[\begin{array}{l}
\mu_{21} \\
\mu_{22}
\end{array}\right] \neq\left[\begin{array}{l}
\delta_1 \\
\delta_2
\end{array}\right] \neq\left[\begin{array}{l}
0 \\
0
\end{array}\right]
\end{array}$$


The corresponding test statistics is given by
$$F = \frac{(n_1 + n_2 − p − 1) × T^2}{(n_1 + n_2 − 2) × p}∼ F_{p,n_1+n_2−p−1}$$

where
$$(1) T^2 = (\bar x −\bar y − δ)′ ×  S_p^{-1} × (\bar x −\bar y − δ)$$
$$(2) S_p = \left( \frac {1}{n_1} + \frac {1}{n_2} \right)× S_{pooled}$$
$$(3) S_{pooled} = \frac {(n_1−1)S_x+(n_2−1)S_y}
{n_1+n_2−2}$$



hote_Two <- function(data1, data2, delta) {
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

#1st Population
x11 <- c(6, 6, 18, 8, 11, 34, 28, 71, 43, 33, 20)
x12 <- c(27, 23, 64, 44, 30, 75, 26, 124, 54, 30, 14)
#2nd Population
x21 <- c(25, 28, 36, 35, 15, 44, 42, 54, 34, 29, 39)
x22 <- c(15, 13, 22, 29, 31, 64, 30, 64, 56, 20, 21)
data1 <- data.frame(x11, x12)
data2 <- data.frame(x21, x22)
delta <- rep(0, 2)

hote_Two(data1 = data1,data2 = data2,delta = delta)
rrcov::T2.test(x = data1,y = data2)


# When the data is in proper (long) format
dat <- data.frame(group=factor(rep(1:2, each=11)),
                  x1= c(x11, x21),
                  x2= c(x12, x22))
rrcov::T2.test(cbind(x1,x2)~group, data=dat)


#### Ques 3
Test whether there exists mean difference between two species of flower:‘setosa’ & ‘versicolor’for built in data set ‘iris’.



data("iris")
data1 <- iris[iris$Species == "setosa", ]
data1
data1 <- data1[,-5]
data2 <- iris[iris$Species == "versicolor", ]
data2 <- data2[,-5]
delta <- rep(0, 4)

hote_Two(data1, data2, delta)
rrcov::T2.test(x = data1, y=data2)


{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

