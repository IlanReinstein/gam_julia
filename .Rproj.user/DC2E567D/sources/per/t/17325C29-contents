size<-c(1.42,1.58,1.78,1.99,1.99,1.99,2.13,2.13,2.13, 2.32,2.32,2.32,2.32,2.32,2.43,2.43,2.78,2.98,2.98)
wear<-c(4.0,4.2,2.5,2.6,2.8,2.4,3.2,2.4,2.6,4.8,2.9, 3.8,3.0,2.7,3.1,3.3,3.0,2.8,1.7) 
x<-size-min(size);x<-x/max(x)
plot(x,wear,xlab="Scaled engine size",ylab="Wear index")

rk<-function(x,z) # R(x,z) for cubic spline on [0,1] 
  { ((z-0.5)^2-1/12)*((x-0.5)^2-1/12)/4-
  ((abs(x-z)-0.5)^4-(abs(x-z)-0.5)^2/2+7/240)/24 }

spl.X<-function(x,xk)
  # set up model matrix for cubic penalized regression spline
  { q<-length(xk)+2
  n<-length(x)
  X<-matrix(1,n,q)
  X[,2]<-x
  X[,3:q]<-outer(x,xk,FUN=rk) # and remaining to R(x,xk) 
  X
  }
xk<-1:4/5 # choose some knots 
X<-spl.X(x,xk) # generate model matrix
mod.1 <- lm(wear ~ X - 1) # fit model
xp<-0:100/100 # x values for prediction 
Xp<-spl.X(xp,xk) # prediction matrix 
lines(xp,Xp%*%coef(mod.1)) #



1:7/8
spl.S<-function(xk)
  # set up the penalized regression spline penalty matrix,
  # given knot sequence xk
{ q<-length(xk)+2;S<-matrix(0,q,q) # initialize matrix to 0
S[3:q,3:q]<-outer(xk,xk,FUN=rk) # fill in non-zero part
S }
mat.sqrt<-function(S) # A simple matrix square root 
  { d<-eigen(S,symmetric=TRUE)
  rS<-d$vectors%*%diag(d$values^0.5)%*%t(d$vectors) }

prs.fit<-function(y,x,xk,lambda)
  # function to fit penalized regression spline to x,y data, # with knots xk, given smoothing parameter, lambda
  { q<-length(xk)+2 # dimension of basis 
  n<-length(x) # number of data
  # create augmented model matrix ....
  Xa <- rbind(spl.X(x,xk),mat.sqrt(spl.S(xk))*sqrt(lambda)) 
  y[(n+1):(n+q)]<-0 # augment the data vector
  lm(y ~ Xa-1) # fit and return the penalized regression spline 
  }

xk<-1:7/8 # choose some knots 
mod.2<-prs.fit(wear,x,xk,0.0001) # fit pen. reg. spline 
Xp<-spl.X(xp,xk) # matrix to map params to fitted values at xp 
plot(x,wear);lines(xp,Xp%*%coef(mod.2)) # plot data & spl. fit



lambda<-1e-8
n<-length(wear)
V<-0


mod<-prs.fit(wear,x,xk,lambda) # fit model, given lambda

# X = model.matrix(mod) # X model matrix
# hat_matrix = X%*%(solve(t(X)%*%X)%*%t(X)) # Hat matrix
# diag(hat_matrix)[1] # First diagonal point in Hat matrix
# fitwithout1 = lm(mpg ~ wt, mtcars[-1,]) # OLS excluding first data point.
# new = data.frame(wt=mtcars[1,'wt']) # Predicting y hat in this OLS w/o first point.
# y_hat_without = predict(fitwithout1, newdata=new) # ... here it is.
# residuals(fit)[1] # The residual when OLS includes data point.
# lev = 1 - (residuals(fit)[1]/(mtcars[1,'mpg'] -  y_hat_without)) # Leverage
# all.equal(diag(hat_matrix)[1],lev) #TRUE


trA<-sum(influence(mod)$hat[1:n]) # find tr(A) 
rss<-sum((wear-fitted(mod)[1:n])^2)# residual sum of squares 
V[i]<-n*rss/(n-trA)^2 # obtain GCV score 
lambda<-lambda*1.5 # increase lambda

plot(1:60,V,type="l",main="GCV score",xlab="i") # plot score


######################
q <- 10
data(trees)
quantile(unique(trees$Girth), 1:(q - 2)/(q - 1))




