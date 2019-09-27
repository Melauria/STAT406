dat <- read.csv("forestfires.csv", header = TRUE)
dat1 = dat[dat$area > 0,]
dat1 = dat1[dat1$area < 1000,]
dat1$area <- log(dat1$area)
dat1$Feb<-as.numeric(dat1$month == "feb")
dat1$Mar<-as.numeric(dat1$month == "mar")
dat1$Apr<-as.numeric(dat1$month == "apr")
dat1$May<-as.numeric(dat1$month == "may")
dat1$Jun<-as.numeric(dat1$month == "jun")
dat1$Jul<-as.numeric(dat1$month == "jul")
dat1$Aug<-as.numeric(dat1$month == "aug")
dat1$Sep<-as.numeric(dat1$month == "sep")
dat1$Oct<-as.numeric(dat1$month == "oct")
dat1$day <- NULL
dat1$month <- NULL
# dat1$area <- log(dat1$area)
# the summary of the full model
summary(lm(area ~., data = dat1))
############################################################################
# using stepwise + AIC to select a model
library(MASS)
null <- lm(area ~ 1, data=dat1)
full <- lm(area ~ ., data=dat1)
st <- stepAIC(null, scope=list(lower=null, upper=full), trace=FALSE)
# using 5-fold CV to check its MSPE
k <- 5
n <- nrow(dat1)
ii <- (1:n) %% k + 1
set.seed(123)
ii <- sample(ii)
pr.st <- rep(0, n)
for(j in 1:k) {
  x0 <- dat1[ii != j, ]
  null0 <- lm(area ~ 1, data=x0)
  full0 <- lm(area ~ ., data=x0) # needed for stepwise
  step.lm0 <- stepAIC(null0, scope=list(lower=null0, upper=full0), trace=FALSE)
  pr.st[ ii == j ] <- predict(step.lm0, newdata=dat1[ii==j,])
}
mspe.st <- mean( (dat1$area - pr.st)^2 )
mspe.st
# mspe for AIC stepwise is 3986.537

############################################################################
# Using Ridge Regression to select a model
library(glmnet)
y <- as.vector(dat1$area)
xm <- as.matrix(dat1[, -11])
lambdas <- exp( seq(-3, 10, length=100))
fittedobject <- glmnet(x=xm, y=y, lambda=rev(lambdas), family='gaussian', alpha=0)
plot(fittedobject, xvar='lambda', label=TRUE)
tmp <- cv.glmnet(x=xm, y=y, lambda=lambdas, nfolds=5, alpha=0, family='gaussian')
ridge.lambda = tmp$lambda.min
# optimal lambda = 170.95736
op.la <- 0
for(j in 1:20) {
  tmp <- cv.glmnet(x=xm, y=y, lambda=lambdas, nfolds=5, alpha=0, family='gaussian')
  op.la <- op.la + tmp$lambda.min # tmp$lambda.1se
}
op.la <- op.la / 20
xm <- scale(as.matrix(dat1[, -11]), scale=FALSE)
xm.svd <- svd(xm)
(est.edf <- sum( xm.svd$d^2 / ( xm.svd$d^2 + op.la ) ))
# edf = 8.55
############################################################################
# Using LASSO Regression to select a model
library(glmnet)
x <- as.matrix(dat1[,-11])
y <- as.vector(dat1$area)

set.seed(800)
fittedobject <- glmnet(x=x, y=y, alpha=1)
plot(fittedobject, xvar ='lambda', label = TRUE)
model = cv.glmnet(x=x,y=y, nfolds = 5, alpha=1)
op.la = model$lambda.min
coef(model, s=op.la)

############################################################################
# Compare MSPEs of Full model, stepAIC, Ridge, and LASSO model
library(MASS)
n <- nrow(xm)
k <- 5
ii <- (1:n) %% k + 1
set.seed(123)
N <- 100
mspe.la <- mspe.st <- mspe.ri <- mspe.f <- rep(0, N)
for(i in 1:N) {
  ii <- sample(ii)
  pr.la <- pr.f <- pr.ri <- pr.st <- rep(0, n)
  for(j in 1:k) {
    tmp.ri <- cv.glmnet(x=xm[ii != j, ], y=y[ii != j], lambda=lambdas,
                        nfolds=5, alpha=0, family='gaussian')
    tmp.la <- cv.glmnet(x = xm[ii != j, ], y = y[ii != j], lambda = lambdas, 
                        nfolds = 5, alpha = 1, family = "gaussian")
    null <- lm(area ~ 1, data=dat1[ii != j, ])
    full <- lm(area ~ ., data=dat1[ii != j, ])
    tmp.st <- stepAIC(null, scope=list(lower=null, upper=full), trace=0)
    
    pr.ri[ ii == j ] <- predict(tmp.ri, s='lambda.min', newx=xm[ii==j,])
    pr.la[ ii == j ] <- predict(tmp.la, s = "lambda.min", newx = xm[ii ==j, ])
    pr.st[ ii == j ] <- predict(tmp.st, newdata=dat1[ii==j,])
    pr.f[ ii == j ] <- predict(full, newdata=dat1[ii==j,])
  }
  mspe.ri[i] <- mean( (dat1$area - pr.ri)^2 )
  mspe.la[i] <- mean( (dat1$area - pr.la)^2)
  mspe.st[i] <- mean( (dat1$area - pr.st)^2 )
  mspe.f[i] <- mean( (dat1$area - pr.f)^2 )
}
boxplot(mspe.la, mspe.ri, mspe.st, mspe.f, names = c("LASSO", "Ridge", "Stepwise", 
      "Full"), col = c("steelblue", "gray80", "tomato", "springgreen"), cex.axis = 1, 
        cex.lab = 1, cex.main = 2)
mtext(expression(hat(MSPE)), side = 2, line = 2.5)



############################################################################
# Use Regression tree
library(rpart)
myc <- rpart.control(minsplit=1, cp = 1e-8)
set.seed(800)
tree.to <- rpart(area~., data = dat1, method = 'class', control = myc)
plot(tree.to,uniform=T,margin=0.1)
c <- tree.to$cptable[which.min(tree.to$cptable[,"xerror"]),"CP"]
tree.t3 <- prune(tree.to, cp=c)
plot(tree.t3, uniform=FALSE)
text(tree.t3, pretty=TRUE)


############################################################################
# Use Random Forest
library("randomForest")
set.seed(123)
dat.rf = randomForest(area~., data = dat1, ntree = 1000, ype='response')
dat.rf
plot(dat.rf)

k <- 5
n <- nrow(dat1)
ii <- (1:n) %% k + 1
set.seed(123)
b<-predict(dat.rf, type='response')
mean((dat1$area - b)^2)
############################################################################
# Use Regression Tree
# library(rpart)
# myc <- rpart.control(minsplit=3, cp = 1e-8)
# set.seed(800)
# tree.to <- rpart(area~., data = dat1, method = 'class', control = myc)
# plot(tree.to,uniform=T,margin=0.1)
# c <- tree.to$cptable[which.min(tree.to$cptable[,"xerror"]),"CP"]
# tree.t3 <- prune(tree.to, cp=c)
# plot(tree.t3, uniform=FALSE)














