library(aciccomp2016)
library(tidyverse)
library(doParallel)
library(foreach)
library(tictoc)
library(comprehenr)
library(MASS)
library(quadprog)
library(glmnet)
library(evalITR)
library(ggplot2)
library(reshape2)
library(ggthemes)
library(latex2exp)
library(ggtext)


dgp_2016_sample <- function(truth, n, random.seed) {
  random.seed <- list(seed = random.seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
  suppressWarnings(do.call("set.seed", random.seed))
  
  selected_sample = sample(1:nrow(truth$x),n, replace = FALSE)
  mu.0 = truth$mu.0[selected_sample]
  mu.1 = truth$mu.1[selected_sample]
  x = truth$x[selected_sample,]
  sigma_y = 1
  
  n1=round(n/2)
  n0=n-n1
  zind=sample(1:n,size=n1)
  z=numeric(n)
  z[zind]=1
  mu = z * mu.1 + (1-z) * mu.0
  y <- mu + sigma_y*rnorm(n)
  result <- list(x = x, z = z, y =y, mu.1 = mu.1, mu.0 = mu.0)
  
}

ntrials = 10000
nseq = 100
no_cores <- 5
cl <- makeCluster(no_cores)
registerDoParallel(cl)
m = 1
input_x = input_2016
final_output8 = as.list(numeric(121))
dim(final_output8) <- c(11,11)
truth = dgp_2016(input_x, 28, 1,extraInfo=TRUE)
train_master <- dgp_2016_sample(truth, nrow(input_x), 2021)
std_dev = sd(truth$y)
y_train = train_master$y
x_train = train_master$x
x_train_expand = model.matrix(~.-1,data=x_train)
z_train = train_master$z
y_train_transformed = y_train * (z_train - mean(z_train)) / (mean(z_train) * (1- mean(z_train)))
ls1<-glmnet(x_train_expand,y_train_transformed,alpha = 1, lambda= 0.05)
tau_full1<-predict(ls1,x_train_expand)
Thatfull = as.numeric(tau_full1 >0)
shift_factor = 1/4*(mean(train_master$mu.1[Thatfull==1]) + mean(train_master$mu.0[Thatfull==0]))
for (s in -5:5) {
train <- train_master
train$y <- train$y - shift_factor + s
train$mu.0 <- train$mu.0 - shift_factor + s
train$mu.1 <- train$mu.1 - shift_factor + s
truePAPE <- function(train, That) {
  n = length(That)
  phat = mean(That)
  pape = (sum(train$mu.1[That==1]) + sum(train$mu.0[That==0])) / n - (phat * mean(train$mu.1) + (1-phat) * mean(train$mu.0))
  return(pape)
}
truePAV <- function(train, That) {
  n = length(That)
  pav = (sum(train$mu.1[That==1]) + sum(train$mu.0[That==0])) / n
  return(pav)
}
PAPEavg<- truePAPE(train ,Thatfull)
PAVavg <- truePAV(train, Thatfull)
tic()
expresult=foreach (l=1:ntrials) %dopar% {
  library(glmnet)
  library(dplyr)
  library(evalITR)
  test = dgp_2016_sample(train, nseq, l)
  y_test = test$y
  x_test = test$x
  z_test = test$z
  x_test_expand = model.matrix(~.-1,data=x_test)
  tau_test1<-predict(ls1,x_test_expand)
  That = as.numeric(tau_test1 >0)
  output = PAPE(z_test,That,y_test, budget= NA, centered = FALSE)
  output2 = PAV(z_test,That,y_test,centered = FALSE)
  pape = output$pape
  sd1 = output$sd
  pav = output2$pav
  sd2 = output2$sd
  return(list(pape+1.96*sd1, pape, pape-1.96*sd1, sd1,pav+1.96*sd2, pav, pav-1.96*sd2, sd2))
}
toc()
PAPEu=matrix(unlist(map(expresult,1)),nrow = ntrials, byrow = TRUE)
PAPEm=matrix(unlist(map(expresult,2)),nrow = ntrials, byrow = TRUE)
PAPEl=matrix(unlist(map(expresult,3)),nrow = ntrials, byrow = TRUE)
sd1=matrix(unlist(map(expresult,4)),nrow = ntrials, byrow = TRUE)
PAVu=matrix(unlist(map(expresult,5)),nrow = ntrials, byrow = TRUE)
PAVm=matrix(unlist(map(expresult,6)),nrow = ntrials, byrow = TRUE)
PAVl=matrix(unlist(map(expresult,7)),nrow = ntrials, byrow = TRUE)
sd2=matrix(unlist(map(expresult,8)),nrow = ntrials, byrow = TRUE)
coverage1=colSums(sweep(PAPEu,2,PAPEavg,">") & sweep(PAPEl,2,PAPEavg,"<"))/ntrials
coverage2=colSums(sweep(PAVu,2,PAVavg,">") & sweep(PAVl,2,PAVavg,"<"))/ntrials

rmse = sqrt(colMeans((sweep(PAPEm,2,PAPEavg))*(sweep(PAPEm,2,PAPEavg))))
mad = colMeans(abs(sweep(PAPEm,2,PAPEavg)))
asd = apply(PAPEm,2,sd)
rmsesd= sqrt(colMeans((sweep(sdm,2,asd))*(sweep(sdm,2,asd))))
madsd = colMeans(abs(sweep(sdm,2,asd)))
final_output8[[s+6,1]]=nseq
final_output8[[s+6,2]]=coverage1
final_output8[[s+6,3]]=coverage2
final_output8[[s+6,4]]=colMeans(PAPEm)
final_output8[[s+6,5]]=PAPEavg
final_output8[[s+6,6]]=colMeans(PAVm)
final_output8[[s+6,7]]=PAVavg
final_output8[[s+6,8]]=apply(PAPEm,2,sd)
final_output8[[s+6,9]]=apply(PAVm,2,sd)
final_output8[[s+6,10]]=colMeans(sd1)
final_output8[[s+6,11]]=colMeans(sd2)
print(s)
}
output_df1 = data.frame(index = seq(-5,5), pav_var = unlist(final_output8[,9]))
stopCluster(cl)
cairo_pdf("var_shift.pdf",width=4,height=4)
ggplot() + 
  geom_line(data = output_df1,aes(x=index, 
                                  y=pav_var),
            alpha=1, size = 1.5) + 
  xlab(TeX("$\\Delta$")) +
  ylab(TeX("$V(\\hat{\\tau}_f(Z))$")) + 
  theme_few() +
  coord_cartesian(ylim = c(0.74, 0.85), xlim = c(-5,5)) + 
  theme(axis.text=element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14))
dev.off()




ntrials = 10000
nseqs=c(100,500,1000,2000)
no_cores <- 5
cl <- makeCluster(no_cores)
registerDoParallel(cl)
input_x = input_2016
final_output9 = as.list(numeric(length(nseqs) * 5))
dim(final_output9) <- c(length(nseqs),5)
truth = dgp_2016(input_x, 28, 1,extraInfo=TRUE)
train_master <- dgp_2016_sample(truth, nrow(input_x), 2021)
std_dev = sd(truth$y)
y_train = train_master$y
x_train = train_master$x
x_train_expand = model.matrix(~.-1,data=x_train)
z_train = train_master$z
y_train_transformed = y_train * (z_train - mean(z_train)) / (mean(z_train) * (1- mean(z_train)))
ls1<-glmnet(x_train_expand,y_train_transformed,alpha = 1, lambda= 0.05)
tau_full1<-predict(ls1,x_train_expand)
Thatfull = as.numeric(tau_full1 >0)
dgp_2016_ante <- function(truth, n, random.seed) {
  random.seed <- list(seed = random.seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
  suppressWarnings(do.call("set.seed", random.seed))
  
  selected_sample = sample(1:nrow(truth$x),n, replace = FALSE)
  mu.0 = truth$mu.0[selected_sample]
  mu.1 = truth$mu.1[selected_sample]
  x = truth$x[selected_sample,]
  sigma_y = 1
  nf=round(n/2)
  nr=n-nf
  zfind=sample(1:n,size=nf)
  z=numeric(n)
  f =numeric(n)
  f[zfind]=1
  x_expand = model.matrix(~.-1,data=x[zfind,])
  tau_f<-predict(ls1,x_expand)
  Thatf = as.numeric(tau_f >0)
  z[zfind[Thatf==1]]=1
  zr = setdiff(1:n, zfind)
  nr1 = round(nr/2)
  nr0 = nr - nr1
  zrind=sample(zr,size=nr1)
  z[zrind]=1
  mu = z * mu.1 + (1-z) * mu.0
  y <- mu + sigma_y*rnorm(n)
  result <- list(x = x, z = z, f=f, y =y, mu.1 = mu.1, mu.0 = mu.0)
}
shift_factor = 1/4*(mean(train_master$mu.1[Thatfull==1]) + mean(train_master$mu.0[Thatfull==0]))
train <- train_master
train$y <- train$y - shift_factor
train$mu.0 <- train$mu.0 - shift_factor
train$mu.1 <- train$mu.1 - shift_factor
truePAPE <- function(train, That) {
  n = length(That)
  phat = mean(That)
  pape = (sum(train$mu.1[That==1]) + sum(train$mu.0[That==0])) / n - (phat * mean(train$mu.1) + (1-phat) * mean(train$mu.0))
  return(pape)
}
PAPEavg<- truePAPE(train ,Thatfull)
for (m in 1:length(nseqs)) {
nseq = nseqs[m]
tic()
expresult=foreach (l=1:ntrials) %dopar% {
  library(glmnet)
  library(dplyr)
  library(evalITR)
  test = dgp_2016_sample(train, nseq, l)
  testa = dgp_2016_ante(train, nseq, l)
  y_test = test$y
  x_test = test$x
  z_test = test$z
  x_test_expand = model.matrix(~.-1,data=x_test)
  tau_test1<-predict(ls1,x_test_expand)
  That = as.numeric(tau_test1 >0)
  output = PAPE(z_test,That,y_test, budget= NA, centered = FALSE)
  pape = output$pape
  n = length(testa$z)
  phat = mean(That)
  x_testa_expand = model.matrix(~.-1,data=testa$x)
  tau_testa1<-predict(ls1,x_testa_expand)
  Thata = as.numeric(tau_testa1 >0)
  phata = mean(Thata)
  pape_ante = n/(n-1)* (1/sum(testa$f) * sum(testa$y[testa$f==1])- (phata/sum(testa$z[testa$f==0]) * sum((1-testa$f)* testa$z * testa$y) + (1-phat)/(sum(1-testa$f)-sum(testa$z[testa$f==0])) * sum((1-testa$f)* (1-testa$z) * testa$y)))
  return(list(pape, pape_ante))
}
toc()
PAPEm=matrix(unlist(map(expresult,1)),nrow = ntrials, byrow = TRUE)
PAPEam=matrix(unlist(map(expresult,2)),nrow = ntrials, byrow = TRUE)
papesd = apply(PAPEm,2,sd)
papeasd = apply(PAPEam,2,sd)
final_output9[[m,1]]=nseq
final_output9[[m,2]]=colMeans(PAPEm)
final_output9[[m,3]]=colMeans(PAPEam)
final_output9[[m,4]]=papesd
final_output9[[m,5]]=papeasd
}
stopCluster(cl)
