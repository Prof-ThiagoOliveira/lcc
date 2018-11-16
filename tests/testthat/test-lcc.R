stopifnot(require("testthat"), require("lcc"))

context("fitting lcc models")

data(hue)

test_that("Name of variables", {
  expect_error(lcc(dataset = ddd, subject = "Fruit", resp = "H_mean",
    method = "Method", time = "Time", qf = 2, qr = 2),"object 'ddd' not found")
  expect_error(lcc(dataset = hue, subject = "AAAAA", resp = "H_mean",
    method = "Method", time = "Time", qf = 2, qr = 2),"Please, verify the name of 'resp', 'subject', 'method', and 'time' variables")
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "AAAAA",
    method = "Method", time = "Time", qf = 2, qr = 2),"Please, verify the name of 'resp', 'subject', 'method', and 'time' variables")
  expect_error(lcc(dataset = hue, subject = "Fruits", resp = "H_mean",
    method = "AAAAA", time = "Time", qf = 2, qr = 2),"Please, verify the name of 'resp', 'subject', 'method', and 'time' variables")
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean",
    method = "Method", time = "AAAAA", qf = 2, qr = 2),"Please, verify the name of 'resp', 'subject', 'method', and 'time' variables")
  })

test_that("If qr>qf should result in error",{
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 3),"'qr' should be less or equal 'qf'")
 })

test_that("pdmat" ,{
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, pdmat = pdIdent()),"Do not include brackets after the pdmat function, e.g. pdSymm()")
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, pdmat = AAAA),"object 'AAAA' not found")

# compatibility names in pdmat
  expect_that(fm1<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, pdmat = "pdIdent()"), is_a("lcc"))
  
  expect_that(fm2<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, pdmat = pdIdent), is_a("lcc"))
  
  expect_that(fm3<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, pdmat = "pdIdent"), is_a("lcc"))
  
  expect_equivalent(fm1,fm2)
  expect_equivalent(fm1,fm3)
  expect_equivalent(fm2,fm3)
})



test_that("var.class and weights.form" ,{
  # without declare weights.forms
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, var.class = varIdent),"Please specify the 'weights.form' argument.")
  
  # wrong name for var.class
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, var.class = varIdent(), weights.form = "method"),"Do not include brackets after the var.class function, e.g. varExp()")
  
  # wrong name for weights.form
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, var.class = varIdent, weights.form = "AAAA"),"The weights.form argument are \"time\", \"method\", \"time.ident\", or \"both\".")
  
  # wrong name for weights.form
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, var.class = varIdent, weights.form = "both"),"Please specify the 'weight.form' correctly for varIdent class")
  
  # Wrong name for var.class and weights.form
  expect_error(lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, var.class = AAA, weights.form = "AAA"),"object 'AAA' not found")
  })

test_that("time_lcc" ,{
  # Testing regular sequence
  expect_that(fm1<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, time_lcc =list(from=min(hue$Time), to=max(hue$Time), n=30)),is_a("lcc"))
  expect_that(fm2<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 2, qr = 2, time_lcc =list(time=seq(0,14,1), from=0, to=14, n=30)),is_a("lcc"))
  expect_equivalent(fm1,fm2)
})

# Testing estimates

## Simulate dataset
Data<-function(N,time){
  Beta01<- 120 # Method A and Region A
  Beta02<- 130 # Method A and Region B
  Beta03<- 120 # Method B and Region A
  Beta04<- 105 # Method B and Region B
  Beta11<--3 # Method A and Region A
  Beta12<--3.5 # Method A and Region B
  Beta13<--3.01 # Method B and Region A
  Beta14<--2.5 # Method B and Region B
  sigmab0<-10
  sigmab1<-0.5
  corb01<-0.8
  covb0b1<-corb01*sqrt(sigmab0*sigmab1)
  error<-c(0.5,0.5)
  par<-data.frame("Par"=c(Beta01,Beta02,Beta03,Beta04,Beta11,Beta12,Beta13,Beta14,
    sigmab0,sigmab1,corb01,error[1],error[2]))
  par$Parameters<-c("Beta01","Beta02","Beta03","Beta04",
    "Beta11","Beta12","Beta13","Beta14",
    "sigmab0","sigmab1","corb01",
    "error[1]","error[2]")
  
  VAR1<-function(time) {
    sigmab0+sigmab1*time^2+2*covb0b1*time+error[1]
  }
  
  VAR2<-function(time) {
    sigmab0+sigmab1*time^2+2*covb0b1*time+error[2]
  }
  
  COV<-function(time) {
    sigmab0+sigmab1*time^2+2*covb0b1*time
  }
  
  E1<-function(time) sum(c(1,time)*c(Beta01,Beta11)) 
  E2<-function(time) sum(c(1,time)*c(Beta02,Beta12)) 
  E3<-function(time) sum(c(1,time)*c(Beta03,Beta13)) 
  E4<-function(time) sum(c(1,time)*c(Beta04,Beta14)) 
  
  mu1T<-Vectorize(function(t) E1(t), "t")
  mu2T<-Vectorize(function(t) E2(t), "t")
  mu3T<-Vectorize(function(t) E3(t), "t")
  mu4T<-Vectorize(function(t) E4(t), "t")
  
  S12T<-Vectorize(function(t) E1(t)-E2(t), "t")
  S13T<-Vectorize(function(t) E1(t)-E3(t), "t")
  S14T<-Vectorize(function(t) E1(t)-E4(t), "t")

  LCC12T <- Vectorize(function(t) 2*COV(t)/(VAR1(t)+VAR1(t)+(E1(t)-E2(t))^2), "t")
  LCC13T <- Vectorize(function(t) 2*COV(t)/(VAR2(t)+VAR2(t)+(E1(t)-E3(t))^2), "t")
  LCC14T <- Vectorize(function(t) 2*COV(t)/(VAR1(t)+VAR2(t)+(E1(t)-E4(t))^2), "t")

  N<-N
  Fruit<-gl(N,length(time),4*N*length(time))
  Method<-gl(4,k=length(time)*N)
  Time<-rep(time, N*4)
  
  E1y<-Vectorize(function(t) E1(t), "t")
  E2y<-Vectorize(function(t) E2(t), "t")
  E3y<-Vectorize(function(t) E3(t), "t")
  E4y<-Vectorize(function(t) E4(t), "t")
  
  Evalue<-c(rep(E1y(time),N),
    rep(E2y(time),N),
    rep(E3y(time),N),
    rep(E4y(time),N))
  
  require(MASS)
  bi<-mvrnorm(N,mu=c(0,0), Sigma = matrix(c(sigmab0,covb0b1,
    covb0b1,sigmab1),
    byrow=TRUE,ncol=2), empirical=FALSE)
  
  Zi<-model.matrix(~time)
  Zb.<-list(NA)
  for(i in 1:N){
    Zb.[[i]]<-Zi%*%bi[i,]
  }
  Zb<-rep(unlist(Zb.),4)
  
  residual<-c(rnorm(length(Zb)/2,mean=0,sd=sqrt(error[1])),
    rnorm(length(Zb)/2,mean=0,sd=sqrt(error[2])))
  
  Response=Evalue+Zb+residual
  
  
  dataset<-data.frame(Fruit,Method,Response, Time)
  
  Time2<-time
  return(list("data"=dataset, "par"=par, 
    "LCC12T"=LCC12T(Time2),
    "LCC13T"=LCC13T(Time2),
    "LCC14T"=LCC14T(Time2),
    "mu1T"=mu1T(Time2),
    "mu2T"=mu2T(Time2),
    "mu3T"=mu3T(Time2),
    "mu4T"=mu4T(Time2),
    "S12T"=S12T(Time2),
    "S13T"=S13T(Time2),
    "S14T"=S14T(Time2)))
}

set.seed(5925670)
dataset<-Data(N=30,time=seq(0,15,1))
test_that("Testing LCC estimates", {
  expect_that(fme1<-lcc(dataset = dataset$data, subject = "Fruit", resp = "Response", method = "Method", time = "Time", qf = 1, qr = 1),is_a("lcc"))
  expect_equal(fme1$Summary.lcc$fitted[[1]][,2], dataset$LCC12T, tolerance = .01)
  expect_equal(fme1$Summary.lcc$fitted[[2]][,2], dataset$LCC13T, tolerance = .01)
  expect_equal(fme1$Summary.lcc$fitted[[3]][,2], dataset$LCC14T, tolerance = .01)
})

# Test for confidence intervals
test_that("Test if confidence interval works",{
  expect_that(fme2<-lcc(dataset = dataset$data, subject = "Fruit", resp = "Response", method = "Method", time = "Time", qf = 1, qr = 1, ci=TRUE, nboot = 500, components = TRUE),is_a("lcc"))
  expect_that(fme3<-lcc(dataset = dataset$data, subject = "Fruit", resp = "Response", method = "Method", time = "Time", qf = 1, qr = 1, ci=TRUE, nboot = 500, components = TRUE, percentileMet = TRUE),is_a("lcc"))
  expect_that(fme4<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 1, qr = 1, ci=TRUE, nboot = 500, components = TRUE),is_a("lcc"))
  expect_that(fme5<-lcc(dataset = hue, subject = "Fruit", resp = "H_mean", method = "Method", time = "Time", qf = 1, qr = 1, ci=TRUE, nboot = 500, components = TRUE, percentileMet = TRUE),is_a("lcc"))
  expect_equal(fme2$Summary.lcc$fitted$LCC, fme3$Summary.lcc$fitted$LCC, tolerance = 0.05)
  expect_equal(fme4$Summary.lcc$fitted$LCC, fme5$Summary.lcc$fitted$LCC, tolerance = 0.05)
})

