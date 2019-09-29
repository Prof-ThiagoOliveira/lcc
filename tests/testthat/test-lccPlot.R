# Testing lcc plot

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
  par<-data.frame("Par"=c(Beta01,Beta02,Beta03,Beta04,Beta11,Beta12,
                          Beta13,Beta14,
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

  LCC12T <-
    Vectorize(function(t) 2*COV(t)/(VAR1(t)+VAR1(t)+
                                      (E1(t)-E2(t))^2), "t")
  LCC13T <-
    Vectorize(function(t) 2*COV(t)/(VAR2(t)+VAR2(t)+
                                      (E1(t)-E3(t))^2), "t")
  LCC14T <-
    Vectorize(function(t) 2*COV(t)/(VAR1(t)+VAR2(t)+
                                      (E1(t)-E4(t))^2), "t")

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
test_that("Object does not inherit from class lcc", {
  aaa<-lm(rnorm(100,10,1)~rnorm(100,50,3))
  expect_error(lccPlot(aaa), "Object must inherit from class \"lcc\"")
})
data(hue)
test_that("LCC, LPC and LA plot test",{
  # Two methods
  expect_that(fme1<-lcc(dataset = hue, subject = "Fruit",
                        resp = "H_mean", method = "Method",
                        time = "Time", qf = 1, qr = 1),is_a("lcc"))
  tmp1<-tempfile()
  expect_known_output(lccPlot(fme1), tmp1)
  ## Components TRUE
  expect_that(fme2<-lcc(dataset = hue, subject = "Fruit",
                        resp = "H_mean", method = "Method",
                        time = "Time", qf = 1, qr = 1,
                        components = TRUE),is_a("lcc"))
  tmp2<-tempfile()
  expect_known_output(lccPlot(fme2), tmp2)
  # More than two methods
  expect_that(fme3<-lcc(dataset = dataset$data, subject = "Fruit",
                        resp = "Response", method = "Method",
                        time = "Time", qf = 1, qr = 1),is_a("lcc"))
  tmp3<-tempfile()
  expect_known_output(lccPlot(fme3), tmp3)
  ## Components TRUE
  expect_that(fme4<-lcc(dataset = dataset$data, subject = "Fruit",
                        resp = "Response", method = "Method",
                        time = "Time", qf = 1, qr = 1,
                        components = TRUE),is_a("lcc"))
  tmp4<-tempfile()
  expect_known_output(lccPlot(fme4), tmp4)
 })

test_that("Confidence intervals plot",{
  # Two methods
  expect_that(fme5<-lcc(dataset = hue, subject = "Fruit",
                        resp = "H_mean", method = "Method",
                        time = "Time", qf = 1, qr = 1, ci = TRUE,
                        nboot = 100),is_a("lcc"))
  tmp5<-tempfile()
  expect_known_output(lccPlot(fme5), tmp5)
  ## Components TRUE
  expect_that(fme6<-lcc(dataset = hue, subject = "Fruit",
                        resp = "H_mean", method = "Method",
                        time = "Time", qf = 1, qr = 1,
                        components = TRUE, ci = TRUE,
                        nboot = 100),is_a("lcc"))
  tmp6<-tempfile()
  expect_known_output(lccPlot(fme6), tmp6)
  # More than two methods
  expect_that(fme7<-lcc(dataset = dataset$data, subject = "Fruit",
                        resp = "Response", method = "Method",
                        time = "Time", qf = 1, qr = 1, ci=TRUE,
                        nboot = 100),is_a("lcc"))
  tmp7<-tempfile()
  expect_known_output(lccPlot(fme7), tmp7)
  # Components TRUE
  expect_that(fme8<-lcc(dataset = dataset$data, subject = "Fruit",
                        resp = "Response", method = "Method",
                        time = "Time", qf = 1, qr = 1,
                        components = TRUE, ci=TRUE, nboot = 100),
              is_a("lcc"))
  tmp8<-tempfile()
  expect_known_output(lccPlot(fme8), tmp8)
})

test_that("labels, shape and colour",{
  expect_that(fm<-lcc(dataset = dataset$data, subject = "Fruit",
                      resp = "Response", method = "Method",
                      time = "Time", qf = 1, qr = 1,
                      components = TRUE),is_a("lcc"))
  tmp<-tempfile()
  expect_known_output(lccPlot(fm, type = "lcc", control=list(
    shape = 2, colour = "red", size = 2,
    xlab = "Time (hours)", ylab = "Longitudinal CC")), tmp)
  tmp2<-tempfile()
  expect_known_output(lccPlot(fm, type = "lpc", control=list(
    shape = 2, colour = "red", size = 2,
    xlab = "Time (hours)", ylab = "Longitudinal PC")), tmp2)
  tmp3<-tempfile()
  expect_known_output(lccPlot(fm, type = "la", control=list(
    shape = 2, colour = "red", size = 2,
    xlab = "Time (hours)", ylab = "Longitudinal Accuracy")), tmp3)
})

test_that("Scales",{
  a<-c(0,1)
  b<-c(-0.5,1)
  c<-c(-0.2,1)
  expect_that(fm<-lcc(dataset = dataset$data, subject = "Fruit",
                      resp = "Response", method = "Method",
                      time = "Time", qf = 1, qr = 1,
                      components = TRUE),is_a("lcc"))
  tmp<-tempfile()
  expect_known_output(lccPlot(fm,
                              control=list(scale_y_continuous = a)),
                      tmp)
  tmp2<-tempfile()
  expect_known_output(lccPlot(fm, type = "lpc",
                              control=list(scale_y_continuous = b)),
                      tmp2)
  tmp3<-tempfile()
  expect_known_output(lccPlot(fm, type = "la",
                              control=list(scale_y_continuous = c)),
                      tmp3)
})

