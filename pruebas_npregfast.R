
# allotest listo
# frfast listo
# interpret.frfastformula listo
# maxp listo
# maxp.diff listo
# print.frfast listo (no lleva ayuda)
# summary.frfast listo
# predict.frfast listo
# plot.frfast listo
# plot.diff listo
# localtest listo
# globaltest listo


# R CMD Rd2pdf npregfast


library(formatR)
tidy_source("globaltest.R",arrow = TRUE, width.cutoff = 80)



library(NPRegfast)
data(barnacle)
barnacle2 <-barnacle
barnacle2$F <- as.factor(barnacle2$F)
levels(barnacle2$F) <- c("oia","muros")

fit <- frfast(DW~RC : F, data = barnacle2, model = "np", kbin = 100)
fit
plot(fit, der = c(0,1), fac = "oia")
summary(fit)



#predict
nd2 <- data.frame(RC = c(10, 15, 20))
pred <- predict.frfast(fit, newdata = nd2)







quartz()
plot(fit$x,fit$p[,1,1])
matplot(fit$x,fit$repboot[,1,1,],type="l")
lines(fit$x,fit$p[,1,1],lwd=4.5)

maxp(fit)


names(fit)









allotest(DW ~ RC : F, data = barnacle, nboot = 300, kbin = 150, seed = 1353)


x <- runif(1000, 0, 2)
a <- 1
b <- 3
y <- a + b * x 
plot(x, y)



x=exp(x)
y=exp(y)


plot(x, y )
data <- data.frame(x, y)

allotest(y ~ x, data = data, nboot = 200, kbin = 150)






setwd("/Users/Sestelo/Dropbox/NO_ACABADO/Marco/nuevo")


datos=read.table("marco.txt",header=T)
datos

strata=as.numeric(datos$strata)
quartz()
plot(datos$Longitud,datos$Peso,col=2)

allotest(Peso ~ Longitud, data = datos, nboot = 500, kbin = 300)
allo <- frfast(Peso~Longitud, data = datos, model = "allo")
plot(allo, der = 0)

np <- frfast(Peso~Longitud, data = datos)
plot(np, der = 0)
