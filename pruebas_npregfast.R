
# allotest listo
# frfast listo


library(NPRegfast)
data(barnacle)
fit <- frfast(DW~RC:F, data = barnacle, model = "np", h=c(0.3))
fit
plot(fit, der = c(0))
summary(fit)



quartz()
plot(fit$x,fit$p[,1,1])
matplot(fit$x,fit$repboot[,1,1,],type="l")
lines(fit$x,fit$p[,1,1],lwd=4.5)

maxp(fit)


names(fit)






fit2 <- frfast(DW~RC:F, data = barnacle, model = "np")
fit2
plot(fit2, fac = 2, der = 0)



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
