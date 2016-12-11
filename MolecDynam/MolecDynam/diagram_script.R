energy <- read.csv("data/energy.csv", header=FALSE)
kinetic <- read.csv("data/kinetic.csv", header=FALSE)
poten <- read.csv("data/poten.csv", header=FALSE)
temp <- read.csv("data/temp.csv", header=FALSE)
vel <- read.csv("data/velocity.csv", header=FALSE)

png('images/energy.png')
plot(energy$V2, energy$V1, type="l", main=abs(max(abs(energy$V1-mean(energy$V1)))/mean(energy$V1)))
abline(a=mean(energy$V1),b=0)
dev.off()

png('images/temp.png')
#plot(temp$V2, temp$V1, type="l", main=abs(max(abs(temp$V1-mean(temp$V1)))/mean(temp$V1)))
step <- 0
plot(temp$V2[step:nrow(temp)], temp$V1[step:nrow(temp)], type="l", main=mean(temp$V1[step:nrow(temp)]), xlab=max(abs(temp$V1[step:nrow(temp)]-mean(temp$V1[step:nrow(temp)])))/mean(temp$V1[step:nrow(temp)]))
abline(a=mean(temp$V1[step:nrow(temp)]),b=0, col="red")
dev.off()

png('images/temp2.png')
step <- 0
plot(temp$V2[step:nrow(temp)], temp$V1[step:nrow(temp)], type="l", main=mean(temp$V1[step:nrow(temp)]), xlab=max(abs(temp$V1[step:nrow(temp)]-mean(temp$V1[step:nrow(temp)])))/mean(temp$V1[step:nrow(temp)]))
abline(a=mean(temp$V1[step:nrow(temp)]),b=0, col="red")
dev.off()

png('images/kinetic.png')
plot(kinetic$V2, kinetic$V1, type="l", main=abs(max(abs(kinetic$V1-mean(kinetic$V1)))/mean(kinetic$V1)))
abline(a=mean(kinetic$V1),b=0)
dev.off()

png('images/poten.png')
plot(poten$V2, poten$V1, type="l", main=abs(max(abs(poten$V1-mean(poten$V1)))/mean(poten$V1)))
abline(a=mean(poten$V1),b=0)
dev.off()

png('images/compare.png')
plot(energy$V2, energy$V1, type="l", xlab="", ylab="", ylim=range(c(energy$V1,kinetic$V1,poten$V1)))
par(new=TRUE)
plot(kinetic$V2, kinetic$V1, type="l", xlab="", ylab="", axes=FALSE, ylim=range(c(energy$V1,kinetic$V1,poten$V1)))
par(new=TRUE)
plot(poten$V2, poten$V1, type="l", xlab="", ylab="", axes=FALSE, ylim=range(c(energy$V1,kinetic$V1,poten$V1)))
dev.off()

png('images/vel.png')
Temp <- mean(temp$V1[(nrow(temp)/2):nrow(temp)])
h <- hist(vel$V1)
par(new=TRUE)
curve(((1/(2*Temp*pi))^(3/2))*4*pi*(x^2)*exp(-(x^2)/(2*Temp)), 0, max(h$breaks), ylab='', yaxt='n', col='red', axes=FALSE)
dev.off()
