temp <- read.csv("temp.csv", header=FALSE)
plot(temp$V2, temp$V1, type="l", main=mean(temp$V1), xlab=max(abs(temp$V1-mean(temp$V1)))/mean(temp$V1))
abline(a=mean(temp$V1),b=0, col="red")

temp <- read.csv("temp.csv", header=FALSE)
vel <- read.csv("velocity.csv", header=FALSE)
Temp <- mean(temp$V1[(nrow(temp)/2):nrow(temp)])
h <- hist(vel$V1)
par(new=TRUE)
curve(((1/(2*Temp*pi))^(3/2))*4*pi*(x^2)*exp(-(x^2)/(2*Temp)),
        0, max(h$breaks), ylab='', yaxt='n', col='red')
axis(1, col='red', col.axis='red')

plot(density(vel$V1))
