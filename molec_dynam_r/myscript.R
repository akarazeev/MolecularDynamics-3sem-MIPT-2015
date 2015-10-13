#plot3d()
#scatterplot3d()
#energy <- read.table("energy.csv", quote="\"")
library("scatterplot3d", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("rgl", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")

energy <- read.csv("energy.csv", header=FALSE)
kinetic <- read.csv("kinetic.csv", header=FALSE)
poten <- read.csv("poten.csv", header=FALSE)

plot(energy$V2, energy$V1, type="l")
par(new=TRUE)
plot(kinetic$V2, kinetic$V1, type="l", axes=FALSE)
par(new=TRUE)
plot(poten$V2, poten$V1, type="l")

temp <- read.csv("temp.csv", header=FALSE)
plot(temp$V2[100:nrow(temp)], temp$V1[100:nrow(temp)], type="l", main=mean(temp$V1[100:nrow(temp)]), xlab=max(abs(temp$V1[100:nrow(temp)]-mean(temp$V1[100:nrow(temp)])))/mean(temp$V1[100:nrow(temp)]))
abline(a=mean(temp$V1[100:nrow(temp)]),b=0, col="red")

plot(energy$V2, energy$V1, type="l", xlab="", ylab="", ylim=range(c(energy$V1,kinetic$V1,poten$V1)))
par(new=TRUE)
plot(kinetic$V2, kinetic$V1, type="l", xlab="", ylab="", axes=FALSE, ylim=range(c(energy$V1,kinetic$V1,poten$V1)))
par(new=TRUE)
plot(poten$V2, poten$V1, type="l", xlab="", ylab="", axes=FALSE, ylim=range(c(energy$V1,kinetic$V1,poten$V1)))

#plot()
plot(energy$V2, energy$V1, type="l", main=max(abs(energy$V1))/mean(energy$V1))
abline(a=mean(energy$V1),b=0)
print(max(energy$V1)/mean(energy$V1))

#scatterplot3d()
scatterplot3d(coords0$V1,coords0$V2,coords0$V3)
par(new=TRUE)
scatterplot3d(coords1$V1,coords1$V2,coords1$V3)
par(new=TRUE)
scatterplot3d(coords2$V1,coords2$V2,coords2$V3)

#plot3d()
#open3d()
open3d()
plot3d(coords0$V1,coords0$V2,coords0$V3, add=TRUE, col="red", decorate3d(box=TRUE))
plot3d(coords1$V1,coords1$V2,coords1$V3, add=TRUE)
plot3d(coords2$V1,coords2$V2,coords2$V3, add=TRUE, col="green")

plot3d(coords0$V1,coords0$V2,coords0$V3, col="red")
plot3d(coords1$V1,coords1$V2,coords1$V3)
plot3d(coords2$V1,coords2$V2,coords2$V3, col="green")

plot3d(allcoords$V1,allcoords$V2,allcoords$V3)

# Velocity
vel <- read.csv("velocity.csv", header=FALSE)
hist(vel$V1, breaks=10)
par(new=TRUE)
plot(density(vel$V1))

