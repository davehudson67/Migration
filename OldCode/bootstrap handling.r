###handling bootstrapped phylogenetic analysis of Pr(decline)~migratory

###which model to use?
m1<-readRDS("FitPhyloHemFlyEOO_5000.RDS")
m2<-readRDS("fit_phyloHemFlyEOO_H1F1_500.RDS")
m1.glm<-readRDS("glm_binomial_same_structure.RDS")

#trying to get the names of the parameters, to guide the arithmetic below
colnames(m1$bootstrap)
colnames(m2$bootstrap)

#checking associations between alpha, EOO and um something else
plot(m1$bootstrap[,32]~m1$bootstrap[,17],xlab=colnames(m1$bootstrap)[17],ylab=colnames(m1$bootstrap)[32])
plot(m1$bootstrap[,32]~m1$bootstrap[,2],xlab=colnames(m1$bootstrap)[2],ylab=colnames(m1$bootstrap)[32])
plot(m1$bootstrap[,2]~m1$bootstrap[,17],xlab=colnames(m1$bootstrap)[17],ylab=colnames(m1$bootstrap)[2])

#failed attempt to do this using matrix multiplication
#GLM.binary<-numeric(32)
#names(GLM.binary)<-colnames(m1$bootstrap)
#m1$bootstrap[1,]%*%GLM.binary

###so, doing it by hand. Get bootstrap density for each combination of Hemisphere, Flyway, Migratory
###Note that these estimates are ALL for a species with log(EOO)_centred of ZERO (i.e. an average breeding range)
###This first set is for 5000 bootstrap of H2F1 reference
#reference category
H2F1M0<-m1$bootstrap[,1]
#reference plus effect of being H1F6 
H1F6M0<-apply(m1$bootstrap[,c(1,3)],1,sum)
#reference plus effect of being H3F6
H3F6M0<-apply(m1$bootstrap[,c(1,4)],1,sum)
#etc...
H1F2M0<-apply(m1$bootstrap[,c(1,5)],1,sum)
H2F2M0<-apply(m1$bootstrap[,c(1,6)],1,sum)
H3F2M0<-apply(m1$bootstrap[,c(1,7)],1,sum)
H1F7M0<-apply(m1$bootstrap[,c(1,8)],1,sum)
H3F7M0<-apply(m1$bootstrap[,c(1,9)],1,sum)
H1F5M0<-apply(m1$bootstrap[,c(1,10)],1,sum)
H3F5M0<-apply(m1$bootstrap[,c(1,11)],1,sum)
H1F1M0<-apply(m1$bootstrap[,c(1,12)],1,sum)
H3F1M0<-apply(m1$bootstrap[,c(1,13)],1,sum)
H1F4M0<-apply(m1$bootstrap[,c(1,14)],1,sum)
H2F4M0<-apply(m1$bootstrap[,c(1,15)],1,sum)
H3F4M0<-apply(m1$bootstrap[,c(1,16)],1,sum)
#now add effect of being migratory (for reference level)
H2F1M1<-apply(m1$bootstrap[,c(1,2)],1,sum)
#add effect of being migratory AND effect of being migratory AND of being migratory in H1F6
H1F6M1<-apply(m1$bootstrap[,c(1,2,3,18)],1,sum)
#etc...
H3F6M1<-apply(m1$bootstrap[,c(1,2,4,19)],1,sum)
H1F2M1<-apply(m1$bootstrap[,c(1,2,5,20)],1,sum)
H2F2M1<-apply(m1$bootstrap[,c(1,2,6,21)],1,sum)
H3F2M1<-apply(m1$bootstrap[,c(1,2,7,22)],1,sum)
H1F7M1<-apply(m1$bootstrap[,c(1,2,8,23)],1,sum)
H3F7M1<-apply(m1$bootstrap[,c(1,2,9,24)],1,sum)
H1F5M1<-apply(m1$bootstrap[,c(1,2,10,25)],1,sum)
H3F5M1<-apply(m1$bootstrap[,c(1,2,11,26)],1,sum)
H1F1M1<-apply(m1$bootstrap[,c(1,2,12,27)],1,sum)
H3F1M1<-apply(m1$bootstrap[,c(1,2,13,28)],1,sum)
H1F4M1<-apply(m1$bootstrap[,c(1,2,14,29)],1,sum)
H2F4M1<-apply(m1$bootstrap[,c(1,2,15,30)],1,sum)
H3F4M1<-apply(m1$bootstrap[,c(1,2,16,31)],1,sum)

#and this second set is 1 500 bootstrap of H1F1 reference
H1F1M0b<-m2$bootstrap[,1]
H1F6M0b<-apply(m2$bootstrap[,c(1,3)],1,sum)
H3F6M0b<-apply(m2$bootstrap[,c(1,4)],1,sum)
H1F2M0b<-apply(m2$bootstrap[,c(1,5)],1,sum)
H2F2M0b<-apply(m2$bootstrap[,c(1,6)],1,sum)
H3F2M0b<-apply(m2$bootstrap[,c(1,7)],1,sum)
H1F7M0b<-apply(m2$bootstrap[,c(1,8)],1,sum)
H3F7M0b<-apply(m2$bootstrap[,c(1,9)],1,sum)
H1F5M0b<-apply(m2$bootstrap[,c(1,10)],1,sum)
H3F5M0b<-apply(m2$bootstrap[,c(1,11)],1,sum)
H2F1M0b<-apply(m2$bootstrap[,c(1,12)],1,sum)
H3F1M0b<-apply(m2$bootstrap[,c(1,13)],1,sum)
H1F4M0b<-apply(m2$bootstrap[,c(1,14)],1,sum)
H2F4M0b<-apply(m2$bootstrap[,c(1,15)],1,sum)
H3F4M0b<-apply(m2$bootstrap[,c(1,16)],1,sum)
H1F1M1b<-apply(m2$bootstrap[,c(1,2)],1,sum)
H1F6M1b<-apply(m2$bootstrap[,c(1,2,3,18)],1,sum)
H3F6M1b<-apply(m2$bootstrap[,c(1,2,4,19)],1,sum)
H1F2M1b<-apply(m2$bootstrap[,c(1,2,5,20)],1,sum)
H2F2M1b<-apply(m2$bootstrap[,c(1,2,6,21)],1,sum)
H3F2M1b<-apply(m2$bootstrap[,c(1,2,7,22)],1,sum)
H1F7M1b<-apply(m2$bootstrap[,c(1,2,8,23)],1,sum)
H3F7M1b<-apply(m2$bootstrap[,c(1,2,9,24)],1,sum)
H1F5M1b<-apply(m2$bootstrap[,c(1,2,10,25)],1,sum)
H3F5M1b<-apply(m2$bootstrap[,c(1,2,11,26)],1,sum)
H2F1M1b<-apply(m2$bootstrap[,c(1,2,12,27)],1,sum)
H3F1M1b<-apply(m2$bootstrap[,c(1,2,13,28)],1,sum)
H1F4M1b<-apply(m2$bootstrap[,c(1,2,14,29)],1,sum)
H2F4M1b<-apply(m2$bootstrap[,c(1,2,15,30)],1,sum)
H3F4M1b<-apply(m2$bootstrap[,c(1,2,16,31)],1,sum)

#and for the GLM version
H2F1M0.glm<-sum(m1.glm$coefficients[c(1)])
H1F6M0.glm<-sum(m1.glm$coefficients[c(1,3)])
H3F6M0.glm<-sum(m1.glm$coefficients[c(1,4)])
H1F2M0.glm<-sum(m1.glm$coefficients[c(1,5)])
H2F2M0.glm<-sum(m1.glm$coefficients[c(1,6)])
H3F2M0.glm<-sum(m1.glm$coefficients[c(1,7)])
H1F7M0.glm<-sum(m1.glm$coefficients[c(1,8)])
H3F7M0.glm<-sum(m1.glm$coefficients[c(1,9)])
H1F5M0.glm<-sum(m1.glm$coefficients[c(1,10)])
H3F5M0.glm<-sum(m1.glm$coefficients[c(1,11)])
H1F1M0.glm<-sum(m1.glm$coefficients[c(1,12)])
H3F1M0.glm<-sum(m1.glm$coefficients[c(1,13)])
H1F4M0.glm<-sum(m1.glm$coefficients[c(1,14)])
H2F4M0.glm<-sum(m1.glm$coefficients[c(1,15)])
H3F4M0.glm<-sum(m1.glm$coefficients[c(1,16)])
H2F1M1.glm<-sum(m1.glm$coefficients[c(1,2)])
H1F6M1.glm<-sum(m1.glm$coefficients[c(1,2,3,18)])
H3F6M1.glm<-sum(m1.glm$coefficients[c(1,4,2,19)])
H1F2M1.glm<-sum(m1.glm$coefficients[c(1,5,2,20)])
H2F2M1.glm<-sum(m1.glm$coefficients[c(1,6,2,21)])
H3F2M1.glm<-sum(m1.glm$coefficients[c(1,7,2,22)])
H1F7M1.glm<-sum(m1.glm$coefficients[c(1,8,2,23)])
H3F7M1.glm<-sum(m1.glm$coefficients[c(1,9,2,24)])
H1F5M1.glm<-sum(m1.glm$coefficients[c(1,10,2,25)])
H3F5M1.glm<-sum(m1.glm$coefficients[c(1,11,2,26)])
H1F1M1.glm<-sum(m1.glm$coefficients[c(1,12,2,27)])
H3F1M1.glm<-sum(m1.glm$coefficients[c(1,13,2,28)])
H1F4M1.glm<-sum(m1.glm$coefficients[c(1,14,2,29)])
H2F4M1.glm<-sum(m1.glm$coefficients[c(1,15,2,30)])
H3F4M1.glm<-sum(m1.glm$coefficients[c(1,16,2,31)])


#calculate the migration-effect for each Hem-Fly category...this is if we just want to plot migration differences
H1F6diff<-H1F6M1-H1F6M0
H3F6diff<-H3F6M1-H3F6M0
H1F2diff<-H1F2M1-H1F2M0
H2F2diff<-H2F2M1-H2F2M0
H3F2diff<-H3F2M1-H3F2M0
H1F7diff<-H1F7M1-H1F7M0
H3F7diff<-H3F7M1-H3F7M0
H1F5diff<-H1F5M1-H1F5M0
H3F5diff<-H3F5M1-H3F5M0
H1F1diff<-H1F1M1-H1F1M0
H3F1diff<-H3F1M1-H3F1M0
H1F4diff<-H1F4M1-H1F4M0
H2F4diff<-H2F4M1-H2F4M0
H3F4diff<-H3F4M1-H3F4M0
H2F1diff<-H2F1M1-H2F1M0

###Now the big narly code to plot the bootstrap posteriors for each Hem-Fly combination
opar<-par(mfcol=c(3,6))
residents<-H1F1M0
migrants<-H1F1M1
plot(density(residents),main="American",xlab="H1F1",ylab="Northern",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H1F1M0b),col=rgb(0,0,1,0.5))
polygon(density(H1F1M1b),col=rgb(1,0,0,0.5))
lines(c(H1F1M0.glm,H1F1M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H1F1M1.glm,H1F1M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

residents<-H3F1M0
migrants<-H3F1M1
plot(density(residents),main="",xlab="H3F1",ylab="Equatorial",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H3F1M0b),col=rgb(0,0,1,0.5))
polygon(density(H3F1M1b),col=rgb(1,0,0,0.5))
lines(c(H3F1M0.glm,H3F1M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H3F1M1.glm,H3F1M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

residents<-H2F1M0
migrants<-H2F1M1
plot(density(residents),main="",xlab="H2F1",ylab="Southern",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H2F1M0b),col=rgb(0,0,1,0.5))
polygon(density(H2F1M1b),col=rgb(1,0,0,0.5))
lines(c(H2F1M0.glm,H2F1M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H2F1M1.glm,H2F1M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

residents<-H1F2M0
migrants<-H1F2M1
plot(density(residents),main="African",xlab="H1F2",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H1F2M0b),col=rgb(0,0,1,0.5))
polygon(density(H1F2M1b),col=rgb(1,0,0,0.5))
lines(c(H1F2M0.glm,H1F2M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H1F2M1.glm,H1F2M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

residents<-H3F2M0
migrants<-H3F2M1
plot(density(residents),main="",xlab="H3F2",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H3F2M0b),col=rgb(0,0,1,0.5))
polygon(density(H3F2M1b),col=rgb(1,0,0,0.5))
lines(c(H3F2M0.glm,H3F2M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H3F2M1.glm,H3F2M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

residents<-H2F2M0
migrants<-H2F2M1
plot(density(residents),main="",xlab="H2F2",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H2F2M0b),col=rgb(0,0,1,0.5))
polygon(density(H2F2M1b),col=rgb(1,0,0,0.5))
lines(c(H2F2M0.glm,H2F2M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H2F2M1.glm,H2F2M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

residents<-H1F4M0
migrants<-H1F4M1
plot(density(residents),main="Asian",xlab="H1F4",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H1F4M0b),col=rgb(0,0,1,0.5))
polygon(density(H1F4M1b),col=rgb(1,0,0,0.5))
lines(c(H1F4M0.glm,H1F4M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H1F4M1.glm,H1F4M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

residents<-H3F4M0
migrants<-H3F4M1
plot(density(residents),main="",xlab="H3F4",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H3F4M0b),col=rgb(0,0,1,0.5))
polygon(density(H3F4M1b),col=rgb(1,0,0,0.5))
lines(c(H3F4M0.glm,H3F4M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H3F4M1.glm,H3F4M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

residents<-H2F4M0
migrants<-H2F4M1
plot(density(residents),main="",xlab="H2F4",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H2F4M0b),col=rgb(0,0,1,0.5))
polygon(density(H2F4M1b),col=rgb(1,0,0,0.5))
lines(c(H2F4M0.glm,H2F4M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H2F4M1.glm,H2F4M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

residents<-H1F5M0
migrants<-H1F5M1
plot(density(residents),main="Amer-Asian",xlab="H1F5",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H1F5M0b),col=rgb(0,0,1,0.5))
polygon(density(H1F5M1b),col=rgb(1,0,0,0.5))
lines(c(H1F5M0.glm,H1F5M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H1F5M1.glm,H1F5M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

residents<-H3F5M0
migrants<-H3F5M1
plot(density(residents),main="",xlab="H3F5",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H3F5M0b),col=rgb(0,0,1,0.5))
polygon(density(H3F5M1b),col=rgb(1,0,0,0.5))
lines(c(H3F5M0.glm,H3F5M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H3F5M1.glm,H3F5M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

#blank plots for missing regions
plot.new()

residents<-H1F6M0
migrants<-H1F6M1
plot(density(residents),main="Afri-Asian",xlab="H1F6",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H1F6M0b),col=rgb(0,0,1,0.5))
polygon(density(H1F6M1b),col=rgb(1,0,0,0.5))
lines(c(H1F6M0.glm,H1F6M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H1F6M1.glm,H1F6M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

residents<-H3F6M0
migrants<-H3F6M1
plot(density(residents),main="",xlab="H3F6",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H3F6M0b),col=rgb(0,0,1,0.5))
polygon(density(H3F6M1b),col=rgb(1,0,0,0.5))
lines(c(H3F6M0.glm,H3F6M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H3F6M1.glm,H3F6M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

plot.new()

residents<-H1F7M0
migrants<-H1F7M1
plot(density(residents),main="3 Flyways",xlab="H1F7",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H1F7M0b),col=rgb(0,0,1,0.5))
polygon(density(H1F7M1b),col=rgb(1,0,0,0.5))
lines(c(H1F7M0.glm,H1F7M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H1F7M1.glm,H1F7M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

residents<-H3F7M0
migrants<-H3F7M1
plot(density(residents),main="",xlab="H3F7",xlim=c(min(c(residents,migrants)),max(c(residents,migrants))))
polygon(density(residents),col=rgb(0,0,1,0.5))
polygon(density(migrants),col=rgb(1,0,0,0.5))
polygon(density(H3F7M0b),col=rgb(0,0,1,0.5))
polygon(density(H3F7M1b),col=rgb(1,0,0,0.5))
lines(c(H3F7M0.glm,H3F7M0.glm),c(0,10),col=rgb(0,0,1,0.5),lwd=3)
lines(c(H3F7M1.glm,H3F7M1.glm),c(0,10),col=rgb(1,0,0,0.5),lwd=3)

plot.new()

par(opar)

####bit of a look at oddball bootstraps
which(m1.boot100.H2F1S$bootstrap[,32]<0.1)
ddd<-apply(m1.boot100.H2F1S$bootdata[,which(m1.boot100.H2F1S$bootstrap[,32]<0.1)],1,sum)
plot(density(resid(m1.boot100.H2F1S)))
dddd<-resid(m1.boot100.H2F1S)[which(ddd==0)]
points(rep(0.2,length(dddd))~dddd)


