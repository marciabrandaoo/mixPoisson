#################################################
# Breast cancer - cure rate MPIGcr and PTcr     #
#         Tables: 5, 6 and 7                    #
#       Figures: 2, 3, 4 and 5                  #
#################################################


#database-----------------------------------------
data_BC <- read.csv("data_breast_cancer.csv",sep = ",")  
head(data_BC)

#===================================================
# Breast cancer  
# Figure 2: Estimated survival function (SF) obtained from the KM estimator for overall patients diagnosed with breast
#           cancer, by clinical stage, surgery, radiotherapy, chemotherapy, and combinations of treatments
#===================================================
pdf("Figure2.pdf", width = 16, height = 8.5) 
par(mfrow=c(2,3))
KM0 <- survfit(Surv(tempo, status) ~ 1, data=data_BC, conf.int=FALSE)
plot(KM0 ,xlab="Time (years)", las=1,cex.axis=1.2,cex.lab=1.2,ylab="Survival Probabily",conf.int=FALSE)
KM1 <- survfit(Surv(tempo, status) ~ ECGRUP, data=data_BC, conf.int=FALSE)
plot(KM1 ,lty = 1:4, col=1:4, las=1,cex.axis=1.2,cex.lab=1.2,xlab="Time (years)",ylab="Survival Probabily",conf.int=FALSE)
legend("bottomleft", title="Stages" ,legend=c("Stage I", "Stage II","Stage III","Stage IV"),lty = 1:4, cex=1.2,  col=1:4, bty="n")
KM2 <- survfit(Surv(tempo, status) ~ CIRURGIA, data=data_BC)
plot(KM2, lty = c(1,5), col=c(2,3), las=1,cex.axis=1.2,cex.lab=1.2,xlab="Time (years)",ylab="Survival Probabily")
legend("bottomleft", title="Surgery" ,legend=c("No", "Yes"),lty = c(1,5), cex=1.2,col=c(2,3),bty="n")
KM3 <- survfit(Surv(tempo, status) ~ RADIO, data=data_BC)
plot(KM3, lty = c(1,5),col=c(2,3),las=1,cex.axis=1.2,cex.lab=1.2,xlab="Time (years)",ylab="Survival Probabily")
legend("bottomleft",  title="Radiotherapy" ,legend=c("No", "Yes"),lty = c(1,5),cex=1.2, col=c(2,3), bty="n")
KM4 <- survfit(Surv(tempo, status) ~ QUIMIO, data=data_BC)
plot(KM4, lty = c(1,5),col=c(2,3),las=1,cex.axis=1.2,cex.lab=1.2,xlab="Time (years)",ylab="Survival Probabily")
legend("bottomleft",  title="Chemotherapy" ,legend=c("No", "Yes"),lty = c(1,5),cex=1.2,col=c(2,3), bty="n")
# iterations of treatments
vari1 = data_BC$CIRURGIA*data_BC$RADIO
vari2 = data_BC$CIRURGIA*data_BC$QUIMIO
vari3 = data_BC$RADIO*data_BC$QUIMIO
vari4 = data_BC$CIRURGIA*data_BC$RADIO*data_BC$QUIMIO
KM50 <- survfit(Surv(tempo, status) ~ vari1+vari2+vari3+vari4, data=data_BC)
plot(KM50, lty = c(1:5),col=c(2,3,4,6,7,8),las=1,cex.axis=1.2,cex.lab=1.2,xlab="Time (years)",ylab="Survival Probabily")
legend("bottomleft",legend=c('Overall', 'Surg+Radio','Surg+Chemo','Radio+Chemo','Surg+Radio+Chemo'),lty = c(1,5), cex=1.2,col=c(2,3,4,6,7,8),bty="n")
dev.off()


 
#===================================================
# Breast cancer  
# Table 5: ML estimate, standard error, and respective of z-value and p-value obtained by 
#           fitting of mixture of the  MPIGcr/WEI.
#===================================================

tab_application_MIG = cbind(round(estimates.psi.MIG,4),c(round(matcoef_MIG[,3],4),"-","-","-"),c(round(matcoef_MIG[,4],4),"-","-","-"))
tab_application_MIG = cbind(
  round(estimates.psi.MIG, 4),
  c(round(matcoef_MIG[,3], 4), NA, NA, NA),
  c(round(matcoef_MIG[,4], 4), NA, NA, NA)
)
colnames(tab_application_MIG) <- c("Estimate_MIG", "SE_hessian_MIG", "z_value", "p_value")
tab_print <- apply(tab_application_MIG, 2, function(x) ifelse(is.na(x), "-", format(x, digits = 4)))
# Mostrar sem aspas
print(tab_print, quote = FALSE)

name2 = paste("Table5")
write.table(tab_print,file=paste(name2,".txt"))
write.csv(tab_print,file=paste(name2,".csv"))



#===================================================
# Breast cancer  
# Figure 3:  Normalized randomized quantile residuals for MPIGcr case applied to the breast cancer dataset. Estimated SF
#       and Hazard function for the MPIGcr model for patients 20, 60 years old who underwent radiotherapy and chemotherapy
#       through different stages of the disease: Stage I (black), Stage II (red), Stage III (green), and Stage IV (blue).
#===================================================


#-----------------------------------
# Normalized randomized quantile
#-----------------------------------
auxiliar <-S_f_hpop(psi.MIG,Data,dist=1,model=1,real=FALSE)$Spop
##-------------
Spop1 <- auxiliar
rq=RQ.res(Spop1)
xy.min=min(c(rq))
xy.max=max(c(rq))
ks_test  <- ks.test(rq, "pnorm")
ad_test  <- ad.test(rq)
cvm_test <- cvm.test(rq)
envelopeDS(rq)
#-----------------------------------
# Estimated SF and Hazard function for the MPIGcr model 
#-----------------------------------
#profile 1
#Estad3
#------------------------Stages of disease
# ECGRUP = 1 (beta0 intercepto)
# ECGRUP = 2 (beta1)
# ECGRUP = 3 (betdata_BC)
# ECGRUP = 4 (beta3)
#-----------------------Treatments
# Surgery                     = A (beta0 intercepto)
# Radiotherapy                = B (beta4)
# Quimiotherapy               = C (beta5)
# Surgery + Radiotherapy      = D (beta6)
# Surgery + Quimiotherapy     = E (beta7)
# Radiotherapy + Quimiotherapy= F (beta8)

#--------------------------------
# Graphich for Estimated Suvivel fuction 
#----------------------------------
psi1_MIG <- psi.MIG
r= length(psi.MIG[1:8])
n=length(t)
SpopMIG11<-c();SpopMIG12<-c();SpopMIG13<-c();SpopMIG14<-c();
SpopMIG21<-c();SpopMIG22<-c();SpopMIG23<-c();SpopMIG24<-c();
SpopMIG31<-c();SpopMIG32<-c();SpopMIG33<-c();SpopMIG34<-c();

hpopMIG11<-c();hpopMIG12<-c();hpopMIG13<-c();hpopMIG14<-c();
hpopMIG21<-c();hpopMIG22<-c();hpopMIG23<-c();hpopMIG24<-c();
hpopMIG31<-c();hpopMIG32<-c();hpopMIG33<-c();hpopMIG34<-c();
#-----------------------------
# 20 years
#----------------------------------
# scenario 1.1 -  20 years; Stage 1, radio and quimio=YES
x11    = c(1,0,0,0,0,1,1,20)
t.seq<-seq(0,15,length=n)
SpopMIG1<-c() 
for(i in 1:length(t.seq))
{
  x = x11
  data.t<- list(t=t.seq[i], delta=1,x=x11)
  SpopMIG11[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$Spop
  hpopMIG11[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$hpop
}
#----------------------------------
# scenario 1.2 - 20 years; Stage 2, radio and quimio=YES
x12    = c(1,1,0,0,0,1,1,20)
SpopMIG<-c(); 
for(i in 1:length(t.seq))
{
  x = x12
  data.t<- list(t=t.seq[i], delta=1,x=x12)
  SpopMIG12[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$Spop
  hpopMIG12[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$hpop
}
#----------------------------------
# scenario 1.3
# 20 years; Stage 3, radio and quimio=YES
x13    = c(1,0,1,0,0,1,1,20)
for(i in 1:length(t.seq))
{
  x = x13
  data.t<- list(t=t.seq[i], delta=1,x=x13)
  SpopMIG13[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$Spop
  hpopMIG13[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$hpop
}
# scenario 1.4 - 20 years; Stage 4, Radio and quimio=YES
x14    = c(1,0,0,1,0,1,1,20)
for(i in 1:length(t.seq))
{
  x = x14
  data.t<- list(t=t.seq[i], delta=1,x=x14)
  SpopMIG14[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$Spop
  hpopMIG14[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$hpop
}

#---------------------------------------
# 60 years
#----------------------------------
# scenario 3.1  -  60years; Stage 1, radio and quimio=YES
x31    = c(1,0,0,0,0,1,1,60)
for(i in 1:length(t.seq))
{
  x = x31
  data.t<- list(t=t.seq[i], delta=1,x=x31)
  SpopMIG31[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$Spop
  hpopMIG31[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$hpop
}
#----------------------------------
# scenario 3.2  -  60 years; Stage 2, radio and quimio=YES
x32    = c(1,1,0,0,0,1,1,60,psi.MIG[9:11])
for(i in 1:length(t.seq))
{
  x = x32
  data.t<- list(t=t.seq[i], delta=1,x=x32)
  SpopMIG32[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$Spop
  hpopMIG32[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$hpop
}
#----------------------------------
# scenario 3.3 - 60 years; Stage 3, radio and quimio=YES
x33    = c(1,0,1,0,0,1,1,60,psi.MIG[9:11])
for(i in 1:length(t.seq))
{
  x = x33
  data.t<- list(t=t.seq[i], delta=1,x=x33)
  SpopMIG33[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$Spop
  hpopMIG33[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$hpop
}
#----------------------------------
# scenario 3.4 -  60 years; Stage 4, radio and quimio=YES
x34    = c(1,0,0,1,0,1,1,60,psi.MIG[9:11])
SpopMIG<-c();SpopPO<-c();SpopBN<-c() 
for(i in 1:length(t.seq))
{
  x = x34
  data.t<- list(t=t.seq[i], delta=1,x=x34)
  SpopMIG34[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$Spop
  hpopMIG34[i] <-S_f_hpop(psi1_MIG,data.t,dist=1,model=1)$hpop
}
#----------------------------------
w_in <- 1018/96
h_in <- 350/96
cairo_pdf("Figure3.pdf", width = w_in, height=h_in)
par(mfrow=c(1,3))
# graphic. ormalized randomized quantile residuals 
envelopeDS(rq)
# graphic. estimated survival function profile 1 and 2
plot(t.seq,SpopMIG11,ylim=c(0,1),type="l",lwd=2,las=1,lty=1,cex.axis=1,cex.lab=1.2,xlab="Time (years)",ylab="Estimated survival function")
points(t.seq,SpopMIG12,type="l",lwd=2,lty=1,col=2)
points(t.seq,SpopMIG13,type="l",lwd=2,lty=1,col=3)
points(t.seq,SpopMIG14,type="l",lwd=2,lty=1,col=4)
points(t.seq,SpopMIG31,type="l",lwd=2,lty=3,col=1)
points(t.seq,SpopMIG32,type="l",lwd=2,lty=3,col=2)
points(t.seq,SpopMIG33,type="l",lwd=2,lty=3,col=3)
points(t.seq,SpopMIG34,type="l",lwd=2,lty=3,col=4)
legend("bottomleft",c(expression(paste("Age 20"))
                      ,expression(paste("Age 60"))  
)
,pch="",lty=c(1,2,3),lwd=c(2,2,2),cex=1,inset=0.01,box.col="white",xjust=1)
#graphic. estimated hazard function profile 1 and 2
plot(t.seq,hpopMIG11,ylim=c(0,0.4),type="l",lwd=2,las=1,lty=1,cex.axis=1,cex.lab=1.2,xlab="Time (years)",ylab="Estimated hazard function")
points(t.seq,hpopMIG12,type="l",lwd=2,lty=1,col=2)
points(t.seq,hpopMIG13,type="l",lwd=2,lty=1,col=3)
points(t.seq,hpopMIG14,type="l",lwd=2,lty=1,col=4)
points(t.seq,hpopMIG31,type="l",lwd=2,lty=3,col=1)
points(t.seq,hpopMIG32,type="l",lwd=2,lty=3,col=2)
points(t.seq,hpopMIG33,type="l",lwd=2,lty=3,col=3)
points(t.seq,hpopMIG34,type="l",lwd=2,lty=3,col=4)
legend("topright",c(expression(paste("Age 20"))
                    ,expression(paste("Age 60"))  
)
,pch="",lty=c(1,2,3),lwd=c(2,2,2),cex=1,inset=0.01,box.col="white",xjust=1)

dev.off()





#===================================================
# Breast cancer  
# Estimated values of cure rate  for patients 
#===================================================
# matrix with estimated cure fractions 
Varcov = solve(hessian(llikeobserved,x0=psi,Data=Data,dist=1,model=1))[1:r,1:r]
psi.MIG     = matrix(psi,nrow=1);colnames(psi.MIG)= c(paste("beta",0:(r-1),sep=""),"alpha","nu","m");round(psi.MIG,6)
beta.est = c(psi.MIG[1:r])
m.est   <- psi.MIG[r+3]
#------------------------Stages of disease
# stage = 1 (beta0 intercept)
# stage = 2 (beta1) 
# stage = 3 (beta2)
# stage = 4 (beta3)
#-----------------------Treatments
# Surgery       = beta4
# radiotherapy  = beta5
# chemotherapy  = beta6
# Age           = beta 7 (20,56 or 60 years)
#---------------------------------------------------
# data vector
# x  = c(1,1,0,0,0,0,0,20) # 0 or 1 in beta0 to beta 6 to select 
# stages and treatments and select beta7 in c(20,56,60)
#---------------------------------------------------
# 20 years
#---------------------------------------------------
#1) No treatment  - Stage 1 - 20 years
x1    = c(1,0,0,0,0,0,0,20)  
theta1 <- exp(t(x1)%*%beta.est) 
p1 = m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta1)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta1)/(1-m.est^2)) ))

varp1 = p1*(1-p1)*sqrt(t(x1)%*%Varcov%*%(x1))
ICp1  = rbind(round(c( p1 - 1.96*varp1, p1 + 1.96*varp1),3));
round(p1,3);ICp1


#---------------------------------------------------
#2) No treatment  - Stage 2 - 20 years
x2    = c(1,1,0,0,0,0,0,20) 
theta2 <- exp(t(x2)%*%beta.est) 
p2 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta2)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta2)/(1-m.est^2)) ))
varp2 = p2*(1-p2)*sqrt(t(x2)%*%Varcov%*%(x2))
ICp2  = rbind(round(c( p2 - 1.96*varp2, p2 + 1.96*varp2),3));
round(p2,3);ICp2
#---------------------------------------------------
#3) No treatment  - Stage 3 - 20 years
x3    = c(1,0,1,0,0,0,0,20) 
theta3 <- exp(t(x3)%*%beta.est) 
p3 =   m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta3)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta3)/(1-m.est^2)) ))

varp3 = p3*(1-p3)*sqrt(t(x3)%*%Varcov%*%(x3))
ICp3  = rbind(round(c( p3 - 1.96*varp3, p3 + 1.96*varp3),3));
round(p3,3);ICp3
#---------------------------------------------------
#4) No treatment  - Stage 4 - 20 years
x4    = c(1,0,0,1,0,0,0,20) 
theta4 <- exp(t(x4)%*%beta.est) 
p4 =   m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta4)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta4)/(1-m.est^2)) ))

varp4 = p4*(1-p4)*sqrt(t(x4)%*%Varcov%*%(x4))
ICp4  = rbind(round(c( p4 - 1.96*varp4, p4 + 1.96*varp4),3));
round(p4,3);ICp4
#---------------------------------------------------
#5) Surgery  - Stage 1 - 20 years
x5    = c(1,0,0,0,1,0,0,20) 
theta5 <- exp(t(x5)%*%beta.est) 
p5 = m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta5)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta5)/(1-m.est^2)) ))

varp5 = p5*(1-p5)*sqrt(t(x5)%*%Varcov%*%(x5))
ICp5  = rbind(round(c( p5 - 1.96*varp5, p5 + 1.96*varp5),3));
round(p5,3);ICp5
#---------------------------------------------------
#6) Surgery  - Stage 2 - 20 years
x6    = c(1,1,0,0,1,0,0,20) 
theta6 <- exp(t(x6)%*%beta.est) 
p6 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta6)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta6)/(1-m.est^2)) ))

varp6 = p6*(1-p6)*sqrt(t(x6)%*%Varcov%*%(x6))
ICp6  = rbind(round(c( p6 - 1.96*varp6, p6 + 1.96*varp6),3));
round(p6,3);ICp6
#---------------------------------------------------
# 7) Surgery  - Stage 3 - 20 years
x7    = c(1,0,1,0,1,0,0,20)  
theta7 <- exp(t(x7)%*%beta.est) 
p7 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta7)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta7)/(1-m.est^2)) ))

varp7 = p7*(1-p7)*sqrt(t(x7)%*%Varcov%*%(x7))
ICp7  = rbind(round(c( p7 - 1.96*varp7, p7 + 1.96*varp7),3));
round(p7,3);ICp7
#---------------------------------------------------
#8) Surgery  - Stage 4 - 20 years
x8    = c(1,0,0,1,1,0,0,20)  
theta8 <- exp(t(x8)%*%beta.est) 
p8 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta8)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta8)/(1-m.est^2)) ))

varp8 = p8*(1-p8)*sqrt(t(x8)%*%Varcov%*%(x8))
ICp8  = rbind(round(c( p8 - 1.96*varp8, p8 + 1.96*varp8),3));
round(p8,3);ICp8
#---------------------------------------------------
# 9) radiotheraphy  - Stage 1 - 20 years
x9    = c(1,0,0,0,0,1,0,20)  
theta9 <- exp(t(x9)%*%beta.est) 
p9 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta9)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta9)/(1-m.est^2)) ))

varp9 = p9*(1-p9)*sqrt(t(x9)%*%Varcov%*%(x9))
ICp9  = rbind(round(c( p9 - 1.96*varp9, p9 + 1.96*varp9),3));
round(p9,3);ICp9
#---------------------------------------------------
#10 radiotheraphy  - Stage 2 - 20 years
x10    = c(1,1,0,0,0,1,0,20)  
theta10 <- exp(t(x10)%*%beta.est) 
p10 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta10)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta10)/(1-m.est^2)) ))

varp10 = p10*(1-p10)*sqrt(t(x10)%*%Varcov%*%(x10))
ICp10  = rbind(round(c( p10 - 1.96*varp10, p10 + 1.96*varp10),3));
round(p10,3);ICp10
#---------------------------------------------------
#11) radiotheraphy  - Stage 3 - 20 years
x11    = c(1,0,1,0,0,1,0,20)  
theta11 <- exp(t(x11)%*%beta.est) 
p11 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta11)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta11)/(1-m.est^2)) ))

varp11 = p11*(1-p11)*sqrt(t(x11)%*%Varcov%*%(x11))
ICp11  = rbind(round(c( p11 - 1.96*varp11, p11 + 1.96*varp11),3));
round(p11,3);ICp11
#---------------------------------------------------
#12) radiotheraphy  - Stage 4 - 20 years
x12    = c(1,0,0,1,0,1,0,20)  
theta12 <- exp(t(x12)%*%beta.est) 
p12 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta12)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta12)/(1-m.est^2)) ))

varp12 = p12*(1-p12)*sqrt(t(x12)%*%Varcov%*%(x12))
ICp12  = rbind(round(c( p12 - 1.96*varp12, p12 + 1.96*varp12),3));
round(p12,3);ICp12
#---------------------------------------------------
#13 Chemotheraphy  - Stage 1 - 20 years
x13    = c(1,0,0,0,0,0,1,20)  
theta13 <- exp(t(x13)%*%beta.est) 
p13 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta13)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta13)/(1-m.est^2)) ))

varp13 = p13*(1-p13)*sqrt(t(x13)%*%Varcov%*%(x13))
ICp13  = rbind(round(c( p13 - 1.96*varp13, p13 + 1.96*varp13),3));
round(p13,3);ICp13
#---------------------------------------------------
#14 Chemotheraphy  - Stage 2 - 20 years
x14    = c(1,1,0,0,0,0,1,20)  
theta14 <- exp(t(x14)%*%beta.est) 
p14 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta14)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta14)/(1-m.est^2)) ))

varp14 = p14*(1-p14)*sqrt(t(x14)%*%Varcov%*%(x14))
ICp14  = rbind(round(c( p14 - 1.96*varp14, p14 + 1.96*varp14),3));
round(p14,3);ICp14
#---------------------------------------------------
#15) Chemotheraphy  - Stage 3 - 20 years
x15    = c(1,0,1,0,0,0,1,20)  
theta15 <- exp(t(x15)%*%beta.est) 
p15 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta15)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta15)/(1-m.est^2)) ))

varp15 = p15*(1-p15)*sqrt(t(x15)%*%Varcov%*%(x15))
ICp15  = rbind(round(c( p15 - 1.96*varp15, p15 + 1.96*varp15),3));
round(p15,3);ICp15
#---------------------------------------------------
#16) Chemotheraphy  - Stage 4 - 20 years
x16    = c(1,0,0,1,0,0,1,20)  
theta16 <- exp(t(x16)%*%beta.est) 
p16 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta16)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta16)/(1-m.est^2)) ))

varp16 = p16*(1-p16)*sqrt(t(x16)%*%Varcov%*%(x16))
ICp16  = rbind(round(c( p16 - 1.96*varp16, p16 + 1.96*varp16),3));
round(p16,3);ICp16

#---------------------------------------------------
#17) Surgery + radiotheraphy  - Stage 1 - 20 years
x17    = c(1,0,0,0,1,1,0,20)  
theta17 <- exp(t(x17)%*%beta.est) 
p17 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta17)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta17)/(1-m.est^2)) ))

varp17 = p17*(1-p17)*sqrt(t(x17)%*%Varcov%*%(x17))
ICp17  = rbind(round(c( p17 - 1.96*varp17, p17 + 1.96*varp17),3));
round(p17,3);ICp17
#---------------------------------------------------
#18) Surgery + radiotheraphy  - Stage 2 - 20 years
x18    = c(1,1,0,0,1,1,0,20)  
theta18 <- exp(t(x18)%*%beta.est) 
p18 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta18)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta18)/(1-m.est^2)) ))

varp18 = p18*(1-p18)*sqrt(t(x18)%*%Varcov%*%(x18))
ICp18  = rbind(round(c( p18 - 1.96*varp18, p18 + 1.96*varp18),3));
round(p18,3);ICp18
#---------------------------------------------------
#19) Surgery + radiotheraphy  - Stage 3 - 20 years
x19    = c(1,0,1,0,1,1,0,20)  
theta19 <- exp(t(x19)%*%beta.est) 
p19 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta19)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta19)/(1-m.est^2)) ))

varp19 = p19*(1-p19)*sqrt(t(x19)%*%Varcov%*%(x19))
ICp19  = rbind(round(c( p19 - 1.96*varp19, p19 + 1.96*varp19),3));
round(p19,3);ICp19
#---------------------------------------------------
#20) Surgery + radiotheraphy  - Stage 4 - 20 years
x20    = c(1,0,0,1,1,1,0,20)  
theta20 <- exp(t(x20)%*%beta.est) 
p20 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta20)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta20)/(1-m.est^2)) ))

varp20 = p20*(1-p20)*sqrt(t(x20)%*%Varcov%*%(x20))
ICp20  = rbind(round(c( p20 - 1.96*varp20, p20 + 1.96*varp20),3));
round(p20,3);ICp20

#---------------------------------------------------
#21) Surgery + Chemotheraphy  - Stage 1 - 20 years
x21    = c(1,0,0,0,1,0,1,20)  
theta21 <- exp(t(x21)%*%beta.est) 
p21 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta21)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta21)/(1-m.est^2)) ))

varp21 = p21*(1-p21)*sqrt(t(x21)%*%Varcov%*%(x21))
ICp21  = rbind(round(c( p21 - 1.96*varp21, p21 + 1.96*varp21),3));
round(p21,3);ICp21
#---------------------------------------------------
#22) Surgery + Chemotheraphy  - Stage 2 - 20 years
x22    = c(1,1,0,0,1,0,1,20)  
theta22 <- exp(t(x22)%*%beta.est) 
p22 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta22)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta22)/(1-m.est^2)) ))

varp22 = p22*(1-p22)*sqrt(t(x22)%*%Varcov%*%(x22))
ICp22  = rbind(round(c( p22 - 1.96*varp22, p22 + 1.96*varp22),3));
round(p22,3);ICp22
#---------------------------------------------------
#23) Surgery + Chemotheraphy  - Stage 3 - 20 years
x23    = c(1,0,1,0,1,0,1,20)  
theta23 <- exp(t(x23)%*%beta.est) 
p23 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta23)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta23)/(1-m.est^2)) ))

varp23 = p23*(1-p23)*sqrt(t(x23)%*%Varcov%*%(x23))
ICp23  = rbind(round(c( p23 - 1.96*varp23, p23 + 1.96*varp23),3));
round(p23,3);ICp23
#---------------------------------------------------
#24) Surgery + Chemotheraphy  - Stage 4 - 20 years
x24    = c(1,0,0,1,1,0,1,20)  
theta24 <- exp(t(x24)%*%beta.est) 
p24 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta24)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta24)/(1-m.est^2)) ))

varp24 = p24*(1-p24)*sqrt(t(x24)%*%Varcov%*%(x24))
ICp24  = rbind(round(c( p24 - 1.96*varp24, p24 + 1.96*varp24),3));
round(p24,3);ICp24
#---------------------------------------------------
#25) radiotherapy + Chemotheraphy  - Stage 1 - 20 years
x25    = c(1,0,0,0,0,1,1,20)  
theta25 <- exp(t(x25)%*%beta.est) 
p25 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta25)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta25)/(1-m.est^2)) ))

varp25 = p25*(1-p25)*sqrt(t(x25)%*%Varcov%*%(x25))
ICp25  = rbind(round(c( p25 - 1.96*varp25, p25 + 1.96*varp25),3));
round(p25,3);ICp25
#---------------------------------------------------
#26) radiotherapy + Chemotheraphy  - Stage 2 - 20 years
x26    = c(1,1,0,0,0,1,1,20)  
theta26 <- exp(t(x26)%*%beta.est) 
p26 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta26)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta26)/(1-m.est^2)) ))

varp26 = p26*(1-p26)*sqrt(t(x26)%*%Varcov%*%(x26))
ICp26  = rbind(round(c( p26 - 1.96*varp26, p26 + 1.96*varp26),3));
round(p26,3);ICp26
#---------------------------------------------------
#27) radiotherapy + Chemotheraphy  - Stage 3 - 20 years
x27    = c(1,0,1,0,0,1,1,20)  
theta27 <- exp(t(x27)%*%beta.est) 
p27 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta27)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta27)/(1-m.est^2)) ))

varp27 = p27*(1-p27)*sqrt(t(x27)%*%Varcov%*%(x27))
ICp27  = rbind(round(c( p27 - 1.96*varp27, p27 + 1.96*varp27),3));
round(p27,3);ICp27
#---------------------------------------------------
#28) radiotherapy + Chemotheraphy  - Stage 4 - 20 years
x28    = c(1,0,0,1,0,1,1,20)  
theta28 <- exp(t(x28)%*%beta.est) 
p28 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta28)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta28)/(1-m.est^2)) ))

varp28 = p28*(1-p28)*sqrt(t(x28)%*%Varcov%*%(x28))
ICp28  = rbind(round(c( p28 - 1.96*varp28, p28 + 1.96*varp28),3));
round(p28,3);ICp28
#---------------------------------------------------
#29) Surgery + radiotherapy + Chemotheraphy  - Stage 1 - 20 years
x29    = c(1,0,0,0,1,1,1,20)  
theta29 <- exp(t(x29)%*%beta.est) 
p29 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta29)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta29)/(1-m.est^2)) ))

varp29 = p29*(1-p29)*sqrt(t(x29)%*%Varcov%*%(x29))
ICp29  = rbind(round(c( p29 - 1.96*varp29, p29 + 1.96*varp29),3));
round(p29,3);ICp29
#---------------------------------------------------
#30) Surgery +  radiotherapy + Chemotheraphy  - Stage 2 - 20 years
x30    = c(1,1,0,0,1,1,1,20)  
theta30 <- exp(t(x30)%*%beta.est) 
p30 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta30)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta30)/(1-m.est^2)) ))

varp30 = p30*(1-p30)*sqrt(t(x30)%*%Varcov%*%(x30))
ICp30  = rbind(round(c( p30 - 1.96*varp30, p30 + 1.96*varp30),3));
round(p30,3);ICp30
#---------------------------------------------------
#31) Surgery +  radiotherapy + Chemotheraphy  - Stage 3 - 20 years
x31    = c(1,0,1,0,1,1,1,20)  
theta31 <- exp(t(x31)%*%beta.est) 
p31 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta31)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta31)/(1-m.est^2)) ))

varp31 = p31*(1-p31)*sqrt(t(x31)%*%Varcov%*%(x31))
ICp31  = rbind(round(c( p31 - 1.96*varp31, p31 + 1.96*varp31),3));
round(p31,3);ICp31
#---------------------------------------------------
#32) Surgery +  radiotherapy + Chemotheraphy  - Stage 4 - 20 years
x32    = c(1,0,0,1,1,1,1,20)  
theta32 <- exp(t(x32)%*%beta.est) 
p32 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta32)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta32)/(1-m.est^2)) ))

varp32 = p32*(1-p32)*sqrt(t(x32)%*%Varcov%*%(x32))
ICp32  = rbind(round(c( p32 - 1.96*varp32, p32 + 1.96*varp32),3));
round(p32,3);ICp32




#===========================================================
# 60 years
#----------------------------------------------------------
#1) No treatment  - Stage 1 - 60 years
x1_60    = c(1,0,0,0,0,0,0,60)  
theta1_60 <- exp(t(x1_60)%*%beta.est) 
p1_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta1_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta1_60)/(1-m.est^2)) ))

varp1_60 = p1_60*(1-p1_60)*sqrt(t(x1_60)%*%Varcov%*%(x1_60))
ICp1_60  = rbind(round(c( p1_60 - 1.96*varp1_60, p1_60 + 1.96*varp1_60),3));
round(p1_60,3);ICp1_60
#---------------------------------------------------
#2) No treatment  - Stage 2 - 60 years
x2_60    = c(1,1,0,0,0,0,0,60) 
theta2_60 <- exp(t(x2_60)%*%beta.est) 
p2_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta2_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta2_60)/(1-m.est^2)) ))
varp2_60 = p2*(1-p2_60)*sqrt(t(x2_60)%*%Varcov%*%(x2_60))
ICp2_60  = rbind(round(c( p2_60 - 1.96*varp2_60, p2_60 + 1.96*varp2_60),3));
round(p2_60,3);ICp2_60
#---------------------------------------------------
#3) No treatment  - Stage 3 - 60 years
x3_60    = c(1,0,1,0,0,0,0,60) 
theta3_60 <- exp(t(x3_60)%*%beta.est) 
p3_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta3_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta3_60)/(1-m.est^2)) ))
varp3_60 = p3_60*(1-p3_60)*sqrt(t(x3_60)%*%Varcov%*%(x3_60))
ICp3_60  = rbind(round(c( p3_60 - 1.96*varp3_60, p3_60 + 1.96*varp3_60),3));
round(p3_60,3);ICp3_60
#---------------------------------------------------
#4) No treatment  - Stage 4 - 60 years
x4_60    = c(1,0,0,1,0,0,0,60) 
theta4_60 <- exp(t(x4_60)%*%beta.est) 
p4_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta4_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta4_60)/(1-m.est^2)) ))
varp4_60 = p4_60*(1-p4_60)*sqrt(t(x4_60)%*%Varcov%*%(x4_60))
ICp4_60  = rbind(round(c( p4_60 - 1.96*varp4_60, p4_60 + 1.96*varp4_60),3));
round(p4_60,3);ICp4_60
#---------------------------------------------------
#5) Surgery  - Stage 1 - 60 years
x5_60    = c(1,0,0,0,1,0,0,60) 
theta5_60 <- exp(t(x5_60)%*%beta.est) 
p5_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta5_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta5_60)/(1-m.est^2)) ))
varp5_60 = p5_60*(1-p5_60)*sqrt(t(x5_60)%*%Varcov%*%(x5_60))
ICp5_60  = rbind(round(c( p5_60 - 1.96*varp5_60, p5_60 + 1.96*varp5_60),3));
round(p5_60,3);ICp5_60
#---------------------------------------------------
#6) Surgery  - Stage 2 - 60 years
x6_60    = c(1,1,0,0,1,0,0,60) 
theta6_60 <- exp(t(x6_60)%*%beta.est) 
p6_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta6_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta6_60)/(1-m.est^2)) ))
varp6_60 = p6_60*(1-p6_60)*sqrt(t(x6_60)%*%Varcov%*%(x6_60))
ICp6_60  = rbind(round(c( p6_60 - 1.96*varp6_60, p6_60 + 1.96*varp6_60),3));
round(p6_60,3);ICp6_60
#---------------------------------------------------
# 7) Surgery  - Stage 3 - 60 years
x7_60    = c(1,0,1,0,1,0,0,60)  
theta7_60 <- exp(t(x7_60)%*%beta.est) 
p7_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta7_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta7_60)/(1-m.est^2)) ))
varp7_60 = p7_60*(1-p7_60)*sqrt(t(x7_60)%*%Varcov%*%(x7_60))
ICp7_60  = rbind(round(c( p7_60 - 1.96*varp7_60, p7_60 + 1.96*varp7_60),3));
round(p7_60,3);ICp7_60
#---------------------------------------------------
#8) Surgery  - Stage 4 - 60 years
x8_60    = c(1,0,0,1,1,0,0,60)  
theta8_60 <- exp(t(x8_60)%*%beta.est) 
p8_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta8_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta8_60)/(1-m.est^2)) ))
varp8_60 = p8_60*(1-p8_60)*sqrt(t(x8_60)%*%Varcov%*%(x8_60))
ICp8_60  = rbind(round(c( p8_60 - 1.96*varp8_60, p8_60 + 1.96*varp8_60),3));
round(p8_60,3);ICp8_60
#---------------------------------------------------
# 9) radiotheraphy  - Stage 1 - 60 years
x9_60    = c(1,0,0,0,0,1,0,60)  
theta9_60 <- exp(t(x9_60)%*%beta.est) 
p9_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta9_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta9_60)/(1-m.est^2)) ))
varp9_60 = p9_60*(1-p9_60)*sqrt(t(x9_60)%*%Varcov%*%(x9_60))
ICp9_60  = rbind(round(c( p9_60 - 1.96*varp9_60, p9_60 + 1.96*varp9_60),3));
round(p9_60,3);ICp9_60
#---------------------------------------------------
#10 radiotheraphy  - Stage 2 - 60 years
x10_60    = c(1,1,0,0,0,1,0,60)  
theta10_60 <- exp(t(x10_60)%*%beta.est) 
p10_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta10_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta10_60)/(1-m.est^2)) ))
varp10_60 = p10_60*(1-p10_60)*sqrt(t(x10_60)%*%Varcov%*%(x10_60))
ICp10_60  = rbind(round(c( p10_60 - 1.96*varp10_60, p10_60 + 1.96*varp10_60),3));
round(p10_60,3);ICp10_60
#---------------------------------------------------
#11) radiotheraphy  - Stage 3 - 60 years
x11_60    = c(1,0,1,0,0,1,0,60)  
theta11_60 <- exp(t(x11_60)%*%beta.est) 
p11_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta11_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta11_60)/(1-m.est^2)) ))
varp11_60 = p11_60*(1-p11_60)*sqrt(t(x11_60)%*%Varcov%*%(x11_60))
ICp11_60  = rbind(round(c( p11_60 - 1.96*varp11_60, p11_60 + 1.96*varp11_60),3));
round(p11_60,3);ICp11_60
#---------------------------------------------------
#12) radiotheraphy  - Stage 4 - 60 years
x12_60    = c(1,0,0,1,0,1,0,60)  
theta12_60 <- exp(t(x12_60)%*%beta.est) 
p12_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta12_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta12_60)/(1-m.est^2)) ))
varp12_60 = p12_60*(1-p12_60)*sqrt(t(x12_60)%*%Varcov%*%(x12_60))
ICp12_60  = rbind(round(c( p12_60 - 1.96*varp12_60, p12_60 + 1.96*varp12_60),3));
round(p12_60,3);ICp12_60
#---------------------------------------------------
#13 Chemotheraphy  - Stage 1 - 60 years
x13_60    = c(1,0,0,0,0,0,1,60)  
theta13_60 <- exp(t(x13_60)%*%beta.est) 
p13_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta13_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta13_60)/(1-m.est^2)) ))
varp13_60 = p13_60*(1-p13_60)*sqrt(t(x13_60)%*%Varcov%*%(x13_60))
ICp13_60  = rbind(round(c( p13_60 - 1.96*varp13_60, p13_60 + 1.96*varp13_60),3));
round(p13_60,3);ICp13_60
#---------------------------------------------------
#14 Chemotheraphy  - Stage 2 - 60 years
x14_60    = c(1,1,0,0,0,0,1,60)  
theta14_60 <- exp(t(x14_60)%*%beta.est) 
p14_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta14_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta14_60)/(1-m.est^2)) ))
varp14_60 = p14_60*(1-p14_60)*sqrt(t(x14_60)%*%Varcov%*%(x14_60))
ICp14_60  = rbind(round(c( p14_60 - 1.96*varp14_60, p14_60 + 1.96*varp14_60),3));
round(p14_60,3);ICp14_60
#---------------------------------------------------
#15) Chemotheraphy  - Stage 3 - 60 years
x15_60    = c(1,0,1,0,0,0,1,60)  
theta15_60 <- exp(t(x15_60)%*%beta.est) 
p15_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta15_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta15_60)/(1-m.est^2)) ))
varp15_60 = p15_60*(1-p15_60)*sqrt(t(x15_60)%*%Varcov%*%(x15_60))
ICp15_60  = rbind(round(c( p15_60 - 1.96*varp15_60, p15_60 + 1.96*varp15_60),3));
round(p15_60,3);ICp15_60
#---------------------------------------------------
#16) Chemotheraphy  - Stage 4 - 60 years
x16_60    = c(1,0,0,1,0,0,1,60)  
theta16_60 <- exp(t(x16_60)%*%beta.est) 
p16_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta16_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta16_60)/(1-m.est^2)) ))
varp16_60 = p16_60*(1-p16_60)*sqrt(t(x16_60)%*%Varcov%*%(x16_60))
ICp16_60  = rbind(round(c( p16_60 - 1.96*varp16_60, p16_60 + 1.96*varp16_60),3));
round(p16_60,3);ICp16_60

#---------------------------------------------------
#17) Surgery + radiotheraphy  - Stage 1 - 60 years
x17_60    = c(1,0,0,0,1,1,0,60)  
theta17_60 <- exp(t(x17_60)%*%beta.est) 
p17_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta17_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta17_60)/(1-m.est^2)) ))
varp17_60 = p17_60*(1-p17_60)*sqrt(t(x17_60)%*%Varcov%*%(x17_60))
ICp17_60  = rbind(round(c( p17_60 - 1.96*varp17_60, p17_60 + 1.96*varp17_60),3));
round(p17_60,3);ICp17_60
#---------------------------------------------------
#18) Surgery + radiotheraphy  - Stage 2 - 60 years
x18_60    = c(1,1,0,0,1,1,0,60)  
theta18_60 <- exp(t(x18_60)%*%beta.est) 
p18_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta18_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta18_60)/(1-m.est^2)) ))
varp18_60 = p18_60*(1-p18_60)*sqrt(t(x18_60)%*%Varcov%*%(x18_60))
ICp18_60  = rbind(round(c( p18_60 - 1.96*varp18_60, p18_60 + 1.96*varp18_60),3));
round(p18_60,3);ICp18_60
#---------------------------------------------------
#19) Surgery + radiotheraphy  - Stage 3 - 60 years
x19_60    = c(1,0,1,0,1,1,0,60)  
theta19_60 <- exp(t(x19_60)%*%beta.est) 
p19_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta19_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta19_60)/(1-m.est^2)) ))
varp19_60 = p19_60*(1-p19_60)*sqrt(t(x19_60)%*%Varcov%*%(x19_60))
ICp19_60  = rbind(round(c( p19_60 - 1.96*varp19_60, p19_60 + 1.96*varp19_60),3));
round(p19_60,3);ICp19_60
#---------------------------------------------------
#20) Surgery + radiotheraphy  - Stage 4 - 60 years
x20_60    = c(1,0,0,1,1,1,0,60)  
theta20_60 <- exp(t(x20_60)%*%beta.est) 
p20_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta20_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta20_60)/(1-m.est^2)) ))
varp20_60 = p20_60*(1-p20_60)*sqrt(t(x20_60)%*%Varcov%*%(x20_60))
ICp20_60  = rbind(round(c( p20_60 - 1.96*varp20_60, p20_60 + 1.96*varp20_60),3));
round(p20_60,3);ICp20_60

#---------------------------------------------------
#21) Surgery + Chemotheraphy  - Stage 1 - 60 years
x21_60    = c(1,0,0,0,1,0,1,60)  
theta21_60 <- exp(t(x21_60)%*%beta.est) 
p21_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta21_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta21_60)/(1-m.est^2)) ))
varp21_60 = p21_60*(1-p21_60)*sqrt(t(x21_60)%*%Varcov%*%(x21_60))
ICp21_60  = rbind(round(c( p21_60 - 1.96*varp21_60, p21_60 + 1.96*varp21_60),3));
round(p21_60,3);ICp21_60
#---------------------------------------------------
#22) Surgery + Chemotheraphy  - Stage 2 - 60 years
x22_60    = c(1,1,0,0,1,0,1,60)  
theta22_60 <- exp(t(x22_60)%*%beta.est) 
p22_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta22_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta22_60)/(1-m.est^2)) ))
varp22_60 = p22_60*(1-p22_60)*sqrt(t(x22_60)%*%Varcov%*%(x22_60))
ICp22_60  = rbind(round(c( p22_60 - 1.96*varp22_60, p22_60 + 1.96*varp22_60),3));
round(p22_60,3);ICp22_60
#---------------------------------------------------
#23) Surgery + Chemotheraphy  - Stage 3 - 60 years
x23_60    = c(1,0,1,0,1,0,1,60)  
theta23_60 <- exp(t(x23_60)%*%beta.est) 
p23_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta23_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta23_60)/(1-m.est^2)) ))
varp23_60 = p23_60*(1-p23_60)*sqrt(t(x23_60)%*%Varcov%*%(x23_60))
ICp23_60  = rbind(round(c( p23_60 - 1.96*varp23_60, p23_60 + 1.96*varp23_60),3));
round(p23_60,3);ICp23_60
#---------------------------------------------------
#24) Surgery + Chemotheraphy  - Stage 4 - 60 years
x24_60    = c(1,0,0,1,1,0,1,60)  
theta24_60 <- exp(t(x24_60)%*%beta.est) 
p24_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta24_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta24_60)/(1-m.est^2)) ))
varp24_60 = p24_60*(1-p24_60)*sqrt(t(x24_60)%*%Varcov%*%(x24_60))
ICp24_60  = rbind(round(c( p24_60 - 1.96*varp24_60, p24_60 + 1.96*varp24_60),3));
round(p24_60,3);ICp24_60
#---------------------------------------------------
#25) radiotherapy + Chemotheraphy  - Stage 1 - 60 years
x25_60    = c(1,0,0,0,0,1,1,60)  
theta25_60 <- exp(t(x25_60)%*%beta.est) 
p25_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta25_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta25_60)/(1-m.est^2)) ))
varp25_60 = p25_60*(1-p25_60)*sqrt(t(x25_60)%*%Varcov%*%(x25_60))
ICp25_60  = rbind(round(c( p25_60 - 1.96*varp25_60, p25_60 + 1.96*varp25_60),3));
round(p25_60,3);ICp25_60
#---------------------------------------------------
#26) radiotherapy + Chemotheraphy  - Stage 2 - 60 years
x26_60    = c(1,1,0,0,0,1,1,60)  
theta26_60 <- exp(t(x26_60)%*%beta.est) 
p26_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta26_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta26_60)/(1-m.est^2)) ))
varp26_60 = p26_60*(1-p26_60)*sqrt(t(x26_60)%*%Varcov%*%(x26_60))
ICp26_60  = rbind(round(c( p26_60 - 1.96*varp26_60, p26_60 + 1.96*varp26_60),3));
round(p26_60,3);ICp26_60
#---------------------------------------------------
#27) radiotherapy + Chemotheraphy  - Stage 3 - 60 years
x27_60    = c(1,0,1,0,0,1,1,60)  
theta27_60 <- exp(t(x27_60)%*%beta.est) 
p27_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta27_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta27_60)/(1-m.est^2)) ))
varp27_60 = p27_60*(1-p27_60)*sqrt(t(x27_60)%*%Varcov%*%(x27_60))
ICp27_60  = rbind(round(c( p27_60 - 1.96*varp27_60, p27_60 + 1.96*varp27_60),3));
round(p27_60,3);ICp27_60
#---------------------------------------------------
#28) radiotherapy + Chemotheraphy  - Stage 4 - 60 years
x28_60    = c(1,0,0,1,0,1,1,60)  
theta28_60 <- exp(t(x28_60)%*%beta.est) 
p28_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta28_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta28_60)/(1-m.est^2)) ))
varp28_60 = p28_60*(1-p28_60)*sqrt(t(x28_60)%*%Varcov%*%(x28_60))
ICp28_60  = rbind(round(c( p28_60 - 1.96*varp28_60, p28_60 + 1.96*varp28_60),3));
round(p28_60,3);ICp28_60
#---------------------------------------------------
#29) Surgery + radiotherapy + Chemotheraphy  - Stage 1 - 60 years
x29_60    = c(1,0,0,0,1,1,1,60)  
theta29_60 <- exp(t(x29_60)%*%beta.est) 
p29_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta29_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta29_60)/(1-m.est^2)) ))
varp29_60 = p29_60*(1-p29_60)*sqrt(t(x29_60)%*%Varcov%*%(x29_60))
ICp29_60  = rbind(round(c( p29_60 - 1.96*varp29_60, p29_60 + 1.96*varp29_60),3));
round(p29_60,3);ICp29_60
#---------------------------------------------------
#30) Surgery +  radiotherapy + Chemotheraphy  - Stage 2 - 60 years
x30_60    = c(1,1,0,0,1,1,1,60)  
theta30_60 <- exp(t(x30_60)%*%beta.est) 
p30_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta30_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta30_60)/(1-m.est^2)) ))
varp30_60 = p30_60*(1-p30_60)*sqrt(t(x30_60)%*%Varcov%*%(x30_60))
ICp30_60  = rbind(round(c( p30_60 - 1.96*varp30_60, p30_60 + 1.96*varp30_60),3));
round(p30_60,3);ICp30_60
#---------------------------------------------------
#31) Surgery +  radiotherapy + Chemotheraphy  - Stage 3 - 60 years
x31_60    = c(1,0,1,0,1,1,1,60)  
theta31_60 <- exp(t(x31_60)%*%beta.est) 
p31_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta31_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta31_60)/(1-m.est^2)) ))
varp31_60 = p31_60*(1-p31_60)*sqrt(t(x31_60)%*%Varcov%*%(x31_60))
ICp31_60  = rbind(round(c( p31_60 - 1.96*varp31_60, p31_60 + 1.96*varp31_60),3));
round(p31_60,3);ICp31_60
#---------------------------------------------------
#32) Surgery +  radiotherapy + Chemotheraphy  - Stage 4 - 60 years
x32_60    = c(1,0,0,1,1,1,1,60)  
theta32_60 <- exp(t(x32_60)%*%beta.est) 
p32_60 =  m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta32_60)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta32_60)/(1-m.est^2)) ))
varp32_60 = p32_60*(1-p32_60)*sqrt(t(x32_60)%*%Varcov%*%(x32_60))
ICp32_60  = rbind(round(c( p32_60 - 1.96*varp32_60, p32_60 + 1.96*varp32_60),3));
round(p32_60,3);ICp32_60


#===================================================
# Breast cancer  
# Table 6: Estimated values of cure rate 
#         for patients 
#===================================================
p_surg_20 = c(p5,p6,p7,p8)
p_radio_20 = c(p9,p10,p11,p12)
p_chemo_20 = c(p13,p14,p15,p16)
p_surg_radio_20 = c(p17,p18,p19,p20)
p_surg_chemo_20 = c(p21,p22,p23,p24)
p_radio_chemo_20 = c(p25,p26,p27,p28)
p_surg_radio_chemo_20 = c(p29,p30,p31,p32)
#summary 60 years
p_surg_60 = c(p5_60,p6_60,p7_60,p8_60)
p_radio_60 = c(p9_60,p10_60,p11_60,p12_60)
p_chemo_60 = c(p13_60,p14_60,p15_60,p16_60)
p_surg_radio_60 = c(p17_60,p18_60,p19_60,p20_60)
p_surg_chemo_60 = c(p21_60,p22_60,p23_60,p24_60)
p_radio_chemo_60 = c(p25_60,p26_60,p27_60,p28_60)
p_surg_radio_chemo_60 = c(p29_60,p30_60,p31_60,p32_60)

pcure_rates = round(rbind( p_surg_20 , p_surg_60 , 
                           p_radio_20 , p_radio_60 ,
                           p_chemo_20, p_chemo_60,
                           p_surg_radio_20, p_surg_radio_60,
                           p_surg_chemo_20, p_surg_chemo_60,
                           p_radio_chemo_20,p_radio_chemo_60,
                           p_surg_radio_chemo_20, p_surg_radio_chemo_60),3)
dimnames(pcure_rates) = list( c("20y_Surgery", "60y_Surgery",
                                "20y_Radiotherapy", "60y_Radiotherapy",
                                "20y_Chemotherapy",  "60y_Chemotherapy",
                                "20y_Surg+Radio", "60y_Surg+Radio  ",
                                "20y_Surg+Chemo", "60y_Surg+Chemo", 
                                "20y_Radio+Chemo",  "60y_Radio+Chemo",
                                "20y_Surg+Radio+Chemo", "60y_Surg+Radio+Chemo"),
                              c("Stage 1"," Stage 2", "Stage 3", "Stage 4"))
pcure_rates
name13 = paste("Table6_p")
write.table(pcure_rates,file=paste(name13,".txt"))
write.csv(pcure_rates,file=paste(name13,".csv"))



#===================================================
# Breast cancer  
# Table 6: Estimated values of IC cure rate part
# 
#===================================================
IC_surg_20 = c(ICp5,ICp6,ICp7,ICp8)
IC_radio_20 = c(ICp9,ICp10,ICp11,ICp12)
IC_chemo_20 = c(ICp13,ICp14,ICp15,ICp16)
IC_surg_radio_20 = c(ICp17,ICp18,ICp19,ICp20)
IC_surg_chemo_20 = c(ICp21,ICp22,ICp23,ICp24)
IC_radio_chemo_20 = c(ICp25,ICp26,ICp27,ICp28)
IC_surg_radio_chemo_20 = c(ICp29,ICp30,ICp31,ICp32)
# summary IC 56 years
IC_surg_60 = c(ICp5_60,ICp6_60,ICp7_60,ICp8_60)
IC_radio_60 = c(ICp9_60,ICp10_60,ICp11_60,ICp12_60)
IC_chemo_60 = c(ICp13_60,ICp14_60,ICp15_60,ICp16_60)
IC_surg_radio_60 = c(ICp17_60,ICp18_60,ICp19_60,ICp20_60)
IC_surg_chemo_60 = c(ICp21_60,ICp22_60,ICp23_60,ICp24_60)
IC_radio_chemo_60 = c(ICp25_60,ICp26_60,ICp27_60,ICp28_60)
IC_surg_radio_chemo_60 = c(ICp29_60,ICp30_60,ICp31_60,ICp32_60)
ICcure_rates = round(rbind(IC_surg_20 ,   IC_surg_60 ,IC_radio_20 ,   IC_radio_60 ,  IC_chemo_20, IC_chemo_60,
                           IC_surg_radio_20,  IC_surg_radio_60,  IC_surg_chemo_20, IC_surg_chemo_60,
                           IC_radio_chemo_20, IC_radio_chemo_60,IC_surg_radio_chemo_20, IC_surg_radio_chemo_60),3)
dimnames(ICcure_rates) = list( c("20y_Surgery", "60y_Surgery",  "20y_Radiotherapy", "60y_Radiotherapy",
                                 "20y_Chemotherapy", "60y_Chemotherapy", "20y_Surg+Radio","60y_Surg+Radio",
                                 "20y_Surg+Chemo",   "60y_Surg+Chemo","20y_Radio+Chemo",  "60y_Radio+Chemo",
                                 "20y_Surg+Radio+Chemo","60y_Surg+Radio+Chemo"),
                               c("St1-LInf","St1-LSup","St2-LInf","St2-LSup","St3-LInf","St3-LSup","St4-LInf","Stage 4-LSup"))
ICcure_rates
name14 = paste("Table6__IC")
write.table(ICcure_rates,file=paste(name14,".txt"))
write.csv(ICcure_rates,file=paste(name14,".csv"))






#===================================================
# Breast cancer  
# Figure 4: Estimated cure rates for the MPIGcr mixture 
#          for patients with 20 and 60 years through treatments
#===================================================

w_in <- 1600/96
h_in <- 850/96

cairo_pdf("Figure4.pdf", width = w_in, height = h_in)

df_long <- pcure_rates %>%
  as.data.frame() %>%
  tibble::rownames_to_column("TreatmentFull") %>%
  rename_with(trimws) %>%
  separate(    TreatmentFull,    into = c("Age", "Treatment"),
               sep = "_",    extra = "merge"  ) %>%
  pivot_longer(    cols = starts_with("Stage"),
                   names_to = "Stage",    values_to = "CureRate"
  )

ggplot(
  df_long,
  aes(     x = Stage,    y = CureRate,    color = Treatment,    linetype = Age,
           group = interaction(Treatment, Age)
  )
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Stages of disease",    y = "Estimated cure rate",
    color = "Treatment",    linetype = "Age"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_blank(),   
    legend.key = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    panel.border = element_rect(
      colour = "black",   fill = NA, linewidth = 1
    )
  )

dev.off()




#===================================================
# Breast cancer  
# Table 7: Estimated values of cure rate and differences
#          percetage points for patients MPIG x PTcr 
#===================================================

psi.MIG     = matrix(psi,nrow=1);colnames(psi.MIG)= c(paste("beta",0:(r-1),sep=""),"alpha","nu","m");round(psi.MIG,4)
Varcov_mig = solve(hessian(llikeobserved,x0=psi.MIG,Data=Data,dist=1,model=1))[1:r,1:r]
#------
Varcov_po = solve(hessian(llikeobserved_PO,x0=psi_PO,data=Data,dist=1))[1:r,1:r]
psi_PO = matrix(EM_PO$estimate[,1],nrow=1);colnames(psi_PO)= c(paste("beta",0:(r-1),sep=""),"alpha","sigma");round(psi_PO,4)
#------
beta.est_MIG = c(psi.MIG[1:r])
beta.est_PO = c(psi_PO[1:r])
m.est   <- psi.MIG[r+3]
#------------------------Stages of disease
# stage = 1 (beta0 intercept)
# stage = 2 (beta1) 
# stage = 3 (betdata_BC)
# stage = 4 (beta3)
#-----------------------Treatments
# Surgery       = beta4
# radiotherapy  = beta5
# chemotherapy  = beta6
# Age           = beta 7 (20,56 or 60 years)
#---------------------------------------------------
# data vector
# x  = c(1,1,0,0,0,0,0,20) # 0 or 1 in beta0 to beta 6 to select 
# stages and treatments and select beta7 in c(20,56,60)
#---------------------------------------------------
#1)  Stage 1 - Radio + chemotherapy  -  56 years
#---------------------------------------------------
x1    = c(1,0,0,0, 0,1,1,56)  
theta1_mig <- exp(t(x1)%*%beta.est_MIG) 
p1_mig = m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta1_mig)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta1_mig)/(1-m.est^2)) ))

varp1_mig = p1_mig*(1-p1_mig)*sqrt(t(x1)%*%Varcov_mig%*%(x1))
ICp1_mig  = rbind(round(c( p1_mig - 1.96*varp1_mig, p1_mig + 1.96*varp1_mig),3));
round(p1_mig,3);ICp1_mig

theta1_po <- exp(t(x1)%*%beta.est_PO) 
p1_po = 1/exp(theta1_po)
varp1_po = p1_po*(1-p1_po)*sqrt(t(x1)%*%Varcov_po%*%(x1))
ICp1_po  = rbind(round(c( p1_po - 1.96*varp1_po, p1_po + 1.96*varp1_po),3));
round(p1_po,3);ICp1_po


#---------------------------------------------------
#2)  Stage 2 - Radio + chemotherapy  -  56 years
#---------------------------------------------------
x2    = c(1,1,0,0, 0,1,1,56)  
thetdata_BC_mig <- exp(t(x2)%*%beta.est_MIG) 
p2_mig = m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*thetdata_BC_mig)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*thetdata_BC_mig)/(1-m.est^2)) ))

varp2_mig = p2_mig*(1-p2_mig)*sqrt(t(x2)%*%Varcov_mig%*%(x2))
ICp2_mig  = rbind(round(c( p2_mig - 1.96*varp2_mig, p2_mig + 1.96*varp2_mig),3));
round(p2_mig,3);ICp2_mig

thetdata_BC_po <- exp(t(x2)%*%beta.est_PO) 
p2_po = 1/exp(thetdata_BC_po)
varp2_po = p2_po*(1-p2_po)*sqrt(t(x2)%*%Varcov_po%*%(x2))
ICp2_po  = rbind(round(c( p2_po - 1.96*varp2_po, p2_po + 1.96*varp2_po),3));
round(p2_po,3);ICp2_po


#---------------------------------------------------
#3)  Stage 3 - Radio + chemotherapy  -  56 years
#---------------------------------------------------
x3    = c(1,0,1,0, 0,1,1,56)  
theta3_mig <- exp(t(x3)%*%beta.est_MIG) 
p3_mig = m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta3_mig)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta3_mig)/(1-m.est^2)) ))

varp3_mig = p3_mig*(1-p3_mig)*sqrt(t(x3)%*%Varcov_mig%*%(x3))
ICp3_mig  = rbind(round(c( p3_mig - 1.96*varp3_mig, p3_mig + 1.96*varp3_mig),3));
round(p3_mig,3);ICp3_mig

theta3_po <- exp(t(x3)%*%beta.est_PO) 
p3_po = 1/exp(theta3_po)
varp3_po = p3_po*(1-p3_po)*sqrt(t(x3)%*%Varcov_po%*%(x3))
ICp3_po  = rbind(round(c( p3_po - 1.96*varp3_po, p3_po + 1.96*varp3_po),3));
round(p3_po,3);ICp3_po

#---------------------------------------------------
#3)  Stage 4 - Radio + chemotherapy  -  56 years
#---------------------------------------------------
x4    = c(1,0,0,1, 0,1,1,56)  
theta4_mig <- exp(t(x4)%*%beta.est_MIG) 
p4_mig = m.est*exp((1-m.est)*(1-sqrt(1+(2*m.est*theta4_mig)/(1-m.est)))) + 
  (1-m.est)*exp((((m.est+1)*(1-m.est^2))/(m.est^2))*(1-sqrt(1+(2*m.est^2*theta4_mig)/(1-m.est^2)) ))

varp4_mig = p4_mig*(1-p4_mig)*sqrt(t(x4)%*%Varcov_mig%*%(x4))
ICp4_mig  = rbind(round(c( p4_mig - 1.96*varp4_mig, p4_mig + 1.96*varp4_mig),3));
round(p4_mig,3);ICp4_mig

theta4_po <- exp(t(x4)%*%beta.est_PO) 
p4_po = 1/exp(theta4_po)
varp4_po = p4_po*(1-p4_po)*sqrt(t(x4)%*%Varcov_po%*%(x4))
ICp4_po  = rbind(round(c( p4_po - 1.96*varp4_po, p4_po + 1.96*varp4_po),3));
round(p4_po,3);ICp4_po


pcure_PO = c(p1_po,p2_po,p3_po,p4_po);pcure_PO
pcure_MIG = c(p1_mig,p2_mig,p3_mig ,p4_mig);pcure_MIG 

IC_pcure_PO = c(ICp1_po,ICp2_po,ICp3_po,ICp4_po);IC_pcure_PO
IC_pcure_MIG = c(ICp1_mig,ICp2_mig,ICp3_mig, ICp4_mig);IC_pcure_MIG

Model = c("PO","MIG")
Stages = c(1,2,3,4)



#=======================================
# Table 7: Estimated values of cure rate and differences
#          percetage points for patients MPIG x PTcr 
#=========================================
pcure_PO  <- c(p1_po, p2_po, p3_po, p4_po)
pcure_MIG <- c(p1_mig, p2_mig, p3_mig, p4_mig)
IC_pcure_PO  <- c(ICp1_po, ICp2_po, ICp3_po, ICp4_po)
IC_pcure_MIG <- c(ICp1_mig, ICp2_mig, ICp3_mig, ICp4_mig)
Stages <- c(1,2,3,4)
# IC (lwr,upr) by stage
IC_PO  <- matrix(IC_pcure_PO,  ncol = 2, byrow = TRUE)
IC_MIG <- matrix(IC_pcure_MIG, ncol = 2, byrow = TRUE)
tab_out <- data.frame(
  Stage   = Stages,
  PO_hat  = pcure_PO,  PO_lwr = IC_PO[,1],  PO_upr = IC_PO[,2],
  MIG_hat = pcure_MIG, MIG_lwr = IC_MIG[,1], MIG_upr = IC_MIG[,2]
) %>%
  mutate(
    `PO (95% CI)`  = sprintf("%.3f [%.3f, %.3f]", PO_hat,  PO_lwr,  PO_upr),
    `MIG (95% CI)` = sprintf("%.3f [%.3f, %.3f]", MIG_hat, MIG_lwr, MIG_upr),
    ` (pp)`       = sprintf("%+.1f", 100*(MIG_hat - PO_hat))
  ) %>%
  select(Stage, `PO (95% CI)`, `MIG (95% CI)`, ` (pp)`)


write.table(tab_out,
            file = "Table7.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.csv(tab_out,
          file = "Table7.csv",
          row.names = FALSE)







#===================================================
# Breast cancer  
# Figure 5: Estimated cure rates with 95% CI and percentage-point difference for 
#            the PTcr and MPIGcr mixture for patient aged 56 years through treatments. 
#===================================================
pcure_PO  <- c(p1_po, p2_po, p3_po, p4_po)
pcure_MIG <- c(p1_mig, p2_mig, p3_mig, p4_mig)
IC_pcure_PO  <- c(ICp1_po, ICp2_po, ICp3_po, ICp4_po)
IC_pcure_MIG <- c(ICp1_mig, ICp2_mig, ICp3_mig, ICp4_mig)
Stages <- c(1,2,3,4)
IC_PO  <- matrix(IC_pcure_PO,  ncol = 2, byrow = TRUE, dimnames = list(NULL, c("lwr","upr")))
IC_MIG <- matrix(IC_pcure_MIG, ncol = 2, byrow = TRUE, dimnames = list(NULL, c("lwr","upr")))
#----------------------------------
tab <- bind_rows(
  data.frame(Stage = Stages, Model = "PO",
             perc = pcure_PO,  lwr = IC_PO[,1],  upr = IC_PO[,2]),
  data.frame(Stage = Stages, Model = "MIG",
             perc = pcure_MIG, lwr = IC_MIG[,1], upr = IC_MIG[,2])
) %>%
  mutate(
    Stage = factor(Stage, levels = c(4,3,2,1)),
    Model = factor(Model, levels = c("PO","MIG"))
  )
#-------------------------------
wide <- tab %>%
  select(Stage, Model, perc) %>%
  pivot_wider(names_from = Model, values_from = perc)
#-------------------------------
wide2 <- wide %>%
  mutate(
    diff_abs_pct = 100 * (MIG - PO),
    lab_abs      = sprintf("%+.1f%%", diff_abs_pct),
    y_mid        = (PO + MIG)/2
  )


w_in <- 1018/96
h_in <- 300/96
cairo_pdf("Figure5.pdf", width = w_in, height=h_in)

ggplot() +
  geom_segment(data = wide2,
               aes(y = Stage, x = PO, xend = MIG, yend = Stage),
               linewidth = 0.9, alpha = 0.5, colour = "grey50") +
  geom_errorbarh(data = tab,
                 aes(y = Stage, xmin = lwr, xmax = upr, colour = Model),
                 height = 0.18, alpha = 0.35, linewidth = 2) +
  geom_point(data = tab,
             aes(y = Stage, x = perc, colour = Model, shape = Model),
             size = 3.2) +
  geom_text(data = wide2,
            aes(y = Stage, x = y_mid, label = lab_abs),
            vjust = -1.0, size = 3.5) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  
  scale_colour_discrete(labels = c("PTcr", "MPIGcr")) +
  scale_shape_discrete(labels = c("PTcr", "MPIGcr")) +
  labs(x = "Estimated cure rate", y = "Stage", colour = "", shape = "") +
  theme_minimal()+ theme(legend.position = "bottom")
dev.off()


