
library(MASS)
library(tidyverse)
library(ISLR)
library(tree)
library(ggcorrplot)
library(visdat)
library(mice)
library(naniar)
library(UpSetR)
library(VIM)
library(GGally)
library(cowplot)
library(car)
library(gridExtra)
library(heplots)
library(mclust)
library(class)
library(caret)
library(ROCR)
select=dplyr::select
heart<-read.csv("C:/Users/fccat/Documents/Projects/Heart_Failure_Analysis//heart.csv",stringsAsFactors= T)
heart$FastingBS=as.factor(heart$FastingBS)
heart$HeartDisease=as.factor(heart$HeartDisease)
anyNA(heart)


#sostituz zeri con NA#####
heart$Cholesterol[heart$Cholesterol==0]=NA
#strategia passiva su RestingBP
heart=heart%>%filter(!is.na(RestingBP))

#divisione test training e validation####
set.seed(1)
test.index=sample(c(0:nrow(heart)),size = 0.2*nrow(heart),F)
test=heart[test.index,]
train_valid=anti_join(heart,test)
validation=train_valid%>%slice_sample(prop=0.2/0.8)
training=anti_join(train_valid,validation)

#verifica bilanciamento target
table(training$HeartDisease)/nrow(training)
table(validation$HeartDisease)/nrow(validation)
table(test$HeartDisease)/nrow(test)
#bilanciata (set seed perchè lo sia sempre)

###metriche valutazione (funzione)#####
metriche=function(tab){
  accuracy=sum(diag(tab))/sum(tab)
  sensitivity=tab[2,2]/sum(tab[,2])
  #possiamo aggiungerne
  return(c("Accuracy"=accuracy,"Sensitivity"=sensitivity))
}


md.pattern(training)

g1=ggplot(training,aes(x=Age,y=after_stat(density)))+
  geom_histogram(aes(fill=is.na(Cholesterol)),col="black",bins=15)+
  theme_minimal()
g2=ggplot(training,aes(x=RestingBP,y=after_stat(density)))+
  geom_histogram(aes(fill=is.na(Cholesterol)),col="black",bins=18)+
  theme_minimal()
g3=ggplot(training,aes(x=MaxHR,y=after_stat(density)))+
  geom_histogram(aes(fill=is.na(Cholesterol)),col="black",bins=15)+
  theme_minimal()
g4=ggplot(training,aes(x=Oldpeak,y=after_stat(density)))+
  geom_histogram(aes(fill=is.na(Cholesterol)),col="black")+
  theme_minimal()
plot_grid(g1,g2,g3,g4)

#correlazioni
ggcorrplot(cor(train_no_NA),show.diag = F,lab = T,type = "upper")
#nessuna correlazione particolarmente alta

#outlier####
boxplot(training%>%select(where(is.numeric)))

#standardizzazione ####

means=colMeans(training%>%select(where(is.numeric)),na.rm=T) #memorizzo medie
variances=apply(training%>%select(where(is.numeric)),2,function(x){var(x,na.rm=T)})  #memorizzo varianze
training_non_stand=training
validation_non_stand=validation
test_non_stand=test
heart_non_stand=heart
training[,c("Age","RestingBP","Cholesterol","MaxHR","Oldpeak")]=scale(training%>%select(where(is.numeric)),center=means,scale = sqrt(variances))
validation[,c("Age","RestingBP","Cholesterol","MaxHR","Oldpeak")]=scale(validation%>%select(where(is.numeric)),center=means,scale = sqrt(variances))
test[,c("Age","RestingBP","Cholesterol","MaxHR","Oldpeak")]=scale(test%>%select(where(is.numeric)),center=means,scale=sqrt(variances))
heart[,c("Age","RestingBP","Cholesterol","MaxHR","Oldpeak")]=scale(heart%>%select(where(is.numeric)),center=means,scale=sqrt(variances))
training_and_valid=rbind(training,validation)

######imputation####
matrix=matrix(1,nrow=12,ncol=12)
matrix[,12]=0
tempDatapmm <- mice(training_and_valid,m=1,maxit=50,meth='pmm',seed=1,predictorMatrix = matrix,printFlag = F)
tempDatamean=mice(training_and_valid,m=1,maxit=50,meth='mean',seed=1,predictorMatrix = matrix,printFlag = F)
tempDatanorm=mice(training_and_valid,m=1,maxit=50,meth='norm',seed=1,predictorMatrix = matrix,printFlag = F)
tempDatanorm.nob=mice(training_and_valid,m=1,maxit=50,meth='norm.nob',seed=1,predictorMatrix = matrix,printFlag = F)
tempDatanormpred=mice(training_and_valid,m=1,maxit=50,meth='norm.predict',seed=1,predictorMatrix = matrix,printFlag = F)
tempDataLasso=mice(training_and_valid,m=1,maxit=50,meth='lasso.norm',seed=1,predictorMatrix = matrix,printFlag = F)
tempDatacart=mice(training_and_valid,m=1,maxit=50,meth='cart',seed=1,predictorMatrix = matrix,printFlag = F)

mice_imputed <- data.frame(
  original = training_and_valid$Cholesterol,
  imputed_pmm = complete(tempDatapmm)$Cholesterol,
  imputed_cart = complete(tempDatacart)$Cholesterol,
  imputed_lasso = complete(tempDataLasso)$Cholesterol,
  imputed_mean = complete(tempDatamean)$Cholesterol,
  imputed_norm = complete(tempDatanorm)$Cholesterol,
  imputed_normpred = complete(tempDatanormpred)$Cholesterol
)

bins_FD=round(apply(mice_imputed,2,function(x) (max(x,na.rm = T)-min(x,na.rm = T))/(2*IQR(x,na.rm = T))*length(x)^(1/3)))
h1_mice <- ggplot(mice_imputed, aes(x = original,y=after_stat(density))) +
  geom_histogram(fill = "#ad1538", color = "#000000", position = "identity",bins=bins_FD[1]) +
  ggtitle("Original distribution") +
  theme_classic()
h2_mice <- ggplot(mice_imputed, aes(x = imputed_pmm,y=after_stat(density))) +
  geom_histogram(fill = "#15ad4f", color = "#000000", position = "identity",bins=bins_FD[2]) +
  ggtitle("PMM-imputed distribution") +
  theme_classic()
h3_mice <- ggplot(mice_imputed, aes(x = imputed_cart,y=after_stat(density))) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity",bins=bins_FD[3]) +
  ggtitle("CART-imputed distribution") +
  theme_classic()
h4_mice <- ggplot(mice_imputed, aes(x = imputed_lasso,y=after_stat(density))) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity",bins=bins_FD[4]) +
  ggtitle("LR-imputed distribution") +
  theme_classic()
h5_mice <- ggplot(mice_imputed, aes(x = imputed_mean,y=after_stat(density))) +
  geom_histogram(fill = "magenta", color = "#000000", position = "identity",bins=bins_FD[5]) +
  ggtitle("LR-imputed distribution") +
  theme_classic()
h6_mice <- ggplot(mice_imputed, aes(x = imputed_norm,y=after_stat(density))) +
  geom_histogram(fill = "gold", color = "#000000", position = "identity",bins=bins_FD[6]) +
  ggtitle("LR-imputed distribution") +
  theme_classic()
h7_mice <- ggplot(mice_imputed, aes(x = imputed_normpred,y=after_stat(density))) +
  geom_histogram(fill = "purple", color = "#000000", position = "identity",bins=bins_FD[7]) +
  ggtitle("LR-imputed distribution") +
  theme_classic()
plot_grid(h1_mice, h2_mice, h3_mice, h4_mice, h5_mice,h6_mice,h7_mice)

densityplot(tempDatapmm,main="PMM",thicker = 0.7)
densityplot(tempDataLasso,main="Lasso",thicker = 0.7)
densityplot(tempDatacart,main="Cart",thicker = 0.7)

#####dataset con imputazione####
training_original=training
validation_original=validation
heart_imputed=mice(heart,predictorMatrix = matrix,meth="cart",m=1,printFlag = F,maxit = 50,ignore=c(1:nrow(heart))%in%test.index)
test$Cholesterol=complete(heart_imputed)[test.index,"Cholesterol"]
training=inner_join(training_original,complete(heart_imputed),join_by(Age,Sex,ChestPainType,RestingBP,FastingBS,RestingECG,MaxHR,Oldpeak,ExerciseAngina,HeartDisease,ST_Slope))%>%
  mutate(Cholesterol=Cholesterol.y)%>%
  select(!contains("."))
validation=inner_join(validation_original,complete(heart_imputed),join_by(Age,Sex,ChestPainType,RestingBP,FastingBS,RestingECG,MaxHR,Oldpeak,ExerciseAngina,HeartDisease,ST_Slope))%>%
  mutate(Cholesterol=Cholesterol.y)%>%
  select(!contains("."))


#verifica assunzioni####
covEllipses(training%>%select(where(is.numeric)), 
            factor(training$HeartDisease), 
            fill = TRUE, 
            pooled = FALSE, 
            col = c("blue", "red"),
            variables=c("Age","RestingBP","Cholesterol","MaxHR","Oldpeak"),
            fill.alpha = 0.05)

#boxplot condizionati
box1=ggplot(training,aes(y=Age,fill=HeartDisease))+
  geom_boxplot(show.legend = F)+
  theme_minimal()
box2=ggplot(training,aes(y=RestingBP,fill=HeartDisease))+
  geom_boxplot(show.legend = F)+
  theme_minimal()
box3=ggplot(training,aes(y=Cholesterol,fill=HeartDisease))+
  geom_boxplot(show.legend = F)+
  theme_minimal()
box4=ggplot(training,aes(y=MaxHR,fill=HeartDisease))+
  geom_boxplot(show.legend = F)+
  theme_minimal()
box5=ggplot(training,aes(y=Oldpeak,fill=HeartDisease))+
  geom_boxplot(show.legend = F)+
  theme_minimal()
box_legend=get_legend(ggplot(training,aes(y=Age,fill=HeartDisease))+
                        geom_boxplot(show.legend = T)+
                        theme_minimal())
plot_grid(box1,box2,box3,box4,box5,box_legend,ncol=3)

par(mfrow=c(2,4))
qqnorm(training$Age[training$HeartDisease==1],main = "Age| Y=1");qqline(training$Age[training$HeartDisease==1],col=2,lwd=2)
qqnorm(training$Age[training$HeartDisease==0],main = "Age| Y=0");qqline(training$Age[training$HeartDisease==0],col=2,lwd=2)
qqnorm(training$RestingBP[training$HeartDisease==1],main = "RestingBP| Y=1");qqline(training$RestingBP[training$HeartDisease==1],col=2,lwd=2)
qqnorm(training$RestingBP[training$HeartDisease==0],main = "RestingBP| Y=0");qqline(training$RestingBP[training$HeartDisease==0],col=2,lwd=2)
qqnorm(training$MaxHR[training$HeartDisease==1],main = "MaxHR| Y=1");qqline(training$MaxHR[training$HeartDisease==1],col=2,lwd=2)
qqnorm(training$MaxHR[training$HeartDisease==0],main = "MaxHR| Y=0");qqline(training$MaxHR[training$HeartDisease==0],col=2,lwd=2)
qqnorm(training$Oldpeak[training$HeartDisease==1],main = "Oldpeak| Y=1");qqline(training$Oldpeak[training$HeartDisease==1],col=2,lwd=2)
qqnorm(training$Oldpeak[training$HeartDisease==0],main = "Oldpeak| Y=0");qqline(training$Oldpeak[training$HeartDisease==0],col=2,lwd=2)
par(mfrow=c(1,1))


apply(training%>%filter(HeartDisease==1)%>%select(where(is.numeric)),2,function(x) ks.test(x,y="dnorm")$p.value)
apply(training%>%filter(HeartDisease==1)%>%select(where(is.numeric)),2,function(x) shapiro.test(x)$p.value)
apply(training%>%filter(HeartDisease==0)%>%select(where(is.numeric)),2,function(x) ks.test(x,y="dnorm")$p.value)
apply(training%>%filter(HeartDisease==0)%>%select(where(is.numeric)),2,function(x) shapiro.test(x)$p.value)
ggplot(training,aes(x=Age,col=HeartDisease))+
  geom_density()

age_dens=ggplot(training,aes(x=Age,col=HeartDisease))+
  geom_density(linewidth=1,show.legend = F)+
  theme_minimal()
resting_dens=ggplot(training,aes(x=RestingBP,col=HeartDisease))+
  geom_density(linewidth=1,show.legend = F)+
  theme_minimal()
maxhr_dens=ggplot(training,aes(x=MaxHR,col=HeartDisease))+
  geom_density(linewidth=1,show.legend = F)+
  theme_minimal()
oldpeak_dens=ggplot(training,aes(x=Oldpeak,col=HeartDisease))+
  geom_density(linewidth=1,show.legend = F)+
  theme_minimal()
choles_dens=ggplot(training,aes(x=Cholesterol,col=HeartDisease))+
  geom_density(linewidth=1,show.legend = F)+
  theme_minimal()
legend=get_legend(ggplot(training,aes(x=Age,col=HeartDisease))+
                    geom_density(linewidth=1)+
                    theme_minimal())
plot_grid(age_dens,resting_dens,maxhr_dens,oldpeak_dens,choles_dens,legend,ncol=3)


#MDA#####
mod_mda=MclustDA(training[,-12]%>%select(where(is.numeric)),class = training$HeartDisease,G=c(1:9))
pred_mda=predict(mod_mda,validation[,-12]%>%select(where(is.numeric)))
confmat_mda=table(pred_mda$classification,validation$HeartDisease)
metriche(confmat_mda)

####classification trees####
#trattano direttamente gli NA
#sul trainning
set.seed(8)
tree.Heart <- tree(HeartDisease ~.,data=training_original)
summary(tree.Heart) #MER: 10,78%

plot(tree.Heart)
text(tree.Heart,pretty = 0) #ST_slope più significativa

tree.pred <- predict(tree.Heart , validation_original , type = "class")
confmat_tree<-table(tree.pred , validation_original$HeartDisease)
metriche(confmat_tree)
#Ci serve la sensitivity poichè problema medico(è più importante classificare correttamente i malati 
#piuttosto che i non malati )
 
#pruning###
cv.Heart <- cv.tree(tree.Heart , FUN = prune.misclass)
par(mfrow = c(1, 2))
plot(cv.Heart$size , cv.Heart$dev, type = "b")
plot(cv.Heart$k, cv.Heart$dev, type = "b")

par(mfrow = c(1, 1))
best_size<-cv.Heart$size[which.min(cv.Heart$dev)]
prune.Heart <- prune.misclass(tree.Heart , best = best_size)

plot(prune.Heart)
text(prune.Heart , pretty = 0)

# validation
tree.pred_pruning <- predict(prune.Heart , validation_original , type = "class")
confmat_pruning<-table(tree.pred_pruning , validation_original$HeartDisease)
metriche(confmat_pruning)

#KNN####
train_dummy=training%>%mutate(Sex=ifelse(Sex=="M",1,0),ChestPainTypeATA=ifelse(ChestPainType=="ATA",1,0),ChestPainTypeNAP=ifelse(ChestPainType=="NAP",1,0),
                           ChestPainTypeTA=ifelse(ChestPainType=="TA",1,0),RestingECGNormal=ifelse(RestingECG=="Normal",1,0),
                           RestingECGST=ifelse(RestingECG=="ST",1,0),ExerciseAngina=ifelse(ExerciseAngina=="Y",1,0),
                           ST_SlopeFlat=ifelse(ST_Slope=="Flat",1,0),ST_SlopeUp=ifelse(ST_Slope=="Up",1,0))
train_dummy=train_dummy[,-c(3,6,10)]
valid_dummy=validation%>%mutate(Sex=ifelse(Sex=="M",1,0),ChestPainTypeATA=ifelse(ChestPainType=="ATA",1,0),ChestPainTypeNAP=ifelse(ChestPainType=="NAP",1,0),
                                ChestPainTypeTA=ifelse(ChestPainType=="TA",1,0),RestingECGNormal=ifelse(RestingECG=="Normal",1,0),
                                RestingECGST=ifelse(RestingECG=="ST",1,0),ExerciseAngina=ifelse(ExerciseAngina=="Y",1,0),
                                ST_SlopeFlat=ifelse(ST_Slope=="Flat",1,0),ST_SlopeUp=ifelse(ST_Slope=="Up",1,0))
valid_dummy=valid_dummy[,-c(3,6,10)]

ks=seq(1,21,by=2)
knn_cv=apply(as.matrix(ks),1,function(x){mod_knn=knn(train_dummy[,-8],k=x,test=valid_dummy[,-8],cl=train_dummy$HeartDisease)
confmat_knn_cv=table(mod_knn,valid_dummy$HeartDisease)
c(x,confmat_knn_cv[2,2]/sum(confmat_knn_cv[,2]),sum(diag(confmat_knn_cv))/nrow(validation))})
k_cv=knn_cv[1,which.max(knn_cv[2,])]
mod_knn_cv=knn(train_dummy[,-8],k=k_cv,test=valid_dummy[,-8],cl=training$HeartDisease)
confmat_knn=table(mod_knn_cv,validation$HeartDisease)
metriche(confmat_knn)


#Reg logistica####

logistic_mod<-glm(HeartDisease~. ,data = training, family="binomial")
summary(logistic_mod)

#per step
step_logistic_mod<-step(logistic_mod, direction = "both", trace = F)
#glm(formula = HeartDisease ~ Sex + ChestPainType + Cholesterol + FastingBS + MaxHR + ExerciseAngina + Oldpeak + ST_Slope, 
#family = "binomial", data = training)
summary(step_logistic_mod)
anova(step_logistic_mod,test = "Chisq")

#outliers and influence points####
outlierTest(step_logistic_mod)
infplot=influencePlot(step_logistic_mod)
obs_inf=as.numeric(rownames(infplot))
#provo a togliere i punti influenti e vedere se AIC migliora
#questi punti sono quelli individuati da influenceplot
logistic_mod_no_inf<-glm(HeartDisease~. ,data = training[-obs_inf,], family="binomial")
summary(logistic_mod_no_inf)
#noto un miglioramenti in termini di AIC rispetto al logistic_mod

#per step
step_logistic_mod_no_inf<-step(logistic_mod_no_inf, direction = "both", trace = F)
summary(step_logistic_mod_no_inf)
#ulteriore miglioramento in termini di AIC rispetto al step_logistic_mod
anova(step_logistic_mod_no_inf,test = "Chisq")


#step_logistic_mod_no_inf
step_logistic_mod_no_inf_pred <- round(predict(step_logistic_mod_no_inf, validation, type = "response"),0)
confmat_step_logistic_mod_no_inf<-table(step_logistic_mod_no_inf_pred , validation$HeartDisease)
metriche(confmat_step_logistic_mod_no_inf)
#step_logistic_mod
step_logistic_mod_pred <- round(predict(step_logistic_mod, validation, type = "response"),0)
confmat_step_logistic_mod<-table(step_logistic_mod_pred , validation$HeartDisease)
metriche(confmat_step_logistic_mod)
#logistic_mod
logistic_mod_pred <- round(predict(logistic_mod, validation, type = "response"),0)
confmat_logistic_mod<-table(logistic_mod_pred , validation$HeartDisease)
metriche(confmat_logistic_mod)
#logistic_mod_no_inf
logistic_mod_no_inf_pred <- round(predict(logistic_mod_no_inf, validation, type = "response"),0)
confmat_logistic_mod_no_inf<-table(logistic_mod_no_inf_pred , validation$HeartDisease)
metriche(confmat_logistic_mod_no_inf)
#i migliori sono il logistic_mod_no_inf e step_logistic_mod_no_inf sia in accuracy che in sesitivity

#diagnostica assunzioni del step_logistic_mod_no_inf ####

probabilities <- predict(step_logistic_mod_no_inf, type = "response")
predictors <- c("Sex","ChestPainType","Cholesterol","FastingBS","MaxHR","ExerciseAngina","Oldpeak","ST_Slope")

logit_training<-training[-obs_inf,]%>%
  select(all_of(predictors))

logit_training <- logit_training %>%
  mutate(logit = log(probabilities/(1-probabilities)))

text = paste("\n Linearità logit-covariate")

plot_grid(ggplot() + 
            annotate("text", x = 4, y = 25, size=8, label = text) + 
            theme_void(),
          ggplot(data=logit_training, aes(Cholesterol,logit))+
            geom_point(alpha=.5,col="orange")+
            geom_smooth(method="loess",col="red")+
            theme_minimal(),
          ggplot(data=logit_training, aes(MaxHR,logit))+
            geom_point(alpha=.5,col="cornflowerblue")+
            geom_smooth(method="loess",col="blue")+
            theme_minimal(),
          ggplot(data=logit_training, aes(Oldpeak,logit))+
            geom_point(alpha=.5,col="green")+
            geom_smooth(method="loess",col="darkgreen")+
            theme_minimal()
)
#Selezione del modello#####
confmats=list(mda=confmat_mda,tree=confmat_tree,pruning=confmat_pruning
              ,knn=confmat_knn,reglog=confmat_logistic_mod,reglognf=confmat_logistic_mod_no_inf,
              regstep=confmat_step_logistic_mod,regstepnf=confmat_step_logistic_mod_no_inf)
sapply(confmats,metriche)
#f1 score
f1_score=function(tab){
  Recall=tab[2,2]/sum(tab[,2])
  Precision=tab[2,2]/sum(tab[2,])
  f1=2*Recall*Precision/(Recall+Precision)
  return(f1)
}
sapply(confmats,f1_score)
#curve roc
#possiamo già escludere l'albero potato e non e la mda

prediction_reglog=prediction(step_logistic_mod_pred,validation$HeartDisease)
performance(prediction_reglog,measure="auc")@y.values
knn.pred=as.numeric(as.vector(mod_knn_cv))
prediction_knn=prediction(knn.pred,validation$HeartDisease)
performance(prediction_knn,measure="auc")@y.values
#par(mfrow=c(2,1))
plot(performance(prediction_knn,"tpr","fpr"),main="ROC knn")
plot(performance(prediction_reglog,"tpr","fpr"),main="ROC logistic regression")
#par(mfrow=c(1,1))

#####Testing######
#standardizzazione 
big_training=rbind(training_non_stand,validation_non_stand)
medie=colMeans(big_training%>%select(where(is.numeric)),na.rm=T) 
varianze=apply(big_training%>%select(where(is.numeric)),2,function(x){var(x,na.rm=T)})
big_training[,c("Age","RestingBP","Cholesterol","MaxHR","Oldpeak")]=scale(big_training[,c("Age","RestingBP","Cholesterol","MaxHR","Oldpeak")],center=medie,scale=sqrt(varianze))
test_non_stand[,c("Age","RestingBP","Cholesterol","MaxHR","Oldpeak")]=scale(test_non_stand[,c("Age","RestingBP","Cholesterol","MaxHR","Oldpeak")],center=medie,scale=sqrt(varianze))
final.test=test_non_stand
heart_non_stand[,c("Age","RestingBP","Cholesterol","MaxHR","Oldpeak")]=scale(heart_non_stand[,c("Age","RestingBP","Cholesterol","MaxHR","Oldpeak")],center=medie,scale=sqrt(varianze))
final.heart=heart_non_stand

#imputazione
heart_imputato=mice(final.heart,predictorMatrix = matrix,meth="cart",m=1,printFlag = F,maxit = 50,ignore=c(1:nrow(heart))%in%test.index)
final.test$Cholesterol=complete(heart_imputato)[test.index,"Cholesterol"]
big_training=inner_join(big_training,complete(heart_imputato),join_by(Age,Sex,ChestPainType,RestingBP,FastingBS,RestingECG,MaxHR,Oldpeak,ExerciseAngina,HeartDisease,ST_Slope))%>%
  mutate(Cholesterol=Cholesterol.y)%>%
  select(!contains("."))

#portiamo avanti il glm e knn 

reg_log<-glm(HeartDisease~. ,data = big_training, family="binomial")
summary(reg_log)

step_reg_log<-step(reg_log, direction = "both", trace = F)
summary(step_reg_log)
outlierTest(step_reg_log)
confmat_test_stepreglog=table(round(predict(step_reg_log,final.test,type="response")),final.test$HeartDisease)
metriche(confmat_test_stepreglog)
influence_plot=influencePlot(step_reg_log)
obs_influenti=as.numeric(rownames(influence_plot))
step_reg_log_no_inf<-glm(HeartDisease~. ,data = training[-obs_influenti,], family="binomial")
confmat_test_stepreglog_noinf=table(round(predict(step_reg_log_no_inf,final.test,type="response")),final.test$HeartDisease)
metriche(confmat_test_stepreglog_noinf)
#Knn
big_training_dummy=big_training%>%mutate(Sex=ifelse(Sex=="M",1,0),ChestPainTypeATA=ifelse(ChestPainType=="ATA",1,0),ChestPainTypeNAP=ifelse(ChestPainType=="NAP",1,0),
       ChestPainTypeTA=ifelse(ChestPainType=="TA",1,0),RestingECGNormal=ifelse(RestingECG=="Normal",1,0),
       RestingECGST=ifelse(RestingECG=="ST",1,0),ExerciseAngina=ifelse(ExerciseAngina=="Y",1,0),
       ST_SlopeFlat=ifelse(ST_Slope=="Flat",1,0),ST_SlopeUp=ifelse(ST_Slope=="Up",1,0))
big_training_dummy=big_training_dummy[,-c(3,6,10)]

final_test_dummy=final.test%>%mutate(Sex=ifelse(Sex=="M",1,0),ChestPainTypeATA=ifelse(ChestPainType=="ATA",1,0),ChestPainTypeNAP=ifelse(ChestPainType=="NAP",1,0),
                                  ChestPainTypeTA=ifelse(ChestPainType=="TA",1,0),RestingECGNormal=ifelse(RestingECG=="Normal",1,0),
                                  RestingECGST=ifelse(RestingECG=="ST",1,0),ExerciseAngina=ifelse(ExerciseAngina=="Y",1,0),
                                  ST_SlopeFlat=ifelse(ST_Slope=="Flat",1,0),ST_SlopeUp=ifelse(ST_Slope=="Up",1,0))
final_test_dummy=final_test_dummy[,-c(3,7,11)]
k_cv
#k già selezionato prima
knn_cv=knn(big_training_dummy[,-8],k=k_cv,test=final_test_dummy[,-9],cl=big_training_dummy$HeartDisease)
confmat_knn_cv=table(knn_cv,final_test_dummy$HeartDisease)
metriche(confmat_knn_cv)




#accuracy test vs accuracy training per vedere se overfitting

knn.train.pred=knn(big_training_dummy[,-8],k=k_cv,test=big_training_dummy[,-8],cl=big_training_dummy$HeartDisease)
confmat_knn_train<-table(knn.train.pred , big_training_dummy$HeartDisease)
metriche(confmat_knn_train)[1]
metriche(confmat_knn_cv)[1]
confmat_train_stepreglog=table(round(predict(step_reg_log,big_training,type="response")),big_training$HeartDisease)
metriche(confmat_train_stepreglog)[1]
metriche(confmat_test_stepreglog)[1]

#knn forse in overfit, comunque reg log migliore quindi porto avanti solo quella
#metriche valutaz sul test

prediction_test_reglog=prediction(predict(step_reg_log,final.test,type="response"),final.test$HeartDisease)
performance(prediction_test_reglog,measure="auc")@y.values
plot(performance(prediction_test_reglog,"tpr","fpr"),main="ROC logistic regression")

prediction_knn=prediction(as.numeric(knn_cv)-1,final.test$HeartDisease)
performance(prediction_knn,measure="auc")@y.values
plot(performance(prediction_knn,"tpr","fpr"),main="ROC knn")

#il modello migliore è la regressione logistica, anche se entrambi sono validi


#variaz test e training error al variare soglia
test_train_error=function(soglia,modello=step_reg_log){
  trpred=predict(modello,big_training,type="response")
  testpred=predict(modello,final.test,type="response")
  ptr=ifelse(trpred>soglia,yes = 1,no=0)
  ptest=ifelse(testpred>soglia,1,0)
  met1=1-metriche(table(ptr,big_training$HeartDisease))[1]
  met2=1-metriche(table(ptest,final.test$HeartDisease))[1]
  met3=1-metriche(table(ptr,big_training$HeartDisease))[2]
  met4=1-metriche(table(ptest,final.test$HeartDisease))[2]
  met=rbind(met1,met2,met3,met4)
  rownames(met)=c("training error","test error", "training fp rate","test fp rate")
  colnames(met)="Metriche"
  return(met)
}

test_train_error=Vectorize(test_train_error,vectorize.args = "soglia")
soglie=seq(0.1,0.9,by=0.05)
errors_soglie=as.data.frame(cbind(soglie,t(test_train_error(soglie))))
colnames(errors_soglie)=c("soglie","training_error","test_error","training_fp_rate","test_fp_rate")
ggplot(errors_soglie,aes(x=soglie))+
  geom_point(aes(y=training_error,col="Training Error Rate"))+
  geom_line(aes(y=training_error,col="Training Error Rate"))+
  geom_point(aes(y=test_error,col="Test Error Rate"))+
  geom_point(aes(y=training_fp_rate,col="Training False Positive Rate"))+
  geom_point(aes(y=test_fp_rate,col="Test False Positive Rate"))+
  geom_line(aes(y=test_error,col="Test Error Rate"))+
  geom_line(aes(y=training_fp_rate,col="Training False Positive Rate"))+
  geom_line(aes(y=test_fp_rate,col="Test False Positive Rate"))+
  ylab("Error Rate")+
  theme_minimal()+
  scale_color_manual(name="",breaks=c("Training Error Rate","Test Error Rate","Training False Positive Rate","Test False Positive Rate"),values=c("Training Error Rate"="red","Test Error Rate"="blue",
                                                                  "Training False Positive Rate"="orange","Test False Positive Rate"="cornflowerblue"))

metriche_soglia=cbind(errors_soglie[c(6,9),1],1-errors_soglie[c(6,9),-1])
colnames(metriche_soglia)=c("soglia","training_accuracy","test_accuracy","training_sensitivity","test_sensitivity")
metriche_soglia
#soglia migliore risulta 0.35, guardando dal grafico si nota
#che abbiamo un'accuray migliore nel test ma anche 
#sensitivty migliore
prob_reg_log=predict(step_reg_log,final.test,type="response")
pred_reg_log_soglia=ifelse(prob_reg_log>0.35,yes = 1,no=0)
confmat_best_soglia=table(pred_reg_log_soglia,final.test$HeartDisease)
f1_score(confmat_best_soglia)
f1_score(confmat_test_stepreglog)
#pure f1 score migliora
prob_reglog=predict(step_logistic_mod,validation,type = "response")
logloss_reglog=-mean(step_logistic_mod_pred*log(prob_reglog)+(1-step_logistic_mod_pred)*log(1-prob_reglog))
logloss_reglog_soglia=-mean(pred_reg_log_soglia*log(prob_reg_log)+(1-pred_reg_log_soglia)*log(1-prob_reg_log))
logloss_reglog
logloss_reglog_soglia
#pure log-loss migliora
#risultato finale mod
summary(step_reg_log)




