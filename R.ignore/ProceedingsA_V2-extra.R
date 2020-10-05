########################################################################################################
##################################    Define some functions     ########################################
########################################################################################################
#Get and plot the training results -------------------------
train.sig<-function(modern_climate,fit,nsig, fit_t,nsig_t, fit_f,nsig_f, fit_tf,nsig_tf){
  train<-cbind.data.frame(modern_climate,fit[["fit"]][,nsig],fit_t[["fit"]][,nsig_t],
                          fit_f[["fit"]][,nsig_f],fit_tf[["fit"]][,nsig_tf])
  colnames(train)<-c("x","WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")
  return(train)
}
plot.train.sig<-function(train_sig_Tmin,train_sig_gdd,train_sig_alpha){
  if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
  if(!require(egg)){ install.packages("egg");library(egg)}
  
  #MTCO
  max_Tmin<-max(train_sig_Tmin);min_Tmin<-min(train_sig_Tmin)
  p1_Tmin<-ggplot(train_sig_Tmin,aes(x,WAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_Tmin, x =min_Tmin,label="(a)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression(paste("WA-PLS  MTCO"," ", (degree~C))), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))
  
  p2_Tmin<-ggplot(train_sig_Tmin,aes(x,TWAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_Tmin, x =min_Tmin,label="(b)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression(paste("TWA-PLS  MTCO"," ", (degree~C))), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))
  
  p3_Tmin<-ggplot(train_sig_Tmin,aes(x,WAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_Tmin, x =min_Tmin,label="(c)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression(paste("WA-PLS with fx  MTCO"," ", (degree~C))), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))
  
  p4_Tmin<-ggplot(train_sig_Tmin,aes(x,TWAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_Tmin, x =min_Tmin,label="(d)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression(paste("TWA-PLS with fx  MTCO"," ", (degree~C))), x = expression(paste("MTCO"," ", (degree~C))))+
    scale_y_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))
  
  
  
  #GDD0
  max_gdd<-max(train_sig_gdd);min_gdd<-min(train_sig_gdd)
  p1_gdd<-ggplot(train_sig_gdd,aes(x,WAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_gdd, x =min_gdd,label="(e)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= bquote('WA-PLS  '~ GDD[0]), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))
  
  p2_gdd<-ggplot(train_sig_gdd,aes(x,TWAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_gdd, x =min_gdd,label="(f)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= bquote('TWA-PLS  '~ GDD[0]), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))
  
  p3_gdd<-ggplot(train_sig_gdd,aes(x,WAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_gdd, x =min_gdd,label="(g)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= bquote('WA-PLS with fx  '~ GDD[0]), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))
  
  p4_gdd<-ggplot(train_sig_gdd,aes(x,TWAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_gdd, x =min_gdd,label="(h)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= bquote('TWA-PLS with fx  '~ GDD[0]), x = bquote(GDD[0]))+
    scale_y_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))
  
  
  #alpha
  max_alpha<-max(train_sig_alpha);min_alpha<-min(train_sig_alpha)
  p1_alpha<-ggplot(train_sig_alpha,aes(x,WAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_alpha, x =min_alpha,label="(i)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression("WA-PLS   "*alpha), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))
  
  p2_alpha<-ggplot(train_sig_alpha,aes(x,TWAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_alpha, x =min_alpha,label="(j)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression("TWA-PLS   "*alpha), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))
  
  p3_alpha<-ggplot(train_sig_alpha,aes(x,WAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_alpha, x =min_alpha,label="(k)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression("WA-PLS with fx   "*alpha), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))
  
  p4_alpha<-ggplot(train_sig_alpha,aes(x,TWAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_alpha, x =min_alpha,label="(l)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression("TWA-PLS with fx   "*alpha), x = expression(alpha))+
    scale_y_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))
  
  p<-ggarrange(p1_Tmin,p1_gdd,p1_alpha,
               p2_Tmin,p2_gdd,p2_alpha,
               p3_Tmin,p3_gdd,p3_alpha,
               p4_Tmin,p4_gdd,p4_alpha,  ncol = 3)
  ggsave(file="Training results.jpeg",p,width=9.5,height=10)
  
}
plot.train.residual.sig<-function(train_sig_Tmin,train_sig_gdd,train_sig_alpha){
  if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
  if(!require(egg)){ install.packages("egg");library(egg)}
  
  #Get the residuals
  train_residual_sig_Tmin<-cbind.data.frame(train_sig_Tmin[,"x"],train_sig_Tmin[,-1]-train_sig_Tmin[,"x"])
  colnames(train_residual_sig_Tmin)<-c("x","WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")
  train_residual_sig_gdd<-cbind.data.frame(train_sig_gdd[,"x"],train_sig_gdd[,-1]-train_sig_gdd[,"x"])
  colnames(train_residual_sig_gdd)<-c("x","WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")
  train_residual_sig_alpha<-cbind.data.frame(train_sig_alpha[,"x"],train_sig_alpha[,-1]-train_sig_alpha[,"x"])
  colnames(train_residual_sig_alpha)<-c("x","WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")
  
  #MTCO
  max_x_Tmin<-max(train_residual_sig_Tmin[,1]);min_x_Tmin<-min(train_residual_sig_Tmin[,1])
  max_y_Tmin<-max(train_residual_sig_Tmin[,-1]);min_y_Tmin<-min(train_residual_sig_Tmin[,-1])
  
  p1_Tmin<-ggplot(train_residual_sig_Tmin,aes(x,WAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_Tmin, x = min_x_Tmin,label="(a)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression(paste("WA-PLS  MTCO"," ", (degree~C))), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_Tmin,max_y_Tmin),breaks=c(-30,-20,-10,0,10,20,30),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_Tmin,max_x_Tmin),breaks=c(-40,-30,-20,-10,0,10,20),labels = function(x) sprintf("%g", x))
  
  p2_Tmin<-ggplot(train_residual_sig_Tmin,aes(x,TWAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_Tmin, x = min_x_Tmin,label="(b)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression(paste("TWA-PLS  MTCO"," ", (degree~C))), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_Tmin,max_y_Tmin),breaks=c(-30,-20,-10,0,10,20,30),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_Tmin,max_x_Tmin),breaks=c(-40,-30,-20,-10,0,10,20),labels = function(x) sprintf("%g", x))
  
  p3_Tmin<-ggplot(train_residual_sig_Tmin,aes(x,WAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_Tmin, x = min_x_Tmin,label="(c)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression(paste("WA-PLS with fx  MTCO"," ", (degree~C))), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_Tmin,max_y_Tmin),breaks=c(-30,-20,-10,0,10,20,30),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_Tmin,max_x_Tmin),breaks=c(-40,-30,-20,-10,0,10,20),labels = function(x) sprintf("%g", x))
  
  p4_Tmin<-ggplot(train_residual_sig_Tmin,aes(x,TWAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_Tmin, x = min_x_Tmin,label="(d)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression(paste("TWA-PLS with fx  MTCO"," ", (degree~C))), x = expression(paste("MTCO"," ", (degree~C))))+
    scale_y_continuous(limits=c(min_y_Tmin,max_y_Tmin),breaks=c(-30,-20,-10,0,10,20,30),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_Tmin,max_x_Tmin),breaks=c(-40,-30,-20,-10,0,10,20),labels = function(x) sprintf("%g", x))
  
  
  
  #GDD0
  max_x_gdd<-max(train_residual_sig_gdd[,1]);min_x_gdd<-min(train_residual_sig_gdd[,1])
  max_y_gdd<-max(train_residual_sig_gdd[,-1]);min_y_gdd<-min(train_residual_sig_gdd[,-1])
  
  p1_gdd<-ggplot(train_residual_sig_gdd,aes(x,WAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_gdd, x =min_x_gdd,label="(e)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= bquote('WA-PLS  '~ GDD[0]), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_gdd,max_y_gdd),breaks=c(-4000,-3000,-2000,-1000,0,1000,2000,3000,4000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_gdd,max_x_gdd),breaks=c(0,2000,4000,6000,8000),labels = function(x) sprintf("%g", x))
  
  p2_gdd<-ggplot(train_residual_sig_gdd,aes(x,TWAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_gdd, x =min_x_gdd,label="(f)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= bquote('TWA-PLS  '~ GDD[0]), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_gdd,max_y_gdd),breaks=c(-4000,-3000,-2000,-1000,0,1000,2000,3000,4000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_gdd,max_x_gdd),breaks=c(0,2000,4000,6000,8000),labels = function(x) sprintf("%g", x))
  
  p3_gdd<-ggplot(train_residual_sig_gdd,aes(x,WAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_gdd, x =min_x_gdd,label="(g)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= bquote('WA-PLS with fx  '~ GDD[0]), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_gdd,max_y_gdd),breaks=c(-4000,-3000,-2000,-1000,0,1000,2000,3000,4000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_gdd,max_x_gdd),breaks=c(0,2000,4000,6000,8000),labels = function(x) sprintf("%g", x))
  
  p4_gdd<-ggplot(train_residual_sig_gdd,aes(x,TWAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_gdd, x =min_x_gdd,label="(h)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= bquote('TWA-PLS with fx  '~ GDD[0]), x = bquote(GDD[0]))+
    scale_y_continuous(limits=c(min_y_gdd,max_y_gdd),breaks=c(-4000,-3000,-2000,-1000,0,1000,2000,3000,4000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_gdd,max_x_gdd),breaks=c(0,2000,4000,6000,8000),labels = function(x) sprintf("%g", x))
  
  
  #alpha
  max_x_alpha<-max(train_residual_sig_alpha[,1]);min_x_alpha<-min(train_residual_sig_alpha[,1])
  max_y_alpha<-max(train_residual_sig_alpha[,-1]);min_y_alpha<-min(train_residual_sig_alpha[,-1])
  
  p1_alpha<-ggplot(train_residual_sig_alpha,aes(x,WAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text",y= max_y_alpha, x =min_x_alpha,label="(i)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression("WA-PLS   "*alpha), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_alpha,max_y_alpha),breaks=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_x_alpha,max_x_alpha),breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2),labels = function(x) sprintf("%0.2f", x))
  
  p2_alpha<-ggplot(train_residual_sig_alpha,aes(x,TWAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_alpha, x =min_x_alpha,label="(j)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression("TWA-PLS   "*alpha), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_alpha,max_y_alpha),breaks=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_x_alpha,max_x_alpha),breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2),labels = function(x) sprintf("%0.2f", x))
  
  p3_alpha<-ggplot(train_residual_sig_alpha,aes(x,WAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_alpha, x =min_x_alpha,label="(k)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression("WA-PLS with fx   "*alpha), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_alpha,max_y_alpha),breaks=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_x_alpha,max_x_alpha),breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2),labels = function(x) sprintf("%0.2f", x))
  
  p4_alpha<-ggplot(train_residual_sig_alpha,aes(x,TWAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_alpha, x =min_x_alpha,label="(l)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression("TWA-PLS with fx   "*alpha), x = expression(alpha))+
    scale_y_continuous(limits=c(min_y_alpha,max_y_alpha),breaks=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_x_alpha,max_x_alpha),breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2),labels = function(x) sprintf("%0.2f", x))
  
  p<-ggarrange(p1_Tmin,p1_gdd,p1_alpha,
               p2_Tmin,p2_gdd,p2_alpha,
               p3_Tmin,p3_gdd,p3_alpha,
               p4_Tmin,p4_gdd,p4_alpha,  ncol = 3)
  ggsave(file="Training residuals.jpeg",p,width=9.5,height=10)
  
}

#Get and plot the core results -------------------------
core.sig<-function(fossil,nsig, fossil_t,nsig_t, fossil_f,nsig_f, fossil_tf,nsig_tf){
  core<-cbind.data.frame(fossil[["fit"]][,nsig],fossil_t[["fit"]][,nsig_t],
                         fossil_f[["fit"]][,nsig_f],fossil_tf[["fit"]][,nsig_tf])
  colnames(core)<-c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")
  return(core)
}
plot.core.sig.sse<-function(sitename, core_sig_Tmin,core_modern_Tmin, core_sig_gdd,core_modern_gdd, core_sig_alpha,core_modern_alpha,zero_Tmin,zero_gdd,zero_alpha,sse_core_sig_Tmin,sse_core_sig_gdd,sse_core_sig_alpha){
  if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
  if(!require(egg)){ install.packages("egg");library(egg)}
  if(!require(gridExtra)){install.packages("gridExtra");library(gridExtra)}
  
  xbreak<-2000*c(seq(0,6))
  
  #MTCO
  plotsite_Tmin<-cbind.data.frame(core_sig_Tmin[which(core_sig_Tmin$Site==sitename),],sse_core_sig_Tmin[which(core_sig_Tmin$Site==sitename),])
  max_Tmin<-max(plotsite_Tmin[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);max_Tmin<-max(max_Tmin,core_modern_Tmin,zero_Tmin)
  min_Tmin<-min(plotsite_Tmin[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);min_Tmin<-min(min_Tmin,core_modern_Tmin,zero_Tmin)
  
  p_Tmin<-ggplot(plotsite_Tmin)+
    geom_point(size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(aes(Age.cal.BP,WAPLS))+
    geom_ribbon(aes(ymin=WAPLS-1.96*sse_WAPLS, ymax=WAPLS+1.96*sse_WAPLS,x=Age.cal.BP),alpha=0.22,fill="black")+
    geom_point(col="red",size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(col="red",aes(Age.cal.BP,TWAPLS.fx))+
    geom_ribbon(aes(ymin=TWAPLS.fx-sse_TWAPLS.fx, ymax=TWAPLS.fx+sse_TWAPLS.fx,x=Age.cal.BP),alpha=0.36,fill="red")+
    labs(y= expression(paste("Reconstructed MTCO"," ", (degree~C))), x = NULL)+
    annotate("text", y= max_Tmin, x = (-1000),label="(a)")+ 
    geom_segment(x = (-100),xend=(-500), y = core_modern_Tmin,yend = core_modern_Tmin)+
    theme_bw()+theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_Tmin,max_Tmin),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))+
    geom_hline(yintercept = zero_Tmin,linetype="dashed")
  
  #GDD0
  plotsite_gdd<-cbind.data.frame(core_sig_gdd[which(core_sig_gdd$Site==sitename),],sse_core_sig_gdd[which(core_sig_gdd$Site==sitename),])
  max_gdd<-max(plotsite_gdd[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);max_gdd<-max(max_gdd,core_modern_gdd,zero_gdd)
  min_gdd<-min(plotsite_gdd[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);min_gdd<-min(min_gdd,core_modern_gdd,zero_gdd)
  
  p_gdd<-ggplot(plotsite_gdd)+
    geom_point(size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(aes(Age.cal.BP,WAPLS))+
    geom_ribbon(aes(ymin=WAPLS-1.96*sse_WAPLS, ymax=WAPLS+1.96*sse_WAPLS,x=Age.cal.BP),alpha=0.22,fill="black")+
    geom_point(col="red",size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(col="red",aes(Age.cal.BP,TWAPLS.fx))+
    geom_ribbon(aes(ymin=TWAPLS.fx-sse_TWAPLS.fx, ymax=TWAPLS.fx+sse_TWAPLS.fx,x=Age.cal.BP),alpha=0.36,fill="red")+
    labs(y= bquote('Reconstructed '~ GDD[0]), x = NULL)+
    annotate("text", y= max_gdd, x = (-1000),label="(b)")+ 
    geom_segment(x = (-100),xend=(-500), y = core_modern_gdd,yend = core_modern_gdd)+
    theme_bw()+theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_gdd,max_gdd),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))+
    geom_hline(yintercept = zero_gdd,linetype="dashed")
  
  #alpha
  plotsite_alpha<-cbind.data.frame(core_sig_alpha[which(core_sig_alpha$Site==sitename),],sse_core_sig_alpha[which(core_sig_alpha$Site==sitename),])
  max_alpha<-max(plotsite_alpha[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);max_alpha<-max(max_alpha,core_modern_alpha,zero_alpha)
  min_alpha<-min(plotsite_alpha[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);min_alpha<-min(min_alpha,core_modern_alpha,zero_alpha)
  
  p_alpha<-ggplot(plotsite_alpha)+
    geom_point(size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(aes(Age.cal.BP,WAPLS))+
    geom_ribbon(aes(ymin=WAPLS-1.96*sse_WAPLS, ymax=WAPLS+1.96*sse_WAPLS,x=Age.cal.BP),alpha=0.22,fill="black")+
    geom_point(col="red",size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(col="red",aes(Age.cal.BP,TWAPLS.fx))+
    geom_ribbon(aes(ymin=TWAPLS.fx-sse_TWAPLS.fx, ymax=TWAPLS.fx+sse_TWAPLS.fx,x=Age.cal.BP),alpha=0.36,fill="red")+
    labs(y= expression("Reconstructed "*alpha), x = "Age (yr BP)")+
    annotate("text", y= max_alpha, x = (-1000),label="(c)")+ 
    geom_segment(x = (-100),xend=(-500), y = core_modern_alpha,yend = core_modern_alpha)+
    theme_bw()+theme(axis.title.x=element_blank())+
    scale_y_continuous(limits=c(min_alpha,max_alpha),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))+
    geom_hline(yintercept = zero_alpha,linetype="dashed")
  
  #Range of reconstruction
  range_WAPLS<-c(max(plotsite_Tmin[,"WAPLS"],na.rm=TRUE)-min(plotsite_Tmin[,"WAPLS"],na.rm=TRUE),
                 max(plotsite_gdd[,"WAPLS"],na.rm=TRUE)-min(plotsite_gdd[,"WAPLS"],na.rm=TRUE),
                 max(plotsite_alpha[,"WAPLS"],na.rm=TRUE)-min(plotsite_alpha[,"WAPLS"],na.rm=TRUE))
  range_TWAPLS.fx<-c(max(plotsite_Tmin[,"TWAPLS.fx"],na.rm=TRUE)-min(plotsite_Tmin[,"TWAPLS.fx"],na.rm=TRUE),
                     max(plotsite_gdd[,"TWAPLS.fx"],na.rm=TRUE)-min(plotsite_gdd[,"TWAPLS.fx"],na.rm=TRUE),
                     max(plotsite_alpha[,"TWAPLS.fx"],na.rm=TRUE)-min(plotsite_alpha[,"TWAPLS.fx"],na.rm=TRUE))
  range<-rbind.data.frame(range_WAPLS,range_TWAPLS.fx)
  colnames(range)<-c("range_Tmin","range_gdd","range_alpha");  
  range$method<-c("WA-PLS","TWA-PLS with fx");
  range$method<-factor(range$method,levels=c("WA-PLS","TWA-PLS with fx"))
  
  p_Tmin_range<-ggplot(range)+geom_col(aes(method,range_Tmin),width=0.2)+ylim(0,1.2*max(range$range_Tmin))+
    geom_text(aes(method,range_Tmin),position = "identity",vjust=(-1),label=sprintf("%0.2f", round(range$range_Tmin,digits = 2)))+
    ylab(expression(paste("Range of MTCO"," ", (degree~C))))+
    annotate("text", y= 1.2*max(range$range_Tmin), x = 0.6,label="(d)")+theme_bw()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())
  
  p_gdd_range<-ggplot(range)+geom_col(aes(method,range_gdd),width=0.2)+ylim(0,1.2*max(range$range_gdd))+
    geom_text(aes(method,range_gdd),position = "identity",vjust=(-1),label=sprintf("%0.2f", round(range$range_gdd,digits = 0)))+
    ylab(bquote('Range of'~ GDD[0]))+annotate("text", y= 1.2*max(range$range_gdd), x = 0.6,label="(e)")+theme_bw()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())
  
  p_alpha_range<-ggplot(range)+geom_col(aes(method,range_alpha),width=0.2)+xlab("Method")+ylim(0,1.2*max(range$range_alpha))+
    geom_text(aes(method,range_alpha),position = "identity",vjust=(-1),label=sprintf("%0.2f", round(range$range_alpha,digits = 2)))+
    labs(y=expression("Range of "*alpha),x=NULL)+annotate("text", y= 1.2*max(range$range_alpha), x = 0.6,label="(f)")+theme_bw()
  
  p1<-ggarrange(p_Tmin,p_gdd,p_alpha,ncol=1)
  p2<-ggarrange(p_Tmin_range,p_gdd_range,p_alpha_range,ncol=1)
  
  p<-arrangeGrob(p1 ,p2, ncol = 7,nrow=1,layout_matrix = cbind(1,1,1,1,1,2,2))
  ggsave(file=paste(sitename,"reconstruction results using the last significant number of components.jpeg"),p,width=10,height=9)
  
}


########################################################################################################
#######################################    Training  ###################################################
########################################################################################################
# install.packages(c("hexSticker", "remotes", "usethis"))
install.packages("remotes")
remotes::install_github("special-uor/fxTWAPLS@v0.0.2")
if(!require(here)){install.packages("here");library(here)}
if(!require(tictoc)){install.packages("tictoc");library(tictoc)}

setwd(here::here("/"))
modern_pollen <- read.csv(system.file("extdata", 
                                      "Modern_Pollen_gdd_alpha_Tmin.csv", 
                                      package = "fxTWAPLS", 
                                      mustWork = TRUE))

taxaColMin <- which(colnames(modern_pollen) == "Abies")
taxaColMax <- which(colnames(modern_pollen) == "Zygophyllaceae")
taxa <- modern_pollen[, taxaColMin:taxaColMax]

# Training ---------------------------------------------

# Get the frequency of each climate variable fx
fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)

# In step 7 of training, climate variable is regressed to the components obtained,
# MTCO
fit_Tmin <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
fit_t_Tmin <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
fit_f_Tmin <- 
  fxTWAPLS::WAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5, usefx = TRUE, fx = fx_Tmin)
fit_tf_Tmin <- 
  fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5, usefx = TRUE, fx = fx_Tmin)

# GDD0
fit_gdd <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$gdd, nPLS = 5)
fit_t_gdd <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$gdd, nPLS = 5)
fit_f_gdd <- 
  fxTWAPLS::WAPLS.w(taxa, modern_pollen$gdd, nPLS = 4, usefx = TRUE, fx = fx_gdd)
fit_tf_gdd <- 
  fxTWAPLS::TWAPLS.w(taxa, modern_pollen$gdd, nPLS = 5, usefx = TRUE, fx = fx_gdd)

# alpha
fit_alpha <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$alpha, nPLS = 5)
fit_t_alpha <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$alpha, nPLS = 5)
fit_f_alpha <- 
  fxTWAPLS::WAPLS.w(taxa, modern_pollen$alpha, nPLS = 4, usefx = TRUE, fx = fx_alpha)
fit_tf_alpha <- 
  fxTWAPLS::TWAPLS.w(taxa, modern_pollen$alpha, nPLS = 5, usefx = TRUE, fx = fx_alpha)

#######################################################################################################################
########################## Plot the training results (Figure 2 and 3) #################################################
#######################################################################################################################

#Get and plot the results using the last significant number of components obtained for cross validation with pseudo sites removed from the training set
train_sig_Tmin <- train.sig(modern_pollen$Tmin, 
                            fit_Tmin, 3, 
                            fit_t_Tmin, 3, 
                            fit_f_Tmin, 3, 
                            fit_tf_Tmin, 4)
train_sig_gdd <- train.sig(modern_pollen$gdd,
                           fit_gdd, 2, 
                           fit_t_gdd, 2, 
                           fit_f_gdd, 2, 
                           fit_tf_gdd, 3)
train_sig_alpha<-train.sig(modern_pollen$alpha, fit_alpha,3, fit_t_alpha,4, fit_f_alpha,2, fit_tf_alpha,3)
write.csv(train_sig_Tmin,"C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Training/train_sig_Tmin.csv")
write.csv(train_sig_gdd,"C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Training/train_sig_gdd.csv")
write.csv(train_sig_alpha,"C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Training/train_sig_alpha.csv")

setwd("C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Training plots")
plot.train.sig(train_sig_Tmin,train_sig_gdd,train_sig_alpha )#Figure2
plot.train.residual.sig(train_sig_Tmin,train_sig_gdd,train_sig_alpha) #Figure3

####################Supplementary Materail plots###########################

#Get and plot the results using the same number of components obtained for cross validation with pseudo sites removed from the training set
train_same_Tmin<-train.sig(modern_pollen$Tmin, fit_Tmin,2, fit_t_Tmin,2, fit_f_Tmin,2, fit_tf_Tmin,2)
train_same_gdd<-train.sig(modern_pollen$gdd, fit_gdd,2, fit_t_gdd,2, fit_f_gdd,2, fit_tf_gdd,2)
train_same_alpha<-train.sig(modern_pollen$alpha, fit_alpha,2, fit_t_alpha,2, fit_f_alpha,2, fit_tf_alpha,2)

setwd("C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Training plots/2 components")
plot.train.sig(train_same_Tmin,train_same_gdd,train_same_alpha)
plot.train.residual.sig(train_same_Tmin,train_same_gdd,train_same_alpha)

#Get and plot the results using the same number of components obtained for cross validation with pseudo sites removed from the training set
train_same_Tmin<-train.sig(modern_pollen$Tmin, fit_Tmin,3, fit_t_Tmin,3, fit_f_Tmin,3, fit_tf_Tmin,3)
train_same_gdd<-train.sig(modern_pollen$gdd, fit_gdd,3, fit_t_gdd,3, fit_f_gdd,3, fit_tf_gdd,3)
train_same_alpha<-train.sig(modern_pollen$alpha, fit_alpha,3, fit_t_alpha,3, fit_f_alpha,3, fit_tf_alpha,3)

setwd("C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Training plots/3 components")
plot.train.sig(train_same_Tmin,train_same_gdd,train_same_alpha)
plot.train.residual.sig(train_same_Tmin,train_same_gdd,train_same_alpha)

#Get and plot the results using the same number of components obtained for cross validation with pseudo sites removed from the training set
train_same_Tmin <- train.sig(modern_pollen$Tmin, fit_Tmin,4, fit_t_Tmin,4, fit_f_Tmin,4, fit_tf_Tmin,4)
train_same_gdd<-train.sig(modern_pollen$gdd, fit_gdd,4, fit_t_gdd,4, fit_f_gdd,4, fit_tf_gdd,4)
train_same_alpha<-train.sig(modern_pollen$alpha, fit_alpha,4, fit_t_alpha,4, fit_f_alpha,4, fit_tf_alpha,4)
setwd("C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Training plots/4 components")
plot.train.sig(train_same_Tmin,train_same_gdd,train_same_alpha)
plot.train.residual.sig(train_same_Tmin,train_same_gdd,train_same_alpha)

########################################################################################################
####################################  Reconstruction  ##################################################
########################################################################################################
setwd(here::here("/"))
Holocene <- read.csv(system.file("extdata", 
                                 "Holocene.csv", 
                                 package = "fxTWAPLS", 
                                 mustWork = TRUE))
taxaColMin <- which(colnames(Holocene) == "Abies")
taxaColMax <- which(colnames(Holocene) == "Zygophyllaceae")
core<-Holocene[,taxaColMin:taxaColMax]

setwd("C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Core")

#MTCO
fossil_Tmin<-fxTWAPLS::WAPLS.predict.w(fit_Tmin,core)
fossil_t_Tmin<-fxTWAPLS::TWAPLS.predict.w(fit_t_Tmin,core)
fossil_f_Tmin<-fxTWAPLS::WAPLS.predict.w(fit_f_Tmin,core)
fossil_tf_Tmin<-fxTWAPLS::TWAPLS.predict.w(fit_tf_Tmin,core)

core_sig_Tmin<-core.sig(fossil_Tmin,3, fossil_t_Tmin,3, fossil_f_Tmin,3, fossil_tf_Tmin,4)
core_sig_Tmin<-cbind.data.frame(Holocene[,c("Site","Age.cal.BP")],core_sig_Tmin)

#GDD0
fossil_gdd<-fxTWAPLS::WAPLS.predict.w(fit_gdd,core)
fossil_t_gdd<-fxTWAPLS::TWAPLS.predict.w(fit_t_gdd,core)
fossil_f_gdd<-fxTWAPLS::WAPLS.predict.w(fit_f_gdd,core)
fossil_tf_gdd<-fxTWAPLS::TWAPLS.predict.w(fit_tf_gdd,core)

core_sig_gdd<-core.sig(fossil_gdd,2, fossil_t_gdd,2, fossil_f_gdd,2, fossil_tf_gdd,3)
core_sig_gdd<-cbind.data.frame(Holocene[,c("Site","Age.cal.BP")],core_sig_gdd)

#alpha
fossil_alpha<-fxTWAPLS::WAPLS.predict.w(fit_alpha,core)
fossil_t_alpha<-fxTWAPLS::TWAPLS.predict.w(fit_t_alpha,core)
fossil_f_alpha<-fxTWAPLS::WAPLS.predict.w(fit_f_alpha,core)
fossil_tf_alpha<-fxTWAPLS::TWAPLS.predict.w(fit_tf_alpha,core)

core_sig_alpha<-core.sig(fossil_alpha,3, fossil_t_alpha,4, fossil_f_alpha,2, fossil_tf_alpha,3)
core_sig_alpha<-cbind.data.frame(Holocene[,c("Site","Age.cal.BP")],core_sig_alpha)

write.csv(core_sig_Tmin,"C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Core/core_sig_Tmin.csv")
write.csv(core_sig_gdd,"C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Core/core_sig_gdd.csv")
write.csv(core_sig_alpha,"C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Core/core_sig_alpha.csv")

#Get the sample specific errors, use nboot=1000
#MTCO
sse_Tmin_WAPLS<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$Tmin,fossil_taxa=core,trainfun=fxTWAPLS::WAPLS.w,predictfun=fxTWAPLS::WAPLS.predict.w,nboot=1000,nPLS=5,nsig=3,usefx=FALSE,fx=NA)
sse_Tmin_TWAPLS<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$Tmin,fossil_taxa=core,trainfun=fxTWAPLS::TWAPLS.w,predictfun=fxTWAPLS::TWAPLS.predict.w,nboot=1000,nPLS=5,nsig=3,usefx=FALSE,fx=NA)
sse_Tmin_WAPLS.fx<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$Tmin,fossil_taxa=core,trainfun=fxTWAPLS::WAPLS.w,predictfun=fxTWAPLS::WAPLS.predict.w,nboot=1000,nPLS=5,nsig=3,usefx=TRUE,fx=fx_Tmin)
sse_Tmin_TWAPLS.fx<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$Tmin,fossil_taxa=core,trainfun=fxTWAPLS::TWAPLS.w,predictfun=fxTWAPLS::TWAPLS.predict.w,nboot=1000,nPLS=5,nsig=4,usefx=TRUE,fx=fx_Tmin)

sse_core_sig_Tmin<-cbind.data.frame(sse_Tmin_WAPLS,sse_Tmin_TWAPLS,sse_Tmin_WAPLS.fx,sse_Tmin_TWAPLS.fx)
colnames(sse_core_sig_Tmin)<-c("sse_WAPLS","sse_TWAPLS","sse_WAPLS.fx","sse_TWAPLS.fx")
write.csv(sse_core_sig_Tmin,"C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Core/sse_core_sig_Tmin.csv")

#GDD0
sse_gdd_WAPLS<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$gdd,fossil_taxa=core,trainfun=fxTWAPLS::WAPLS.w,predictfun=fxTWAPLS::WAPLS.predict.w,nboot=1000,nPLS=5,nsig=2,usefx=FALSE,fx=NA)
sse_gdd_TWAPLS<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$gdd,fossil_taxa=core,trainfun=fxTWAPLS::TWAPLS.w,predictfun=fxTWAPLS::TWAPLS.predict.w,nboot=1000,nPLS=5,nsig=2,usefx=FALSE,fx=NA)
sse_gdd_WAPLS.fx<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$gdd,fossil_taxa=core,trainfun=fxTWAPLS::WAPLS.w,predictfun=fxTWAPLS::WAPLS.predict.w,nboot=1000,nPLS=4,nsig=2,usefx=TRUE,fx=fx_gdd)
sse_gdd_TWAPLS.fx<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$gdd,fossil_taxa=core,trainfun=fxTWAPLS::TWAPLS.w,predictfun=fxTWAPLS::TWAPLS.predict.w,nboot=1000,nPLS=5,nsig=3,usefx=TRUE,fx=fx_gdd)

sse_core_sig_gdd<-cbind.data.frame(sse_gdd_WAPLS,sse_gdd_TWAPLS,sse_gdd_WAPLS.fx,sse_gdd_TWAPLS.fx)
colnames(sse_core_sig_gdd)<-c("sse_WAPLS","sse_TWAPLS","sse_WAPLS.fx","sse_TWAPLS.fx")
write.csv(sse_core_sig_gdd,"C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Core/sse_core_sig_gdd.csv")

#alpha
sse_alpha_WAPLS<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$alpha,fossil_taxa=core,trainfun=fxTWAPLS::WAPLS.w,predictfun=fxTWAPLS::WAPLS.predict.w,nboot=1000,nPLS=5,nsig=3,usefx=FALSE,fx=NA)
sse_alpha_TWAPLS<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$alpha,fossil_taxa=core,trainfun=fxTWAPLS::TWAPLS.w,predictfun=fxTWAPLS::TWAPLS.predict.w,nboot=1000,nPLS=5,nsig=4,usefx=FALSE,fx=NA)
sse_alpha_WAPLS.fx<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$alpha,fossil_taxa=core,trainfun=fxTWAPLS::WAPLS.w,predictfun=fxTWAPLS::WAPLS.predict.w,nboot=1000,nPLS=4,nsig=2,usefx=TRUE,fx=fx_alpha)
sse_alpha_TWAPLS.fx<-fxTWAPLS::sse.sample(modern_taxa=taxa,modern_climate=modern_pollen$alpha,fossil_taxa=core,trainfun=fxTWAPLS::TWAPLS.w,predictfun=fxTWAPLS::TWAPLS.predict.w,nboot=1000,nPLS=5,nsig=3,usefx=TRUE,fx=fx_alpha)

sse_core_sig_alpha<-cbind.data.frame(sse_alpha_WAPLS,sse_alpha_TWAPLS,sse_alpha_WAPLS.fx,sse_alpha_TWAPLS.fx)
colnames(sse_core_sig_alpha)<-c("sse_WAPLS","sse_TWAPLS","sse_WAPLS.fx","sse_TWAPLS.fx")
write.csv(sse_core_sig_alpha,"C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Core/sse_core_sig_alpha.csv")


#######################################################################################################################
########################## Plot the reconstruction results (Figure 5 and 6) ###########################################
#######################################################################################################################

setwd("C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Core plots")
#WA-PLS zero compression points
zero_Tmin<-(-2.1);zero_gdd<-3291;zero_alpha<-0.79

#Plot Basa de la Mora (Figure 5)
sitename<-"Basa de la Mora"
core_modern_Tmin<-(-2.4);core_modern_gdd<-1849;core_modern_alpha<-1.04
plot.core.sig.sse(sitename, core_sig_Tmin,core_modern_Tmin, core_sig_gdd,core_modern_gdd, core_sig_alpha,core_modern_alpha,zero_Tmin,zero_gdd,zero_alpha,sse_core_sig_Tmin,sse_core_sig_gdd,sse_core_sig_alpha)

#Plot Estanya (Figure 6)
sitename<-"Estanya"
core_modern_Tmin<-4.2;core_modern_gdd<-4352;core_modern_alpha<-0.71
plot.core.sig.sse(sitename, core_sig_Tmin,core_modern_Tmin, core_sig_gdd,core_modern_gdd, core_sig_alpha,core_modern_alpha,zero_Tmin,zero_gdd,zero_alpha,sse_core_sig_Tmin,sse_core_sig_gdd,sse_core_sig_alpha)

########################################################################################################
########################    Maps to show the training data  (Figure 1)  ################################
########################################################################################################
# Maps -------------------------------------------------------------------
#Modern sites
if(!require(ggmap)){ install.packages("ggmap");library(ggmap)}
if(!require(ggsn)){ install.packages("ggsn");library(ggsn)}
if(!require(maps)){ install.packages("maps");library(maps)}
if(!require(mapdata)){ install.packages("mapdata");library(mapdata)}
if(!require(sf)){ install.packages("sf");library(sf)}

setwd("C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Training plots")
modern_pollen <- read.csv(system.file("extdata", 
                                      "Modern_Pollen_gdd_alpha_Tmin.csv", 
                                      package = "fxTWAPLS", 
                                      mustWork = TRUE),
                          row.names = 1)

world <- map_data("world") 
minLong<-min(modern_pollen$Long);maxLong<-max(modern_pollen$Long);
minLat<-min(modern_pollen$Lat);maxLat<-max(modern_pollen$Lat);

region<-world[which(world$long>minLong&world$long<maxLong&world$lat>minLat&world$lat<maxLat),]

xat <- pretty(region$long)
yat <- pretty(region$lat)

if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(egg)){ install.packages("egg");library(egg)}

p1<-ggplot() + geom_polygon(data = region,aes(x=long, y = lat, group = group),fill='gray90',color='gray50') + theme_bw()+
  geom_point(data = modern_pollen, aes(x = Long, y = Lat), color = "black", size = 1)+
  annotate("text", y= maxLat, x = minLong,label="(a)",size=8)+
  scalebar(region, dist = 1000, dist_unit = "km",transform = TRUE, model = "WGS84",st.size = 3)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x=element_blank(),axis.text=element_text(size=18))+
  scale_y_continuous(breaks = yat, labels = paste0(yat,'?N')) 

p2<-ggplot() + geom_polygon(data = region, aes(x=long, y = lat, group = group),fill='gray90',color='gray50') + theme_bw()+
  geom_point(data = modern_pollen, aes(x = Long, y = Lat, color= Tmin), size = 1)+
  scale_colour_gradientn(colours=c("blue","dodgerblue2","khaki","orangered"))+
  annotate("text", y= maxLat, x = minLong,label="(b)",size=8)+
  scalebar(region, dist = 1000, dist_unit = "km",transform = TRUE, model = "WGS84",st.size = 3)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x=element_blank(),axis.text.y=element_blank(),axis.text=element_text(size=18))+
  theme(legend.position = c(0.85,0.85),legend.title = element_blank(),legend.text = element_text(size=12))

p3<-ggplot() + geom_polygon(data = region, aes(x=long, y = lat, group = group),fill='gray90',color='gray50') + theme_bw()+
  geom_point(data = modern_pollen, aes(x = Long, y = Lat, color= gdd), size = 1)+
  scale_colour_gradientn(colours=c("dodgerblue2","khaki","orangered"))+
  annotate("text", y= maxLat, x = minLong,label="(c)",size=8)+
  scalebar(region, dist = 1000, dist_unit = "km",transform = TRUE, model = "WGS84",st.size = 3)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text=element_text(size=18))+
  scale_y_continuous(breaks = yat, labels = paste0(yat,'?N')) +
  scale_x_continuous(breaks = xat, labels = paste0(xat,'?E')) +
  theme(legend.position = c(0.85,0.85),legend.title = element_blank(),legend.text = element_text(size=12))

p4<-ggplot() + geom_polygon(data = region, aes(x=long, y = lat, group = group),fill='gray90',color='gray50') + theme_bw()+
  geom_point(data = modern_pollen, aes(x = Long, y = Lat, color= alpha), size = 1)+
  scale_colour_gradientn(colours=c("orangered","khaki","dodgerblue2"))+
  annotate("text", y= maxLat, x = minLong,label="(d)",size=8)+
  scalebar(region, dist = 1000, dist_unit = "km",transform = TRUE, model = "WGS84",st.size = 3)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.text=element_text(size=18))+
  scale_x_continuous(breaks = xat, labels = paste0(xat,'?E')) +
  theme(legend.position = c(0.85,0.85),legend.title = element_blank(),legend.text = element_text(size=12))

p<-ggarrange(p1,p2,p3,p4, ncol = 2)

ggsave(file=paste("Map.jpeg"),p,width=14.5,height=11.5)

########################################################################################################
########################    Plot to show the compression principle (Figure 4)    #######################
########################################################################################################
setwd("C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Principle")
#Plots to show the principle

#u1>u2, t1>t2
a1_a<-log(0.6);u1_a<-10;t1_a<-8
a2_a<-log(0.6);u2_a<-(-6);t2_a<-4

x<-2;P1_a<-exp(a1_a-(x-u1_a)^2/(2*t1_a^2));P2_a<-exp(a2_a-(x-u2_a)^2/(2*t2_a^2))

pred_x<-(P1_a*u1_a+P2_a*u2_a)/(P1_a+P2_a);pred_xt<-(P1_a*u1_a/(t1_a^2)+P2_a*u2_a/(t2_a^2))/(P1_a/(t1_a^2)+P2_a/(t2_a^2))
SE<-sqrt(1/sum(P1_a/(t1_a^2)+P2_a/(t2_a^2)))

p1<-ggplot(data.frame(x=c(-30,30)), aes(x)) + ylim(0,0.7)+
  stat_function(fun=function(x) exp(a1_a-(x-u1_a)^2/(2*t1_a^2)))+
  stat_function(fun=function(x) exp(a2_a-(x-u2_a)^2/(2*t2_a^2)))+
  geom_vline(xintercept=pred_x,linetype="dashed")+geom_vline(xintercept=pred_xt,linetype="dashed")+
  annotate("text", x=u1_a,y=0.02,label = deparse(bquote(u[1])),parse = T)+
  annotate("text", x=u2_a,y=0.02,label = deparse(bquote(u[2])),parse = T)+
  annotate("text", x=pred_x+4.5,y=0.65,label = deparse(bquote(hat(x)[iWA])),parse = T)+
  annotate("text", x=pred_xt-3,y=0.65,label = deparse(bquote(hat(x)[i])),parse = T)+
  labs(subtitle = expression(paste("(a) ", u[1]," > ",u[2]," , ",t[1]," > ",t[2])))+theme_void()

#u1<u2, t1<t2
a2_b<-log(0.6);u2_b<-10;t2_b<-8
a1_b<-log(0.6);u1_b<-(-6);t1_b<-4

x<-2;P1_b<-exp(a1_b-(x-u1_b)^2/(2*t1_b^2));P2_b<-exp(a2_b-(x-u2_b)^2/(2*t2_b^2))

pred_x<-(P1_b*u1_b+P2_b*u2_b)/(P1_b+P2_b);pred_xt<-(P1_b*u1_b/(t1_b^2)+P2_b*u2_b/(t2_b^2))/(P1_b/(t1_b^2)+P2_b/(t2_b^2))
SE<-sqrt(1/sum(P1_b/(t1_b^2)+P2_b/(t2_b^2)))

p2<-ggplot(data.frame(x=c(-30,30)), aes(x)) + ylim(0,0.7)+
  stat_function(fun=function(x) exp(a1_b-(x-u1_b)^2/(2*t1_b^2)))+
  stat_function(fun=function(x) exp(a2_b-(x-u2_b)^2/(2*t2_b^2)))+
  geom_vline(xintercept=pred_x,linetype="dashed")+geom_vline(xintercept=pred_xt,linetype="dashed")+
  annotate("text", x=u1_b,y=0.02,label = deparse(bquote(u[1])),parse = T)+
  annotate("text", x=u2_b,y=0.02,label = deparse(bquote(u[2])),parse = T)+
  annotate("text", x=pred_x+4.5,y=0.65,label = deparse(bquote(hat(x)[iWA])),parse = T)+
  annotate("text", x=pred_xt-3,y=0.65,label = deparse(bquote(hat(x)[i])),parse = T)+
  labs(subtitle = expression(paste("(b) ", u[1]," < ",u[2]," , ",t[1]," < ",t[2])))+theme_void()

#u1>u2, t1<t2
a1_d<-log(0.6);u1_d<-10;t1_d<-4
a2_d<-log(0.6);u2_d<-(-6);t2_d<-8

x<-4;P1_d<-exp(a1_d-(x-u1_d)^2/(2*t1_d^2));P2_d<-exp(a2_d-(x-u2_d)^2/(2*t2_d^2))

pred_x<-(P1_d*u1_d+P2_d*u2_d)/(P1_d+P2_d);pred_xt<-(P1_d*u1_d/(t1_d^2)+P2_d*u2_d/(t2_d^2))/(P1_d/(t1_d^2)+P2_d/(t2_d^2))
SE<-sqrt(1/sum(P1_d/(t1_d^2)+P2_d/(t2_d^2)))

p3<-ggplot(data.frame(x=c(-30,30)), aes(x)) + ylim(0,0.7)+
  stat_function(fun=function(x) exp(a1_d-(x-u1_d)^2/(2*t1_d^2)))+
  stat_function(fun=function(x) exp(a2_d-(x-u2_d)^2/(2*t2_d^2)))+
  geom_vline(xintercept=pred_x,linetype="dashed")+geom_vline(xintercept=pred_xt,linetype="dashed")+
  annotate("text", x=u1_d,y=0.02,label = deparse(bquote(u[1])),parse = T)+
  annotate("text", x=u2_d,y=0.02,label = deparse(bquote(u[2])),parse = T)+
  annotate("text", x=pred_x-4.5,y=0.65,label = deparse(bquote(hat(x)[iWA])),parse = T)+
  annotate("text", x=pred_xt+3,y=0.65,label = deparse(bquote(hat(x)[i])),parse = T)+
  labs(subtitle = expression(paste("(c) ", u[1]," > ",u[2]," , ",t[1]," < ",t[2])))+theme_void()


#u1<u2, t1>t2
a2_c<-log(0.6);u2_c<-10;t2_c<-4
a1_c<-log(0.6);u1_c<-(-6);t1_c<-8

x<-4;P1_c<-exp(a1_c-(x-u1_c)^2/(2*t1_c^2));P2_c<-exp(a2_c-(x-u2_c)^2/(2*t2_c^2))

pred_x<-(P1_c*u1_c+P2_c*u2_c)/(P1_c+P2_c);pred_xt<-(P1_c*u1_c/(t1_c^2)+P2_c*u2_c/(t2_c^2))/(P1_c/(t1_c^2)+P2_c/(t2_c^2))
SE<-sqrt(1/sum(P1_c/(t1_c^2)+P2_c/(t2_c^2)))

p4<-ggplot(data.frame(x=c(-30,30)), aes(x)) + ylim(0,0.7)+
  stat_function(fun=function(x) exp(a1_c-(x-u1_c)^2/(2*t1_c^2)))+
  stat_function(fun=function(x) exp(a2_c-(x-u2_c)^2/(2*t2_c^2)))+
  geom_vline(xintercept=pred_x,linetype="dashed")+geom_vline(xintercept=pred_xt,linetype="dashed")+
  annotate("text", x=u1_c,y=0.02,label = deparse(bquote(u[1])),parse = T)+
  annotate("text", x=u2_c,y=0.02,label = deparse(bquote(u[2])),parse = T)+
  annotate("text", x=pred_x-4.5,y=0.65,label = deparse(bquote(hat(x)[iWA])),parse = T)+
  annotate("text", x=pred_xt+3,y=0.65,label = deparse(bquote(hat(x)[i])),parse = T)+
  labs(subtitle = expression(paste("(d) ", u[1]," < ",u[2]," , ",t[1]," > ",t[2])))+theme_void()

p_principle<-ggarrange(p1,p2,p3,p4,ncol = 2)

# Optimum and tolerance ---------------------------------------------------
getut<-function(x,y){
  t<-matrix(NA,ncol(y),1); #tolerance
  
  sumk_yik<-rowSums(y)
  sumi_yik<-colSums(y)
  
  u = t(y)%*%x / sumi_yik;
  n2<-matrix(NA,ncol(y),1)
  for(k in 1:ncol(y)){
    t[k] = sqrt(sum(y[,k]*(x-u[k])^2)/sumi_yik[k])
    n2[k]<-1/sum((y[,k]/sum(y[,k]))^2)
    t[k]<-t[k]/sqrt(1-1/n2[k])
  }
  
  ut<-cbind.data.frame(u,t)
  return(ut)
}# get u and t
ut_Tmin<-getut(modern_pollen$Tmin,taxa)
ut_gdd<-getut(modern_pollen$gdd,taxa)
ut_alpha<-getut(modern_pollen$alpha,taxa)

p_Tmin<-ggplot(ut_Tmin,aes(u,t))+geom_point(size=0.8)+annotate("text", y= max(ut_Tmin$t), x = min(ut_Tmin$u),label="(e)")+ theme_bw()+
  labs(y= expression(paste("Tolerance of MTCO"," ", (degree~C))), x =expression(paste("Optimum of MTCO"," ",(degree~C))))
p_gdd<-ggplot(ut_gdd,aes(u,t))+geom_point(size=0.8)+annotate("text", y= max(ut_gdd$t), x = min(ut_gdd$u),label="(f)")+theme_bw()+
  labs(y= bquote('Tolerance of'~ GDD[0]), x = bquote('Optimum of'~ GDD[0]))
p_alpha<-ggplot(ut_alpha,aes(u,t))+geom_point(size=0.8)+annotate("text", y= max(ut_alpha$t), x = min(ut_alpha$u),label="(g)")+theme_bw()+
  labs(y= expression("Tolerance of "*alpha), x = expression("Optimum of "*alpha))

p_optimum_tolerance<-ggarrange(p_Tmin,p_gdd,p_alpha,ncol = 1)


# fx ----------------------------------------------------------------------
h_Tmin<-ggplot(modern_pollen,aes(x=Tmin))+geom_histogram(binwidth=0.02)+theme_bw()+
  labs(y= "Frequency", x =expression(paste("MTCO"," ",(degree~C))))+
  annotate("text", y= 40, x = min(modern_pollen$Tmin),label="(h)")
h_gdd<-ggplot(modern_pollen,aes(x=gdd))+geom_histogram(binwidth=20)+theme_bw()+
  labs(y= "Frequency", x =bquote(GDD[0]))+
  annotate("text", y= 110, x = min(modern_pollen$gdd),label="(i)")
h_alpha<-ggplot(modern_pollen,aes(x=alpha))+geom_histogram(binwidth=0.002)+theme_bw()+
  labs(y= "Frequency", x = expression(alpha))+
  annotate("text", y= 60 , x = min(modern_pollen$alpha),label="(j)")

h<-ggarrange(h_Tmin,h_gdd,h_alpha,ncol = 1)

# Put them together -------------------------------------------------------
if(!require(gridExtra)){install.packages("gridExtra");library(gridExtra)}
p<-arrangeGrob(p_principle,p_optimum_tolerance,h,ncol = 7,nrow=1,layout_matrix = cbind(1,1,1,2,2,3,3),padding=unit(2,"line"))

ggsave(file="The principle of compression in WA-PLS.jpeg",p,width=11,height=7)

########################################################################################################
############################################    Table 3    #############################################
########################################################################################################
sum(taxa>0)/ncol(taxa)

narrow_Tmin<-which(ut_Tmin$t<3);narrow_gdd<-which(ut_gdd$t<500);narrow_alpha<-which(ut_alpha$t<0.1)
sum(taxa[,narrow_Tmin]>0)/length(narrow_Tmin)
sum(taxa[,narrow_gdd]>0)/length(narrow_gdd)
sum(taxa[,narrow_alpha]>0)/length(narrow_alpha)


#######################################################################################################################
######################################## Residuals and elevations #####################################################
#######################################################################################################################
# Residual and elevation --------------------------------------------------
Elv<-modern_pollen$Elv
train_sig_Tmin<-train.sig(modern_pollen$Tmin, fit_Tmin,3, fit_t_Tmin,3, fit_f_Tmin,3, fit_tf_Tmin,4)
train_sig_gdd<-train.sig(modern_pollen$gdd, fit_gdd,2, fit_t_gdd,2, fit_f_gdd,2, fit_tf_gdd,3)
train_sig_alpha<-train.sig(modern_pollen$alpha, fit_alpha,3, fit_t_alpha,4, fit_f_alpha,2, fit_tf_alpha,3)

train_sig_Tmin_elv<-cbind.data.frame(train_sig_Tmin,Elv)
train_sig_gdd_elv<-cbind.data.frame(train_sig_gdd,Elv)
train_sig_alpha_elv<-cbind.data.frame(train_sig_alpha,Elv)

setwd("C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Training plots")
if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(egg)){ install.packages("egg");library(egg)}
#MTCO
p1_Tmin<-ggplot(data=train_sig_Tmin_elv)+geom_point(aes(Elv,WAPLS-x),size=0.8)+theme_bw()+
  annotate("text", y= 40, x = -500,label="(a)")+ylim(-40,40)+
  labs(y= expression(paste("WA-PLS MTCO"," ", (degree~C))), x = NULL)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())
p2_Tmin<-ggplot(data=train_sig_Tmin_elv)+geom_point(aes(Elv,TWAPLS-x),size=0.8)+theme_bw()+
  annotate("text", y= 40, x = -500,label="(b)")+ylim(-40,40)+
  labs(y= expression(paste("TWA-PLS MTCO"," ", (degree~C))), x = NULL)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())
p3_Tmin<-ggplot(data=train_sig_Tmin_elv)+geom_point(aes(Elv,WAPLS.fx-x),size=0.8)+theme_bw()+
  annotate("text", y= 40, x = -500,label="(c)")+ylim(-40,40)+
  labs(y= expression(paste("WA-PLS with fx MTCO"," ", (degree~C))), x = NULL)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())
p4_Tmin<-ggplot(data=train_sig_Tmin_elv)+geom_point(aes(Elv,TWAPLS.fx-x),size=0.8)+theme_bw()+
  annotate("text", y= 40, x = -500,label="(d)")+ylim(-40,40)+
  labs(y= expression(paste("TWA-PLS with fx MTCO"," ", (degree~C))), x = "Elevation (m)")

#gdd
p1_gdd<-ggplot(data=train_sig_gdd_elv)+geom_point(aes(Elv,WAPLS-x),size=0.8)+theme_bw()+
  annotate("text", y= 6000, x = -500,label="(e)")+ylim(-6000,6000)+
  labs(y= bquote('WA-PLS'~ GDD[0]), x = NULL)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())
p2_gdd<-ggplot(data=train_sig_gdd_elv)+geom_point(aes(Elv,TWAPLS-x),size=0.8)+theme_bw()+
  annotate("text", y= 6000, x = -500,label="(f)")+ylim(-6000,6000)+
  labs(y= bquote('TWA-PLS'~ GDD[0]), x = NULL)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())
p3_gdd<-ggplot(data=train_sig_gdd_elv)+geom_point(aes(Elv,WAPLS.fx-x),size=0.8)+theme_bw()+
  annotate("text", y= 6000, x = -500,label="(g)")+ylim(-6000,6000)+
  labs(y= bquote('WA-PLS with fx'~ GDD[0]), x = NULL)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())
p4_gdd<-ggplot(data=train_sig_gdd_elv)+geom_point(aes(Elv,TWAPLS.fx-x),size=0.8)+theme_bw()+
  annotate("text", y= 6000, x = -500,label="(h)")+ylim(-6000,6000)+
  labs(y= bquote('TWA-PLS with fx'~ GDD[0]), x = "Elevation (m)")

#alpha
p1_alpha<-ggplot(data=train_sig_alpha_elv)+geom_point(aes(Elv,WAPLS-x),size=0.8)+theme_bw()+
  annotate("text", y= 1, x = -500,label="(i)")+ylim(-1,1)+
  labs(y= expression("WA-PLS "*alpha), x = NULL)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())
p2_alpha<-ggplot(data=train_sig_alpha_elv)+geom_point(aes(Elv,TWAPLS-x),size=0.8)+theme_bw()+
  annotate("text", y= 1, x = -500,label="(j)")+ylim(-1,1)+
  labs(y= expression("TWA-PLS "*alpha), x = NULL)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())
p3_alpha<-ggplot(data=train_sig_alpha_elv)+geom_point(aes(Elv,WAPLS.fx-x),size=0.8)+theme_bw()+
  annotate("text", y= 1, x = -500,label="(k)")+ylim(-1,1)+
  labs(y= expression("WA-PLS with fx "*alpha), x = NULL)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())
p4_alpha<-ggplot(data=train_sig_alpha_elv)+geom_point(aes(Elv,TWAPLS.fx-x),size=0.8)+theme_bw()+
  annotate("text", y= 1, x = -500,label="(l)")+ylim(-1,1)+
  labs(y= expression("TWA-PLS with fx "*alpha), x = "Elevation (m)")


p<-ggarrange(p1_Tmin,p1_gdd,p1_alpha,
             p2_Tmin,p2_gdd,p2_alpha,
             p3_Tmin,p3_gdd,p3_alpha,
             p4_Tmin,p4_gdd,p4_alpha,  ncol = 3)
ggsave(file="Residual and elevation.jpeg",p,width=9.5,height=10)

#######################################################################################################################
############################################### non_unimodal taxa #######################################################
#######################################################################################################################

non_unimodal_taxon<-"Artemisia"

abun<-modern_pollen[,non_unimodal_taxon]

plot(abun~modern_pollen$Tmin)
plot(abun~modern_pollen$gdd)
plot(abun~modern_pollen$alpha)

if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(egg)){ install.packages("egg");library(egg)}

setwd("C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Non_unimodality")

p1<-ggplot(modern_pollen)+geom_point(size=0.8,aes(Tmin,Artemisia))+theme_bw()+
  labs(y= "Abundance of Artemisa",  x = expression(paste("MTCO"," ", (degree~C))))  
ggsave(file="Abundance of Artemisia to MTCO.jpeg",p1,width=6,height=4)



# Training with non_unimodal taxa removed ---------------------------------------------
setwd("C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Non-unimodality")

#Get the frequency of each climate variable fx
fx_Tmin<-fxTWAPLS::fx(modern_pollen$Tmin,bin=0.02)
fx_gdd<-fxTWAPLS::fx(modern_pollen$gdd,bin=20)
fx_alpha<-fxTWAPLS::fx(modern_pollen$alpha,bin=0.002)

#In step 7 of training, climate variable is regressed to the components obtained,
taxa1<-taxa[,-which(colnames(taxa)==non_unimodal_taxon)]
#MTCO
fit_Tmin1<-fxTWAPLS::WAPLS.w(taxa1,modern_pollen$Tmin,nPLS=5)
fit_t_Tmin1<-fxTWAPLS::TWAPLS.w(taxa1,modern_pollen$Tmin,nPLS=5)
fit_f_Tmin1<-fxTWAPLS::WAPLS.w(taxa1,modern_pollen$Tmin,nPLS=5,usefx=TRUE,fx=fx_Tmin)
fit_tf_Tmin1<-fxTWAPLS::TWAPLS.w(taxa1,modern_pollen$Tmin,nPLS=5,usefx=TRUE,fx=fx_Tmin)

#GDD0
fit_gdd1<-fxTWAPLS::WAPLS.w(taxa1,modern_pollen$gdd,nPLS=5)
fit_t_gdd1<-fxTWAPLS::TWAPLS.w(taxa1,modern_pollen$gdd,nPLS=5)
fit_f_gdd1<-fxTWAPLS::WAPLS.w(taxa1,modern_pollen$gdd,nPLS=4,usefx=TRUE,fx=fx_gdd)
fit_tf_gdd1<-fxTWAPLS::TWAPLS.w(taxa1,modern_pollen$gdd,nPLS=5,usefx=TRUE,fx=fx_gdd)

#alpha
fit_alpha1<-fxTWAPLS::WAPLS.w(taxa1,modern_pollen$alpha,nPLS=5)
fit_t_alpha1<-fxTWAPLS::TWAPLS.w(taxa1,modern_pollen$alpha,nPLS=5)
fit_f_alpha1<-fxTWAPLS::WAPLS.w(taxa1,modern_pollen$alpha,nPLS=4,usefx=TRUE,fx=fx_alpha)
fit_tf_alpha1<-fxTWAPLS::TWAPLS.w(taxa1,modern_pollen$alpha,nPLS=5,usefx=TRUE,fx=fx_alpha)

# Reconstruction with non_unimodal taxa removed ---------------------------------------------

core1<-core[,-which(colnames(core)==non_unimodal_taxon)]

#MTCO
fossil_Tmin1<-fxTWAPLS::WAPLS.predict.w(fit_Tmin1,core1)
fossil_t_Tmin1<-fxTWAPLS::TWAPLS.predict.w(fit_t_Tmin1,core1)
fossil_f_Tmin1<-fxTWAPLS::WAPLS.predict.w(fit_f_Tmin1,core1)
fossil_tf_Tmin1<-fxTWAPLS::TWAPLS.predict.w(fit_tf_Tmin1,core1)

core_sig_Tmin1<-core.sig(fossil_Tmin1,3, fossil_t_Tmin1,3, fossil_f_Tmin1,3, fossil_tf_Tmin1,4)
core_sig_Tmin1<-cbind.data.frame(Holocene[,c("Site","Age.cal.BP")],core_sig_Tmin1)

#GDD0
fossil_gdd1<-fxTWAPLS::WAPLS.predict.w(fit_gdd1,core1)
fossil_t_gdd1<-fxTWAPLS::TWAPLS.predict.w(fit_t_gdd1,core1)
fossil_f_gdd1<-fxTWAPLS::WAPLS.predict.w(fit_f_gdd1,core1)
fossil_tf_gdd1<-fxTWAPLS::TWAPLS.predict.w(fit_tf_gdd1,core1)

core_sig_gdd1<-core.sig(fossil_gdd1,2, fossil_t_gdd1,2, fossil_f_gdd1,2, fossil_tf_gdd1,3)
core_sig_gdd1<-cbind.data.frame(Holocene[,c("Site","Age.cal.BP")],core_sig_gdd1)

#alpha
fossil_alpha1<-fxTWAPLS::WAPLS.predict.w(fit_alpha1,core1)
fossil_t_alpha1<-fxTWAPLS::TWAPLS.predict.w(fit_t_alpha1,core1)
fossil_f_alpha1<-fxTWAPLS::WAPLS.predict.w(fit_f_alpha1,core1)
fossil_tf_alpha1<-fxTWAPLS::TWAPLS.predict.w(fit_tf_alpha1,core1)

core_sig_alpha1<-core.sig(fossil_alpha1,3, fossil_t_alpha1,4, fossil_f_alpha1,2, fossil_tf_alpha1,3)
core_sig_alpha1<-cbind.data.frame(Holocene[,c("Site","Age.cal.BP")],core_sig_alpha1)


#Get the sample specific error
#MTCO
sse_Tmin_WAPLS1<-fxTWAPLS::sse.sample(modern_taxa=taxa1,modern_climate=modern_pollen$Tmin,fossil_taxa=core1,trainfun=fxTWAPLS::WAPLS.w,predictfun=fxTWAPLS::WAPLS.predict.w,nboot=1000,nPLS=5,nsig=3,usefx=FALSE,fx=NA)
sse_Tmin_TWAPLS1<-fxTWAPLS::sse.sample(modern_taxa=taxa1,modern_climate=modern_pollen$Tmin,fossil_taxa=core1,trainfun=fxTWAPLS::TWAPLS.w,predictfun=fxTWAPLS::TWAPLS.predict.w,nboot=1000,nPLS=5,nsig=3,usefx=FALSE,fx=NA)
sse_Tmin_WAPLS.fx1<-fxTWAPLS::sse.sample(modern_taxa=taxa1,modern_climate=modern_pollen$Tmin,fossil_taxa=core1,trainfun=fxTWAPLS::WAPLS.w,predictfun=fxTWAPLS::WAPLS.predict.w,nboot=1000,nPLS=5,nsig=3,usefx=TRUE,fx=fx_Tmin)
sse_Tmin_TWAPLS.fx1<-fxTWAPLS::sse.sample(modern_taxa=taxa1,modern_climate=modern_pollen$Tmin,fossil_taxa=core1,trainfun=fxTWAPLS::TWAPLS.w,predictfun=fxTWAPLS::TWAPLS.predict.w,nboot=1000,nPLS=5,nsig=4,usefx=TRUE,fx=fx_Tmin)

sse_core_sig_Tmin1<-cbind.data.frame(sse_Tmin_WAPLS1,sse_Tmin_TWAPLS1,sse_Tmin_WAPLS.fx1,sse_Tmin_TWAPLS.fx1)
colnames(sse_core_sig_Tmin1)<-c("sse_WAPLS","sse_TWAPLS","sse_WAPLS.fx","sse_TWAPLS.fx")

#GDD0
sse_gdd_WAPLS1<-fxTWAPLS::sse.sample(modern_taxa=taxa1,modern_climate=modern_pollen$gdd,fossil_taxa=core1,trainfun=fxTWAPLS::WAPLS.w,predictfun=fxTWAPLS::WAPLS.predict.w,nboot=1000,nPLS=5,nsig=2,usefx=FALSE,fx=NA)
sse_gdd_TWAPLS1<-fxTWAPLS::sse.sample(modern_taxa=taxa1,modern_climate=modern_pollen$gdd,fossil_taxa=core1,trainfun=fxTWAPLS::TWAPLS.w,predictfun=fxTWAPLS::TWAPLS.predict.w,nboot=1000,nPLS=5,nsig=2,usefx=FALSE,fx=NA)
sse_gdd_WAPLS.fx1<-fxTWAPLS::sse.sample(modern_taxa=taxa1,modern_climate=modern_pollen$gdd,fossil_taxa=core1,trainfun=fxTWAPLS::WAPLS.w,predictfun=fxTWAPLS::WAPLS.predict.w,nboot=1000,nPLS=4,nsig=2,usefx=TRUE,fx=fx_gdd)
sse_gdd_TWAPLS.fx1<-fxTWAPLS::sse.sample(modern_taxa=taxa1,modern_climate=modern_pollen$gdd,fossil_taxa=core1,trainfun=fxTWAPLS::TWAPLS.w,predictfun=fxTWAPLS::TWAPLS.predict.w,nboot=1000,nPLS=5,nsig=3,usefx=TRUE,fx=fx_gdd)

sse_core_sig_gdd1<-cbind.data.frame(sse_gdd_WAPLS1,sse_gdd_TWAPLS1,sse_gdd_WAPLS.fx1,sse_gdd_TWAPLS.fx1)
colnames(sse_core_sig_gdd1)<-c("sse_WAPLS","sse_TWAPLS","sse_WAPLS.fx","sse_TWAPLS.fx")

#alpha
sse_alpha_WAPLS1<-fxTWAPLS::sse.sample(modern_taxa=taxa1,modern_climate=modern_pollen$alpha,fossil_taxa=core1,trainfun=fxTWAPLS::WAPLS.w,predictfun=fxTWAPLS::WAPLS.predict.w,nboot=1000,nPLS=5,nsig=3,usefx=FALSE,fx=NA)
sse_alpha_TWAPLS1<-fxTWAPLS::sse.sample(modern_taxa=taxa1,modern_climate=modern_pollen$alpha,fossil_taxa=core1,trainfun=fxTWAPLS::TWAPLS.w,predictfun=fxTWAPLS::TWAPLS.predict.w,nboot=1000,nPLS=5,nsig=4,usefx=FALSE,fx=NA)
sse_alpha_WAPLS.fx1<-fxTWAPLS::sse.sample(modern_taxa=taxa1,modern_climate=modern_pollen$alpha,fossil_taxa=core1,trainfun=fxTWAPLS::WAPLS.w,predictfun=fxTWAPLS::WAPLS.predict.w,nboot=1000,nPLS=4,nsig=2,usefx=TRUE,fx=fx_alpha)
sse_alpha_TWAPLS.fx1<-fxTWAPLS::sse.sample(modern_taxa=taxa1,modern_climate=modern_pollen$alpha,fossil_taxa=core1,trainfun=fxTWAPLS::TWAPLS.w,predictfun=fxTWAPLS::TWAPLS.predict.w,nboot=1000,nPLS=5,nsig=3,usefx=TRUE,fx=fx_alpha)

sse_core_sig_alpha1<-cbind.data.frame(sse_alpha_WAPLS1,sse_alpha_TWAPLS1,sse_alpha_WAPLS.fx1,sse_alpha_TWAPLS.fx1)
colnames(sse_core_sig_alpha1)<-c("sse_WAPLS","sse_TWAPLS","sse_WAPLS.fx","sse_TWAPLS.fx")

#plot them
setwd("C:/Users/ml4418.SPHB-LT-069/Desktop/Master Project/Data/Output data/Non-unimodality")
#WA-PLS zero compression points
if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
xbreak<-2000*c(seq(0,6))

#Plot Basa de la Mora
sitename<-"Basa de la Mora"
sitename<-"Estanya"

if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
xbreak<-2000*c(seq(0,6))

#MTCO
plotsite_Tmin<-cbind.data.frame(core_sig_Tmin[which(core_sig_Tmin$Site==sitename),],sse_core_sig_Tmin[which(core_sig_Tmin$Site==sitename),])
plotsite_Tmin1<-cbind.data.frame(core_sig_Tmin1[which(core_sig_Tmin1$Site==sitename),],sse_core_sig_Tmin1[which(core_sig_Tmin1$Site==sitename),])
max_Tmin<-max(plotsite_Tmin[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);max_Tmin<-max(max_Tmin,core_modern_Tmin)
min_Tmin<-min(plotsite_Tmin[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);min_Tmin<-min(min_Tmin,core_modern_Tmin)

p_Tmin<-ggplot()+
  geom_point(data=plotsite_Tmin,size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(data=plotsite_Tmin,aes(Age.cal.BP,WAPLS))+
  geom_ribbon(data=plotsite_Tmin,aes(ymin=WAPLS-1.96*sse_WAPLS, ymax=WAPLS+1.96*sse_WAPLS,x=Age.cal.BP),alpha=0.22,fill="black")+
  geom_point(data=plotsite_Tmin1,size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(linetype="dashed",data=plotsite_Tmin1,aes(Age.cal.BP,WAPLS))+
  geom_ribbon(data=plotsite_Tmin1,aes(ymin=WAPLS-1.96*sse_WAPLS, ymax=WAPLS+1.96*sse_WAPLS,x=Age.cal.BP),alpha=0.22,fill="black")+
  geom_point(data=plotsite_Tmin,size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(data=plotsite_Tmin,col="red",aes(Age.cal.BP,TWAPLS.fx))+
  geom_ribbon(data=plotsite_Tmin,aes(ymin=TWAPLS.fx-1.96*sse_TWAPLS.fx, ymax=TWAPLS.fx+1.96*sse_TWAPLS.fx,x=Age.cal.BP),alpha=0.22,fill="red")+
  geom_point(data=plotsite_Tmin1,size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(linetype="dashed",col="red",data=plotsite_Tmin1,aes(Age.cal.BP,TWAPLS.fx))+
  geom_ribbon(data=plotsite_Tmin1,aes(ymin=TWAPLS.fx-1.96*sse_TWAPLS.fx, ymax=TWAPLS.fx+1.96*sse_TWAPLS.fx,x=Age.cal.BP),alpha=0.22,fill="red")+
  labs(y= expression(paste("Reconstructed MTCO"," ", (degree~C))),  x = "Age (yr BP)")+
  theme_bw()+
  scale_y_continuous(limits=c(min_Tmin,max_Tmin),labels = function(x) sprintf("%g", x))+
  scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))
ggsave(file=paste(sitename,"Comparision of reconstructions with and without non_unimodal taxa.jpeg"),p_Tmin,width=10,height=6)


# Principle of observed non-unimodality -------------------------------------
x<-seq(-80,80)
a1<-0.3;u1<-(-40);t1<-5
f1<-function(x){exp(a1-(x-u1)^2/(2*t1^2))}

a2<-0.5;u2<-0;t2<-10
f2<-function(x){exp(a2-(x-u2)^2/(2*t2^2))}

a3<-0.25;u3<-0;t3<-20
f3<-function(x){exp(a3-(x-u3)^2/(2*t3^2))}

a4<-0.3;u4<-40;t4<-5
f4<-function(x){exp(a4-(x-u4)^2/(2*t4^2))}

f<-function(x){exp(a1-(x-u1)^2/(2*t1^2))+exp(a2-(x-u2)^2/(2*t2^2))+exp(a3-(x-u3)^2/(2*t3^2)+exp(a4-(x-u4)^2/(2*t4^2)))}
f1s<-function(x){ (exp(a1-(x-u1)^2/(2*t1^2)))/(exp(a1-(x-u1)^2/(2*t1^2))+exp(a2-(x-u2)^2/(2*t2^2))+exp(a3-(x-u3)^2/(2*t3^2))+exp(a4-(x-u4)^2/(2*t4^2)))}
f2s<-function(x){ (exp(a2-(x-u2)^2/(2*t2^2)))/(exp(a1-(x-u1)^2/(2*t1^2))+exp(a2-(x-u2)^2/(2*t2^2))+exp(a3-(x-u3)^2/(2*t3^2))+exp(a4-(x-u4)^2/(2*t4^2)))}
f3s<-function(x){ (exp(a3-(x-u3)^2/(2*t3^2)))/(exp(a1-(x-u1)^2/(2*t1^2))+exp(a2-(x-u2)^2/(2*t2^2))+exp(a3-(x-u3)^2/(2*t3^2))+exp(a4-(x-u4)^2/(2*t4^2)))}
f4s<-function(x){ (exp(a4-(x-u4)^2/(2*t4^2)))/(exp(a1-(x-u1)^2/(2*t1^2))+exp(a2-(x-u2)^2/(2*t2^2))+exp(a3-(x-u3)^2/(2*t3^2))+exp(a4-(x-u4)^2/(2*t4^2)))}

curve(f,from=-80,to=80,lty="dashed",xaxt='n', yaxt='n',xlab="Climate",ylab="Abundance");text(x=x[which.max(f(x))],y=max(f(x)),"sum")
curve(f1,from=-80,to=80,add=TRUE);text(x=x[which.max(f1(x))],y=max(f1(x)),"1")
curve(f2,from=-80,to=80,add=TRUE);text(x=x[which.max(f2(x))],y=max(f2(x)),"2")
curve(f3,from=-80,to=80,add=TRUE);text(x=x[which.max(f3(x))],y=max(f3(x)),"3")
curve(f4,from=-80,to=80,add=TRUE);text(x=x[which.max(f4(x))],y=max(f4(x)),"4")
curve(f3s,from=-80,to=80,add=TRUE,col="red");text(x=x[which.max(f3s(x))],y=max(f3s(x)),"3",col="red")



