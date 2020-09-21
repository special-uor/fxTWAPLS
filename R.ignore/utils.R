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
plot.core.sig<-function(sitename, core_sig_Tmin,core_modern_Tmin, core_sig_gdd,core_modern_gdd, core_sig_alpha,core_modern_alpha,zero_Tmin,zero_gdd,zero_alpha){
  if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
  if(!require(egg)){ install.packages("egg");library(egg)}
  if(!require(gridExtra)){install.packages("gridExtra");library(gridExtra)}
  
  xbreak<-2000*c(seq(0,6))
  
  #MTCO
  plotsite_Tmin<-core_sig_Tmin[which(core_sig_Tmin$Site==sitename),]
  max_Tmin<-max(plotsite_Tmin[,-c(1,2)],na.rm=TRUE);max_Tmin<-max(max_Tmin,core_modern_Tmin,zero_Tmin)
  min_Tmin<-min(plotsite_Tmin[,-c(1,2)],na.rm=TRUE);min_Tmin<-min(min_Tmin,core_modern_Tmin,zero_Tmin)
  
  p_Tmin<-ggplot(plotsite_Tmin)+
    geom_point(size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(aes(Age.cal.BP,WAPLS))+
    geom_point(col="red",size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(col="red",aes(Age.cal.BP,TWAPLS.fx))+
    labs(y= expression(paste("Reconstructed MTCO"," ", (degree~C))), x = NULL)+
    annotate("text", y= max_Tmin, x = (-1000),label="(a)")+ 
    geom_segment(x = (-100),xend=(-500), y = core_modern_Tmin,yend = core_modern_Tmin)+
    theme_bw()+theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_Tmin,max_Tmin),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))+
    geom_hline(yintercept = zero_Tmin,linetype="dashed")
  
  #GDD0
  plotsite_gdd<-core_sig_gdd[which(core_sig_gdd$Site==sitename),]
  max_gdd<-max(plotsite_gdd[,-c(1,2)],na.rm=TRUE);max_gdd<-max(max_gdd,core_modern_gdd,zero_gdd)
  min_gdd<-min(plotsite_gdd[,-c(1,2)],na.rm=TRUE);min_gdd<-min(min_gdd,core_modern_gdd,zero_gdd)
  
  p_gdd<-ggplot(plotsite_gdd)+
    geom_point(size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(aes(Age.cal.BP,WAPLS))+
    geom_point(col="red",size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(col="red",aes(Age.cal.BP,TWAPLS.fx))+
    labs(y= bquote('Reconstructed '~ GDD[0]), x = NULL)+
    annotate("text", y= max_gdd, x = (-1000),label="(b)")+ 
    geom_segment(x = (-100),xend=(-500), y = core_modern_gdd,yend = core_modern_gdd)+
    theme_bw()+theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_gdd,max_gdd),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))+
    geom_hline(yintercept = zero_gdd,linetype="dashed")
  
  #alpha
  plotsite_alpha<-core_sig_alpha[which(core_sig_alpha$Site==sitename),]
  max_alpha<-max(plotsite_alpha[,-c(1,2)],na.rm=TRUE);max_alpha<-max(max_alpha,core_modern_alpha,zero_alpha)
  min_alpha<-min(plotsite_alpha[,-c(1,2)],na.rm=TRUE);min_alpha<-min(min_alpha,core_modern_alpha,zero_alpha)
  
  p_alpha<-ggplot(plotsite_alpha)+
    geom_point(size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(aes(Age.cal.BP,WAPLS))+
    geom_point(col="red",size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(col="red",aes(Age.cal.BP,TWAPLS.fx))+
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
  #range_TWAPLS<-c(max(plotsite_Tmin[,"TWAPLS"],na.rm=TRUE)-min(plotsite_Tmin[,"TWAPLS"],na.rm=TRUE),
  #  max(plotsite_gdd[,"TWAPLS"],na.rm=TRUE)-min(plotsite_gdd[,"TWAPLS"],na.rm=TRUE),
  # max(plotsite_alpha[,"TWAPLS"],na.rm=TRUE)-min(plotsite_alpha[,"TWAPLS"],na.rm=TRUE))
  # range_WAPLS.fx<-c(max(plotsite_Tmin[,"WAPLS.fx"],na.rm=TRUE)-min(plotsite_Tmin[,"WAPLS.fx"],na.rm=TRUE),
  # max(plotsite_gdd[,"WAPLS.fx"],na.rm=TRUE)-min(plotsite_gdd[,"WAPLS.fx"],na.rm=TRUE),
  # max(plotsite_alpha[,"WAPLS.fx"],na.rm=TRUE)-min(plotsite_alpha[,"WAPLS.fx"],na.rm=TRUE))
  range_TWAPLS.fx<-c(max(plotsite_Tmin[,"TWAPLS.fx"],na.rm=TRUE)-min(plotsite_Tmin[,"TWAPLS.fx"],na.rm=TRUE),
                     max(plotsite_gdd[,"TWAPLS.fx"],na.rm=TRUE)-min(plotsite_gdd[,"TWAPLS.fx"],na.rm=TRUE),
                     max(plotsite_alpha[,"TWAPLS.fx"],na.rm=TRUE)-min(plotsite_alpha[,"TWAPLS.fx"],na.rm=TRUE))
  #range<-rbind.data.frame(range_WAPLS,range_TWAPLS,range_WAPLS.fx,range_TWAPLS.fx)
  #colnames(range)<-c("range_Tmin","range_gdd","range_alpha");  
  #range$method<-c("WA-PLS","TWA-PLS","WA-PLS with fx","TWA-PLS with fx");
  #range$method<-factor(range$method,levels=c("WA-PLS","TWA-PLS","WA-PLS with fx","TWA-PLS with fx"))
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
  
  pleft<-ggarrange(p_Tmin,p_gdd,p_alpha,ncol=1)
  pright<-ggarrange(p_Tmin_range,p_gdd_range,p_alpha_range,ncol=1)
  
  p<-arrangeGrob(pleft ,pright, ncol = 7,nrow=1,layout_matrix = cbind(1,1,1,1,1,2,2))
  ggsave(file=paste(sitename,"reconstruction results using the last significant number of components.jpeg"),p,width=10,height=9)
  
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