

########################################################################################################################################
# correlation between different variable types

mixed_assoc = function(df, cor_method="spearman", adjust_cramersv_bias=TRUE){
  
  df_comb = data.frame(t(combn(sort(names(df)), 2)), stringsAsFactors = F) %>% set_names("X1", "X2")
  
  is_nominal = function(x) class(x) %in% c("factor", "character")
  is_numeric <- function(x) { is.integer(x) || is_double(x)}
  
  f = function(xName,yName) {
    x =  pull(df, xName)
    y =  pull(df, yName)
    
    result = if(is_nominal(x) && is_nominal(y)){
      # use bias corrected cramersV as described in https://rdrr.io/cran/rcompanion/man/cramerV.html
      cv = cramerV(as.character(x), as.character(y), bias.correct = adjust_cramersv_bias)
      data.frame(xName, yName, assoc=cv, type="cramersV")
      
    }else if(is_numeric(x) && is_numeric(y)){
      correlation = cor(x, y, method=cor_method, use="complete.obs")
      data.frame(xName, yName, assoc=correlation, type="correlation")
      
    }else if(is_numeric(x) && is_nominal(y)){
      r_squared = summary(lm(x ~ y))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else if(is_nominal(x) && is_numeric(y)){
      r_squared = summary(lm(y ~x))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else {
      warning(paste("unmatched column type combination: ", class(x), class(y)))
    }
    
    # finally add complete obs number and ratio to table
    result %>% mutate(complete_obs_pairs=sum(!is.na(x) & !is.na(y)), complete_obs_ratio=complete_obs_pairs/length(x)) %>% rename(x=xName, y=yName)
  }
  
  # apply function to each variable combination
  map2_df(df_comb$X1, df_comb$X2, f)
}

########################################################################################################################################

#This is a function calculating the pure effects plus overlaps in lm()-type models, anovas and dbRDAs via variation partitioning

#Y=Target variable OR distance matrix (dbRDA variance partitioning) OR lm output
#p=predictor variables (character vector, size n, NOT needed if Y is lm output) 
#data=Data.frame containing predictor variables (not needed if Y is lm output)
#scale (logical): should the data be z-transformed?
#nested: is any factor nested within another? write as "outer:inner" (not needed if Y is lm output)
#minus.null (logical): Convert negative variance partitions to zero?

VarPart<-function(Y,p,data,adjust=TRUE,scale=FALSE,nested=NULL,minus.null=FALSE){
  #should the data be z-transformed?
  if(scale==TRUE){
    data<-data.frame(scale(data[,c(Y,p)]))
  }
  if(missing(data)){
    ifelse(class(Y)=="lm",data<-Y[["model"]],warning("source data not defined"))
  }
  if(missing(p)){
    ifelse(class(Y)=="lm",{
      p<-colnames(data)[-1]
      if(grepl(":",Y[["call"]][2])){
        inputs<-strsplit(as.character(Y[["call"]][2])," ")[[1]]
        nested<-inputs[grepl(":",inputs)]
        response<-inputs[1]
      }
    },warning("predictor variables not defined"))
  }
  
  # combination of all explanatory variables
  combo<-lapply(1:length(p),function(m){combn(p,m)})
  
  #calculation of Rsquared and adjusted Rsquared for any combination of explanatory variables
  r2.res<-sapply(combo,function(a){
    apply(a,2,function(b){
      #If there are nested effects they will be added here
      if(!is.null(nested)){
        nest_in<-strsplit(nested,":")[[1]][2]
        if(nest_in%in%b){
          b<-c(b[b!=nest_in],nested)
        }}
      
      #formula for the individual models
      form<-as.formula(paste(ifelse(
        class(Y)=="lm",response,ifelse(
          is.dist(Y),"Y",Y)),"~",paste(b,collapse="+")))
      
      #type of response data (dbRDA vs linear model)
      ifelse(is.dist(Y),
             {mod.res<-dbrda(form,data,na.action = na.omit)
             },{mod.res<-lm(form,data,na.action = na.omit)
             })
      unlist(c(paste(b,collapse="+"),RsquareAdj(mod.res)))
    })})
  r2.res<-do.call(cbind,r2.res)
  r2.res<-data.frame(predictor=r2.res[1,],expl.Var=as.numeric(r2.res[2,]),adjust=as.numeric(r2.res[3,]))
  ifelse(adjust==TRUE,pre<-r2.res$adjust,pre<-r2.res$expl.Var)
  r2.ttl<-pre[nrow(r2.res)]
  prd.all<-strsplit(as.character(r2.res$predictor),"\\+")
  prd.n<-prd.all[[nrow(r2.res)]]
  
  expl<-vector(mode="numeric",length=nrow(r2.res))
  nexpl<-vector(mode="numeric",length=nrow(r2.res))
  varP<-vector(mode="numeric",length=nrow(r2.res))
  
  #Calculate the pure partitions and overlaps
  for(o in r2.res$predictor){
    prd.o<-strsplit(as.character(o),"\\+")[[1]]
    combo.o<-unlist(sapply(1:length(prd.o),function(m){
      apply(combn(prd.o,m),2,paste,collapse = "+")
    }))
    ttl.p<-prd.n[!prd.n%in%prd.o]
    expl[r2.res$predictor==o]<-sum(varP[r2.res$predictor%in%combo.o&r2.res$predictor!=as.character(o)])
    nexpl[r2.res$predictor==o]<-sum(pre[r2.res$predictor==paste(ttl.p,collapse = "+")])
    varP<-r2.ttl-expl-nexpl
  }
  
  #If there are negative partitions show them as zeroes?
  #Legendre (2008) [doi: 10.1093/jpe/rtm001] "Negative values of R2 are interpreted as zeros; they correspond to cases where the explanatory variables explain less variation than random normal variables would."
  if(minus.null==TRUE){varP[varP<0]<-0}
  r2.res$varP<-varP
  
  #with negative partitions, the total sum of partitions might not equal the total explained variance of the full model
  if(sum(r2.res[,3][r2.res[,3]>0])!=r2.ttl){warning("total explained variance != sum of partions. check for negative variances")}
  return(r2.res)
}
