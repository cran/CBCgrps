twogrps <-
function(df,gvar,varlist = NULL,
    p.rd=3,
    skewvar=NULL,
    norm.rd=2,
    sk.rd=2,
    tabNA="no",#need to replace NaN with NA for all factors
    cat.rd=0, pnormtest=0.05,
    maxfactorlevels=30,
    minfactorlevels=10,
    sim = FALSE,#to use simulated p value
    workspace=2e5){
##group varibale must be a factor
df[,gvar]<-as.factor(df[,gvar])
#NaN is forced to be NA, NaN can cause problem
df<-replace(df,is.na(df),NA)
g1<-levels(df[,gvar])[1]
g2<-levels(df[,gvar])[2]
if(is.null(varlist)){
  varlist<-names(df)[!names(df)%in%gvar]
}else if(sum(!(varlist%in%names(df)))>0){
  stop("varlist contains variables not in the data frame")
} else {
  varlist = varlist
}
if(sum(!(skewvar%in%names(df)))>0){
  stop("skewvar contains variables not in the data frame")
}


  
table.norm<-data.frame(matrix(rep(999,7),nrow=1))
names(table.norm)<-c("mean","sd","mean.1",
    "sd.1","mean.2","sd.2","p")
table.skew<-data.frame(matrix(rep(999,10),nrow=1))
names(table.skew)<-c("median","IQR1","IQR3",
    "median.1","IQR1.1","IQR3.1",
    "median.2","IQR1.2","IQR3.2",
    "p")
table.cat<-data.frame(matrix(rep(999,7),nrow=1))
names(table.cat)<-c("No.tot","per.tot","No.1",
    "per.1","No.2","per.2","p")  
for (var in varlist){
if(class(df[,var])=="factor"&
length(levels(factor(df[,var])))>maxfactorlevels)
     {print(paste("the factor variable",var,
     	"contains more than",
     	maxfactorlevels,"levels",sep=' '))
     next }else{	
 if(class(df[,var])=="factor"|
 length(levels(factor(df[,var])))<=minfactorlevels){
   if(var%in%skewvar){
     stop("skewvar contains categorical variables")
   }
   if(tabNA=="no"){
		df[,var]<-factor(df[,var])
		}else{
		  df[,var]<-factor(df[,var],exclude = NULL)
	    }
   table<-table(df[,var],useNA=tabNA)
	 per<-prop.table(table)
	 table.sub<-table(df[,var],df[,gvar],useNA=tabNA)
	 per.sub<-prop.table(table.sub,2)
	 p<-tryCatch({#using fisher's test when scarce data
          chisq.test(table.sub)$p.value
       }, warning = function(w) {
          fisher.test(table.sub,
          workspace = workspace,
          simulate.p.value = sim)$p.value
       })
       	frame<-data.frame(No.tot=as.data.frame(table)[,"Freq"],
	     per.tot=round(as.data.frame(per)[,"Freq"]*100,cat.rd),
	     No.1=as.data.frame.matrix(table.sub)[,g1],
	     per.1=round(as.data.frame.matrix(per.sub)[,g1]*100,cat.rd),
	     No.2=as.data.frame.matrix(table.sub)[,g2],
	     per.2=round(as.data.frame.matrix(per.sub)[,g2]*100,cat.rd),
	     p=round(p,p.rd))
	rownames(frame)<-paste(var,levels(df[,var]),sep="_.._")
	table.cat<-rbind(table.cat,frame)
 }else{
 	if((ad.test(df[,var])$p.value>=pnormtest&is.null(skewvar))|(!(var%in%skewvar)&!is.null(skewvar))){
     mean<-round(mean(df[,var],na.rm=T),norm.rd)
	 sd<-round(sd(df[,var],na.rm=T),norm.rd)
	 mean.1<-round(mean(df[df[,gvar]==g1,var],na.rm=T),norm.rd)
	 sd.1<-round(sd(df[df[,gvar]==g1,var],na.rm=T),norm.rd)
	 mean.2<-round(mean(df[df[,gvar]==g2,var],na.rm=T),norm.rd)
	 sd.2<-round(sd(df[df[,gvar]==g2,var],na.rm=T),norm.rd)
	 p<-round(t.test(df[,var]~df[,gvar])$p.value,p.rd)
	 output1<-data.frame(mean,sd,mean.1,
	 sd.1,mean.2,sd.2,p)
	 rownames(output1)<-var
	table.norm<-rbind(table.norm,output1)
	}else{
		median<-as.numeric(summary(df[,var])[3])
		IQR1<-as.numeric(summary(df[,var])[2])
		IQR3<-as.numeric(summary(df[,var])[5])
		median.1<-as.numeric(summary(df[df[,gvar]==g1,var])[3])
		IQR1.1<-as.numeric(summary(df[df[,gvar]==g1,var])[2])
		IQR3.1<-as.numeric(summary(df[df[,gvar]==g1,var])[5])
		median.2<-as.numeric(summary(df[df[,gvar]==g2,var])[3])
		IQR1.2<-as.numeric(summary(df[df[,gvar]==g2,var])[2])
		IQR3.2<-as.numeric(summary(df[df[,gvar]==g2,var])[5])
		p<-wilcox.test(df[,var]~df[,gvar])$p.value
		output2<-cbind(round(data.frame(median,
		 IQR1,IQR3,
         median.1,IQR1.1,IQR3.1,
         median.2,IQR1.2,IQR3.2),sk.rd),
         p=round(p,p.rd))
        rownames(output2)<-var
        table.skew<-rbind(table.skew,output2)
	}}
 		}
}
#paste all varibles into one table
##skewed data
tot<-paste(table.skew$median,"(",
    table.skew$IQR1,",",
    table.skew$IQR3,")",sep="")
grp1<-paste(table.skew$median.1,"(",
    table.skew$IQR1.1,",",
    table.skew$IQR3.1,")",sep="")
grp2<-paste(table.skew$median.2,"(",
    table.skew$IQR1.2,",",
    table.skew$IQR3.2,")",sep="")
tmsk<-cbind(tot,grp1,grp2,
            p=ifelse(table.skew$p<1*10^(-p.rd),
                     paste("<",1*10^(-p.rd),sep = ""),
                     table.skew$p))
rownames(tmsk)<-rownames(table.skew)
##norm data
tot<-paste(table.norm$mean,"\U00B1",
    table.norm$sd,sep="")
grp1<-paste(table.norm$mean.1,"\U00B1",
    table.norm$sd.1,sep="")
grp2<-paste(table.norm$mean.2,"\U00B1",
    table.norm$sd.2,sep="")
tmnorm<-cbind(tot,grp1,grp2,
              p=ifelse(table.norm$p<1*10^(-p.rd),
                       paste("<",1*10^(-p.rd),sep = ""),
                       table.norm$p))
rownames(tmnorm)<-rownames(table.norm)
##categorical data
tot<-paste(table.cat$No.tot,"(",
    table.cat$per.tot,")",sep="")
grp1<-paste(table.cat$No.1,"(",
    table.cat$per.1,")",sep="")
grp2<-paste(table.cat$No.2,"(",
    table.cat$per.2,")",sep="")
tmcat<-cbind(tot,grp1,grp2,
             p=ifelse(table.cat$p<1*10^(-p.rd),
                      paste("<",1*10^(-p.rd),sep = ""),
                      table.cat$p))
rownames(tmcat)<-rownames(table.cat) 
table<-rbind(tmnorm,tmsk,tmcat)
table<-table[!table[,4]=="999",]
tabdup <- as.data.frame(apply(table,MARGIN = 2, as.character),
                        stringsAsFactors = F)
rownames(tabdup) <- rownames(table)
tabdup[,"var"] <- gsub("_.._.*","",rownames(tabdup))
tabAddhead <- NULL;
for (var in varlist) {
  tabsub <- tabdup[tabdup$var==var,]
  if(nrow(tabsub)>1){
    newline <- data.frame(matrix(c(rep("",3),tabsub$p[1],tabsub$var[1]),nrow = 1),
                          stringsAsFactors = F)
    names(newline) <- names(tabsub)
    rownames(newline) <- tabsub$var[1]
    tab1line <- rbind.data.frame(newline,tabsub)
    tab1line[2:(nrow(tabsub)+1),"p"] <- ""
    rownames(tab1line)[2:(nrow(tabsub)+1)] <- sub(paste(var,"_.._",sep = ""),
                                                   rownames(tab1line)[2:(nrow(tabsub)+1)],
                                                  replacement = "    ")
    tabAddhead <- rbind(tabAddhead,tab1line)
  } else{
    tabAddhead <- rbind(tabAddhead,tabsub)
  }
}
tabAddhead <- tabAddhead[,-5]
tabAddhead<-rbind(c(paste("Total (n=",nrow(df),")",sep = ""),
                    paste(g1,"(n=",nrow(df[df[,gvar]==g1,]),")",sep = ""),
                    paste(g2,"(n=",nrow(df[df[,gvar]==g2,]),")",sep = ""),"p"),
                  tabAddhead)
rownames(tabAddhead)[1] <- "Variables"
colnames(tabAddhead)<-NULL
##combine them together
table.norm<-table.norm[-1,]
table.skew<-table.skew[-1,]
table.cat<-table.cat[-1,]
results<-list(table=tabAddhead,
  g1=g1,g2=g2,
  table.cat=table.cat,
  table.norm=table.norm,
  table.skew=table.skew)
return(results)       	      	
#the end of the function      	
   }
