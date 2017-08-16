multigrps <-
function(df,gvar,
    p.rd=3,
    normtest='yes',
    norm.rd=2,
    sk.rd=2,
    tabNA="no",#need to replace NaN with NA for all factors
    cat.rd=2,
    maxfactorlevels=30,
    minfactorlevels=10,
    sim=FALSE,
    workspace=2e5){
##group varibale must be a factor
df[,gvar]<-as.factor(df[,gvar])
#NaN is forced to be NA, NaN can cause problem
df<-replace(df,is.na(df),NA)
for(i in 1:length(levels(df[,gvar]))){
  assign(paste("g", i, sep = ""), levels(df[,gvar])[i])    
}
varlist<-names(df)[!names(df)%in%gvar]      	
table.norm<-data.frame(matrix(rep(999,length(levels(df[,gvar]))*2+3),nrow=1))
nor.names<-c("mean","sd")
for(i in 1:length(levels(df[,gvar]))){
	addmean<-paste("mean", i, sep = ".")
	addsd<-paste("sd", i, sep = ".")
	nor.names<-c(nor.names,addmean,addsd)
}
names(table.norm)<-c(nor.names,"p")
table.skew<-data.frame(matrix(rep(999,length(levels(df[,gvar]))*3+4),nrow=1))
skew.names<-c("median","IQR1","IQR3")
for(i in 1:length(levels(df[,gvar]))){
	addmedian<-paste("median", i, sep = ".")
	addiqr1<-paste("IQR1", i, sep = ".")
	addiqr3<-paste("IQR3", i, sep = ".")
	skew.names<-c(skew.names,addmedian,addiqr1,addiqr3)
}
names(table.skew)<-c(skew.names,"p")
table.cat<-data.frame(matrix(rep(999,length(levels(df[,gvar]))*2+3),nrow=1))
cat.names<-c("No.tot","per.tot")
for(i in 1:length(levels(df[,gvar]))){
	addno<-paste("No", i, sep = ".")
	addper<-paste("per", i, sep = ".")
	cat.names<-c(cat.names,addno,addper)
}
names(table.cat)<-c(cat.names,"p")  
for (var in varlist){
if(class(df[,var])=="factor"&
length(levels(factor(df[,var])))>maxfactorlevels)
     {print(paste("the factor variable",var,
     	"contains more than",
     	maxfactorlevels,"levels",sep=' '))
     next }else{	
 if(class(df[,var])=="factor"|
 length(levels(factor(df[,var])))<=minfactorlevels){
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
	     per.tot=round(as.data.frame(per)[,"Freq"],cat.rd))
	for(i in 1:length(levels(df[,gvar]))){
		assign(paste("No",i,sep="."),as.data.frame.matrix(table.sub)[,get(paste("g",i,sep=''))])
		assign(paste("per",i,sep="."),round(as.data.frame.matrix(per.sub)[,get(paste("g",i,sep=''))],cat.rd))
		frame<-data.frame(frame,
		   get(paste("No",i,sep=".")),
		   get(paste("per",i,sep=".")))
	}
	frame<-data.frame(frame,p=round(p,p.rd))
	colnames(frame)<-names(table.cat)	
	rownames(frame)<-paste(var,levels(df[,var]),sep="_")
	table.cat<-rbind(table.cat,frame)
 }else{
 	if(ad.test(df[,var])$p.value>=0.05|normtest=="no"){
     mean<-round(mean(df[,var],na.rm=T),norm.rd)
	 sd<-round(sd(df[,var],na.rm=T),norm.rd)
	 output<-data.frame(mean,sd)
	 for(i in 1:length(levels(df[,gvar]))){
	 	assign(paste("mean",i,sep="."),round(mean(df[df[,gvar]==get(paste("g",i,sep='')),var],na.rm=T),norm.rd))
		assign(paste("sd",i,sep="."),round(sd(df[df[,gvar]==get(paste("g",i,sep='')),var],na.rm=T),norm.rd))
		output<-data.frame(output,
		   get(paste("mean",i,sep=".")),
		   get(paste("sd",i,sep=".")))
	 }	 
	 p<-summary(aov(df[,var]~df[,gvar]))[[1]][1,"Pr(>F)"]
	 output<-data.frame(output,p=round(p,p.rd))
	 colnames(output)<-names(table.norm) 
	 rownames(output)<-var
	table.norm<-rbind(table.norm,output)
	}else{
		median<-as.numeric(summary(df[,var])[3])
		IQR1<-as.numeric(summary(df[,var])[2])
		IQR3<-as.numeric(summary(df[,var])[5])
		output2<-round(data.frame(median,
		 IQR1,IQR3),sk.rd)
		for(i in 1:length(levels(df[,gvar]))){
			assign(paste("median",i,sep="."),round(as.numeric(summary(df[df[,gvar]==get(paste("g",i,sep='')),var])[3]),sk.rd))
		assign(paste("IQR1",i,sep="."),round(as.numeric(summary(df[df[,gvar]==get(paste("g",i,sep='')),var])[2]),sk.rd))
		assign(paste("IQR3",i,sep="."),round(as.numeric(summary(df[df[,gvar]==get(paste("g",i,sep='')),var])[5]),sk.rd))
		output2<-data.frame(output2,
		   get(paste("median",i,sep=".")),
		   get(paste("IQR1",i,sep=".")),
		   get(paste("IQR3",i,sep=".")))
		}		
		p<-kruskal.test(df[,var]~df[,gvar])$p.value
		output2<-data.frame(output2,p=round(p,p.rd))
		colnames(output2)<-names(table.skew) 
        rownames(output2)<-var
        table.skew<-rbind(table.skew,output2)
	}}
 		}
}
#paste all varibles into one table
tm.names<-"tot"
tot<-paste(table.skew$median,"(",
    table.skew$IQR1,",",
    table.skew$IQR3,")",sep="")
tmsk<-tot
for(i in 1:length(levels(df[,gvar]))){
	assign(paste("grp",i,sep=''),
	paste(table.skew[,paste("median",i,sep=".")],"(",
    table.skew[,paste("IQR1",i,sep=".")],",",
    table.skew[,paste("IQR3",i,sep=".")],")",sep=""))
    tm.names<-c(tm.names,paste("grp",i,sep=''))
    tmsk<-cbind(tmsk,
	get(paste('grp',i,sep="")))
}
tm.names<-c(tm.names,"p")    
tmsk<-cbind(tmsk,p=table.skew$p)
colnames(tmsk)<-tm.names
rownames(tmsk)<-rownames(table.skew)
#normal data
tot<-paste(table.norm$mean,"\U00B1",
    table.norm$sd,sep="")
tmnorm<-tot    
for(i in 1:length(levels(df[,gvar]))){
	assign(paste("grp",i,sep=''),
	paste(table.norm[,paste("mean",i,sep=".")],
	 "\U00B1",
     table.norm[,paste("sd",i,sep=".")],sep=""))
    tmnorm<-cbind(tmnorm,
	get(paste('grp',i,sep="")))
}
tmnorm<-cbind(tmnorm,p=table.norm$p)
colnames(tmnorm)<-tm.names
rownames(tmnorm)<-rownames(table.norm)
#categorical data
tot<-paste(table.cat$No.tot,"(",
    table.cat$per.tot,")",sep="")
tmcat<-tot    
for(i in 1:length(levels(df[,gvar]))){
	assign(paste("grp",i,sep=''),
	paste(table.cat[,paste("No",i,sep=".")],
	 "(",
     table.cat[,paste("per",i,sep=".")],")",sep=""))
    tmcat<-cbind(tmcat,
	get(paste('grp',i,sep="")))
}
tmcat<-cbind(tmcat,p=table.cat$p)
colnames(tmcat)<-tm.names
rownames(tmcat)<-rownames(table.cat)
#combine all variable types  
table<-rbind(tmnorm,tmsk,tmcat)
colnames(table)<-NULL
table<-rbind(tm.names,table)
table<-table[!table[,length(levels(df[,gvar]))+2]=="999",]
table.norm<-table.norm[-1,]
table.skew<-table.skew[-1,]
table.cat<-table.cat[-1,]
results<-list(table=table,
  table.cat=table.cat,
  table.norm=table.norm,
  table.skew=table.skew)
for(i in 1:length(levels(df[,gvar]))){
  	results<-c(results,
  	setNames(as.list(levels(df[,gvar])[i]),paste("g", i, sep = "")))
  }
return(results)       	      	     	
   }
