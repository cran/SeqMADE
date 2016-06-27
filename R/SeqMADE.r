# copyright: Mingli Lei<leimingli2013@sjtu.edu.cn>
# .packageName <- "SeqMADE"

#############################################################################################
##
##'exprs' is a data frame or matrix for two groups or conditions, with rows as variables (genes) and columns as samples
##gA is the sample list in group A, for example gA = c("A1","A2","A3","A4","A5")
##gB is the sample list in group B, for example gB = c("B1","B2","B3","B4","B5")
###output:
##Count (factor1): A numeric expression value for all genes in group A and group B
##Group (factor2) indicator factor for the genes when Group=1 represent that the expression value is generated from group A, otherwise gB
##Direction (factor3) indicator factor for the genes when Direction=1 represent that gene is up-regulated, otherwise down-regulated
##
#############################################################################################################

Factor <- function(exprs,case,control) {

  #Direction factor
   Direction <- matrix(nrow=nrow(exprs),ncol=1)
	for (i in 1:nrow(exprs)) {
		groupA <- sum(as.numeric(exprs[i,case]))
		groupB <- sum(as.numeric(exprs[i,control]))
		if(groupA > groupB) {
			Direction[i]=1
		}
		if(groupA <= groupB) {
			Direction[i]=-1
		}
	}
	Direction <- Direction[,rep(1,length(case)+length(control))]
	Direction <- as.factor(Direction) 
	
  #Count value
   exprs <- exprs[,c(case,control)] 
   exprs <- as.matrix(exprs)
   Count <- as.numeric(exprs)
  
  #Group factor
   Group <- matrix(rep(1:0,c(length(case)*nrow(exprs),length(control)*nrow(exprs))),nrow=nrow(exprs),ncol=length(case)+length(control))
   Group <- as.factor(Group)
  
  #result
   factors <- list(Count=Count,Group=Group,Direction=Direction)
   return(factors)
}


#########################################################################################################
##'exprs' is a data frame or matrix for two groups or conditions, with rows as variables (genes) and columns as samples
##'networkModule' is the gene sets or modules in the biological network or metabolic pathway, with the 1th column as the module names and the 2th columnn as the gene symbol constituting the module
##output: modulematrix is a matrix, in which the indicator variables 1 or 0 represent whether a gene belong to a given module
#########################################################################################################
"moduleMatrix" <- function(exprs,networkModule){
	modulematrix <- matrix(nrow=nrow(exprs), ncol=nrow(networkModule)+1)
	for (i in 1:nrow(exprs)) {
		gene <- as.character(exprs[i,1])
    
		for (j in 1:nrow(networkModule)) {
			genelist <- as.character(networkModule[j,2])
			genes <- strsplit(genelist, ";")
			sum <- 0
      
			for (k in 1:length(genes[[1]])) {
				if (genes[[1]][k] == gene)
				{
					sum <- sum + 1 
				}
				if (sum == 1) 
				{
					next
				}
			}	
      
			modulematrix[i,j] <- sum
      
		}
    
	}
  
	for (i in 1:nrow(exprs)) {
		if (sum(modulematrix[i,1:nrow(networkModule)]) == 0) {
			modulematrix[i,nrow(networkModule)+1] <- 1
		} else {
			modulematrix[i,nrow(networkModule)+1] <- 0
		}  
	}
	return(modulematrix)
}


############################################################################################################################
##N represents the sample size.
############################################################################################################################
"nbGLM" <- function(factors ,N ,networkModule, modulematrix, distribution = c("poisson", "NB")[1]) {

    Count <- factors$Count
	Group <- factors$Group
	p.nominal <- matrix(nrow=ncol(modulematrix)-1,ncol=1)
	for (i in 1:(ncol(modulematrix)-1)) {
		Module <- as.matrix(modulematrix[,i])
		Module <- Module[,rep(1,N)]
		Module <- as.factor(Module)  
		if (distribution=='NB'){
		fm <- glm.nb(Count ~ Group * Module )
		}
		if (distribution=='poisson'){
		fm <- glm(Count ~ Group * Module ,family=poisson())
		}
		p.nominal[i] <- summary(fm)$coefficients[4,4]
	}
	fdr <- p.adjust(p.nominal,method="fdr")
  	Result <- data.frame(p.value=p.nominal,fdr=fdr,Gene.Symbol=networkModule[,2]);
  	row.names(Result) <-networkModule[,1];
	colnames(Result) <- c( "p.nominal", "fdr", "Gene.Symbol")
	return(Result)
} 

############################################################################################################################
##nbGLMdir is the generalized linear model with three indicator factors Group, Module and Direction
############################################################################################################################
#'nbGLM is the generalized linear model with two indicator factors Group and Module
"nbGLMdir" <- function(factors,N,networkModule, modulematrix, distribution = c("poisson", "NB")[1]){
	Count <- factors$Count
	Group <- factors$Group
	Direction <- factors$Direction
	p.nominal <- matrix(nrow=ncol(modulematrix)-1,ncol=1)
for (i in 1:(ncol(modulematrix)-1)) {
		Module <- as.matrix(modulematrix[,i])
		Module <- Module[,rep(1,N)]
		Module <- as.factor(Module)  
		if (distribution=='NB'){
			fm <- glm.nb(Count ~ Group * Module * Direction)
			if (length(summary(fm)$coefficients)==32){
				p.nominal[i] <- summary(fm)$coefficients[8,4]
			}else{
				fm1 <- glm.nb(Count ~ Group * Module)
				p.nominal[i] <- summary(fm1)$coefficients[4,4]
		    }
		}
		if (distribution=='poisson'){
			fm <- glm(Count ~ Group * Module * Direction,family=poisson() )
			if (length(summary(fm)$coefficients)==32){
				p.nominal[i] <- summary(fm)$coefficients[8,4]
			}else{
				next
		    }
		}
		
		
	}
	
		
	fdr <- p.adjust(p.nominal,method="fdr")
  	Result <- data.frame(p.value=p.nominal,fdr=fdr,Gene.Symbol=networkModule[,2]);
  	row.names(Result) <-networkModule[,1];
	colnames(Result) <- c( "p.nominal", "fdr", "Gene.Symbol")
	return(Result)
} 




 
########################################################################################################################
#'nbGLMdir is the generalized linear model with three indicator factors Group, Module and Direction using permutation
## Disturb the sample labels for the two conditions and re-assign the samples to two datasets,then calculate the 'z-score' for 
## N times and then pool all the z-score together to construct a 'NULL' distribution.
#########################################################################################################################
"nbGLMdirperm" <- function(exprs,case,control,factors,networkModule,modulematrix,N,distribution = c("poisson", "NB")[1]) {
	Count <- factors$Count
	Group <- factors$Group
	Direction <- factors$Direction
	exprsvalue <- exprs[,c(case,control)]
	z.score <- matrix(nrow=ncol(modulematrix)-1,ncol=1)
	for (i in 1:(ncol(modulematrix)-1)) {
		Module <- as.matrix(modulematrix[,i])
		Module <- Module[,rep(1,ncol(exprsvalue))]
		Module <- as.factor(Module)  
		if (distribution=='NB'){
		fm <- glm.nb(Count ~ Group * Module * Direction)
		}
		if (distribution=='poisson'){
		fm <- glm(Count ~ Group * Module * Direction, family=poisson() )
		}
		z.score[i] <- summary(fm)$coefficients[8,3]
	}
	if(N>0) {
		ZS <- matrix(nrow=ncol(modulematrix)-1,ncol=N)
		n.1 = length(case)
		n.2 = length(control)
		cat.j = 0
		for(j in 1:N) {
			if ( (j*100/N)%/%10>cat.j) {
				cat.j = cat.j+1
				cat(cat.j*10,'%','\n')
			}
			seq <- sample(n.1+n.2)
			exprs.1 <- exprs[,seq[1:n.1]];
			exprs.2 <- exprs[,seq[(n.1+1):(n.1+n.2)]];
			fm <- glm.nb(Count ~ Group * Module * Direction)
  		ZS[,j] = summary(fm)$coefficients[8,3]
		}

		p.value <- sapply(z.score,function(x) sum(as.vector(ZS)>x,na.rm=T)/length(!is.na(as.vector(ZS))))
  		fdr <- p.adjust(p.value,method="fdr")
  		Result <- data.frame(z.score=z.score,p.value=p.value,fdr=fdr,Gene.Symbol=networkModule[,2]);
  		row.names(Result) <-networkModule[,1];
		colnames(Result) <- c("z.score", "pvalue.CRC", "pvalue.fdr", "Gene.Symbol")
	}
	else { 
  		Result<- data.frame(z.score=z.score,p.value=p.value,fdr=fdr,networkModule[,2]);
  		row.names(Result) <- networkModule[,1];
		colnames(Result) <- c("z.score", "pvalue.CRC", "pvalue.fdr", "Gene.Symbol")
  	} 
  	return(Result)
	
}