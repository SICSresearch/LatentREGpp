

#'@name inivals_MultiPoly
#'@title Initial values for multidimensional polytomous data 
#'@description Get initial values for polytomous data according dataset categories and
#'size cluster for dimensionality. It has no 3PL model. Then, no initial values for c parameter.
#'@param data_poly Polytomous data to get initial values.
#'@param size.cluster A vector with dimensionality test.
#'@param verbose True for get information about process in runtime. False in otherwise. False by default.
#'@export
inivals_MultiPoly<- function(data_poly, size.cluster, verbose=F)
{

#
#	SEPARATE CLUSTERS (each cluster mean one dim)
#
if(verbose){print("Starting Process to Compute Multidimensional Polytomous Initial Values ")}
ncat<- apply(data_poly,2, max)

data_f<- matrix(NA, ncol= ncol(data_poly), nrow=nrow(data_poly))
for(i in 1:ncol(data_poly))	
{
	data_f[,i]<- as.factor(data_poly[,i])
}
 #separate cluster
  cluster<- list()
  start = 1
  for(i in 1:length(size.cluster)){
    cluster[[i]]<- data_poly[,start:(start + (size.cluster[i] - 1 ))]
    start = start + size.cluster[i] 
  }
	
#
#	ESTIMATE TRAITS VIA PCA
#

traits<- matrix(NA, ncol= length(size.cluster), nrow=nrow(data_poly))
for(i in 1:length(size.cluster))
{
ACP<- FactoMineR::PCA(cluster[[i]], graph=F)
traits[,i]<-ACP$ind$coor[,1]
}

C<- cov(traits)
E<- eigen(C); V<- E$vectors; D<- diag(E$values)
Sqrt<- V%*%sqrt(D)%*%solve(V)
std.traits<- t(solve(Sqrt)%*%t(traits))

traits<- std.traits



#
#	ESTIMATE ITEM PARAMETER VIA Ordinal Logistic Regression
#


MODEL<- list()
COEF<- list()
for(i in 1:ncol(data_poly))
{
## fit ordered logit model and store results 'm'
MODEL[[i]] <- polr(as.factor(data_poly[,i]) ~ traits,
	method="logistic",Hess=TRUE, model=T)

d<- -MODEL[[i]]$zeta
a<- MODEL[[i]]$coef
M<- matrix(c(a,d), nrow=1)
colnames(M)<- c(paste("a",seq(1,length(a)), sep="") ,paste("d",seq(1,length(d)), sep="") )
rownames(M)<- paste("Item", i)

COEF[[i]]<- M
}


coef<- matrix(NA, nrow=length(COEF), ncol= ((max(ncat)-1)+length(size.cluster)))

for(i in 1:ncol(data_poly))
{
coef[i,1:length(COEF[[i]])] <- COEF[[i]]
}
COEF[[1]]
id_aux<- unlist(lapply(COEF, ncol))
id_max_aux<-which(id_aux == max(id_aux))[1]
colnames(coef)<- colnames(COEF[[id_max_aux]])


#
#	Model$zeta: Is the threshold parameter of the polytomous model
#	Model$coef: Is the discrimination parameter
#


if(verbose){print("Done")}

res<- coef
return(res)

}

