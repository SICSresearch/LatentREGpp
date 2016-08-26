#######################################################################
#' @name LatentREGpp
#' @docType package
#' @title LatentREGpp : Item Response Theory Implemented in R and Cpp
#' @description latentregpp is a c++ implementation of the Multidimensional Item Respone Theory (MIRT)
#' cappable of performing parameter and traits estimations. It also provides a list of options to 
#' perform optiman analysis and provides usefull information about the obtained model.
#' @details
#' \tabular{ll}{
#'Package: \tab LatentREGpp\cr
#'Type: \tab Package\cr
#'Version: \tab 0.0.5\cr
#'Date: \tab 2016-08-22\cr
#'License: \tab MIT + file LICENSE \cr
#'}
#'@author SICS Research Team
#'@keywords IRT MIRT Psychometry 
#'@useDynLib LatentREGpp
#'@importFrom Rcpp sourceCpp
#'@importFrom randtoolbox sobol
#'@importFrom fastGHQuad gaussHermiteData
#'@importFrom RSpectra eigs
#'@importFrom sirt noharm.sirt
#'@importFrom IRTpp irtpp
#'@importFrom IRTpp parameter.matrix
#'@importFrom IRTpp individual.traits
#'@importFrom FactoMineR PCA
#'@importFrom FactoMineR HCPC
#'@importFrom MASS polr
#'@section Getting Started:
#'Get started with the LatentREGpp package browsing the index of this documentation
#'if you need help the vignettes should be helpful.
#'@section Getting Started:
#'The LatentREGpp package allows you to use the LatentREGpp methodology for simulating, analyzing and scoring tests \cr
#'You can browse the package vignettes to get started.
#'
NULL

#'@name latentreg
#'@title Parameter estimation of a test
#'@description Estimates the test parameters according to the Multidimensional Item Response Theory
#'@param data The matrix containing the answers of tested individuals
#'@param dim The Dimensionality of the test
#'@param model 1PL, 2PL or 3PL
#'@param EMepsilon Convergence value to determine the accuracy of the test
#'@param clusters Clusters per dimention
#'@param quadratura_technique Quasi-Monte Carlo or Gaussian
#'@param quad_points Amount of quadrature points
#'@param individual_weights Weights of the quadrature points
#'@param initial_values Initial Values of the estimation
#'@export
latentreg = function(data, dim, model = "2PL", EMepsilon = 1e-4, clusters = NULL,
				  quad_tech = NULL, quad_points = NULL, 
				  individual_weights = as.integer(c()),
				  initial_values = NULL,
				  verbose = TRUE, save_time = TRUE ) {
	# Quadrature technique
	if ( is.null(quad_tech) ) {
		if ( dim < 5 ) quad_tech = "Gaussian"
		else quad_tech = "QMCEM"
	} else {
		if ( dim >= 5 && quad_tech == "Gaussian" ) {
			message("For dim >= 5 QMCEM quadrature technique is recommended")
			input = readline(prompt = "Are you sure you want continue with Gaussian quadrature? [y/n]: ")
			if ( input == "n" || input == "N" )
				quad_tech = "QMCEM"
		}
	}

	if ( save_time ) first_time = Sys.time()

	# Asserting matrix type of data
	data = data.matrix(data)

	# Model id
	if ( model == "1PL" ) m = 1
	else if ( model == "2PL" ) m = 2
	else if ( model == "3PL" ) m = 3

	q = quadpoints(dim = dim, quad_tech = quad_tech, 
				   quad_points = quad_points)
	theta = q$theta
	weights = q$weights

	# Type of data
	dichotomous_data = is_data_dicho(data)
	if ( dichotomous_data == -1 )
		stop("Inconsistent data")

	if ( dim == 1 ) {
		# Item parameters estimation
		obj_return = (latentregcpp(Rdata = data, dim = dim, model = m, EMepsilon = EMepsilon, 
							Rtheta = theta, Rweights = weights, 
							Rindividual_weights = individual_weights,
							dichotomous_data = dichotomous_data))
	} else {
		if ( dichotomous_data ) {

			if ( is.null(clusters) ) {
				#1. Find a temporal cluster
				p = ncol(data)
				d = dim
				clusters = find_temporal_cluster(p=p,d=d)

				#2. Find initial values
				III = inivals_MultiUni_NOHARM(data, clusters, model=model, 
				                            find.restrictions=TRUE, verbose=FALSE, probit=FALSE)

				#3. do something weird :v 
				# Normalize Discrimination vectors per item
				#A_matrix<- III$coefs[,1:d]
				#Sigma<- cov(A_matrix)
				#CholSigm<- chol(Sigma)
				#A_asterisco<- A_matrix%*%solve(CholSigm)
				#Betas<- normalize(t(A_asterisco))
				#cov(t(Betas))

				#4. Clustering and find clustering
				#CLUST<- find_cluster(data=t(A_matrix), ang=22.5, h7=0.9, q_proy= 0.6)
				#unlist(lapply(CLUST,ncol))
				#paste("N clust",length(CLUST))
				#acpc<-PCA(CLUST[[1]])


				#BONUS. Another Reclassify
				acp = FactoMineR::PCA(X = III$coefs,graph = FALSE)
				hcpc = FactoMineR::HCPC(acp,nb.clust = d,graph = FALSE)
				CLUST_final<- list()
				for (i in 1:length(table(hcpc$data.clust$clust))) {
					CLUST_final[[i]]<- data[,which(hcpc$data.clust$clust==i)]
				}

				#5. Find Reclassify
				#P_axes_filter<- Principal_axes(CLUST)
				#CLUSTB<- reclass_data_Praxes(data=t(A_matrix), Principal_axes=t(P_axes_filter))
				#unlist(lapply(CLUSTB, ncol)); paste("N clust",length(CLUSTB));

				#BONUS. Printing the clusters
				#print("New Clusters are: ")
				#print(CLUST_final)

				clusters = unlist(lapply(CLUST_final, ncol))
			}

			#Initial values
			if ( is.null(initial_values) ) {
				#Find initial values again
				initial_values = inivals_MultiUni_NOHARM(data, clusters, model=model, 
				                          find.restrictions=FALSE, verbose=FALSE, probit=FALSE)$coefs
			} else 
				initial_values = data.matrix(initial_values)
		  
			# Item parameters estimation
			obj_return = (latentregcpp(Rdata = data, dim = dim, model = m, EMepsilon = EMepsilon, 
								Rtheta = theta, Rweights = weights, 
								Rindividual_weights = individual_weights,
								dichotomous_data = dichotomous_data,
								Rclusters = clusters,
								Rinitial_values = initial_values ))
		} else {
			# TODO find clusters
			if ( is.null(clusters) )
				stop("You must specify clusters")

			#Initial values
			if ( is.null(initial_values) ) {
				#Find initial values again
				initial_values = inivals_MultiPoly(data_poly = data, size.cluster = clusters, verbose=F)
			} else 
				initial_values = data.matrix(initial_values)

			# Item parameters estimation
			obj_return = latentregcpp(Rdata = data, dim = dim, model = m, EMepsilon = EMepsilon, 
									Rtheta = theta, Rweights = weights, 
									Rindividual_weights = individual_weights,
									dichotomous_data = dichotomous_data,
									Rclusters = clusters,
									Rinitial_values = initial_values )
		}
	}

	if ( save_time ) second_time = Sys.time()
	obj_return$time = second_time - first_time	

	return (obj_return)
}