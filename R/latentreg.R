#######################################################################
#' @name latentregpp
#' @docType package
#' @title latentregpp : Item Response Theory Implemented in R and Cpp
#' @description latentregpp is a c++ implementation of the Multidimensional Item Respone Theory (MIRT)
#' cappable of performing parameter and traits estimations. It also provides a list of options to 
#' perform optiman analysis and provides usefull information about the obtained model.
#' @details
#' \tabular{ll}{
#'Package: \tab latentregpp\cr
#'Type: \tab Package\cr
#'Version: \tab 0.0.5\cr
#'Date: \tab 2016-08-22\cr
#'License: \tab MIT + file LICENSE \cr
#'}
#'@author SICS Research Team
#'@keywords IRT MIRT Psychometry 
#'@useDynLib LatentRegpp
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
#'@section Getting Started:
#'Get started with the latentregpp package browsing the index of this documentation
#'if you need help the vignettes should be helpful.
#'@section Getting Started:
#'The latentregpp package allows you to use the latentregpp methodology for simulating, analyzing and scoring tests \cr
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
#'@param quadrature_points Amount of quadrature points
#'@param individual_weights Weights of the quadrature points
#'@param initial_values Initial Values of the estimation
#'@export
latentreg = function(data, dim, model = "2PL", EMepsilon = 1e-4, clusters = NULL,
				  quadrature_technique = NULL, quadrature_points = NULL, 
				  individual_weights = as.integer(c()),
				  initial_values = NULL,
				  verbose = TRUE ) {
	# Asserting matrix type of data
	data = data.matrix(data)

	# Model id
	if ( model == "1PL" ) m = 1
	else if ( model == "2PL" ) m = 2
	else if ( model == "3PL" ) m = 3

	# Quadrature technique
	if ( is.null(quadrature_technique) ) {
		if ( dim < 5 ) quadrature_technique = "Gaussian"
		else quadrature_technique = "QMCEM"
	} else {
		if ( dim >= 5 && quadrature_technique == "Gaussian" ) {
			print("Better use QMCEM")
			# Try to change the qudrature technique
		}
	}

	q = quadpoints(dim = dim, quadrature_technique = quadrature_technique, 
				   quadrature_points = quadrature_points)
	theta = q$theta
	weights = q$weights

	# Type of data
	dichotomous_data = is_data_dicho(data)
	if ( dichotomous_data == -1 ) {
		print("Inconsistent data")
		return -1
	}

	if ( dim == 1 ) {
		# Item parameters estimation
		latentregcpp(Rdata = data, dim = dim, model = m, EMepsilon = EMepsilon, 
					Rtheta = theta, Rweights = weights, 
					Rindividual_weights = individual_weights,
					dichotomous_data = dichotomous_data)
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
				list_initial_values = inivals_MultiUni_NOHARM(data, clusters, model=model, 
				                          find.restrictions=FALSE, verbose=FALSE, probit=FALSE)
			}
		  
			# Item parameters estimation
			latentregcpp(Rdata = data, dim = dim, model = m, EMepsilon = EMepsilon, 
						Rtheta = theta, Rweights = weights, 
						Rindividual_weights = individual_weights,
						dichotomous_data = dichotomous_data,
						Rclusters = clusters,
						Rinitial_values = list_initial_values$coefs )
		} else {
			# TODO find clusters
			# TODO find initial values

			# Item parameters estimation
			#latentregppcpp(Rdata = data, dim = dim, model = m, EMepsilon = EMepsilon, 
			#			Rtheta = theta, Rweights = weights, 
			#			Rindividual_weights = individual_weights,
			#			dichotomous_data = dichotomous_data,
			#			Rclusters = clusters,
			#			Rinitial_values = list_initial_values$coefs )
		}
	}
}