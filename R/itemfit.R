#######################################################################
#' @name LatentREGpp
#' @docType package
#' @title LatentREGpp : Item Response Theory Implemented in R and Cpp
#' @description This Package is a c++ implementation of the Multidimensional Item Response Theory (MIRT) capable of performing parameter and traits estimations. It also provides a list of options to perform an optimal analysis and provides useful information about the obtained model.\cr
#' This package is a work of SICS Research Group, Universidad Nacional de Colombia.\cr
#' @details
#' \tabular{ll}{
#'Package: \tab LatentREGpp\cr
#'Type: \tab Package\cr
#'Version: \tab 1.0\cr
#'Date: \tab 2016-09-07\cr
#'License: \tab MIT + file LICENSE \cr
#'}
#'@author Milder Hernandez Cagua <milhernandezcag@unal.edu.co>
#'@author Jhonatan Javier Guzman del Rio <jhjguzmanri@unal.edu.co>
#'@author Camilo Antonio Suarez Bolanos <caasuarezbo@unal.edu.co>
#'@author Alvaro Mauricio Montenegro <ammontenegrod@unal.edu.co>
#'@author Sergio Paez Moncaleano <spaezm@unal.edu.co>
#'@author Juan David Cortes Castillo <jdcortesc@unal.edu.co>
#'@author Campo Elias Pardo Turriago <cepardot@unal.edu.co>
#'@author Luisa Fernanda Parra Arboleda <lfparraar@unal.edu.co>
#'@author Emilio Pablo Berdugo Camacho <epberdugoc@unal.edu.co>
#'@keywords IRT MIRT Psychometry 
#'@useDynLib LatentREGpp
#'@importFrom Rcpp sourceCpp
#'@importFrom randtoolbox sobol
#'@importFrom fastGHQuad gaussHermiteData
#'@importFrom RSpectra eigs
#'@importFrom sirt noharm.sirt
#'@importFrom FactoMineR PCA
#'@importFrom FactoMineR HCPC
#'@importFrom MASS polr
#'@importFrom ade4 dudi.pca
#'@importFrom stats cor cov dist plogis quantile rexp rnorm runif
#'@importFrom utils head
#'@section Getting Started:
#'Get started with the LatentREGpp package browsing the index of this documentation
#'if you need help the vignettes should be helpful.
#'@section Getting Started:
#'The LatentREGpp package allows you to use the LatentREGpp methodology for simulating, analyzing and scoring tests \cr
#'You can browse the package vignettes to get started.
#'@section Acknowledgment:
#'This work was Supported by Colciencias Research Grant 0039-2013 and SICS Research Group, Universidad Nacional de Colombia.
NULL


#'@name itemfit
#'@title Parameter estimation of a test
#'@description Estimates the test parameters according to the Multidimensional Item Response Theory
#'@param data The matrix containing the answers of tested individuals
#'@param dim The dimensionality of the test
#'@param model "1PL", "2PL" or "3PL"
#'@param EMepsilon Convergence value to determine the accuracy of the test
#'@param clusters A vector with cluster per dimension 
#'@param quad_tech A string with technique. "Gaussian" for Gaussian quadrature 
#'or "QMCEM" for Quasi-Monte Carlo quadrature
#'@param quad_points Amount of quadrature points. If quadratura_technique is "Gaussian". It can be NULL
#'@param individual_weights A vector with Weights of the quadrature points.
#'@param initial_values A matrix with initial values for estimation process. Be sure about
#'dimension, model and consistency with data. 
#'@param verbose True for get information about estimation process in runtime. False in otherwise. 
#'@param save_time True for save estimation time. False otherwise.
#'@examples
#' #Example 1
#'
#' dir = system.file(package="LatentREGpp")
#' folder = "/dataset/1D/dicho/"
#' file = "1000x50-1.csv"
#' data_dir = paste(c(dir, folder, file), collapse = "")
#' data = read.table(file = data_dir, sep = ";")
#' est <- itemfit(data = data, dim = 1)
#'
#' #Example 2
#'
#' #Dichotomous and multidimensional data
#' dir = system.file(package="LatentREGpp")
#' folder = "/dataset/3D/dicho/"
#' file = "1000x55-1.csv"
#' data_dir = paste(c(dir, folder, file), collapse = "")
#' data = read.table(file = data_dir, sep = ";")
#' clust <- c(20,20,15)
#' st <- itemfit(data = data, model = "2PL",dim = 3, EMepsilon = 1e-03, clusters = clust, quad_tech = "Gaussian")
#'@export
itemfit = function(data, dim, model = "2PL", EMepsilon = 1e-4, clusters = NULL,
				  quad_tech = NULL, quad_points = NULL, 
				  individual_weights = as.integer(c()),
				  initial_values = NULL,
				  verbose = TRUE, save_time = TRUE ) {
	
	# Quadrature technique
	if ( is.null(quad_tech) ) {
		if ( dim < 5 ) quad_tech = "Gaussian"
		else quad_tech = "QMCEM"
	} else {
		if ( quad_tech != "Gaussian" && quad_tech != "QMCEM" )
			stop("Quadrature technique not found")

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

	# Type of data
	dichotomous_data = is_data_dicho(data)
	if ( dichotomous_data == -1 )
		stop("Inconsistent data")

	# Model id
	if ( model == "1PL" ) m = 1
	else if ( model == "2PL" ) m = 2
	else if ( model == "3PL" ) m = 3

	q = quadpoints(dim = dim, quad_tech = quad_tech, 
				   quad_points = quad_points)
	theta = q$theta
	weights = q$weights

	if ( dim == 1 ) {
		# Item parameters estimation
		obj_return = itemfitcpp(Rdata = data, dim = dim, model = m, EMepsilon = EMepsilon, 
							Rtheta = theta, Rweights = weights, 
							Rindividual_weights = individual_weights,
							dichotomous_data = dichotomous_data, 
							verbose = verbose )
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

				acp = PCA(X = III$coefs,graph = FALSE)
				hcpc = HCPC(acp,nb.clust = d,graph = FALSE)
				CLUST_final<- list()
				for (i in 1:length(table(hcpc$data.clust$clust))) {
					CLUST_final[[i]]<- data[,which(hcpc$data.clust$clust==i)]
				}

				clusters = unlist(lapply(CLUST_final, ncol))
			} else
				if ( length(clusters) != dim )
					stop("Clusters length must be equal to the number of dimensions")

			#Initial values
			if ( is.null(initial_values) ) {
				#Find initial values again
				initial_values = inivals_MultiUni_NOHARM(data, clusters, model=model, 
				                          find.restrictions=FALSE, verbose=FALSE, probit=FALSE)$coefs
			} else {
				initial_values = data.matrix(initial_values)
				if ( nrow(initial_values) != ncol(data) )
					stop("Inconsistent initial_values. Number of rows must be equal to the number of items")
			}
		  
			# Item parameters estimation
			obj_return = itemfitcpp(Rdata = data, dim = dim, model = m, EMepsilon = EMepsilon, 
								Rtheta = theta, Rweights = weights, 
								Rindividual_weights = individual_weights,
								dichotomous_data = dichotomous_data,
								Rclusters = clusters,
								Rinitial_values = initial_values, 
								verbose = verbose)
		} else {
			if ( is.null(clusters) )
				stop("You must specify clusters")

			#Initial values
			if ( is.null(initial_values) ) {
				#Find initial values again
				initial_values = inivals_MultiPoly(data_poly = data, size.cluster = clusters, verbose=F)
			} else {
				initial_values = data.matrix(initial_values)
				if ( nrow(initial_values) != ncol(data) )
					stop("Inconsistent initial_values. Number of rows must be equal to the number of items")
			}

			# Item parameters estimation
			obj_return = itemfitcpp(Rdata = data, dim = dim, model = m, EMepsilon = EMepsilon, 
									Rtheta = theta, Rweights = weights, 
									Rindividual_weights = individual_weights,
									dichotomous_data = dichotomous_data,
									Rclusters = clusters,
									Rinitial_values = initial_values, 
									verbose = verbose )
		}
	}

	if ( save_time ) {
		second_time = Sys.time()
		obj_return$time = second_time - first_time	
	}

	obj_return$dimension = dim
	obj_return$model = model
	obj_return$clusters = clusters
	obj_return$convergence = obj_return$iterations < 500
	obj_return$epsilon = EMepsilon
  	
  	if ( dichotomous_data )
		colnames(obj_return$zetas) = c(paste("a",c(1:obj_return$dimension),sep = ""),"d","c")
	else {
		max_ncat = ncol(obj_return$zetas) - obj_return$dimension
		colnames(obj_return$zetas)<- c(paste("a",c(1:obj_return$dimension), sep=""),
			paste("d",c(1:max_ncat), sep=""))
	}
  
	return (obj_return)
}