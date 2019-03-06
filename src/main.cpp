////////////////////////////////////////////////////////////////////////////////////////
///
///	Created by Thomas Alauzet on November 12, 2018.
/// Copyright 2018. All rights reserved.
///
////////////////////////////////////////////////////////////////////////////////////////

#include <omp.h>
#include <cstdio>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "bfgs.h"
#include "class.h"
#include "parser.h"
#include "recuit.h"
#include "neldermead.h"
#include "objectivefunction.h"

using namespace std;


// Input File :
// ndim
// npoint
// optim_method 		NM	BFGS	RS
// estimation_method 	MC	NN	MM	AE	KL	MST

// Or arguments :
// ./spacefilling ndim ndpoint RS MST

// export OMP_NUM_THREADS=4


void launchOptimization(optim_method optim_method, estimation_method estimation_method,
						int ndim, int npoint, double q, double* param, double* param_min,
						double* param_max, double& f, double* fVector)
{
	////////////////////////////////////////////////////////////////////////////////////////////////    
	// Create the function to optimize
	////////////////////////////////////////////////////////////////////////////////////////////////
	
	Function* function = new Function(estimation_method, ndim, npoint, q);


	////////////////////////////////////////////////////////////////////////////////////////////////
	// Create the optimization class and minimize
	////////////////////////////////////////////////////////////////////////////////////////////////


	if (optim_method == P_BFGS)
	{
		Bfgs* optimization = new Bfgs(function, ndim, npoint);
		
		optimization->minimize(param, param_min, param_max);
		
		delete optimization;
	}
	else if (optim_method == P_NELDER_MEAD)
	{
		NelderMead* optimization = new NelderMead(function, ndim, npoint);
		
		optimization->minimize(param, param_min, param_max);
		
		delete optimization;
	}
	else if (optim_method == P_RECUIT_SIMULE)
	{
		Recuit* optimization = new Recuit(function, ndim, npoint);
		
		optimization->minimize(param, param_min, param_max);
		
		delete optimization;
	}


	f = function->calcf(param);
	function->calcFVector(param, fVector);

	delete function;
}


int main(int argc, char *argv[])
{
	stdout = freopen("../output.txt", "w", stdout);

	cout << "##############################################################################################" << endl;
	cout << "##############################################################################################" << endl;
	cout << "#################################### Space-Filling Design ####################################" << endl;
	cout << "##############################################################################################" << endl;
	cout << "##############################################################################################" << endl;

	int					sys					= 0;
	int					dim					= (argc > 1) ? atoi(argv[1]) : getDim();
	int					npoint				= (argc > 2) ? atoi(argv[2]) : getPoint();
	double				q					= 0.9;
	double				f					= 0.0;
	double*				a					= new double[npoint * dim];
	double*				b					= new double[npoint * dim];
	double*				param				= new double[npoint * dim];
	double*				fVector				= new double[6];
	optim_method		optim_method		= P_RECUIT_SIMULE;
	estimation_method	estimation_method	= P_MINIMAL_SPANNING_TREE;


	parse(dim, npoint, a, b, optim_method, estimation_method, argc, argv);
	initializeCenter(dim, npoint, param);
	double time = omp_get_wtime();
	launchOptimization(optim_method, estimation_method, dim, npoint, q, param, a, b, f, fVector);
	time = omp_get_wtime() - time;
	printResults(dim, npoint, param, f);


	// Print Time
	stdout = freopen ("/dev/tty", "a", stdout);
	cout << "Time : " << time << endl;


	// Python Plot
	string command = "python ../python/plot.py "  + intToString(dim) + " " + intToString(npoint);
	sys = system(command.c_str());


	// Benchmark Output
	string fileName	= "../Results/";
	fileName += argv[4];
	fileName += "_";
	fileName += doubleToString(q);
	fileName += "_";
	fileName += intToString(dim);
	fileName += "_";
	fileName += intToString(npoint);
	fileName += "_";
	fileName += argv[5];
	printFinalResults(dim, npoint, 6, param, fVector, fileName);


	delete[] a;
	delete[] b;
	delete[] param;
	delete[] fVector;

	return(sys);
}