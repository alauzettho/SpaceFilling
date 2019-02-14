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
// estimation_method 	MC	MST		NN	ROS

// export OMP_NUM_THREADS=4

// Complexite : maxiter * npoint * npoint C(num_threads, ndim)


void launchOptimization(optim_method optim_method, estimation_method estimation_method,
						int ndim, int npoint, double* param, double* param_min, double* param_max, double& f)
{
	////////////////////////////////////////////////////////////////////////////////////////////////    
	// Create the function to optimize
	////////////////////////////////////////////////////////////////////////////////////////////////
	
	Function* function = new Function(estimation_method, ndim, npoint);


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

	int					sys		= 0;
	int					dim		= getDim();
	int					npoint	= getPoint();
	double				f		= 0.0;
	double*				a		= new double[npoint * dim];
	double*				b		= new double[npoint * dim];
	double*				param	= new double[npoint * dim];
	optim_method		optim_method;
	estimation_method	estimation_method;


	parse(dim, npoint, a, b, optim_method, estimation_method);

	initializeCenter(dim, npoint, param);

	double time = omp_get_wtime();

	launchOptimization(optim_method, estimation_method, dim, npoint, param, a, b, f);

	time = omp_get_wtime() - time;

	printResults(dim, npoint, param, f);


	stdout = freopen ("/dev/tty", "a", stdout);
	cout << "Time : " << time << endl;


	if (dim == 2)
	{
		string command = "python ../python/plot2D.py";
		sys = system(command.c_str());
	}
	else if (dim == 3)
	{
		string command = "python ../python/plot3D.py";
		sys = system(command.c_str());
	}
	else
	{
		string command = "python ../python/plotND.py";
		sys = system(command.c_str());
	}

	return(sys);
}