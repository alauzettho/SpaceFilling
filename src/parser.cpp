////////////////////////////////////////////////////////////////////////////////////////
///
///	Created by Thomas Alauzet on November 12, 2018.
/// Copyright 2018. All rights reserved.
///
////////////////////////////////////////////////////////////////////////////////////////

#include "parser.h"

using namespace std;


string intToString(int number)
{
	stringstream stream;
	stream << setprecision(17) << number;
	string s = stream.str();
	return(s);
}


string doubleToString(double number)
{
	stringstream stream;
	stream << setprecision(17) << number;
	string s = stream.str();
	return(s);
}


void parse(int ndim, int npoint, double* a, double* b, optim_method& optim_method, estimation_method& estimation_method, int argc, char* argv[])
{
	optim_method		= P_RECUIT_SIMULE;
	estimation_method	= P_MINIMAL_SPANNING_TREE;

	if (argc > 4)
	{
		string arg4 (argv[4]);

		if		(arg4 == "MC")	estimation_method	= P_MONTE_CARLO;
		else if	(arg4 == "NN")	estimation_method	= P_NEAREST_NEIGHBOR;
		else if	(arg4 == "MM")	estimation_method	= P_MAXIMIN;
		else if	(arg4 == "AE")	estimation_method	= P_AE;
		else if	(arg4 == "KL")	estimation_method	= P_KL;
		else if	(arg4 == "MST")	estimation_method	= P_MINIMAL_SPANNING_TREE;
		else if	(arg4 == "ROS")	estimation_method	= P_ROSENBROCK;
	}

	if (argc > 3)
	{
		string arg3 (argv[3]);

		if		(arg3 == "NM")		optim_method	= P_NELDER_MEAD;
		else if	(arg3 == "BFGS")	optim_method	= P_BFGS;
		else if	(arg3 == "RS")		optim_method	= P_RECUIT_SIMULE;
	}
}


void initializeCenter(int ndim, int npoint, double* param)
{
	time_t timer;
	time(&timer);
	srand(timer);

	for (int i = 0; i < ndim * npoint; i++)
	{
		param[i] = (rand() % 1000000) / 1000000.0;
	}
}


void printResults(int ndim, int npoint, double* param, double f)
{
	cout << endl << "F :	" + doubleToString(f) << endl;
	cout << endl << "Param :" << endl;

	for (int i = 0; i < ndim; i++)
	{
		for (int j = 0; j < npoint; j++)
		{
			cout << doubleToString(param[i * npoint + j]) + "	";
		}
		cout << endl;
	}
}


void printFinalResults(int ndim, int npoint, int numberEstim, double* param, double* fVector, string fileName1)
{
	string fileName2	= fileName1 + "_F";
	const char* temp1	= fileName2.c_str();
	const char* temp2	= fileName1.c_str();
	
	stdout = freopen(temp1, "w", stdout);

	for (int i = 0; i < numberEstim; i++)
	{
		cout << doubleToString(fVector[i]) + "	";
	}

	stdout = freopen(temp2, "w", stdout);

	for (int i = 0; i < ndim; i++)
	{
		for (int j = 0; j < npoint; j++)
		{
			cout << doubleToString(param[i * npoint + j]) + "	";
		}
		cout << endl;
	}
}