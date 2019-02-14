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
	stream << setprecision(13) << number;
	string s = stream.str();
	return(s);
}


string doubleToString(double number)
{
	stringstream stream;
	stream << setprecision(13) << number;
	string s = stream.str();
	return(s);
}


int getDim()
{
	int ndim = 0;

	ifstream file("../input.txt");

	if (file)
	{
		file >> ndim;
	}
	else
	{
		cerr << "Error reading filename" << endl;
	}

	assert (ndim != 0 && "consider changing input file");

	return(ndim);
}


int getPoint()
{
	int ndim	= 0;
	int npoint	= 0;

	ifstream file("../input.txt");

	if (file)
	{
		file >> ndim >> npoint;
	}
	else
	{
		cerr << "Error reading filename" << endl;
	}

	assert (npoint != 0 && "consider changing input file");

	return(npoint);
}


void parse(int ndim, int npoint, double* a, double* b, optim_method& optim_method, estimation_method& estimation_method)
{
	// Defaults values
	for (int i = 0; i < ndim * npoint; i++)
	{
		a[i] = 0.0;
		b[i] = 1.0;
	}

	optim_method		= P_BFGS;
	estimation_method	= P_MONTE_CARLO;


	// Import values
	ifstream file("../input.txt");

	if (file)
	{
		string temp1;
		string temp2;

		file >> ndim >> npoint >> temp1 >> temp2;

		if		(temp1 == "NM")			optim_method		= P_NELDER_MEAD;
		else if	(temp1 == "BFGS")		optim_method		= P_BFGS;
		else if	(temp1 == "RS")			optim_method		= P_RECUIT_SIMULE;


		if		(temp2 == "MC")			estimation_method	= P_MONTE_CARLO;
		else if	(temp2 == "NN")			estimation_method	= P_NEAREST_NEIGHBOR;
		else if	(temp2 == "MST")		estimation_method	= P_MINIMAL_SPANNING_TREE;
		else if	(temp2 == "ROS")		estimation_method	= P_ROSENBROCK;


		file.close();
	}
	else
	{
		cerr << "Error reading filename" << endl;
	}
}


void initializeCenter(int ndim, int npoint, double* param)
{
	time_t timer;
	time(&timer);
	srand(timer);

	for (int i = 0; i < ndim * npoint; i++)
	{
		param[i] = (rand() % 100) / 100.0;
	}

	// param[0] = 0.001;
	// param[1] = 0.002;
	// param[2] = 0.999;
	// param[3] = 0.998;
	// param[4] = 0.003;
	// param[5] = 0.997;
	// param[6] = 0.996;
	// param[7] = 0.004;
	// param[8] = 0.01; 
	// param[9] = 0.05;
	// param[10] = 0.51;
	// param[11] = 0.01;

}


void printResults(int ndim, int npoint, double* param, double f)
{
	cout << "##############################################################################################" << endl;

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