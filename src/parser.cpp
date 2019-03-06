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


void parse(int ndim, int npoint, double* a, double* b, optim_method& optim_method, estimation_method& estimation_method, int argc, char* argv[])
{
	for (int i = 0; i < ndim * npoint; i++)
	{
		a[i] = 0.0;
		b[i] = 1.0;
	}

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
	else
	{
		ifstream file("../input.txt");
	
		if (file)
		{
			int		ndim0;
			int		npoint0;
			string	temp1;
			string	temp2;
	
			file >> ndim0 >> npoint0 >> temp1 >> temp2;
	
			if		(temp1 == "NM")			optim_method		= P_NELDER_MEAD;
			else if	(temp1 == "BFGS")		optim_method		= P_BFGS;
			else if	(temp1 == "RS")			optim_method		= P_RECUIT_SIMULE;
	
			if		(temp2 == "MC")		estimation_method	= P_MONTE_CARLO;
			else if	(temp2 == "NN")		estimation_method	= P_NEAREST_NEIGHBOR;
			else if	(temp2 == "MM")		estimation_method	= P_MAXIMIN;
			else if	(temp2 == "AE")		estimation_method	= P_AE;
			else if	(temp2 == "KL")		estimation_method	= P_KL;
			else if	(temp2 == "MST")	estimation_method	= P_MINIMAL_SPANNING_TREE;
			else if	(temp2 == "ROS")	estimation_method	= P_ROSENBROCK;

			file.close();
		}
		else
		{
			cerr << "Error reading filename" << endl;
		}
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