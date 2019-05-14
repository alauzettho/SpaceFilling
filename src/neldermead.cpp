////////////////////////////////////////////////////////////////////////////////////////
///
/// Created by Thomas Alauzet on November 14, 2018.
/// Copyright 2018. All rights reserved.
///
////////////////////////////////////////////////////////////////////////////////////////

#include "neldermead.h"

using namespace std;


NelderMead::NelderMead(Function* function, int dim, int npoint) : m_function (function), m_ndim (dim), m_npoint (npoint)
{
	m_nparam		= m_ndim * m_npoint;
	m_simplex_size	= m_nparam + 1;
	m_nitermax		= 10000000;
	m_nfunevalmax	= 10000000;
	m_tolf			= 1e-15;
	m_alpha			= 1.0;
	m_beta			= 1.00 + 2.0 / double(m_nparam);
	m_gamma			= 0.75 - 0.5 / double(m_nparam);
	m_delta			= 1.00 - 1.0 / double(m_nparam);
	m_tau1			= 0.05;
	m_tau2			= 0.00025;
	m_tolx			= new double[m_nparam];

	for (int i = 0; i < m_nparam; i++)
	{
		m_tolx[i]	= m_tolf;
	}
}


NelderMead::~NelderMead()
{
	delete[] m_tolx;
}


void NelderMead::minimize(double* param, const double* min_param, const double* max_param)
{
	// cout << "##############################################################################################" << endl;

	int		iter			= 0;
	int		test			= -1;	// Value corresponding to the Nelder_Mead test that occurs at each iteration
	int		nfuneval		= 0;	// Number of function computation
	double	fr				= 0.0;	// Reflexion value
	double	f_new			= 0.0;	// Other function computation
	string	message;

	// cout << "##################################### Starting Nelder Mead ###################################" << endl;

	// Array to save function values
	vector<double> f (m_simplex_size, 1e+15);

	// Nelder-Mead initial simplex
	vector <vector<double> > e; 
	vector <vector<double> > x;

	// Base
	for (int i = 0; i < m_nparam; i++)
	{
		e.push_back(vector<double>(m_nparam, 0.0));
		e[i][i] = 1.0;
	}

	x.push_back(vector<double>(m_nparam, 0.0));
	
	for (int j = 0; j < m_nparam; j++)
	{
		x[0][j] = param[j];
	}

	for (int i = 1; i < m_nparam + 1; i++)
	{
		x.push_back(vector<double>(m_nparam, 0.0));
		for (int j = 0; j < m_nparam; j++)
		{
			x[i][j] = param[j] * (1.0 + m_tau1 * e[i - 1][j]);
		}
		
		if (fabs(x[i][i - 1]) < 1e-15)
		{
			x[i][i-1] = m_tau2;
		}
	}


	// Evaluate and sort initial simplex
	evalSimplex(x, &f[0], nfuneval);
	sortSimplex(x, &f[0]);


	while (iter < m_nitermax && nfuneval < m_nfunevalmax && !testF(&f[0]) && !testX(x))
	{
		iter++;

		// INIT
		vector<vector<double> > new_points;
		for (int i = 0; i < 4; i++)
		{
			new_points.push_back(vector<double>(m_nparam, 0.0));
		}

		// CENTROID (excluding the worst point at position m_simplex_size)
		vector<double> centroid(m_nparam, 0.0);
		for (int i = 0; i < m_nparam; i++)
		{
			for (int j = 0; j < m_nparam; j++)
			{
				centroid[i] += x[j][i];
			}
			centroid[i] *= 1.0 / double(m_nparam);
		}

		// REFLECTION
		for (int i = 0; i < m_nparam; i++)
		{
			new_points[REFLECTION][i] = centroid[i] + m_alpha * (centroid[i] - x[m_nparam][i]);
		}

		// EXPANSION
		for (int i = 0; i < m_nparam; i++)
		{
			new_points[EXPANSION][i] = centroid[i] + m_beta * (new_points[REFLECTION][i] - centroid[i]);
		}

		// OUTSIDE CONTRACTION
		for (int i = 0; i < m_nparam; i++)
		{
			new_points[CONTRACTION_OUT][i] = centroid[i] + m_gamma * (new_points[REFLECTION][i] - centroid[i]);
		}

		// INSIDE CONTRACTION
		for (int i = 0; i < m_nparam; i++)
		{
			new_points[CONTRACTION_IN][i] = centroid[i] - m_gamma * (centroid[i] - x[m_nparam][i]);
		}

		// Compute f at reflexion
		fr = m_function->calcf(&new_points[REFLECTION][0]);
		nfuneval++;


		// Testing which point will be selected
		test = ERROR;
		if		(fr >= f[m_nparam])							test = CONTRACTION_IN;
		else if	(fr >= f[m_nparam - 1] && fr < f[m_nparam])	test = CONTRACTION_OUT;
		else if	(fr < f[0])									test = EXPANSION;
		else if	(fr >= f[0] && fr < f[m_nparam - 1])		test = REFLECTION;
		assert(test != -1);


		// Compute new value (Not fr)
		if (test == REFLECTION)
		{
			f_new = fr;
		}
		else
		{
			f_new = m_function->calcf(&new_points[test][0]);
			nfuneval++;
		}


		// Next part updates Simplexe
		if (test == REFLECTION)
		{
			// message = "Iter: " + intToString(iter) + " REFLECTION with param: ";
			// for (int i = 0; i < m_nparam; i++) message += doubleToString(new_points[EXPANSION][i]) + " ";
			// message += " F: " + doubleToString(f_new) + " Funeval: " + intToString(nfuneval);
			// cout << endl << message << endl;

			updateSimplex(new_points[REFLECTION], f_new, x, &f[0]);
		}

		else if (test == EXPANSION)
		{
			// message = "Iter: " + intToString(iter) + " EXPANSION with param: ";
			// for (int i = 0; i < m_nparam; i++) message += doubleToString(new_points[EXPANSION][i]) + " ";
			// message += " F: " + doubleToString(f_new) + " Funeval: " + intToString(nfuneval);
			// cout << endl << message << endl;

			if (f_new < fr)
			{
				updateSimplex(new_points[EXPANSION], f_new, x, &f[0]);
			}
			else
			{
				updateSimplex(new_points[REFLECTION], fr, x, &f[0]);
			}
		}

		else if (test == CONTRACTION_OUT)
		{
			// message = "Iter: " + intToString(iter) + " CONTRACTION_OUT with param: ";
			// for (int i = 0; i < m_nparam; i++) message += doubleToString(new_points[EXPANSION][i]) + " ";
			// message += " F: " + doubleToString(f_new) + " Funeval: " + intToString(nfuneval);
			// cout << endl << message << endl;

			if (f_new <= fr)
			{
				updateSimplex(new_points[CONTRACTION_OUT], f_new, x, &f[0]);
			}
			else
			{
				shrinkSimplex(x);					// x[0] do not change but other points get closer to minimum
				evalSimplex(x, &f[0], nfuneval);	// Re-evaluate simplexe
				sortSimplex(x, &f[0]);				// Sort simplexe values
			}
		}

		else if (test == CONTRACTION_IN)
		{
			// message = "Iter: " + intToString(iter) + " CONTRACTION_IN with param: ";
			// for (int i = 0; i < m_nparam; i++) message += doubleToString(new_points[EXPANSION][i]) + " ";
			// message += " F: " + doubleToString(f_new) + " Funeval: " + intToString(nfuneval);
			// cout << endl << message << endl;

			if (f_new < f[m_nparam])
			{
				updateSimplex(new_points[CONTRACTION_IN], f_new, x, &f[0]);
			}
			else
			{
				shrinkSimplex(x);					// x[0] do not change but other points get closer to minimum
				evalSimplex(x, &f[0], nfuneval);	// Re-evaluate simplexe
				sortSimplex(x, &f[0]);				// Sort simplexe values
			}
		}

		// Update param with the best point of the final simplex
		for (int j = 0; j < m_nparam; j++) param[j] = x[0][j];

		double fval = f[0];

		// Updating Output
		printResults(m_ndim, m_npoint, param, fval);
	}


	// if	(nfuneval >= m_nfunevalmax)	cout << "nfuneval = " 	+ doubleToString(nfuneval)	+ " > maxfuneval = "	+ doubleToString(m_nfunevalmax)	<< endl;
	// if	(iter >= m_nitermax)		cout << "iter = "		+ doubleToString(iter)		+ " > maxiter = "		+ doubleToString(m_nitermax)	<< endl;

	// cout << "###################################### Ending Nelder Mead ####################################" << endl;
	// cout << "##############################################################################################" << endl;
}


bool NelderMead::testX(const vector<vector<double> >& x)
{
	bool isfinished = true;

	for (int isimplex = 1; isimplex < m_simplex_size ; isimplex++)
	{
		for (int iparam = 0; iparam < m_nparam; iparam++)
		{
			if (fabs(x[isimplex][iparam] - x[0][iparam]) > m_tolx[iparam]) isfinished = false;
		}
	}

	if (isfinished)
	{
		cout << endl << "Stopping optimization: simplex dimensions < user defined tolx" << endl;
	}

	return isfinished;
}


bool NelderMead::testF(double* f)
{
	double diff;
	double maxf = fabs(f[1] - f[0]);

	for (int i = 2; i < m_simplex_size; i++)
	{
		diff = fabs(f[i] - f[0]);
		if (diff > maxf) maxf = diff;
	}
	if (maxf <= m_tolf)
	{
		cout << endl << "Stopping optimization: deltaf = " << doubleToString(maxf) << " < tolf = " << doubleToString(m_tolf) << endl;
		return true;
	}
	else
	{
		return false;
	}
}


void NelderMead::evalSimplex(vector<vector<double> >& x, double* f, int& nfuneval)
{
	for (int i = 0; i < m_simplex_size; i++)
	{
		f[i] = m_function->calcf(&x[i][0]);
		nfuneval++;
	}
}


void NelderMead::updateSimplex(const vector<double>& new_point, const double new_f, vector<vector<double> >& x, double* f)
{
	for (int i = 0; i < m_nparam; i++)
	{
		x[m_nparam][i] = new_point[i];
	}

	f[m_nparam] = new_f;

	insertion(x, f);
}


void NelderMead::shrinkSimplex(vector<vector<double> >& x)
{
	for (int i = 1; i < m_simplex_size; i++)
	{
		for (int j = 0; j < m_nparam; j++)
		{
			x[i][j] = x[0][j] + m_delta * (x[i][j] - x[0][j]);
		}
	}
}


void NelderMead::insertion(vector<vector<double> >& x, double* f)
{
	bool	is_finished	= false;
	double	fn			= f[m_nparam];
	double*	temp		= new double[m_nparam];

	for (int i = 0; i < m_nparam; i++)
	{
		temp[i] = x[m_nparam][i];
	}

	int i = m_nparam - 1;

	while (i >= 0 && !is_finished)    
	{
		if (fn < f[i] && i!=0)
		{
			f[i + 1] = f[i];
			for (int j = 0; j < m_nparam; j++)
			{
				x[i + 1][j] = x[i][j];
			}
		}
		else if (fn < f[i] && i == 0)
		{
			f[i + 1] = f[i];
			for (int j = 0; j < m_nparam; j++)
			{
				x[i + 1][j] = x[i][j];
			}
			f[0] = fn;
			for (int j = 0; j < m_nparam; j++)
			{
				x[0][j] = temp[j];
			}
		}
		else
		{
			f[i + 1] = fn;
			for (int j = 0; j < m_nparam; j++)
			{
				x[i + 1][j] = temp[j];
			}
			is_finished = true;
		}
		i--;
	}

	delete [] temp;
}


void NelderMead::sortSimplex(vector<vector<double> >& x, double *f)
{
	int		ntemp2;
	double	ntemp1;
	int*	count = new int[m_nparam + 1];
	bool	order = true;

	for (int i = 0; i < m_simplex_size; i++)
	{
		count[i] = i;
	}

	for (int i = 0; i < m_simplex_size && order; i++)
	{
		order = false;
		for (int j = 1; j < m_simplex_size- i; j++)
		{
			if (f[j - 1] > f[j])
			{
				ntemp1			= f[j];
				f[j]			= f[j - 1];
				f[j - 1]		= ntemp1;
				ntemp2			= count[j];
				count[j]		= count[j-1];
				count[j-1]		= ntemp2;
				order			= true;
			}
		}
	}

	vector<vector<double> > copyx;
	for (int i = 0; i < m_simplex_size; i++)
	{
		copyx.push_back(vector<double>(m_nparam,0));
		for (int j = 0; j < m_nparam; j++)
		{
			copyx[i][j] = x[i][j];
		}
	}

	for (int i = 0; i < m_simplex_size; i++)
	{
		for (int j = 0; j < m_nparam; j++)
		{
			x[i][j] = copyx[count[i]][j];
		}
	}

	delete[] count;
}