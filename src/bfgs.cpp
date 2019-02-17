////////////////////////////////////////////////////////////////////////////////////////
///
/// Created by Thomas Alauzet on November 14, 2018.
/// Copyright 2018. All rights reserved.
///
////////////////////////////////////////////////////////////////////////////////////////

#include "bfgs.h"

using namespace std;


Bfgs::Bfgs(Function* function, int ndim, int npoint) : m_function (function), m_ndim (ndim), m_npoint (npoint)
{
	m_maxiter		= 10000;
	m_maxfuneval	= 10000000;
	m_tolg			= 0.0000001;
	m_nparam		= m_npoint * m_ndim;
}


Bfgs::~Bfgs(){}


void Bfgs::minimize(double* param, const double* min_param, const double* max_param)
{
	cout << "##############################################################################################" << endl;
	int		iter			= 0;
	int		nfuneval		= 0;
	double  f				= 0.0;
	double  ys				= 1.0;
	double  alpha			= 0.0;
	double* g				= new double[m_nparam];
	double* d				= new double[m_nparam];
	double* s				= new double[m_nparam];
	double* y				= new double[m_nparam];
	double* gold			= new double[m_nparam];
	vector<double>			line(m_nparam, 0.0);
	vector<vector<double> >	H;
	vector<vector<double> >	P;
	vector<vector<double> >	sst;
	vector<vector<double> >	mat1;
	vector<vector<double> >	mat2;

	for (int i = 0; i < m_nparam; i++)
	{
		H.push_back(line);
		sst.push_back(line);
		mat1.push_back(line);
		mat2.push_back(line);
		P.push_back(line);
	}

	cout << "######################################## Starting BFGS #######################################" << endl;

	f = m_function->calcf(param);
	nfuneval++;

	// Updating Output
	// printResults(m_ndim, m_npoint, param, f);

	m_function->calcGradient(param, f, g, nfuneval);


	while (iter < m_maxiter && nfuneval < m_maxfuneval && normInf(m_nparam, g) > m_tolg && ys > 0)
	{
		bool temp = true;

		if (temp) //iter == 0)
		{
			// Uses identity matrice has hessian model

			for (int i = 0; i < m_nparam; i++)
			{
				for (int j = 0; j < m_nparam; j++) H[i][j] = 0.0;
				H[i][i] = 1.0;
				d[i] = -g[i];
			}
		}
		else
		{
			// Update hessian has BFGS states
			
			if (iter == 1)
			{
				for (int i = 0; i < m_nparam; i++)
				{
					H[i][i] = scalarProduct(m_nparam, y, s) / scalarProduct(m_nparam, y, y);
				}
			}

			for (int i = 0; i < m_nparam; i++)
			{
				for (int j = 0; j < m_nparam; j++)
				{
					sst[i][j]	=  s[i] * s[j] / scalarProduct(m_nparam, y, s);
					mat1[i][j]	= -s[i] * y[j] / scalarProduct(m_nparam, y, s);
					mat2[i][j]	= -y[i] * s[j] / scalarProduct(m_nparam, y, s);
				}
			}

			for (int i = 0; i < m_nparam; i++)
			{
				mat1[i][i] += 1.0;
				mat2[i][i] += 1.0;
			}

			// Now Update: H = mat1 * H * mat2 + sst

			for (int i = 0; i < m_nparam; i++)
			{
				for (int j = 0; j < m_nparam; j++)
				{
					P[i][j] = 0.0;
					for (int k = 0; k < m_nparam; k++)
					{
						P[i][j] += mat1[i][k] * H[k][j];
					}
				}
			}

			for (int i = 0; i < m_nparam; i++)
			{
				for (int j = 0; j < m_nparam; j++)
				{
					H[i][j] = 0.0;
					for (int k = 0; k < m_nparam; k++)
					{
						H[i][j] += P[i][k] * mat2[k][j];
					}
				}
			}

			for (int i = 0; i < m_nparam; i++)
			{
				for (int j = 0; j < m_nparam; j++)
				{
					H[i][j] += sst[i][j];
				}
			}

			// Update descent direction : d = - H * g
			for (int i = 0; i < m_nparam; i++)
			{
				d[i] = 0.0;
				for (int j = 0; j < m_nparam; j++) d[i] -= H[i][j] * g[j];
			}
		}


		getStep(iter, param, min_param, max_param, d, alpha, f, nfuneval);

		//if (alpha < 1e-15) break;

		// Updating values for next Iteration
		for (int j = 0; j < m_nparam; j++)
		{
			s[j]		= alpha * d[j];
			param[j]	+= s[j];
		}


		// Updating Output
		// printResults(m_ndim, m_npoint, param, f);


		if (iter + 1 < m_maxiter && nfuneval + m_nparam < m_maxfuneval)
		{
			for (int j = 0; j < m_nparam; j++) gold[j] = g[j];
			m_function->calcGradient(param, f, g, nfuneval);
			for (int j = 0; j < m_nparam; j++) y[j] = g[j] - gold[j];
		}

		ys = scalarProduct(m_nparam, y, s);

		ys = 1.0;

		iter++;
	}

	cout << "######################################### Ending BFGS ########################################" << endl;
	cout << "##############################################################################################" << endl;
	if		(normInf(m_nparam, g) <= m_tolg)	cout << "||g|| = "		+ doubleToString(normInf(m_nparam, g)) + " < tolg = " + doubleToString(m_tolg) << endl;
	else if	(nfuneval >= m_maxfuneval)			cout << "nfuneval = " 	+ doubleToString(nfuneval) + " > maxfuneval = " + doubleToString(m_maxfuneval) << endl;
	else if	(iter >= m_maxiter)					cout << "iter = "		+ doubleToString(iter) + " > maxiter = " + doubleToString(m_maxiter) << endl;
	else if	(ys <= 0)							cout << "Y'S = "		+ doubleToString(ys) + " < 0" << endl;
	else										cout << "Unknown stop" << endl;


	delete[] s;
	delete[] g;
	delete[] y;
	delete[] d;
	delete[] gold;
}


void Bfgs::getStep(int iter, double* param, const double* min_param, const double* max_param, double* d, double& alpha, double& f, int& nfuneval)
{
	// TODO : Upgrade to Wolfe Armijo Goldstein

	if (iter == 0)
	{
		alpha = 1.0 / normInf(m_nparam, d);
		if (alpha > 1.0) alpha = 1.0;
	}
	else
	{
		alpha = 2.0;
	}

	for (int i = 0; i < m_nparam; i++)
	{
		if (param[i] - alpha * d[i] > max_param[i % m_ndim])
		{
			alpha = fabs((param[i] - max_param[i % m_ndim]) / d[i]);
		}
		if (param[i] - alpha * d[i] < min_param[i % m_ndim])
		{
			alpha = fabs((param[i] - min_param[i % m_ndim]) / d[i]);
		}
	}
	
	double	fnew = f;
	double* new_param = new double[m_nparam];

	while (fnew >= f)
	{
		if (alpha < 1e-15) break;

		alpha /= 2.0;
		for (int i = 0; i < m_nparam; i++)
		{
			new_param[i] = param[i] + alpha * d[i];
		}

		fnew = m_function->calcf(new_param);
		nfuneval++;
	}

	f = fnew;
}


double Bfgs::scalarProduct(int size, double* x, double* y)
{
	double scalar = 0.0;

	for (int i = 0; i < size; i++)
	{
		scalar += x[i] * y[i];
	}

	return(scalar);
}


double Bfgs::normInf(int size, double* vector)
{
	double norm = fabs(vector[0]);

	for (int i = 1; i < size; i++)
	{
		if (fabs(vector[i]) > norm) norm = fabs(vector[i]);
	}

	return(norm);
}