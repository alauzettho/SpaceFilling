////////////////////////////////////////////////////////////////////////////////////////
///
/// Created by Thomas Alauzet on November 28, 2018.
/// Copyright 2018. All rights reserved.
///
////////////////////////////////////////////////////////////////////////////////////////

#include "recuit.h"

using namespace std;


Recuit::Recuit(Function* function, int ndim, int npoint) : m_function (function), m_ndim (ndim), m_npoint (npoint)
{
	// TODO : find better param
	m_temp_init		= 10000;
	m_lambda		= 0.999;
	m_maxiter		= 5000;
	m_nparam		= m_npoint * m_ndim;
}


Recuit::~Recuit(){}


void Recuit::minimize(double* param, const double* min_param, const double* max_param)
{
	// cout << "##############################################################################################" << endl;

	int		iter			= 0;
	double	e				= m_function->calcf(param);
	double	m				= e;
	double	en				= e;
	double	temp			= m_temp_init;
	double*	s				= new double[m_nparam];
	double*	g				= new double[m_nparam];
	double*	sn				= new double[m_nparam];

	
	#pragma omp parallel for
	for (int i = 0; i < m_nparam; i++)
	{
		s[i] = param[i];
		g[i] = param[i];
	}


	// cout << "####################################### Starting Recuit ######################################" << endl;


	while (iter < m_maxiter)
	{
		voisin(iter, s, sn);

		en = m_function->calcf(sn);

		// Metropolis
		if (en < e || (rand() / double(RAND_MAX)) < exp((e - en) / temp))
		{
			for (int i = 0; i < m_nparam; i++)
			{
				s[i] = sn[i];
			}

			e = en;
		}

		// Global Optimum
		if (e < m)
		{
			for (int i = 0; i < m_nparam; i++)
			{
				g[i] = s[i];
			}

			m = e;
		}

		// Updating Output
		// printResults(m_ndim, m_npoint, s, e);

		temp *= m_lambda; 
		iter++;
	}


	// cout << "######################################## Ending Recuit #######################################" << endl;

	for (int i = 0; i < m_nparam; i++)
	{
		param[i] = g[i];
	}

	delete[] s;
	delete[] g;

	// cout << "##############################################################################################" << endl;
}


void Recuit::voisin(int iter, double* param, double* paramVoisin)
{
	double product = 1.0 / sqrt(m_npoint);

	for (int i = 0; i < m_nparam; i++)
	{
		paramVoisin[i] = param[i];
	}

	for (int i = 0; i < m_npoint; i++)
	{
		int		x		= 0;
		double	dist1	= 1e+14;
		double	dist2	= 1e+14;

		for (int j = 0; j < m_npoint; j++)
		{
			if (j != i)
			{
				dist2 = 0.00;

				for (int k = 0; k < m_ndim; k++)
				{
					dist2 += abs(param[k * m_npoint + i] - param[k * m_npoint + j]);
				}

				if (dist2 < dist1)
				{
					x		= j;
					dist1	= dist2;
				}
			}
		}

		for (int k = 0; k < m_ndim; k++)
		{
			paramVoisin[k * m_npoint + i] += product * (param[k * m_npoint + i] - param[k * m_npoint + x]);

			if (paramVoisin[k * m_npoint + i] > 1.0)
			{
				paramVoisin[k * m_npoint + i] = 0.99;
			}

			if (paramVoisin[k * m_npoint + i] < 0.0)
			{
				paramVoisin[k * m_npoint + i] = 0.01;
			}
		}
	}
}


void Recuit::voisinRandom(int iter, double* param, double* paramVoisin)
{
	int r = rand() % m_nparam;

	for (int i = 0; i < m_nparam; i++)
	{
		paramVoisin[i] = param[i];
	}

	paramVoisin[r] = rand() / double(RAND_MAX);
}