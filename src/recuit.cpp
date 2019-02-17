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
	m_temp_init		= 10;
	m_lambda		= 0.5;
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
		voisin(s, sn);

		en = m_function->calcf(sn);

		// Metropolis
		if (en < e || (rand() / double(RAND_MAX)) < exp((e - en) / temp))
		{
			#pragma omp parallel for
			for (int i = 0; i < m_nparam; i++)
			{
				s[i] = sn[i];
			}

			e = en;
		}

		// Global Optimum
		if (e < m)
		{
			#pragma omp parallel for
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


void Recuit::voisinAll(double* param, double* paramVoisin)
{
	int*	itArray	= new int[m_npoint];
	double	product	= 1.0 / sqrt(m_npoint);

	for (int i = 0; i < m_npoint; i++)
	{
		itArray[i] = 0;
	}

	for (int i = 0; i < m_npoint; i++)
	{
		if (itArray[i] == 0)
		{
			int		x		= 0;
			double	dist1	= 1e+14;
			double	dist2	= 1e+15;
			double	dist3	= 0.0;

			for (int j = 0; j < m_npoint; j++)
			{
				if (j != i)
				{
					dist2 = 0.00;
	
					for (int k = 0; k < m_ndim; k++)
					{
						dist2 += (param[k * m_npoint + i] - param[k * m_npoint + j]) * (param[k * m_npoint + i] - param[k * m_npoint + j]);
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
				dist3 = param[k * m_npoint + i] - param[k * m_npoint + x];

				paramVoisin[k * m_npoint + i] = param[k * m_npoint + i] + product * dist3;
				paramVoisin[k * m_npoint + x] = param[k * m_npoint + x] - product * dist3;
	
				if (paramVoisin[k * m_npoint + i] > 1.0)
				{
					paramVoisin[k * m_npoint + i] = 1.0;
				}
				if (paramVoisin[k * m_npoint + i] < 0.0)
				{
					paramVoisin[k * m_npoint + i] = 0.0;
				}
				if (paramVoisin[k * m_npoint + x] > 1.0)
				{
					paramVoisin[k * m_npoint + x] = 1.0;
				}
				if (paramVoisin[k * m_npoint + x] < 0.0)
				{
					paramVoisin[k * m_npoint + x] = 0.0;
				}
			}

			itArray[i] = 1;
			itArray[x] = 1;
		}
	}
}


void Recuit::voisin(double* param, double* paramVoisin)
{
	double product = 1.0 / sqrt(m_npoint);

	#pragma omp parallel for
	for (int i = 0; i < m_npoint; i++)
	{
		int		x		= 0;
		double	dist1	= 1e+14;
		double	dist2	= 1e+15;

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


void Recuit::voisinRandom(double* param, double* paramVoisin)
{
	for (int i = 0; i < m_nparam; i++)
	{
		paramVoisin[i] = param[i];
	}

	int r = rand() % m_nparam;

	paramVoisin[r] = rand() / double(RAND_MAX);
}