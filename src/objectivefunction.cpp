////////////////////////////////////////////////////////////////////////////////////////
///
///	Created by Thomas Alauzet on November 13, 2018.
/// Copyright 2018. All rights reserved.
/// Able to add pragma OMP.
///
////////////////////////////////////////////////////////////////////////////////////////

#include "objectivefunction.h"

using namespace std;


Function::Function(estimation_method estimation_method, int ndim, int npoint) : m_estimation_method (estimation_method), m_ndim (ndim), m_npoint (npoint)
{
	m_q				= 0.5;
	m_eps_diff		= 0.00001;
	m_nparam 		= m_npoint * m_ndim;
}


Function::~Function(){}


double Function::squareNorm(int size, double* vector)
{
	double norm = 0.0;

	for (int i = 0; i < size; i++)
	{
		norm += vector[i] * vector[i];
	}

	return(norm);
}


void Function::calcGradient(double* param, const double f, double* g, int& nfuneval)
{
	double	fp;
	double	fm;
	double* new_param = new double[m_nparam];

	for (int i = 0; i < m_nparam; i++)
	{
		new_param[i] = param[i];
	}

	for (int i = 0; i < m_nparam; i++)
	{
		new_param[i] += m_eps_diff;

		fp = calcf(new_param);
		nfuneval++;

		new_param[i] -= 2 * m_eps_diff;

		fm = calcf(new_param);
		nfuneval++;

		g[i] = (fp - fm) / (2 * m_eps_diff);

		new_param[i] += m_eps_diff;

		// cout << "g[" << i << "] = (" << fp << " - " << fm << ") / " << m_eps_diff << " = " << g[i] << endl;
	}
}


double Function::calcf(double* param)
{
	double fvalue = 1e+15;

	if (m_estimation_method == P_MONTE_CARLO)
	{
		fvalue = monteCarlo(param);
	}
	else if (m_estimation_method == P_NEAREST_NEIGHBOR)
	{
		fvalue = nearestNeighbor(param);
	}
	else if (m_estimation_method == P_MINIMAL_SPANNING_TREE)
	{
		fvalue = minimalSpanningTree(param);
	}
	else if (m_estimation_method == P_ROSENBROCK)
	{
		fvalue = rosenbrock(param);
	}
	else if (m_estimation_method == P_MAXIMIN)
	{
		fvalue = maximin(param);
	}
	else if (m_estimation_method == P_AE)
	{
		fvalue = aeMaximizer(param);
	}
	else if (m_estimation_method == P_KL)
	{
		fvalue = shannon(param);
	}


	// Checking if parameters are within boudaries
	for (int i = 0; i < m_nparam; i++)
	{
		if (param[i] > 1.0 || param[i] < 0.0) fvalue = 1e+15;
	}

	return(fvalue);
}


double Function::rosenbrock(double* param)
{
	double fvalue = 0.0;

	for (int i = 0; i < m_nparam - 1; i++)
	{
		fvalue += 100 * (param[i + 1] - param[i] * param[i]) * (param[i + 1] - param[i] * param[i]) + (param[i] - 1) * (param[i] - 1);
	}

	return(fvalue);
}


double Function::monteCarlo(double* param)
{
	double fvalue = 0.0;
	double square = m_ndim / 12.0;
	double hquare = 1.0 / (12.0 * pow(m_npoint, 2.0 / (m_ndim + 4)));
	double coeff1 = - 1.0 / (2.0 * square * hquare);
	double coeff2 = 0.0;
	double* diffv = new double[m_ndim];


	for (int i = 0; i < m_npoint; i++)
	{
		coeff2 = 0.0;

		for (int j = i; j < m_npoint; j++)
		{
			for (int k = 0; k < m_ndim; k++)
			{
				diffv[k] = param[k * m_npoint + i] - param[k * m_npoint + j];
			}

			coeff2 += exp(coeff1 * squareNorm(m_ndim, diffv));
		}

		fvalue += pow(coeff2, m_q - 1);
	}

	delete[] diffv;

	return(-fvalue);
}


double Function::nearestNeighbor(double* param)
{
	double	dist	= 0.0;
	double	norm	= 0.0;
	double	fvalue	= 0.0;
	double	coeff1	= m_ndim * (1 - m_q);
	double* diffv	= new double[m_ndim];


	for (int i = 0; i < m_npoint; i++)
	{
		dist = 1e+15;
		for (int j = 0; j < m_npoint; j++)
		{
			if (i != j)
			{
				for (int k = 0; k < m_ndim; k++)
				{
					diffv[k] = param[k * m_npoint + i] - param[k * m_npoint + j];
				}
				
				norm = sqrt(squareNorm(m_ndim, diffv));

				if (norm < dist)
				{
					dist = norm;
				}
			}
		}

		fvalue += pow(dist, coeff1);
	}

	delete[] diffv;

	return(-fvalue);
}


double Function::minimalSpanningTree(double* param)
{
	int		V		= m_npoint;
	int		E		= V * (V - 1) / 2;
	double	dist	= 0.0;
	double	fvalue	= 0.0;
	double* diffv	= new double[m_ndim];


	Graph g(V, E);

	// Create Complete Graph
	for (int i = 0; i < V; i++)
	{
		for (int j = i + 1; j < V; j++)
		{
			for (int k = 0; k < m_ndim; k++)
			{
				diffv[k] = param[k * m_npoint + i] - param[k * m_npoint + j];
			}

			dist = sqrt(squareNorm(m_ndim, diffv));

			g.addEdge(i, j, dist);
		}
	}

	fvalue = g.kruskalMST(m_ndim, m_q);

	delete[] diffv;

	return(-fvalue);
}


void Graph::addEdge(double u, double v, double w)
{
	edges.push_back({w, {u, v}});
}


double Graph::kruskalMST(int ndim, double q)
{
	double mst_wt	= 0.0;
	double gamma	= ndim * (1 - q);


	sort(edges.begin(), edges.end());

	DisjointSets ds(V);

	vector< pair<double, iPair> >::iterator it;


	for (it = edges.begin(); it != edges.end(); it++)
	{
		int u = it->second.first;
		int v = it->second.second;

		int set_u = ds.find(u);
		int set_v = ds.find(v);

		if (set_u != set_v)
		{
			// cout << u << " - " << v << endl;

			mst_wt += pow(it->first, gamma);

			ds.merge(set_u, set_v);
		}
	}

	return(mst_wt);
}


double Function::maximin(double* param)
{
	double	dist	= 0.0;
	double	norm	= 0.0;
	double	fvalue	= 1e+15;
	double* diffv	= new double[m_ndim];


	for (int i = 0; i < m_npoint; i++)
	{
		dist = 1e+15;
		for (int j = 0; j < m_npoint; j++)
		{
			if (i != j)
			{
				for (int k = 0; k < m_ndim; k++)
				{
					diffv[k] = param[k * m_npoint + i] - param[k * m_npoint + j];
				}
				
				norm = sqrt(squareNorm(m_ndim, diffv));

				if (norm < dist)
				{
					dist = norm;
				}
			}
		}

		if (dist < fvalue)
		{
			fvalue = dist;
		}
	}

	delete[] diffv;

	return(-fvalue);
}


double Function::aeMaximizer(double* param)
{
	double	dist	= 0.0;
	double	fvalue	= 0.0;
	double* diffv	= new double[m_ndim];


	for (int i = 0; i < m_npoint - 1; i++)
	{
		for (int j = i + 1; j < m_npoint; j++)
		{
			for (int k = 0; k < m_ndim; k++)
			{
				diffv[k] = param[k * m_npoint + i] - param[k * m_npoint + j];
			}

			dist = squareNorm(m_ndim, diffv);

			fvalue += 1 / dist;
		}
	}

	delete[] diffv;

	return(fvalue);
}


double Function::shannon(double* param)
{
	double fvalue = 0.0;
	double square = m_ndim / 12.0;
	double hquare = 1.0 / (12.0 * pow(m_npoint, 2.0 / (m_ndim + 4)));
	double coeff1 = - 1.0 / (2.0 * square * hquare);
	double coeff2 = 0.0;
	double* diffv = new double[m_ndim];


	for (int i = 0; i < m_npoint; i++)
	{
		coeff2 = 0.0;

		for (int j = i; j < m_npoint; j++)
		{
			for (int k = 0; k < m_ndim; k++)
			{
				diffv[k] = param[k * m_npoint + i] - param[k * m_npoint + j];
			}

			coeff2 += exp(coeff1 * squareNorm(m_ndim, diffv));
		}

		fvalue += log(coeff2);
	}

	delete[] diffv;

	return(fvalue);
}