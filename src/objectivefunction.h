////////////////////////////////////////////////////////////////////////////////////////
///
///	Created by Thomas Alauzet on November 13, 2018.
/// Copyright 2018. All rights reserved.
///
////////////////////////////////////////////////////////////////////////////////////////

#ifndef OBJECTIVEFUNCTION_H
#define OBJECTIVEFUNCTION_H

#include <iostream>
#include <bits/stdc++.h>
#include "math.h"
#include "class.h"


typedef pair<int, int> iPair;


class Function
{
	int					m_ndim;
	int					m_nparam;
	int					m_npoint;
	double				m_q;
	double				m_eps_diff;
	estimation_method	m_estimation_method;


	public:
		Function(estimation_method estimation_method, int ndim, int npoint);

		~Function();

		void	calcGradient(double* param, const double f, double* g, int& nfuneval);

		double	squareNorm(int size, double* vector);

		double	calcf(double* param);

		double	rosenbrock(double* param);

		double	monteCarlo(double* param);

		double	nearestNeighbor(double* param);

		double	minimalSpanningTree(double* param);
};


struct Graph
{
	int V;
	int E;

	vector<pair<double, iPair> > edges;

	Graph(int V, int E)
	{
		this->V = V;
		this->E = E;
	}

	void addEdge(double u, double v, double w);

	double kruskalMST(int ndim, double q);
};


struct DisjointSets
{
	int *parent;
	int *rnk;
	int n;

	DisjointSets(int n)
	{
		this->n	= n;
		parent	= new int[n + 1];
		rnk		= new int[n + 1];

		for (int i = 0; i <= n; i++)
		{
			rnk[i] = 0;
			parent[i] = i;
		}
	}

	int find(int u)
	{
		if (u != parent[u])
		{
			parent[u] = find(parent[u]); 
		}
		return parent[u];
	}

	void merge(int x, int y)
	{
		x = find(x), y = find(y);
  
		if (rnk[x] > rnk[y])
		{
			parent[y] = x;
		}
		else
		{
			parent[x] = y;
		}
		if (rnk[x] == rnk[y])
		{
			rnk[y]++;
		}
	}
};


#endif 
