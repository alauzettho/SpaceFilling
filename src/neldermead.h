////////////////////////////////////////////////////////////////////////////////////////
///
/// Created by Thomas Alauzet on November 14, 2018.
/// Copyright 2018. All rights reserved.
///
////////////////////////////////////////////////////////////////////////////////////////

#ifndef NELDERMEAD_H
#define NELDERMEAD_H


#include <vector>
#include <iostream>
#include <assert.h>
#include "class.h"
#include "parser.h"
#include "objectivefunction.h"

using namespace std;


class NelderMead
{
	int			m_ndim;
	int			m_npoint;
	int			m_nparam;
	int			m_nitermax;
	int			m_nfunevalmax;
	int			m_simplex_size;
	double		m_beta;
	double		m_tolf;
	double		m_tau1;
	double		m_tau2;
	double		m_alpha;
	double		m_delta;
	double		m_gamma;
	double*		m_tolx;
	Function*	m_function;

	public:

		NelderMead(Function* function, int dim, int npoint);

		~NelderMead();

		void minimize(double* param, const double* min_param, const double* max_param);

		bool testX(const vector<vector<double> >& x);

		bool testF(double* f);

		void insertion(vector<vector<double> >& x, double* f);

		void evalSimplex(vector<vector<double> >& x, double* f, int& nfuneval);

		void sortSimplex(vector<vector<double> >& x, double *f);

		void updateSimplex(const vector<double>& new_point, const double new_f, vector<vector<double> >& x, double* f);

		void shrinkSimplex(vector<vector<double> >& x);
};


#endif