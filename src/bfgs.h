////////////////////////////////////////////////////////////////////////////////////////
///
/// Created by Thomas Alauzet on November 14, 2018.
/// Copyright 2018. All rights reserved.
///
////////////////////////////////////////////////////////////////////////////////////////

#ifndef BFGS_H
#define BFGS_H


#include <string>
#include <vector>
#include "parser.h"
#include "objectivefunction.h"

using namespace std;


class Bfgs
{
	int			m_maxiter;
	int			m_maxfuneval;
	int			m_ndim;
	int			m_npoint;
	int			m_nparam;
	double		m_tolg;
	Function*	m_function;	

	public:

		Bfgs(Function* function, int ndim, int npoint);

		~Bfgs();

		void	minimize(double* param, const double* min_param, const double* max_param);

		void	getStep(int iter, double* param, const double* min_param, const double* max_param, double* d, double& alpha, double& f, int& nfuneval);
		
		double	scalarProduct(int size, double* x, double* y);

		double	normInf(int size, double* vector);
};


#endif