////////////////////////////////////////////////////////////////////////////////////////
///
/// Created by Thomas Alauzet on November 28, 2018.
/// Copyright 2018. All rights reserved.
///
////////////////////////////////////////////////////////////////////////////////////////

#ifndef RECUIT_H
#define RECUIT_H

#include <omp.h>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include "parser.h"
#include "objectivefunction.h"

using namespace std;


class Recuit
{
	int			m_maxiter;
	int			m_maxfuneval;
	int			m_ndim;
	int			m_npoint;
	int			m_nparam;
	double		m_lambda;
	double		m_temp_init;
	Function*	m_function;	

	public:

		Recuit(Function* function, int ndim, int npoint);

		~Recuit();

		void	minimize(double* param, const double* min_param, const double* max_param);

		void	voisin(double* param, double* paramVoisin);

		void	voisinAll(double* param, double* paramVoisin);

		void	voisinRandom(double* param, double* paramVoisin);
};


#endif
