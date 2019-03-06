////////////////////////////////////////////////////////////////////////////////////////
///
///	Created by Thomas Alauzet on November 14, 2018.
/// Copyright 2018. All rights reserved.
///
////////////////////////////////////////////////////////////////////////////////////////

#ifndef CLASS_H
#define CLASS_H


using namespace std;


enum optim_method
{
	P_NELDER_MEAD,
	P_BFGS,
	P_RECUIT_SIMULE
};


enum estimation_method
{
	P_MONTE_CARLO,
	P_NEAREST_NEIGHBOR,
	P_MINIMAL_SPANNING_TREE,
	P_MAXIMIN,
	P_AE,
	P_KL,
	P_ROSENBROCK
};


enum nelder_mead_test
{
	REFLECTION,
	EXPANSION,
	CONTRACTION_OUT,
	CONTRACTION_IN,
	ERROR = -1
};


#endif