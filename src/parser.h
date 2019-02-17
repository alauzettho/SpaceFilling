////////////////////////////////////////////////////////////////////////////////////////
///
///	Created by Thomas Alauzet on November 12, 2018.
/// Copyright 2018. All rights reserved.
///
////////////////////////////////////////////////////////////////////////////////////////

#ifndef PARSER_H
#define PARSER_H

#include <ctime>
#include <string>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include "class.h"


int		getDim();
int		getPoint();
void	initializeCenter(int ndim, int npoint, double* param);
void	printResults(int ndim, int npoint, double* param, double f);
void	parse(int ndim, int npoint, double* a, double* b, optim_method& optim_method, estimation_method& estimation_method, int argc, char *argv[]);
string	intToString(int number);
string	doubleToString(double number);


#endif 