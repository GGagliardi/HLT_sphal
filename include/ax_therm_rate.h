#ifndef __ax_therm_rate__
#define __ax_therm_rate__

#include "numerics.h"
#include "highPrec.h"
#include "Spectral.h"
#include "stat.h"

using namespace std;

void get_axion_therm_rate(double m, double sigma, int T, int Nboots, const distr_t_list& C, string out_path, int INCLUDE_ERRORS);


#endif
