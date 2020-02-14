#include <iostream>
#include <stdio.h>
#include <math.h>
#include "camd.h"

using namespace std;

extern "C" void reordering(int * nin, int * Mp, int * Mi, int * reo, int * cind)
{
  int n = nin[0];
  camd_order(n,Mp,Mi,reo,(double *) NULL, (double *) NULL, cind);
}

				
