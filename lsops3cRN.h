#ifndef __LSOPS3CRN_H
#define __LSOPS3CRN_H
#include "llistRN.h"
#include "mem_seg_levelset030414sf.h"
#include <math.h>
void ls_iteration(double *F, double *phi, double* label, long* dims, 
                  LL* Lz, LL* Ln1, LL* Lp1, LL *Ln2, LL *Lp2, 
                  LL *Lin2out, LL* Lout2in);
void ls_iteration_RN(double *phi, double* label, long* dims, 
                  LL* Lz, LL* Ln1, LL* Lp1, LL *Ln2, LL *Lp2, 
		     LL *Lin2out, LL* Lout2in);


void ls_mask2phi3c(double* mask, double* phi, double* label, long* dims, 
                   LL* Lz, LL* Ln1, LL* Ln2, LL* Lp1, LL* Lp2);

double ls_min_hood_onlevel(int idx, long x, long y, long z, long *dims, 
                           double *phi,double *label, double level);

double ls_max_hood_onlevel(int idx, long x, long y, long z, long *dims, 
                           double *phi,double *label, double level);
#endif
