/* Copyright (C) 2024, Murad Banaji
 *
 * This file is part of EPIcomp, for compartmental models in epidemiology
 *
 * EPIcomp is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2, 
 * or (at your option) any later version.
 *
 * EPIcomp is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with EPIcomp: see the file COPYING.  If not, write to 
 * the Free Software Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA. 

 */

#include "Epi.h"



//Ii+Si-->(ki)2Ii, infection process
//Ii+Sj->Ii+Ij, cross infection (occurs in compartment j)
//Ii-->(kr) Ri (death or recovery with immunity)
//We compare theoretical and simulated mean outbreak sizes in a
//complete symmetric network as we vary the leak
int main(int argc, char *argv[]){

  int mainargs=0;

  double tmax=500.0;//cut off by time tmax if extinction isn't achieved
  char *outfile;
  int sims_per_param;  
  double RepNo;//set R0, infer ki
  int Nc;//Total population in each compartment
  int ncomp;//number of compartments

  double leak;
  double eps;//This is the *total* out-infection fraction (so eps_{ii} = 1-eps)
  
  double compV;
  double ki;//infection rate constant
  double kr=1.0;//recovery rate constant
  double **epsM;

  int *S0, *R0;
  double *V;

  int i;
  int totparams=21;
 
  double det_sz;//deterministic expected epidemic size
  double prior=0.0;//no prior immunity
  FILE *fd;  
  
  static std::mt19937 rng(std::random_device{}()); // Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> runif(0.0, 1.0);


  while(optind < argc){
    switch (mainargs) {
    case 0:
      outfile=argv[optind];
      break;
    case 1:
      sims_per_param=atoi(argv[optind]);
      break;
    case 2:
      RepNo=atof(argv[optind]);
      break;
    case 3:
      Nc=atoi(argv[optind]);
      break;
    case 4:
      ncomp=atoi(argv[optind]);
      break;
    }
    mainargs++;
    optind++;
  }

  if(mainargs<5){
    fprintf(stderr, "Insufficient parameters. These should be (1) number of compartments, (2) compartment population, (3) R0, (4) simulations per parameter value, and (5) output file. EXITING.\n");exit(0);
  }


  fd=fopen(outfile,"w");
  compV=(double)Nc;//set population density to 1
  ki=RepNo*kr*compV/Nc;
  epsM=dmatrix(ncomp,ncomp);
  S0=(int*) malloc(sizeof(int) * (ncomp+1));
  R0=(int*) malloc(sizeof(int) * (ncomp+1));
  V=(double*) malloc(sizeof(double) * (ncomp+1));

  
  for(i=0;i<totparams;i++){
    leak = i*0.25;
    eps=leak/((double)Nc);
    sym_coupling(eps, epsM, ncomp);
    inittoconst(S0,ncomp,Nc);
    inittozero(R0,ncomp);
    inittoconst(V,ncomp,compV);
    det_sz=(1.0+boost::math::lambert_w0(-1.0*RepNo*exp(-RepNo))/RepNo);
    //a=pow((1.0-eps*(RepNo-1.0)/((double)(ncomp-1))),det_sz);
    //fprintf(stderr, "detsz=%.4f, a=%.4f\n",det_sz, a);
    fprintf(fd, "%.2f %.4f %.4f\n", leak, mean_outbreak_sym(leak, RepNo, Nc, ncomp, prior), mean_SIR_outbreak(rng, ncomp, V, S0, R0, epsM, ki, kr, tmax, sims_per_param, 0, NULL)/det_sz/(1.0-1.0/RepNo));
    //fprintf(stderr, "%.2f, %.4f, %.4f\n", leak, mean_outbreak_sym(leak, RepNo, Nc, ncomp));
  }

  free_dmatrix(epsM,ncomp,ncomp);
  fclose(fd);

  
  return 0;
}
