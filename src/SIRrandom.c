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


//random network. To explore how increasing the mean degree increases the expected outbreak size
//Ii+Si-->(ki)2Ii, infection process
//Ii+Sj->Ii+Ij, cross infection (occurs in compartment j)
//Ii-->(kr) Ri (death or recovery with immunity)
int main(int argc, char *argv[]){

  int opt, mainargs=0;
  double tmax=500.0;//cut off by time tmax if extinction isn't achieved

  int grids_per_param;//number of grids
  int sims_per_grid;//number of experiments per grid
  double RepNo;//set R0, infer ki
  double leak;
  int Nc;//Total population in each compartment
  int ncomp;//number of compartments
  
  double compV;
  double eps;
  double kr=1.0;//recovery rate constant
  double ki;//infection rate constant

  //constout=0: constant total infectivity, variable in-infection,
  // and constant *total* out-infection;
  //constout=1: constant total infectivity, constant in-infection
  // and constant *total* out-infection;
  //constout=2: variable total infectivity, constant in-infection and
  //constant out-infection on each edge
  //constout=3: variable total infectivity, constant in-infection and
  //constant susceptibility
  int constout=2;
  double meandeg;//shape parameter for gamma distributed coupling strengths

  double **epsM; //The coupling matrix
  int **epsI; //The coupling matrix structure
  double mean_eps;//for variable leak simulations
  double epsvar;

  int i,j;

  int *S0, *R0;
  double *V;

  FILE *fd=fopen("SIRgeneral/outfiles/random_meandeg", "w");
  FILE *fd1=fopen("SIRgeneral/outfiles/random_meandeg_mean", "w");
  FILE *fd2=fopen("SIRgeneral/outfiles/network_properties", "w");

  static std::mt19937 rng(std::random_device{}()); // Standard mersenne_twister_engine seeded with rd()
  
  int **sparse=NULL;
  double sz_natural, sz_natural_tot;

  //Set this to one to choose introductions according to the susceptible population rather than the total population in each compartment
  int choose_by_S=0;

  //reset constout if necessary
  while ((opt = getopt(argc, argv, "c:")) != -1) {
    switch (opt) {
    case 'c': //reset constout
      constout = atoi(optarg);
      break;
    default: /* '?' */
      fprintf(stderr, "Option not recognised");
      exit(EXIT_FAILURE);
    }
  }

  
  while(optind < argc){
    switch (mainargs) {
    case 0:
      grids_per_param=atoi(argv[optind]);
      break;
    case 1:
      sims_per_grid=atoi(argv[optind]);
      break;
    case 2:
      RepNo=atof(argv[optind]);
      break;
    case 3:
      leak=atof(argv[optind]);
      break;
    case 4:
      Nc=atoi(argv[optind]);
      break;
    case 5:
      ncomp=atoi(argv[optind]);
      break;
    }
    mainargs++;
    optind++;
  }

  if(mainargs<6){
    fprintf(stderr, "Insufficient parameters. These should be:\n\t(1) grids per parameter value\n\t(2) simulations per grid\n\t(3) basic reproduction number (R0)\n\t(4) \"leak\" (fraction of out-infections X compartment population)\n\t(5) \"Nc\" (compartment population) and\n\t(6) number of compartments.\nEXITING.\n");exit(0);
  }
  
  compV=(double)Nc;//set population density to 1
  eps=leak/((double)Nc);//strength of total cross infection
  ki=RepNo*kr*compV/(double)Nc;//infection rate constant
  epsM=dmatrix(ncomp,ncomp); //The coupling matrix
  epsI=imatrix(ncomp,ncomp); //The coupling matrix structure
  S0=(int*) malloc(sizeof(int) * (ncomp+1));
  R0=(int*) malloc(sizeof(int) * (ncomp+1));
  V=(double*) malloc(sizeof(double) * (ncomp+1));
   
  //Initialisation
  inittoconst(V,ncomp,compV);inittoconst(S0,ncomp,Nc);inittozero(R0,ncomp);

  //varying the mean degree
  for(i=0;i<=10;i++){
    meandeg=4+i;
    //set the coupling matrix (gridtype=6)
    sz_natural_tot=0.0;
    for(j=0;j<grids_per_param;j++){
      fprintf(stderr, "%d/%d", j+1,grids_per_param);
      setgrid(rng, 6, epsM, epsI, ncomp, 0, 0, Nc, &sparse, RepNo, leak, eps, meandeg, 0, 0, 0, 0, &mean_eps, &epsvar, constout, fd2);
      sz_natural=mean_SIR_outbreak(rng, ncomp, V, S0, R0, epsM, ki, kr, tmax, sims_per_grid, choose_by_S, sparse);
      sz_natural_tot+=sz_natural;
      fprintf(fd, "%.4f %.4e %.4e %.4f\n", meandeg, mean_eps/(double)(ncomp-1), sqrt(epsvar), sz_natural*100.0);
      //fprintf(stderr, "mean degree=%.4f, mean (SD) eps=%.4e (%.4e), mean epi size=%.4f%%\n", meandeg, mean_eps/(double)(ncomp-1), sqrt(epsvar), sz_natural*100.0);
      
    }
    sz_natural_tot/=(double)grids_per_param;

    fprintf(fd1, "%.4f %.4f\n", meandeg, sz_natural_tot*100.0);
    fprintf(stderr, "mean degree=%.4f, mean epi size=%.4f%%\n", meandeg, sz_natural_tot*100.0);
      

  }

  fclose(fd);
  fclose(fd1);
  fclose(fd2);

  free_dmatrix(epsM,ncomp,ncomp);
  free_imatrix(epsI,ncomp,ncomp);
  free((char*)S0);free((char*)R0);free((char*)V);
  if(sparse)
    free_imat(sparse,ncomp);

  
  return 0;
}
