/* Copyright (C) 2024-5, Murad Banaji
 *
 * This file is part of EPIcode, for compartmental models in epidemiology
 *
 * EPIcode is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2, 
 * or (at your option) any later version.
 *
 * EPIcode is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with EPIcode: see the file COPYING.  If not, write to 
 * the Free Software Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA. 

 */

#include "Epi.h"

// Ii+Sj->Ii+Ij (i,j=1,...,n) (infection, but no migration)
// Ii-->Ri (death or recovery with immunity)

// To explore how increasing the variance of the coupling weights
// reduces the expected outbreak size in Gamma networks

int main(int argc, char *argv[]){

  int opt, mainargs=0;
  double tmax=500.0;//cut off by time tmax if extinction isn't achieved

  int grids_per_param=100;//number of grids
  int sims_per_grid=1000;//number of experiments per grid
  double RepNo=2.0;//set R0, infer ki
  double leak=2.5;
  int Nc=50;//Total population in each compartment
  int ncomp=30;//number of compartments
  double shp;//shape parameter for gamma distributed coupling strengths

  //constout=0: constant total infectivity, variable in-infection,
  // and constant *total* out-infection;
  //constout=1: constant total infectivity, constant in-infection
  // and constant *total* out-infection;
  //constout=2: variable total infectivity, constant in-infection and
  //constant out-infection on each edge
  //constout=3: variable total infectivity, constant in-infection and
  //constant susceptibility
  int constout=2;
  double compV;//set population density to 1
  double eps;//strength of total cross infection

  double kr=1.0;//recovery rate constant
  double ki;//infection rate constant

  double **epsM; //The coupling matrix
  int **epsI; //The coupling matrix structure
  double mean_eps;//for variable leak simulations
  double epsvar;
  int sym=0;//symmetric coupling or not?

  int i,j;

  int *S0, *R0;
  double *V;

  FILE *fd=fopen("SIRgeneral/outfiles/gamma_effect", "w");
  FILE *fd1=fopen("SIRgeneral/outfiles/gamma_effect_mean", "w");
  FILE *fd2=fopen("SIRgeneral/outfiles/network_properties", "w");

  int **sparse=NULL;
  double sz_natural, sz_natural_tot;

  //Set this to one to choose introductions according to the susceptible population rather than the total population in each compartment
  int choose_by_S=0;

  static std::mt19937 rng(std::random_device{}()); // Standard mersenne_twister_engine seeded with rd()

  /////////////////
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
    case 6:
      sym=atoi(argv[optind]);
      break;
    }
    mainargs++;
    optind++;
  }

  if(mainargs<7){
    fprintf(stderr, "Insufficient parameters. These should be:\n\t(1) grids per parameter value\n\t(2) simulations per grid\n\t(3) basic reproduction number (R0)\n\t(4) \"leak\" (fraction of out-infections X compartment population)\n\t(5) \"Nc\" (compartment population)\n\t(6) number of compartments\n\t(7) symmetric coupling?\nEXITING.\n");exit(0);
  }

  
  compV=(double)Nc;//set population density to 1
  eps=leak/((double)Nc);
  ki=RepNo*kr*compV/(double)Nc;
  epsM=dmatrix(ncomp,ncomp);
  epsI=imatrix(ncomp,ncomp);
  S0=(int*) malloc(sizeof(int) * (ncomp+1));
  R0=(int*) malloc(sizeof(int) * (ncomp+1));
  V=(double*) malloc(sizeof(double) * (ncomp+1));

  //Initialisation
  inittoconst(V,ncomp,compV);inittoconst(S0,ncomp,Nc);inittozero(R0,ncomp);

  //varying the shape parameter
  for(i=0;i<=10;i++){
    shp=pow(10.0, -(double)i/5.0);
    //set the coupling matrix (gridtype=2)
    sz_natural_tot=0.0;
    for(j=0;j<grids_per_param;j++){
      fprintf(stderr, "%d/%d", j+1,grids_per_param);
      setgrid(rng, 2, epsM, epsI, ncomp, 0, 0, Nc, &sparse, RepNo, leak, eps, shp, 0, sym, 0, 0, &mean_eps, &epsvar, constout, fd2);
      sz_natural=mean_SIR_outbreak(rng, ncomp, V, S0, R0, epsM, ki, kr, tmax, sims_per_grid, choose_by_S, sparse);
      sz_natural_tot+=sz_natural;
      fprintf(fd, "%.4f %.4e %.4e %.4f\n", shp, mean_eps/(double)(ncomp-1), sqrt(epsvar), sz_natural*100.0);
      fprintf(stderr, "shape=%.4f, mean (SD) eps=%.4e (%.4e), mean epi size=%.4f%%\n", shp, mean_eps/(double)(ncomp-1), sqrt(epsvar), sz_natural*100.0);
      
    }
    sz_natural_tot/=(double)grids_per_param;

    fprintf(fd1, "%.4f %.4f\n", shp, sz_natural_tot*100.0);
    fprintf(stderr, "shape=%.4f, mean epi size=%.4f%%\n", shp, sz_natural_tot*100.0);
      

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
