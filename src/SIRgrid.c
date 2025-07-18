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
//We compute mean outbreak sizes in a 2D lattice as we increase the rewiring probability
int main(int argc, char *argv[]){

  int opt, mainargs=0;
  double tmax=500.0;//cut off by time tmax if extinction isn't achieved
  
  int grids_per_param;
  int sims_per_grid;
  double RepNo;
  double leak;
  int Nc;
  int Nx, Ny;

  int ncomp;
  double compV;//set population density to 1
  double eps;
  double ki;//infection rate constant
  
  double kr=1.0;//recovery rate constant


  double **epsM;
  int **epsI;
  int **sparse=NULL;
  double meaneps, epsvar;

  int i, j;
  double rewiring_prob;

  int *S0, *R0;
  double *V;

  //constout=0: constant total infectivity, variable in-infection,
  // and constant *total* out-infection;
  //constout=1: constant total infectivity, constant in-infection
  // and constant *total* out-infection;
  //constout=2: variable total infectivity, constant in-infection and
  //constant out-infection on each edge
  int constout=2;

 
  double epi_sz, epi_sz_tot;
  FILE *fd=fopen("SIRgeneral/outfiles/network_properties", "w");
  FILE *fd1=fopen("SIRgeneral/outfiles/grid_rewire","w");
  FILE *fd2=fopen("SIRgeneral/outfiles/grid_rewire_mean","w");

  static std::mt19937 rng(std::random_device{}()); // Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> runif(0.0, 1.0);

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
      Nx=atoi(argv[optind]);
      break;
    case 6:
      Ny=atoi(argv[optind]);
      break;
    }
    mainargs++;
    optind++;
  }

  if(mainargs<7){
    fprintf(stderr, "Insufficient parameters. These should be:\n\t(1) grids per parameter value\n\t(2) simulations per grid\n\t(3) basic reproduction number (R0)\n\t(4) \"leak\" (fraction of out-infections X compartment population)\n\t(5) \"Nc\" (compartment population)\n\t(6) \"Nx\" the x-dimension of the grid, and\n\t(7) \"Ny\", the y dimension of the grid\nEXITING.\n");exit(0);
  }

  ncomp=Nx*Ny;
  compV=(double)Nc;//set population density to 1
  eps=leak/((double)Nc);
  ki=RepNo*kr*compV/(double)Nc;
  epsM=dmatrix(ncomp,ncomp); //The coupling matrix
  epsI=imatrix(ncomp,ncomp);
  S0=(int*) malloc(sizeof(int) * (ncomp+1));
  R0=(int*) malloc(sizeof(int) * (ncomp+1));
  V=(double*) malloc(sizeof(double) * (ncomp+1));

  
  for(i=0;i<21;i+=2){
    rewiring_prob = (double)i*0.01;
    inittoconst(S0,ncomp,Nc);inittozero(R0,ncomp);inittoconst(V,ncomp,compV);
    epi_sz_tot=0.0;
    for(j=0;j<grids_per_param;j++){
      setgrid(rng, 4, epsM, epsI, ncomp, Nx, Ny, Nc, &sparse, RepNo, leak, eps, 0, rewiring_prob, 0, 0, 0, &meaneps, &epsvar, constout, fd);
      epi_sz=mean_SIR_outbreak(rng, ncomp, V, S0, R0, epsM, ki, kr, tmax, sims_per_grid, 0, sparse);
      epi_sz_tot+=epi_sz;
      fprintf(fd1, "%.2f %.4f\n", rewiring_prob, epi_sz*100.0);
      //fprintf(stderr, "%.2f: Mean epi sz = %.4f%%\n", rewiring_prob, epi_sz*100);
    }
    epi_sz_tot/=(double)grids_per_param;
    fprintf(stderr, "rewiring: %.2f%%, mean epi sz: %.4f%%\n", rewiring_prob*100.0, epi_sz_tot*100.0);
    fprintf(fd2, "%.2f %.4f\n", rewiring_prob, epi_sz_tot*100.0);
    //fprintf(stderr, "%.2f, %.4f, %.4f\n", leak, mean_outbreak_sym(leak, RepNo, Nc, ncomp));
  }

  free_dmatrix(epsM,ncomp,ncomp);
  free_imatrix(epsI,ncomp,ncomp);
  free((char*)S0);free((char*)R0);free((char*)V);
  if(sparse)
    free_imat(sparse,ncomp);
  fclose(fd);
  fclose(fd1);
  fclose(fd2);
  
  return 0;
}
