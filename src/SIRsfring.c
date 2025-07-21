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

// We compute mean outbreak sizes in ring-preferential networks
// as we increase the initial ring size

int main(int argc, char *argv[]){


  int opt, mainargs=0;
  double tmax=500.0;//cut off by time tmax if extinction isn't achieved
  
  int grids_per_param;
  int sims_per_grid;
  double RepNo;
  double leak;
  int Nc;
  int ncomp;
  int clq_sz=4;
  int new_edges;
  
  int init_ring;
  double kr=1.0;//recovery rate constant
  double ki;//infection rate constant

  //constout=0: constant total infectivity, variable in-infection,
  // and constant *total* out-infection;
  //constout=1: constant total infectivity, constant in-infection
  // and constant *total* out-infection;
  //constout=2: variable total infectivity, constant in-infection and
  //constant out-infection on each edge
  int constout=2;

  double compV;
  double eps;

  double **epsM;
  int **epsI;
  int **sparse=NULL;
  double meaneps, epsvar;

  int i, j;
  int *S0, *R0;
  double *V;

 
  double epi_sz, epi_sz_tot;
  FILE *fd=fopen("SIRgeneral/outfiles/network_properties", "w");
  FILE *fd1=fopen("SIRgeneral/outfiles/outsize_ringsf_rewire","w");
  FILE *fd2=fopen("SIRgeneral/outfiles/outsize_ringsf_rewire_mean","w");

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
      ncomp=atoi(argv[optind]);
      break;
    case 6:
      clq_sz=atoi(argv[optind]);
      break;
    case 7:
      new_edges=atoi(argv[optind]);
      break;
    }
    mainargs++;
    optind++;
  }


  if(mainargs<8){
    fprintf(stderr, "Insufficient parameters. These should be:\n\t(1) grids per parameter value\n\t(2) simulations per grid\n\t(3) basic reproduction number (R0)\n\t(4) \"leak\" (fraction of out-infections X compartment population)\n\t(5) \"Nc\" (compartment population),\n\t(6) number of compartments,\n\t(7) clique size, and\n\t(7) new edges per step.\n\t(initial clique size is always 2)\n\t(initial ring size is increased from 3)\nEXITING.\n");exit(0);
  }


  compV=(double)Nc;//set population density to 1
  eps=leak/((double)Nc);  
  ki=RepNo*kr*compV/(double)Nc;//infection rate constant
  epsM=dmatrix(ncomp,ncomp); //The coupling matrix
  epsI=imatrix(ncomp,ncomp);
  S0=(int*) malloc(sizeof(int) * (ncomp+1));
  R0=(int*) malloc(sizeof(int) * (ncomp+1));
  V=(double*) malloc(sizeof(double) * (ncomp+1));
  
  
  for(i=0;i<11;i++){
    init_ring=6+2*i;
    inittoconst(S0,ncomp,Nc);inittozero(R0,ncomp);inittoconst(V,ncomp,compV);
    epi_sz_tot=0.0;
    for(j=0;j<grids_per_param;j++){
      setgrid(rng, 7, epsM, epsI, ncomp, 0, 0, Nc, &sparse, RepNo, leak, eps, 0, 0, init_ring, clq_sz, new_edges, &meaneps, &epsvar, constout, fd);
      epi_sz=mean_SIR_outbreak(rng, ncomp, V, S0, R0, epsM, ki, kr, tmax, sims_per_grid, 0, sparse);
      epi_sz_tot+=epi_sz;
      fprintf(fd1, "%d %.4f\n", init_ring, epi_sz*100.0);
      //fprintf(stderr, "%.2f: Mean epi sz = %.4f%%\n", rewiring_prob, epi_sz*100);
    }
    epi_sz_tot/=(double)grids_per_param;
    fprintf(stderr, "initial ring: %d, mean epi sz: %.4f%%\n", init_ring, epi_sz_tot*100.0);
    fprintf(fd2, "%d %.4f\n", init_ring, epi_sz_tot*100.0);
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
