/* Copyright (C) 2024-6, Murad Banaji
 *
 * This file is part of EPIcode, for metapopulation models in epidemiology
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


// Variable graph topology, but identical compartments
// Ii+Sj->Ii+Ij (i,j=1,...,n) (infection, but no migration)
// Ii-->Ri (death or recovery with immunity)

// To compare primary and secondary outbreak sizes
int main(int argc, char *argv[]){

  int opt, mainargs=0;
  char *outfile, *outfile1; // for data and summary
  char *arg8, *arg9, *arg10, *arg11;

  double tmax=500.0;//cut off by time tmax if extinction isn't achieved

  int maxruns;//maximum total runs
  int epimax=20;//desired number of primary epidemics
  double threshold=0.1;//An epidemic if 10% or more of mean size
  int secondepi=0;//count secondary epidemics (mean)
  int largersecondepi=0;//count secondary epidemics larger than first
  int n_sim;//n_sim: simulations used to compute mean outbreak size
  double RepNo;//set R0, infer ki
  double leak;
  int Nc;//Total population in each compartment
  int gridtype;
  
  double compV;
  double eps;
  double kr=1.0;//recovery rate constant
  double ki;//infection rate constant

  //network type and coupling
  //gridtype=1:complete; 2:gamma, 3:ring small_world, 4:2d grid small world, 5:scale-free, 6:star-shaped

  int Nx,Ny,ncomp;
  double rewiring_prob;//in [0,1]: probability of rewiring for small-world networks
  //shape parameter for gamma distributed coupling strengths;
  //mean degree for random networks
  double param1;
  int clq_sz;//initial clique size for ring small-world (must be even)
  int scalefree_clq;//initial clique size for scale free networks
  //initial ring size for ring-based preferential attachement networks
  //symmetric or not for gamma coupling
  int param2=0;
  

  //Set this to one to choose introductions according to the susceptible population rather than the total population in each compartment
  int choose_by_S=0;
  int Iinit=1;//initial infections
  int compinit;
  double r1;

  //constout=0: constant total infectivity, variable in-infection,
  // and constant *total* out-infection;
  //constout=1: constant total infectivity, constant in-infection
  // and constant *total* out-infection;
  //constout=2: variable total infectivity, constant in-infection and
  //constant out-infection on each edge
  int constout=2;
  
  int totpop;

  double **epsM;
  int **epsI;
  double **NextGen;
  int **sparse=NULL;

  double mean_eps, epsvar;//for variable leak simulations

  int i,jrun,epicount=0;
  double t=0.0;

  int *S0, *I0, *R0, *S, *I, *R, *S1, *R1;
  double *V;

  FILE *fd;//for plotting
  FILE *fd1;
  FILE *fd2=fopen("SIRgeneral/outfiles/network_properties", "w");//for plotting purposes

  double sz_natural, sz_natural1, immune_frac;
  int numbins=14;
  int count[numbins];

  static std::mt19937 rng(std::random_device{}()); // Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> runif(0.0, 1.0);


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
      outfile=argv[optind];
      break;
    case 1:
      outfile1=argv[optind];
      break;
    case 2:
      maxruns=atoi(argv[optind]);
      break;
    case 3:
      epimax=atoi(argv[optind]);
      break;
    case 4:
      n_sim=atoi(argv[optind]);
      break;
    case 5:
      RepNo=atof(argv[optind]);
      break;
    case 6:
      leak=atof(argv[optind]);
      break;
    case 7:
      Nc=atoi(argv[optind]);
      break;
    case 8:
      gridtype=atoi(argv[optind]);
      break;
    case 9:
      arg8=argv[optind];
      break;
    case 10:
      arg9=argv[optind];
      break;
    case 11:
      arg10=argv[optind];
      break;
    case 12:
      arg11=argv[optind];
      break;
    }
    mainargs++;
    optind++;
  }
  fd=fopen(outfile, "w");//for plotting
  fd1=fopen(outfile1, "w");//summary

  
  if(mainargs<10){
    fprintf(stderr, "ERROR: expecting at least the followng arguments:\n\t\"outfile\" (output file for primary/secondary outbreak sizes)\n\t\"outfile1\" (output file with summary statistics)\n\t\"maxruns\" (total base runs)\n\t\"epimax\" (desired number of epidemics)\n\t\"n_sim\" (runs per mean value)\n\t\"RepNo\" (basic reproduction number)\n\t\"leak\" (fraction of out-infections X compartment population)\n\t\"Nc\" (compartment population)\n\t\"gridtype\" (see source file for available types)\n\tnumber of compartments\"Nc\", or \"Nx\" and \"Ny\" for a 2D grid\n\nEXITING\n");exit(0);
  }



  if(gridtype==4){//2d grid small world
    if(mainargs<12){
      fprintf(stderr, "for a 2d grid you must specify the dimensions \"Nx\" and \"Ny\", and the rewiring probability \"rewiring_prob\". EXITING\n");exit(0);
    }
    Nx=atoi(arg8);Ny=atoi(arg9);ncomp=Nx*Ny;rewiring_prob=atof(arg10);
  }
  else{
    ncomp=atoi(arg8);//number of compartments
    if(gridtype==2){//gamma
      if(mainargs<12){
	fprintf(stderr, "for gamma distributed weights, you must specify the shape parameter, and whether coupling is symmetric (0=not/1=symmetric). EXITING\n");exit(0);
      }
      param1=atof(arg9);
      param2=atoi(arg10);
    }
    else if(gridtype==3){//ring small world
      if(mainargs<12){
	fprintf(stderr, "for a (ring) small-world network you must specify:\n\tthe initial clique size \"clq_sz\" (must be even), and\n\tthe rewiring probability \"rewiring_prob\".\n\nEXITING\n");exit(0);
      }
      clq_sz=atoi(arg9);
      rewiring_prob=atof(arg10);

    }
    else if(gridtype==5){//scale free
      if(mainargs<11){
	fprintf(stderr, "for a scale free network, you must specify the initial clique size (a positive integer)\n\nEXITING\n");exit(0);
      }
      scalefree_clq=atoi(arg9);
      if(scalefree_clq<=0){
	fprintf(stderr, "for a scale-free network, the initial clique size must be a positive integer (currently %d)\n\nEXITING\n", scalefree_clq);exit(0);
      }
    }
    else if(gridtype==6){//random
      if(mainargs<11){
	fprintf(stderr, "for a random network, you must specify the mean degree.\nEXITING\n");exit(0);
      }
      param1=atof(arg9);
    }
    else if(gridtype==7){//ring scale-free
      if(mainargs<13){
	fprintf(stderr, "For a (ring) preferential attachment network you must specify:\n\tthe ring size \"init_ring\", \n\tthe initial clique size \"clq_sz\" (must be even), and\n\tnumber of new edges at each step.\n\nEXITING\n");exit(0);
      }
      param2=atoi(arg9);
      clq_sz=atoi(arg10);
      scalefree_clq=atoi(arg11);
    }
    else if(gridtype!=1 &&gridtype!=8){
      fprintf(stderr, "\"gridtype\" must be in the range 1-8. (Currently %d).\nEXITING.\n", gridtype);exit(0);
    }

  }

  
  compV=(double)Nc;//set population density to 1
  eps=leak/((double)Nc);
  ki=RepNo*kr*compV/Nc;
  totpop=ncomp*Nc;

  epsM=dmatrix(ncomp,ncomp); //The coupling matrix
  epsI=imatrix(ncomp,ncomp); //The coupling matrix structure
  NextGen=dmatrix(ncomp,ncomp); //The next generation matrix
  S0=(int*) malloc(sizeof(int) * (ncomp+1));
  I0=(int*) malloc(sizeof(int) * (ncomp+1));
  R0=(int*) malloc(sizeof(int) * (ncomp+1));
  S=(int*) malloc(sizeof(int) * (ncomp+1));
  I=(int*) malloc(sizeof(int) * (ncomp+1));
  R=(int*) malloc(sizeof(int) * (ncomp+1));
  S1=(int*) malloc(sizeof(int) * (ncomp+1));
  R1=(int*) malloc(sizeof(int) * (ncomp+1));
  V=(double*) malloc(sizeof(double) * (ncomp+1));
  

  if(n_sim<1 || n_sim>100000){
    fprintf(stderr, "numvax must be in the range 0-100000. EXITING.\n");exit(0);
  }


  inittoconst(V,ncomp,compV);//Volumes all equal

  //set the coupling matrix
  setgrid(rng, gridtype, epsM, epsI, ncomp, Nx, Ny, Nc, &sparse, RepNo, leak, eps, param1, rewiring_prob, param2, clq_sz, scalefree_clq, &mean_eps, &epsvar, constout, fd2);
  fclose(fd2);

  inittozero(count,numbins);
  
  //initial run to find expected initial outbreak size
  //set n_sim to 10000 for this initial run
  inittoconst(S1,ncomp,Nc);inittozero(R1,ncomp);
  sz_natural=mean_SIR_outbreak(rng, ncomp, V, S1, R1, epsM, ki, kr, tmax, 10000, choose_by_S, sparse);
  fprintf(stderr, "Mean outbreak: %.4f%%\n", sz_natural*100.0);

  for(jrun=0;jrun<maxruns;jrun++){
    fprintf(stderr, "%d/%d\n", jrun+1, maxruns);

    //"Iinit" infections in a randomly chosen compartment
    //randomness is important in the case of non-regular networks
    do{r1 = runif(rng);}while(r1==1.0);
    compinit=(int)(r1*ncomp);
    for(i=0;i<ncomp;i++){
      if(i==compinit){I0[i]=Iinit;}else{I0[i]=0;}
      S0[i]=Nc-I0[i];R0[i]=0;
    }

    //Run the epidemic simulation
    if((t=runSIRepi(rng, ncomp, V, I0, S0, R0, I, S, R, NULL, NULL, NULL, epsM, ki, kr, tmax, 0, 0, 0, 0, NULL, sparse))<0.0){
      fprintf(stderr,"Epidemic duration longer than tmax. EXITING.\n");exit(0);
    }
    

    immune_frac=((double)vsum(R,ncomp))/((double)totpop);
    if(immune_frac>threshold*sz_natural){//only significant epidemics
      epicount++;
      sz_natural1=mean_SIR_outbreak(rng, ncomp, V, S, R, epsM, ki, kr, tmax, n_sim, choose_by_S, sparse);
      fprintf(stderr, "primary: %.4f%%, secondary: %.4f%%\n", immune_frac*100, 100.0*sz_natural1);
      fprintf(fd, "%.4f %.4f\n", immune_frac*100.0, sz_natural1*100.0);
      if(sz_natural1>threshold*sz_natural){
	secondepi++;
	if(sz_natural1>immune_frac)
	  largersecondepi++;
      }
    }
    if(epicount>=epimax)
      break;
  }

  fprintf(fd1, "Herd immunity: %.4f%%, %.4f%%\n", 100-(double)secondepi/(double)epicount*100, (double)largersecondepi/(double)epicount*100);
  
  fclose(fd);
  fclose(fd1);
  free_dmatrix(NextGen,ncomp,ncomp);
  free_dmatrix(epsM,ncomp,ncomp);
  free_imatrix(epsI,ncomp,ncomp);
  free((char*)S0);free((char*)I0);free((char*)R0);free((char*)V);
  free((char*)S);free((char*)I);free((char*)R);
  free((char*)S1);free((char*)R1);
  if(sparse)
    free_imat(sparse,ncomp);

  
  return 0;
}
