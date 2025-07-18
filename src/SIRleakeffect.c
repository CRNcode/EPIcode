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


//Variable graph topology, but identical compartments
//Ii+Si-->(ki)2Ii, infection process
//Ii+Sj->Ii+Ij, cross infection (occurs in compartment j)
//Ii-->(kr) Ri (death or recovery with immunity)
int main(int argc, char *argv[]){

  int opt, mainargs=0;
  char *outfile;
  int sims_per_param;//number of simulations per leak value
  double RepNo;//set R0, infer ki
  double leakmin, leakmax; //min and max leak
  int Nc;//Total population in each compartment
  int gridtype;//1:complete; 2:gamma, 3:ring small_world, 4:2d grid small world, 5:scale-free, 6:star-shaped
  char *arg8, *arg9, *arg10, *arg11;

  double tmax=500.0;//cut off by time tmax if extinction isn't achieved
  double compV;
  double kr=1.0;//recovery rate constant
  double ki;//infection rate constant
  double leak, eps;
  int Nx,Ny;
  int ncomp;//number of compartments
  //shape parameter for gamma distributed coupling strengths;
  //mean degree for random networks
  double param1;
  double rewiring_prob;//in [0,1]: probability of rewiring for small-world networks
  int clq_sz;//initial clique size for small-world networks (must be even)
  int scalefree_clq;//initial clique size for scale free networks
  //initial ring size for ring-based preferential attachement networks
  //(potentially) symmetric or not for gamma coupling
  int param2=0;

  double meandeg;
  double maxleak;

  //constout=0: constant total infectivity, variable in-infection,
  // and constant *total* out-infection;
  //constout=1: constant total infectivity, constant in-infection
  // and constant *total* out-infection;
  //constout=2: variable total infectivity, constant in-infection and
  //constant out-infection on each edge
  int constout=2;

  int i;
  FILE *fd, *fd1;
  int **sparse=NULL;
  int totparams=21;
  double stepsz;
  double meansz;
  double det_sz;//deterministic expected epidemic size

  double mean_eps, epsvar;//for variable leak simulations
  double **epsM; //The coupling matrix
  int **epsI; //The network structure
  
  int *S0, *R0;
  double *V;

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
      outfile=argv[optind];
      break;
    case 1:
      sims_per_param=atoi(argv[optind]);
      break;
    case 2:
      RepNo=atof(argv[optind]);
      break;
    case 3:
      leakmin=atof(argv[optind]);
      break;
    case 4:
      leakmax=atof(argv[optind]);
      break;
    case 5:
      Nc=atoi(argv[optind]);
      break;
    case 6:
      gridtype=atoi(argv[optind]);
      break;
    case 7:
      arg8=argv[optind];
      break;
    case 8:
      arg9=argv[optind];
      break;
    case 9:
      arg10=argv[optind];
      break;
    case 10:
      arg11=argv[optind];
      break;
    }
    mainargs++;
    optind++;
  }

  ////////////////////
  if(mainargs<8){
    fprintf(stderr, "ERROR: expecting at least the following arguments:\n\toutput file\n\t\"totruns\" (number of runs per parameter)\n\t\"RepNo\" (basic reproduction number)\n\tminimum \"leak\"\n\tmaximum \"leak\"\n\t\"Nc\" (compartment population)\n\t\"gridtype\" (see source file for available types)\n\tnumber of compartments\"Nc\", or \"Nx\" and \"Ny\" for a 2D grid\n\nEXITING\n");exit(0);
  }

  
  if(gridtype==4){//2d grid small world
    if(mainargs<10){
      fprintf(stderr, "for a 2d grid you must specify the dimensions \"Nx\" and \"Ny\", and the rewiring probability \"rewiring_prob\". EXITING\n");exit(0);
    }
    Nx=atoi(arg8);Ny=atoi(arg9);ncomp=Nx*Ny;rewiring_prob=atof(arg10);
    meandeg=4;
  }
  else{
    ncomp=atoi(arg8);//number of compartments
    if(gridtype==1)
      meandeg=(double)(ncomp-1);
    else if(gridtype==3){//ring small world
      if(mainargs<10){
	fprintf(stderr, "for a (ring) small-world network you must specify:\n\tthe initial clique size \"clq_sz\" (must be even), and\n\tthe rewiring probability \"rewiring_prob\".\n\nEXITING\n");exit(0);
      }
      clq_sz=atoi(arg9);
      rewiring_prob=atof(arg10);
      meandeg=clq_sz;
    }
    else if(gridtype==2){//gamma
      fprintf(stderr, "Leak effect routine currently not suitable for gamma-distributed coupling strengths. (Resetting the grid each time could be misleading). EXITING.\n"); exit(0);
      /* if(argc<10){ */
      /* 	fprintf(stderr, "for gamma distributed weights, you must specify the shape parameter \"shp\". EXITING\n");exit(0); */
      /* } */
      /* shp=atof(argv[9]); */
    }
    else if(gridtype==5){//scale free
      if(mainargs<9){
	fprintf(stderr, "for a scale free network, you must specify the initial clique size (a positive integer)\n\nEXITING\n");exit(0);
      }
      scalefree_clq=atoi(arg9);
      if(scalefree_clq<=0){
	fprintf(stderr, "for a scale-free network, the initial clique size must be a positive integer (currently %d)\n\nEXITING\n", scalefree_clq);exit(0);
      }
      meandeg=(double)(scalefree_clq*(2*ncomp-scalefree_clq-1))/(double)ncomp;
    }
    else if(gridtype==6){//random
      if(mainargs<9){
	fprintf(stderr, "for a random network, you must specify the mean degree.\nEXITING\n");exit(0);
      }
      param1=atof(arg9);
      meandeg=param1;
    }
    else if(gridtype==7){//ring scale-free
      if(mainargs<11){
	fprintf(stderr, "For a (ring) preferential attachment network you must specify:\n\tthe ring size \"init_ring\", \n\tthe initial clique size \"clq_sz\" (must be even), and\n\tnumber of new edges at each step.\n\nEXITING\n");exit(0);
      }
      param2=atoi(arg9);
      clq_sz=atoi(arg10);
      scalefree_clq=atoi(arg11);
      meandeg=(double)(param2*clq_sz+2*scalefree_clq*(ncomp-param2))/(double)ncomp;
    }
    else if(gridtype==8){//star
      meandeg=(double)(2*(ncomp-1))/(double)ncomp;
    }
    else{
      fprintf(stderr, "\"gridtype\" must be in the range 1-8. (Currently %d).\nEXITING.\n", gridtype);exit(0);
    }

  }


  S0=(int*) malloc(sizeof(int) * (ncomp+1));
  R0=(int*) malloc(sizeof(int) * (ncomp+1));
  V=(double*) malloc(sizeof(double) * (ncomp+1));

  epsM=dmatrix(ncomp,ncomp); //The coupling matrix
  epsI=imatrix(ncomp,ncomp); //The coupling matrix structure

  stepsz=(leakmax-leakmin)/(double)(totparams-1);
  fd=fopen(outfile, "w");
  fd1=fopen("SIRgeneral/outfiles/network_properties", "w");


  compV=(double)Nc;//set population density to 1
  ki=RepNo*kr*compV/Nc;//infection rate constant
    

  maxleak=meandeg*(double)Nc/(meandeg+1.0);
  //Initialisation
  inittoconst(V,ncomp,compV);//Volumes all equal
  inittoconst(S0,ncomp,Nc);
  inittozero(R0,ncomp);
  
  //set up grid
  setgrid(rng, gridtype, epsM, epsI, ncomp, Nx, Ny, Nc, &sparse, RepNo, leakmin, leakmin/((double)Nc), param1, rewiring_prob, param2, clq_sz, scalefree_clq, &mean_eps, &epsvar, constout, fd1);


  for(i=0;i<totparams;i++){
    leak = leakmin+i*stepsz;
    eps=leak/((double)Nc);

    if(gridtype==1){//complete
      sym_coupling(eps, epsM, ncomp);
    }
    else{//all others: reset coupling matrix, but not structure
      coupling_from_epsI(eps, epsI, epsM, ncomp, constout);
    }

    //printmat(epsM,ncomp,ncomp);fflush(stderr);
    det_sz=(1.0+boost::math::lambert_w0(-1.0*RepNo*exp(-RepNo))/RepNo);
    meansz=mean_SIR_outbreak(rng, ncomp, V, S0, R0, epsM, ki, kr, tmax, sims_per_param, 0, sparse);
    fprintf(stderr, "leak=%.4f%% size=%.4f%%\n", leak/maxleak*100.0, meansz*100/det_sz/(1.0-1.0/RepNo));
    fprintf(fd, "%.4f %.4f\n", leak/maxleak, meansz/det_sz/(1.0-1.0/RepNo));
    //fprintf(stderr, "%.2f, %.4f, %.4f\n", leak, mean_outbreak_sym(leak, RepNo, Nc, ncomp));
  }

  
  fclose(fd);fclose(fd1);
  free_dmatrix(epsM,ncomp,ncomp);
  free_imatrix(epsI,ncomp,ncomp);
  free((char*)S0);free((char*)R0);free((char*)V);

  if(sparse)
    free_imat(sparse,ncomp);

  return 0;
}
