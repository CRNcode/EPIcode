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

// Variable graph topology, but identical compartments
// Ii+Sj->Ii+Ij (i,j=1,...,n) (infection, but no migration)
// Ii-->Ri (death or recovery with immunity)

// To compute theoretical and simulated expected outbreak sizes
// as a function of prior random immunity: illustrates that
// even when compartmentalisation does not reduce expected
// outbreak size, some random immunity can have a greater protective
// effect in a compartmental model than a homogeneous model. 

int main(int argc, char *argv[]){

  int opt, mainargs=0;
  char *arg8, *arg9, *arg10, *arg11;
  
  double tmax=500.0;//cut off by time tmax if extinction isn't achieved
  
  char *fout=argv[1];
  int numvax=atoi(argv[2]);//number of experiments
  int runs_per_vax=atoi(argv[3]); //simulations after each outbreak
  double RepNo=atof(argv[4]);//set R0, infer ki
  double leak=atof(argv[5]);
  int Nc=atoi(argv[6]);//Total population in each compartment
  int gridtype=atoi(argv[7]);
  
  double compV;
  double eps;
  double kr=1.0;//recovery rate constant
  double ki;//infection rate constant

  //network type and coupling
  //gridtype=1:complete; 2:gamma, 3:ring small_world, 4:2d grid small world, 5:scale-free, 6:star-shaped, 7: ring scale-free

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

  //constout=0: constant total infectivity, variable in-infection,
  // and constant *total* out-infection;
  //constout=1: constant total infectivity, constant in-infection
  // and constant *total* out-infection;
  //constout=2: variable total infectivity, constant in-infection and
  //constant out-infection on each edge
  int constout=2;

  double **epsM; //The coupling matrix
  int **epsI; //The coupling matrix structure
  double **NextGen; //The next generation matrix
  int **sparse=NULL;

  double mean_eps, epsvar;//for variable leak simulations

  int i,j;

  int *S1, *R1;
  double *V;
  
  FILE *fd1;//for plotting
  FILE *fd2=fopen("SIRgeneral/outfiles/network_properties", "w");//for plotting purposes

  double sz_vax=0, sz_vax_max;
  double S0frac, R0frac, det_sz_frac, det_sz_frac_cond, det_sz_max;
  double scalefac;
  int theory_only=0;//set to 1 to do no simulations

  static std::mt19937 rng(std::random_device{}()); 
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
      fout=argv[optind];
      break;
    case 1:
      numvax=atoi(argv[optind]);
      break;
    case 2:
      runs_per_vax=atoi(argv[optind]);
      break;
    case 3:
      RepNo=atof(argv[optind]);
      break;
    case 4:
      leak=atof(argv[optind]);
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


  if(mainargs<8){
    fprintf(stderr, "ERROR: expecting at least the followng arguments:\n\toutput file name\n\t\"numvax\" (number of population vaccinations)\n\t\"runs_per_vax\" (runs per vaccination state)\n\t\"RepNo\" (basic reproduction number)\n\t\"leak\" (fraction of out-infections X compartment population)\n\t\"Nc\" (compartment population)\n\t\"gridtype\" (see source file for available types)\n\tnumber of compartments\"Nc\", or \"Nx\" and \"Ny\" for a 2D grid\n\nEXITING\n");exit(0);
  }
    
  if(gridtype==4){//2d grid small world
    if(mainargs<10){
      fprintf(stderr, "for a 2d grid you must specify the dimensions \"Nx\" and \"Ny\", and the rewiring probability \"rewiring_prob\". EXITING\n");exit(0);
    }
    Nx=atoi(arg8);Ny=atoi(arg9);ncomp=Nx*Ny;rewiring_prob=atof(arg10);
  }
  else{
    ncomp=atoi(arg8);//number of compartments
    if(gridtype==2){//gamma
      if(mainargs<10){
	fprintf(stderr, "for gamma distributed weights, you must specify the shape parameter, and whether coupling is symmetric (0=not/1=symmetric). EXITING\n");exit(0);
      }
      param1=atof(arg9);
      param2=atoi(arg10);
    }
    if(gridtype==3){//ring small world
      if(mainargs<10){
	fprintf(stderr, "for a (ring) small-world network you must specify:\n\tthe initial clique size \"clq_sz\" (must be even), and\n\tthe rewiring probability \"rewiring_prob\".\n\nEXITING\n");exit(0);
      }
      clq_sz=atoi(arg9);
      rewiring_prob=atof(arg10);
    }
    else if(gridtype==5){//scale free
      if(mainargs<9){
	fprintf(stderr, "for a scale free network, you must specify the initial clique size (a positive integer)\n\nEXITING\n");exit(0);
      }
      scalefree_clq=atoi(arg9);
      if(scalefree_clq<=0){
	fprintf(stderr, "for a scale-free network, the initial clique size must be a positive integer (currently %d)\n\nEXITING\n", scalefree_clq);exit(0);
      }
    }
    else if(gridtype==6){//random
      if(mainargs<9){
	fprintf(stderr, "for a random network, you must specify the mean degree.\nEXITING\n");exit(0);
      }
      param1=atof(arg9);
    }
    else if(gridtype==7){//ring scale-free
      if(mainargs<11){
	fprintf(stderr, "For a (ring) preferential attachment network you must specify:\n\tthe ring size \"init_ring\", \n\tthe initial clique size \"clq_sz\" (must be even), and\n\tnumber of new edges at each step.\n\nEXITING\n");exit(0);
      }
      param2=atoi(arg9);
      clq_sz=atoi(arg10);
      scalefree_clq=atoi(arg11);
    }
    else if(gridtype!=1){
      fprintf(stderr, "\"gridtype\" must be 1 or 3--7. (Currently %d).\nEXITING.\n", gridtype);exit(0);
    }
  }


  
  compV=(double)Nc;//set population density to 1
  eps=leak/((double)Nc);
  ki=RepNo*kr*compV/(double)Nc;
  epsM=dmatrix(ncomp,ncomp);
  epsI=imatrix(ncomp,ncomp);
  NextGen=dmatrix(ncomp,ncomp);
  S1=(int*) malloc(sizeof(int) * (ncomp+1));
  R1=(int*) malloc(sizeof(int) * (ncomp+1));
  V=(double*) malloc(sizeof(double) * (ncomp+1));
  fd1=fopen(fout, "w");//for plotting


  if(numvax<0 || numvax>1000){
    fprintf(stderr, "numvax must be in the range 0-1000. EXITING.\n");exit(0);
  }
  if(runs_per_vax<0 || runs_per_vax>100000){
    fprintf(stderr, "runs_per_vax must be in the range 0-100000. EXITING.\n");exit(0);
  }


  inittoconst(V,ncomp,compV);//Volumes all equal

  //set the coupling matrix
  setgrid(rng, gridtype, epsM, epsI, ncomp, Nx, Ny, Nc, &sparse, RepNo, leak, eps, param1, rewiring_prob, param2, clq_sz, scalefree_clq, &mean_eps, &epsvar, constout, fd2);
  fclose(fd2);


  //for(i=100;i>=0;i--){
  for(i=100;i>=0;i-=5){
    S0frac=(double)i/double(100);R0frac=1.0-S0frac;
    scalefac=S0frac>0?max((1.0-1.0/(S0frac*RepNo)),0):0;
    det_sz_frac_cond=(S0frac+boost::math::lambert_w0(-S0frac*RepNo*exp(-RepNo*S0frac))/RepNo);
    det_sz_frac=scalefac*det_sz_frac_cond;//homogeneous case expected size
    //fprintf(stderr, "%.4f\n", mean_outbreak_sym(leak, RepNo, Nc, ncomp, 1.0-S0frac));
    if(!theory_only){//run simulations
      sz_vax=0;
      for(j=0;j<numvax;j++){
	inittoconst(S1,ncomp,Nc);inittozero(R1,ncomp);
	vaccinate(rng, ncomp, S1, R1, (double)(100-i)/(double)100);
	sz_vax+=mean_SIR_outbreak(rng, ncomp, V, S1, R1, epsM, ki, kr, tmax, runs_per_vax, choose_by_S, sparse);
      }
      sz_vax/=(double)numvax;
      if(i==100){
	sz_vax_max=sz_vax;det_sz_max=det_sz_frac;
      }
      if(gridtype==1 && ncomp<=12)//complete, symmetric: also print the expected sizes
	fprintf(fd1, "%.2f %.4f %.4f %.4f %.4f %.4f\n", R0frac, det_sz_frac, scalefac*mean_outbreak_sym(leak, RepNo, Nc, ncomp, 1.0-S0frac)*det_sz_frac_cond, det_sz_max,sz_vax,sz_vax_max);
      else
	fprintf(fd1, "%.2f %.4f %.4f %.4f %.4f\n", R0frac, det_sz_frac,det_sz_max,sz_vax,sz_vax_max);
      fprintf(stderr, "%.2f %.4f\n", R0frac, sz_vax);
    }
    else
      fprintf(fd1, "%.2f %.4f\n", R0frac, scalefac*mean_outbreak_sym(leak, RepNo, Nc, ncomp, 1.0-S0frac)*det_sz_frac_cond);

    

  }
  
  fclose(fd1);
  free_dmatrix(NextGen,ncomp,ncomp);
  free_dmatrix(epsM,ncomp,ncomp);
  free_imatrix(epsI,ncomp,ncomp);
  free((char*)S1);free((char*)R1);
  if(sparse)
    free_imat(sparse,ncomp);

  
  return 0;
}
