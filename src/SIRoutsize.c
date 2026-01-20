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

// Get the expected outbreak sizes along a trajectory from data
// Need to input the correct network type.
// Starts by projecting the state onto the disease-free states

int main(int argc, char *argv[]){
  if(argc<9){
    fprintf(stderr, "ERROR: expecting at least the followng arguments:\n\tinput file\n\toutput file\n\t\"secondruns\" (to find expected outbreak size)\n\t\"RepNo\" (basic reproduction number)\n\t\"leak\" (fraction of out-infections X compartment population)\n\t\"Nc\" (compartment population)\n\t\"gridtype\" (see source file for available types)\n\tnumber of compartments\"Nc\", or \"Nx\" and \"Ny\" for a 2D grid\n\nEXITING\n");exit(0);
  }

  int opt;
  double tmax=500.0;//cut off by time tmax if extinction isn't achieved
  char *infile=argv[1];//input file
  char *outfile=argv[2];//input file
  int secondruns=atoi(argv[3]); //simulations after each outbreak
  double RepNo=atof(argv[4]);//set R0, infer ki
  double leak=atof(argv[5]);
  int Nc=atoi(argv[6]);//Total population in each compartment
  double compV=(double)Nc;//set population density to 1
  double eps=leak/((double)Nc);//strength of total cross infection
  double kr=1.0;//recovery rate constant
  double ki=RepNo*kr*compV/Nc;//infection rate constant

  //network type and coupling
  //gridtype=1:complete; 2:gamma, 3:ring small_world, 4:2d grid small world, 5:scale-free, 6:star-shaped
  int gridtype=atoi(argv[7]);
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
  
  if(gridtype==4){//2d grid small world
    if(argc<11){
      fprintf(stderr, "for a 2d grid you must specify the dimensions \"Nx\" and \"Ny\", and the rewiring probability \"rewiring_prob\". EXITING\n");exit(0);
    }
    Nx=atoi(argv[8]);Ny=atoi(argv[9]);ncomp=Nx*Ny;rewiring_prob=atof(argv[10]);
  }
  else{
    ncomp=atoi(argv[8]);//number of compartments

    if(gridtype==2){//gamma
      if(argc<11){
	fprintf(stderr, "for gamma distributed weights, you must specify the shape parameter, and whether coupling is symmetric (0=not/1=symmetric). EXITING\n");exit(0);
      }
      param1=atof(argv[9]);
      param2=atoi(argv[10]);
    }
    else if(gridtype==3){//ring small world
      if(argc<11){
	fprintf(stderr, "for a (ring) small-world network you must specify:\n\tthe initial clique size \"clq_sz\" (must be even), and\n\tthe rewiring probability \"rewiring_prob\".\n\nEXITING\n");exit(0);
      }
      clq_sz=atoi(argv[9]);
      rewiring_prob=atof(argv[10]);
    }
    else if(gridtype==5){//scale free
      if(argc<10){
	fprintf(stderr, "for a scale free network, you must specify the initial clique size (a positive integer)\n\nEXITING\n");exit(0);
      }
      scalefree_clq=atoi(argv[9]);
      if(scalefree_clq<=0){
	fprintf(stderr, "for a scale-free network, the initial clique size must be a positive integer (currently %d)\n\nEXITING\n", scalefree_clq);exit(0);
      }
    }
    else if(gridtype==6){//random
      if(argc<10){
	fprintf(stderr, "for a random network, you must specify the mean degree.\nEXITING\n");exit(0);
      }
      param1=atof(argv[9]);
    }
    else if(gridtype==7){//ring scale-free
      if(argc<12){
	fprintf(stderr, "For a (ring) preferential attachment network you must specify:\n\tthe ring size \"init_ring\", \n\tthe initial clique size \"clq_sz\" (must be even), and\n\tnumber of new edges at each step.\n\nEXITING\n");exit(0);
      }
      param2=atoi(argv[9]);
      clq_sz=atoi(argv[10]);
      scalefree_clq=atoi(argv[11]);
    }
    else if(gridtype!=1 &&gridtype!=8){
      fprintf(stderr, "\"gridtype\" must be in the range 1-8. (Currently %d).\nEXITING.\n", gridtype);exit(0);
    }

  }

  //Set this to one to choose introductions according to the susceptible population rather than the total population in each compartment
  int choose_by_S=0;

  //constout=0: constant total infectivity, variable in-infection,
  // and constant *total* out-infection;
  //constout=1: constant total infectivity, constant in-infection
  // and constant *total* out-infection;
  //constout=2: variable total infectivity, constant in-infection and
  //constant out-infection on each edge
  int constout=2;
  

  double **epsM=dmatrix(ncomp,ncomp); //The coupling matrix
  int **epsI=imatrix(ncomp,ncomp); //The coupling matrix structure
  int **sparse=NULL;

  double mean_eps, epsvar;//for variable leak simulations

  int i;

  int S0[ncomp], R0[ncomp];
  double V[ncomp];


  FILE *fdin, *fdout;
  FILE *fd1=fopen("SIRgeneral/outfiles/vax_protect", "w");//for plotting
  FILE *fd2=fopen("SIRgeneral/outfiles/network_properties", "w");//for plotting purposes

  double sz_natural;
  unsigned long numl=0;
  int maxl=0, linelen;
  char *linenew;
  char *wd, *wdtime;

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
  
  if(secondruns<0 || secondruns>100000){
    fprintf(stderr, "secondruns must be in the range 0-100000. EXITING.\n");exit(0);
  }

  numl = numflines(infile, &maxl); // number of lines and max line length
  if(numl==0){
    fprintf(stderr, "ERROR: \"%s\" is an empty file. EXITING.\n", infile);exit(0);
  }
  linenew = (char*) malloc(sizeof(char) * (maxl));

  inittoconst(V,ncomp,compV);//Volumes all equal

  //set the coupling matrix
  setgrid(rng, gridtype, epsM, epsI, ncomp, Nx, Ny, Nc, &sparse, RepNo, leak, eps, param1, rewiring_prob, param2, clq_sz, scalefree_clq, &mean_eps, &epsvar, constout, fd2);
  fclose(fd2);

  //open th2 files
  fdin=fopen(infile, "r");
  if(!(fdout=fopen(outfile, "w"))){
    fprintf(stderr, "ERROR in resample - File: %s could not be opened...\n", outfile);exit(0);
  }

  //run simulations
  //set all infecteds to be recovered (projection onto disease-free states
  //which preserves total population in each compartment
  while((linelen = gtline(fdin, linenew, maxl)) > 0){
    wdtime=getnthwd(linenew, 1);
    for(i=0;i<ncomp;i++){
      wd=getnthwd(linenew, 3*i+7);S0[i]=atoi(wd);free(wd);
      wd=getnthwd(linenew, 3*i+8);R0[i]=atoi(wd);free(wd);
      wd=getnthwd(linenew, 3*i+9);R0[i]+=atoi(wd);free(wd);
    }
    printvec(S0,ncomp);printvec(R0,ncomp);

    sz_natural=mean_SIR_outbreak(rng, ncomp, V, S0, R0, epsM, ki, kr, tmax, secondruns, choose_by_S, sparse);
    fprintf(stderr, "%s %.4f\n", wdtime, sz_natural*100.0);
    fprintf(fdout, "%s %.4f\n", wdtime, sz_natural);
    free(wdtime);
  }
  

  fclose(fd1);
  fclose(fdin);fclose(fdout);
  free((char*)linenew);
  free_dmatrix(epsM,ncomp,ncomp);
  free_imatrix(epsI,ncomp,ncomp);
  if(sparse)
    free_imat(sparse,ncomp);

  
  return 0;
}
