/* Copyright (C) 2025, Murad Banaji
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
  char *arg6, *arg7, *arg8, *arg9;
  
  double tmax=500.0;//cut off by time tmax if extinction isn't achieved
  
  int totruns;//number of experiments
  double RepNo;//set R0, infer ki
  double leak;
  int Nc;//Total population in each compartment
  int gridtype;//1:complete; 2:gamma, 3:ring small_world, 4:2d grid small world, 5:scale-free, 6:random, 7:ring based scale free, 8: star-shaped

  double compV;
  double eps;//strength of total cross infection
  double kr=1.0;//recovery rate constant
  double ki;//infection rate constant
  int Iinit=1;//initial infections
  int compinit;
  double r1;

  int Nx,Ny;
  int ncomp;//number of compartments
  //param1: shape parameter for gamma distributed coupling strengths;
  //mean degree for random networks
  double param1;
  double rewiring_prob;//in [0,1]: probability of rewiring for small-world networks
  int clq_sz;//initial clique size for small-world networks (must be even)
  int scalefree_clq;//initial clique size for scale free networks
  //initial ring size for ring-based preferential attachement networks
  //symmetric or not for gamma coupling
  int param2=0;

  //constout=0: constant total infectivity, variable in-infection,
  // and constant *total* out-infection;
  //constout=1: constant total infectivity, constant in-infection
  // and constant *total* out-infection;
  //constout=2: variable total infectivity, constant in-infection and
  //constant out-infection on each edge
  int constout=2;

  double **epsM;
  int **epsI;
  double **NextGen; //The next generation matrix
  double mean_eps, epsvar;//for variable leak simulations

  int totpop;
  
  int printfull=1;//full output: long file; safest set to 0
  int printmean=0;//do we want the mean output?

  int i,i1,j,jrun,iter,totepis;

  int timeint;
  time_t timepoint;
  double t=0.0;

  int *S0, *I0, *R0, *S, *I, *R;
  double *V;

  int S0tot=0,R0tot=0, Stot,Rtot;
  int *Episz;//number of comps. infected at the end of Epi.
  int *Episzhist;//histogram of outbreak sizes
  double *outfrac;//fraction infected during outbreak
  double *Rt_classical, *Rt_nextgen;//reproduction numbers
  double *q;//vector of extinction probabilities
  double outbreak_P;//probability of outbreak, given introduction
 
  FILE *fd=fopen("SIRgeneral/outfiles/Epincompout", "w");
  FILE *fd1=fopen("SIRgeneral/outfiles/Epincompmean", "w");
  FILE *fd2=fopen("SIRgeneral/outfiles/Epincompfinal", "w");
  FILE *fd3=fopen("SIRgeneral/outfiles/Epincompeps", "w");
  FILE *fd4=fopen("SIRgeneral/outfiles/Epincompduration.csv", "w");
  FILE *fd6=fopen("SIRgeneral/outfiles/vax_protect", "w");//for plotting purposes
  FILE *fd7=fopen("SIRgeneral/outfiles/network_properties", "w");


  //mean of the simulations
  int meantot=1000;//total time points to output
  double **meanS, **meanI, **meanR;
  double dt=(double)meantot/tmax;
  
  //outbreak size: raw, and conditioned on at least one compartment being hit.
  double *mean_outbreak_sz, mean_outbreak_sztot;
  double *mean_outbreak_sz_cond, mean_outbreak_sz_condtot;
  double meanfinalT=0.0;
  double det_sz;//deterministic expected epidemic size

  double **alphap;
  int **sparse=NULL;
  int brief=1;//avoid calculating Rt_nextgen and outbreak_P
  int Epihistnum=20;
  int EpiHist[Epihistnum+1];
  int Epihistmax=0;

  static std::mt19937 rng(std::random_device{}()); // Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> runif(0.0, 1.0);



  //reset constout, and brief
  while ((opt = getopt(argc, argv, "c:b:")) != -1) {
    switch (opt) {
    case 'c': //reset constout
      constout = atoi(optarg);
      break;
    case 'b': //reset constout
      brief = atoi(optarg);
      break;
    default: /* '?' */
      fprintf(stderr, "Option not recognised");
      exit(EXIT_FAILURE);
    }
  }

  
  while(optind < argc){
    switch (mainargs) {
    case 0:
      totruns=atoi(argv[optind]);
      break;
    case 1:
      RepNo=atof(argv[optind]);
      break;
    case 2:
      leak=atof(argv[optind]);
      break;
    case 3:
      Nc=atoi(argv[optind]);
      break;
    case 4:
      gridtype=atoi(argv[optind]);
      break;
    case 5:
      arg6=argv[optind];
      break;
    case 6:
      arg7=argv[optind];
      break;
    case 7:
      arg8=argv[optind];
      break;
    case 8:
      arg9=argv[optind];
      break;
    }
    mainargs++;
    optind++;
  }



  if(mainargs<6){
    fprintf(stderr, "ERROR: expecting at least the following arguments:\n\t\"totruns\" (number of base-runs)\n\t\"RepNo\" (basic reproduction number)\n\t\"leak\" (fraction of out-infections X compartment population)\n\t\"Nc\" (compartment population)\n\t\"gridtype\" (see source file for available types)\n\tnumber of compartments\"Nc\", or \"Nx\" and \"Ny\" for a 2D grid\n\nEXITING\n");exit(0);
  }
  if(gridtype==4){//2d grid small world
    if(mainargs<8){
      fprintf(stderr, "for a 2d grid you must specify the dimensions \"Nx\" and \"Ny\", and the rewiring probability \"rewiring_prob\". EXITING\n");exit(0);
    }
    Nx=atoi(arg6);Ny=atoi(arg7);ncomp=Nx*Ny;rewiring_prob=atof(arg8);
  }
  else{
    ncomp=atoi(arg6);//number of compartments
    if(gridtype==2){//gamma
      if(mainargs<8){
	fprintf(stderr, "for gamma distributed weights, you must specify the shape parameter, and whether coupling is symmetric (0=not/1=symmetric). EXITING\n");exit(0);
      }
      param1=atof(arg7);param2=atoi(arg8);
    }
    else if(gridtype==3){//ring small world
      if(mainargs<8){
	fprintf(stderr, "for a (ring) small-world network you must specify:\n\tthe initial clique size \"clq_sz\" (must be even), and\n\tthe rewiring probability \"rewiring_prob\".\n\nEXITING\n");exit(0);
      }
      clq_sz=atoi(arg7);rewiring_prob=atof(arg8);

    }
    else if(gridtype==5){//scale free
      if(mainargs<7){
	fprintf(stderr, "for a scale free network, you must specify the initial clique size (a positive integer)\n\nEXITING\n");exit(0);
      }
      scalefree_clq=atoi(arg7);
    }
    else if(gridtype==6){//random
      if(mainargs<7){
	fprintf(stderr, "for a random network, you must specify the mean degree.\nEXITING\n");exit(0);
      }
      param1=atof(arg7);
    }
    else if(gridtype==7){//ring scale-free
      if(mainargs<9){
	fprintf(stderr, "For a (ring) preferential attachment network you must specify:\n\tthe ring size \"init_ring\", \n\tthe initial clique size \"clq_sz\" (must be even), and\n\tnumber of new edges at each step.\n\nEXITING\n");exit(0);
      }
      param2=atoi(arg7);clq_sz=atoi(arg8);scalefree_clq=atoi(arg9);      
    }
    else if(gridtype!=1 &&gridtype!=8){
      fprintf(stderr, "\"gridtype\" must be in the range 1-8. (Currently %d).\nEXITING.\n", gridtype);exit(0);
    }

  }


  compV=(double)Nc;//set population density to 1
  eps=leak/((double)Nc);
  ki=RepNo*kr*compV/Nc;
  epsM=dmatrix(ncomp,ncomp);epsI=imatrix(ncomp,ncomp);
  NextGen=dmatrix(ncomp,ncomp);
  totpop=ncomp*Nc;
  S0=(int*) malloc(sizeof(int) * (ncomp+1));
  I0=(int*) malloc(sizeof(int) * (ncomp+1));
  R0=(int*) malloc(sizeof(int) * (ncomp+1));
  S=(int*) malloc(sizeof(int) * (ncomp+1));
  I=(int*) malloc(sizeof(int) * (ncomp+1));
  R=(int*) malloc(sizeof(int) * (ncomp+1));
  V=(double*) malloc(sizeof(double) * (ncomp+1));
  meanS=dmatrix(ncomp,meantot);
  meanI=dmatrix(ncomp,meantot);
  meanR=dmatrix(ncomp,meantot);
  Episz=(int*) malloc(sizeof(int) * (totruns+1));
  Episzhist=(int*) malloc(sizeof(int) * (ncomp+2));
  Rt_classical=(double*) malloc(sizeof(double) * (totruns+1));
  Rt_nextgen=(double*) malloc(sizeof(double) * (totruns+1));
  q=(double*) malloc(sizeof(double) * (ncomp+1));
  outfrac=(double*) malloc(sizeof(double) * (totruns+1));
  mean_outbreak_sz=(double*) malloc(sizeof(double) * (ncomp+1));
  mean_outbreak_sz_cond=(double*) malloc(sizeof(double) * (ncomp+1));
  alphap=alphapoly(eps, RepNo, Nc, ncomp, 1);
  
  inittozero(EpiHist, Epihistnum+1);

  //random seeding
  timeint = time(&timepoint); /* set and convert time to an integer */
  srand(timeint);

  fprintf(fd3, "epsilon\ttheory\tsimulation\t(C1\tC2)\n");
  
  //Initialisation
  inittoconst(V,ncomp,compV);//Volumes all equal

  //set the coupling matrix
  setgrid(rng, gridtype, epsM, epsI, ncomp, Nx, Ny, Nc, &sparse, RepNo, leak, eps, param1, rewiring_prob, param2, clq_sz, scalefree_clq, &mean_eps, &epsvar, constout, fd7);
  fclose(fd7);

  if(Iinit>Nc){
    fprintf(stderr, "ERROR in SIRgeneral: the initial infections \"Iinit\" (%d) are all put in one compartment, and so must be fewer than or equal to the compartment population \"Nc\" (%d). EXITING.\n", Iinit, Nc);exit(0);
  }
  //exit(0);
  
  //printmat(epsM,ncomp,ncomp);
  /* exit(0); */

  //
  // Start of simulation
  //
  for(i1=0;i1<1;i1++){//varying epsilon?
    inittozero(mean_outbreak_sz,ncomp);inittozero(mean_outbreak_sz_cond,ncomp);
    inittozero(Episz,totruns);
    inittozero(Episzhist,ncomp+1);
    meanfinalT=0.0;
    totepis=0;
    
    //eps=i1*0.0005;
    //eps=0.003;
    for(j=0;j<meantot;j++){
      for(i=0;i<ncomp;i++){
	meanS[i][j]=0.0;meanI[i][j]=0.0;meanR[i][j]=0.0;
      }
    }

    R0tot=0;S0tot=totpop-Iinit;

    //for each run at current epsilon
    for(jrun=0;jrun<totruns;jrun++){
      fprintf(stderr, "%d/%d\n", jrun+1, totruns);

      //"Iinit" infections in a randomly chosen compartment
      //randomness is important in the case of non-regular networks
      do{r1 = runif(rng);}while(r1==1.0);
      compinit=(int)(r1*ncomp);
      for(i=0;i<ncomp;i++){
	if(i==compinit){I0[i]=Iinit;}else{I0[i]=0;}
	S0[i]=Nc-I0[i];R0[i]=0;
      }

      //Run the epidemic simulation
      if((t=runSIRepi(rng, ncomp, V, I0, S0, R0, I, S, R, meanI, meanS, meanR, epsM, ki, kr, tmax, printfull, printmean, meantot, dt, fd, sparse))<0.0){
	fprintf(stderr,"Epidemic duration longer than tmax. EXITING.\n");exit(0);
      }
      //Epidemic over
      //get details of the epidemic: number of compartments hit (heuristic)
      Episz[jrun]=episz(R, R0, ncomp, Nc, mean_outbreak_sz, mean_outbreak_sz_cond, &totepis);
      Episzhist[Episz[jrun]]++;meanfinalT+=t;
      //fprintf(stderr, "Epi size: %d\n", Episz[jrun]);

      Stot=vsum(S,ncomp);Rtot=vsum(R,ncomp);
      //Fraction of pop. infected
      outfrac[jrun]=((double)(Rtot-R0tot))/((double)(totpop));
      EpiHist[(int)(Epihistnum*outfrac[jrun])]++;
      Epihistmax=max(Epihistmax,EpiHist[(int)(Epihistnum*outfrac[jrun])]);
	  
      //classical and next-gen effective Rep No
      Rt_classical[jrun]=Rt_clas(Stot, ki, kr, totpop);
      if(!brief){
	NextG(NextGen, S, V, epsM, ki, kr, ncomp);
	Rt_nextgen[jrun]=SpecRad(NextGen,ncomp,1e-7);	  	  
	//printmat(NextGen,ncomp,ncomp);
	//outbreak probability at end of epidemic (branching process approx.)
	//Athreya and Ney, V.3, Thm 2. "q"
	outbreak_P=outbreak_prob(ki, kr, S, V, epsM, ncomp, q, 1e-6, &iter);
      }


      if(outfrac[jrun]>0.05){//5% or more infected
	fprintf(fd2, "%.4f", t);
	/* for(i=0;i<ncomp;i++) */
	/*   fprintf(fd2, ", %d", R[i]); */
	if(!brief)
	  fprintf(fd2, ", %d, %.4f%%, %.4f%%, %.4f, %.4f\n", Episz[jrun], outfrac[jrun]*100.0, outbreak_P*100.0, Rt_classical[jrun], Rt_nextgen[jrun]);
	else
	  fprintf(fd2, ", %d, %.4f%%, %.4f\n", Episz[jrun], outfrac[jrun]*100.0, Rt_classical[jrun]);
      }

      if(Episz[jrun]==ncomp)
	fprintf(fd4, "%d, %.4f\n", Episz[jrun], t);

      fprintf(fd,"\n\n");
    }//End of loop to run epis

    //print out the mean evolution
    if(printmean)
      printSIRmean(fd1, ncomp, meantot, totruns, dt, S0, I0, R0, meanS, meanI, meanR);
 
    mean_outbreak_sztot=vsum(mean_outbreak_sz,ncomp);
    mean_outbreak_sz_condtot=vsum(mean_outbreak_sz_cond,ncomp);

    //printmat(epsM,ncomp,ncomp);
    
    fprintf(stderr, "\nepsilon = %.4f, \"R0\" = %.4f, mean eps = %.4f, mean path=%.4f\n", eps, RepNo, mean_eps, mean_path_length(epsI,ncomp));

    //The deterministic epidemic size in the one-compartment case
    //See https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology
    det_sz=(1.0+boost::math::lambert_w0(-(double)S0tot/(double)totpop*RepNo*exp(-RepNo*(1.0-(double)R0tot/(double)totpop)))/RepNo);
    fprintf(stderr, "Deterministic final outbreak size (%%) = %.4f\n", det_sz*100.0);
    fprintf(stderr, "Stochastic mean final outbreak size (%%) = %.4f\n", mean_outbreak_sztot/(double)totruns/(double)totpop*100.0);

    if((gridtype==3 || gridtype==5 || gridtype==6) && (constout==0 || constout==2)){//variable out infection
      for(i=0;i<ncomp;i++){
	fprintf(stderr, "%d: %.4f, %.4f (%d), %.4f\n", i, epsM[i][i], vsum(epsM[i],ncomp)-epsM[i][i], sparse[i][0]-1, mean_outbreak_sz[i]/(double)totruns/((double)Nc)*100.0);
      }
    }
    else{
      for(i=0;i<ncomp;i++){
	fprintf(stderr, "%.4f ", mean_outbreak_sz[i]/(double)totruns/((double)Nc)*100.0);
      }
    }

    fprintf(stderr, "\n");

    if(totepis>0){
      fprintf(stderr, "Conditioned mean final outbreak size (%%) = %.4f ( ", mean_outbreak_sz_condtot/(double)totepis/((double)(totpop))*100.0);
      for(i=0;i<ncomp;i++){
	fprintf(stderr, "%.4f ", mean_outbreak_sz_cond[i]/((double)totepis)/((double)Nc)*100.0);
      }
      fprintf(stderr, ")\n");
    }
    
    fprintf(stderr, "Stochastic mean final outbreak time = %.4f\n", meanfinalT/totruns);
    fprintf(stderr, "Episzhist: "); printvec(Episzhist,ncomp+1);

    fprintf(fd3, "%.4f\t%.4f\t%.4f\t(", eps, det_sz*(2.0-pow(1.0-eps,det_sz*Nc))/(2.0)*100.0, mean_outbreak_sztot/(double)totruns/((double)2*Nc)*100.0);
    for(i=0;i<ncomp;i++){
      fprintf(fd3, "%.4f\t", mean_outbreak_sz[i]/(double)totruns/((double)Nc)*100.0);
    }
    fprintf(fd3, ")\n");

    if(gridtype==1){//theoretical results in the symmetric, complete, case
      fprintf(stderr, "Episz_expect: ");
      for(i=0;i<=ncomp;i++){
	fprintf(stderr, "%d ", (int)(alphap[ncomp][i]*totepis));
      }
      fprintf(stderr, "\n");
    }

    //visualise outbreak sizes
    printstarhist(EpiHist, Epihistnum, Epihistmax, 100);

    
  }


  fclose(fd);
  fclose(fd1);
  fclose(fd2);
  fclose(fd3);
  fclose(fd4);
  fclose(fd6);
  free_dmatrix(NextGen,ncomp,ncomp);
  free_dmatrix(epsM,ncomp,ncomp);
  free_imatrix(epsI,ncomp,ncomp);
  free_dmatrix(meanS,ncomp,meantot);
  free_dmatrix(meanI,ncomp,meantot);
  free_dmatrix(meanR,ncomp,meantot);
  free((char*)S0);free((char*)I0);free((char*)R0);free((char*)V);
  free((char*)S);free((char*)I);free((char*)R);
  free((char*)Episz);free((char*)Episzhist);
  free((char*)Rt_classical);free((char*)Rt_nextgen);free((char*)q);
  free((char*)outfrac);free((char*)mean_outbreak_sz);
  free((char*)mean_outbreak_sz_cond);
  
  free_dmatrix(alphap,ncomp+1,ncomp+1);
  if(sparse)
    free_imat(sparse,ncomp);

  
  return 0;
}
