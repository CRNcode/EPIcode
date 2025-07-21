
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

// SIR epidemics on grids: specialised for generating videos
// We can select an epidemic hitting roughly a given proportion of
// compartments through "size_select"

int main(int argc, char *argv[]){

  if(argc<8){
    fprintf(stderr, "Insufficient parameters. These should be:\n\t(1) number of simulations\n\t(2) basic reproduction number (R0)\n\t(3) \"leak\" (fraction of out-infections X compartment population)\n\t(4) \"Nc\" (compartment population)\n\t(5) \"Nx\" the x-dimension of the grid, and\n\t(6) \"Ny\", the y dimension of the grid\n\t(7) rewiring probability.\nEXITING.\n");exit(0);
  }

  int opt;
  double tmax=500.0;//cut off by time tmax if extinction isn't achieved
  int totruns=atoi(argv[1]);
  double RepNo=atof(argv[2]);//set R0, infer ki
  double leak=atof(argv[3]);//Given Nc introductions, how many are external?
  int Nc=atoi(argv[4]);//Total population in each compartment
  int Nx=atoi(argv[5]);//horizontal dimension of grid (>=3)
  int Ny=atoi(argv[6]);//vertical dimension of grid (>=3)
  double rewiring_prob=atof(argv[7]);
  double compV=(double)Nc;//set population density to 1
  int ncomp=Nx*Ny;
  double eps=leak/((double)Nc);//strength of total cross infection
  double kr=1.0;
  double ki=RepNo*kr*compV/Nc;
  int Iinit=5;//initial infections

  int i,j,jrun,iter,totepis=0;
  int printall=1;//set to 0 if we don't want to store the evolution of first run

  double t;
  double vaxfrac=0.0;

  int **I0, **S0, **R0;
  int Rtot,S0tot=0,R0tot=0;
  int Episz[totruns];
  int Episzhist[ncomp+1];
  double Rt_classical[totruns], Rt_nextgen[totruns];//reproduction numbers
  double q[ncomp];//vector of extinction probabilities
  double outbreak_P;//probability of outbreak, given introduction
  double outfrac[totruns];//fraction infected during outbreak
  
  double **VM=dmatrix(Nx,Ny); //compartment volumes
  double **epsM=dmatrix(ncomp,ncomp); //The coupling matrix
  int **epsI=imatrix(ncomp,ncomp); //The coupling matrix
  double **NextGen=dmatrix(ncomp,ncomp);//Next generation matrix
  double mean_eps, epsvar;
  
  int **I,**S,**R;
  double **mean_outbreak_sz, mean_outbreak_sztot;
  double **mean_outbreak_sz_cond, mean_outbreak_sz_condtot;

  FILE *fd1,*fd4,*fd5;
  FILE *fd=fopen("SIRgeneral/outfiles/network_properties", "w");
  FILE *fd2=fopen("SIRtorus/outfiles/SIRtorusfinal", "w");
  FILE *fd3=fopen("SIRtorus/outfiles/SIRtoruseps", "w");
 
  double meanfinalT=0.0, det_sz;
  int sampling=10;//for even time steps (number per second)
  double Icomp=(double)(Nc/5);
  double Scomp=(double)(Nc);
  int **sparse;

  //constout=0: constant total infectivity, variable in-infection,
  // and constant *total* out-infection;
  //constout=1: constant total infectivity, constant in-infection
  // and constant *total* out-infection;
  //constout=2: variable total infectivity, constant in-infection and
  //constant out-infection on each edge
  int constout=2;
  int size_select=2;//2:try to find a medium sized epi
  int epimin=0, epimax=Nx*Ny;
  int found=0; // searching for an epi of roughly a given size?

  static std::mt19937 rng(std::random_device{}()); // Standard mersenne_twister_engine seeded with rd()

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


  if(size_select==1){
    epimin=(int)(0.02*Nx*Ny);epimax=(int)(0.1*Nx*Ny);
  }
  else if(size_select==2){
    epimin=(int)(0.08*Nx*Ny);epimax=(int)(0.18*Nx*Ny);
  }
  else if(size_select==3){
    epimin=(int)(0.18*Nx*Ny);epimax=(int)(0.38*Nx*Ny);
  }
  else if(size_select==4){
    epimin=(int)(0.38*Nx*Ny);epimax=(int)(0.58*Nx*Ny);
  }

  S=imatrix(Nx, Ny);I=imatrix(Nx, Ny);R=imatrix(Nx, Ny);
  S0=imatrix(Nx, Ny);I0=imatrix(Nx, Ny);R0=imatrix(Nx, Ny);
  mean_outbreak_sz=dmatrix(Nx, Ny);mean_outbreak_sz_cond=dmatrix(Nx, Ny);

  inittoconst(VM,Nx,Ny,compV);
  inittozero(epsI, ncomp, ncomp);inittozero(epsM, ncomp, ncomp);
  setgrid(rng, 4, epsM, epsI, ncomp, Nx, Ny, Nc, &sparse, RepNo, leak, eps, 0, rewiring_prob, 0, 0, 0, &mean_eps, &epsvar, constout, fd);

  
  //initial populations
  inittozero(I0, Nx, Ny);
  I0[Nx/2][Ny/2]=Iinit;// initial infections
  
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      S0[i][j]=Nc-I0[i][j];R0[i][j]=0;
    }
  }

  inittozero(mean_outbreak_sz,Nx,Ny);inittozero(mean_outbreak_sz_cond,Nx,Ny);
  inittozero(Episz,totruns);inittozero(Episzhist,ncomp+1);

  //vaccinate some fraction of the population?
  vaccinate(rng, Nx, Ny, S0, R0, vaxfrac);
  //printmat(S0, Nx, Ny);exit(0);

  S0tot=msum(S0,Nx,Ny);R0tot=msum(R0,Nx,Ny);
  fprintf(fd3, "epsilon\ttheory\tsimulation\t(C1\tC2)\n");
  for(jrun=0;jrun<totruns;jrun++){//main loop

    fd1=fopen("SIRtorus/outfiles/SIRtorusmean", "w");
    fd4=fopen("SIRtorus/outfiles/SIRtorusout", "w");

    fprintf(stderr, "%d/%d\n",jrun+1,totruns);

    if((t=runSIRepi(rng, Nx, Ny, VM, I0, S0, R0, I, S, R, Iinit, Icomp, Scomp, epsM, ki, kr, tmax, sampling, printall, fd4, fd1, sparse))<0.0){
      fprintf(stderr,"Epidemic duration longer than tmax. EXITING.\n");exit(0);
    }

     
    //fprintf(stderr, "got here1\n");fflush(stderr);

    //Epidemic over
    Episz[jrun]=episz(R, R0, Nx, Ny, Nc, mean_outbreak_sz, mean_outbreak_sz_cond, &totepis);
    Episzhist[Episz[jrun]]++;meanfinalT+=t;
    //fprintf(stderr, "Epi size: %d\n", Episz[jrun]);

    Rtot=msum(R,Nx,Ny);
    if(size_select && Episz[jrun]>=epimin && Episz[jrun]<=epimax){
      fprintf(stderr, "Epi size: %d\n", Episz[jrun]);
      found=1;
    }
       
    //Fraction of pop. infected
    outfrac[jrun]=((double)(Rtot-R0tot))/((double)(ncomp*Nc));
    //classical and next-gen effective Rep No
    Rt_classical[jrun]=Rt_clas_nn(S,VM,epsM,ki,kr,Nx,Ny);
    NextG(NextGen, S, VM, epsM, ki, kr, Nx, Ny);
    Rt_nextgen[jrun]=SpecRad(NextGen,ncomp,1e-7);
    //fprintf(stderr, "got here2\n");fflush(stderr);
  
    //outbreak probability at end of epidemic (branching process approx.)
    //Athreya and Ney, V.3, Thm 2. "q"
    outbreak_P=outbreak_prob(ki, kr, S, VM, epsM, Nx, Ny, q, 1e-6, &iter);
	
    fprintf(fd2, "%.4f", t);
    /* for(i=0;i<Nx;i++) */
    /*   for(j=0;j<Ny;j++) */
    /*     fprintf(fd2, ", %d", R[i][j]); */
    fprintf(fd2, ", %d, %.4f%%, %.4f%%, %.4f, %.4f\n", Episz[jrun], outfrac[jrun]*100.0, outbreak_P*100.0, Rt_classical[jrun], Rt_nextgen[jrun]);
    fclose(fd1);
    fclose(fd4);

    if(found==1){
      fd5=fopen("SIRtorus/gnuargs", "w");//for plotting purposes
      fprintf(fd5, "%d %.1f %d %d %.1f %.1f %d %.2f\n", (int)(sampling*t + 100.0), RepNo, Nx, Ny, ((double)sampling), leak, Nc, rewiring_prob);
      fclose(fd5);
      break;
    }

  }
  //fprintf(stderr, "got here\n");fflush(stderr);

  mean_outbreak_sztot=msum(mean_outbreak_sz,Nx,Ny);
  mean_outbreak_sz_condtot=msum(mean_outbreak_sz_cond,Nx,Ny);

  fprintf(stderr, "\nepsilon = %.4f, \"R0\" = %.4f, mean eps = %.4f, mean path=%.4f\n", eps, RepNo, mean_eps, mean_path_length(epsI,ncomp));

  //Epidemic sizes (i.e., number of compartments hit)
  fprintf(stderr, "Episzhist: ");
  for(i=0;i<=ncomp;i++)
    fprintf(stderr, "%d ", Episzhist[i]);
  fprintf(stderr, "\n");

  //The deterministic epidemic size in the one-compartment case
  //See https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology
  det_sz=(Nc+Nc*boost::math::lambert_w0(-(double)(S0tot)/(double)(ncomp*Nc)*RepNo*exp(-RepNo*(1.0-(double)(R0tot)/(double)(ncomp*Nc))))/RepNo)/((double)Nc);
  fprintf(stderr, "Deterministic final outbreak size (%%) = %.4f\n", det_sz*100.0);

  fprintf(stderr, "Stochastic mean final outbreak size (%%) = %.4f\n", mean_outbreak_sztot/(double)jrun/((double)(ncomp*Nc))*100.0);
  for(i=0;i<Ny;i++){
    for(j=0;j<Nx;j++){
      fprintf(stderr, "%.4f ", mean_outbreak_sz[j][i]/(double)jrun/((double)Nc)*100.0);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");

  if(totepis>0){
    fprintf(stderr, "Conditioned mean final outbreak size (%%) = %.4f\n", mean_outbreak_sz_condtot/(double)totepis/((double)(ncomp*Nc))*100.0);
    for(i=0;i<Ny;i++){
      for(j=0;j<Nx;j++){
	fprintf(stderr, "%.4f ", mean_outbreak_sz_cond[j][i]/((double)totepis)/((double)Nc)*100.0);
      }
      fprintf(stderr, "\n");
    }
  }
  
  fprintf(stderr, "Stochastic mean final outbreak time = %.4f\n", meanfinalT/jrun);


  fprintf(fd3, "%.4f\t%.4f\t%.4f\t(", eps, det_sz*(2.0-pow(1.0-eps,det_sz*Nc))/(2.0)*100.0, mean_outbreak_sztot/(double)jrun/((double)2*Nc)*100.0);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      fprintf(fd3, "%.4f\t", mean_outbreak_sz[i][j]/(double)jrun/((double)Nc)*100.0);
    }
  }
  fprintf(fd3, ")\n");


  free_imatrix(I, Nx, Ny);
  free_imatrix(S, Nx, Ny);
  free_imatrix(R, Nx, Ny);
  free_imatrix(I0, Nx, Ny);
  free_imatrix(S0, Nx, Ny);
  free_imatrix(R0, Nx, Ny);
  free_dmatrix(mean_outbreak_sz, Nx, Ny);
  free_dmatrix(mean_outbreak_sz_cond, Nx, Ny);
  free_dmatrix(VM, Nx, Ny);
  free_dmatrix(epsM,ncomp,ncomp); //The coupling matrix
  free_imatrix(epsI,ncomp,ncomp);
  free_dmatrix(NextGen,ncomp,ncomp);//Next generation matrix
  free_imat(sparse,ncomp);

  fclose(fd);
  fclose(fd2);
  fclose(fd3);

  
  return 0;
}
