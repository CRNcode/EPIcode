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

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

//overloading - zero offset
double **dmatrix(long nrh, long nch){
  double **m=dmatrix(0, nrh-1, 0, nch-1);
  return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

//overloading - zero offset
void free_dmatrix(double **m, long nrh, long nch){
  free_dmatrix(m,0,nrh-1,0,nch-1);
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;

  /* allocate pointers to rows */
  m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

//overloading - zero offset
int **imatrix(long nrh, long nch){
  int **m=imatrix(0, nrh-1, 0, nch-1);
  return m;
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

//overloading - zero offset
void free_imatrix(int **m, long nrh, long nch){
  free_imatrix(m,0,nrh-1,0,nch-1);
}

//simple version
void free_imat(int **A, int numA){
  int i;
  for(i=0;i<numA;i++)
    free ((char *)(A[i]));
  if(numA && A)
    free((char *) A);
}

//output to file (double matrix)
void printmat(FILE *fd, double **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      fprintf(fd, "%.4f ", imat[i][j]);
    fprintf(fd, "\n");
  }
  fprintf(fd, "\n");
  return;
}

//output to stderr/cerr (double matrix)
void printmat(double **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      fprintf(stderr, "%.4f ", imat[i][j]);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

//output to stderr/cerr (integer matrix)
void printmat(int **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      fprintf(stderr, "%d ", imat[i][j]);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

void printsparse(int **sparse, int n){
  int i;
  for(i=0;i<n;i++){
    fprintf(stderr, "%d: ", i);
    printvec(sparse[i],sparse[i][0]+1);
  }
  return;
}

//overloading: output to file
void printsparse(FILE *fd, int **sparse, int n){
  int i;
  for(i=0;i<n;i++){
    fprintf(fd, "%d: ", i);
    fprintvec(fd, sparse[i],sparse[i][0]+1);
  }
  return;
}

void printvec(int *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    fprintf(stderr, "%d ", ivec[i]);
  fprintf(stderr, "\n");
  return;
}

void printvec(double *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    fprintf(stderr, "%.5f ", ivec[i]);
  fprintf(stderr, "\n");
  return;
}

void fprintvec(FILE *fd, int *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    fprintf(fd, "%d ", ivec[i]);
  fprintf(fd, "\n");
  return;
}


double **cpmat(double **mat, int n, int m){
  int i,j;
  double **tmp=dmatrix(0, n-1, 0, m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){tmp[i][j]=mat[i][j];}
  }
  return tmp;
}

//overloading the nonallocating version
void cpmat(double **mat, double **out, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){out[i][j]=mat[i][j];}
  }
  return;
}

int **cpmat(int **mat, int n, int m){
  int i,j;
  int **tmp=imatrix(0, n-1, 0, m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){tmp[i][j]=mat[i][j];}
  }
  return tmp;
}

//overloading the nonallocating version
void cpmat(int **mat, int **out, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){out[i][j]=mat[i][j];}
  }
  return;
}

//copy vecin to vecout
void cpvec(double *vecin, double *vecout, int n){
  int i;
  for(i=0;i<n;i++)
    vecout[i]=vecin[i];
  return;
}

//copy vecin to vecout
void cpvec(int *vecin, int *vecout, int n){
  int i;
  for(i=0;i<n;i++)
    vecout[i]=vecin[i];
  return;
}


//row total
double rowsum(double **mat, int n, int m, int rowk){
  int i;
  double tot=0.0;
  if(rowk>=n){
    fprintf(stderr, "ERROR in rowsum: j out of range. Exiting.\n");
    exit(0);
  }
  for(i=0;i<m;i++)
    tot+=mat[rowk][i];
  return tot;
}

//row total of vector
double vsum(double *v, int n){
  int i;
  double tot=0.0;
  for(i=0;i<n;i++)
    tot+=v[i];
  return tot;
}

//row total of vector
int vsum(int *v, int n){
  int i;
  int tot=0;
  for(i=0;i<n;i++)
    tot+=v[i];
  return tot;
}

//sum total of matrix entries
int msum(int **v, int n, int m){
  int i,j;
  int tot=0.0;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
    tot+=v[i][j];
  return tot;
}

double msum(double **v, int n, int m){
  int i,j;
  double tot=0.0;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
    tot+=v[i][j];
  return tot;
}

//row maximum: for nonnegative vectors
double vmax(double v[], int n){
  int i;
  double tot=v[0];
  for(i=1;i<n;i++)
    tot=max(tot,v[i]);
  return tot;
}

//row minimum: assume entries no larger than 1e6
double vmin(double v[], int n){
  int i;
  double tot=v[0];
  for(i=1;i<n;i++)
    tot=min(tot,v[i]);
  return tot;
}

//row minimum: assume entries no larger than 1e6
int vmin(int v[], int n){
  int i;
  int tot=v[0];
  for(i=1;i<n;i++)
    tot=min(tot,v[i]);
  return tot;
}

//position, not value of minimum entry
int vminpos(int v[], int n){
  int i,pos=0;
  int tot=v[0];
  for(i=1;i<n;i++){
    if(v[i]<tot){
      tot=v[i];
      pos=i;
    }
  }
  return pos;
}

//position of first nonzero entry (int vector)
int firstnonz(int *vec, int n){
  int i;
  for(i=0;i<n;i++){
    if(vec[i])
      return i;
  }
  return -1;
}


//overloading: also require to be in list Q
int vminpos(int v[], int Q[], int n){
  int i;
  int tot;
  int pos=firstnonz(Q,n);
  if(pos<0){fprintf(stderr, "ERROR in vminpos: empty list. EXITING.\n");exit(0);}
  tot=v[pos];
  for(i=pos+1;i<n;i++){
    if(v[i]<tot && Q[i]){tot=v[i];pos=i;}
  }
  return pos;
}

void inittozero(int *vec, int len){
  int k;for(k=0;k<len;k++){vec[k]=0;}
} 

void inittozero(double *vec, int len){
  int k;for(k=0;k<len;k++){vec[k]=0.0;}
}

void inittozero(int **mat, int n, int m){
  int j,k;
  for(j=0;j<n;j++){
    for(k=0;k<m;k++){
      mat[j][k]=0;
    }
  }
}

void inittozero(double **mat, int n, int m){
  int j,k;
  for(j=0;j<n;j++){
    for(k=0;k<m;k++){
      mat[j][k]=0.0;
    }
  }
}

void inittoconst(int *vec, int n, int c){
  int j;
  for(j=0;j<n;j++)
    vec[j]=c;
  return;
}

void inittoconst(double *vec, int n, double c){
  int j;
  for(j=0;j<n;j++)
    vec[j]=c;
  return;
}

void inittoconst(int **mat, int n, int m, int c){
  int j,k;
  for(j=0;j<n;j++){
    for(k=0;k<m;k++){
      mat[j][k]=c;
    }
  }
}

void inittoconst(double **mat, int n, int m, double c){
  int j,k;
  for(j=0;j<n;j++){
    for(k=0;k<m;k++){
      mat[j][k]=c;
    }
  }
}

//count the nonzero entries in an integer vector
int nonzentries(int *ivec, int n){
  int i,tot=0;
  for(i=0;i<n;i++){
    if(ivec[i])
      tot++;
  }
  return tot;
}

//overloading: count the nonzero entries in the ith column of a matrix
int nonzentries(int **imat, int i, int n, int m){
  int j,tot=0;
  if(i>=m){
    fprintf(stderr, "ERROR in \"nonzentries\": column index too large. EXITING.\n");exit(0);
  }
  for(j=0;j<n;j++){
    if(imat[j][i])
      tot++;
  }
  return tot;
}


//reindexing routine for 2d arrays (lexicographic):
//no error checking
int d2tod1(int x, int y, int Nx, int Ny){
  return y*Nx+x;
}

void d1tod2(int t, int Nx, int Ny, int *x, int *y){
  (*x)= t%Nx;
  (*y)=(t-(*x))/Nx;
  return;
}

//nearest neighbours in a 2D array
int neighbour_r(int i, int j, int Nx, int Ny){
  return d2tod1((i+1)%Nx, j, Nx, Ny);
}
//1d version
int neighbour_r(int t, int Nx, int Ny){
  int i, j;
  d1tod2(t, Nx, Ny, &i, &j);
  return neighbour_r(i, j, Nx, Ny);
}
int neighbour_l(int i, int j, int Nx, int Ny){
  return d2tod1((Nx+i-1)%Nx, j, Nx, Ny);
}
//1d version
int neighbour_l(int t, int Nx, int Ny){
  int i, j;
  d1tod2(t, Nx, Ny, &i, &j);
  return neighbour_l(i, j, Nx, Ny);
}

int neighbour_u(int i, int j, int Nx, int Ny){
  return d2tod1(i, (j+1)%Ny, Nx, Ny);
}
//1d version
int neighbour_u(int t, int Nx, int Ny){
  int i, j;
  d1tod2(t, Nx, Ny, &i, &j);
  return neighbour_u(i, j, Nx, Ny);
}

int neighbour_d(int i, int j, int Nx, int Ny){
  return d2tod1(i, (Ny+j-1)%Ny, Nx, Ny);
}
//1d version
int neighbour_d(int t, int Nx, int Ny){
  int i, j;
  d1tod2(t, Nx, Ny, &i, &j);
  return neighbour_d(i, j, Nx, Ny);
}



// Spectral radius of an n X n matrix, assumed primitive
// Following: https://arxiv.org/pdf/1907.04175
double SpecRad(double **M, int n, double tol){
  double r[n];//row sums
  int i,j,k=0;
  int kmax=1000;//maximum iterations
  double **Mt=cpmat(M,n,n);
  double **Mt1=dmatrix(0,n-1,0,n-1);
  double error;
  
  //printmat(Mt,n,n);

  for(i=0;i<n;i++){
    r[i]=rowsum(Mt,n,n,i);
  }
  //printvec(r,n);
  error=vmax(r,n)-vmin(r,n);
  //fprintf(stderr, "Error = %.4f\n", error);

  while(error>tol && k<kmax){
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
	Mt1[i][j]=r[j]/r[i]*Mt[i][j];
      }
    }
    //printmat(Mt1,n,n);
    cpmat(Mt1,Mt,n,n);

    for(i=0;i<n;i++){
      r[i]=rowsum(Mt,n,n,i);
    }
    //printvec(r,n);
    error=vmax(r,n)-vmin(r,n);
    //fprintf(stderr, "Error = %.4f\n", error);
    k++;
  }
  if(k>=kmax && error>1e-3){
    fprintf(stderr, "max iterations reached in spectral radius calculation: error: %.4f.\n",error);
    free_dmatrix(Mt,0,n-1,0,n-1);
    free_dmatrix(Mt1,0,n-1,0,n-1);
    return -1;
  }
  
  free_dmatrix(Mt,0,n-1,0,n-1);
  free_dmatrix(Mt1,0,n-1,0,n-1);
  return 0.5*(vmax(r,n)+vmin(r,n));
}

//Compute the vector of extinction probabilities
//Athreya and Ney, V.3, Thm 2. (output in double s[])
//The Euclidean max-norm distance between two iterates must be < tol
void qvec(double ki, double kr, int *S, double *V, double **eps, int ncomp, double *s, double tol, int *iter){
  double L[ncomp],stmp[ncomp];
  double t1,t2,r1,r2,err;
  int i,j,flag=1,maxit=10000,c=0;
  for(i=0;i<ncomp;i++)//initialise s
    s[i]=0.1;

  while(flag){
    c++;
    //printvec(s,ncomp);
    if(c>maxit){
      fprintf(stderr, "ERROR in \"qvec\": took more than %d iterates to converge to the vector of extinction probabilities with tol %.7f. EXITING.\n",maxit,tol);
      exit(0);
    }
    for(i=0;i<ncomp;i++){
      t1=0.0;t2=0.0;for(j=0;j<ncomp;j++){r1=eps[i][j]*((double)(S[j]))/V[j];r2=r1*s[j];t1+=r1;t2+=r2;}
      L[i]=kr+ki*t1;
      stmp[i]=1.0/L[i]*(kr+ki*s[i]*t2);
    }
    //printvec(L,ncomp);
    err=0.0;
    for(i=0;i<ncomp;i++){
      err=max(err,abs(stmp[i]-s[i]));
    }
    cpvec(stmp,s,ncomp);
    if(err<tol){
      flag=0;
      *iter=c;
    }
    
  }
  return;
}

//overloading to accept a 2D array as input
//Not the most efficient for a sparsely coupled array, e.g., nn
//But easiest to code. Reconsider for large arrays.
void qvec(double ki, double kr, int **S, double **V, double **eps, int Nx, int Ny, double *s, double tol, int *iter){
  int Svec[Nx*Ny];
  double Vvec[Nx*Ny];
  int i,j,t;

  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      t=d2tod1(i,j,Nx,Ny);Svec[t]=S[i][j];Vvec[t]=V[i][j];
    }
  }

  qvec(ki, kr, Svec, Vvec, eps, Nx*Ny, s, tol, iter);
  return;
}

//outbreak probability assuming that the probability of introduction
//in a compartment is proportional to compartment volume. Equivalently,
//assuming all compartments have population density 1, so V stands in for
//compartment population
double outbreak_prob(double ki, double kr, int *S, double *V, double **eps, int ncomp, double *q, double tol, int *iter){
  int j;
  double Vtot=0, outbreak_prob=0.0;
  qvec(ki, kr, S, V, eps, ncomp, q, tol, iter);
  for(j=0;j<ncomp;j++){outbreak_prob+=V[j]*(1.0-q[j]);Vtot+=V[j];}
  outbreak_prob/=Vtot;
  return outbreak_prob;
}

//overloading to take matrix input
double outbreak_prob(double ki, double kr, int **S, double **V, double **eps, int Nx, int Ny, double *q, double tol, int *iter){
  int i,j,t;
  int ncomp=Nx*Ny;
  int Svec[ncomp];
  double Vvec[ncomp];
  double Vtot=0, outbreak_prob=0.0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      t=d2tod1(i,j,Nx,Ny);Svec[t]=S[i][j];Vvec[t]=V[i][j];
    }
  }
  qvec(ki, kr, Svec, Vvec, eps, ncomp, q, tol, iter);
  for(j=0;j<ncomp;j++){outbreak_prob+=Vvec[j]*(1.0-q[j]);Vtot+=Vvec[j];}
  outbreak_prob/=Vtot;
  return outbreak_prob;
}


//infection propensity: infector from comp. i, infectee in comp. j
double inf_prop(int *S, int *I, double *V, double **eps, double ki, int i, int j){
  return (double)(I[i]*S[j])*eps[i][j]*ki/V[j];
}

//infection propensities (overloading)
//2d arrays for S,I,V, but 1D indices i, j
double inf_prop(int **S, int **I, double **V, double **eps, double ki, int i, int j, int Nx, int Ny){
  int x1,y1,x2,y2;
  d1tod2(i, Nx, Ny, &x1, &y1);
  d1tod2(j, Nx, Ny, &x2, &y2);
  return (double)(I[x1][y1]*S[x2][y2])*eps[i][j]*ki/V[x2][y2];
}

//total infection propensity: infections occurring in comp. j
double inf_prop_tot(int *S, int *I, double *V, double **eps, double ki, int j, int ncomp){
  double tot=0.0;
  int i;
  for(i=0;i<ncomp;i++)
    tot+=(double)(I[i])*eps[i][j];
  tot*=ki*(double)(S[j])/V[j];
  return tot;
}

//compartmental SIR model: set the infection and recovery
//intensities in each compartment and return the total intensity
double SIRintensity(int *S, int *I, double *V, double **epsM, double ki, double kr, double *vr, double *vitot, int ncomp, int **sparse){
  int i,j;
  double vc,a0=0.0;
  if(sparse==NULL){
    for(i=0;i<ncomp;i++){
      vr[i]=(double)(I[i])*kr;a0+=vr[i];//recovery in comp. i
      for(j=0;j<ncomp;j++){//total infection intensity *into* comp. i
	vc=inf_prop(S, I, V, epsM, ki, j, i);//infection *into* compartment i, *from* compartment j: S_i+I_j = I_i + I_j
	a0+=vc;vitot[i]+=vc;
      }
    }
  }
  else{//sparse matrix
    for(i=0;i<ncomp;i++){
      vr[i]=(double)(I[i])*kr;a0+=vr[i];//recovery in comp. i
      for(j=1;j<1+sparse[i][0];j++){//total infection intensity *into* comp. i
	vc=inf_prop(S, I, V, epsM, ki, sparse[i][j], i);//infection *into* compartment i, *from* compartment j: S_i+I_j = I_i + I_j
	a0+=vc;vitot[i]+=vc;
      }
    }
  }
  return a0;
}

//compartmental SEIRS model: set the exposure, E-->I, I-->R and R-->S
//intensities in each compartment and return the total intensity
double SEIRSintensity(int *S, int *E, int *I, int *R, double *V, double **epsM, double ki, double kei, double kir, double krs, double *vei, double *vir, double *vrs, double *vitot, int ncomp, int **sparse){
  int i,j;
  double vc,a0=0.0;
  for(i=0;i<ncomp;i++){
    vir[i]=(double)(I[i])*kir;a0+=vir[i];//I-->R
    vei[i]=(double)(E[i])*kei;a0+=vei[i];//E-->I
    vrs[i]=(double)(R[i])*krs;a0+=vrs[i];//R-->S
  }
  if(sparse==NULL){
    for(i=0;i<ncomp;i++){
      for(j=0;j<ncomp;j++){//total infection intensity *into* comp. i
	vc=inf_prop(S, I, V, epsM, ki, j, i);//infection *into* compartment i, *from* compartment j: S_i+I_j = I_i + I_j
	a0+=vc;vitot[i]+=vc;
      }
    }
  }
  else{//sparse matrix
    for(i=0;i<ncomp;i++){
      for(j=1;j<1+sparse[i][0];j++){//total infection intensity *into* comp. i
	vc=inf_prop(S, I, V, epsM, ki, sparse[i][j], i);//infection *into* compartment i, *from* compartment j: S_i+I_j = I_i + I_j
	a0+=vc;vitot[i]+=vc;
      }
    }
  }
  return a0;
}

//The 2d grid version: because of the sparsity, it seems more
//efficient to do it this way
//Set the infection and recovery intensities in each compartment
//return the total intensity
double SIRintensity(int **S, int **I, double **V, double **epsM, double ki, double kr, double **vr, double **vitot, int Nx, int Ny, int **sparse){
  int i,j,k,i1;
  double vc,a0=0.0;
  inittozero(vitot,Nx,Ny);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      i1=d2tod1(i,j,Nx,Ny);
      vr[i][j]=(double)(I[i][j])*kr;a0+=vr[i][j];
      for(k=1;k<1+sparse[i1][0];k++){//total infection intensity *into* comp. i1
	vc=inf_prop(S, I, V, epsM, ki, sparse[i1][k], i1, Nx, Ny);
	a0+=vc;vitot[i][j]+=vc;
      }
    }
  }

  return a0;
}

//Get the approximate number of compartments affected
int episz(int *R, int *R0, int ncomp, int Nc, double *mean_outbreak_sz, double *mean_outbreak_sz_cond, int *totepis){
  int i,epsz=0;
  for(i=0;i<ncomp;i++){
    mean_outbreak_sz[i]+=(double)(R[i]-R0[i]);
    if(R[i]-R0[i]>Nc/10)//heuristic threshold for outbreak
      epsz++;
  }
  if(epsz>0){//conditioned Epi size
    (*totepis)++;
    for(i=0;i<ncomp;i++){
      mean_outbreak_sz_cond[i]+=(double)(R[i]-R0[i]);
    }
  }
  return epsz;
}

//overloading: 2d input
int episz(int **R, int **R0, int Nx, int Ny, int Nc, double **mean_outbreak_sz, double **mean_outbreak_sz_cond, int *totepis){
  int i,j,epsz=0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      mean_outbreak_sz[i][j]+=(double)(R[i][j]-R0[i][j]);
      if(R[i][j]-R0[i][j]>Nc/10)//heuristic threshold for outbreak
	epsz++;
    }
  }
  if(epsz>0){//conditioned Epi size
    (*totepis)++;
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	mean_outbreak_sz_cond[i][j]+=(double)(R[i][j]-R0[i][j]);
      }
    }
  }
  return epsz;
}



//Next Generation matrix
void NextG(double **NextGen, int *S, double *V, double **epsM, double ki, double kr, int ncomp){
  int i,j;
  for(i=0;i<ncomp;i++){
    for(j=0;j<ncomp;j++){
      NextGen[i][j]=epsM[i][j]*ki/kr*((double)S[j])/V[j];
    }
  }
}

//Overloading for 2d array: Next Generation matrix
void NextG(double **NextGen, int **S, double **V, double **epsM, double ki, double kr, int Nx, int Ny){
  int i,j;
  int x,y;
  int ncomp=Nx*Ny;
  for(i=0;i<ncomp;i++){
    for(j=0;j<ncomp;j++){
      d1tod2(j, Nx, Ny, &x, &y);
      //fprintf(stderr,"(x,y)=(%d,%d)\n",x,y);
      NextGen[i][j]=epsM[i][j]*ki/kr*((double)(S[x][y]))/V[x][y];
    }
  }
}

//classical Rt
double Rt_clas(int *S, double *V, int ncomp, int *Nc, double **epsM, double ki, double kr, int totpop){
  int i,j;
  double tot=0.0;
  for(i=0;i<ncomp;i++){
    for(j=0;j<ncomp;j++){
      tot+=epsM[i][j]*((double)Nc[i]*S[j])/V[j];
    }
  }
  tot*=ki/kr/((double)totpop);
  return tot;
}

//classical Rt - overloading for equal compartment populations
double Rt_clas(int *S, double *V, int ncomp, double **epsM, double ki, double kr){
  int i,j;
  double tot=0.0;
  for(i=0;i<ncomp;i++){
    for(j=0;j<ncomp;j++){
      tot+=epsM[i][j]*((double)S[j])/V[j];
    }
  }
  tot*=ki/kr/((double)ncomp);
  return tot;
}

//classical Rt for identical compartments
double Rt_clas(int Stot, double ki, double kr, int Ntot){
  return ki*((double)Stot)/kr/(double)Ntot;
}


//classical Rt - overloading for 2D array (identical compartments)
//with only nearest neighbour coupling
double Rt_clas_nn(int **S, double **V, double **epsM, double ki, double kr, int Nx, int Ny){
  int i,j,jj;
  double tot=0.0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      jj=d2tod1(i,j,Nx,Ny);
      tot+=epsM[jj][jj]*((double)S[i][j])/V[i][j];
      tot+=epsM[neighbour_r(i,j,Nx,Ny)][jj]*((double)S[i][j])/V[i][j];
      tot+=epsM[neighbour_l(i,j,Nx,Ny)][jj]*((double)S[i][j])/V[i][j];
      tot+=epsM[neighbour_u(i,j,Nx,Ny)][jj]*((double)S[i][j])/V[i][j];
      tot+=epsM[neighbour_d(i,j,Nx,Ny)][jj]*((double)S[i][j])/V[i][j];
    }
  }
  tot*=ki/kr/((double)(Nx*Ny));
  return tot;
}



// Basic output
void printtSIR(FILE *fd, double t, int *S, int *I, int *R, int ncomp){
  int i;
  fprintf(fd, "%.4f", t);
  for(i=0;i<ncomp;i++)
    fprintf(fd, ", %d, %d, %d", S[i],I[i],R[i]);
  fprintf(fd, "\n");
}

// overloading to get the Rt values
void printtSIR(FILE *fd, double t, double Rt_classical, double Rt_nextgen, int *S, int *I, int *R, int ncomp){
  int i,Stot=vsum(S,ncomp), Itot=vsum(I,ncomp), Rtot=vsum(R,ncomp);
  //to quietly produce nan if Rt_nextgen couldn't be calculated
  fprintf(fd, "%.4f %.4f %.4f %d %d %d", t, Rt_classical, Rt_nextgen>0?Rt_nextgen:0.0/0.0, Stot, Itot, Rtot);
  for(i=0;i<ncomp;i++)
    fprintf(fd, " %d %d %d", S[i],I[i],R[i]);
  fprintf(fd, "\n");
}

// A total of meantot time-points; dt is the dime increment
void printSIRmean(FILE *fd, int ncomp, int meantot, int totruns, double dt, int *S0, int *I0, int *R0, double **meanS, double **meanI, double **meanR){
  int i,j;
  //mean evolution
  fprintf(fd, "0.0");
  for(i=0;i<ncomp;i++)
    fprintf(fd, ", %.4f, %.4f, %.4f", (double)S0[i],(double)I0[i],(double)R0[i]);
  fprintf(fd, "\n");
	  
  for(j=1;j<meantot;j++){
    fprintf(fd, "%.4f", j/dt);
    for(i=0;i<ncomp;i++){
      fprintf(fd, ", %.4f, %.4f, %.4f", meanS[i][j-1]/((double)totruns),meanI[i][j-1]/((double)totruns), meanR[i][j-1]/((double)totruns));
    }
    fprintf(fd, "\n");
  }
  return;
}

//Quick and rough visual output of a histogram. Epihistnum is assumed to be a
//divisor of 100
void printstarhist(int *Hist, int histnum, int histmax, int maxlen){
  int i,j;
  for(i=0;i<=histnum;i++){
    if(i<histnum)
      fprintf(stderr, "%3d - %3d%%:", i*100/histnum, (i+1)*100/histnum);
    else
      fprintf(stderr, "      100%%:");
    for(j=0;j<Hist[i]*maxlen/histmax;j++){
      fprintf(stderr, "*");
    }
    fprintf(stderr, "\n");
  }
  return;
}

int SIRmeanupdate(int ncomp, int jstart, int jend, int *S, int *I, int *R, double **meanS, double **meanI, double **meanR){
  int i,j;
  for(j=jstart;j<jend;j++){
    for(i=0;i<ncomp;i++){
      meanS[i][j]+=S[i];meanI[i][j]+=I[i];meanR[i][j]+=R[i];
    }
  }
  return j;
}

int SEIRmeanupdate(int ncomp, int jstart, int jend, int *S, int *E, int *I, int *R, double **meanS, double **meanE, double **meanI, double **meanR){
  int i,j;
  for(j=jstart;j<jend;j++){
    for(i=0;i<ncomp;i++){
      meanS[i][j]+=S[i];meanE[i][j]+=E[i];meanI[i][j]+=I[i];meanR[i][j]+=R[i];
    }
  }
  return j;
}

//complete, symmetric, network
void sym_coupling(double eps, double **epsM, int ncomp){
  int i,j;
  for(i=0;i<ncomp;i++){
    for(j=0;j<ncomp;j++){
      if(i==j)
	epsM[i][j]=(1.0-eps);
      else
	epsM[i][j]=eps/((double)(ncomp-1));
    }
  }
}


//complete graph with gamma distributed coupling strengths with mean "leak"
//If "shp" is an integer, then sum of "shp" exponentials with rate "rate"
//With a fixed mean, decreasing shp (and hence increaseing scl) increases
//variance
//sym controls whether we start with symmetric coupling or not
//constout=0: constant total infectivity, variable in-infection,
// and constant *total* out-infection;
//constout=1: constant total infectivity, constant in-infection
// and constant *total* out-infection;
//constout=2: variable total infectivity, constant in-infection and
//constant out-infection on each edge
//constout=3: variable total infectivity, constant in-infection and
//constant susceptibility
void gamma_coupling(RNG &rng, double leak, double **epsM, int ncomp, int Nc, double shp, int sym, int constout){
  //double rate = shp/leak;
  double scl = leak/shp;
  double eps = leak/((double)Nc);
  double totepsrow[ncomp], totepscol[ncomp];
  int i, j;
  //static std::mt19937 rng(std::random_device{}()); 
  std::gamma_distribution<double> dist(shp, scl);
  inittoconst(totepsrow,ncomp,1.0-eps);
  inittoconst(totepscol,ncomp,1.0-eps);
  if(constout==0){
    for(i=0;i<ncomp;i++)
      epsM[i][i]=1.0;
    if(!sym){
      for(i=0;i<ncomp;i++){
	for(j=0;j<ncomp;j++){
	  if(i!=j){
	    epsM[i][j]=dist(rng)/(double)Nc/(double)(ncomp-1);
	    epsM[i][i]-=epsM[i][j];
	    if(epsM[i][i]<0.1){
	      fprintf(stderr, "Coupling too strong in \"gamma_coupling\". EXITING.\n");exit(0);
	    }
	  }
	}
      }
    }
    else{
      for(i=0;i<ncomp;i++){
	for(j=i+1;j<ncomp;j++){
	  epsM[i][j]=dist(rng)/(double)Nc/(double)(ncomp-1);
	  epsM[j][i]=epsM[i][j];
	  epsM[i][i]-=epsM[i][j];epsM[j][j]-=epsM[i][j];
	  if(epsM[i][i]<0.1 || epsM[j][j]<0.1){
	    fprintf(stderr, "Coupling too strong in \"gamma_coupling\". EXITING.\n");exit(0);
	  }
	}
      }
    }
  }
  else{//constout = 1 or 2 or 3
    for(i=0;i<ncomp;i++)
      epsM[i][i]=1.0-eps;
    if(!sym){
      for(i=0;i<ncomp;i++){
	for(j=0;j<ncomp;j++){
	  if(i!=j){
	    epsM[i][j]=dist(rng)/(double)Nc/(double)(ncomp-1);
	    totepsrow[i]+=epsM[i][j];
	    totepscol[j]+=epsM[i][j];
	  }
	}
      }
    }
    else{
      for(i=0;i<ncomp;i++){
	for(j=i+1;j<ncomp;j++){
	  epsM[i][j]=dist(rng)/(double)Nc/(double)(ncomp-1);
	  epsM[j][i]=epsM[i][j];
	  totepsrow[i]+=epsM[i][j];totepsrow[j]+=epsM[i][j];
	  totepscol[j]+=epsM[i][j];totepscol[i]+=epsM[i][j];
	}
      }
    }
  }

  if(constout==1){//rescale leaks to maintain constant total infectivity
    for(i=0;i<ncomp;i++){
      for(j=0;j<ncomp;j++){
	epsM[i][j]/=totepsrow[i];
      }
    }
  }
  else if(constout==3){//rescale leaks to maintain constant total susceptibility
    for(i=0;i<ncomp;i++){
      for(j=0;j<ncomp;j++){
	epsM[i][j]/=totepscol[i];
      }
    }
  }
  

  //printmat(epsM,ncomp,ncomp);exit(0);
  return;
}



//create epsM, given epsI and constout
//constout=0: constant total infectivity, variable in-infection,
// and constant *total* out-infection;
//constout=1: constant total infectivity, constant in-infection
// and constant *total* out-infection;
//constout=2: variable total infectivity, constant in-infection and
//constant out-infection on each edge
//constout=3: variable total infectivity, constant in-infection and
//constant susceptibility
void coupling_from_epsI(double eps, int **epsI, double **epsM, int ncomp, int constout){
  int tot=msum(epsI,ncomp,ncomp)-ncomp;//total (directed) edges
  double theta=(double)ncomp*eps/(double)tot;//eps per edge (eps/(mean degree))
  int i,j,t;
  //fprintf(stderr, "theta = %.4e\n", theta);

  if(constout<0 || constout > 3){
    fprintf(stderr, "ERROR in \"coupling_from_epsI\": constout must be an integer in [0,3]. EXITING.\n");exit(0);
  }
  
  if(constout==0){//variable in-infection
    for(i=0;i<ncomp;i++){//rows
      t=nonzentries(epsI[i],ncomp);//out-degree
      for(j=0;j<ncomp;j++){//cols
	if(epsI[i][j]){
	  if(i==j){
	    epsM[i][j]=1.0-(double)(t-1)*theta;
	    if(epsM[i][j]<0){
	      fprintf(stderr, "ERROR in \"coupling_from_epsI\". Inwards infection cannot be less than zero. EXITING.\n"); exit(0);
	    } 
	  }
	  else
	    epsM[i][j]=theta;
	}
      }
    }
  }
  else if(constout==1){//constant infectivity
    for(i=0;i<ncomp;i++){//rows
      t=nonzentries(epsI[i],ncomp);//out-degree
      for(j=0;j<ncomp;j++){//cols
	if(epsI[i][j]){
	  if(i==j)
	    epsM[i][j]=1.0-eps;
	  else
	    epsM[i][j]=eps/(double)(t-1);
	}
      }
    }
  }
  else if(constout==2){//variable susceptibility and infectivity
    for(i=0;i<ncomp;i++){//rows
      t=nonzentries(epsI[i],ncomp);//out-degree
      for(j=0;j<ncomp;j++){//cols
	if(epsI[i][j]){
	  if(i==j)
	    epsM[i][j]=1.0-eps;
	  else
	    epsM[i][j]=theta;
	}
      }
    }
  }
  else if(constout==3){//constant susceptibility
    for(i=0;i<ncomp;i++){//cols
      t=nonzentries(epsI,i,ncomp,ncomp);//in-degree
      for(j=0;j<ncomp;j++){//rows
	if(epsI[i][j]){//adjacent vertices
	  if(i==j)
	    epsM[j][i]=1.0-eps;
	  else
	    epsM[j][i]=eps/(double)(t-1);
	}
      }
    }
  }

    
  return;
}

//create sparse, given epsI
void sparse_from_epsI(int **epsI, int ncomp, int ***sparse){
  int i,j,t,t1;
  (*sparse) = (int**) malloc(sizeof(int*) * (ncomp));
  for(i=0;i<ncomp;i++){
    t=nonzentries(epsI[i],ncomp);
    (*sparse)[i]=(int*) malloc(sizeof(int) * (t+1));
    (*sparse)[i][0]=t;t1=1;
    for(j=0;j<ncomp;j++){
      if(epsI[i][j]){
	(*sparse)[i][t1++]=j;
      }
    }
  }


  return;
}



//star shaped network (node zero is the hub)
void star_coupling(int **epsI, int ncomp){
  int i,j;

  for(i=0;i<ncomp;i++){
    for(j=0;j<ncomp;j++){
      if(i==j){
	epsI[i][j]=1;
      }
      else if(i==0 || j==0){
	epsI[i][j]=1;
      }
    }
  }

  /* //set epsM and sparse */
  /* coupling_from_epsI(eps, epsI, epsM, ncomp, sparse, constout); */
  return;
}



//2D grid with some small-world rewiring
//For meaning of constout, see above
int nn_smallworld(RNG &rng, int **epsI, int Nx, int Ny, double rewiring_prob){
  int i,j,t1,t2;
  int ncomp=Nx*Ny;
  double r1,r2;
  //static std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> runif(0.0, 1.0);
  inittozero(epsI,Nx*Ny,Nx*Ny);

  //basic coupling structure
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      t1=d2tod1(i, j, Nx, Ny);epsI[t1][t1]=1;
      t2=neighbour_r(i,j,Nx,Ny);epsI[t1][t2]=1;
      t2=neighbour_l(i,j,Nx,Ny);epsI[t1][t2]=1;
      t2=neighbour_u(i,j,Nx,Ny);epsI[t1][t2]=1;
      t2=neighbour_d(i,j,Nx,Ny);epsI[t1][t2]=1;
    }
  }

  //now do the rewiring
  for(i=0;i<ncomp;i++){
    r1 = runif(rng);
    if(r1<rewiring_prob){//rewire right neighbour
      do{r2 = runif(rng);}while((int)(ncomp*r2)>=ncomp || epsI[i][(int)(ncomp*r2)]);//find a non-edge
      epsI[i][(int)(ncomp*r2)]=1;epsI[(int)(ncomp*r2)][i]=1;//add
      t2=neighbour_r(i,Nx,Ny);
      epsI[i][t2]=0;epsI[t2][i]=0;//delete
    }
    r1 = runif(rng);
    if(r1<rewiring_prob){//rewire down neighbour
      do{r2 = runif(rng);}while((int)(ncomp*r2)>=ncomp || epsI[i][(int)(ncomp*r2)]);//find a non-edge
      epsI[i][(int)(ncomp*r2)]=1;epsI[(int)(ncomp*r2)][i]=1;//add
      t2=neighbour_d(i,Nx,Ny);
      epsI[i][t2]=0;epsI[t2][i]=0;//delete
    }
  }

  /* //set epsM and sparse */
  /* coupling_from_epsI(eps, epsI, epsM, ncomp, sparse, constout); */

  return is_SC(epsI,ncomp);//strongly connected?
}



//coupling matrix: total out-infection fraction = eps
//K is the mean degree, assumed even. We expect K << ncomp
//rewiring_prob is the "randomness" parameter - likelihood of rewiring
//see https://en.wikipedia.org/wiki/Watts%E2%80%93Strogatz_model
//For meaning of constout, see above
int small_world(RNG &rng, int **epsI, int ncomp, int K, double rewiring_prob){
  int i,j;
  double r1, r2;
  //static std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> runif(0.0, 1.0);
  //fprintf(stderr, "entering small_world\n");

  if(K<=0 || K%2!=0 || K>ncomp-2){
    fprintf(stderr, "for a (ring) small-world network the initial clique size must be a positive, even integer (currently %d) and must be less than the number of compartments (currently %d) minus 1\n\nEXITING\n", K, ncomp);exit(0);
  }

  //Set up regular ring lattice
  inittozero(epsI,ncomp,ncomp);
  for(i=0;i<ncomp;i++){
    for(j=0;j<ncomp;j++){
      if(i==j)
	epsI[i][j]=1;//(1.0-eps);
      else if(abs(j-i)%(ncomp-1-K/2) > 0 && abs(j-i)%(ncomp-1-K/2) <= K/2)
	epsI[i][j]=1;//eps/((double)(K));
    }
  }

  //now do the rewiring
  for(i=0;i<ncomp;i++){
    for(j=i+1;j<i+1+K/2;j++){
      r1 = runif(rng);
      if(r1<rewiring_prob){//rewire
	do{r2 = runif(rng);}while((int)(ncomp*r2)>=ncomp || epsI[i][(int)(ncomp*r2)]);//find a non-edge
	epsI[i][(int)(ncomp*r2)]=1;//epsM[i][j%ncomp];
	epsI[(int)(ncomp*r2)][i]=1;//epsM[i][j%ncomp];
	epsI[i][j%ncomp]=0;epsI[j%ncomp][i]=0;
      }
    }
  }
  
  return is_SC(epsI,ncomp);//strongly connected?
}

//Random (symmetric) graph G(n,p)
//https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model
int randomgraph(RNG &rng, int **epsI, int ncomp, double meandeg){
  int i,j;
  double r1;
  double prob = double(meandeg)/(double)(ncomp-1);
  std::uniform_real_distribution<> runif(0.0, 1.0);
  
  if(prob < 0 || prob>1){
    fprintf(stderr, "ERROR in randomgraph: edges must be chosen with probability in [0,1]. EXITING.\n");exit(0);
  }


  inittozero(epsI,ncomp,ncomp);
  for(i=0;i<ncomp;i++){
    for(j=i;j<ncomp;j++){
      if(i==j)
	epsI[i][j]=1;//(1.0-eps);
      else{
	r1 = runif(rng);
	if(r1<prob){
	  epsI[i][j]=1;
	  epsI[j][i]=1;
	}
      }
    }
  }
  
  return is_SC(epsI,ncomp);//strongly connected?
}


//preferential attachment starting with a regular ring lattice of size m0 and clique-size K
//after this we attach m1<=K<=m0-1 edges at each stage. K must be even
void ring_scalefree(RNG &rng, int **epsI, int ncomp, int m0, int K, int m1){
  int i,j,R1,cumul,totdegree=0,totd,newedges;
  double r1;
  std::uniform_real_distribution<> runif(0.0, 1.0);

  if(K<=0 || K%2!=0){
    fprintf(stderr, "For a (ring) preferential attachment network the initial clique size must be a positive, even integer (currently %d)\n\nEXITING\n", K);exit(0);
  }
  if(m1>K || K>m0-2){
    fprintf(stderr, "For a (ring) preferential attachment network the edges to attach at each stage (currently %d) must be less than or equal to the initial clique size (currently %d), which must be less than the initial ring (currently %d) minus 1\n\nEXITING\n", m1, K, m0);exit(0);
  }

  //Set up regular ring lattice (without loops)
  inittozero(epsI,ncomp,ncomp);
  for(i=0;i<m0;i++){
    for(j=i+1;j<m0;j++){
      if(abs(j-i)%(m0-1-K/2) > 0 && abs(j-i)%(m0-1-K/2) <= K/2){
	epsI[i][j]=1;epsI[j][i]=1;totdegree+=2;
      }
    }
  }

  for(i=m0;i<ncomp;i++){
    //fprintf(stderr, "i=%d\n", i);
    totd=totdegree;newedges=0;
    while(newedges<m1){//m1 new edges
      do{r1 = runif(rng);}while(r1==1);
      R1=(int)(r1*totd);cumul=0;
      for(j=0;j<i;j++){
	cumul+=vsum(epsI[j],ncomp);
	if(R1<=cumul && epsI[i][j]==0){//nonedge: attach to j
	  epsI[i][j]=1;totdegree+=2;newedges++;break;
	}
      }
    }
    //attach the m1 reverse edges (don't do this earlier - messes up cumul)
    for(j=0;j<i;j++){
      if(epsI[i][j]==1)
	epsI[j][i]=1;
    }
  }

  for(i=0;i<ncomp;i++)//loops at the end
    epsI[i][i]=1;
  
  return; //strongly connected?
}



//Preferential attachment algorithm to create a network with parameter m0
//(initial clique of m0 nodes, and each new node attaches to m0 existing ones)
//For meaning of constout, see above
void scale_free(RNG &rng, int **epsI, int ncomp, int m0){
  int i,j,R1,cumul,totdegree=0,totd,newedges;
  double r1;
  //static std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> runif(0.0, 1.0);

  if(m0>ncomp){
    fprintf(stderr, "ERROR in scale_free: initial clique (%d) cannot be larger than the number of compartments (%d). EXITING.\n",m0,ncomp);exit(0);
  }
  
  //fprintf(stderr, "entering scale_free\n");

  //leave loops for the end
  inittozero(epsI,ncomp,ncomp);
  for(i=0;i<m0;i++){//initial clique
    for(j=i+1;j<m0;j++){
      epsI[i][j]=1;epsI[j][i]=1;totdegree+=2;
    }
  }

  
  for(i=m0;i<ncomp;i++){
    //fprintf(stderr, "i=%d\n", i);
    totd=totdegree;newedges=0;
    while(newedges<m0){//m0 new edges
      do{r1 = runif(rng);}while(r1==1);
      R1=(int)(r1*totd);cumul=0;
      for(j=0;j<i;j++){
	cumul+=vsum(epsI[j],ncomp);
	if(R1<=cumul && epsI[i][j]==0){//nonedge: attach to j
	  epsI[i][j]=1;totdegree+=2;newedges++;break;
	}
      }
    }
    //attach the m0 reverse edges (don't do this earlier - messes up cumul)
    for(j=0;j<i;j++){
      if(epsI[i][j]==1)
	epsI[j][i]=1;
    }
  }
  
  for(i=0;i<ncomp;i++)//loops at the end
    epsI[i][i]=1;
    
  return; //definitely strongly connected
}



//Run an SIR epidemic. Terminate when there are no more infecteds
//or time tmax is reached
double runSIRepi(RNG &rng, int ncomp, double *V, int *I0, int *S0, int *R0, int *I, int *S, int *R, double **meanI, double **meanS, double **meanR, double **epsM, double ki, double kr, double tmax, int printfull, int printmean, int meantot, double dt, FILE *fd, int **sparse){
  int i;
  double t=0.0,tau=0.0;
  int jlast=0;
  double vr[ncomp],vitot[ncomp];//intensities
  double a0,r1, r2, cumul;
  int Ntot;
  double Rt_classical, Rt_nextgen;
  double **NextGen=dmatrix(ncomp,ncomp);
  std::uniform_real_distribution<> runif(0.0, 1.0);

  //initialise
  for(i=0;i<ncomp;i++){I[i]=I0[i];S[i]=S0[i];R[i]=R0[i];}

  if(printfull){
    //assuming identical compartments
    Ntot=vsum(I0,ncomp)+vsum(S0,ncomp)+vsum(R0,ncomp);
    Rt_classical=Rt_clas(vsum(S,ncomp), ki, kr, Ntot);
    NextG(NextGen, S, V, epsM, ki, kr, ncomp);
    Rt_nextgen=SpecRad(NextGen,ncomp,1e-7);	
    printtSIR(fd, t, Rt_classical, Rt_nextgen, S, I, R, ncomp);
  }

      
  while(t<(int)tmax+10){
    inittozero(vitot,ncomp);
    a0=SIRintensity(S, I, V, epsM, ki, kr, vr, vitot, ncomp, sparse);

    do{r1 = runif(rng);}while(r1==0);tau=1.0/a0*log(1.0/r1);t+=tau;
    r2 = runif(rng);
    if(printfull){
      //assuming identical compartments
      Rt_classical=Rt_clas(vsum(S,ncomp), ki, kr, Ntot);
      NextG(NextGen, S, V, epsM, ki, kr, ncomp);
      Rt_nextgen=SpecRad(NextGen,ncomp,1e-7);	
      printtSIR(fd, t, Rt_classical, Rt_nextgen, S, I, R, ncomp);
    }
    
    if(printmean)
      jlast=SIRmeanupdate(ncomp, min(jlast,meantot), min(meantot-1,(int)(dt*t)), S, I, R, meanS, meanI, meanR);

    //Gillespie step
    cumul=0.0;
    for(i=0;i<ncomp;i++){//compartment i
      if(r2>=cumul/a0){
	if(r2<(cumul+vitot[i])/a0){//infection
	  S[i]--;I[i]++;break;
	}
	else if(r2<(cumul+vitot[i]+vr[i])/a0){//recovery
	  I[i]--;R[i]++;break;
	}
      }
      cumul+=vitot[i]+vr[i];
    }

    if(printfull){//assuming identical compartments
      Rt_classical=Rt_clas(vsum(S,ncomp), ki, kr, Ntot);
      NextG(NextGen, S, V, epsM, ki, kr, ncomp);
      Rt_nextgen=SpecRad(NextGen,ncomp,1e-7);	
      printtSIR(fd, t, Rt_classical, Rt_nextgen, S, I, R, ncomp);
    }
	
    if(vsum(I,ncomp)==0){//epidemic over
      if(printmean)
	SIRmeanupdate(ncomp, min(jlast,meantot), meantot-1, S, I, R, meanS, meanI, meanR);
      free_dmatrix(NextGen,ncomp,ncomp);
      return t;
    }
  }
  
  free_dmatrix(NextGen,ncomp,ncomp);
  return -1;
}

//Run an SEIRS epidemic. Terminate when there are no more infecteds
//or time tmax is reached
//simplified from runSIRepi
double runSEIRSepi(RNG &rng, int ncomp, double *V, int *I0, int *S0, int *E0, int *R0, int *I, int *S, int *E, int *R, double **meanI, double **meanS, double **meanE, double **meanR, double **epsM, double ki, double kei, double kir, double krs, double tmax, int printmean, int meantot, double dt, FILE *fd, int **sparse){

  int i;
  double t=0.0,tau=0.0;
  int jlast=0;
  double vei[ncomp],vir[ncomp],vrs[ncomp],vitot[ncomp];//intensities
  double a0,r1, r2, cumul;
  /* int Ntot; */
  /* double Rt_classical, Rt_nextgen; */
  /* double **NextGen=dmatrix(ncomp,ncomp); */
  std::uniform_real_distribution<> runif(0.0, 1.0);

  //initialise
  for(i=0;i<ncomp;i++){S[i]=S0[i];E[i]=E0[i];I[i]=I0[i];R[i]=R0[i];}

  /* if(printfull){ */
  /*   //assuming identical compartments */
  /*   Ntot=vsum(I0,ncomp)+vsum(S0,ncomp)+vsum(R0,ncomp); */
  /*   Rt_classical=Rt_clas(vsum(S,ncomp), ki, kr, Ntot); */
  /*   NextG(NextGen, S, V, epsM, ki, kr, ncomp); */
  /*   Rt_nextgen=SpecRad(NextGen,ncomp,1e-7);	 */
  /*   printtSIR(fd, t, Rt_classical, Rt_nextgen, S, I, R, ncomp); */
  /* } */

      
  while(t<(int)tmax+10){
    inittozero(vitot,ncomp);
    a0=SEIRSintensity(S, E, I, R, V, epsM, ki, kei, kir, krs, vei, vir, vrs, vitot, ncomp, sparse);

    do{r1 = runif(rng);}while(r1==0);tau=1.0/a0*log(1.0/r1);t+=tau;
    r2 = runif(rng);
    /* if(printfull){ */
    /*   //assuming identical compartments */
    /*   Rt_classical=Rt_clas(vsum(S,ncomp), ki, kr, Ntot); */
    /*   NextG(NextGen, S, V, epsM, ki, kr, ncomp); */
    /*   Rt_nextgen=SpecRad(NextGen,ncomp,1e-7);	 */
    /*   printtSIR(fd, t, Rt_classical, Rt_nextgen, S, I, R, ncomp); */
    /* } */
    
    if(printmean)
      jlast=SEIRmeanupdate(ncomp, min(jlast,meantot), min(meantot-1,(int)(dt*t)), S, E, I, R, meanS, meanE, meanI, meanR);


    //Gillespie step
    cumul=0.0;
    for(i=0;i<ncomp;i++){//compartment i
      if(r2>=cumul/a0){
	if(r2<(cumul+vitot[i])/a0){//exposure
	  S[i]--;E[i]++;break;
	}
	else if(r2<(cumul+vitot[i]+vei[i])/a0){//E->I
	  E[i]--;I[i]++;break;
	}
	else if(r2<(cumul+vitot[i]+vei[i]+vir[i])/a0){//I-->R
	  I[i]--;R[i]++;break;
	}
	else if(r2<(cumul+vitot[i]+vei[i]+vir[i]+vrs[i])/a0){//R-->S
	  R[i]--;S[i]++;break;
	}
      }
      cumul+=vitot[i]+vei[i]+vir[i]+vrs[i];
    }

    /* if(printfull){//assuming identical compartments */
    /*   Rt_classical=Rt_clas(vsum(S,ncomp), ki, kr, Ntot); */
    /*   NextG(NextGen, S, V, epsM, ki, kr, ncomp); */
    /*   Rt_nextgen=SpecRad(NextGen,ncomp,1e-7);	 */
    /*   printtSIR(fd, t, Rt_classical, Rt_nextgen, S, I, R, ncomp); */
    /* } */
	
    if(vsum(I,ncomp)==0){//epidemic over
      if(printmean)
	jlast=SEIRmeanupdate(ncomp, min(jlast,meantot), meantot-1, S, E, I, R, meanS, meanE, meanI, meanR);
      /* free_dmatrix(NextGen,ncomp,ncomp); */
      return t;
    }
  }
  
  /* free_dmatrix(NextGen,ncomp,ncomp); */
  return -1;
}


//print a matrix scaled by some expected maximum
void printscaledmat(FILE *fd, int **A, int N, int Nx, int Ny){
  int i,j;
  for(i=0;i<Ny;i++){
    for(j=0;j<Nx;j++){
      fprintf(fd, "%.4f ", (double)(A[j][i])/(double)(N));
    }
    fprintf(fd, "\n");
  }
  fprintf(fd, "\n");
}

//overloading: vector representation of the matrix A
void printscaledmat(FILE *fd, int *A, int N, int Nx, int Ny){
  int i,j;
  for(i=0;i<Ny;i++){
    for(j=0;j<Nx;j++){
      fprintf(fd, "%.4f ", (double)(A[d2tod1(j, i, Nx, Ny)])/(double)(N));
    }
    fprintf(fd, "\n");
  }
  fprintf(fd, "\n");
}



void printSImatrices(const char basedir[], int **I, int **S, int Nx, int Ny, double Icomp, double Scomp, int tlast){
  FILE *fd;
  char fname[strlen(basedir)+10];
  sprintf(fname, "%s/dataI/m%04d",basedir,tlast);
  fd=fopen(fname, "w");
  printscaledmat(fd, I, Icomp, Nx, Ny);//print histogram (matrix)
  fclose(fd);
  sprintf(fname, "%s/dataS/m%04d",basedir,tlast);
  fd=fopen(fname, "w");
  printscaledmat(fd, S, Scomp, Nx, Ny);//print histogram (matrix)
  fclose(fd);
}

//overloading: vector representation of the matrices I and S
void printSImatrices(const char basedir[], int *I, int *S, int Nx, int Ny, double Icomp, double Scomp, int tlast){
  FILE *fd;
  char fname[strlen(basedir)+10];
  sprintf(fname, "%s/dataI/m%04d",basedir,tlast);
  fd=fopen(fname, "w");
  printscaledmat(fd, I, Icomp, Nx, Ny);//print histogram (matrix)
  fclose(fd);
  sprintf(fname, "%s/dataS/m%04d",basedir,tlast);
  fd=fopen(fname, "w");
  printscaledmat(fd, S, Scomp, Nx, Ny);//print histogram (matrix)
  fclose(fd);
}

//overloading
//The 2d grid case
double runSIRepi(RNG &rng, int Nx, int Ny, double **V, int **I0, int **S0, int **R0, int **I, int **S, int **R, int Iinit, double Icomp, double Scomp, double **epsM, double ki, double kr, double tmax, int sampling, int printall, FILE *fd4, FILE *fd1, int **sparse){
  int i,j,tlast,flag,Isum,Imax;
  double t=0.0, tau=0.0, cumul, r1, r2, a0;
  double **vr, **vitot;//intensities
  int maxpts=(int)((tmax+10)*sampling);
  int Isumtot[maxpts];
  std::uniform_real_distribution<> runif(0.0, 1.0);
  tlast=0;

  vitot=dmatrix(0, Nx-1, 0, Ny-1);
  vr=dmatrix(0, Nx-1, 0, Ny-1);

  //initialise
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      I[i][j]=I0[i][j];S[i][j]=S0[i][j];R[i][j]=R0[i][j];
    }
  }

  if(printall){
    printscaledmat(fd4, I, Icomp, Nx, Ny);
    printSImatrices("SIRtorus", I, S, Nx, Ny, Icomp, Scomp, tlast);
    inittozero(Isumtot,maxpts);Isumtot[0]=Iinit;Imax=Iinit;
  }

    
  while(t<tmax+10){
    while(tlast <(int)(t*sampling + 0.0001)){
      tlast++;
      if(printall){
	printSImatrices("SIRtorus", I, S, Nx, Ny, Icomp, Scomp, tlast);
	Isum=msum(I,Nx,Ny);Isumtot[tlast]=Isum;Imax=max(Imax,Isum);
      }
    }

    inittozero(vitot,Nx,Ny);
    a0=SIRintensity(S, I, V, epsM, ki, kr, vr, vitot, Nx, Ny, sparse);

    do{r1 = runif(rng);}while(r1==0);tau=1.0/a0*log(1.0/r1);t+=tau;
    r2 = runif(rng);

    if(printall)
      printscaledmat(fd4, I, Icomp, Nx, Ny);
	
    cumul=0.0;flag=0;
    for(i=0;i<Nx;i++){//compartment i,j
      for(j=0;j<Ny;j++){
	//fprintf(stderr, "r2=%.4f, cumul = %.4f, a0: %.4f, vitot: %.4f, vr: %.4f\n", r2, cumul/a0,a0,vitot[i][j]/a0,vr[i][j]/a0);
	if(r2>=cumul/a0){
	  if(r2<(cumul+vitot[i][j])/a0){//infection
	    S[i][j]--;I[i][j]++;flag=1;break;
	    //fprintf(stderr, "I:r2=%.4f, cumul = %.4f, a0: %.4f, vitot: %.4f, vr: %.4f\n", r2, cumul/a0,a0,vitot[i][j]/a0,vr[i][j]/a0);
	  }
	  else if(r2<(cumul+vitot[i][j]+vr[i][j])/a0){//recovery
	    I[i][j]--;R[i][j]++;flag=1;break;
	    //fprintf(stderr, "R:r2=%.4f, cumul = %.4f, a0: %.4f, vitot: %.4f, vr: %.4f\n", r2, cumul/a0,a0,vitot[i][j]/a0,vr[i][j]/a0);
	  }
	}
	cumul+=vitot[i][j]+vr[i][j];
      }
      if(flag==1)
	break;
    }
      
    if(printall)
      printscaledmat(fd4, I, Icomp, Nx, Ny);    
           
    Isum=msum(I,Nx,Ny);
    
    //if(msum(R,Nx,Ny)>10000){exit(0);}
    if(Isum==0){//epidemic over
      //fprintf(stderr, "got hereC\n");fflush(stderr);
      //fprintf(stderr, "Rsum=%d\n", msum(R,Nx,Ny));fflush(stderr);
      //if(msum(R,Nx,Ny)<2000){exit(0);}
      if(printall){
	//padding out
	while(tlast<(int)(sampling*t + 110)){
	  tlast++;
	  printSImatrices("SIRtorus", I, S, Nx, Ny, Icomp, Scomp, tlast);
	}
	fprintf(fd4,"\n\n");
	for(i=0;i<tlast;i++)
	  fprintf(fd1, "%.4f %.4f\n", (double)i/(double)sampling, (double)Isumtot[i]/(double)Imax);
	fprintf(fd1, "\n");
      }

      free_dmatrix(vitot, 0, Nx-1, 0, Ny-1);
      free_dmatrix(vr, 0, Nx-1, 0, Ny-1);
      return t;
    }
  }

  free_dmatrix(vitot, 0, Nx-1, 0, Ny-1);
  free_dmatrix(vr, 0, Nx-1, 0, Ny-1);

  return -1.0; //Epidemic didn't end

}


//numsims is the number of simulations
//Assume no infecteds, assume identical compartments
//choose_by_S non-zero means choose introductions according to
//the *susceptible* population (rather than the total population)
double mean_SIR_outbreak(RNG &rng, int ncomp, double *V, int *S0, int *R0, double **epsM, double ki, double kr, double tmax, int numsims, int choose_by_S, int **sparse){
  int I[ncomp],S[ncomp],R[ncomp],I0[ncomp];
  double t;
  int i,i0,j,r2;
  double mean_sz=0;
  int Nc[ncomp];
  int pop_tot=0;
  int totsims=0;
  std::uniform_real_distribution<> runif(0.0, 1.0);
  int Scum[ncomp+1];
  int alleq=1;//compartments of equal size
  int S0tot=0, R0tot=0;
  //fprintf(stderr, "initial S state: "); printvec(S0,ncomp);
  //fprintf(stderr, "initial R state: "); printvec(R0,ncomp);

  Scum[0]=0;

  for(i=0;i<ncomp;i++){
    Nc[i]=S0[i]+R0[i];S0tot+=S0[i];R0tot+=R0[i];pop_tot+=Nc[i];
    if(i>0&&Nc[i]!=Nc[0])
      alleq=0;
    if(choose_by_S)
      Scum[i+1]=Scum[i]+S0[i];
    else
      Scum[i+1]=Scum[i]+Nc[i];
  }
  
  if(S0tot==0){//no susceptibles
    fprintf(stderr,"WARNING in \"mean_SIR_outbreak\": empty susceptible population, so no epidemic is possible.\n");
    return 0.0;
  }
  
  inittozero(I0,ncomp);
  j=0;
  while(j<numsims){
    if(!choose_by_S){//choose an individual
      do{r2 = (int)(runif(rng)*pop_tot);}while(r2==pop_tot);
    }
    else{
      do{r2 = (int)(runif(rng)*S0tot);}while(r2==S0tot);
    }

    
    //in which compartment do they reside? (i0)
    if(!choose_by_S && alleq)
      i0=r2%ncomp;
    else{
      for(i=0;i<ncomp;i++){
	if(r2>Scum[i]&&r2<=Scum[i+1] &&S0[i]>0){//initial infection in comp. i
	  i0=i;break;
	}
      }
    }


    //introduce infection into i0
    S0[i0]--;I0[i0]++;totsims++;
    
    if(totsims%100==0)//to track progress
      fprintf(stderr, "*");
    if((t=runSIRepi(rng, ncomp, V, I0, S0, R0, I, S, R, NULL, NULL, NULL, epsM, ki, kr, tmax, 0, 0, 0, 0, NULL, sparse))<0.0){
      fprintf(stderr,"ERROR in \"mean_SIR_outbreak\": epidemic duration longer than tmax. EXITING.\n");exit(0);
    }
    mean_sz+=vsum(R,ncomp)-R0tot;
    //fprintf(stderr, "%.4f\n", (double)(vsum(R,ncomp)-R0tot)/(double)pop_tot);
    j++;S0[i0]++;I0[i0]--;
    
  }

  fprintf(stderr, "\n");
  //fprintf(stderr, "totsims = %d\n", totsims);
  return ((double)mean_sz)/((double)totsims)/((double)pop_tot);
}

//2d grid version
double mean_SIR_outbreak(RNG &rng, int Nx, int Ny, double **V, int **S0, int **R0, double **epsM, double ki, double kr, double tmax, int numsims, int **sparse){
  int **I, **S, **R, **I0;
  double t;
  int i,j,jsim;
  double mean_sz=0;
  int Nc[Nx][Ny];
  int pop_tot=0;
  int sims,totsims=0;

  I=imatrix(0, Nx-1, 0, Ny-1);
  S=imatrix(0, Nx-1, 0, Ny-1);
  R=imatrix(0, Nx-1, 0, Ny-1);
  I0=imatrix(0, Nx-1, 0, Ny-1);

  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      Nc[i][j]=S0[i][j]+R0[i][j];
      pop_tot+=Nc[i][j];
    }
  }

  inittozero(I0,Nx,Ny);

  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      sims=(int)(((double)Nc[i][j])/((double)pop_tot)*((double)numsims));
      I0[i][j]++;S0[i][j]--;
      jsim=0;
      //fprintf(stderr, "i=%d\n", i);
      while(jsim<sims){
	totsims++;
	if(totsims%100==0)//to track progress
	  fprintf(stderr, "*");
	//run the epidemic
	if((t=runSIRepi(rng, Nx, Ny, V, I0, S0, R0, I, S, R, 1, 0, 0, epsM, ki, kr, tmax, 1, 0, NULL, NULL, sparse))<0.0){
	  fprintf(stderr,"Epidemic duration longer than tmax. EXITING.\n");exit(0);
	}
	mean_sz+=msum(R,Nx,Ny)-msum(R0,Nx,Ny);
	jsim++;
      }
      I0[i][j]--;S0[i][j]++;
    }
  }
  fprintf(stderr, "\n");
  //fprintf(stderr, "totsims = %d\n", totsims);
  free_imatrix(I, 0, Nx-1, 0, Ny-1);
  free_imatrix(S, 0, Nx-1, 0, Ny-1);
  free_imatrix(R, 0, Nx-1, 0, Ny-1);
  free_imatrix(I0, 0, Nx-1, 0, Ny-1);
  return ((double)mean_sz)/((double)totsims)/((double)pop_tot);
}

//numsims is the number of simulations
//Assume no E or I, assume identical compartments
//choose_by_S non-zero means choose introductions according to
//the *susceptible* population (rather than the total population)
 double mean_SEIRS_outbreak(RNG &rng, int ncomp, double *V, int *S0, int *R0, double **epsM, double ki, double kei, double kir, double krs, double tmax, int numsims, int choose_by_S, int **sparse){
  int I[ncomp],E[ncomp],S[ncomp],R[ncomp],I0[ncomp],E0[ncomp];
  double t;
  int i,i0,j,r2;
  double mean_sz=0;
  int Nc[ncomp];
  int pop_tot=0;
  int totsims=0;
  std::uniform_real_distribution<> runif(0.0, 1.0);
  int Scum[ncomp+1];
  int alleq=1;//compartments of equal size
  int S0tot=0, R0tot=0;
  //fprintf(stderr, "initial S state: "); printvec(S0,ncomp);
  //fprintf(stderr, "initial R state: "); printvec(R0,ncomp);

  Scum[0]=0;

  for(i=0;i<ncomp;i++){
    Nc[i]=S0[i]+R0[i];S0tot+=S0[i];R0tot+=R0[i];pop_tot+=Nc[i];
    if(i>0&&Nc[i]!=Nc[0])
      alleq=0;
    if(choose_by_S)
      Scum[i+1]=Scum[i]+S0[i];
    else
      Scum[i+1]=Scum[i]+Nc[i];
  }
  
  if(S0tot==0){//no susceptibles
    fprintf(stderr,"WARNING in \"mean_SIR_outbreak\": empty susceptible population, so no epidemic is possible.\n");
    return 0.0;
  }
  
  inittozero(I0,ncomp);inittozero(E0,ncomp);
  j=0;
  while(j<numsims){
    if(!choose_by_S){//choose an individual
      do{r2 = (int)(runif(rng)*pop_tot);}while(r2==pop_tot);
    }
    else{
      do{r2 = (int)(runif(rng)*S0tot);}while(r2==S0tot);
    }
    //in which compartment do they reside? (i0)
    if(!choose_by_S && alleq)
      i0=r2%ncomp;
    else{
      for(i=0;i<ncomp;i++){
	if(r2>Scum[i]&&r2<=Scum[i+1] &&S0[i]>0){//initial infection in comp. i
	  i0=i;break;
	}
      }
    }
    //introduce infection into i0
    S0[i0]--;I0[i0]++;totsims++;
    if(totsims%100==0)//to track progress
      fprintf(stderr, "*");
    if((t=runSEIRSepi(rng, ncomp, V, I0, S0, E0, R0, I, S, E, R, NULL, NULL, NULL, NULL, epsM, ki, kei, kir, krs, tmax, 0, 0, 0, NULL, sparse))<0.0){
      fprintf(stderr,"ERROR in \"mean_SEIRS_outbreak\": epidemic duration longer than tmax. EXITING.\n");exit(0);
    }
    mean_sz+=vsum(R,ncomp)-R0tot;
    //fprintf(stderr, "%.4f\n", (double)(vsum(R,ncomp)-R0tot)/(double)pop_tot);
    j++;S0[i0]++;I0[i0]--;
    
  }

  fprintf(stderr, "\n");
  //fprintf(stderr, "totsims = %d\n", totsims);
  return ((double)mean_sz)/((double)totsims)/((double)pop_tot);
}

 

//randomly choose susceptible individuals to vaccinate with probabilities
//proportional to susceptible populations
//Note: this is different from random vaccination
//of possibly recovered individuals
//Note that frac refers to the fraction of *susceptible* individuals
void vaccinate(RNG &rng, int ncomp, int *S, int *R, double frac){
  int i,vaxed,r2;
  int Stot=vsum(S,ncomp);
  int tovax=(int)(frac*((double)Stot));
  std::uniform_real_distribution<> runif(0.0, 1.0);
  int Scum[ncomp+1];

  Scum[0]=0;
  for(i=0;i<ncomp;i++)
    Scum[i+1]=Scum[i]+S[i];

  vaxed=0;
  while(vaxed<tovax){
    r2 = (int)(runif(rng)*Stot);
    for(i=0;i<ncomp;i++){
      if(r2>Scum[i]&&r2<=Scum[i+1] &&S[i]>0){
	S[i]--;R[i]++;vaxed++;
	break;
      }
    }
    
  }
  return;

}

//overloading: vaccination on a 2d grid
void vaccinate(RNG &rng, int Nx, int Ny, int **S, int **R, double frac){
  int i,j,t,vaxed,r2;
  int Stot=msum(S,Nx, Ny);
  int tovax=(int)(frac*((double)Stot));
  std::uniform_real_distribution<> runif(0.0, 1.0);
  int Scum[Nx*Ny+1];
  int flag=0;

  Scum[0]=0;t=1;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      Scum[t]=Scum[t-1]+S[i][j];
      t++;
    }
  }
  //fprintf(stderr, "Stot=%d\n", Stot);

  vaxed=0;
  while(vaxed<tovax){
    t=0;flag=0;
    r2 = (int)(runif(rng)*Stot);
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	if(r2>Scum[t]&&r2<=Scum[t+1] &&S[i][j]>0){
	  S[i][j]--;R[i][j]++;vaxed++;
	  flag=1;break;
	}
	t++;
      }
      if(flag)
	break;
    }
  }
  return;
}

//Choose a grid type and set the coupling matrices; 
//1:complete; 2:gamma, 3:ring small_world, 4:2d grid small world,
//5:clique preferential attachment, 6:random graph G(n,p),
//7: ring preferential attachment, 8: star shaped
//output network data to fd
//param1: shape for gamma distribution/ mean degree for random network
//For meaning of constout, see above
void setgrid(RNG &rng, int gridtype, double **epsM, int **epsI, int ncomp, int Nx, int Ny, int Nc, int ***sparse, double RepNo, double leak, double eps, double param1, double rewiring_prob, int param2, int clq_sz, int scalefree_clq, double *mean_eps, double *epsvar, int constout, FILE *fd){
  int i,j,flag,tot;
  double t,meaneps,RepNoA;

  if(constout<0 || constout>3){
    fprintf(stderr, "ERROR in \"setgrid\": constout must be 0, 1, 2 or 3. Exiting.\n");exit(0);
  }
  
  inittozero(epsM,ncomp,ncomp);
  if(gridtype==1){//complete
    sym_coupling(eps, epsM, ncomp);
    inittoconst(epsI, ncomp, ncomp, 1);
    fprintf(fd, "complete, symmetric network with %d compartments, each of population %d.\n", ncomp, Nc);
  }
  else if(gridtype==2){//gamma distributed
    gamma_coupling(rng, leak, epsM, ncomp, Nc, param1, param2, constout);
    inittoconst(epsI, ncomp, ncomp, 1);
    meaneps=leak/((double)Nc)/((double)(ncomp-1));
    fprintf(fd, "complete network with %d compartments, each of population %d, and coupling-strengths set from a Gamma distribution with mean = %.4e, SD = %.4e and shape = %.4f.\n", ncomp, Nc, meaneps, sqrt(meaneps*meaneps/param1), param1);
  }
  else if(gridtype==3){//small world
    flag=0;tot=0;
    while(!flag && tot<10){//get a strongly connected network
      flag=small_world(rng, epsI, ncomp, clq_sz, rewiring_prob);
      tot++;
    }
    if(tot==10){
      fprintf(stderr, "Could not get a strongly connected (small world) network. Exiting.\n");exit(0);
    }
    //now set up sparse and epsM
    sparse_from_epsI(epsI, ncomp, sparse);
    coupling_from_epsI(eps, epsI, epsM, ncomp, constout);
  
    fprintf(fd, "Regular ring lattice small world network with %d compartments, each of population %d, original neighbourhood size %d, and rewiring probability %.4f.\n", ncomp, Nc, clq_sz, rewiring_prob);
    //fprintf(stderr, "Mean path=%.4f\n", mean_path_length(epsI,ncomp));
    //printsparse(*sparse,ncomp);
  }
  else if(gridtype==4){//2d grid small world
    flag=0;tot=0;
    while(!flag && tot<10){//ensure strongly connected network
      flag=nn_smallworld(rng, epsI, Nx, Ny, rewiring_prob);
      tot++;
    }
    if(tot==10){
      fprintf(stderr, "Could not get a strongly connected (2d grid small world) network. Exiting.\n");exit(0);
    }
    //set epsM and sparse
    sparse_from_epsI(epsI, ncomp, sparse);
    coupling_from_epsI(eps, epsI, epsM, ncomp, constout);
    fprintf(fd, "2D grid small world network with %d X %d compartments, each of population %d, and with rewiring probability %.4f.\n", Nx, Ny, Nc, rewiring_prob);

    /* fprintf(stderr, "Mean path=%.4f\n", mean_path_length(epsI,ncomp)); */
    /* printsparse(*sparse,ncomp); */
    /* exit(0); */
  }
  else if(gridtype==5){//scale free network with parameter scalefree_clq
    scale_free(rng, epsI, ncomp, scalefree_clq);
    //now set up sparse and epsM
    sparse_from_epsI(epsI, ncomp, sparse);
    coupling_from_epsI(eps, epsI, epsM, ncomp, constout);
    fprintf(fd, "Preferential attachment network with %d compartments, each of population %d, initial clique of size %d.\n", ncomp, Nc, scalefree_clq);
  }
  else if(gridtype==6){//random network
    flag=0;tot=0;
    while(!flag && tot<100){//ensure strongly connected network
      flag=randomgraph(rng, epsI, ncomp, param1);
      tot++;
    }
    if(tot==100){
      fprintf(stderr, "Could not get a strongly connected random network. Exiting.\n");exit(0);
    }
    sparse_from_epsI(epsI, ncomp, sparse);
    coupling_from_epsI(eps, epsI, epsM, ncomp, constout);
    fprintf(fd, "random graph with %d compartments, each of population %d, and edge probability %.4f.\n", ncomp, Nc, (double)param1/(double)(ncomp-1));
  }
  else if(gridtype==7){
    ring_scalefree(rng, epsI, ncomp, param2, clq_sz, scalefree_clq);
    //now set up sparse and epsM
    sparse_from_epsI(epsI, ncomp, sparse);
    coupling_from_epsI(eps, epsI, epsM, ncomp, constout);
    fprintf(fd, "Ring-based preferential attachment network with %d compartments, each of population %d. Initial ring of size %d, initial cliques of size %d, and %d edges attach at each subsequent stage.\n", ncomp, Nc, param2, clq_sz, scalefree_clq);
  }
  else if(gridtype==8){
    star_coupling(epsI, ncomp);
    sparse_from_epsI(epsI, ncomp, sparse);
    coupling_from_epsI(eps, epsI, epsM, ncomp, constout);
    fprintf(fd, "Star shaped network with %d compartments, each of population %d.\n", ncomp, Nc);
  }
  (*mean_eps)=0.0;
  for(i=0;i<ncomp;i++){
    for(j=0;j<ncomp;j++){
      if(j != i)
	(*mean_eps)+=epsM[i][j];
    }
  }
  (*mean_eps)/=(double)ncomp;
  (*epsvar)=0.0;
  for(i=0;i<ncomp;i++){
    for(j=0;j<ncomp;j++){
      if(j != i)
	(*epsvar)+=pow(epsM[i][j]-(*mean_eps)/(double)(ncomp-1),2.0);
    }
  }
  (*epsvar)/=(double)(ncomp*(ncomp-1));

  RepNoA=RepNo*msum(epsM,ncomp,ncomp)/(double)ncomp;
  fprintf(fd, "\nset eps = %.4f, mean eps = %.4f, set \"R0\" = %.4f, actual \"R0\" = %.4f, eps_ij mean (SD) = %.4e (%.4e), mean path=%.4f, degree SD=%.4f, constout=%d\n", eps, (*mean_eps), RepNo, RepNoA, (*mean_eps)/(double)(ncomp-1), sqrt((*epsvar)), mean_path_length(epsI,ncomp), sqrt(degree_var(*sparse,ncomp)),constout);
  

  if(*sparse){
    fprintf(fd, "Network structure:\n");
    printsparse(fd, *sparse, ncomp);
  }
  if((gridtype==2 || gridtype==3 || gridtype==4 || gridtype==5 || gridtype==6 || gridtype==7)){//variable out infection
    fprintf(fd, "\nin/out/tot-infection values:\n");
    for(i=0;i<ncomp;i++){
      t=vsum(epsM[i],ncomp);
      fprintf(fd, "%d: %.4f, %.4f, %.4f (%d)\n", i, epsM[i][i], t-epsM[i][i], t, (*sparse)?(*sparse)[i][0]-1:ncomp-1);
    }
    //fprintf(fd, "\nFull epsM:\n");
    //printmat(fd, epsM, ncomp, ncomp);
  }
  
  return;
}





//Dijkstra's algorithm for finding path lengths from initial node k
//Based on pseudocode at en.wikipedia.org/wiki/Dijkstra%27s_algorithm
//adj is the adjacency matrix, Q is vertices not yet dealt with
//dists stores the minimum distances
void Dijkstra(int k, int **adj, int n, int *dists){
  int j,pos;
  int Q[n],tot=n;
  for(j=0;j<n;j++){
    dists[j]=n;Q[j]=1;
  }
  dists[k]=0;
  while(tot){
    pos=vminpos(dists,Q,n);
    Q[pos]=0;tot--;
    for(j=0;j<n;j++){//neighbours
      if(adj[pos][j] && Q[j]){
	if(dists[pos]+adj[pos][j]<dists[j])
	  dists[j]=dists[pos]+adj[pos][j];
      }
    }
  }
  return;
}

//Mean path length between vertices of a (possibly integer weighted) digraph
double mean_path_length(int **adj, int n){
  int i,j;
  double meandist=0.0;
  int **dists=imatrix(n,n);
  for(i=0;i<n;i++)
    Dijkstra(i, adj, n, dists[i]);

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      meandist+=dists[i][j];
    }
  }

  free_imatrix(dists,n,n);
  return meandist/((double)(n*(n-1)));
}

//Variance in the degrees for a non-regular graph
double degree_var(int **sparse, int n){
  int i;
  double meandeg=0;
  double degvar=0;
  if(!sparse)//regular graph
    return 0.0;

  for(i=0;i<n;i++)
    meandeg+=sparse[i][0]-1;

  meandeg/=(double)n;
  for(i=0;i<n;i++)
    degvar+=pow((double)(sparse[i][0]-1)-meandeg,2.0);

  return degvar/=(double)n;
}

// is the integer "i" in the list lst?
int isinlist(int i, int ilst[], int tot){
  int j;
  for(j=0;j<tot;j++){
    if(i==ilst[j])
      return 1;
  }
  return 0;
}



void strongconnect_light(int i, int **mat, int n, int *inds, int *lowlink, int *ind, int *totS, int *Sc, int *Scpos){
  int j,k;

  /* fprintf(stderr, "current: "); */
  /* printvec1(Sc, *Scpos); */
  // Set the depth index for v to the smallest unused index
  inds[i]=*ind;lowlink[i]=*ind;(*ind)++;
  Sc[(*Scpos)++]=i;//add to stack

  // Consider successors of v
  for (j=0;j<n;j++){
    if(mat[i][j]){// each out-edge
      if(inds[j]<0){
	strongconnect_light(j, mat, n, inds, lowlink, ind, totS, Sc, Scpos);
	lowlink[i]=min(lowlink[i],lowlink[j]);
      }
      else if(isinlist(j, Sc, *Scpos))
	lowlink[i]=min(lowlink[i],inds[j]);
    }
  }
  // If node i is a root node, pop the stack
  if (lowlink[i]==inds[i]){
    k=(*Scpos)-1;
    while(Sc[k]!=i)
      k--;
    (*totS)++;(*Scpos)=k;//reset stack
  }

}

int Tarjan_light(int **mat, int n){
  int i, ind=0; //ind increments as we encounter vertices
  int *inds=(int *)malloc((size_t) ((n)*sizeof(int)));
  int *lowlink=(int *)malloc((size_t) ((n)*sizeof(int)));
  int *Sc= (int *)malloc((size_t) ((n)*sizeof(int)));
  int Scpos=0;//stack-counter
  int totS=0;

  for(i=0;i<n;i++)
    inds[i]=-1;

  for(i=0;i<n;i++){
    if(inds[i]<0)
      strongconnect_light(i,mat,n,inds,lowlink,&ind, &totS, Sc, &Scpos);
  }

  free((char *)Sc);
  free((char *)inds);
  free((char *)lowlink);

  return totS;
}

//Is the digraph with adjacency matrix AM strongly connected?
int is_SC(int **AM, int n){
  int comps=Tarjan_light(AM, n);
  if(comps==1)//1 SCC
    return 1;
  return 0;
}

//The contagion probabilities (complete, symmetric graph):
//p[j][i] is the probability that given j compartments in total
//and given an outbreak in 1, the outbreak will eventually
//hit a total of i compartments.
//"a" is the probility of *not* having contagion to a neighbour
//"prior" is prior immunity as a fraction
double **alphapoly(double eps, double RepNo, int Nc, int maxcomp, double prior){
  double **p=dmatrix(maxcomp+1,maxcomp+1);
  double a;
  int i,ncomp;
  double det_sz;
  double s0=1.0-prior;//initial susceptible fraction

  inittozero(p,maxcomp,maxcomp);

  //absolute epi size in each compartment in the deterministic case,
  //via Lambert W function
  det_sz=(double)Nc*(s0+boost::math::lambert_w0(-s0*RepNo*exp(-s0*RepNo))/RepNo);

  ncomp=2;
  a=pow((1.0-eps*(s0*RepNo-1.0)/((double)(maxcomp-1))),det_sz);
  //fprintf(stderr,"eps=%.4f, det_sz=%.4f, a=%.4f\n", eps, det_sz, a);//exit(0);

 
  p[ncomp][1]=a;
  p[ncomp][2]=1.0-a;

  for(ncomp=3;ncomp<=maxcomp;ncomp++){
    for(i=1;i<ncomp;i++)
      p[ncomp][i]=((double)(ncomp-1))/((double)(ncomp-i))*pow(a,(double)i)*p[ncomp-1][i];
    p[ncomp][ncomp]=1.0-vsum(p[ncomp],ncomp);
  }

  return p;
}

//Here "eps" is Theta/Nc, namely (n-1)*eps in the write-up
double **alphapolyA(double eps, double RepNo, int Nc, int maxcomp, double prior){
  double **p=dmatrix(maxcomp+1,maxcomp+1);
  double a;
  int m,ncomp;
  double det_sz;
  double s0=1.0-prior;//initial susceptible fraction

  inittozero(p,maxcomp+1,maxcomp+1);

  //absolute epi size in each compartment in the deterministic case,
  //via Lambert W function
  //det_sz=(double)Nc*(s0+boost::math::lambert_w0(-(1.0-eps)*s0*RepNo*exp(-(1.0-eps)*s0*RepNo))/((1.0-eps)*RepNo));
  //uncorrect R0
  det_sz=(double)Nc*(s0+boost::math::lambert_w0(-s0*RepNo*exp(-s0*RepNo))/(RepNo));

  if(maxcomp==1)
    a=1.0;
  else
    a=pow((1.0-eps*(s0*RepNo-1.0/(1.0-eps))/(double)(maxcomp-1)),det_sz);
  //fprintf(stderr,"eps=%.4f, det_sz=%.4f, a=%.4f\n", eps, det_sz, a);//exit(0);

  p[1][1]=1.0;
  for(ncomp=2;ncomp<=maxcomp;ncomp++){
    for(m=1;m<ncomp;m++)
      p[ncomp][m]=((double)(ncomp-1))/((double)(ncomp-m))*pow(a,(double)m)*p[ncomp-1][m];
    p[ncomp][ncomp]=1.0-vsum(p[ncomp],ncomp);
  }

  return p;
}

//Theoretical mean outbreak size in the complete, symmetrical network,
//as a fraction of the theoretical size in the strongly coupled case
double mean_outbreak_sym(double leak, double RepNo, int Nc, int ncomp, double prior){
  double eps=leak/((double)Nc);
  double **p=alphapolyA(eps, RepNo, Nc, ncomp, prior);
  int m;
  double sz=0;

  for(m=1;m<=ncomp;m++)
    sz+=(double)m*p[ncomp][m];

  sz/=(double)ncomp;
  free_dmatrix(p,ncomp+1,ncomp+1);
  return sz;
}

//
// get the number of non-empty lines and 
// maximum linelength in a file
//
unsigned long numflines(const char *fname, int *maxline){
  FILE *fd;
  int i=0, c=0;
  unsigned long numl=0;
  (*maxline)=0;

  if(!(fd= fopen(fname, "r"))){
    fprintf(stderr, "WARNING in numflines: \"%s\" could not be opened for reading.\n", fname);
    exit(0);
  }

  while((c=getc(fd)) != EOF){
    i++;
    if (c=='\n'){
      numl++;
      (*maxline)=max(i, (*maxline));
      i=0;
    }
  }
  if(i!=0){ // final line does not end in newline
    numl++;
    (*maxline)=max(i, (*maxline));
  }
  (*maxline)+=2; // to hold the terminal character plus possible EOF

  fclose(fd);
  return numl;
}

int gtline(FILE *fp, char s[], int lim){
  /* store a line as a string, including the terminal newline character */
  int c=0, i;

  for(i=0; i<lim-1 && (c=getc(fp))!=EOF && c!='\n';++i){s[i] = c;}
  if (c == '\n') {s[i++] = c;}
  s[i] = '\0';
  if(i==lim-1 && c!='\n' && c!=0)
    fprintf(stderr, "WARNING in gtline: line length exceeded. *%s*\n", s);
  return i;
}

// chops a string (allocates memory for the output)
char *strchop2(char *t, int n, int n1){
  char *s;
  int i=0;
  if((n<= (int)strlen(t)) && n1>0){
    s = (char*) malloc(sizeof(char) * (n1+1));
    while (t[n+i] && (i< n1)){s[i] = t[n+i];i++;}
    s[i] = '\0';
  }
  else
    s = strdup("");
  
  return s;
}
 
 //nth "word" from a string (any characters allowed)
char *getnthwd(char *s, int n){
  int i, j, k=0;
  char *v=NULL;
  while(s[k]){
    while(s[k] && isspace(s[k])){k++;} // skip space
    for(i=0;i<n-1;i++){
      while(s[k] && !isspace(s[k])){k++;} // skip word
      while(s[k] && isspace(s[k])){k++;} // skip nonword
    }
    j=0;
    while(!isspace(s[k])){j++;k++;} // get the word
    v = strchop2(s, k-j, j);
    return v;
  }
  if(!v)
    v=strdup("");
  return v;
}


//assume first column is time, and resample to get the state at intervals
//of "sampling". Return number of lines
int resample(const char infile[], const char outfile[], double sampling){
  FILE *fdin, *fdout;
  unsigned long numl=0;
  int maxl=0, linelen;
  int index, indexnew;
  double t;
  char *wdnew, *wdold;
  char *lineold, *linenew;
  int tot=0;
  numl = numflines(infile, &maxl); // number of lines and max line length
  if(numl==0)
    return 0;
  
  lineold = (char*) malloc(sizeof(char) * (maxl));
  linenew = (char*) malloc(sizeof(char) * (maxl));

  fdin=fopen(infile, "r");

  if(!(fdout=fopen(outfile, "w"))){
    fprintf(stderr, "ERROR in resample - File: %s could not be opened...\n", outfile);
    exit(0);
  }

  gtline(fdin, lineold, maxl);
  wdold=getnthwd(lineold, 1);
  index=(int)((atof(wdold))/sampling);
  fprintf(fdout, "%.4f%s", atof(wdold),lineold+strlen(wdold));
  tot++;
  while((linelen = gtline(fdin, linenew, maxl)) > 0){
    wdnew=getnthwd(linenew, 1);
    t=atof(wdnew);indexnew=(int)(t/sampling);
    while(index<indexnew){//gone past, print last
      index++;
      fprintf(fdout, "%.4f%s", (double)(index)*sampling,lineold+strlen(wdold)); 
      tot++;
    }
    strcpy(lineold, linenew);
    free(wdold);wdold=strdup(wdnew);free(wdnew);
    
  }


  free(wdold);
  free(lineold);free(linenew);
  fclose(fdin);fclose(fdout);
  return tot;
  
}
