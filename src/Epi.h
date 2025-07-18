#include <stdio.h>
#include <math.h>
#include <stdlib.h> /* for srand() */
#include <string.h>
#include <time.h> /* for random seeding */
#include <ctype.h>
#include <unistd.h>
#include <glpk.h>
#include <dirent.h>
#include <sys/stat.h>
#include <random>
#include <boost/math/special_functions/lambert_w.hpp>
//#define _XOPEN_SOURCE
#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

typedef std::mt19937 RNG;
#define NR_END 1
#define FREE_ARG char*

double **dmatrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrh, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrh, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrh, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrh, long nch);
void free_imat(int **A, int numA);
void printmat(FILE *fd, double **imat, int n, int m);
void printmat(double **imat, int n, int m);
void printmat(int **imat, int n, int m);
void printsparse(int **sparse, int n);
void printsparse(FILE *fd, int **sparse, int n);
void printvec(int *ivec, int n);
void printvec(double *ivec, int n);
void fprintvec(FILE *fd, int *ivec, int n);
double **cpmat(double **mat, int n, int m);
void cpmat(double **mat, double **out, int n, int m);
int **cpmat(int **mat, int n, int m);
void cpmat(int **mat, int **out, int n, int m);
void cpvec(double *vecin, double *vecout, int n);
void cpvec(int *vecin, int *vecout, int n);
double rowsum(double **mat, int n, int m, int rowk);
double vsum(double *v, int n);
int vsum(int *v, int n);
double msum(double **v, int n, int m);
int msum(int **v, int n, int m);
double vmax(double v[], int n);
double vmin(double v[], int n);
int d2tod1(int x, int y, int Nx, int Ny);
void d1tod2(int t, int Nx, int Ny, int *x, int *y);
int neighbour_r(int i, int j, int Nx, int Ny);
int neighbour_l(int i, int j, int Nx, int Ny);
int neighbour_u(int i, int j, int Nx, int Ny);
int neighbour_d(int i, int j, int Nx, int Ny);
int neighbour_r(int t, int Nx, int Ny);
int neighbour_l(int t, int Nx, int Ny);
int neighbour_u(int t, int Nx, int Ny);
int neighbour_d(int t, int Nx, int Ny);
void inittozero(int *vec, int len);
void inittozero(double *vec, int len);
void inittozero(int **mat, int n, int m);
void inittozero(double **mat, int n, int m);
void inittoconst(int *vec, int n, int c);
void inittoconst(double *vec, int n, double c);
void inittoconst(int **mat, int n, int m, int c);
void inittoconst(double **mat, int n, int m, double c);
double Rt_clas(int *S, double *V, int ncomp, int *comp_pop, double **epsM, double ki, double kr, int totpop);
double Rt_clas(int *S, double *V, int ncomp, double **epsM, double ki, double kr);
double Rt_clas(int Stot, double ki, double kr, int Ntot);
double Rt_clas_nn(int **S, double **V, double **epsM, double ki, double kr, int Nx, int Ny);
double SpecRad(double **M, int n, double tol);
void qvec(double ki, double kr, int *S, double *V, double **eps, int n, double *s, double tol, int *iter);
void qvec(double ki, double kr, int **S, double **V, double **eps, int Nx, int Ny, double *s, double tol, int *iter);
double outbreak_prob(double ki, double kr, int *S, double *V, double **eps, int ncomp, double *s, double tol, int *iter);
double outbreak_prob(double ki, double kr, int **S, double **V, double **eps, int Nx, int Ny, double *s, double tol, int *iter);
double inf_prop(int *S, int *I, double *V, double **eps, double ki, int i, int j);
double inf_prop(int **S, int **I, double **V, double **eps, double ki, int i, int j, int Nx, int Ny);
double inf_prop_tot(int *S, int *I, double *V, double **eps, double ki, int j, int ncomp);
double intensity(int *S, int *I, double *V, double **epsM, double ki, double kr, double *vr, double *vitot, int ncomp, int **sparse);
double intensity(int **S, int **I, double **V, double **epsM, double ki, double kr, double **vr, double **vitot, int Nx, int Ny);
int episz(int *R, int *R0, int ncomp, int comp_pop, double *mean_outbreak_sz, double *mean_outbreak_sz_cond, int *totepis);
int episz(int **R, int **R0, int Nx, int Ny, int comp_pop, double **mean_outbreak_sz, double **mean_outbreak_sz_cond, int *totepis);
void NextG(double **NextGen, int *S, double *V, double **epsM, double ki, double kr, int ncomp);
void NextG(double **NextGen, int **S, double **V, double **epsM, double ki, double kr, int Nx, int Ny);
void printtSIR(FILE *fd, double t, int *S, int *I, int *R, int ncomp);
// A total of meantot time-points; dt is the dime increment
void printSIRmean(FILE *fd, int ncomp, int meantot, int totruns, double dt, int *S0, int *I0, int *R0, double **meanS, double **meanI, double **meanR);
void printscaledmat(FILE *fd, int **A, int N, int Nx, int Ny);
void printscaledmat(FILE *fd, int *A, int N, int Nx, int Ny);
void printSImatrices(const char basedir[], int **I, int **S, int Nx, int Ny, double Icomp, double Scomp, int tlast);
void printstarhist(int *Hist, int histnum, int histmax, int maxlen);
int SIRmeanupdate(int ncomp, int jstart, int jend, int *S, int *I, int *R, double **meanS, double **meanI, double **meanR);
void coupling_from_epsI(double eps, int **epsI, double **epsM, int ncomp, int constout);
void sparse_from_epsI(int **epsI, int ncomp, int ***sparse);
void sym_coupling(double eps, double **epsM, int ncomp);
void gamma_coupling(RNG &rng, double leak, double **epsM, int ncomp, int comp_pop, double shp, int sym, int constout);
int small_world(RNG &rng, double eps, double **epsM, int **epsI, int ncomp, int K, double beta, int ***sparse, int constout);
int nn_smallworld(RNG &rng, int **epsI, int Nx, int Ny, double rewiring_prob);
void scale_free(RNG &rng, double eps, double **epsM, int **epsI, int ncomp, int m0, int ***sparse, int constout);
double runSIRepi(RNG &rng, int ncomp, double *V, int *I0, int *S0, int *R0, int *I, int *S, int *R, double **meanI, double **meanS, double **meanR, double **epsM, double ki, double kr, double tmax, int printfull, int printmean, int meantot, double dt, FILE *fd, int **sparse);
double runSIRepi(RNG &rng, int Nx, int Ny, double **V, int **I0, int **S0, int **R0, int **I, int **S, int **R, int Iinit, double Icomp, double Scomp, double **epsM, double ki, double kr, double tmax, int sampling, int printall, FILE *fd4, FILE *fd1, int **sparse);

double mean_SIR_outbreak(RNG &rng, int ncomp, double *V, int *S0, int *R0, double **epsM, double ki, double kr, double tmax, int numsims, int choose_by_S, int **sparse);
double mean_SIR_outbreak(RNG &rng, int Nx, int Ny, double **V, int **S0, int **R0, double **epsM, double ki, double kr, double tmax, int numsims, int **sparse);
void setgrid(RNG &rng, int gridtype, double **epsM, int **epsI, int ncomp, int Nx, int Ny, int comp_pop, int ***sparse, double RepNo, double leak, double eps, double shp, double rewiring_prob, int init_ring, int clq_sz, int scalefree_clq, double *mean_eps, double *epsvar, int constout, FILE *fd);

void vaccinate(RNG &rng, int ncomp, int *S, int *R, double frac);
void vaccinate(RNG &rng, int Nx, int Ny, int **S, int **R, double frac);


int is_SC(int **AM, int n);
double mean_path_length(int **adj, int n);
double degree_var(int **sparse, int n);
double **alphapoly(double eps, double RepNo, int Nc, int maxcomp, double prior);double **alphapolyA(double eps, double RepNo, int Nc, int maxcomp, double prior);double mean_outbreak_sym(double leak, double RepNo, int Nc, int ncomp, double prior);
unsigned long numflines(const char *fname, int *maxline);
int gtline(FILE *fp, char s[], int lim);
char *strchop2(char *t, int n, int n1);
char *getnthwd(char *s, int n);
int resample(const char infile[], const char outfile[], double sampling);
