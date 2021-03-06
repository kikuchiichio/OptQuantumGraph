/************************************/
/*PARTICLE SWARM OPTIMIZATION */
/************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <iostream>
/*CONSTANTS */
#define NOPS 10    /*NUMBER OF PARTICLES*/
#define LIMITL 256 /*MAXIMUM SIZE OF THE ARRAY FOR PARTICLE DATA*/
#define ILIMIT 150 /*MAXIMUM NUMBER OF ITERATION*/
#define SEED 32767 /*INITIALIZATION OF RANDOM NUMBER*/
#define W 0.7      /*INERTIA CONSTANT*/
#define C1 1.4     /*CONSTANT OF ATTRACTION FROM PERSONAL BEST*/
#define C2 1.4     /*CONSTANT OF ATTRACTION FROM GROUP BEST*/

#define DMAX 1 /*STRUCTURE FOR ARBITRARY COORDINATE*/
struct point
{
  double x[DMAX];
};

typedef std::complex<double> Complex;

extern "C"
{
  void zheev_(const char &JOBZ, const char &UPLO,
              const int &N, Complex *A, const int &LDA,
              double *WW, Complex *WORK, const int &LWORK,
              double *RWORK,
              int &INFO, int JOBZlen, int UPLOlen);
  void zgetrf_(int *M, int *N, Complex *A, int *LDA, int *IPIV, int *INFO);
  void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
  void zhegv_(const int &ITYPE, const char &JOBZ, const char &UPLO,
              const int &N,
              Complex **A, const int &LDA,
              Complex **B, const int &LDB,
              double *WW, Complex *WORK, const int &LWORK,
              double *RWORK,
              int &INFO, int JOBZlen, int UPLOlen);
};

using namespace std;

/*STRUCTURE FOR PARTILE*/
struct particle
{
  struct point pos;     /*POSITION*/
  double value;         /*EVALUATED VALUE*/
  struct point v;       /*VELOCITY*/
  struct point bestpos; /*THE BEST OF EVER-REACHED POSITIONS*/
  double bestval;       /*THE BEST OF EVER-EVALUATED VALUES*/
};

/*PROTOTYPES OF FUNCTIONS */
void initps(struct particle ps[]);   /*INITIALIZATION OF PARTICLE SWARM*/
void printps(struct particle ps[]);  /*PRINT OUT OF PARTICLE SWARM DATA*/
double frand();                      /*UNIFORM REAL RANDOM NUMBER FROM ZERO TO ONE*/
double calcval(double *x);           /*THE OBJECTIVE FUNCTION*/
double calcresidue(double *x);       /*THE RESIDUE OF THE OPTIMIZATION*/
void optimize(struct particle ps[]); /*OPTIMIZATION*/
void setgbest(struct particle ps[]); /*FIND THE GROUP BEST*/

/*GLOBAL VARIABLES*/
struct point gbestpos; /*THE BEST POSITION IN THE GROUP*/
double gbestval;       /*THE BEST VALUE IN THE GROUP*/

/**********************************************/
/*UNIFORM REAL RANDOM NUMBER FROM ZERO TO ONE */
/**********************************************/
double doublerand(void)
{
  double result;
  while ((result = (double)rand() / RAND_MAX) >= 1)
    ;
  return result;
}

/*VARIABLES TO CONSTRUCT A QUANTUM GRAPH*/
double lambda[10];
double L[10];
double A[10];
int B = 5;
double k;

double h00(double k)
{
  double cc = 0.;
  cc += lambda[0] / k;
  for (int i = 1; i <= B; i++)
  {
    cc += 1. / tan(k * L[i]);
  }
  return cc;
}

double hjj(double k, int j)
{
  double cc = 0.0;
  cc += lambda[j] / k;
  cc += 1. / tan(k * L[j]);
  return cc;
}

complex<double> h0j(double k, int j)
{
  complex<double>
      cc;
  cc = complex<double>(cos(A[j] * L[j]), -sin(A[j] * L[j]));
  cc /= sin(k * L[j]);
  cc *= -1;
  return cc;
}

double zetah(double k)
{
  complex<double> A[(B + 1) * (B + 1)];
  A[0] = h00(k);
  for (int i = 1; i <= B; i++)
  {
    A[0 + (B + 1) * i] = h0j(k, i);
    A[i + (B + 1) * 0] = conj(A[0 + (B + 1) * i]);
    A[i + (B + 1) * i] = hjj(k, i);
  }

  int m = B + 1;   // ???????????????
  int n = B + 1;   // ???????????????
  int lda = B + 1; // m????????????

  int info;                    // ????????????????????????0?????????
  int ipiv[B + 1];             // ????????????m,n??????????????????????????????
  int lwork = B + 1;           // n????????????
  complex<double> work[B + 1]; // ????????????lwork????????????

  zgetrf_(&m, &n, A, &lda, ipiv, &info);

  double det = 1;
  for (int i = 0; i <= B; i++)
  {
    int j = i + (B + 1) * i;
    det *= real(A[j]);
    if (ipiv[i] != i + 1)
      det = -det;
  }
  return det;
}

double zetak(double k)
{
  double cc = h00(k);

  for (int i = 1; i <= B; i++)
  {
    double absh0j = abs(h0j(k, i));
    cc -= absh0j * absh0j / hjj(k, i);
  }
  for (int i = 1; i <= B; i++)
  {
    cc *= hjj(k, i);
  }
  return cc;
}

int main()
{


  for (int i = 0; i <= 10; i++)
  {
    L[i] = 1;
    A[i] = 0;
    lambda[i] = 1.;
  }

  for (double k = 0.05; k <= 3; k += 0.01)
  {
    cout << "K " << k << " " << zetak(k) << " :" << zetah(k) << endl;
  }

  struct particle ps[LIMITL]; /*THE SWARM OF PARTICLES*/
  int i;                      /*ITERATOR OF THE REPETITION*/

  /*INITIALIZE THE RANDOM NUMBER*/
  srand(SEED);

  /*INITIALIZATION OF PARTICLE SWARM*/
  initps(ps);

  printps(ps);

  /*OPTIMIZATION*/
  for (i = 0; i < ILIMIT; ++i)
  {
    optimize(ps);
    printf("\n No. %d\n", i);
    printps(ps);
  }

  /*COMPUTE RESIDUE ( THE REMMNANT OF GRAD)*/
  //calcresidue(gbestpos.x);

  return 0;
}

/*************************************/
/*PARTICLE SWARM OPTIMIZATION        */
/*************************************/

void optimize(struct particle ps[])
{
  int i, j;
  double r1, r2; /*RANDOM NUMBER*/

  for (i = 0; i < NOPS; ++i)
  {

    /*SET RANDOM NUMBERS*/
    r1 = frand();
    r2 = frand();

    /*VELOCITY UPDATE*/
    for (j = 0; j < DMAX; j++)
    {
      ps[i].v.x[j] = W * ps[i].v.x[j] + C1 * r1 * (ps[i].bestpos.x[j] - ps[i].pos.x[j]) + C2 * r2 * (gbestpos.x[j] - ps[i].pos.x[j]);
    }

    /*POSITION UPDATE*/
    for (j = 0; j < DMAX; j++)
    {
      ps[i].pos.x[j] += ps[i].v.x[j];
    }

    /*UPDATE THE BEST FOR EACH PARTICLE*/
    ps[i].value = calcval(ps[i].pos.x);
    //printf("%lf %lf ",ps[i].pos.x[0],ps[i].value);

    if (ps[i].value < ps[i].bestval)
    {
      ps[i].bestval = ps[i].value;
      ps[i].bestpos = ps[i].pos;
    }
  }

  /*UPDATE THE BEST IN THE GROUP*/
  setgbest(ps);
}

double calceigen(double *XIN)
{
  complex<double> U[(B + 1) * (B + 1)];
  int N = B + 1;
  double E[B + 1];
  int i, j, info;
  complex<double> cwork[4 * (B + 1)];
  double rwork[4 * (B + 1)];
  U[0] = complex<double>(h00(XIN[0]), 0.0);
  for (i = 1; i <= B; i++)
  {

    U[0 + (B + 1) * i] = h0j(XIN[0], i);
    U[i + (B + 1) * 0] = conj(h0j(XIN[0], i));
    U[i + (B + 1) * i] = Complex(hjj(XIN[0], i));
  }
 
  zheev_('V', 'U', N, (Complex *)U, N, E, cwork, 4 * N, rwork, info, 1, 1);
  double det = 1;
  for (i = 0; i <= B; i++)
  {
    det *= E[i];
  }

  return 0;
}

double calcval(double *XIN)
{
  double zeta_ = zetak(XIN[0]);
  calceigen(XIN);
  return zeta_ * zeta_;
}


/****************************/
/*FIND THE BEST IN THE GROUP*/
/****************************/
void setgbest(struct particle ps[])
{
  int i, j;
  double besti;
  double x[DMAX];
  besti = ps[0].value;
  for (j = 0; j < DMAX; j++)
  {
    x[j] = ps[0].pos.x[j];
  }

  /*FIND THE CURRENT BEST*/
  for (i = 0; i < NOPS; ++i)
  {
    //printf("---%d %lf\n",i,ps[i].value);
    if (ps[i].value < besti)
    {
      besti = ps[i].value;
      for (j = 0; j < DMAX; j++)
      {
        x[j] = ps[i].pos.x[j];
      }
      //printf("%d besti=%lf\n ",i,besti);
    }
    /*IF POSSIBLE, UPDATE THE BEST*/
    if (besti < gbestval)
    {
      gbestval = besti;
      for (j = 0; j < DMAX; j++)
      {
        gbestpos.x[j] = x[j];
      }
    }
  }
//  printf("gbestval=%lf\n", gbestval);
}

/**********************************/
/*INITIALIZATION OF PARTICLE SWARM*/
/**********************************/
void initps(struct particle ps[])
{

  int i, j;
  double x[DMAX];

  for (i = 0; i < NOPS; ++i)
  {
    /*POSITION*/
    for (j = 0; j < DMAX; j++)
    {
      x[j] = ps[i].pos.x[j] = 0.01 * (frand());
    }
    /*EVALUATED VALUE*/
    ps[i].value = calcval(x);
    /*VELOCITY*/
    for (j = 0; j < DMAX; j++)
    {
      ps[i].v.x[j] = 0.000001 * (frand() * 2 - 1.0);
    }
    /*THE BEST POSITION*/
    for (j = 0; j < DMAX; j++)
    {
      ps[i].bestpos.x[j] = ps[i].pos.x[j];
    }
    /*THE BEST EVALUATED VALUE*/
    ps[i].bestval = ps[i].value;
    printf("%d %lf %lf %lf\n", i, ps[i].pos.x[0], ps[i].v.x[0], ps[i].value);
  }
  /*THE BEST IN THE GROUP*/
  gbestval = ps[0].value;
  printf("gbestval=%lf", gbestval);
  for (j = 0; j < DMAX; j++)
  {
    gbestpos.x[j] = ps[0].pos.x[j];
  }
  setgbest(ps);
}

/***********************************/
/*PRINT OUT OF PARTICLE SWARM DATA */
/***********************************/
void printps(struct particle ps[])
{
  double average[DMAX], sqaverage[DMAX];
  int debug = 1;
  for (int j = 0; j < DMAX; j++)
  {
    average[j] = 0;
    sqaverage[j] = 0;
  }
  if (debug == 1)
  {
    for (int i = 0; i < NOPS; ++i)
    {
      printf("%d ", i);
      for (int j = 0; j < DMAX; j++)
      {
        average[j] += ps[i].pos.x[j];
        sqaverage[j] += ps[i].pos.x[j] * ps[i].pos.x[j];
        printf("%lf ", ps[i].pos.x[j]);
      }
      for (int j = 0; j < DMAX; j++)
      {
        printf("%lf ", ps[i].v.x[j]);
      }
      printf("%lf ", ps[i].value);
      printf("\n");
    }
    printf("\nSTATISTICS: AVR= ");
    for (int j = 0; j < DMAX; j++)
    {
      printf(" %lf ", average[j] / NOPS);
    }
    printf(" VAR=");
    for (int j = 0; j < DMAX; j++)
    {
      printf(" %lf ", sqaverage[j] / NOPS - average[j] * average[j] / NOPS / NOPS);
    }
  }
  printf("  BEST: ");
  for (int j = 0; j < DMAX; j++)
  {
    printf("%lf ", gbestpos.x[j]);
  }

  printf(" -> %lf\n", gbestval);
  //printf("gbestval=%lf %lf\n", gbestval, calcval(gbestpos.x));

  double xxx[1];

  for (double rx = 0.5; rx <= 1.5; rx += 0.1)
  {
    xxx[0] = rx * gbestpos.x[0];
  }
}

void printps___(struct particle ps[])
{
  int debug = 1;
  if (debug == 1)
  {
    for (int i = 0; i < NOPS; ++i)
    {
      printf("%d ", i);
      for (int j = 0; j < DMAX; j++)
      {
        printf("%lf ", ps[i].pos.x[j]);
      }
      for (int j = 0; j < DMAX; j++)
      {
        printf("%lf ", ps[i].v.x[j]);
      }
      printf("%lf ", ps[i].value);

      printf("\n");
    }
  }
  printf("\n BEST: ");
  printf(" -> %lf\n", gbestval);
}

/*********************************************/
/*UNIFORM REAL RANDOM NUMBER FROM ZERO TO ONE*/
/*********************************************/
double frand(void)
{
  double result;
  while ((result = (double)rand() / RAND_MAX) >= 1)
    ;
  return result;
}
