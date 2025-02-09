/*Part of structure.c.

  This bit of the program is involved in collecting data (DataCollection)
  and printing results (OutPutResults). */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structure.h"
#include "params.h"
#include "mymath.h"

void UpdateSums (double *Q, double *QSum, int *Z, double *P, double *PSum,
                 double *Fst, double *FstSum, int *NumAlleles,
                 int *AncestDist, double *Epsilon, double *SumEpsilon,
                 double *lambda, double *sumlambda,
                 double *LocPrior, double *sumLocPrior, int LocPriorLen);
double CalcLike (int *Geno, int *PreGeno, double *Q, double *P, int *Recessive,
                 double *sumindlike, double *indlike_norm);
double EstLogProb (double sumlikes, double sumsqlikes, int reps);
double KLDiv (int pop1, int pop2, double *P, double *LogP, int *NumAlleles, int reps);
void PrintKLD (FILE * file, double *P, double *LogP, int *NumAlleles, int reps, int format);
void PrintNET (FILE * file, double *P, int *NumAlleles, int reps, int format);

void PrintBanner (int rep, double *Alpha, double *Fst, double like, double *lambda);
void PrintMainParams (FILE * file, int rep, int argc, char *argv[]);
void PrintQ (FILE * file, int *Geno, int rep, double *QSum, struct IND *Individual,
             int *AncestDist, double *UsePopProbs,double *sumR);
void PrintP (FILE * file, int rep, int *Geno, double *PSum, int *Translation,
             int *NumAlleles, double *SumEpsilon, char *Markername);
void PrintMembership (FILE * file, double *QSum, struct IND *Individual);
void PrintSiteBySite (FILE * file, double *SiteBySiteSum, int rep,
                      char *Markername, double *PSum,
                      int *NumAlleles, int *Geno, int *Translation,double *SumEpsilon, struct IND *Individual);
void PrintSequences (FILE * file, int *Geno, char *Markername, double *SiteBySiteSum, int rep, int *Translation);
void PrintSums (FILE * file, int rep, double sumlikes,
                double sumsqlikes, double *FstSum, double *sumAlpha, double *sumlambda,
                double *sumR,double *varR, struct IND *Individual,
                double *sumLocPriors, int LocPriorLen, double DIC);
void PrintGeneName(FILE * file, int loc, char *Markername);
int EqualGeneNames(int loc1,int loc2,char *Markername);


/*=================================================*/
void
DataCollection (int *Geno, int *PreGeno,
                double *Q, double *QSum, int *Z, int *Z1,
                double *SiteBySiteSum, double *P, double *PSum,
                double *Fst, double *FstSum, int *NumAlleles,
                int *AncestDist, double *Alpha, double *sumAlpha,
                double *sumR, double *varR, double *like,
                double *sumlikes, double *sumsqlikes, double *R,
                double *Epsilon, double *SumEpsilon, double recomblikelihood,
                double *lambda, double *sumlambda, int *Recessive,
                double *LocPrior, double *sumLocPrior, int LocPriorLen,
                double *sumindlikes, double *indlikes_norm, int rep)
{
  int ind, pop;
  UpdateSums (Q, QSum, Z, P, PSum, Fst, FstSum, NumAlleles, AncestDist,
              Epsilon, SumEpsilon, lambda, sumlambda,
              LocPrior, sumLocPrior, LocPriorLen);
  if (LINKAGE) {
    for (ind = 0; ind < NUMINDS; ind++) {
      sumR[ind] += R[ind];
      varR[ind] += R[ind] * R[ind];
    }
  }

  if (LOCPRIOR && NOADMIX==0) {
    for (pop=0; pop<MAXPOPS; pop++)
      for (int loc=0; loc<=NUMLOCATIONS; loc++) {
        int pos = AlphaPos(loc, pop);
        sumAlpha[pos] += Alpha[pos];
      }
  } else if (!(NOADMIX) && (!(NOALPHA))) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      sumAlpha[pop] += Alpha[pop];
    }
  }


  if (COMPUTEPROB) {
    if (LINKAGE) {
      *like = recomblikelihood;
    }

    if (rep < BURNIN) {
      *like = CalcLike (Geno, PreGeno, Q, P, Recessive, NULL, NULL);
    } else {
      *like = CalcLike (Geno, PreGeno, Q, P, Recessive,
                        sumindlikes, indlikes_norm);
    }
    *sumlikes += *like;
    *sumsqlikes += (*like) * (*like);
  }
  /*printf("%f %f %f\n", *like, *sumlikes, *sumsqlikes); */
}

/*---------------------------------------------------*/
void
PrintLike (double like, int rep, int *Geno, int *PreGeno, double *Q,
           double *P,double recomblikelihood,
           int *Recessive)
{
  if (rep + 1 > BURNIN) {        /*already calculated */
    printf ("%6.0f\n", like);
  } else {
    if (LINKAGE) {
      printf("%6.0f\n",recomblikelihood);
    } else {
      printf ("%6.0f\n", CalcLike (Geno, PreGeno, Q, P, Recessive, NULL, NULL));
    }
  }
}


/*---------------------------------------------------*/
double
CalculateRAverage (double *R)
{
  int ind;
  double temp;
  temp = 0.0;
  for (ind = 0; ind < NUMINDS; ind++) {
    temp += R[ind];
  }
  return temp / (double) NUMINDS;
}

/*----------------------------------------------------*/
void
PrintUpdate (int rep, int *Geno, int *PreGeno, double *Alpha, double *Fst, double *P, double *Q,
             double like, double sumlikes, double sumsqlikes, int *NumAlleles,
             double *R, double *lambda, struct IND *Individual,
             double recomblikelihood, int *Recessive,
             double *LocPrior, int LocPriorLen)
{
  /*print a bunch of stuff to screen during run: rep, alpha, f, KLD, likelihood...
    Also occasionally print a header banner to define the variables. */

  double logprob=0.0;
  /*  int i;
   *  int printalign; */
  int pop;

  if ((rep < BURNIN + UPDATEFREQ) && (rep > BURNIN)) {
    printf ("\nBURNIN completed");
  }

  if (((NOADMIX) && (ADMBURNIN > 0)) && !LOCPRIOR ) {
    if ((rep < ADMBURNIN + UPDATEFREQ) && (rep >= ADMBURNIN)) {
      printf ("\nAdmixture Burnin complete.  Current alpha");
      if (POPALPHAS) {
        printf ("s = ");
        for (pop = 0; pop < MAXPOPS; pop++) {
          printf ("%1.3f ", Alpha[pop]);
        }
      } else {
        printf (" = %1.3f ", Alpha[0]);
      }
      printf ("\n\n");
    }
  }
  
  if ((LINKAGE) && (ADMBURNIN> 0)) {
    if ((rep < ADMBURNIN+UPDATEFREQ) && (rep >= ADMBURNIN)) {
      printf ("\nNo recombination Burnin complete.  Current rec = %1.3f",
              CalculateRAverage (R));
      PrintBanner (rep, Alpha, Fst, like, lambda);
    }
  }

  /*calculate some stuff */

  if (LINKAGE) {
    if (rep <= ADMBURNIN) {
      like = CalcLike (Geno, PreGeno, Q, P, Recessive, NULL, NULL); 
    } else {
      like = recomblikelihood;
    }

    if (rep >= BURNIN + 2) { /* +2 because need 2 observations for variance*/
      logprob = EstLogProb (sumlikes, sumsqlikes, rep - BURNIN);
    } else if (COMPUTEPROB) { /*not linkage model*/
      if (rep <= BURNIN + 2) {
        like = CalcLike (Geno, PreGeno, Q, P, Recessive, NULL, NULL);
      } else {
        logprob = EstLogProb (sumlikes, sumsqlikes, rep - BURNIN);
      }
    }
  }

  /*possibly print banner to indicate what the numbers are */
  if ((rep == UPDATEFREQ) ||
      ((rep < BURNIN) && ((rep % (BANNERFREQ * UPDATEFREQ)) == 0)) ||
      ((rep > BURNIN) && (((rep - BURNIN) % (BANNERFREQ * UPDATEFREQ)) == 0))) {
    if (rep != NUMREPS + BURNIN) {        /*don't bother for last line of output */
      if ((rep > 1) && (PRINTQSUM) && (MAXPOPS > 1)) {
        PrintMembership (stdout, Q, Individual);
      }
      PrintBanner (rep, Alpha, Fst, like, lambda);
    }
  }
  
  /*print current values to screen */
  printf ("%5d:    ", rep);

  if (LINKAGE && rep >= ADMBURNIN) {
    if (!INDIVIDUALR) {
      printf ("%1.09f  ", R[0]);
    } else {
      printf ("%1.09f   ", CalculateRAverage (R));
    }
  }

  if (PRINTLAMBDA) {
    if (POPSPECIFICLAMBDA) {
      for (pop=0;pop<MAXPOPS;pop++) {
        printf ("%1.2f    ", lambda[pop]);
      }
    } else {
      printf ("%1.2f    ", lambda[0]);
    }
  }

  if ((!(NOADMIX)) && (!(NOALPHA))) {
    if (POPALPHAS) {
      for (pop = 0; pop < MAXPOPS; pop++) {
        printf ("%1.3f  ", Alpha[pop]);
        if (pop > 8) {
          printf (" "); /*extra space for number */
        }
      }
    } else {
      printf ("%1.3f  ", Alpha[0]);
    }
  }
  
  if (FREQSCORR) {
    printf ("  ");
    if (ONEFST) {
      printf ("%1.3f ", Fst[0]);
    } else {
      for (pop = 0; pop < MAXPOPS; pop++) {
        printf ("%1.3f ", Fst[pop]);
      }
    }
    printf ("   ");
  } else {
    printf ("  ");
  }

  if (LOCPRIOR) {
    printf ("%1.3f ", LocPrior[0]);
    printf ("   ");
  }

  /*currently it only net distances not KLD ones */
  if (PRINTKLD || PRINTNET) {
    PrintNET (stdout, P, NumAlleles, 1, 0);
  }
  
  if (COMPUTEPROB) {
    if (rep > BURNIN + 2) {
      printf ("  %.0f  ", like);
      printf ("  %.0f ", logprob);
    } else {
      printf ("  --  ");
    }
  }

  printf ("\n");

  if (rep == BURNIN) {
    printf ("\nBURNIN completed");
    PrintBanner (rep, Alpha, Fst, like, lambda);
  }
  fflush(stdout);

}
/*----------------------------------------------------*/
void
PrintBanner (int rep, double *Alpha, double *Fst, double like, double *lambda)
    /*print banner to screen during run giving variable names */
{
  int i, j, k;
  int pop;

  printf ("\n");
  for (i = 4; i < ((int) log10 (rep)); i++) {
    printf (" ");
  }
  printf (" Rep#:   ");

  if ((LINKAGE) && (rep >= ADMBURNIN)) {
    printf (" r           ");
  }

  if (PRINTLAMBDA) {
    if (POPSPECIFICLAMBDA) {
      for (pop=0;pop<MAXPOPS;pop++) {
        printf("Lambda%d ",pop+1);
      }
    } else {
      printf ("Lambda  ");
    }
  }

  if (((!(NOADMIX)) && (!(NOALPHA)))) {
    printf (" ");
    if (POPALPHAS) {
      for (pop = 0; pop < MAXPOPS; pop++) {
        printf ("Alpha%d ", pop + 1);
      }
    } else {
      printf ("Alpha  ");
    }
  }

  if (FREQSCORR) {
    printf ("   ");
    if (ONEFST)
      printf ("Fst   ");
    else
      for (pop = 0; pop < MAXPOPS; pop++)
        printf (" F%d   ", pop + 1);
    printf (" ");
  }
  else
    printf (" ");

  if (LOCPRIOR)
  {
    printf ("  r     ");
  }

  if (PRINTKLD || PRINTNET)
  {
    for (j = 0; j < MAXPOPS - 1; j++)
      for (k = j + 1; k < MAXPOPS; k++)
        printf (" D%d,%d ", j + 1, k + 1);
  }

  if (COMPUTEPROB)
  {
    printf ("   Ln Like ");
    for (i = 8; i < ((int) log10 (rep)); i++)
      printf (" ");

    if (rep >= BURNIN)
      printf (" Est Ln P(D)");
  }
  printf ("\n");


}

/*----------------------------------------------------*/
double
EstLogProb (double sumlikes, double sumsqlikes, int reps)
    /*returns the current estimated Prob of Data.  Reps is the
      number of reps, not including burnin */
{
  double mean = sumlikes / reps;
  double var = SampleVar (sumsqlikes, sumlikes, reps);

  return (mean - var / 2.0);

}

/*-------------------------------------------------------*/
double
NETiv (int pop1, int pop2, double *P, int *NumAlleles, int reps)
{
  /* This function returns the current estimated average net nucleotide
     distance between the allele frequencies in populations 1 and 2.
     Here reps is the number of reps over which the P is an average, rather
     than the number of reps since the start of the program, as elsewhere */
  double sum, d1, d2,d12, norm;
  int loc,allele;
  norm = (double) reps;
  norm *= norm;
  sum=0.0;
  for (loc=0; loc<NUMLOCI;loc++) {
    d1=0.0;
    d2=0.0;
    d12=0.0;
    for (allele=0;allele<NumAlleles[loc];allele++) {
      d1+=P[PPos(loc,pop1,allele)]*P[PPos(loc,pop1,allele)]/norm;
      d2+=P[PPos(loc,pop2,allele)]*P[PPos(loc,pop2,allele)]/norm;
      d12+=P[PPos(loc,pop1,allele)]*P[PPos(loc,pop2,allele)]/norm;
    }
    sum+= 0.5*(d1+d2)-d12;
  }
  return sum/NUMLOCI;
}


/*-------------------------------------------------------*/
double
GROSSiv (int pop1, int pop2, double *P, int *NumAlleles, int reps)
{
  /* This function returns the current estimated average gross nucleotide
     distance between the allele frequencies in populations 1 and 2.
     Here reps is the number of reps over which the P is an average, rather
     than the number of reps since the start of the program, as elsewhere */
  double sum, d12, norm;
  int loc,allele;
  norm = (double) reps;
  norm *= norm;
  sum=0.0;
  for (loc=0; loc<NUMLOCI;loc++)
  {
    d12=0.0;
    for (allele=0;allele<NumAlleles[loc];allele++)
    {
      d12+=P[PPos(loc,pop1,allele)]*P[PPos(loc,pop2,allele)]/norm;
    }
    sum+= 1.0-d12;
  }
  return sum/NUMLOCI;
}

/*----------------------------------------------------*/
double
KLDiv (int pop1, int pop2, double *P, double *LogP, int *NumAlleles, int reps)
    /*This function returns the current (average) estimated
      Kullback-Leibler divergence between the allele frequencies in
      pops 1 and 2.  Here reps is the number of reps over which the P
      is an average, rather than the number of reps since the start of
      the program, as elsewhere. */
{
  double sum = 0.0;
  int allele, loc;

  for (loc = 0; loc < NUMLOCI; loc++)
    for (allele = 0; allele < NumAlleles[loc]; allele++)
      sum += ((double) P[PPos (loc, pop1, allele)] / reps)
          * log (P[PPos (loc, pop1, allele)] / P[PPos (loc, pop2, allele)]);

  return sum / NUMLOCI;

}

/*----------------------------------------------------*/
void
PrintNET (FILE * file, double *P, int *NumAlleles, int reps, int format)
    /*This function prints the current (average) estimated
      Net-nucleotide
