/*=======================================================

  STRUCTURE.C

  Program for inferring population structure using multilocus
  genotype data.

  Code written by Daniel Falush, Melissa Hubisz, and Jonathan Pritchard

  See additional details in README file.

  =========================================================*/
#define VERSION "2.3.4 (Jul 2012)"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "mymath.h"
#include "ran.h"
#include "params.h"
#include "datain.h"
#include "output.h"

void InitializeZ (int *Geno, struct IND *Individual, int *Z);
void UpdateQAdmixture (double *Q, int *Z, double *Alpha, struct IND *Individual);

/*========================================================
  ==========================================================*/
void
Welcome (FILE * file)
{
  fprintf (file, "\n\n");
  fprintf (file, "----------------------------------------------------\n");
  fprintf (file, "STRUCTURE by Pritchard, Stephens and Donnelly (2000)\n");
  fprintf (file, "     and Falush, Stephens and Pritchard (2003)\n");
  fprintf (file, "       Code by Pritchard, Falush and Hubisz\n");
  fprintf (file, "             Version %s\n", VERSION);
  fprintf (file, "----------------------------------------------------\n");

  fprintf (file, "\n\n");
  fflush(file);
}
/*-------------------------------------------------------*/
void Kill ()                            /*exit program */
{
  printf ("\nExiting the program due to error(s) listed above.\n\n");
  exit (1);
}
/*---------------------------------------*/
void CheckParamCombinations()
{
  if ((LINES>2) && (USEPOPINFO==1))
  {
    printf(" USEPOPINFO does not work for ploidy >2\n");
    Kill();
  }
  if ((PHASEINFO) && (LINES!=2))
  {
    printf("phase info is only applicable to diploid data!! \n");
    Kill();
  }
  if (LINES ==2 && PHASEINFO==1 && PHASED ==0 && MARKOVPHASE==-1 && LINKAGE)
  {
    printf("You need to specify a phase model using the parameter MARKOVPHASE!! \n");
    Kill();
  }
  if (LINKAGE && !MAPDISTANCES)
  {
    printf("Map distance information is required in order to run the linkage model. \n");
    Kill();
  }

  if ((LINKAGE) && (!PHASED) && (LINES!=2))
  {
    printf("unphased data only permitted for diploids!! \n");
    Kill();
  }

  if ((LINKAGE) && (NOADMIX))
  {
    printf("please choose the LINKAGE option or the NOADMIX option, but not both \n");
    Kill();
  }

  if ((INFERLAMBDA) && (FREQSCORR))
  {
    printf("Warning!, choosing both INFERLAMBDA and FREQSCORR parameters may leave the computer with insufficient information to estimate either parameter accurately \n");
  }


  if ((LINKAGE==0) && (SITEBYSITE))
  {
    printf("SITEBYSITE is not currently implemented for no-linkage model\n");
    SITEBYSITE=0;
  }


  if (((NOADMIX) || (LINKAGE)) && ADMBURNIN >= BURNIN)
  {
    printf("The ADMBURNIN should finish before the BURNIN!! \n");
    Kill();
  }

  if ((PFROMPOPFLAGONLY) && (!(POPFLAG)))
  {
    printf("PFROMPOPFLAGONLY can only be turned on when the data file contains POPFLAG data\n");
    Kill();
  }
  if (LINES>2 && RECESSIVEALLELES && MISSING==NOTAMBIGUOUS)
  {
    printf("The code for missingdata (MISSING) should be set differently to the code (NOTAMBIGUOUS) for loci whose genotypes are known");
    Kill();
  }
  if (LOCPRIOR && LINKAGE) {
    printf("LOCPRIOR model is not compatible with linkage model\n");
    Kill();
  }
  if (LOCPRIOR && USEPOPINFO) {
    printf("LOCPRIOR model is not compatible with USEPOPINFO option\n");
    Kill();
  }
  /*  if (RANDOMIZE && SEED!=-1) {
      printf("Warning: Seed from parameter file will not be used as RANDOMIZE is set to 1.  SEED in output file will be random seed drawn from time\n");
      }*/  /* modified by JKP */

  if (RANDOMIZE)
    printf("Note: RANDOMIZE is set to 1. The random number generator will be initialized using the system clock, ignoring any specified value of SEED.\n");

}
/*---------------------------------------*/
/*void FreeAll (int *Geno, double *Mapdistance, char *Markername,
              struct IND *Individual, int *Translation,
              int *NumAlleles, int *Z, int *Z1, double *Q, double *P, double *LogP, double *Epsilon,
              double *Fst, int *NumLociPop, double *PSum, double *QSum,
              double *SiteBySiteSum,
              double *FstSum, int *AncestDist, double *UsePopProbs, double *R,
              double *sumR, double *varR, double *LocPrior, double *sumLocPrior)
*/
void FreeAll(double *Mapdistance, double *Phase, int *Phasemodel, double *lambda, double *sumlambda,
             char *Markername, int *Geno, int* PreGeno, int* Recessive, struct IND *Individual,
             int *Translation, int *NumAlleles, int *Z, int *Z1, double *Q, double *P, double *LogP,
             double *R, double *sumR, double *varR, double *Epsilon, double *SumEpsilon, double *Fst,
             double *FstSum, int *NumLociPop, double *PSum, double *QSum, double *SiteBySiteSum,
             int *AncestDist, double *UsePopProbs, double *LocPrior, double *sumLocPrior,
             double *Alpha, double *sumAlpha, double *sumIndLikes, double *indLikesNorm)
{
  /** these variables are calloc'd in main and freed in the same order */

  free (Mapdistance);
  
  free(Phase);
  if (LINES==2 && PHASED == 0) {
    free(Phasemodel);
  }
  


  free(lambda);
  free(sumlambda);

  
  free (Markername);
  free (Geno);

  if (RECESSIVEALLELES) {
    free(PreGeno);
    free(Recessive);
  }


  free (Individual);
  free (Translation);
  free (NumAlleles);
  
  free (Z);
  free (Z1);
  
  free (Q);
  free (P);
  
  free (LogP);
  free (R);
  free (sumR);
  free (varR);

  free (Epsilon);
  
  if (FREQSCORR) {
    free(SumEpsilon);
  }
  
  
  free (Fst);
  free (FstSum);
  free (NumLociPop);
  
  free (PSum);
  free (QSum);

  if ( SITEBYSITE) {
    free (SiteBySiteSum);
  }

  if (ANCESTDIST) {
    free (AncestDist);
  }

  if (USEPOPINFO)  {
    free (UsePopProbs);
  }

  if (LOCPRIOR) {
    free(LocPrior);
    free(sumLocPrior);
  }

  
  free(Alpha);
  free(sumAlpha);

  free(sumIndLikes);
  free(indLikesNorm);
}
/*---------------------------------------*/
void
PrintTranslation (int *Translation, int *NumAlleles)
{
  int loc, allele;
  printf ("Translation matrix:\n");
  for (loc = 0; loc < NUMLOCI; loc++)
  {
    for (allele = 0; allele < NumAlleles[loc]; allele++)
      printf ("%2d ", Translation[TransPos (loc, allele)]);
    printf ("\n");
  }
  printf ("\n");
  for (loc = 0; loc < NUMLOCI; loc++)
    printf ("%d ", NumAlleles[loc]);
  printf ("\n");


}

/*-------------------------------------------*/
void
InitFromGeogPop (int *Geno, struct IND *Individual, int *Z, int verbose)
{
  /*initialize the population of origin of each allele. These are
    started in their given populations (when there is population info).
    It is assumed that the populations are numbered 1..K.  If the
    population identifier is out of range, the allele is assigned to a
    random population.  These is used to check that the MCMC is not
    producing unexpected results because it couldn't find the mode. */

  int ind, line, loc;
  int poperrors = 0;
  int pop;

  if (!(POPDATA)) {
    InitializeZ (Geno, Individual, Z);
    if (verbose) {
      printf ("Starting from a random configuration because POP=0\n");
    }
  } else {
    for (ind = 0; ind < NUMINDS; ind++) {
      for (line = 0; line < LINES; line++) {
        for (loc = 0; loc < NUMLOCI; loc++) {
          if (Geno[GenPos (ind, line, loc)] == MISSING) {
            Z[ZPos (ind, line, loc)] = UNASSIGNED;
          } else {
            pop = Individual[ind].Population;
            if ((pop > 0) && (pop <= MAXPOPS)) {
              Z[ZPos (ind, line, loc)] = pop - 1;
            } else {
              Z[ZPos (ind, line, loc)] = RandomInteger (0, MAXPOPS - 1);
              poperrors++;
            }
          }
        }
      }
    }
    if (verbose) {
      printf ("USING GIVEN POPULATION INFO TO SET INITIAL CONDITION\n");
    }
    if ((verbose) && (poperrors)) {
      printf ("WARNING: unable to initialize %d individuals to the predefined\n", poperrors);
      printf ("populations because their population names were not in the range\n");
      printf ("{1..%d}.  These individuals were initialized at random.\n",MAXPOPS);
    }
  }
}
/*---------------------------------------*/
void
InitializeZ (int *Geno, struct IND *Individual, int *Z)
{
  /*initialize the population of origin of each allele. I pick these at
    random because this seems to produce better behaviour than starting
    out with everybody in one pop, for instance.  I also set missing Data
    to the unassigned pop from the start.  STARTATPOPINFO indicates that
    individuals should be started at their given populations. It is
    assumed that the populations are numbered 1..K.  If the population
    identifier is out of range, the allele is assigned to a random
    population. */

  int ind, line, loc;
  int allele;
  int pop;

  for (ind = 0; ind < NUMINDS; ind++) {
    for (line = 0; line < LINES; line++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        allele = Geno[GenPos (ind, line, loc)]; /*missing data */
        if (allele == MISSING) {
          Z[ZPos (ind, line, loc)] = UNASSIGNED;
        } else {  /*------data present-----------*/
          if ((STARTATPOPINFO) && (POPDATA)) {    /*use given pops as initial Z */
            pop = Individual[ind].Population;
            if ((pop > 0) && (pop <= MAXPOPS)) {
              Z[ZPos (ind, line, loc)] = pop - 1;
            } else {
              Z[ZPos (ind, line, loc)] = RandomInteger (0, MAXPOPS - 1);
            }
          } else {          /*initial Z random */
            Z[ZPos (ind, line, loc)] = RandomInteger (0, MAXPOPS - 1);
          }
        }
        /*printf("%d ",Z[ZPos(ind,line,loc)]); */
      }
      /*printf("\n"); */
    }
  }

  if ((STARTATPOPINFO) && (POPDATA)) {
    printf ("USING GIVEN POPULATION INFO TO SET INITIAL CONDITION\n");
  }
}
/*---------------------------------------*/
void
InitFreqPriors (double *Epsilon, double *Fst, int *Geno, int *NumAlleles)
{
  int ind, line, loc, allele,pop;
  int value;
  int *Count /*[MAXALLELES] */ ;        /*stores number of copies of each allele
                                          at present locus */
  int total;

  if (!(FREQSCORR)) {           /*allele frequencies uncorrelated across populations */
    for (loc = 0; loc < NUMLOCI; loc++) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        Epsilon[EpsPos (loc, allele)] = LAMBDA;
      }
    }
    for (pop = 0; pop< MAXPOPS; pop++) {
      Fst[pop] = 0.5;
    }
  } else {                      /*correlated allele frequencies------------------- */
    Count = calloc (MAXALLELES, sizeof (int));
    if (Count == NULL) {
      printf ("Error in assigning memory, InitFreqPriors\n");
      Kill ();
    }

    for (loc = 0; loc < NUMLOCI; loc++) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        Count[allele] = 0;
      }
      total = 0;
      for (ind = 0; ind < NUMINDS; ind++) {
        for (line = 0; line < LINES; line++) {
          value = Geno[GenPos (ind, line, loc)];
          if (value != MISSING) {
            total++;
            if ((value < 0) || (value >= NumAlleles[loc])) {
              printf ("WARNING: expected allele value, InitFreqPriors: loc %d, allele %d\n", loc, value);
            } else {
              Count[value]++;
            }
          }
        }
      }

      /*Start Epsilon at (mean sample frequencies) */
      /* add lambda to all counts to ensure positive frequencies
       * for recessive model etc */
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        Epsilon[EpsPos (loc, allele)] =
            (((double) LAMBDA +
              (double) Count[allele]) / ((double) NumAlleles[loc] *
                                         (double) LAMBDA +
                                         (double) total));
        /* ((double) Count[allele] / total); */
      }
      /*printf("\n"); */
    }


    for (pop= 0; pop < MAXPOPS; pop++) {
      Fst[pop] = FPRIORMEAN;    /*start Fst at the prior mean */
    }
    if (Count != NULL) {
      free (Count);
    }
  }                             /*end, correlated allele frequencies------------- */

}
/*---------------------------------------*/


void
CheckPopPriors (struct IND *Individual)
    /*This function is called when USEPOPINFO==1 to check that the
      prior information on populations is ok */
{
  int ind;
  int numnopop;

  if ((MIGRPRIOR < 0.0) || (MIGRPRIOR > 1.0)) {
    printf ("MIGRPRIOR (which is currently set to %1.3f) must be in the range [0.0, 1.0]\n", MIGRPRIOR);
    Kill ();
  }

  if (!(POPDATA)) {
    printf ("Can't apply USEPOPINFO because no POPDATA in input data file\n");
    Kill ();
  }

  if (!(POPFLAG)) {              /*if no popflag, then assume that everybody should be used */
    for (ind = 0; ind < NUMINDS; ind++) {
      Individual[ind].PopFlag = 1;
    }
  }

  /*Check that the given population is within range for all individuals, and if
    not, then turn popflag off for that individual. */
  for (ind = 0; ind < NUMINDS; ind++) {
    if (Individual[ind].PopFlag) {
      if ((Individual[ind].Population < 1) || (Individual[ind].Population > MAXPOPS)) {
        printf ("Warning: population prior for individual %d is %d, which is not\n",
                ind + 1, Individual[ind].Population);
        printf ("  in the range 1..%d.  Population prior for this individual will\n",
                MAXPOPS);
        printf ("  be ignored\n");
        Individual[ind].PopFlag = 0;
      }
    }
  }

  if ((INFERALPHA) && (!(NOADMIX))) {     /*check whether alpha is needed at all */
    numnopop = 0;
    for (ind = 0; ind < NUMINDS; ind++) {
      if (Individual[ind].PopFlag == 0) {
        numnopop++;
      }
    }
    if (numnopop == 0) {
      NOALPHA = 1;
      INFERALPHA = 0;
    }
  }
}


/*GetNumLocations: Melissa added 7/12/07.  Sets the variable NUMLOCATIONS and also setse
  all the individuals so that ind[i].myloc is in (0...NUMLOCATIONS).  ind[i].loc is unchanged
  to whatever the input file indicates*/
void GetNumLocations (struct IND *ind) {
  int maxloc=0, i, j, *freq, pos;
  for (i=0; i<NUMINDS; i++) {
    /* for now we're not dealing with unknown location */
    if (ind[i].Location < 0) {
      printf("individual %s has negative location!  locations should be >= 0\n", ind[i].Label);
      Kill();
    }
    if (ind[i].Location > maxloc) {
      maxloc = ind[i].Location;
    }
  }

  freq = malloc((maxloc+1)*sizeof(int));
  for (i=0; i<=maxloc; i++)
