/*
   Part of structure.c.  

   This bit is in charge of reading in the information from the datafile, and 
   preparing it. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structure.h"
#include "params.h"

FILE *INPUT;			/*Data file */

void OpenData ();
void ReadData (int *Geno, double *Mapdistance, char *Markername, struct IND *Individual, double *Phase,
                int *Recessive);
void WarnEOF (int ind);
int CheckIfValidInt (const char intstr[ALLELLEN], int ind, int loc);
void Mismatch (int ind, const char description[30]);

void ExtraData ();
void PrintSomeData (int *Geno, struct IND *Individual, FILE * file);
void CountLineLens ();
void RecessiveDataCheck(int *Geno, int *Recessive);


/*================================================*/
void 
ReadInputFile (int *Geno, double *Mapdistance, char *Markername, struct IND *Individual, double *Phase,
                int *Recessive)
{
  OpenData ();			/*open input file */
  printf ("Reading file \"%s\".\n", DATAFILE);
  ReadData (Geno, Mapdistance, Markername, Individual, Phase,Recessive);
  fclose (INPUT);
  if (ECHODATA)
    PrintSomeData (Geno, Individual, stdout);
  if (RECESSIVEALLELES) RecessiveDataCheck(Geno,Recessive);

}
/*-------------------------------------------*/
void 
OpenData ()
{
  /*  int trouble = 0; */
  INPUT = fopen (DATAFILE, "r");

  if (INPUT == NULL)
    {
      printf ("Unable to open the file %s.\n", DATAFILE);
      Kill ();
    }

}


/*-------------------------------------------*/
void 
ReadData (int *Geno, double *Mapdistance, char *Markername, struct IND *Individual, double *Phase, int *Recessive)
{
  int ind;
  int loc;
  int line;
  int col;
  char label[LABELLEN];
  char intstr[ALLELLEN];
  char gene[GENELEN];
  int next;
  int pop=-1;
  int popflag=-1;
  int phen=-1;
  int strlength;
  int valid = 1;
  /*  int polymorphismcounter; */
  /* int toggle; */
  int rowsperind;  /*the next 2 variables allow single-row input for data*/
  int colsperloc;  /*number of columns per locus = 1 or PLOIDY depending on format*/ 
  int i,j;
  /*  int loc1,loc2,temp; */
  int *inputorder = malloc(NUMLOCI*sizeof(int));
  int location=-1;
  int phenotypecol=-1;
  /*  int numcols; */
  int CheckIfValidDouble (const char intstr[ALLELLEN], int ind,int loc);


  ind = -1;

 /*inputorder is initialised with the loci in the right order */
  for (loc=0;loc<NUMLOCI;loc++)
  inputorder[loc]=loc;

 /* The following part (which should normally be commented out) randomizes the order of the loci*/
 /* while leaving mapdistances, etc, intact. */
 /* this can be useful when testing if there is significant map distance info in the sample. */
 /*
   for (ind=0;ind<NUMLOCI*100;ind++)
 {
 loc1=RandomInteger(0,NUMLOCI-1);
 loc2=RandomInteger(0,NUMLOCI-1);
 temp=inputorder[loc1];
 inputorder[loc1]=inputorder[loc2];
 inputorder[loc2]=temp;
 }
 */
   /*Read in locus information*/
   
     if (MARKERNAMES) {
       for (loc = 0; loc < NUMLOCI; loc++)
       {
 	strlength = ReadString (gene, GENELEN, INPUT);
 	if (strlength == 0) {
 	  printf("%i MARKERNAMES\n", loc);WarnEOF (ind);}
 	for (j=0; j<strlength; j++)
 	  Markername[MarkernamePos(loc,j)] = gene[j];
       }
     }
   if (RECESSIVEALLELES)
     for (loc = 0; loc < NUMLOCI; loc++)
       {
 	strlength = ReadString (intstr, ALLELLEN, INPUT);
 	if (strlength == 0)
 	  WarnEOF (ind);
 	Recessive[loc] = (double) atoi (intstr);
       }
   if (MAPDISTANCES)
     for (loc = 0; loc < NUMLOCI; loc++)
       {
 	strlength = ReadString (intstr, ALLELLEN, INPUT);
 	if (strlength == 0)
 	  WarnEOF (ind);
 	
 	/* daniel made the change: change from atoi to atof on Dec 24, 2002 */
 	Mapdistance[loc] = (double) atof (intstr);
       }
   /*End reading in locus information-------------------------------
     Read in individual/genotype data-------------------------------*/
   
   if (ONEROWPERIND) {rowsperind = 1; colsperloc = LINES;}
   else {rowsperind = LINES; colsperloc = 1;}
   
   for (ind = 0; ind < NUMINDS; ind++)
     {
       for (line = 0; line < rowsperind; line++)
 	{
 	  if (LABEL)		/*read in extra data for each individual */
 	    {
 	      strlength = ReadString (label, LABELLEN, INPUT);
 	      if (strlength == 0)
 		WarnEOF (ind);
 	    }
 	  if (POPDATA)
 	    {
 	      strlength = ReadString (intstr, ALLELLEN, INPUT);
 	      if (strlength == 0)
 		WarnEOF (ind);
 	      valid = valid * CheckIfValidInt (intstr, ind, loc);
 	      pop = atoi (intstr);
 	    }
 	  if (POPFLAG)
 	    {
 	      strlength = ReadString (intstr, ALLELLEN, INPUT);
 	      if (strlength == 0)
 		WarnEOF (ind);
 	      valid = valid * CheckIfValidInt (intstr, ind, loc);
 	      popflag = atoi (intstr);
 	    }
 	  /* melissa added 7/12/07 */
 	  if (LOCDATA==1) {
 	    strlength = ReadString(intstr, ALLELLEN, INPUT);
 	    if (strlength ==0)
 	      WarnEOF(ind);
 	    valid = valid + CheckIfValidInt(intstr, ind, loc);
 	    location = atoi(intstr);
 	  }
 	  if (PHENOTYPE) /*-------------------------*/
 	    {
 	      if (PHENOTYPECOL == -9) phenotypecol = 0;   /*if STRATPARAMS is not read, then
 							    this is given a value of -9*/
 
 	      else if (PHENOTYPECOL <= LABEL + POPDATA + POPFLAG)  /*too small*/
 		{
 		  printf("Error: PHENOTYPECOL (now set in STRATPARAMS as %d) must be at least %d\n",
 			PHENOTYPECOL, LABEL + POPDATA + POPFLAG + 1);
 		  printf("given your current values of LABEL, POPDATA and POPFLAG\n");
 		  Kill();
 		}
 	      else if (PHENOTYPECOL > LABEL + POPDATA + POPFLAG + 1 + EXTRACOLS)  /*too large*/
 		{
 		  printf("Error: PHENOTYPECOL (set in STRATPARAMS as %d) must be at most %d,\n",
 			 PHENOTYPECOL, LABEL + POPDATA + POPFLAG + 1 + EXTRACOLS);
 		  printf("  given your current values of LABEL, POPDATA, POPFLAG and EXTRACOLS\n");
 		  Kill();
 		}
 	      else phenotypecol = PHENOTYPECOL - LABEL - POPDATA - POPFLAG -1;   /*figure out which
 									           column to read*/
 
 	      /*printf("phenotypecol = %d; EXTRACOLS = %d\n",phenotypecol, EXTRACOLS);*/
 	      for (col = 0; col < EXTRACOLS + 1; col++) 
 		{
 		  if (col == phenotypecol)
 		    {
 		      strlength = ReadString (intstr, ALLELLEN, INPUT);
 		      if (strlength == 0)
 			WarnEOF (ind);
 		      valid = valid * CheckIfValidInt (intstr, ind, loc);
 		      phen = atoi (intstr);
 		      /*printf("phenotypecol = %d, col = %d, phen = %d\n",phenotypecol,col,phen);*/
 		    }
 		  else ReadString (intstr, ALLELLEN, INPUT);   /*skip these data*/
 		}
 	    
 	    }          /*---End if (Phenotype)----------------------*/
 
 	  else for (col = 0; col < EXTRACOLS; col++)  /*no phenotypes*/
 	    ReadString (intstr, ALLELLEN, INPUT);
 
 	  
 	  if (line == 0)		/*save the preliminary data */
 	    {
 	      if (LABEL)
 		strcpy (Individual[ind].Label, label);
 	      if (POPDATA)
 		Individual[ind].Population = pop;
 	      if (POPFLAG)
 		Individual[ind].PopFlag = popflag;
 	      if (PHENOTYPE)
 		Individual[ind].Phenotype = phen;
 
 	      /* melissa added 7/12/07 */
 	      if (LOCDATA)
 		Individual[ind].Location = location;
 	      if (LOCISPOP) 
 		Individual[ind].Location = pop;
 	    }
 	  if (line > 1)		/*check for consistency across lines */
 	    {
 	      if (LABEL)
 		if ((strcmp (Individual[ind].Label, label)))
 		  {
 		    Mismatch (ind, "label");
 		    valid = 0;
 		  }
 	      /*printf("The labels are %s and %s\n",Individual[ind].Label,label);} */
 	      if (POPDATA)
 		if (Individual[ind].Population != pop)
 		  {
 		    Mismatch (ind, "pop");
 		    valid = 0;
 		  }
 	      if (POPFLAG)
 		if (Individual[ind].PopFlag != popflag)
 		  {
 		    Mismatch (ind, "popflag");
 		    valid = 0;
 		  }
 	      if (PHENOTYPE)
 		if (Individual[ind].Phenotype != phen)
 		  {
 		    Mismatch (ind, "phenotype");
 		    valid = 0;
 		  }
 	    }
 	  
 	  for (loc = 0; loc < NUMLOCI; loc++)    /*read in genotype data here*/
 	    for (i=0; i<colsperloc; i++)
 	      {
 		strlength = ReadString (intstr, ALLELLEN, INPUT);
 		if (strlength == 0) {printf("readlociEOF\n");
 		WarnEOF (ind);}
 		valid = valid * CheckIfValidInt (intstr, ind, loc);
 		if (ONEROWPERIND) Geno[GenPos (ind, i, inputorder[loc])] = atoi (intstr);
 		else Geno[GenPos (ind, line, inputorder[loc])] = atoi (intstr);
 	      }
 	/* printf(" % 4ld % 4ld % 4ld     .....    % 4ld % 4ld \n",Geno[GenPos(ind,line,0)],Geno[GenPos(ind,line,1)],Geno[GenPos(ind,line,2)],Geno[GenPos(ind,line,NUMLOCI-2)],Geno[GenPos(ind,line,NUMLOC[...]
 	}
 
       if (PHASEINFO) /*possibly read in row of phase information*/
 	for (loc=0;loc<NUMLOCI;loc++)
 	  {
 	    strlength = ReadString (intstr, ALLELLEN, INPUT);
 	    if (strlength == 0) WarnEOF (ind);
 	    valid=valid * CheckIfValidDouble(intstr,ind,loc);
 	    Phase[PhasePos(ind,loc)]=atof(intstr);
 	    /*check that probabilities are in [0,1]*/
 	    if (Phase[PhasePos(ind,loc)]>1.0 || Phase[PhasePos(ind,loc)]<0.0)
 	      {
 		printf("Phase information for individual %d locus %d (%1.3f) is not a real in [0.0,1.0] \n",ind,loc,Phase[PhasePos(ind,loc)]);
 		valid=0;
 	      }
 	  }
     } /*end of loop over individuals*/


   /*check if anything else left in the input file*/
   do
     {
       next = fgetc (INPUT);
       if ((!Whitespace (next)) && (next != EOF))
 	{
 	  ExtraData ();
 	  valid = 0;
 	  break;
 	}
     }
   while (next != EOF);

   if (!(valid))
     {
       CountLineLens ();
       Kill ();
     }
   free(inputorder);
}
/*------------------------------------*/
void RecessiveDataCheck(int *Geno, int *Recessive)
    /* this function checks whether any genotypes have both a recessive and dominant
     * allele.  If any do, it terminates the program.
     * in the polyploid case, it also checks whether the NOTAMBIGUOUS code is anywhere 
     * in the datset. If it is it also terminates the program. */
{
  int ind, line, loc;
  int rec, dom;
  int error=0;

  int recessive_shown = 0;
  int recessive_allele_shown = 0;
  int recessive_allele_not_shown = 0;
  int unambiguous_loci = 0;
  int ambiguous_norecessive = 0;
  for (loc=0; loc<NUMLOCI; loc++)
    {
      
      if (Recessive[loc] != MISSING)
 	{
 	  recessive_shown=0;
 	  for (ind=0; ind<NUMINDS; ind++)
 	    {
 	      rec=0; 
 	      dom=0;
 	      for (line=0; line<LINES; line++)
 		{
 		  if (Geno[GenPos (ind,line,loc)] == Recessive[loc]) {
 		    rec=1;
 		    recessive_shown=1;
 		  }
 		  else if (Geno[GenPos (ind,line,loc)] != MISSING) dom=1;
 		}
 	      if (rec*dom==1 && error<100) 
 		{
 		  printf("WARNING: both recessive and dominant alleles in genotype: ind=%d, loc=%d\n",ind+1,loc+1);
 		  error++;
 		}

 	    }
 	  
 	  if(recessive_shown) 
 	    recessive_allele_shown++;
 	  else
 	    if(Recessive[loc]!= NOTAMBIGUOUS)
 	      recessive_allele_not_shown++;
 	}
    }
  if (error>100) printf("Total of %d such errors\n",error);
  if (error>0) {printf("Terminating program.\n"); Kill();}
  if (LINES>2)
    for (loc=0;loc<NUMLOCI;loc++){
      
      if(Recessive[loc] == NOTAMBIGUOUS)
 	unambiguous_loci++;
      if(Recessive[loc] == MISSING)
 	ambiguous_norecessive++;
      for (ind=0;ind<NUMINDS;ind++)
 	for (line=0;line<LINES;line++)
 	  if (Geno[GenPos(ind,line,loc)]==NOTAMBIGUOUS)
 	    {
 	      printf("WARNING: the code for NOTAMBIGUOUS alleles, %d, appears  at least once in the dataset at : ind=%d, loc=%d\n", NOTAMBIGUOUS, ind+1,loc+1);
 	      Kill(); /* modify by William - kill (small case) is not a standard ANSI C function */
 	    }
   
    }
  printf("\n");
  if(LINES>2){
    printf("Number of loci without genotypic ambiguity: %d \n",unambiguous_loci);
    printf("Number of loci with ambiguous copy numbers but no recessive alleles: %d \n",ambiguous_norecessive);
  }
  printf("Number of loci with recessive allele present in data: %d\n", recessive_allele_shown);
  printf("Number of loci with recessive allele absent in data: %d\n", recessive_allele_not_shown);
  fflush(stdout);

  
}
/*------------------------------------*/
void 
WarnEOF (int ind)
{
  /*This function is called from ReadData if there is an unexpected end of file */

  printf ("\n\nWARNING:  Unexpected end of input file.  The details of the\n");
  printf ("input file are set in mainparams.  I ran out of data while reading\n");
  printf ("the data for individual %d.\n\n", ind + 1);

  CountLineLens ();
  Kill ();
}
/*------------------------------------*/
int 
CheckIfValidInt (const char intstr[ALLELLEN], int ind, int loc)
{
  /*This function checks for non-numeric data in the input file (not
     used when reading the labels).  Returns 1 if valid, otherwise 0. */

  int i;
  int ok = 1;

  for (i = 0; i < ALLELLEN; i++)
    {
      if (intstr[i] == '\0')
 	break;
      if (((intstr[i] < '0') || (intstr[i] > '9')) && (intstr[i] != '-'))
 	ok = 0;
    }

  if (ok == 0)
    {
      printf ("\nWARNING! Probable error in the input file.  \n");
      printf ("Individual %d, locus %d;  ", ind + 1, loc + 1);
      printf ("encountered the following data \n");
      printf
