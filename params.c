/*Part of structure.c.  

This bit is in charge of reading in the parameters from the files
mainparams and extraparams.  It begins by defining default values for
everything (just in case!) and then scans through these files to get
the correct values.  All the parameters are stored in structure.c as
global variables.

Command line flags: enter the appropriate flag followed by new value :

-m mainparams file
-e extraparams file
-s stratparams file
-i input file 
-o output file
-K MAXPOPS 
-L NUMLOCI
-N NUMINDS
-D seeD

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "structure.h"
#include <string.h>

#define MAXNAME 30  /*maximum length of parameter names*/
FILE *PARAMS; /*parameter files*/


void ReadFile();
void GetParamValue();
void SetValue(char Word[]);
void OpenFile(const char input[15]);
void PresetValues();
void CheckIfMissing();
int Whitespace(char next);
void InputFileNames(int argc, const char *argv[], char filename[], const char flag[50]);
void CommandLineValues(int argc, const char *argv[]);


/*=======MAIN FUNCTION FOR THIS PART OF THE PROGRAM============*/
void GetParams(int strat, int argc, const char *argv[])
     /*the value of strat is 1 if this is called from STRAT.c, otherwise 0*/
     /*This function takes in all the values of the parameters.  The default
       is to get these all from main- extra- and STRATparams, but they can
       be taken from other files, depending on command line flags, or from
       the command line for a few special parameters*/
{
  
  /* William Wen made change: length from 50 to 500, avoid buffer overflow
   * in case user choose a very long path/file name. */
  char filename[500]; 
  int SetUpdateFreq();

  PresetValues();

  InputFileNames(argc,argv,filename,"-m");  /*check for command line change to input file*/
  OpenFile(filename);
  ReadFile();
  fclose(PARAMS);

  InputFileNames(argc,argv,filename,"-e");
  OpenFile(filename);
  ReadFile();
  fclose(PARAMS);

  if (strat)  /*called from STRAT.c*/
    {
      InputFileNames(argc,argv,filename,"-s");
      OpenFile(filename);
      ReadFile();
      fclose(PARAMS);
    }

  CommandLineValues(argc,argv);   /*check for command line change to parameters*/
  CheckIfMissing(strat);
  UPDATEFREQ = SetUpdateFreq();
} 
/*-------------------------------------------*/
void CommandLineValues(int argc, const char *argv[])
{
  char value[100];
  int i;
  int CommandLineInts(const char *value, const char parameter[100]);
  /*
-i input file 
-o output file
-K MAXPOPS 
-L NUMLOCI
-N NUMINDS
  */

  for (i=0; i<argc-1; i++)
    {
      if (strcmp(argv[i],"-i")==0) sprintf (DATAFILE,"%s",argv[i+1]); 
      else if (strcmp(argv[i],"-o")==0) sprintf (OUTFILE,"%s",argv[i+1]); 
      else if (strcmp(argv[i],"-K")==0) 
	{sprintf (value,"%s",argv[i+1]); MAXPOPS = CommandLineInts(value,"MAXPOPS");}
      else if (strcmp(argv[i],"-L")==0)  
	{sprintf (value,"%s",argv[i+1]); NUMLOCI = CommandLineInts(value,"NUMLOCI");}
      else if (strcmp(argv[i],"-N")==0)  
	{sprintf (value,"%s",argv[i+1]); NUMINDS = CommandLineInts(value,"NUMINDS");}	
      else if (strcmp(argv[i],"-r")==0)  
	{sprintf (value,"%s",argv[i+1]); NUMSIMSTATS = CommandLineInts(value,"NUMSIMSTATS");}	
      else if (strcmp(argv[i],"-c")==0)  
	{sprintf (value,"%s",argv[i+1]); PHENOTYPECOL = CommandLineInts(value,"PHENOTYPECOL");}
      else if (strcmp(argv[i], "-D")==0) 
	SEED = atoi(argv[++i]);
    }
}
/*-------------------------------------------*/
int CommandLineInts(const char *value, const char parameter[100])
     /*the point of this function is to take the values entered at the command
       line (for MAXPOPS, etc) and make sure they are valid integers*/
{
  int i;
  int ok = 1;

  for (i = 0; i < 100; i++)
    {
      if (value[i] == '\0')
	break;
      if (((value[i] < '0') || (value[i] > '9')) && (value[i] != '-'))
	{
	  ok = 0;
	  printf("Error in command-line data entry for %s\n",parameter);
	  printf("Value input at command line is %s\n",value);
	  Kill();
	}
    }
  
  return atoi(value);
}
/*-------------------------------------------*/
void InputFileNames(int argc, const char *argv[], char filename[], const char flag[50])
     /*This function checks whether there is a command line flag to
     indicate that a different input file should be used instead of
     mainparams, extraparams or STRATparams*/
{
  int found = 0;
  int i,j;

  int shift = 0;

  memset(filename,0,500);
  for (i=0; i<argc-1; i++)
    if (strcmp(argv[i],flag)==0) 
      {
	found = 1;
	
	/* modify by William to enable user to use space character in path/file name
         * Jan 29, 2003 */
	j=i+1;
	while(argv[j][0]!='-'&&j<argc){
	  strcat(filename,argv[j]);
	  strcat(filename," ");
	  j++;
	  if(j==argc)
	    break;
	}
	filename[strlen(filename)-1]='\0';
	
	/* sprintf (filename,"%s",argv[i+1]); */
        /* end of modification */
	
	break;
      }

  if (!(found)) 
    {
      if (strcmp(flag,"-m")==0) sprintf (filename, "%s","mainparams");
      else if (strcmp(flag,"-e")==0) sprintf (filename, "%s","extraparams");
      else if (strcmp(flag,"-s")==0) sprintf (filename, "%s","STRATparams");
    }
  

  /* modify by William to enable user to use space character in path/file name
   * Jan 29, 2003 */
 
  for(i=0;i<strlen(filename);i++){
    if(filename[i]=='\"'){           /*" [fix highlighting*/
      shift++;
    }
    filename[i]=filename[i+shift];
    if(filename[i]=='\"')          /*" [fix highlighting*/
      filename[i]='\0';
    if(filename[i]=='\0')
      break;
    
  }

  

}
/*-------------------------------------------*/
void ReadFile()
{
  /*The files may contain comments etc, but I look for strings in the
    form "#define NUMINDS 424" which define parameter values*/

  int notEOF = 1;  /*switch this to 0 when we hit EOF*/
  int next;

  do         /*scan the file for # until EOF is reached*/
    {
      next = fgetc(PARAMS);  
      if (next==EOF) notEOF=0;      
      else if (next=='#') /*signals parameter value coming*/
	GetParamValue(); 
    } while (notEOF);
}
/*-------------------------------------------*/
void GetParamValue()
{
  char NextWord[MAXNAME];
  int  ReadString(char Word[],int maxlength,FILE *THEFILE);

  if (ReadString(NextWord,MAXNAME,PARAMS))
    if (!(strcmp(NextWord,"define")))  /*check that next word is "define"*/
                                  /*strcmp returns 0 if strings unequal*/ 
      {
	if (ReadString(NextWord,MAXNAME,PARAMS))
	  SetValue(NextWord);
      }
}
/*-------------------------------------------*/
int  ReadSpaceString(char Word[],int maxlength,FILE *THEFILE)
    /* Add in by William Wen, for reading space-allowed path/file name
     * read through the string, includs all intermediate spaces, terminate
     * reading at "//" sign or new line. */
{
  int next;
  int last=0;
  int length = 0;
  /*  int white; */ /* 1 if next is whitespace*/
  int temp_length;
  
  while (length<maxlength-1) 
    {
      next = fgetc(THEFILE);
      if (next==EOF) break;    /*break if EOF*/
      
      if(last=='/'&&next=='/'){
	length--;
	break;
      }
      
      if(last==0&&isspace(next))
	continue;
      
      if (next!='\n'&&next!='\r')
	{
	  Word[length] = next;
	  length++;
	  last = next;
	}
      
      
      else if (length>0) break;  /*break if new whitespace*/
    }  
  
  temp_length=length-1;
  while(isspace(Word[temp_length]))
    temp_length--;
  
  length= temp_length+1;
	
 
  if (length>0)   /*add end position*/
    {
  
      
      Word[length] = '\0';
      length++;
    }
  if (length==maxlength)  /*move cursor to next whitespace*/
    {
      do 
	next = fgetc(THEFILE);
      while ((next != EOF) && (!(Whitespace(next)))); 
    }
    
  return (length);	
}	
/*-------------------------------------------*/

int  ReadString(char Word[],int maxlength,FILE *THEFILE)
     /*returns the length of the string.  The next word (continuous
     non-whitespace characters) is put into Word. Returns 0 if EOF if
     reached before any words.  If string is too long, moves cursor
     to next whitespace anyway.*/
{
  int next;
  int length = 0;
  int white; /* 1 if next is whitespace*/

  while (length<maxlength-1) 
    {
      next = fgetc(THEFILE);
      if (next==EOF) break;    /*break if EOF*/
      white = Whitespace(next);
      if (!(white))
	{
	  Word[length] = next;
	  length++;
	}
      else if (length>0) break;  /*break if new whitespace*/
    }  
    
  if (length>0)   /*add end position*/
    {
      Word[length] = '\0';
      length++;
    }
  if (length==maxlength)  /*move cursor to next whitespace*/
    {
      do 
	next = fgetc(THEFILE);
      while ((next != EOF) && (!(Whitespace(next)))); 
    }
    
  return (length);	
}	
/*-------------------------------------------*/


void SetValue(char Word[])
{
  if (VERBOSE) printf("Reading value of \"%s\"\n",Word);


  /* William modify this to fit the requirement of windows,
   * where space is allowed in path and file name. */
  if (!(strcmp(Word,"INFILE"))) {ReadSpaceString(DATAFILE,STRLEN,PARAMS);
  
  printf("datafile is\n%s\n",DATAFILE);
  }
  else if (!(strcmp(Word,"OUTFILE"))) ReadSpaceString(OUTFILE,STRLEN,PARAMS);
  else if (!(strcmp(Word,"NUMINDS"))) fscanf(PARAMS,"%d",&NUMINDS);
  else if (!(strcmp(Word,"NUMLOCI"))) fscanf(PARAMS,"%d",&NUMLOCI);
  else if (!(strcmp(Word,"MISSING"))) fscanf(PARAMS,"%d",&MISSING);
  else if (!(strcmp(Word,"LABEL"))) fscanf(PARAMS,"%d",&LABEL);
  else if (!(strcmp(Word,"POPDATA"))) fscanf(PARAMS,"%d",&POPDATA);
  else if (!(strcmp(Word,"MARKERNAMES"))) fscanf(PARAMS,"%d",&MARKERNAMES);
  else if (!(strcmp(Word,"PHASEINFO"))) fscanf(PARAMS,"%d",&PHASEINFO);
  else if (!(strcmp(Word,"MAPDISTANCES"))) fscanf(PARAMS,"%d",&MAPDISTANCES);
  else if (!(strcmp(Word,"MARKOVPHASE"))) fscanf(PARAMS,"%d",&MARKOVPHASE);
  else if (!(strcmp(Word,"PLOIDY"))) fscanf(PARAMS,"%d",&LINES);
  else if (!(strcmp(Word,"ONEROWPERIND"))) fscanf(PARAMS,"%d",&ONEROWPERIND);
  else if (!(strcmp(Word,"RECESSIVEALLELES"))) fscanf(PARAMS,"%d",&RECESSIVEALLELES);
  else if (!(strcmp(Word,"NOTAMBIGUOUS"))) fscanf(PARAMS,"%d",&NOTAMBIGUOUS);
  else if (!(strcmp(Word,"POPFLAG"))) fscanf(PARAMS,"%d",&POPFLAG);
  else if (!(strcmp(Word,"PHENOTYPE"))) fscanf(PARAMS,"%d",&PHENOTYPE);
  else if (!(strcmp(Word,"EXTRACOLS"))) fscanf(PARAMS,"%d",&EXTRACOLS);
  else if (!(strcmp(Word,"MAXPOPS"))) fscanf(PARAMS,"%d",&MAXPOPS);
  else if (!(strcmp(Word,"BURNIN"))) fscanf(PARAMS,"%d",&BURNIN);
  else if (!(strcmp(Word,"NUMREPS"))) fscanf(PARAMS,"%d",&NUMREPS);
  else if (!(strcmp(Word,"USEPOPINFO"))) fscanf(PARAMS,"%d",&USEPOPINFO);
  else if (!(strcmp(Word,"PFROMPOPFLAGONLY"))) fscanf(PARAMS,"%d",&PFROMPOPFLAGONLY);
  else if (!(strcmp(Word,"INFERALPHA"))) fscanf(PARAMS,"%d",&INFERALPHA);
  else if (!(strcmp(Word,"INFERLAMBDA"))) fscanf(PARAMS,"%d",&INFERLAMBDA);
  else if (!(strcmp(Word,"POPSPECIFICLAMBDA"))) fscanf(PARAMS,"%d",&POPSPECIFICLAMBDA);
  else if (!(strcmp(Word,"POPALPHAS"))) fscanf(PARAMS,"%d",&POPALPHAS);
  else if (!(strcmp(Word,"COMPUTEPROB"))) fscanf(PARAMS,"%d",&COMPUTEPROB);
  else if (!(strcmp(Word,"NOADMIX"))) fscanf(PARAMS,"%d",&NOADMIX);
  else if (!(strcmp(Word,"ADMBURNIN"))) fscanf(PARAMS,"%d",&ADMBURNIN);
  else if (!(strcmp(Word,"SITEBYSITE"))) fscanf(PARAMS,"%d",&SITEBYSITE);
  else if (!(strcmp(Word,"PHASED"))) fscanf(PARAMS,"%d",&PHASED);
  else if (!(strcmp(Word,"PRINTQHAT"))) fscanf(PARAMS,"%d",&PRINTQHAT);
  else if (!(strcmp(Word,"PRINTQSUM"))) fscanf(PARAMS,"%d",&PRINTQSUM);
  else if (!(strcmp(Word,"UPDATEFREQ"))) fscanf(PARAMS,"%d",&UPDATEFREQ);
  else if (!(strcmp(Word,"PRINTLIKES"))) fscanf(PARAMS,"%d",&PRINTLIKES);
  else if (!(strcmp(Word,"INTERMEDSAVE"))) fscanf(PARAMS,"%d",&INTERMEDSAVE);
  else if (!(strcmp(Word,"PRINTKLD"))) fscanf(PARAMS,"%d",&PRINTKLD);
  else if (!(strcmp(Word,"PRINTNET"))) fscanf(PARAMS,"%d",&PRINTNET);
  else if (!(strcmp(Word,"PRINTLAMBDA"))) fscanf(PARAMS,"%d",&PRINTLAMBDA);
  else if (!(strcmp(Word,"ECHODATA"))) fscanf(PARAMS,"%d",&ECHODATA);
  else if (!(strcmp(Word,"ANCESTDIST"))) fscanf(PARAMS,"%d",&ANCESTDIST);
  else if (!(strcmp(Word,"NUMBOXES"))) fscanf(PARAMS,"%d",&NUMBOXES);
  else if (!(strcmp(Word,"ANCESTPINT"))) fscanf(PARAMS,"%lf",&ANCESTPINT);
  else if (!(strcmp(Word,"GENSBACK"))) fscanf(PARAMS,"%d",&GENSBACK);
  else if (!(strcmp(Word,"MIGRPRIOR"))) fscanf(PARAMS,"%lf",&MIGRPRIOR);
  else if (!(strcmp(Word,"ALPHA"))) fscanf(PARAMS,"%lf",&ALPHA);
  else if (!(strcmp(Word,"LOG10RPROPSD"))) fscanf(PARAMS,"%lf",&LOG10RPROPSD);
  else if (!(strcmp(Word,"LOG10RMIN"))) fscanf(PARAMS,"%lf",&LOG10RMIN);
  else if (!(strcmp(Word,"LOG10RMAX"))) fscanf(PARAMS,"%lf",&LOG10RMAX);
  else if (!(strcmp(Word,"LOG10RSTART"))) fscanf(PARAMS,"%lf",&LOG10RSTART);
  else if (!(strcmp(Word,"FREQSCORR"))) fscanf(PARAMS,"%d",&FREQSCORR);
  else if (!(strcmp(Word,"ONEFST"))) fscanf(PARAMS,"%d",&ONEFST);
  else if (!(strcmp(Word,"FPRIORMEAN"))) fscanf(PARAMS,"%lf",&FPRIORMEAN);
  else if (!(strcmp(Word,"FPRIORSD"))) fscanf(PARAMS,"%lf",&FPRIORSD);
  else if (!(strcmp(Word,"LAMBDA"))) fscanf(PARAMS,"%lf",&LAMBDA);
  else if (!(strcmp(Word,"UNIFPRIORALPHA"))) fscanf(PARAMS,"%d",&UNIFPRIORALPHA);
  else if (!(strcmp(Word,"ALPHAMAX"))) fscanf(PARAMS,"%lf",&ALPHAMAX);
  else if (!(strcmp(Word,"ALPHAPRIORA"))) fscanf(PARAMS,"%lf",&ALPHAPRIORA);
  else if (!(strcmp(Word,"ALPHAPRIORB"))) fscanf(PARAMS,"%lf",&ALPHAPRIORB);
  else if (!(strcmp(Word,"ALPHAPROPSD"))) fscanf(PARAMS,"%lf",&ALPHAPROPSD);
  else if (!(strcmp(Word,"STARTATPOPINFO"))) fscanf(PARAMS,"%d",&STARTATPOPINFO);
  else if (!(strcmp(Word,"RANDOMIZE"))) fscanf(PARAMS,"%d",&RANDOMIZE);
    else if (!(strcmp(Word,"LINKAGE"))) fscanf(PARAMS,"%d",&LINKAGE);
  else if (!(strcmp(Word,"METROFREQ"))) fscanf(PARAMS,"%d",&METROFREQ);
  else if (!(strcmp(Word,"REPORTHITRATE"))) fscanf(PARAMS,"%d",&REPORTHITRATE);

  /* HIDDEN PARAMETERS */
  else if (!(strcmp(Word,"NOQS"))) fscanf(PARAMS,"%
