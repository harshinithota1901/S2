/* cgp_czj.c version 1.03 implements CGP preprocessing for lil-gp */
/* using lil-gp 1.02 */
/* extended to use weights for mutation sets 5/31/96 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lilgp.h" 

#define SMALL 0.000000001                        /* delta for double compare */
#define MINWGHT 0.00001           /* min. maintanble wght for a member of MS */

static int NumF;                                    /* num functions in fset */
static int NumT;                         /* num typeII/III functions in fset */

const int RepeatsSrc_czj=5;    /* max repeats in crossover on infeasible src */
const int RepeatsBad_czj=5;                  /* same on bad (size) offspring */
static double WghtsExt;   /* sum of wghts of cross-feasible leaves of a tree */
static double WghtsInt;                       /* same for the internal nodes */

MS_czj_t MS_czj;                      /* the global table with mutation sets */

int Function_czj; 
        /* the node to be modified has this parent; if the node  is the root 
	   then this number is the number of type I functions */

int Argument_czj; 
        /* the node to be modified is this argument of its parent
           uninitialized if the node is the root */

int MinDepth_czj;
        /* if depth=m-n is given in eithe initialization or mutation
           and the corresponding depth_abs=true, this is used to grow only
           trees at least as deep as 'm' (if possible) */

typedef int *specArray_t; /* dynamic for NumF+NumT; array[i]==1 indicates that
                             function indexed i is present in a set */
typedef struct 
{ int numF;                         /* number of F set elements in specArray */
  int numT;
  specArray_t members; 
} specArrayArg_t;
typedef specArrayArg_t *specArrays_t;     /* dynamic array of specArrayArg_t */
typedef struct 
{ int arity;
  specArrays_t Tspecs;    /* dynamic array of ptrs to specArray_t for Tspecs */
  specArrays_t Fspec;                                       /* just one here */
  specArrays_t Fspecs;
} constraint_t;                     /* NOTE: Root will use Tspecs and Fspecs */
typedef constraint_t *constraints_t; /* dynamic array for functions and Root */

static constraints_t Constraints;

static int funNumber(const char* funName, int treeNumber)
/* given funName, return its index in fset[treeNumber] or -1 if not found */
{ int i;
  for (i=0; i<fset[treeNumber].size; i++)
    if (!strcmp(funName,fset[treeNumber].cset[i].string))
      return(i);
  return(-1);
}

static void displayHeader(void)
{ printf("\n\n\t\tWELCOME TO cgp-lilgp 1.1/1.02\n");
  printf("\n\t\tdeveloped by\n");
  printf("\tCezary Z. Janikow\n");
  printf("\tUniversity of Missouri - St. Louis\n");
  printf("\temailto:janikow@arch.umsl.edu\n");
  printf("\thttp://www.cs.umsl.edu/Faculty/janikow/janikow.html\n");
  printf("\tftp://radom.umsl.edu\n");
  printf("\n\t\timplementation team:\n");
  printf("\tCezary Z. Janikow, leader\n");
  printf("\tGreg Banister, UMR student\n");
  printf("\tScott DeWeese, UMSL student\n");
  printf("\n\n\n\n\tThis is distributed as addition to lil-gp\n");
  printf("\n\tNo explicit/implicit warranty\n");
  printf("\n\n\n");
}

static void displayConstraints(int Ts, int F, int Fs)
	/* for debugging, arguments state which to display */
{ int fun, arg, entry;
  printf("\n\n\t\tCONSTRAINTS\n");
  for (fun=0; fun<NumF; fun++)
  { printf("Function \"%s\" [#%d]:\n",fset[0].cset[fun].string,fun);
    if (F)
    { printf("\tF_%s [#Fs=%d:#Ts=%d] =",fset[0].cset[fun].string,
          Constraints[fun].Fspec[0].numF,Constraints[fun].Fspec[0].numT);
      for (entry=0; entry<NumF; entry++)
        if (Constraints[fun].Fspec[0].members[entry])
           printf(" %s",fset[0].cset[entry].string);
      printf(" ||");
      for (; entry<NumF+NumT; entry++)
        if (Constraints[fun].Fspec[0].members[entry])
           printf(" %s",fset[0].cset[entry].string);
      printf("\n");
    }
    if (Fs)
      for (arg=0; arg<Constraints[fun].arity; arg++)
      { printf("\tF_%s_%d [#Fs=%d:#Ts=%d] =",fset[0].cset[fun].string,arg,
          Constraints[fun].Fspecs[arg].numF,Constraints[fun].Fspecs[arg].numT);
        for (entry=0; entry<NumF; entry++)
          if (Constraints[fun].Fspecs[arg].members[entry])
            printf(" %s",fset[0].cset[entry].string);
        printf(" ||");
        for (; entry<NumF+NumT; entry++)
          if (Constraints[fun].Fspecs[arg].members[entry])
            printf(" %s",fset[0].cset[entry].string);
        printf("\n");
      } 
    if (Ts)
      for (arg=0; arg<Constraints[fun].arity; arg++)
      { printf("\tT_%s_%d [#Fs=%d:#Ts=%d] =",fset[0].cset[fun].string,arg,
          Constraints[fun].Tspecs[arg].numF,Constraints[fun].Tspecs[arg].numT);
        for (entry=0; entry<NumF; entry++)
          if (Constraints[fun].Tspecs[arg].members[entry])
            printf(" %s",fset[0].cset[entry].string);
        printf(" ||");
        for (; entry<NumF+NumT; entry++)
          if (Constraints[fun].Tspecs[arg].members[entry])
            printf(" %s",fset[0].cset[entry].string);
        printf("\n");
      } 
  }
  printf("Root:\n",fun);
  if (Fs)
  { printf("\tF_Root [#Fs=%d:#Ts=%d] = ",
          Constraints[NumF].Fspecs[0].numF,Constraints[NumF].Fspecs[0].numT);
    for (entry=0; entry<NumF; entry++)
      if (Constraints[NumF].Fspecs[0].members[entry])
        printf(" %s",fset[0].cset[entry].string);
    printf(" ||");
    for (; entry<NumF+NumT; entry++)
      if (Constraints[NumF].Fspecs[0].members[entry])
        printf(" %s",fset[0].cset[entry].string);
    printf("\n");
  } 
  if (Ts)
  { printf("\tT_Root [#Fs=%d:#Ts=%d] = ",
          Constraints[NumF].Tspecs[0].numF,Constraints[NumF].Tspecs[0].numT);
    for (entry=0; entry<NumF; entry++)
      if (Constraints[NumF].Tspecs[0].members[entry])
        printf(" %s",fset[0].cset[entry].string);
    printf(" ||");
    for (; entry<NumF+NumT; entry++)
      if (Constraints[NumF].Tspecs[0].members[entry])
        printf(" %s",fset[0].cset[entry].string);
    printf("\n");
  } 
}

static void displayMS(void)
	/* display the mutation sets from MS_czj */
{ int i,j,k;
  printf("\n\nThe following mutation sets were computed...\n\n");
  for (i=0; i<NumF; i++)                      /* First display mutation sets */
  { printf("Function \"%s\" [#%d]: %d mutation set pairs\n",
           fset[0].cset[i].string,i,fset[0].cset[i].arity);
    for (j=0; j<fset[0].cset[i].arity; j++) 
    { printf("\tArgument %d:\n",j);
      printf("\t\tF [%d members] =",MS_czj[i][j].numF);
      for (k=0; k<MS_czj[i][j].numF; k++)
        printf(" %s",fset[0].cset[MS_czj[i][j].members[k]].string);
      printf("\n\t\tT [%d members] =",MS_czj[i][j].numT);
      for (k=0; k<MS_czj[i][j].numT; k++)
        printf(" %s",
          fset[0].cset[MS_czj[i][j].members[MS_czj[i][j].numF+k]].string);
      printf("\n");
    }
    printf("\n");
  }
  printf("Root:\n");
  printf("\t\tF [%d members] =",MS_czj[NumF][0].numF);
  for (k=0; k<MS_czj[NumF][0].numF; k++)
        printf(" %s",fset[0].cset[MS_czj[NumF][0].members[k]].string);
  printf("\n\t\tT [%d members] =",MS_czj[NumF][0].numT);
  for (k=0; k<MS_czj[NumF][0].numT; k++)
    printf(" %s",
         fset[0].cset[MS_czj[NumF][0].members[MS_czj[NumF][0].numF+k]].string);
  printf("\n\n\n");
}

static void displayWeightsWheels(int weights,int wheels)
        /* display weights/wheels from MS_czj */
        /* if weights/wheels==0 then do not display weights/wheels */
{ int i,j,k;
  if (weights==0 && wheels==0)
    return;
  printf("\n\nThese are %s%s...\n\n",weights ? "weights/" : "",
         wheels ? "wheels" : "");
  for (i=0; i<NumF; i++)              
  { printf("Function \"%s\" [#%d]: %d arguments\n",
           fset[0].cset[i].string,i,fset[0].cset[i].arity);
    for (j=0; j<fset[0].cset[i].arity; j++)
    { printf("\tArgument %d: ",j);
      printf("F [%d members, %s]  and T [%d members, %s]\n",MS_czj[i][j].numF,
             MS_czj[i][j].areFs ? "used":"not used", MS_czj[i][j].numT,
             MS_czj[i][j].areTs ? "used":"not used");
      if (weights)
      { printf("\tWeights:");
        for (k=0; k<NumF+NumT; k++)
        /* if (MS_czj[i][j].weights[k]>=MINWGHT) */
            printf(" %.3f",MS_czj[i][j].weights[k]);
        printf("\n");
      }
      if (wheels)
      { printf("\tWheels:");
        for (k=0; k<MS_czj[i][j].numFT; k++)
          printf(" %.3f",MS_czj[i][j].wheel[k]);
        printf("\n");
      }
    }
    printf("\n");
  }
  printf("Root: ");
  printf("F [%d members, %s] and T [%d members, %s]\n",MS_czj[NumF][0].numF,
         MS_czj[NumF][0].areFs ? "used":"not used", MS_czj[NumF][0].numT,
         MS_czj[NumF][0].areTs ? "used":"not used");
  if (weights)
  { printf("\tWeights:");
    for (k=0; k<NumF+NumT; k++)
      /* if (MS_czj[NumF][0].weights[k]>MINWGHT) */
        printf(" %.3f",MS_czj[NumF][0].weights[k]);
    printf("\n");
  }
  if (wheels)
  { printf("\tWheels:");
    for (k=0; k<MS_czj[NumF][0].numFT; k++)
      printf(" %.3f",MS_czj[NumF][0].wheel[k]);
    printf("\n");
  }
  printf("\n");
}

static double readMinWghtto1(const char *prompt)
/* read (0,1], and if entry < MINWGHT then set it as MINWGHT */
{ double what;
  int res;
  printf("%s: ",prompt);
  res=scanf("%lf",&what);
  while (res<1 || what >1 || what <0)
  { if (res==EOF)
      exit(1);
    fflush(stdin);
    printf("\tInvalid weight: %s: ",prompt);
    res=scanf("%lf",&what);
  }
  if (what<MINWGHT)
    what=MINWGHT;                           /* smaller values become minimal */
  return(what);
}

static void readWeightsSetWheels(void)
/* read weights for mutation set entries and construct wheels */
/* assume weights/wheels are set for all members equalweighted */
/* assume weights for non-members are set to -1 */
/* areFs and areTs members of MS_czj are set to true if ate least one $/
/*   members has weight >MINWGHT */
{ int i,j,k;
  double adjWght;
  int areFs, areTs;
  printf("\n\nSetting initial weights for mutation set members...\n\n");
  printf("\nInitial weights are all equal. Do you accept [0 to change]: ");
  scanf("%d",&i);
  if (i) 
    return;                               /* leave inital weights and wheels */
  for (i=0; i<NumF; i++)                   
  { printf("\n");
    printf("Function \"%s\" [#%d]: %d mutation set pairs\n",
           fset[0].cset[i].string,i,fset[0].cset[i].arity);
    for (j=0; j<fset[0].cset[i].arity; j++)
    { areFs=0; areTs=0; 
      printf("Argument %d:\n",j);
      printf("\tF [%d members] =",MS_czj[i][j].numF);
      for (k=0; k<MS_czj[i][j].numF; k++)
        printf(" %s",fset[0].cset[MS_czj[i][j].members[k]].string);
      printf("\n\tT [%d members] =",MS_czj[i][j].numT);
      for (k=0; k<MS_czj[i][j].numT; k++)
        printf(" %s",
               fset[0].cset[MS_czj[i][j].members[MS_czj[i][j].numF+k]].string);
      printf("\n\n\tReading the weights for type I functions...\n");
      for (k=0; k<MS_czj[i][j].numF; k++)
      { printf("\tFunction \"%s\" [%d]: ",
               fset[0].cset[MS_czj[i][j].members[k]].string,
               MS_czj[i][j].members[k]);
        MS_czj[i][j].weights[MS_czj[i][j].members[k]]=
                                            readMinWghtto1("give weight (0,1]");
      }
      printf("\n\tReading the weights for type II/III terminals...\n");
      for (k=0; k<MS_czj[i][j].numT; k++)
      { printf("\tTerminal \"%s\" [%d]: ",
               fset[0].cset[MS_czj[i][j].members[MS_czj[i][j].numF+k]].string,
               MS_czj[i][j].members[MS_czj[i][j].numF+k]);
        MS_czj[i][j].weights[MS_czj[i][j].members[MS_czj[i][j].numF+k]]=
                                            readMinWghtto1("give weight (0,1]");
      }
    /* now all non-memb weights are -1, all memb weights are in [MINWGHT..1] */
                   /* now set mut wheel skipping over weights <MINWGHT+SMALL */
      for (k=0; k<MS_czj[i][j].numF; k++) 
      { if (MS_czj[i][j].weights[MS_czj[i][j].members[k]]<MINWGHT+SMALL)
          adjWght=0;
        else
        { adjWght=MS_czj[i][j].weights[MS_czj[i][j].members[k]];
          areFs=1;
        }
        MS_czj[i][j].wheel[k]= (k==0) ? adjWght:MS_czj[i][j].wheel[k-1]+adjWght;
      }
      for (k=MS_czj[i][j].numF; k<MS_czj[i][j].numFT; k++)
      { if (MS_czj[i][j].weights[MS_czj[i][j].members[k]]<MINWGHT+SMALL)
          adjWght=0;
        else
        { adjWght=MS_czj[i][j].weights[MS_czj[i][j].members[k]];
          areTs=1;
        }
        MS_czj[i][j].wheel[k]= (k==0) ? adjWght:MS_czj[i][j].wheel[k-1]+adjWght;
      }
      MS_czj[i][j].areFs=areFs;
      MS_czj[i][j].areTs=areTs;
      if (!areFs && !areTs)
      { fprintf(stderr,
                "\tno member of f=%d arg=%d has any weight >MINWGHT\n",i,j);
        exit(1);
      }
    }
    printf("\n\n");
  }
  printf("Root:\n");
  areFs=0; areTs=0;
  printf("\t\tF [%d members] =",MS_czj[NumF][0].numF);
  for (k=0; k<MS_czj[NumF][0].numF; k++)
        printf(" %s",fset[0].cset[MS_czj[NumF][0].members[k]].string);
  printf("\n\t\tT [%d members] =",MS_czj[NumF][0].numT);
  for (k=0; k<MS_czj[NumF][0].numT; k++)
    printf(" %s",
         fset[0].cset[MS_czj[NumF][0].members[MS_czj[NumF][0].numF+k]].string);
  printf("\n\tReading the weights for type I functions...\n");
  for (k=0; k<MS_czj[NumF][0].numF; k++)
  { printf("\tFunction \"%s\" [%d]: ",
           fset[0].cset[MS_czj[NumF][0].members[k]].string,
           MS_czj[NumF][0].members[k]);
    MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]]=
                                            readMinWghtto1("give weight (0,1]");
  }
  printf("\n\tReading the weights for type II/III terminals...\n");
  for (k=0; k<MS_czj[NumF][0].numT; k++)
  { printf("\tTerminal \"%s\" [%d]: ",
           fset[0].cset[MS_czj[NumF][0].members[MS_czj[NumF][0].numF+k]].string,
           MS_czj[NumF][0].members[MS_czj[NumF][0].numF+k]);
    MS_czj[NumF][0].weights[MS_czj[NumF][0].members[MS_czj[NumF][0].numF+k]]=
                                        readMinWghtto1("give weight (0,1]");
  }
  for (k=0; k<MS_czj[NumF][0].numF; k++)
  { if (MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]]<MINWGHT+SMALL)
      adjWght=0;
    else
    { adjWght=MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]];
      areFs=1;
    }
    MS_czj[NumF][0].wheel[k]= (k==0) ? adjWght : 
                                       MS_czj[NumF][0].wheel[k-1]+adjWght;
  }
  for (k=MS_czj[NumF][0].numF; k<MS_czj[NumF][0].numFT; k++)
  { if (MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]]<MINWGHT+SMALL)
      adjWght=0;
    else
    { adjWght=MS_czj[NumF][0].weights[MS_czj[NumF][0].members[k]];
      areTs=1;
    }
    MS_czj[NumF][0].wheel[k]= (k==0) ? adjWght : 
                                       MS_czj[NumF][0].wheel[k-1]+adjWght;
  }
  MS_czj[NumF][0].areFs=areFs;
  MS_czj[NumF][0].areTs=areTs;
  if (!areFs && !areTs)
  { fprintf(stderr,"\tno member of Root sets has any weight >MINWGHT\n");
    exit(1);
  }
  printf("\n\n");
}

static void displayFT(void)
{ int i;
  for (i=0; i<25; i++)
    printf("\n");
  printf("%d ordinary functions: \n",NumF);
  for (i=0; i<NumF; i++)
    printf("%5d = %s\n",i,fset[0].cset[i].string); 
  printf("%d terminals (ordinary and ephemeral): \n",NumT);
  for (; i<NumF+NumT; i++)
    printf("%5d = %s\n",i,fset[0].cset[i].string); 
  printf("Separate entries by [ ,;]  Hit <ENTER> for empty set\n");
  printf("Use either function names or numbers, in any order\n\n");
}

static void readOneConstraintSet(const char*prompt, specArrays_t setP, int max)
	/* read one set from one line; max is the highest index allowed */
{ int entry, status;
  char buf[BUFSIZ];
  char sep[]=" ,;\n";
  char *p;
  for (entry=0; entry<NumF+NumT; entry++)
    setP->members[entry]=0;                            /* reset set to empty */
  setP->numF=setP->numT=0;                          /* reset member counters */
  printf("%s [0..%d] = ",prompt,max);
  if (fgets(buf,BUFSIZ,stdin)==NULL)
  { fprintf(stderr,"ERROR: failed reading constrained\n");
    exit(1);
  }
  p=strtok(buf,sep);
  while (p!=NULL)
  { if ((entry=funNumber(p,0))>=0 || (sscanf(p,"%d",&entry)>0))
    { if (entry<0 || entry >max)
        printf("\a\t\t%d entry out of range\n",entry);
      else                                /* entry is a valid function index */
      { setP->members[entry]=1;
        if (entry<NumF) 
          setP->numF++;
        else 
          setP->numT++;
      }
    }
    else               /* failed reading an integer or invalid function name */
      printf("\t\ainvalid entry\n");
    p=strtok(NULL,sep);
  }
}

static void readFTspecs(void)
{ int i, j;
  char prompt[BUFSIZ];
  
  Constraints=(constraints_t)calloc((size_t)(NumF+1),sizeof(constraint_t));
                                                  /* last entry for the Root */
  for (i=0; i<NumF; i++)                   /* first work on type I functions */
  { displayFT();
    printf("Function %d=%s:\n",i,fset[0].cset[i].string);
    Constraints[i].arity=fset[0].cset[i].arity;
    Constraints[i].Fspec=(specArrays_t)calloc((size_t)1,sizeof(specArrayArg_t));
    Constraints[i].Fspec[0].members=(specArray_t)calloc((size_t)(NumF+NumT),
                                                 sizeof(int));
    sprintf(prompt,"\tF_%s (exclusions)",fset[0].cset[i].string);
    readOneConstraintSet(prompt,&(Constraints[i].Fspec[0]),NumF-1);
                                            /* Note: type I only for callers */
    Constraints[i].Tspecs=(specArrays_t)calloc((size_t)Constraints[i].arity,
                                               sizeof(specArrayArg_t));
    Constraints[i].Fspecs=(specArrays_t)calloc((size_t)Constraints[i].arity,
                                               sizeof(specArrayArg_t));
    for (j=0; j<Constraints[i].arity; j++)
    { Constraints[i].Tspecs[j].members=(specArray_t)calloc((size_t)(NumF+NumT),
                                                 sizeof(int));
      Constraints[i].Fspecs[j].members=(specArray_t)calloc((size_t)(NumF+NumT),
                                                 sizeof(int));
      sprintf(prompt,"\tF_%s_%d (exclusions)",fset[0].cset[i].string,j);
      readOneConstraintSet(prompt,&(Constraints[i].Fspecs[j]),NumF+NumT-1);
      sprintf(prompt,"\tT_%s_%d (inclusions)",fset[0].cset[i].string,j);
      readOneConstraintSet(prompt,&(Constraints[i].Tspecs[j]),NumF+NumT-1);
    }
  }
  Constraints[NumF].arity=1;
  Constraints[NumF].Fspec=NULL;
  Constraints[i].Tspecs=(specArrays_t)calloc((size_t)1,sizeof(specArrayArg_t));
  Constraints[i].Fspecs=(specArrays_t)calloc((size_t)1,sizeof(specArrayArg_t));
  Constraints[i].Tspecs[0].members=(specArray_t)calloc((size_t)(NumF+NumT),
                                                 sizeof(int));
  Constraints[i].Fspecs[0].members=(specArray_t)calloc((size_t)(NumF+NumT),
                                                 sizeof(int));
  displayFT();
  printf("Root:");
  readOneConstraintSet("\tF^Root (exclusions)",&(Constraints[NumF].Fspecs[0]),
                       NumF+NumT-1);
  readOneConstraintSet("\tT^Root (inclusions)",&(Constraints[NumF].Tspecs[0]),
                       NumF+NumT-1);
}

static void generateNF(void)
	/* from specs of Constraints generate NF in Fspecs of Constraints */
	/* this involves creating T-extensive Fspecs */
{ int fun, arg, entry;
  for (fun=0; fun<NumF; fun++)           /* first create T-extensive F-specs */
    for (arg=0; arg<Constraints[fun].arity; arg++)
      for (entry=0; entry<NumF+NumT; entry++)
        if (Constraints[fun].Tspecs[arg].members[entry]==0)
          Constraints[fun].Fspecs[arg].members[entry]=1;
  for (entry=0; entry<NumF+NumT; entry++)               /* same for the Root */
    if (Constraints[NumF].Tspecs[0].members[entry]==0)
      Constraints[NumF].Fspecs[0].members[entry]=1;

  for (fun=0; fun<NumF; fun++)              /* now create F-extensive Fspecs */
    for (entry=0; entry<NumF; entry++)
      if (Constraints[fun].Fspec[0].members[entry]!=0)     /* must extend it */
        for (arg=0; arg<Constraints[entry].arity; arg++)
          Constraints[entry].Fspecs[arg].members[fun]=1;

  for (fun=0; fun<NumF+1; fun++)            /* recount set entries in Fspecs */
    for (arg=0; arg<Constraints[fun].arity; arg++)
    { Constraints[fun].Fspecs[arg].numF=0;
      Constraints[fun].Fspecs[arg].numT=0;
      for (entry=0; entry<NumF; entry++)
        if (Constraints[fun].Fspecs[arg].members[entry]!=0)
          Constraints[fun].Fspecs[arg].numF++;
      for (; entry<NumF+NumT; entry++)
        if (Constraints[fun].Fspecs[arg].members[entry]!=0)
          Constraints[fun].Fspecs[arg].numT++;
    }
}

static void generateMS(void)
	/* create MS from the Fspecs part (ie F-intensive) of Constraints */
{ int fun, arg, entry, k;
  MS_czj=(MS_czj_t)calloc((size_t)(NumF+1),sizeof(mutationSets_czj_t));        
                                            /* one (last) entry for the Root */

  for (fun=0; fun<NumF; fun++)                   /* set all type I functions */
  { MS_czj[fun]=(mutationSets_czj_t)calloc((size_t)fset[0].cset[fun].arity,
                                           sizeof(mutationSet_czj_t));
    for (arg=0; arg<fset[0].cset[fun].arity; arg++)
    { MS_czj[fun][arg].numF=NumF-Constraints[fun].Fspecs[arg].numF;
      MS_czj[fun][arg].numT=NumT-Constraints[fun].Fspecs[arg].numT;
      MS_czj[fun][arg].numFT=MS_czj[fun][arg].numF+MS_czj[fun][arg].numT;
      if (MS_czj[fun][arg].numFT==0)
        fprintf(stderr,"\n\tBoth sets empty for function %d=%s argument %d\n\a",
                fun,fset[0].cset[fun].string,arg);
      MS_czj[fun][arg].members=
         (int*)calloc((size_t)(MS_czj[fun][arg].numFT),sizeof(int));
      MS_czj[fun][arg].weights=
         (double*)calloc((size_t)(NumF+NumT),sizeof(double));
      MS_czj[fun][arg].wheel=
         (double*)calloc((size_t)(MS_czj[fun][arg].numFT),sizeof(double));
      for (entry=0,k=0; k<NumF+NumT; k++)
        if (Constraints[fun].Fspecs[arg].members[k]==0)
        { MS_czj[fun][arg].members[entry]=k;
          MS_czj[fun][arg].weights[k]=1.0;
          MS_czj[fun][arg].wheel[entry]= (entry==0) ? 
            MS_czj[fun][arg].weights[k] : 
            MS_czj[fun][arg].wheel[entry-1]+MS_czj[fun][arg].weights[k];
          entry++;
        }
        else
          MS_czj[fun][arg].weights[k]= -1.0;
      MS_czj[fun][arg].areFs= !!MS_czj[fun][arg].numF;
      MS_czj[fun][arg].areTs= !!MS_czj[fun][arg].numT;
    }
  }
  MS_czj[NumF]=(mutationSets_czj_t)calloc((size_t)1,
                    sizeof(mutationSet_czj_t));              /* for the Root */
  MS_czj[NumF][0].numF=NumF-Constraints[NumF].Fspecs[0].numF;
  MS_czj[NumF][0].numT=NumT-Constraints[NumF].Fspecs[0].numT;
  MS_czj[NumF][0].numFT=MS_czj[NumF][0].numF+MS_czj[NumF][0].numT;
  if (MS_czj[NumF][0].numFT==0)
  { printf("\n\tBoth Root sets empty - no valid programs exist\n\a");
    exit(1);
  }
  MS_czj[NumF][0].members=
      (int*)calloc((size_t)(MS_czj[NumF][0].numFT),sizeof(int));
  MS_czj[NumF][0].weights=
      (double*)calloc((size_t)(NumF+NumT),sizeof(double));
  MS_czj[NumF][0].wheel=
      (double*)calloc((size_t)(MS_czj[NumF][0].numFT),sizeof(double));
  for (entry=0,k=0; k<NumF+NumT; k++)
    if (Constraints[NumF].Fspecs[0].members[k]==0)
    { MS_czj[NumF][0].members[entry]=k;
      MS_czj[NumF][0].weights[k]=1.0;
      MS_czj[NumF][0].wheel[entry]=(entry==0) ?
            MS_czj[NumF][0].weights[k] :
            MS_czj[NumF][0].wheel[entry-1]+MS_czj[NumF][0].weights[k];
      entry++;
    }
    else
      MS_czj[NumF][0].weights[k]= -1.0;
  MS_czj[NumF][0].areFs= !!MS_czj[NumF][0].numF;
  MS_czj[NumF][0].areTs= !!MS_czj[NumF][0].numT;
}

void create_MS_czj(void)
        /* will access global fset function table, and will allocate and
           initialize MS_czj; check that no pair is completely empty */
	/* Must be called after fset is created and ordered,
           but before initializing the population */
{ int what=0;
  NumF=fset[0].function_count;                                      /* |F_I| */
  NumT=fset[0].terminal_count;                           /* |F_II| + |F_III| */
  displayHeader();
  readFTspecs();
/*  printf("\nRead the following constraints...\n");
    displayConstraints(1,1,1); 
*/
  generateNF(); 
  printf("\nThe normal constraints are...\n");
  displayConstraints(0,0,1);
  generateMS();
  displayMS();
  readWeightsSetWheels();
  printf("\nWheels are...\n"); 
  displayWeightsWheels(1,1); 
  printf("Send 1 to continue, anything else to quit cgp-lil-gp: ");
  scanf("%d",&what);
  if (what!=1)
    exit(1);
  printf("\n\n");
}

static int spinWheel(int startI, int stopI, double *wheel)
/* spin the wheel returning an index between startI and stopI inclusively,
   with probability proportional to wheel allocation (roulette) */
{ double mark,begining;
  begining=startI ? wheel[startI-1] : 0;
  mark=begining+random_double()*(wheel[stopI]-begining);
  while (mark > wheel[startI]) startI++;
  return(startI);
}

int random_F_czj(void)
/* return a random type I index from fset, but which appear in the */
/*    mutation set for Function_czj/Argument_czj */
/* if the set is empty, call random_FT_czj() */
/* NOTE: set is considered empty if numF==0 or each weight is <MINWGHT+SMALL */
{ int randIndex;
  randIndex=MS_czj[Function_czj][Argument_czj].numF;
  if (randIndex==0 || MS_czj[Function_czj][Argument_czj].areFs==0)
    return(random_FT_czj());    
  randIndex=spinWheel(0,randIndex-1,MS_czj[Function_czj][Argument_czj].wheel);
  return MS_czj[Function_czj][Argument_czj].members[randIndex];
}
 
int random_T_czj(void)
/* as random_F_czj, except that extract members of T */
{ int randIndex;
  if (MS_czj[Function_czj][Argument_czj].numT==0 || 
      MS_czj[Function_czj][Argument_czj].areTs==0)
    return(random_FT_czj());           
  randIndex=spinWheel(MS_czj[Function_czj][Argument_czj].numF,
                      MS_czj[Function_czj][Argument_czj].numFT-1,
                      MS_czj[Function_czj][Argument_czj].wheel); 
  return MS_czj[Function_czj][Argument_czj].members[randIndex]; 
}

int random_FT_czj(void)
        /* return a random type I/II/III index from fset, but which appear in
           the mutation set for Function_czj/Argument_czj */
        /* assume both sets (F and T) are not empty */
{ int randIndex;
  if (MS_czj[Function_czj][Argument_czj].numFT==0)
  { fprintf(stderr,"\nERROR: both sets should not be empty\n");
    exit(1);
  }
  if (MS_czj[Function_czj][Argument_czj].areFs==0 &&
      MS_czj[Function_czj][Argument_czj].areTs==0)
  { fprintf(stderr,"\nERROR: both sets should not have empty weights\n");
    exit(1);
  }
  randIndex=spinWheel(0,MS_czj[Function_czj][Argument_czj].numFT-1,
                      MS_czj[Function_czj][Argument_czj].wheel);
  return MS_czj[Function_czj][Argument_czj].members[randIndex];
} 

static int markXNodes_recurse_czj( lnode **t )
/* assume Function_czj and Argument_czj are set to indicate dest point */
/* mark and count all feasible source nodes for crossover in tree */
/* marking is done with the corresponding weights w/r to dest parent */
/*   and wheel values are accumulated */
/* clear all other marking weights and wheel values to 0 */
/* sum the weights of feasible internal/external nodes in WghtsInt/WghtsWxt */
/* return the number of marked nodes */
/* NOTE: wheel entry may be the same as that of the last node if this node */
/*   is infeasible -> this will ensure that this node is not selected later */
{ function *f = (**t).f;
  double *wheelExt_czj=&((**t).wheelExt_czj);
  double *wheelInt_czj=&((**t).wheelInt_czj);
  int j;
  double wght=0;
  int total;

  ++*t;                              /* step the pointer over the function. */

  if ( f->arity == 0 )                                    /* it is external */
  { if (f->ephem_gen)
               ++*t;                                  /* etra value to skip */
    wght=MS_czj[Function_czj][Argument_czj].weights[f->index];
    if (wght<(MINWGHT+SMALL))     /* not in this mutation set or do not use */
      total=0;
    else
    { WghtsExt+=wght;
      total=1;
    }
    *wheelInt_czj=WghtsInt;
    *wheelExt_czj=WghtsExt;
    return total;
  }
  switch (f->type)                           /* here only for internal nodes */
  { case FUNC_DATA: case EVAL_DATA:
      wght=MS_czj[Function_czj][Argument_czj].weights[f->index];
      if (wght<(MINWGHT+SMALL))    /* not in this mutation set or do not use */
        total=0;
      else
      { WghtsInt+=wght;
        total=1;
      }
      *wheelInt_czj=WghtsInt;
      *wheelExt_czj=WghtsExt;
      for (j=0; j<f->arity; ++j)
        total+=markXNodes_recurse_czj(t);    /* t has already been advanced */
      break;
    case FUNC_EXPR: case EVAL_EXPR:
      wght=MS_czj[Function_czj][Argument_czj].weights[f->index];
      if (wght<(MINWGHT+SMALL))    /* not in this mutation set or do not use */
        total=0;
      else
      { WghtsInt+=wght;
        total=1;
      }
      *wheelInt_czj=WghtsInt;
      *wheelExt_czj=WghtsExt;
      for (j=0; j<f->arity; ++j)
      { ++*t;                   /* skip the pointer over the skipsize node. */
        total+=markXNodes_recurse_czj(t);
      }
      break;
  } /* switch */
  return total;
}

int markXNodes_czj( lnode *data )
/* assume Function_czj and Argument_czj are set to indicate dest point */
/* mark all nodes in tree which are feasible sources with their wghts */
/*   while contructing the wheels for internal and external nodes */
/* accummulate total wght marked in WghtsInt and WghtsExt */
/*   for the internal and the external nodes, respectively */
/* return the number of marked nodes */
{ lnode *t=data;
  WghtsInt=0;
  WghtsExt=0;
  return markXNodes_recurse_czj (&t);
}

static lnode *getSubtreeMarked_recurse_czj(lnode **t, double mark)
/* assume feasible internal and external nodes are marked with wheel values */
/*   and 'mark' determines which wheel entry is used */
/* this function spins the wheel looking for any node */
{ function *f = (**t).f;
  double *wheelExt_czj=&((**t).wheelExt_czj);
  double *wheelInt_czj=&((**t).wheelInt_czj);
  lnode *r;
  int i;

  if (mark < (*wheelExt_czj + *wheelInt_czj))      
    return *t;                 /* this might be either internal or external */
  ++*t;                                              /* move t to next node */
  if (f->arity==0)
  { if (f->ephem_gen)
      ++*t;                                /* skip over the terminal nodes. */
    return NULL;                                          /* not found here */
  }
  switch (f->type)
  { case FUNC_DATA: case EVAL_DATA:
      for (i=0; i<f->arity; i++)
      { r=getSubtreeMarked_recurse_czj(t,mark);
        if (r!=NULL)
          return r;                                      /* subtree found */
      }
      break;
    case FUNC_EXPR: case EVAL_EXPR:
      for (i=0; i<f->arity; i++)
      { ++*t;
        r=getSubtreeMarked_recurse_czj(t,mark);
        if (r!=NULL)
          return r;
      }
      break;
  }
  return NULL;
}

static lnode *getSubtreeMarkedInt_recurse_czj(lnode **t, double mark)
/* same as getSubtreeMarked_recurse_czj except look only internal nodes */
{ function *f = (**t).f;
  double *wheelInt_czj=&((**t).wheelInt_czj);
  lnode *r;
  int i;

  if (f->arity==0)                               /* it is external, skip it */
  { ++*t; 
    if (f->ephem_gen)
      ++*t; 
    return NULL;
  }
  if (mark < *wheelInt_czj)                           /* return this node */
      return *t;
  ++*t;                                              /* move t to next node */
  switch (f->type)
  { case FUNC_DATA: case EVAL_DATA:
      for (i=0; i<f->arity; i++)
      { r=getSubtreeMarkedInt_recurse_czj(t,mark);
        if (r!=NULL)
          return r;                                      /* subtree found */
      }
      break;
    case FUNC_EXPR: case EVAL_EXPR:
      for (i=0; i<f->arity; i++)
      { ++*t;
        r=getSubtreeMarkedInt_recurse_czj(t,mark);
        if (r!=NULL)
          return r;
      }
      break;
  }
  return NULL;
}

static lnode *getSubtreeMarkedExt_recurse_czj(lnode **t, double mark)
/* same as getSubtreeMarked_recurse_czj except look only external nodes */
{ function *f = (**t).f;
  double *wheelExt_czj=&((**t).wheelExt_czj);
  lnode *r;
  int i;

  if (f->arity==0)                           /* it is external, check it out */
  { if (mark<*wheelExt_czj)                            /* return this node */
        return *t;
    ++*t; 
    if (f->ephem_gen)
      ++*t;
    return NULL;
  }
  ++*t;                        /* if here than it is an internal node - skip */
  switch (f->type)
  { case FUNC_DATA: case EVAL_DATA:
      for (i=0; i<f->arity; i++)
      { r=getSubtreeMarkedExt_recurse_czj(t,mark);
        if (r!=NULL)
          return r;                                         /* subtree found */
      }
      break;
    case FUNC_EXPR: case EVAL_EXPR:
      for (i=0; i<f->arity; i++)
      { ++*t;
        r=getSubtreeMarkedExt_recurse_czj(t,mark);
        if (r!=NULL)
          return r;
      }
      break;
  }
  return NULL;
}

lnode *getSubtreeMarked_czj(lnode *data, int intExt)
/* assume tree is filled with both internal and external wheels */
/*   accummulated in WghtsInt and WghtsExt and at least one node is feasible */
/* return a node with selection prob. proportional to its weight */
/* if no nodes found return NULL */
/* if intExt==0 then looking for both internal and external */
/* if intExt==1 then looking for internal, and switch to external only if */
/*   no internal marked */
/* the opposite for intExt==2 */
{ lnode *el = data;
  if (intExt==0 || intExt==1 && WghtsInt<SMALL ||
      intExt==2 && WghtsExt<SMALL)
    return              /* override 'intExt' parameter and look for any node */
      getSubtreeMarked_recurse_czj(&el,(WghtsInt+WghtsExt)*random_double());
  if (intExt==1)
    return
      getSubtreeMarkedInt_recurse_czj(&el,WghtsInt*random_double());
  return 
    getSubtreeMarkedExt_recurse_czj(&el,WghtsExt*random_double());
}

static int verify_tree_czj_recurse ( lnode **t )
/* return #times the tree pointed by t violates MS_czj */
/*   note: *t always points at a function node here: save the function */
{ function *f = (**t).f;
  int i;
  int total=0;

  if (MS_czj[Function_czj][Argument_czj].weights[f->index]<0)
    total++;                                                     /* invalid */
  ++*t;                              /* step the pointer over the function. */
  if (f->arity==0)
  { if (f->ephem_gen)
      ++*t;    /* skip the pointer over the ERC value if this node has one. */
    return total; 
  }
  switch (f->type)
  { case FUNC_DATA: case EVAL_DATA:
      for (i=0; i<f->arity; ++i)
      { Function_czj = f->index;
        Argument_czj = i;
        total+=verify_tree_czj_recurse (t);
      }
      break;
    case FUNC_EXPR: case EVAL_EXPR:
      for (i=0; i<f->arity; ++i)
      { ++*t;                 /* skip the pointer over the skipsize node. */
        Function_czj = f->index;
        Argument_czj = i;
        total+=verify_tree_czj_recurse (t);
      }
      break;
  } /* switch */
  return total;
}

int verify_tree_czj ( lnode *tree )
/* return #times the tree pointed by tree violates MS_czj */
{ lnode *t = tree;
  int f_save=Function_czj;
  int a_save=Argument_czj;
  int total;

  Function_czj = fset[0].function_count;
  Argument_czj = 0;                                  /* start from the Root */
  total=verify_tree_czj_recurse (&t);
  Function_czj=f_save;
  Argument_czj=a_save;
  return total;
}
