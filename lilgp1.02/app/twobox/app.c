/*  lil-gp Genetic Programming System, version 1.0, 11 July 1995
 *  Copyright (C) 1995  Michigan State University
 * 
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of version 2 of the GNU General Public License as
 *  published by the Free Software Foundation.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *  
 *  Douglas Zongker       (zongker@isl.cps.msu.edu)
 *  Dr. Bill Punch        (punch@isl.cps.msu.edu)
 *
 *  Computer Science Department
 *  A-714 Wells Hall
 *  Michigan State University
 *  East Lansing, Michigan  48824
 *  USA
 *  
 */

#include <stdio.h>
#include <math.h>

#include <lilgp.h>

globaldata g;

static double app_fitness_cases[APP_FITNESS_CASES][7];

int app_build_function_sets ( void )
{
     function_set *fset;
     int *tree_map;
     char *tree_name[2];
     int ret;

     function sets[2][20] = 
     { { { f_multiply,   NULL,      NULL,         2, "*",     FUNC_DATA, -1, 0 },
         { f_protdivide, NULL,      NULL,         2, "/",     FUNC_DATA, -1, 0 },
         { f_add,        NULL,      NULL,         2, "+",     FUNC_DATA, -1, 0 },
         { f_subtract,   NULL,      NULL,         2, "-",     FUNC_DATA, -1, 0 },
         { NULL,         NULL,      NULL,        -1, "ADF0",  EVAL_DATA,  1, 0 },
         { f_l0,         NULL,      NULL,         0, "L0",    TERM_NORM, -1, 0 },
         { f_w0,         NULL,      NULL,         0, "W0",    TERM_NORM, -1, 0 },
         { f_h0,         NULL,      NULL,         0, "H0",    TERM_NORM, -1, 0 },
         { f_l1,         NULL,      NULL,         0, "L1",    TERM_NORM, -1, 0 },
         { f_w1,         NULL,      NULL,         0, "W1",    TERM_NORM, -1, 0 },
         { f_h1,         NULL,      NULL,         0, "H1",    TERM_NORM, -1, 0 } },
            
       { { f_multiply,   NULL,      NULL,         2, "*",     FUNC_DATA, -1, 0 },
         { f_protdivide, NULL,      NULL,         2, "/",     FUNC_DATA, -1, 0 },
         { f_add,        NULL,      NULL,         2, "+",     FUNC_DATA, -1, 0 },
         { f_subtract,   NULL,      NULL,         2, "-",     FUNC_DATA, -1, 0 },
         { NULL,         NULL,      NULL,         0, "ARG0",  TERM_ARG,   0, 0 },
         { NULL,         NULL,      NULL,         0, "ARG1",  TERM_ARG,   1, 0 },
         { NULL,         NULL,      NULL,         0, "ARG2",  TERM_ARG,   2, 0 } } };
         
     
     fset = (function_set *)MALLOC ( 2 * sizeof ( function_set ) );

     fset[0].size = 11;
     fset[0].cset = sets[0];

     fset[1].size = 7;
     fset[1].cset = sets[1];

     tree_map = (int *)MALLOC ( 2 * sizeof ( int ) );

     tree_map[0] = 0;
     tree_name[0] = "RPB";
     tree_map[1] = 1;
     tree_name[1] = "ADF0";

     ret = function_sets_init ( fset, 2, tree_map, tree_name, 2 );

     FREE ( fset );
     FREE ( tree_map );

     return ret;
}

void app_eval_fitness ( individual *ind )
{

     int i;
     DATATYPE value;
     double diff;

     set_current_individual ( ind );

     ind->r_fitness = 0.0;
     ind->hits = 0;

     /* loop over all the fitness cases. */
     for ( i = 0; i < APP_FITNESS_CASES; ++i )
     {
          /* fill in global structure according to current fitness case. */
          g.l0 = app_fitness_cases[i][0];
          g.w0 = app_fitness_cases[i][1];
          g.h0 = app_fitness_cases[i][2];
          g.l1 = app_fitness_cases[i][3];
          g.w1 = app_fitness_cases[i][4];
          g.h1 = app_fitness_cases[i][5];

          /* evaluate the tree. */
          value = evaluate_tree ( ind->tr[0].data, 0 );

          /* here you would score the value returned by the individual
           * and update the raw fitness and/or hits. */

          diff = fabs ( value - app_fitness_cases[i][6] );
          if ( diff < 0.01 )
               ++ind->hits;
          ind->r_fitness += diff;
          
     }

     if ( isnan(ind->r_fitness) || !finite(ind->r_fitness) )
          ind->s_fitness = MAXFLOAT;
     else
          ind->s_fitness = ind->r_fitness;
     
     ind->a_fitness = 1/(1+ind->s_fitness);

     /* always leave this line in. */
     ind->evald = EVAL_CACHE_VALID;
}

int app_end_of_evaluation ( int gen, multipop *mpop, int newbest,
                           popstats *gen_stats, popstats *run_stats )
{
     return run_stats[0].besthits == APP_FITNESS_CASES;
}

void app_end_of_breeding ( int gen, multipop *mpop )
{
     return;
}

int app_create_output_streams()
{
     return 0;
}
 

int app_initialize ( int startfromcheckpoint )
{
     int i, j;

     if ( !startfromcheckpoint )
     {
	  oprintf ( OUT_PRG, 20, "generating fitness cases:\n" );

	  for ( i = 0; i < APP_FITNESS_CASES; ++i )
	  {
	       for ( j = 0; j < 6; ++j )
		    app_fitness_cases[i][j] = (double)(random_int(10)+1);
	       app_fitness_cases[i][6] = app_fitness_cases[i][0] *
		    app_fitness_cases[i][1] * app_fitness_cases[i][2] -
			 app_fitness_cases[i][3] * app_fitness_cases[i][4] *
			      app_fitness_cases[i][5];
	  }
     }
     else
     {
	  oprintf ( OUT_PRG, 20, "fitness cases loaded from checkpoint:\n" );
     }
     
     oprintf ( OUT_PRG, 20, " L0  W0  H0    L1  W1  H1      D\n" );
     oprintf ( OUT_PRG, 20, "%3.0lf %3.0lf %3.0lf   %3.0lf %3.0lf %3.0lf  %5.0lf\n",
	      app_fitness_cases[i][0], app_fitness_cases[i][1],
	      app_fitness_cases[i][2], app_fitness_cases[i][3],
	      app_fitness_cases[i][4], app_fitness_cases[i][5],
	      app_fitness_cases[i][6] );

     return 0;
          
}

void app_uninitialize ( void )
{
     return;
}

void app_read_checkpoint ( FILE *f )
{
     int i, j;
     fscanf ( f, "%*s %d\n", &i );
     if ( i != APP_FITNESS_CASES )
          error ( E_FATAL_ERROR, "wrong # of fitness cases in checkpoint file." );
     for ( i = 0; i < APP_FITNESS_CASES; ++i )
     {
	  for ( j = 0; j < 7; ++j )
	  {
	       read_hex_block ( app_fitness_cases[i]+j, sizeof ( double ), f );
	       fscanf ( f, " %*f\n" );
	  }
     }
}

void app_write_checkpoint ( FILE *f )
{
     int i, j;
     fprintf ( f, "fitness-case-count: %d\n", APP_FITNESS_CASES );
     for ( i = 0; i < APP_FITNESS_CASES; ++i )
     {
	  for ( j = 0; j < 7; ++j )
	  {
	       write_hex_block ( app_fitness_cases[i]+j,
				sizeof ( double ), f );
	       fprintf ( f, " %f\n", app_fitness_cases[i][j] );
	  }
     }
}

    
