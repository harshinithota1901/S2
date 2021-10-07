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

int app_build_function_sets ( void )
{
     function_set fset;
     int tree_map;
     char *tree_name;
     function sets[] = 
     { { f_and, NULL, NULL, 2, "and", FUNC_DATA, -1, 0 },
       { f_or,  NULL, NULL, 2, "or",  FUNC_DATA, -1, 0 },
       { f_not, NULL, NULL, 1, "not", FUNC_DATA, -1, 0 },
       { f_if,  NULL, NULL, 3, "if",  FUNC_EXPR, -1, 0 },
       { f_a0,  NULL, NULL, 0, "a0",  TERM_NORM, -1, 0 },
       { f_a1,  NULL, NULL, 0, "a1",  TERM_NORM, -1, 0 },
       { f_a2,  NULL, NULL, 0, "a2",  TERM_NORM, -1, 0 },
       { f_d0,  NULL, NULL, 0, "d0",  TERM_NORM, -1, 0 },
       { f_d1,  NULL, NULL, 0, "d1",  TERM_NORM, -1, 0 },
       { f_d2,  NULL, NULL, 0, "d2",  TERM_NORM, -1, 0 },
       { f_d3,  NULL, NULL, 0, "d3",  TERM_NORM, -1, 0 },
       { f_d4,  NULL, NULL, 0, "d4",  TERM_NORM, -1, 0 },
       { f_d5,  NULL, NULL, 0, "d5",  TERM_NORM, -1, 0 },
       { f_d6,  NULL, NULL, 0, "d6",  TERM_NORM, -1, 0 },
       { f_d7,  NULL, NULL, 0, "d7",  TERM_NORM, -1, 0 } };

     fset.size = 15;
     fset.cset = sets;

     tree_map = 0;
     tree_name = "TREE";

     return function_sets_init ( &fset, 1, &tree_map, &tree_name, 1 );
}

void app_eval_fitness ( individual *ind )
{

     int i;
     int v, dv;

     set_current_individual ( ind );
     
     ind->r_fitness = 0.0;
     ind->hits = 0;
     
     for ( i = 0; i < 2048; ++i )
     {
          g.a = i&7;
          g.d = i>>3;

          v = evaluate_tree ( ind->tr[0].data, 0 );
          dv = !!(g.d&(1<<g.a));
          v = !!v;
          
          if ( v == dv )
          {
               ++ind->r_fitness;
               ++ind->hits;
          }
          
     }
     ind->s_fitness = 2048-ind->r_fitness;
     ind->a_fitness = ind->r_fitness/2048.0;
     ind->evald = EVAL_CACHE_VALID;
}

int app_end_of_evaluation ( int gen, multipop *mpop, int newbest,
                           popstats *gen_stats, popstats *run_stats )
{
     int i, v, dv;
     char strings[2][2] = { ".", "*" };

     oprintf ( OUT_BST, 80, "best-of-run individual's scoreboard:\n" );
     
     for ( i = 0; i < 2048; ++i )
     {
          g.a = i&7;
          g.d = i>>3;
          
          v = evaluate_tree ( run_stats[0].best[0]->ind->tr[0].data, 0 );
          dv = !!(g.d&(1<<g.a));
          v = !!v;
          
          oputs ( OUT_BST, 80, strings[dv==v] );
          if ( (i+1)%64 == 0 )
               oputs ( OUT_BST, 80, "\n" );
     }
     
     return (run_stats[0].best[0]->ind->s_fitness==0);
}

void app_end_of_breeding ( int gen, multipop *mpop )
{
     return;
}

int app_create_output_streams ( void )
{
     return 0;
}

int app_initialize ( int startfromcheckpoint )
{
     return 0;
}

void app_uninitialize ( void )
{
     return;
}

void app_write_checkpoint ( FILE *f )
{
     return;
}

void app_read_checkpoint ( FILE *f )
{
     return;
}
