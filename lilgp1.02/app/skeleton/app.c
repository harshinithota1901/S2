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

/* declare the (globaldata) struct if you use it.  the name "g" is suggested
 * for this instance.
 */

globaldata g;

/* app_build_function_sets()
 *
 * this function should build data structures describing the function
 * set(s) and the trees of the individual.  it should then pass these
 * structures to function_sets_init(), destroy its data structures,
 * and return the (integer) return value from function_sets_init().
 *
 * see the manual for more information on constructing this function.
 */

int app_build_function_sets ( void )
{
     function_set *fset;
     int *tree_map;
     char **tree_name;
     int ret;

     int num_function_sets = 1;
     int num_trees = 1;

     fset = (function_set *)MALLOC ( num_function_sets *
                                    sizeof ( function_set ) );
     fset[0].size = /* # of functions and terminals */;
     fset[0].cset = /* definition of the function set */;

     tree_map = (int *)MALLOC ( num_trees * sizeof ( int ) );
     /* specify, for each tree, what function set to use. */
     tree_map[0] = 0;

     tree_name = (char **)MALLOC ( num_trees * sizeof ( char * ) );
     /* give each tree a name. */
     tree_name[0] = "TREE";

     ret = function_sets_init ( fset, num_function_sets, tree_map,
                               tree_name, num_trees );

     FREE ( fset );
     FREE ( tree_map );
     FREE ( tree_name );

     return ret;
}

/* app_eval_fitness()
 *
 * this function should evaluate the fitness of the individual.  typically
 * this function will loop over all the fitness cases.  the following
 * fields in the (individual *) should be filled in:
 *    r_fitness    (raw fitness)
 *    s_fitness    (standardized fitness)
 *    a_fitness    (adjusted fitness)
 *    hits         (hits)
 *    evald        (always set to EVAL_CACHE_VALID)
 */

void app_eval_fitness ( individual *ind )
{

     int i;
     DATATYPE value;

     set_current_individual ( ind );
     
     ind->r_fitness = 0.0;
     ind->hits = 0;

     /* loop over all the fitness cases. */
     for ( i = 0; i < 10; ++i )
     {
          /* fill in global structure according to current fitness case. */
          g.x = /* something */;

          /* evaluate a tree of the individual.  this example uses tree 0,
           but you can evaluate any tree as many times as you want. */
          value = evaluate_tree ( ind->tr[0].data, 0 );

          /* here you would score the value returned by the individual
           * and update the raw fitness and/or hits. */
        
     }

     /* compute the standardized and raw fitness. */
     ind->s_fitness = ind->r_fitness;
     ind->a_fitness = 1/(1+ind->s_fitness);

     /* always leave this line in. */
     ind->evald = EVAL_CACHE_VALID;
}

/* app_end_of_evaluation()
 *
 * this is called every generation after the entire popualation has
 * been evaluated.  "mpop" is a pointer to basically everything --
 * every individual in every subpopulation, exchange topology for
 * multipop runs, etc.  see the documentation for useful things to
 * find inside it.
 *
 * run_stats[0].best[0]->ind is always a pointer to the best-of-run
 * individual.  if the "newbest" flag is 1, it indicates that this
 * is a NEW best-of-run individual.
 *
 * similarly, gen_stats[0].best[0]->ind points to the best-of-generation
 * individual.
 *
 * DO NOT MODIFY ANY OF THE STRUCTURES YOU ARE PASSED POINTERS TO.
 *
 * should return 1 to indicate that the termination criterion has been
 * met and the run should terminate, 0 otherwise.  you do not need to
 * test the generation number against the maximum -- this is done by
 * the kernel.
 *
 * you can dynamically modify breeding parameters and/or subpop exchange
 * topology by making the appropriate changes to the parameter database.
 * see the manual for more information.  after making changes, you need
 * to call:
 *    rebuild_breeding_table ( mpop );
 * or
 *    rebuild_exchange_topology ( mpop );
 * for your changes to take effect.
 */
 
int app_end_of_evaluation ( int gen, multipop *mpop, int newbest,
                           popstats *gen_stats, popstats *run_stats )
{
     return 0;
}

/* app_end_of_breeding()
 *
 * this is called every generation after the next population is created
 * (but before it is evaluated).
 *
 * DO NOT MODIFY ANY OF THE STRUCTURES YOU ARE PASSED POINTERS TO.
 */

void app_end_of_breeding ( int gen, multipop *mpop )
{
     return;
}

/* app_create_output_streams()
 *
 * if you are going to create any custom output streams, do it here.
 * this is not documented yet, look at the source for the symbolic
 * regression problem for a guide.
 *
 * return 0 if init went OK and the run can proceed, or 1 to abort.
 */

int app_create_output_streams()
{
     return 0;
}
 

/* app_initialize()
 *
 * this should perform any application-specific initialization the user
 * needs to do.
 *
 * return 0 if init went OK and the run can proceed, or 1 to abort.
 */

int app_initialize ( int startfromcheckpoint )
{
     return 0;
}

/* app_uninitialize()
 *
 * perform application cleanup (free memory, etc.)
 */

void app_uninitialize ( void )
{
     return;
}

/* app_read_checkpoint()
 *
 * read state information from a checkpoint file.  it is passed a handle
 * for a file that has been opened for read in text mode.  it should
 * leave the file pointer at the end of the user information.
 */

void app_read_checkpoint ( FILE *f )
{
     return;
}

/* app_write_checkpoint()
 *
 * write state information to a checkpoint file.  it is passed a handle
 * for a file that has been opened for write in text mode.  it should
 * leave the file pointer at the end of the user information.
 */

void app_write_checkpoint ( FILE *f )
{
     return;
}

    

