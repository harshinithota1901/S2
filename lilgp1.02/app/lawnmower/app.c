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
     function_set *fset;
     int *tree_map;
     char *tree_name[3];
     int ret;

     function sets[4][10] =

            /*** RPB ***/
     { { { f_left,  NULL,     NULL,     0, "left",  TERM_NORM, -1, 0 },
         { f_mow,   NULL,     NULL,     0, "mow",   TERM_NORM, -1, 0 },
         { NULL,    f_vecgen, f_vecstr, 0, "Rvm",   TERM_ERC,  -1, 0 },
              
         { f_frog,  NULL,     NULL,     1, "frog",  FUNC_DATA, -1, 0 },
         { f_vma,   NULL,     NULL,     2, "vma",   FUNC_DATA, -1, 0 },
         { f_prog2, NULL,     NULL,     2, "prog2", FUNC_DATA, -1, 0 },
         { NULL,    NULL,     NULL,    -1, "ADF0",  EVAL_DATA,  1, 0 },
         { NULL,    NULL,     NULL,    -1, "ADF1",  EVAL_DATA,  2, 0 } },

            /*** ADF0 ***/
       { { f_vma,   NULL,     NULL,     2, "vma",   FUNC_DATA, -1, 0 },
         { f_prog2, NULL,     NULL,     2, "prog2", FUNC_DATA, -1, 0 },
              
         { f_left,  NULL,     NULL,     0, "left",  TERM_NORM, -1, 0 },
         { f_mow,   NULL,     NULL,     0, "mow",   TERM_NORM, -1, 0 },
         { NULL,    f_vecgen, f_vecstr, 0, "Rvm2",  TERM_ERC,  -1, 0 } },


            /*** ADF1 ***/
       { { f_left,  NULL,     NULL,     0, "left",  TERM_NORM, -1, 0 },
         { f_mow,   NULL,     NULL,     0, "mow",   TERM_NORM, -1, 0 },
         { NULL,    f_vecgen, f_vecstr, 0, "Rvm",   TERM_ERC,  -1, 0 },
         { NULL,    NULL,     NULL,     0, "ARG0",  TERM_ARG,   0, 0 },

         { f_frog,  NULL,     NULL,     1, "frog",  FUNC_DATA, -1, 0 },
         { f_vma,   NULL,     NULL,     2, "vma",   FUNC_DATA, -1, 0 },
         { f_prog2, NULL,     NULL,     2, "prog2", FUNC_DATA, -1, 0 },
         { NULL,    NULL,     NULL,    -1, "ADF0",  EVAL_DATA,  1, 0 } },
            
            /*** non-ADF function set ***/
       { { f_left,  NULL,     NULL,     0, "left",  TERM_NORM, -1, 0 },
         { f_mow,   NULL,     NULL,     0, "mow",   TERM_NORM, -1, 0 },
         { NULL,    f_vecgen, f_vecstr, 0, "Rvm",   TERM_ERC,  -1, 0 },
              
         { f_frog,  NULL,     NULL,     1, "frog",  FUNC_DATA, -1, 0 },
         { f_vma,   NULL,     NULL,     2, "vma",   FUNC_DATA, -1, 0 },
         { f_prog2, NULL,     NULL,     2, "prog2", FUNC_DATA, -1, 0 } } };

     binary_parameter ( "app.use_adfs", 1 );
     if ( atoi ( get_parameter ( "app.use_adfs" ) ) )
     {

          fset = (function_set *)MALLOC ( 3 * sizeof ( function_set ) );
          fset[0].size = 8;
          fset[0].cset = sets[0];
          fset[1].size = 5;
          fset[1].cset = sets[1];
          fset[2].size = 8;
          fset[2].cset = sets[2];
          
          tree_map = (int *)MALLOC ( 3 * sizeof ( int ) );
          tree_map[0] = 0;
          tree_name[0] = "RPB";
          tree_map[1] = 1;
          tree_name[1] = "ADF0";
          tree_map[2] = 2;
          tree_name[2] = "ADF1";
          
          ret = function_sets_init ( fset, 3, tree_map, tree_name, 3 );
     }
     else
     {
          fset = (function_set *)MALLOC ( 1 * sizeof ( function_set ) );
          fset[0].size = 6;
          fset[0].cset = sets[3];
          
          tree_map = (int *)MALLOC ( 1 * sizeof ( int ) );
          tree_map[0] = 0;
          tree_name[0] = "TREE";
          
          ret = function_sets_init ( fset, 1, tree_map, tree_name, 1 );
     }

     FREE ( tree_map );
     FREE ( fset );

     return ret;
}

void app_eval_fitness ( individual *ind )
{

     int i, j, k;
     DATATYPE value;

     set_current_individual ( ind );

     g.xpos = 4;
     g.ypos = 4;
     g.dir = APP_NORTH;
     g.leftcount = 0;
     g.movecount = 0;
     g.abort = 0;

     for ( i = 0; i < g.lawn_width; ++i )
          for ( j = 0; j < g.lawn_height; ++j )
               g.lawn[i][j] = APP_UNMOWN;
     
     /* evaluate the tree. */
     evaluate_tree ( ind->tr[0].data, 0 );

     k = 0;
     for ( i = 0; i < g.lawn_width; ++i )
          for ( j = 0; j < g.lawn_height; ++j )
               k += (g.lawn[i][j] == APP_MOWN);

     ind->hits = k;
     ind->r_fitness = k;
        
     /* compute the standardized and raw fitness. */
     ind->s_fitness = ( g.lawn_width * g.lawn_height ) - ind->r_fitness;
     ind->a_fitness = 1/(1+ind->s_fitness);

     /* always leave this line in. */
     ind->evald = EVAL_CACHE_VALID;
}

int app_end_of_evaluation ( int gen, multipop *mpop, int newbest,
                           popstats *gen_stats, popstats *run_stats )
{
     return run_stats[0].best[0]->ind->hits == g.lawn_width * g.lawn_height;
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
     int i;
     char *param;
     
     param = get_parameter ( "app.lawn_width" );
     if ( param == NULL )
          g.lawn_width = APP_DEFAULT_LAWN_WIDTH;
     else
     {
          g.lawn_width = atoi ( param );
          if ( g.lawn_width < 1 )
               error ( E_FATAL_ERROR, "invalid value for \"app.lawn_width\"." );
     }
     
     param = get_parameter ( "app.lawn_height" );
     if ( param == NULL )
          g.lawn_height = APP_DEFAULT_LAWN_HEIGHT;
     else
     {
          g.lawn_height = atoi ( param );
          if ( g.lawn_height < 1 )
               error ( E_FATAL_ERROR, "invalid value for \"app.lawn_height\"." );
     }

     g.lawn = (unsigned char **)MALLOC ( g.lawn_width *
					sizeof ( unsigned char * ) );
     for ( i = 0; i < g.lawn_width; ++i )
          g.lawn[i] = (unsigned char *)MALLOC ( g.lawn_height *
                                               sizeof ( unsigned char ) );

     return 0;
}

/* app_uninitialize()
 *
 * perform application cleanup (free memory, etc.)
 */

void app_uninitialize ( void )
{
     int i;

     for ( i = 0; i < g.lawn_width; ++i )
          FREE ( g.lawn[i] );
     FREE ( g.lawn );
     
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

    
