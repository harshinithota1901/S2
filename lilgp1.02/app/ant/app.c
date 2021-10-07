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

#include <lilgp.h>

globaldata g;

int app_build_function_sets ( void )
{
     function_set fset;
     int tree_map;
     char *tree_name;
     function sets[] = 
     { { f_if_food_ahead, NULL, NULL, 2, "if-food-ahead", FUNC_EXPR, -1, 0 },
       { f_progn2,        NULL, NULL, 2, "progn2",        FUNC_DATA, -1, 0 },
       { f_progn3,        NULL, NULL, 3, "progn3",        FUNC_DATA, -1, 0 },
       { f_move,          NULL, NULL, 0, "move",          TERM_NORM, -1, 0 },
       { f_right,         NULL, NULL, 0, "right",         TERM_NORM, -1, 0 },
       { f_left,          NULL, NULL, 0, "left",          TERM_NORM, -1, 0 },
       { f_progn4,        NULL, NULL, 4, "progn4",        FUNC_DATA, -1, 0 } };

     binary_parameter ( "app.use_progn4", 0 );
     if ( atoi ( get_parameter ( "app.use_progn4" ) ) )
          fset.size = 7;
     else
          fset.size = 6;
     fset.cset = sets;

     tree_map = 0;
     tree_name = "TREE";

     return function_sets_init ( &fset, 1, &tree_map, &tree_name, 1 );
}
void app_eval_fitness ( individual *ind )
{
     int i, j, k;

     set_current_individual ( ind );
     
     /* reset the map. */
     for ( i = 0; i < g.ysize; ++i )
          memcpy ( g.cmap[i], g.map[i], g.xsize );

     /* reset the simulation variables. */
     g.x = g.y = 0;
     g.dir = APP_EAST;
     g.food = 0;
     g.time = 0;

     /* evaluate the tree until the ant collects all the food
        or runs out of time. */
     while ( g.time < g.maxtime && g.food < g.maxfood )
          evaluate_tree ( ind->tr[0].data, 0 );

     ind->hits = g.food;
     ind->r_fitness = g.food;
     ind->s_fitness = g.maxfood - g.food;
     ind->a_fitness = 1.0/(1.0+ind->s_fitness);
     ind->evald = EVAL_CACHE_VALID;
     
}

int app_end_of_evaluation ( int gen, multipop *mpop, int newbest,
                           popstats *gen_stats, popstats *run_stats )
{
     individual temp;
     
     temp.tr = run_stats[0].best[0]->ind->tr;
     app_eval_fitness ( &temp );
     oprintf ( OUT_BST, 50, "BEST OF RUN INDIVIDUAL:\n" );
     oprintf ( OUT_BST, 50, "   time taken: %d/%d\n", g.time, g.maxtime );
     oprintf ( OUT_BST, 50, "food consumed: %d/%d\n", g.food, g.maxfood );
     oprintf ( OUT_BST, 50, "path: \n" );
     app_print_trail ( OUT_BST, g.cmap, g.xsize, g.ysize );
     
     return (run_stats[0].best[0]->ind->hits==g.maxfood);
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
     char *param;
     FILE *in;
     int i, j, k;

     param = get_parameter ( "app.maxtime" );
     if ( param == NULL )
          g.maxtime = 400;
     else
     {
          g.maxtime = atoi ( param );
          if ( g.maxtime < 1 )
               error ( E_FATAL_ERROR, "\"app.maxtime\" must be positive." );
     }
     
     param = get_parameter ( "app.trail" );
     if ( param == NULL )
          error ( E_FATAL_ERROR, "\"app.trail\" not specified." );

     in = fopen ( param, "r" );
     if ( in == NULL )
          error ( E_FATAL_ERROR, "can't open \"%s\".", param );
     
     fscanf ( in, "%d %d", &g.xsize, &g.ysize );
     if ( g.xsize < 1 || g.ysize < 1 )
          error ( E_FATAL_ERROR, "bad trail size." );

     g.map = (unsigned char **)MALLOC ( g.ysize * sizeof ( unsigned char * ) );
     g.cmap = (unsigned char **)MALLOC ( g.ysize * sizeof ( unsigned char * ));
     for ( i = 0; i < g.ysize; ++i )
     {
          g.map[i] = (unsigned char *)MALLOC ( g.xsize *
                                              sizeof ( unsigned char ) );
          g.cmap[i] = (unsigned char *)MALLOC ( g.xsize *
                                              sizeof ( unsigned char ) );
          memset ( g.map[i], APP_EMPTY, g.xsize );
     }
     
     i = j = 0;
     while ( (k=fgetc(in))!='\n' && k != EOF );
     g.maxfood = 0;

     while ( j < g.ysize )
     {
          k = fgetc(in);
          if ( k == EOF )
               j = g.ysize;
          else if ( k == '\n' )
          {
               i = 0;
               ++j;
          }
          else if ( k == '#' && i < g.xsize )
          {
               g.map[j][i] = APP_FOOD;
               ++i;
               ++g.maxfood;
          }
          else if ( i < g.xsize )
               ++i;
     }

     fclose ( in );

     app_print_trail ( OUT_PRG, g.map, g.xsize, g.ysize );

     return 0;
     
}

void app_print_trail ( int stream, unsigned char **m, int x, int y )
{
     int i, j;
     char string[] = ".#xENWS";

     for ( j = 0; j < y; ++j )
     {
          oprintf ( stream, 50, "     " );
          for ( i = 0; i < x; ++i )
               oprintf ( stream, 50, "%c", string[m[j][i]] );
          oprintf ( stream, 50, "\n" );
     }
}
               

void app_uninitialize ( void )
{
     int i;

     for ( i = 0; i < g.ysize; ++i )
     {
          FREE ( g.map[i] );
          FREE ( g.cmap[i] );
     }

     FREE ( g.map );
     FREE ( g.cmap );
}

void app_write_checkpoint ( FILE *f )
{
     return;
}

void app_read_checkpoint ( FILE *f )
{
     return;
}

               
     
     
     
