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

DATATYPE f_if_food_ahead ( int tree, farg *args )
{
     int spot;

     if ( g.time >= g.maxtime )
          return 0;
     
     switch ( g.dir )
     {
        case APP_EAST:
          spot = g.cmap[g.y][(g.x+1)%g.xsize];
          break;
        case APP_WEST:
          spot = g.cmap[g.y][(g.x-1+g.xsize)%g.xsize];
          break;
        case APP_SOUTH:
          spot = g.cmap[(g.y+1)%g.ysize][g.x];
          break;
        case APP_NORTH:
          spot = g.cmap[(g.y-1+g.ysize)%g.ysize][g.x];
          break;
     }

     if ( spot == APP_FOOD )
          evaluate_tree ( args[0].t, tree );
     else
          evaluate_tree ( args[1].t, tree );

     return 0;
}

DATATYPE f_progn2 ( int tree, farg *args )
{
     return 0;
}

DATATYPE f_progn3 ( int tree, farg *args )
{
     return 0;
}

DATATYPE f_progn4 ( int tree, farg *args )
{
     return 0;
}

DATATYPE f_move ( int tree, farg *args )
{
     if ( g.time >= g.maxtime )
          return 0;

     ++g.time;
     
     g.cmap[g.y][g.x] = APP_PATH;
     
     switch ( g.dir )
     {
        case APP_EAST:
          ++g.x;
          break;
        case APP_WEST:
          g.x += g.xsize-1;
          break;
        case APP_SOUTH:
          ++g.y;
          break;
        case APP_NORTH:
          g.y += g.ysize-1;
          break;
     }

     g.x %= g.xsize;
     g.y %= g.ysize;

     if ( g.cmap[g.y][g.x] == APP_FOOD )
          ++g.food;

     g.cmap[g.y][g.x] = APP_ANT + g.dir;
}

DATATYPE f_right ( int tree, farg *args )
{
     if ( g.time >= g.maxtime )
          return 0;
     ++g.time;
     g.dir = (g.dir+3)%4;
     return 0;
}

DATATYPE f_left ( int tree, farg *args )
{
     if ( g.time >= g.maxtime )
          return 0;
     ++g.time;
     g.dir = (g.dir+1)%4;
     return 0;
}

