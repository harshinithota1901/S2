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

DATATYPE f_vma ( int tree, farg *args )
{
     vector sum;
     sum.x = ( args[0].d.x + args[1].d.x ) % g.lawn_width;
     sum.y = ( args[0].d.y + args[1].d.y ) % g.lawn_height;
     return sum;
}

DATATYPE f_frog ( int tree, farg *args )
{
     vector jump;
     
     if ( g.abort )
          return args[0].d;
     
     switch ( g.dir )
     {
        case APP_NORTH:
          g.xpos += args[0].d.x;
          g.ypos -= args[0].d.y;
          break;
        case APP_SOUTH:
          g.xpos -= args[0].d.x;
          g.ypos += args[0].d.y;
          break;
        case APP_EAST:
          g.xpos += args[0].d.y;
          g.ypos += args[0].d.x;
          break;
        case APP_WEST:
          g.xpos -= args[0].d.y;
          g.ypos -= args[0].d.x;
          break;
     }
     g.xpos = ((g.xpos%g.lawn_width)+g.lawn_width)%g.lawn_width;
     g.ypos = ((g.ypos%g.lawn_height)+g.lawn_height)%g.lawn_height;
     g.lawn[g.xpos][g.ypos] = APP_MOWN;
     ++g.movecount;
     if ( g.movecount > 100 )
          g.abort = 1;

     return args[0].d;
}

DATATYPE f_prog2 ( int tree, farg *args )
{
     return args[1].d;
}

DATATYPE f_left ( int tree, farg *args )
{
     vector bob = {0,0};
     if ( g.abort )
          return bob;
     g.dir = (g.dir+1)%4;
     ++g.leftcount;
     if ( g.leftcount > 100 )
          g.abort = 1;
     return bob;
}

DATATYPE f_mow ( int tree, farg *args )
{
     vector bob = {0,0};
     
     if ( g.abort )
          return bob;
     
     switch ( g.dir )
     {
        case APP_NORTH:
          g.ypos = ( g.ypos + g.lawn_height - 1 ) % g.lawn_height;
          break;
        case APP_SOUTH:
          g.ypos = ( g.ypos + 1 ) % g.lawn_height;
          break;
        case APP_EAST:
          g.xpos = ( g.xpos + 1 ) % g.lawn_width;
          break;
        case APP_WEST:
          g.xpos = ( g.xpos + g.lawn_width - 1 ) % g.lawn_width;
          break;
     }
     g.lawn[g.xpos][g.ypos] = APP_MOWN;
     ++g.movecount;
     if ( g.movecount > 100 )
          g.abort = 1;
     return bob;
}

void f_vecgen ( DATATYPE *r )
{
     r->x = random_int ( g.lawn_width );
     r->y = random_int ( g.lawn_height );
}

char *f_vecstr ( DATATYPE d )
{
     static char buffer[20];

     sprintf ( buffer, "(%d,%d)", d.x, d.y );
     return buffer;
}
          
          
          
          
     
