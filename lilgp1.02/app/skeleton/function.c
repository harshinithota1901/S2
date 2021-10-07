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

#include <math.h>
#include <stdio.h>

#include <lilgp.h>

/* here you need one C function for every ordinary (non-ERC) function
 * and terminal.  each should take an int and array of type (farg) and
 * return a type (DATATYPE).
 *
 * functions declared as type DATA can get the value of child i by
 * reading "args[i].d".
 *
 * functions declared as type EXPR can evaluate the child and gets its
 * value by passing "args[i].t" (and the int argument) to evaluate_tree().
 *
 * terminals should not need to use either argument passed to them.
 *
 */

DATATYPE f_SAMPLE_data_function ( int tree, farg *args )
{
     /* just return the value of my child. */
     return args[0].d;
}

DATATYPE f_SAMPLE_expr_function ( int tree, farg *args )
{
     /* evaluate my child and return its value. */
     return evaluate_tree ( args[0].t, tree );
}

DATATYPE f_SAMPLE_terminal ( int tree, farg *args )
{
     /* just return 0. */
     return (DATATYPE)0;
}

/* here you need for every terminal that is a random constant, a C function
 * to generate a new constant and a C function to print the constant to
 * a string.
 *
 * the generation function should take a (DATATYPE *), into which it places
 *     the generated constant, and should return void.
 * the printing function should take a (DATATYPE), which it should print
 *     into a (char *) that it returns.  the character buffer should be
 *     declared static.
 */

void f_SAMPLE_erc_gen ( DATATYPE *v )
{
     *v = random_int ( 10 );
}

char *f_SAMPLE_erc_str ( DATATYPE v )
{
     static char buffer[20];

     sprintf ( buffer, "%d", v );
     return buffer;
}

     
     

