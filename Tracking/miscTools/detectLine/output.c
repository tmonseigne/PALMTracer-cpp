/*  Output the extracted lines as ASCII file; part of detect-lines.
    Copyright (C) 1996-1998 Carsten Steger

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */



#include "lines.h"



/* Output the results as ASCII table. */
void print_ascii(f, contours, num_cont, sigma, low, high, mode, has_width,
                 correct_pos, infile)
  FILE    *f;
  contour **contours;
  long    num_cont;
  double  sigma;
  double  low;
  double  high;
  long    mode;
  BOOL    has_width;
  BOOL    correct_pos;
  char    *infile;
{
  long    i, j, num_points;
  contour *cont;

  fprintf(f,"# Output of %s ",program_name);
  fprintf(f,"-s %.7g -l %.7g -u %.7g",sigma,low,high);
  if (mode == MODE_LIGHT)
    fprintf(f," -b");
  else
    fprintf(f," -d");
  if (!has_width)
    fprintf(f," -W");
  if (has_width && !correct_pos)
    fprintf(f," -C");
  if (infile != NULL)
    fprintf(f," %s",infile);
  else
    fprintf(f," -");
  fprintf(f,"\n\n");
  if (!has_width)
    fprintf(f,"no-");
  fprintf(f,"width\n");
  if (!correct_pos)
    fprintf(f,"no-");
  fprintf(f,"correct\n\n");
  fprintf(f,"%ld\n\n",num_cont);
  for (i=0; i<num_cont; i++) {
    cont = contours[i];
    num_points = cont->num;
    fprintf(f,"%ld\n",num_points);
    for (j=0; j<num_points; j++) {
      if (has_width && correct_pos)
        fprintf(f,"%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f\n",
                cont->row[j]+0.5,cont->col[j]+0.5,cont->angle[j],
                cont->width_l[j],cont->width_r[j],cont->contrast[j],
                cont->asymmetry[j]);
      else if (has_width && !correct_pos)
        fprintf(f,"%12.7f %12.7f %12.7f %12.7f %12.7f\n",cont->row[j]+0.5,
                cont->col[j]+0.5,cont->angle[j],cont->width_l[j],
                cont->width_r[j]);
      else
        fprintf(f,"%12.7f %12.7f %12.7f\n",cont->row[j]+0.5,cont->col[j]+0.5,
                cont->angle[j]);
    }
    fprintf(f,"\n");
  }
}
