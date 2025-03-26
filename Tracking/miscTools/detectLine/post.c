/*  Output the extracted lines as a PostScript file; part of detect-lines.
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



/* Output results as PostScript file.  The parameter include_image determines
   whether the input image should be output along with the results, while
   as_epsf determines whether an Encapsulated PostScript file should be
   generated. */
void print_postscript(f, image, width, height, include_image, contours,
                      num_cont, as_epsf, sigma, low, high, mode, has_width,
                      correct_pos, infile)
  FILE    *f;
  byte    *image;
  long    width;
  long    height;
  BOOL    include_image;
  contour **contours;
  long    num_cont;
  BOOL    as_epsf;
  double  sigma;
  double  low;
  double  high;
  long    mode;
  BOOL    has_width;
  BOOL    correct_pos;
  char    *infile;
{
  long    i, j, r, c, l, num_points;
  int     g;
  contour *cont;
  double  px, py, nx, ny, px_r, py_r, px_l, py_l;
  double  last_px_r, last_py_r, last_px_l, last_py_l;
  double  last_px, last_py, last_w_r, last_w_l;
  double  w, h, s, sx, sy, tx, ty;
  time_t  current_time;
#ifdef HAVE_UNISTD_H
  char    *login_name, *user_name, *office;
  struct passwd *pwent;
#endif

  sx = 500.0/width;
  sy = 700.0/height;
  if (sx < sy) s = sx;
  else s = sy;
  w = s*width;
  h = s*height;
  if (as_epsf) {
    tx = 0;
    ty = 0;
  } else {
    tx = floor((595-w)/2);
    ty = floor((842-h)/2);
  }
  fprintf(f,"%%!PS-Adobe-3.0");
  if (as_epsf)
     fprintf(f," EPSF-3.0\n");
  else
    fprintf(f,"\n");
  fprintf(f,"%%%%BoundingBox: %.0f %.0f %.0f %.0f\n",tx,ty,tx+floor(w),
          ty+floor(h));
  time(&current_time);
  fprintf(f,"%%%%CreationDate: %s",ctime(&current_time));
  fprintf(f,"%%%%Creator: %s %d.%d\n",program_name,VERSION,REVISION);
  fprintf(f,"%%%%Copyright: %s is Copyright (C) 1996-1998 ",program_name);
  fprintf(f,"Carsten Steger\n");
  fprintf(f,"%%%%Title: Output of %s ",program_name);
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
  fprintf(f,"\n");
#ifdef HAVE_UNISTD_H
  login_name = getlogin();
  if (login_name != NULL) {
    pwent = getpwnam(login_name);
    if (pwent != NULL) {
      user_name = pwent->pw_gecos;
      office = strchr(user_name,',');
      if (office != NULL)
        *office = '\0';
      fprintf(f,"%%%%For: %s (%s)\n",user_name,login_name);
    }
  }
#endif
  if (!as_epsf)
    fprintf(f,"%%%%Pages: 1\n");
  fprintf(f,"%%%%EndComments\n");
  fprintf(f,"%%%%BeginProlog\n");
  fprintf(f,"%%%%BeginResource: procset ");
  fprintf(f,"%s %d %d\n",program_name,VERSION,REVISION);
  fprintf(f,"/LinesDict 10 dict def\n");
  fprintf(f,"LinesDict begin\n");
  fprintf(f,"/L {1 0 0 setrgbcolor newpath moveto lineto stroke} bind def\n");
  fprintf(f,"/E {0 1 0 setrgbcolor newpath moveto lineto stroke} bind def\n");
  fprintf(f,"/W %ld def\n",width);
  fprintf(f,"/H %ld def\n",height);
  fprintf(f,"/B W string def\n");
  fprintf(f,"/I {W H 8 [W 0 0 H neg 0 H] {currentfile B readhexstring pop}");
  fprintf(f," image} bind def\n");
  fprintf(f,"/C {newpath 0 0 moveto 0 1 lineto 1 1 lineto 1 0 lineto\n");
  fprintf(f,"    closepath clip} bind def\n");
  fprintf(f,"end\n");
  fprintf(f,"/LinesBegin {LinesDict begin /LinesState save def} def\n");
  fprintf(f,"/LinesEnd {LinesState restore end} def\n");
  fprintf(f,"%%%%EndResource\n");
  fprintf(f,"%%%%EndProlog\n");
  if (!as_epsf)
    fprintf(f,"\n%%%%Page: 1 1");
  fprintf(f,"\nLinesBegin\n");
  fprintf(f,"%.0f %.0f translate %.7g %.7g scale\n",tx,ty,w,h);
  if (contours != NULL)
    fprintf(f,"1 setlinecap 0.1 setlinewidth\n");
  fprintf(f,"C\n");

  if (include_image) {
    fprintf(f,"I\n");
    for (r=0; r<height; r++) {
      for (c=0; c<width; c++) {
        l = LINCOOR(r,c,width);
        g = image[l];
        fprintf(f,"%02x",g);
        if (c % 36 == 35)
          fprintf(f,"\n");
      }
      fprintf(f,"\n");
    }
  }

  if (contours != NULL) {
    fprintf(f,"0 1 translate %.7g %.7g scale\n",1.0/width,-1.0/height);
    for (i=0; i<num_cont; i++) {
      cont = contours[i];
      num_points = cont->num;
      last_w_r = last_w_l = 0;
      for (j=0; j<num_points; j++) {
        px = cont->row[j];
        py = cont->col[j];
        nx = cos(cont->angle[j]);
        ny = sin(cont->angle[j]);
        if (has_width) {
          px_r = px+cont->width_r[j]*nx;
          py_r = py+cont->width_r[j]*ny;
          px_l = px-cont->width_l[j]*nx;
          py_l = py-cont->width_l[j]*ny;
        }
        if (j > 0) {
          fprintf(f,"%.5g %.5g %.5g %.5g L\n",
                  last_py+0.5,last_px+0.5,py+0.5,px+0.5);
        }
        if (has_width) {
          if (last_w_r > 0 && cont->width_r[j] > 0) {
            fprintf(f,"%.5g %.5g %.5g %.5g E\n",last_py_r+0.5,last_px_r+0.5,
                    py_r+0.5,px_r+0.5);
          }
          if (last_w_l > 0 && cont->width_l[j] > 0) {
            fprintf(f,"%.5g %.5g %.5g %.5g E\n",last_py_l+0.5,last_px_l+0.5,
                    py_l+0.5,px_l+0.5);
          }
        }
        last_px = px;
        last_py = py;
        if (has_width) {
          last_w_r = cont->width_r[j];
          last_px_r = px_r;
          last_py_r = py_r;
          last_w_l = cont->width_l[j];
          last_px_l = px_l;
          last_py_l = py_l;
        }
      }
    }
  }

  fprintf(f,"LinesEnd\n");
  fprintf(f,"showpage\n");
  fprintf(f,"%%%%EOF\n");
}
