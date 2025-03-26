/*  Read a PGM file; part of detect-lines.
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



/* Local function prototypes */

static void skip_whitespace __P((FILE *f));
static void read_non_whitespace __P((FILE *f, char *buf, size_t buf_size));



/* Skip whitespace in the pgm header. */
static void skip_whitespace(f)
  FILE *f;
{
  int c;

  c = fgetc(f);
  while (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '#') {
    while (c == ' ' || c == '\t' || c == '\n' || c == '\r')
      c = fgetc(f);
    if (c == '#') {
      while (c != '\n' && c != '\r')
        c = fgetc(f);
    }
  }
  ungetc(c,f);
}



static void read_non_whitespace(f,buf,buf_size)
  FILE   *f;
  char   *buf;
  size_t buf_size;
{
  int c, i;

  i = 0;
  c = fgetc(f);
  while (c != ' ' && c != '\t' && c != '\n' && c != '\r' && c != '#') {
    if (i >= (long)buf_size-1)
      break;
    buf[i] = c;
    c = fgetc(f);
    i++;
  }
  buf[i] = '\0';
  ungetc(c,f);
}



/* Read a pgm image from file f and return it and its width and height.  The
   caller is responsible for freeing the memory of the image. */
void read_pgm(f,image,width,height)
  FILE *f;
  byte **image;           
  long *width;
  long *height;
{  
  long i, j, k, n, val, max_gray;
  int  c;
  char buf[1024], filetype[2];
  BOOL is_ascii;

  /* Determine pgm format */
  n = fread(filetype,sizeof(char),2,f);
  if (n < 2)
    print_error(ERR_NOPGM,"no header");
  if (!strncmp(filetype,"P2",2))
    is_ascii = TRUE;
  else if (!strncmp(filetype,"P5",2))
    is_ascii = FALSE;
  else
    print_error(ERR_NOPGM,"not a pgm header");

  skip_whitespace(f);

  read_non_whitespace(f,buf,1024);
  sscanf(buf,"%ld",width);
  if (*width <= 0)
    print_error(ERR_NOPGM,"width <= 0");
  skip_whitespace(f);
  read_non_whitespace(f,buf,1024);
  sscanf(buf,"%ld",height);
  if (*height <= 0)
    print_error(ERR_NOPGM,"height <= 0");

  skip_whitespace(f);

  read_non_whitespace(f,buf,1024);
  sscanf(buf,"%ld",&max_gray);
  if (max_gray <= 0)
    print_error(ERR_NOPGM,"maximum gray value <= 0");
  if (max_gray > 255)
    print_error(ERR_NOPGM,"maximum gray value > 255");

  c = fgetc(f); /* Read the newline character. */

  *image = xcalloc((*width)*(*height),sizeof(byte));
  
  if (is_ascii) {
    k = 0;
    for (i=0; i<*height; i++) {
      for (j=0; j<*width; j++) {
        if (fscanf(f,"%ld",&val) == EOF)
          print_error(ERR_NOPGM,"premature end of file");
        (*image)[k++] = (byte)val;
      }
    }
  } else {
    n = fread(*image,1,(*width)*(*height),f);
    if (n != (*width)*(*height))
      print_error(ERR_NOPGM,"premature end of file");
  }
}
