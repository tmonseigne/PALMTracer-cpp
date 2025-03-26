/*  detect-lines, extract lines and their width from images.
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



/* Type declarations */

/* Command line options for the program */
typedef struct {
  double sigma;
  double low;
  double high;
  long   mode;
  BOOL   correct;
  BOOL   width;
  BOOL   extend;
  BOOL   postscript;
  BOOL   encapsulated;
  BOOL   image;
  char   *infile;
  char   *outfile;
} options;



/* Global variables */

/* The basename of the program as called from the shell */
char *program_name;



/* Local function prototypes */

static void parse_args __P((int argc, char **argv, options *opts));
static void print_usage __P((void));
static void print_version __P((void));
static void check_sigma __P((double sigma, long width, long height));



/* Print usage information 
static void print_usage()
{
  fprintf(stderr,"usage: %s -s sigma -l low -u high ",program_name);
  fprintf(stderr,"[-bdwWcCjJpeihv] [pgm-image] [output-file]\n");
  fprintf(stderr," -s --sigma             value of sigma for the Gaussian ");
  fprintf(stderr,"smoothing masks\n");
  fprintf(stderr," -l --low               lower hysteresis threshold\n");
  fprintf(stderr," -u --high              upper hysteresis threshold\n");
  fprintf(stderr," -b --bright            extract bright lines (default)\n");
  fprintf(stderr," -d --dark              extract dark lines\n");
  fprintf(stderr," -w --width             extract line width (default)\n");
  fprintf(stderr," -W --no-width          do not extract line width\n");
  fprintf(stderr," -c --correct           correct line width and position ");
  fprintf(stderr,"(default)\n");
  fprintf(stderr," -C --no-correct        do not correct line width and ");
  fprintf(stderr,"position\n");
  fprintf(stderr," -j --add-junctions     add junctions by extending lines\n");
  fprintf(stderr," -J --no-add-junctions  do not add junctions (default)\n");
  fprintf(stderr," -p --postscript        output results as PostScript ");
  fprintf(stderr,"file\n");
  fprintf(stderr," -e --epsf              output results as Encapsulated ");
  fprintf(stderr,"PostScript file\n");
  fprintf(stderr," -i --include-image     include image in PostScript file\n");
  fprintf(stderr," -h --help              give this help\n");
  fprintf(stderr," -v --version           display version number\n");
  fprintf(stderr," pgm-image              input image (if not given, use ");
  fprintf(stderr,"standard input)\n");
  fprintf(stderr," output-file            output file (if not given, use ");
  fprintf(stderr,"standard output)\n");
  exit(0);
}

*/

/* Print version information */
static void print_version()
{
  fprintf(stderr,"%s %d.%d\n",program_name,VERSION,REVISION);
  exit(0);
}



/* Parse command line arguments and return the selected options in opts.
   This routine is not as smart as the GNU getopt_long since options and
   their arguments must always be separated by whitespace and short options
   must not be clustered together, e.g., -pi is forbidden. 
static void parse_args(argc, argv, opts)
  int     argc;
  char    **argv;
  options *opts;
{
  long i, n;
  bool sigma_found = FALSE;
  bool low_found = FALSE;
  bool high_found = FALSE;
  bool infile_found = FALSE;
  bool outfile_found = FALSE;

  // We need to know at least sigma, low, and high. 
  if (argc == 1)
    print_usage();

  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) {
      print_usage();
    } else if (!strcmp(argv[i],"-v") || !strcmp(argv[i],"--version")) {
      print_version();
    } else if (!strcmp(argv[i],"-s") || !strcmp(argv[i],"--sigma")) {
      i++;
      if (i >= argc)
        print_error(ERR_SNS,NULL);
      n = sscanf(argv[i],"%lf",&opts->sigma);
      if (n < 1)
        print_error(ERR_SNS,NULL);
      sigma_found = TRUE;
    } else if (!strcmp(argv[i],"-l") || !strcmp(argv[i],"--low")) {
      i++;
      if (i >= argc)
        print_error(ERR_LNS,NULL);
      n = sscanf(argv[i],"%lf",&opts->low);
      if (n < 1)
        print_error(ERR_LNS,NULL);
      low_found = TRUE;
    } else if (!strcmp(argv[i],"-u") || !strcmp(argv[i],"--high")) {
      i++;
      if (i >= argc)
        print_error(ERR_HNS,NULL);
      n = sscanf(argv[i],"%lf",&opts->high);
      if (n < 1)
        print_error(ERR_HNS,NULL);
      high_found = TRUE;
    } else if (!strcmp(argv[i],"-b") || !strcmp(argv[i],"--bright")) {
      opts->mode = MODE_LIGHT;
    } else if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"--dark")) {
      opts->mode = MODE_DARK;
    } else if (!strcmp(argv[i],"-c") || !strcmp(argv[i],"--correct")) {
      opts->correct = TRUE;
      opts->width = TRUE;
    } else if (!strcmp(argv[i],"-C") || !strcmp(argv[i],"--no-correct")) {
      opts->correct = FALSE;
    } else if (!strcmp(argv[i],"-w") || !strcmp(argv[i],"--width")) {
      opts->width = TRUE;
    } else if (!strcmp(argv[i],"-W") || !strcmp(argv[i],"--no-width")) {
      opts->width = FALSE;
      opts->correct = FALSE;
    } else if (!strcmp(argv[i],"-j") || !strcmp(argv[i],"--add-junctions")) {
      opts->extend = TRUE;
    } else if (!strcmp(argv[i],"-J") ||
               !strcmp(argv[i],"--no-add-junctions")) {
      opts->extend = FALSE;
    } else if (!strcmp(argv[i],"-p") || !strcmp(argv[i],"--postscript")) {
      opts->postscript = TRUE;
    } else if (!strcmp(argv[i],"-e") || !strcmp(argv[i],"--epsf")) {
      opts->postscript = TRUE;
      opts->encapsulated = TRUE;
    } else if (!strcmp(argv[i],"-i") || !strcmp(argv[i],"--include-image")) {
      opts->image = TRUE;
    } else if (!strcmp(argv[i],"-")) {
      // Make stdin or stdout the respective input or output file.
      if (!infile_found) {
        opts->infile = NULL;
        infile_found = TRUE;
      } else if (!outfile_found) {
        opts->outfile = NULL;
        outfile_found = TRUE;
      } else {
        print_error(ERR_TMF,argv[i]);
      }
    } else if (argv[i][0] == '-') {
      print_error(ERR_UKO,argv[i]);
    } else {
      if (!infile_found) {
        opts->infile = argv[i];
        infile_found = TRUE;
      } else if (!outfile_found) {
        opts->outfile = argv[i];
        outfile_found = TRUE;
      } else {
        print_error(ERR_TMF,argv[i]);
      }
    }
  }
  if (!sigma_found)
    print_error(ERR_SNS,NULL);
  if (!low_found)
    print_error(ERR_LNS,NULL);
  if (!high_found)
    print_error(ERR_HNS,NULL);
  if (opts->sigma < 0)
    print_error(ERR_SOR,"< 0");
  if (opts->low < 0)
    print_error(ERR_LOR,"< 0");
  if (opts->high < 0)
    print_error(ERR_HOR,"< 0");
  if (opts->high < opts->low)
    print_error(ERR_LGH,NULL);
  if (opts->correct && !opts->width)
    print_error(ERR_CNW,NULL);
  if (opts->image && !opts->postscript)
    print_error(ERR_INP,NULL);
}*/



/* Check whether sigma is in the allowed range for the current image size. */
static void check_sigma(sigma, width, height)
  double sigma;
  long width;
  long height;
{
  long min_dim;

  min_dim = width < height ? width : height;
  if (sigma < 0.4)
    print_error(ERR_SOR,"< 0.4");
  if (MASK_SIZE(MAX_SIZE_MASK_2,sigma) >= min_dim)
    print_error(ERR_SOR,"too large for image size");
}



/* The main program 
int main(argc, argv)
  int  argc;
  char **argv;
{
  byte    *image;
  long    width, height;
  contour **contours, *cont;
  long    i, num_cont;
  FILE    *inf, *outf;
  char    err_str[200];
  static options opts = {
    -1.0, -1.0, -1.0, MODE_LIGHT, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
    NULL, NULL
  };

  program_name = strrchr(argv[0],'/');
  if (program_name == NULL)
    program_name = argv[0];
  else
    program_name++;

  parse_args(argc,argv,&opts);

  if (opts.infile != NULL) {
    inf = fopen(opts.infile,"rb");
    if (inf == NULL) {
      sprintf(err_str,"%s: %s",opts.infile,sys_errlist[errno]);
      print_error(ERR_FNF,err_str);
    }
  } else
    inf = stdin;

  if (opts.outfile != NULL) {
    outf = fopen(opts.outfile,"w");
    if (outf == NULL) {
      sprintf(err_str,"%s: %s",opts.outfile,sys_errlist[errno]);
      print_error(ERR_FNF,err_str);
    }
  } else
    outf = stdout;

  read_pgm(inf,&image,&width,&height);

  check_sigma(opts.sigma,width,height);

  detect_lines(image,width,height,&contours,&num_cont,opts.sigma,opts.low,
               opts.high,opts.mode,opts.width,opts.correct,opts.extend);





  if (opts.postscript)
    print_postscript(outf,image,width,height,opts.image,contours,num_cont,
                     opts.encapsulated,opts.sigma,opts.low,opts.high,opts.mode,
                     opts.width,opts.correct,opts.infile);
  else
    print_ascii(outf,contours,num_cont,opts.sigma,opts.low,opts.high,opts.mode,
                opts.width,opts.correct,opts.infile);

  printf("Nombre de contours: %d\n",num_cont);
  for (i=0; i<num_cont; i++) {
    
    cont = contours[i];
    printf("Contour N0 %d :%d points\n",i+1,(int)cont->num);
    if (cont->width_r != NULL)
      free(cont->width_r);
    if (cont->width_l != NULL)
      free(cont->width_l);
    if (cont->angle != NULL)
      free(cont->angle);
    if (cont->col != NULL)
      free(cont->col);
    if (cont->row != NULL)
      free(cont->row);
    free(cont);
  }
  free(contours);
  free(image);

  if (inf != stdin)
    fclose(inf);
  if (outf != stdout)
    fclose(outf);

  return 0;
}
*/




int detectLines(unsigned char *image, long width, long height, contour ***contours, float sigma, float low){


  long    num_cont=0;

 
  /*static options opts = {
    -1.0, -1.0, -1.0, MODE_LIGHT, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
    NULL, NULL
  };*/

  check_sigma(sigma,width,height);

  detect_lines(image,width,height,contours,&num_cont,sigma,low,
               low*3.,MODE_LIGHT,FALSE,FALSE,FALSE);


  return (num_cont);
}



void deleteContours(long num_cont, contour ***contours){
	int i;
	contour *cont;
  for (i=0; i<num_cont; i++) {
    
    cont = (*contours)[i];
    if (cont->width_r != NULL)
      free(cont->width_r);
    if (cont->width_l != NULL)
      free(cont->width_l);
    if (cont->angle != NULL)
      free(cont->angle);
    if (cont->col != NULL)
      free(cont->col);
    if (cont->row != NULL)
      free(cont->row);
    free(cont);
  }
  free(*contours);
}