/* -*- mode: c++; c-basic-offset: 3 -*-
 * @(#)pandore.h		PANDORE		2001-04-03
 *
 * PANTHEON Project
 *
 * GREYC IMAGE
 * 6 Boulevard Maréchal Juin
 * F-14050 Caen Cedex France
 *
 * This file is free software. You can use it, distribute it
 * and/or modify it. However, the entire risk to the quality
 * and performance of this program is with you.
 *
 * Please email any bugs, comments, and/or additions to this file to:
 * Regis.Clouard@greyc.ensicaen.fr
 * http://www.greyc.ensicaen.fr/~regis
 */

/*
 * (C)Régis Clouard - 1999-10-08
 * (R)Régis Clouard - 2001-04-03 (version 3.00)
 */

#ifndef	__PPANDOREH__
#define	__PPANDOREH__

#ifndef _WIN32
#include <unistd.h>
#endif
#include <stdlib.h>
#include <fcntl.h>

#include <string>

/*
 * The general constants.
 */
#ifndef NULL
#define	NULL	(0L)
#endif
#ifndef TRUE
#define TRUE	true
#endif
#ifndef FALSE
#define FALSE	false
#endif
#ifndef RESFIL
#define RESFIL  "/.pandore"
#endif

/*
 * The basic scalar types.
 * (`BASE' is a flag for the
 * preprocessor macro `forallbase')
 * (`IMAGE' is a flag for the
 * preprocessor macro `forallbaseimage')
 * DO NOT DELETE this comment!
 */
typedef char Char;             // BASE
typedef unsigned char Uchar;   // BASE IMAGE
typedef short Short;           // BASE
typedef unsigned short Ushort; // BASE
typedef long Long;             // BASE IMAGE
typedef unsigned long Ulong;   // BASE
typedef float Float;           // BASE IMAGE
typedef double Double;         // BASE

/*
 * Definition des variables de limites de valeurs.
 */
#include <limits.h>
#include <float.h>

#ifndef MINCHAR
#define	MINCHAR SCHAR_MIN
#endif
#ifndef MAXCHAR
#define	MAXCHAR SCHAR_MAX
#endif
#ifndef MAXUCHAR
#define	MAXUCHAR UCHAR_MAX
#endif

#ifndef MINSHORT
#define	MINSHORT SHRT_MIN
#endif
#ifndef MAXSHORT
#define	MAXSHORT SHRT_MAX
#endif
#ifndef MAXUSHORT
#define	MAXUSHORT USHRT_MAX
#endif

#ifndef MINLONG
#define	MINLONG LONG_MIN
#endif
#ifndef MAXLONG
#define	MAXLONG LONG_MAX
#endif
#ifndef MAXULONG
#define	MAXULONG ULONG_MAX
#endif

#ifndef MINFLOAT
#define	MINFLOAT FLT_MIN
#endif
#ifndef MAXFLOAT
#define	MAXFLOAT FLT_MAX
#endif

#ifndef MINDOUBLE
#define	MINDOUBLE DBL_MIN
#endif
#ifndef MAXDOUBLE
#define	MAXDOUBLE DBL_MAX
#endif

// Needed for cerr...
#include <iostream>

//using namespace std;

//#include "errc.h"
/*#include "main.h"
#include "panfile.h"
#include "pobject.h"
#include "point.h"
#include "neighb.h"
#include "image.h"
#include "region.h"
#include "graph.h"
#include "plist.h"
#include "promotions.h"
#include "plimits.h"
#include "collection.h"*/

/*
 *  Some macro functions
 */
#ifndef MIN
#define	MIN(x,y) (((x)<(y))? (x) : (y))
#endif
#ifndef MAX
#define	MAX(x,y) (((x)>(y))? (x) : (y))
#endif
#ifndef SQR
#define	SQR(x) ((x)*(x))
#endif
#ifndef CUB
#define	CUB(x) ((x)*(x)*(x))
#endif
#ifndef ABS
#define	ABS(x) (((x)>0)? (x) : (-1*(x)))
#endif

#define	WHITE 255
#define	BLACK 0

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Compatibility with previous Pandore versions...

/*typedef Pobject Pobject3d;	
typedef PColorSpace EspaceCoul;
#define RVB RGB
typedef class GEdge Lnode;*/

#endif // __PPANDOREH__
