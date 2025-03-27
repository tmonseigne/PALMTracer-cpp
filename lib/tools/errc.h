/* -*- mode: c++; c-basic-offset: 3 -*-
 * @(#)errc.h			PANDORE		2001-04-03
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
 * (C)Alexandre Duret-Lutz - 1999-10-14
 * (R)Régis Clouard - 2001-04-03 (version 3.00)
 */

#ifndef __PERRCH__
#define __PERRCH__

#include <memory.h>
#include "pandore.h"

typedef enum { FAILURE = 0, SUCCESS } FS_t;

class Errc;

void Exit(const Errc&);
void Exit(FS_t);
void Exit(Char);
void Exit(FS_t);
void Exit(Uchar);
void Exit(FS_t);
void Exit(Short);
void Exit(FS_t);
void Exit(Ushort);
void Exit(FS_t);
void Exit(Long);
void Exit(FS_t);
void Exit(Ulong);
void Exit(FS_t);
void Exit(Float);
void Exit(FS_t);
void Exit(Double);

/**
 * The <code>Errc</code> defines a type of anything
 * An Errc instance can be from any type (Char, Uchar, Short, Ushort, 
 * Ulong, Long, Foat, Double, Bool) plus the two predefined values
 * SUCCESS and FAILURE).
 * <i>However, <code>Errc </code>cannot be of type int</i>.
 */
class Errc
{
private:
	enum
	{
		FS_RET,
		Char_RET,
		Uchar_RET,
		Short_RET,
		Ushort_RET,
		Long_RET,
		Ulong_RET,
		Float_RET,
		Double_RET
	} ret;

public :
	/**
	 * Create a return value from
	 * the predefined SUCCESS or FAILURE
	 */
	Errc(const FS_t f = SUCCESS)
	{
		ret = FS_RET;
		val.fs = f;
	}

	/**
	 * Creates a return value from the 
	 * given boolean value.
	 */
	Errc(const bool f)
	{
		ret = FS_RET;
		val.fs = f ? SUCCESS : FAILURE;
	}

	/**
	 * Creates a return value from the 
	 * given Char value.
	 */
	Errc(Char l)
	{
		ret = Char_RET;
		val.v_Char = l;
	}

	/**
	 * Creates a return value from the 
	 * given Uchar value.
	 */
	Errc(Uchar l)
	{
		ret = Uchar_RET;
		val.v_Uchar = l;
	}

	/**
	 * Creates a return value from the 
	 * given Short value.
	 */
	Errc(Short l)
	{
		ret = Short_RET;
		val.v_Short = l;
	}

	/**
	 * Creates a return value from the 
	 * given Ushort value.
	 */
	Errc(Ushort l)
	{
		ret = Ushort_RET;
		val.v_Ushort = l;
	}

	/**
	 * Creates a return value from the 
	 * given Long value.
	 */
	Errc(Long l)
	{
		ret = Long_RET;
		val.v_Long = l;
	}

	/**
	 * Creates a return value from the 
	 * given Ulong value.
	 */
	Errc(Ulong l)
	{
		ret = Ulong_RET;
		val.v_Ulong = l;
	}

	/**
	 * Creates a return value from the 
	 * given Float value.
	 */
	Errc(Float l)
	{
		ret = Float_RET;
		val.v_Float = l;
	}

	/**
	 * Creates a return value from the 
	 * given Double value.
	 */
	Errc(Double l)
	{
		ret = Double_RET;
		val.v_Double = l;
	}

	/**
	 * Creates a return value from the
	 * given Errc value
	 */
	Errc(const Errc& o)
	{
		ret = o.ret;
		memcpy(&val, &o.val, sizeof(val));
	}

	Errc& operator=(const Errc& o)
	{
		ret = o.ret;
		memcpy(&val, &o.val, sizeof(val));
		return *this;
	}

	operator bool() { return ret != FS_RET || val.fs == SUCCESS; }
#ifndef _WIN32
      operator FS_t() { return val.fs ; }
      operator Char() { return val.v_Char; }
      operator Uchar() { return val.v_Uchar; }
      operator Short() { return val.v_Short; }
      operator Ushort() { return val.v_Ushort; }
      operator Long() { return val.v_Long; }
      operator Ulong() { return val.v_Ulong; }
      operator Float() { return val.v_Float; }
      operator Double() { return val.v_Double; }
#endif
	void Exit() const;

private:
	union
	{
		FS_t fs;
		Char v_Char;
		Uchar v_Uchar;
		Short v_Short;
		Ushort v_Ushort;
		Long v_Long;
		Ulong v_Ulong;
		Float v_Float;
		Double v_Double;
	} val;
};

#endif // __PERRCH__
