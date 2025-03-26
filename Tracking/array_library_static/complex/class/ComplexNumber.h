


/*******************************************************************************************************************
								Pactage for complex numbers
									Institut Curie
									UMR - 144
									by Victor Racine
									sept, 23rd 2003
*******************************************************************************************************************/


#ifdef _COMPLEX_NUMBER_

//#ifndef		_Complex_Numbers_
//#define		_Complex_Numbers_

//#define DataTypeComplex	float

#define		DIV_ZERO	101



#include	<math.h>
#include	<iostream.h>

/*#ifdef ComplexNumber
#error ComplexNumber
#endif*/ 

class ComplexNumber{
public:
	DataTypeComplex	im;
	DataTypeComplex	re;
	ComplexNumber();
	ComplexNumber(DataTypeComplex _im, DataTypeComplex _re);
	~ComplexNumber();

	

	DataTypeComplex	getAngle(void);
	DataTypeComplex	getAbs(void);

	//ComplexNumber& operator=(const ComplexNumber &complexNumber);
	ComplexNumber& operator+=(const ComplexNumber &complexNumber);
	ComplexNumber& operator-=(const ComplexNumber &complexNumber);
	ComplexNumber& operator*=(const ComplexNumber &complexNumber);
	ComplexNumber& operator/=(const ComplexNumber &complexNumber);

	ComplexNumber& operator=(char &dataTypeComplex2);
	ComplexNumber& operator=(const char &dataTypeComplex2);
	ComplexNumber& operator=(const unsigned char &dataTypeComplex2);
	ComplexNumber& operator=(const short &dataTypeComplex2);
	ComplexNumber& operator=(const unsigned short &dataTypeComplex2);
	ComplexNumber& operator=(const long &dataTypeComplex2);
	ComplexNumber& operator=(const unsigned long &dataTypeComplex2);
	ComplexNumber& operator=(const float &dataTypeComplex2);
	ComplexNumber& operator=(const double &dataTypeComplex2);
	ComplexNumber& operator=(const bool &dataTypeComplex2);

	ComplexNumber& operator=(const ComplexChr &complexNumber);
	ComplexNumber& operator=(const ComplexSrt &complexNumber);
	ComplexNumber& operator=(const ComplexLng &complexNumber);
	ComplexNumber& operator=(const ComplexFlt &complexNumber);
	ComplexNumber& operator=(const ComplexDbl &complexNumber);



	ComplexNumber& operator++(void){};
	ComplexNumber& operator--(void){};

	
	ComplexNumber& operator+=(const DataTypeComplex &dataTypeComplex);
	ComplexNumber& operator-=(const DataTypeComplex &dataTypeComplex);
	ComplexNumber& operator*=(const DataTypeComplex &dataTypeComplex);
	ComplexNumber& operator/=(const DataTypeComplex &dataTypeComplex);

	

	


	friend ostream& operator<< ( ostream& os, ComplexNumber &complexNumber);
};



//#endif	//	_Complex_Numbers_

#endif //_COMPLEX_NUMBER_