

/*******************************************************************************************************************
								Define operators for all complex data types
									Institut Curie
									UMR - 144
									by Victor Racine
									sept, 24th 2003
*******************************************************************************************************************/

#ifdef ComplexNumber2


ComplexNumber& ComplexNumber::operator=(const ComplexNumber2 &complexNumber){
	im=(DataTypeComplex)complexNumber.im;
	re=(DataTypeComplex)complexNumber.re;
	return(*this);
}


#endif //ComplexNumber2