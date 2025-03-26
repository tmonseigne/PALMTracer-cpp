

/*******************************************************************************************************************
								Define operators for all complex data types
									Institut Curie
									UMR - 144
									by Victor Racine
									sept, 24th 2003
*******************************************************************************************************************/

#ifdef DataTypeComplex2


ComplexNumber& ComplexNumber::operator=(const DataTypeComplex2 &dataTypeComplex2){
	im=(DataTypeComplex)0;
	re=(DataTypeComplex)dataTypeComplex2;
	return(*this);
}


#endif //DataTypeComplex2