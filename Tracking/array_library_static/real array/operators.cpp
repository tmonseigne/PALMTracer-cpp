/*******************************************************************************************************************
								All =, +=, -=, *=, /= Operators 
*******************************************************************************************************************/


#ifdef ArrayType2


//**************		=		*********************
ArrayType& ArrayType::operator=(const ArrayType2& arrayType)
{
	if (nbPoints != arrayType.nbPoints) {
		exit(SIZE_ERR);
	}
	if (type() == arrayType.type()) memcpy(data, arrayType.data, sizeof(DataType) * nbPoints);
	else
		for (unsigned long i = 0; i < nbPoints; i++) {
#ifdef BOOL_TYPE
			if(arrayType.data[i]==(DataType2)0)
				data[i]=false;
			else
				data[i]=true;
#else
			data[i] = (DataType)arrayType.data[i];
#endif
		}
	return (*this);
}


ArrayType& ArrayType::operator=(const DataType2 val)
{
	DataType val2;
#ifdef BOOL_TYPE
	if(val==0) val2=false;
	else val2=true;
#else
	val2 = (DataType)val;
#endif
	for (unsigned long i = 0; i < nbPoints; i++) {
		data[i] = val2;
	}
	return (*this);
}


//**************		+=		*********************
#ifndef BOOL_TYPE
#ifndef BOOL_TYPE1
ArrayType& ArrayType::operator+=(const ArrayType2& arrayType)
{
	if (nbPoints != arrayType.nbPoints) {
		exit(SIZE_ERR);
	}
	for (unsigned long i = 0; i < nbPoints; i++) {
#ifdef BOOL_TYPE
		if (arrayType.data[i]!=(DataType2)0)
			data[i]=true;
#else
		data[i] += (DataType)arrayType.data[i];
#endif
	}
	return (*this);
}


ArrayType& ArrayType::operator+=(const DataType2 val)
{
	DataType val2;
	val2 = (DataType)val;
	for (unsigned long i = 0; i < nbPoints; i++) {
		data[i] += val2;
	}

	return (*this);
}
#endif
#endif

//**************		-=		*********************

#ifndef BOOL_TYPE
#ifndef BOOL_TYPE1
ArrayType& ArrayType::operator-=(const ArrayType2& arrayType)
{
	if (nbPoints != arrayType.nbPoints) {
		exit(SIZE_ERR);
	}
	for (unsigned long i = 0; i < nbPoints; i++) {
		data[i] -= (DataType)arrayType.data[i];
	}
	return (*this);
}

ArrayType& ArrayType::operator-=(const DataType2 val)
{
	DataType val2;
	val2 = (DataType)val;


	for (unsigned long i = 0; i < nbPoints; i++) {
		data[i] -= val2;
	}
	return (*this);
}


//**************		*=		*********************
ArrayType& ArrayType::operator*=(const ArrayType2& arrayType)
{
	if (nbPoints != arrayType.nbPoints) {
		exit(SIZE_ERR);
	}
	for (unsigned long i = 0; i < nbPoints; i++) {
#ifdef BOOL_TYPE
			if(arrayType.data[i]==(DataType2)0)
				data[i]=false;
#else
		data[i] *= (DataType)arrayType.data[i];
#endif
	}
	return (*this);
}
#endif
#endif

#ifndef BOOL_TYPE
#ifndef BOOL_TYPE1
ArrayType& ArrayType::operator*=(const DataType2 val)
{
	DataType val2;
	val2 = (DataType)val;

	for (unsigned long i = 0; i < nbPoints; i++) {
		data[i] *= val2;
	}
	return (*this);
}
#endif
#endif

//**************		/=		*********************

#ifndef BOOL_TYPE
#ifndef BOOL_TYPE1
ArrayType& ArrayType::operator/=(const ArrayType2& arrayType)
{
	if (nbPoints != arrayType.nbPoints) {
		exit(SIZE_ERR);
	}
	for (unsigned long i = 0; i < nbPoints; i++) {
		data[i] /= (DataType)arrayType.data[i];
	}
	return (*this);
}


ArrayType& ArrayType::operator/=(const DataType2 val)
{
	DataType val2;
	val2 = (DataType)val;
	for (unsigned long i = 0; i < nbPoints; i++) {
		data[i] /= val2;
	}
	return (*this);
}
#endif
#endif

#endif
