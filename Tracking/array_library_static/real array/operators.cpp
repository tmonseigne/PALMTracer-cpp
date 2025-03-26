#include "stdafx.h"

/*******************************************************************************************************************
								All =, +=, -=, *=, /= Operators 
*******************************************************************************************************************/


#ifdef ArrayType2




//**************		=		*********************
ArrayType& ArrayType::operator=(const ArrayType2 &arrayType){
	if(nbPoints!=arrayType.nbPoints){
			openDebugStreamALS();
			fprintf(debugStreamALS, "££ £ SIZE_ERR in operator=(const ArrayType2 0x%p &arrayType), this:0x%p £ ££\n",&arrayType, this);
			fprintf(debugStreamALS, "££ £ nbPoints: %d <> %d\n",nbPoints, arrayType.nbPoints);
			closeDebugStreamALS();
			exit(SIZE_ERR);
	}
	if(type()==arrayType.type())
		memcpy(data, arrayType.data, sizeof(DataType)*nbPoints);
	else
		for(unsigned long i=0; i<nbPoints; i++){
			#ifdef BOOL_TYPE
			if(arrayType.data[i]==(DataType2)0)
				data[i]=false;
			else
				data[i]=true;
			#else
				#ifdef COMPLEX
					#ifdef COMPLEX2
						data[i][0]=(DataTypeSimple)arrayType.data[i][0];
						data[i][1]=(DataTypeSimple)arrayType.data[i][1];
					#else
						data[i][0]=(DataTypeSimple)arrayType.data[i];
						data[i][1]=(DataTypeSimple)0;
					#endif
				#else
					data[i]=(DataType)arrayType.data[i];
				#endif
			#endif
		}
	return(*this);
}


ArrayType& ArrayType::operator=(const DataType2 val){
	DataType val2;
	#ifdef BOOL_TYPE
	if(val==0) val2=false;
	else val2=true;
	#else
		#ifdef COMPLEX
			#ifdef COMPLEX2
				val2[0]=(DataTypeSimple)val[0];
				val2[1]=(DataTypeSimple)val[1];
			#else
				val2[0]=(DataTypeSimple)val;
				val2[1]=(DataTypeSimple)0;
			#endif
		#else
			val2=(DataType)val;
		#endif
	#endif
	for(unsigned long i=0; i<nbPoints; i++){
		#ifdef COMPLEX
				data[i][0]=val2[0];
				data[i][1]=val2[1];
		#else
			data[i]=val2;
		#endif
	}
	return(*this);
}



//**************		+=		*********************
#ifndef BOOL_TYPE
#ifndef BOOL_TYPE1
ArrayType& ArrayType::operator+=(const ArrayType2 &arrayType){
	if(nbPoints!=arrayType.nbPoints){
		openDebugStreamALS();
		fprintf(debugStreamALS, "££ £ SIZE_ERR in operator+=(const ArrayType2 0x%p &arrayType), this:0x%p £ ££\n",&arrayType, this);
		closeDebugStreamALS();
		exit(SIZE_ERR);
	}
	for(unsigned long i=0; i<nbPoints; i++){
		#ifdef BOOL_TYPE
		if (arrayType.data[i]!=(DataType2)0)
			data[i]=true;
		#else
			#ifdef COMPLEX
				#ifdef COMPLEX2
					data[i][0]+=(DataTypeSimple)arrayType.data[i][0];
					data[i][1]+=(DataTypeSimple)arrayType.data[i][1];
				#else
					data[i][0]+=(DataTypeSimple)arrayType.data[i];
				#endif
			#else
				data[i]+=(DataType)arrayType.data[i];
			#endif	
		#endif
	}
	return(*this);
}



ArrayType& ArrayType::operator+=(const DataType2 val){
	DataType val2;
	#ifdef COMPLEX
		#ifdef COMPLEX2
			val2[0]=(DataTypeSimple)val[0];
			val2[1]=(DataTypeSimple)val[1];
		#else
			val2[0]=(DataTypeSimple)val;
			val2[1]=(DataTypeSimple)0;
		#endif
	#else
		val2=(DataType)val;
	#endif
	for(unsigned long i=0; i<nbPoints; i++){
	#ifdef COMPLEX
		#ifdef COMPLEX2
			data[i][0]+=val2[0];
			data[i][1]+=val2[1];
		#else
			data[i][0]+=val2[0];
		#endif
	#else
		data[i]+=val2;
	#endif
	}
		
	return(*this);
}
#endif
#endif

//**************		-=		*********************

#ifndef BOOL_TYPE
#ifndef BOOL_TYPE1
ArrayType& ArrayType::operator-=(const ArrayType2 &arrayType){
	if(nbPoints!=arrayType.nbPoints){
		openDebugStreamALS();
		fprintf(debugStreamALS, "££ £ SIZE_ERR in operator-=(const ArrayType2 0x%p &arrayType), this:0x%p £ ££\n",&arrayType, this);
		closeDebugStreamALS();
		exit(SIZE_ERR);
	}
	for(unsigned long i=0; i<nbPoints; i++){
		#ifdef COMPLEX
			#ifdef COMPLEX2
				data[i][0]-=(DataTypeSimple)arrayType.data[i][0];
				data[i][1]-=(DataTypeSimple)arrayType.data[i][1];
			#else
				data[i][0]-=(DataTypeSimple)arrayType.data[i];
			#endif
		#else
			data[i]-=(DataType)arrayType.data[i];
		#endif	
	}
	return(*this);
}

ArrayType& ArrayType::operator-=(const DataType2 val){
	DataType val2;
	#ifdef COMPLEX
		#ifdef COMPLEX2
			val2[0]=(DataTypeSimple)val[0];
			val2[1]=(DataTypeSimple)val[1];
		#else
			val2[0]=(DataTypeSimple)val;
			val2[1]=(DataTypeSimple)0;
		#endif
	#else
		val2=(DataType)val;
	#endif

	
	for(unsigned long i=0; i<nbPoints; i++){
	#ifdef COMPLEX
		#ifdef COMPLEX2
			data[i][0]-=val2[0];
			data[i][1]-=val2[1];
		#else
			data[i][0]-=val2[0];
		#endif
	#else
		data[i]-=val2;
	#endif
	}
	return(*this);
}


//**************		*=		*********************
ArrayType& ArrayType::operator*=(const ArrayType2 &arrayType){
	if(nbPoints!=arrayType.nbPoints){
		openDebugStreamALS();
		fprintf(debugStreamALS, "££ £ SIZE_ERR in operator*=(const ArrayType2 0x%p &arrayType), this:0x%p £ ££\n",&arrayType, this);
		closeDebugStreamALS();
		exit(SIZE_ERR);
	}
	for(unsigned long i=0; i<nbPoints; i++){
		#ifdef BOOL_TYPE
			if(arrayType.data[i]==(DataType2)0)
				data[i]=false;
		#else
			#ifdef COMPLEX
				#ifdef COMPLEX2
					DataTypeSimple data_i__0_=data[i][0];
					data[i][0]=data_i__0_*(DataTypeSimple)arrayType.data[i][0]-data[i][1]*(DataTypeSimple)arrayType.data[i][1];
					data[i][1]=data[i][1]*(DataTypeSimple)arrayType.data[i][0]+data_i__0_*(DataTypeSimple)arrayType.data[i][1];
				#else
					data[i][0]*=(DataTypeSimple)arrayType.data[i];
					data[i][1]*=(DataTypeSimple)arrayType.data[i];
				#endif
			#else
				data[i]*=(DataType)arrayType.data[i];
			#endif			
		#endif
	}
	return(*this);
}
#endif
#endif

#ifndef BOOL_TYPE
#ifndef BOOL_TYPE1
ArrayType& ArrayType::operator*=(const DataType2 val){
	DataType val2;
	#ifdef COMPLEX
		#ifdef COMPLEX2
			val2[0]=(DataTypeSimple)val[0];
			val2[1]=(DataTypeSimple)val[1];
		#else
			val2[0]=(DataTypeSimple)val;
			val2[1]=(DataTypeSimple)0;
		#endif
	#else
		val2=(DataType)val;
	#endif
	
	for(unsigned long i=0; i<nbPoints; i++){
		#ifdef COMPLEX
			#ifdef COMPLEX2
				DataTypeSimple data_i__0_=data[i][0];
				data[i][0]=data_i__0_*val2[0]-data[i][1]*val2[1];
				data[i][1]=data[i][1]*val2[0]+data_i__0_*val2[1];
			#else
				data[i][0]*=val2[0];
				data[i][1]*=val2[0];
			#endif
		#else
			data[i]*=val2;
		#endif
	}
	return(*this);
}
#endif
#endif

//**************		/=		*********************

#ifndef BOOL_TYPE
#ifndef BOOL_TYPE1
ArrayType& ArrayType::operator/=(const ArrayType2 &arrayType){
	if(nbPoints!=arrayType.nbPoints){
		openDebugStreamALS();
		fprintf(debugStreamALS, "££ £ SIZE_ERR in operator/=(const ArrayType2 0x%p &arrayType), this:0x%p £ ££\n",&arrayType, this);
		closeDebugStreamALS();
		exit(SIZE_ERR);
	}
	for(unsigned long i=0; i<nbPoints; i++){
		#ifdef COMPLEX
			#ifdef COMPLEX2
				DataTypeSimple data_i__0_=data[i][0];
				data[i][0]=(data_i__0_*(DataTypeSimple)arrayType.data[i][0]+data[i][1]*(DataTypeSimple)arrayType.data[i][1])
					/((DataTypeSimple)arrayType.data[i][0]*(DataTypeSimple)arrayType.data[i][0]
					+(DataTypeSimple)arrayType.data[i][1]*(DataTypeSimple)arrayType.data[i][1]);
				data[i][1]=(data[i][1]*(DataTypeSimple)arrayType.data[i][0]-data_i__0_*(DataTypeSimple)arrayType.data[i][1])
					/((DataTypeSimple)arrayType.data[i][0]*(DataTypeSimple)arrayType.data[i][0]
					+(DataTypeSimple)arrayType.data[i][1]*(DataTypeSimple)arrayType.data[i][1]);
			#else
				data[i][0]/=(DataTypeSimple)arrayType.data[i];
				data[i][1]/=(DataTypeSimple)arrayType.data[i];
			#endif
		#else
			data[i]/=(DataType)arrayType.data[i];
		#endif
	}
	return(*this);
}


ArrayType& ArrayType::operator/=(const DataType2 val){
	DataType val2;
	#ifdef COMPLEX
		#ifdef COMPLEX2
			val2[0]=(DataTypeSimple)val[0];
			val2[1]=(DataTypeSimple)val[1];
		#else
			val2[0]=(DataTypeSimple)val;
			val2[1]=(DataTypeSimple)0;
		#endif
	#else
		val2=(DataType)val;
	#endif
	for(unsigned long i=0; i<nbPoints; i++){
		#ifdef COMPLEX
			#ifdef COMPLEX2
				DataTypeSimple data_i__0_=data[i][0];
				data[i][0]=(data_i__0_*val2[0]+data[i][1]*val2[1])
					/(val2[0]*val2[0]+val2[1]*val2[1]);
				data[i][1]=(data[i][1]*val2[0]-data_i__0_*val2[1])
					/(val2[0]*val2[0]+val2[1]*val2[1]);
			#else
				data[i][0]/=val2[0];
				data[i][1]/=val2[0];
			#endif
		#else
			data[i]/=val2;;
		#endif
	}
	return(*this);
}
#endif
#endif

#endif