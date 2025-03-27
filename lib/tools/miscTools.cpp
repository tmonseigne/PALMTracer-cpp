/*******************************************************************************************************************
								Miscellaneous tools
*******************************************************************************************************************/


#include "miscTools.h"
extern FILE* debugStreamALS;

double getAutoThreshold(FltArray& im)
{
	double threshold;
	double oldThreshold;
	double m0, m1; //the average of the 2 parts of the histogram
	double count0, count1, count;
	double s0, s1;
	double sigma;


	long cpt = 30;
	threshold = 0.;
	oldThreshold = 1e10;
	count = 0;
	for (long i = 0; i < (long)im.nbPoints; i++) {
		if (im.data[i] > 0.) {
			count++;
			threshold += im.data[i];
		}
	}
	threshold /= count;

	while ((fabs(oldThreshold - threshold) > 1e-12) && (cpt > 0)) {
		count0 = count1 = 0.;
		m0 = m1 = 0.;
		for (long i = 0; i < (long)im.nbPoints; i++) {
			if (im.data[i] > 0.) {
				if (im.data[i] < threshold) {
					m0 += im.data[i];
					count0++;
				}
				else {
					m1 += im.data[i];
					count1++;
				}
			}
		}

		m0 /= count0;
		m1 /= count1;
		s0 = s1 = 0.;

		long i;
		for (i = 0; i < (long)im.nbPoints; i++) {
			if (im.data[i] > 0.) {
				if (im.data[i] < threshold) { s0 += (im.data[i] - m0) * (im.data[i] - m0); }
				else { s1 += (im.data[i] - m1) * (im.data[i] - m1); }
			}
		}

		s0 = sqrt(s0);
		s1 = sqrt(s1);

		sigma = (s0 + s1) / count;

		oldThreshold = threshold;
		threshold = (m0 + m1) / 2. + sigma * sigma * log(count0 / count1) / (m1 - m0);

		cpt--;
	}


	return (threshold);
}

void printErrno(FILE* file, int errno)
{
	switch (errno) {
		case ECHILD: fprintf(file, "errno: No spawned processes.\n");
			break;
		case EAGAIN: fprintf(file, "errno: No more processes. An attempt to create a new process failed because there are no more process slots, or there is not enough memory, or the maximum nesting level has been reached.\n");
			break;
		case E2BIG: fprintf(file, "errno: Argument list too long.\n");
			break;
		case EACCES: fprintf(file, "errno: Permission denied. The file's permission setting does not allow the specified access. This error signifies that an attempt was made to access a file (or, in some cases, a directory) in a way that is incompatible with the file's attributes.\n");
			fprintf(file, "For example, the error can occur when an attempt is made to read from a file that is not open, to open an existing read-only file for writing, or to open a directory instead of a file. Under MS-DOS operating system versions 3.0 and later, EACCES may also indicate a locking or sharing violation.\n");
			fprintf(file, "The error can also occur in an attempt to rename a file or directory or to remove an existing directory.\n");
			break;
		case EBADF: fprintf(file, "errno: Bad file number. There are two possible causes: 1) The specified file handle is not a valid file-handle value or does not refer to an open file. 2) An attempt was made to write to a file or device opened for read-only access.\n");
			break;
		case EDEADLOCK: fprintf(file, "errno: Resource deadlock would occur. The argument to a math function is not in the domain of the function.\n");
			break;
		case EDOM: fprintf(file, "errno: Math argument.\n");
			break;
		case EEXIST: fprintf(file, "errno: Files exist. An attempt has been made to create a file that already exists. For example, the _O_CREAT and _O_EXCL flags are specified in an _open call, but the named file already exists.\n");
			break;
		case EINVAL: fprintf(file, "errno: Invalid argument. An invalid value was given for one of the arguments to a function. For example, the value given for the origin when positioning a file pointer (by means of a call to fseek) is before the beginning of the file.\n");
			break;
		case EMFILE: fprintf(file, "errno: Too many open files. No more file handles are available, so no more files can be opened.\n");
			break;
		case ENOENT: fprintf(file, "errno: No such file or directory. The specified file or directory does not exist or cannot be found. This message can occur whenever a specified file does not exist or a component of a path does not specify an existing directory.\n");
			break;
		case ENOEXEC: fprintf(file, "errno: Exec format error. An attempt was made to execute a file that is not executable or that has an invalid executable-file format.\n");
			break;
		case ENOMEM: fprintf(file, "errno: Not enough core. Not enough memory is available for the attempted operator. For example, this message can occur when insufficient memory is available to execute a child process, or when the allocation request in a _getcwd call cannot be satisfied.\n");
			break;
		case ENOSPC: fprintf(file, "errno: No space left on device. No more space for writing is available on the device (for example, when the disk is full).\n");
			break;
		case ERANGE: fprintf(file, "errno: Result too large. An argument to a math function is too large, resulting in partial or total loss of significance in the result. This error can also occur in other functions when an argument is larger than expected (for example, when the buffer argument to _getcwd is longer than expected).\n");
			break;
		case EXDEV: fprintf(file, "errno: Cross-device link. An attempt was made to move a file to a different device (using the rename function).\n");
			break;
		default: fprintf(file, "errno: Unknown error.\n");
	}
}

double tenTo(double a) { return (exp(a * log(10.))); }

wchar_t* stringReplace(wchar_t s, wchar_t d, wchar_t* st)
{
	wchar_t* ptrSt = st;
	while (*ptrSt != '\0') {
		if (*ptrSt == s) { *ptrSt = d; }
		ptrSt++;
	}
	return (st);
}


unsigned long twoToThe(unsigned long val)
{
	if (val == 0) { return (1); }
	return (((unsigned long)1) << val);
}


double distance(double X1, double Y1, double X2, double Y2) { return (sqrt(distance2(X1, Y1, X2, Y2))); }
double distance(double X1, double Y1, double Z1, double X2, double Y2, double Z2) { return (sqrt(distance2(X1, Y1, Z1, X2, Y2, Z2))); }

double distance2(double X1, double Y1, double X2, double Y2)
{
	double res;
	res = (X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2);
	if (res < 0) { return (0.); }
	return (res);
}

double distance2(double X1, double Y1, double Z1, double X2, double Y2, double Z2)
{
	double res;
	res = (X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2) + (Z1 - Z2) * (Z1 - Z2);
	if (res < 0) { return (0.); }
	return (res);
}

inline long turnClockwiseX(const long x, const long y)
{
	if (x != 0) { return (0); }
	else if (y == 1) { return (-1); }
	else { return (1); }
}

inline long turnClockwiseY(const long x, const long y)
{
	if (y != 0) { return (0); }
	else if (x == 1) { return (1); }
	else { return (-1); }
}

inline void turnClockwise(long& x, long& y)
{
	if (y == 0) {
		if (x == 1) {
			x = 0;
			y = 1;
		}
		else {
			x = 0;
			y = -1;
		}
	}
	else if (y == 1) {
		x = -1;
		y = 0;
	}
	else {
		x = 1;
		y = 0;
	}
}

inline void turnAntiClockwise(long& x, long& y)
{
	if (y == 0) {
		if (x == 1) {
			x = 0;
			y = -1;
		}
		else {
			x = 0;
			y = 1;
		}
	}
	else if (y == 1) {
		x = 1;
		y = 0;
	}
	else {
		x = -1;
		y = 0;
	}
}


#define ISOUT(x,y) (isOutside((mainSegmentation), (val), (x), (y), (xSize), (ySize)))

bool isOutside(USrtArray& mainSegmentation, const long structure, const long x, const long y, const long xSize, const long ySize)
{
	if (x < 0) { return (true); }
	if (y < 0) { return (true); }
	if (x >= xSize) { return (true); }
	if (y >= ySize) { return (true); }
	return (mainSegmentation.data[x * ySize + y] != structure);
}

#undef ISOUT
