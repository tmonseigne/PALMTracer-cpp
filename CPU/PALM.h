// Le bloc ifdef suivant est la façon standard de créer des macros qui facilitent l'exportation 
// à partir d'une DLL. Tous les fichiers contenus dans cette DLL sont compilés avec le symbole PALM_EXPORTS
// défini sur la ligne de commande. Ce symbole ne doit pas être défini dans les projets
// qui utilisent cette DLL. De cette manière, les autres projets dont les fichiers sources comprennent ce fichier considèrent les fonctions 
// PALM_API comme étant importées à partir d'une DLL, tandis que cette DLL considère les symboles
// définis avec cette macro comme étant exportés.
//#ifdef PALM_EXPORTS
//#define PALM_API __declspec(dllexport)
//#else
//#define PALM_API __declspec(dllimport)
//#endif




//extern PALM_API int nPALM;
//PALM_API int fnPALM(void);
//PALM_API void OpenPALMProcessing(unsigned short * image,  double * pointsList, unsigned long maxNbPoints, unsigned long sizeX, unsigned long sizeY, 
//							unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin);
_declspec(dllexport) void _stdcall _OpenPALMProcessing(unsigned short * image,  double * pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY, 
							unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX,
							double sigmaY, double angle, unsigned short size);
//PALM_API long PALMProcessing();
//PALM_API void closePALMProcessing();
_declspec(dllexport) long _stdcall _PALMProcessing();
_declspec(dllexport) void _stdcall _closePALMProcessing();