// Le bloc ifdef suivant est la façon standard de créer des macros qui facilitent l'exportation 
// à partir d'une DLL. Tous les fichiers contenus dans cette DLL sont compilés avec le symbole PALM_EXPORTS
// défini sur la ligne de commande. Ce symbole ne doit pas étre défini dans les projets
// qui utilisent cette DLL. De cette manière, les autres projets dont les fichiers sources comprennent ce fichier considérent les fonctions 
// PALM_API comme étant importées à partir d'une DLL, tandis que cette DLL considére les symboles
// définis avec cette macro comme étant exportés.

_declspec(dllexport) void _stdcall _OpenPALMProcessing(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size);
_declspec(dllexport) long _stdcall _PALMProcessing();
_declspec(dllexport) void _stdcall _closePALMProcessing();
_declspec(dllexport) long _stdcall Process(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size);
