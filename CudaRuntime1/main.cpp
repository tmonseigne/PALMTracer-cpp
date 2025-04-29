#include <windows.h>

extern "C" void LaunchLoadKernel(int durationSeconds);

_declspec(dllexport)
long _stdcall Process() {
    LaunchLoadKernel(10); // Charge le GPU pendant 10 secondes
    return 0;
}