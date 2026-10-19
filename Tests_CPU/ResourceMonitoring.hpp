#pragma once

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <optional>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

// Dépend de Windows.h : conserver cette inclusion après les définitions des types Win32.
#include <psapi.h>
#else
#include <sys/resource.h>
#include <unistd.h>
#if defined(__APPLE__)
#include <mach/mach.h>
#endif
#endif

namespace palm_tracer::monitoring {
/// <summary>Relevé ponctuel des ressources du processus C++, tous ses threads compris, sans ses processus enfants.</summary>
struct Snapshot {
	std::chrono::steady_clock::time_point timestamp{}; ///< Instant monotone associé au compteur CPU.
	std::optional<double> cpu_seconds;                 ///< Temps CPU cumulé utilisateur et noyau, en secondes.
	std::optional<std::uint64_t> resident_bytes;       ///< Mémoire résidente courante ; absence si la lecture échoue.
};

/// <summary>Renvoie le nombre de processeurs logiques de la machine, ou zéro si cette information est indisponible.</summary>
inline unsigned int LogicalProcessorCount()
{
	#if defined(_WIN32)
	return GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);
	#else
	const long Count = sysconf(_SC_NPROCESSORS_ONLN);
	return Count > 0 ? static_cast<unsigned int>(Count) : 0U;
	#endif
}

/// <summary>Lit les compteurs natifs sans attente, allocation de tampon de calcul ni thread de surveillance.</summary>
inline Snapshot Capture()
{
	Snapshot result;
	#if defined(_WIN32)
	FILETIME creation{}, exit{}, kernel{}, user{};
	if (GetProcessTimes(GetCurrentProcess(), &creation, &exit, &kernel, &user)) {
		const auto Kernel  = (static_cast<std::uint64_t>(kernel.dwHighDateTime) << 32U) | kernel.dwLowDateTime;
		const auto User    = (static_cast<std::uint64_t>(user.dwHighDateTime) << 32U) | user.dwLowDateTime;
		result.cpu_seconds = (static_cast<double>(Kernel) + static_cast<double>(User)) * 1e-7;
	}
	#else
	rusage usage{};
	if (getrusage(RUSAGE_SELF, &usage) == 0) {
		result.cpu_seconds = static_cast<double>(usage.ru_utime.tv_sec) + static_cast<double>(usage.ru_stime.tv_sec)
		                     + (static_cast<double>(usage.ru_utime.tv_usec) + static_cast<double>(usage.ru_stime.tv_usec)) * 1e-6;
	}
	#endif
	result.timestamp = std::chrono::steady_clock::now();
	#if defined(_WIN32)
	PROCESS_MEMORY_COUNTERS memory{};
	if (K32GetProcessMemoryInfo(GetCurrentProcess(), &memory, sizeof(memory))) { result.resident_bytes = memory.WorkingSetSize; }
	#elif defined(__APPLE__)
	mach_task_basic_info_data_t memory{};
	mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
	if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&memory), &count) == KERN_SUCCESS) {
		result.resident_bytes = memory.resident_size;
	}
	#else
	std::ifstream stream("/proc/self/statm");
	std::uint64_t virtual_pages = 0, resident_pages = 0;
	const long PageSize         = sysconf(_SC_PAGESIZE);
	if (PageSize > 0 && stream >> virtual_pages >> resident_pages) { result.resident_bytes = resident_pages * static_cast<std::uint64_t>(PageSize); }
	#endif
	return result;
}

/// <summary>Calcule la part moyenne de la capacité CPU totale occupée par le processus, en pourcentage entre 0 et 100.</summary>
/// <remarks>Formule : 100 × temps CPU / (durée écoulée × processeurs logiques). Une information absente reste indisponible.
/// Les tests très courts sont sensibles à la résolution des compteurs CPU du système.</remarks>
inline std::optional<double> CpuPercent(const Snapshot& start, const Snapshot& end, const unsigned int processors)
{
	const double Elapsed = std::chrono::duration<double>(end.timestamp - start.timestamp).count();
	if (!start.cpu_seconds || !end.cpu_seconds || processors == 0 || Elapsed <= 0 || *end.cpu_seconds < *start.cpu_seconds) { return std::nullopt; }
	return std::clamp(100.0 * (*end.cpu_seconds - *start.cpu_seconds) / (Elapsed * processors), 0.0, 100.0);
}
} // namespace palm_tracer::monitoring
