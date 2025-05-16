///-------------------------------------------------------------------------------------------------
/// 
/// \file Constants.hpp
/// \brief  Fichier de constantes necessaire aux algorithmes.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
/// 
///-------------------------------------------------------------------------------------------------

#pragma once

//=================================================
//==================== Threads ====================
//=================================================
constexpr int N_THREADS = 4;


//================================================
//==================== Errors ====================
//================================================
// OK
constexpr int RET_OK = 0;

// MIA
constexpr int SIZE_ERR = -1;
constexpr int ALLOC_ERR = -2;

// Wavelet Basic Static
constexpr int W_SIZE_ERR = -101;

// Segmentation Static
constexpr int S_ALLOC_ERR = -201;


// Tracking
constexpr int T_ERR_OPENNING_FILE = -801;
constexpr int T_ERR_NBDIM = -805;
constexpr int T_ERR_ALLOC = -806;
constexpr int T_ERR_CORRUPT_FILE = -807;
constexpr int T_ERR_BAD_MODEL = -808;
constexpr int T_ERR_OOP = -809;
