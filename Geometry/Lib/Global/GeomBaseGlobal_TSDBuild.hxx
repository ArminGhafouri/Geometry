#pragma once

#ifndef _GeomBaseGlobal_TSDBuild_Header
#define _GeomBaseGlobal_TSDBuild_Header

# if defined(Geometry_LIB)

// Building Library

#  define Geometry_EXPORT __declspec(dllexport)
# else
// Reading Library
#  define Geometry_EXPORT __declspec(dllimport)
# endif

#endif // !_GeomBaseGlobal_TSDBuild_Header