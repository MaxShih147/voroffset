#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "vor3d/Common.h"
#include "vor3d/CompressedVolume.h"
#include <set>
#include <limits>
#include <cassert>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace voroffset3d
{
/**
 * @brief Deprecated: Creates dexel volume from a model file.
 * @note This modification was implemented by max_shih.
 *
 * This function uses geogram (under GEO namespace) routines to generate dexels from a file.
 * Now models are pre-parsed into vertices and facet indices externally, and GEO functions
 * are being replaced by custom code. Use the new (name) function instead.
 *
 */
/**
 * @brief         Creates dexels from a triangle mesh.
 *
 * @param[in]     filename    { Filename of the mesh to load. }
 * @param[in,out] voxel_size  { Voxel size for the 2D dexel grid. }
 * @param[in]     padding     { Additional padding on each side. }
 * @param[in]     num_voxels  { Explicitly set the 2D grid size (max length), before padding. }
 *
 * @return        { The dexelized volume. }
 */
// CompressedVolume create_dexels(const std::string &filename,
// 	double &voxel_size, int padding = 0, int num_voxels = -1);

/**
 * @brief Dumps dexel structure as a mesh file.
 * @note This modification was implemented by max_shih.
 * 
 * This function exports the dexel representation into a complete mesh and writes it to the specified file path.
 * However, for our current requirements, we only need to export data structures suitable for Three.js,
 * namely a vertices array and a facet indices array. Therefore, this function is deprecated,
 * and has been replaced by the new (need a name...) function.
 * 
 */	
// void dexel_dump(const std::string &filename, const CompressedVolume &voxels);

} // namespace voroffset3d
