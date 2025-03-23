#include "vor3d/Dexelize.h"
#include "vor3d/AABB.h"

#include <random>
#include <chrono>
#include <algorithm>
#include <array>
#include <iterator>

#ifndef TEST_WID
#define TEST_WID 25
#endif
#ifndef TEST_LEN
#define TEST_LEN 30
#endif
#ifndef TEST_HEI
#define TEST_HEI 40
#endif

/**
 * @brief Compute the signed area of triangle (0,0)-(x1,y1)-(x2,y2), used to determine orientation.
 * @param x1,y1,x2,y2 Coordinates of triangle vertices
 * @param[out] twiceSignedArea The twice signed area result
 * @return -1, 0, or +1 depending on orientation
 */
int Orientation(double x1, double y1, double x2, double y2, double& twiceSignedArea) {
    twiceSignedArea = y1 * x2 - x1 * y2;
    if (twiceSignedArea > 0) return 1;
    else if (twiceSignedArea < 0) return -1;
    else if (y2 > y1) return 1;
    else if (y2 < y1) return -1;
    else if (x1 > x2) return 1;
    else if (x1 < x2) return -1;
    else return 0;
}

/**
 * @brief Robust 2D point-in-triangle test with barycentric coordinate output.
 */
bool PointInTriangle2D(
    double x0, double y0,
    double x1, double y1,
    double x2, double y2,
    double x3, double y3,
    double& a, double& b, double& c)
{
    x1 -= x0; x2 -= x0; x3 -= x0;
    y1 -= y0; y2 -= y0; y3 -= y0;
    int signA = Orientation(x2, y2, x3, y3, a);
    if (signA == 0) return false;
    int signB = Orientation(x3, y3, x1, y1, b);
    if (signB != signA) return false;
    int signC = Orientation(x1, y1, x2, y2, c);
    if (signC != signA) return false;
    double sum = a + b + c;
    a /= sum;
    b /= sum;
    c /= sum;
    return true;
}

/**
 * @brief Intersect a vertical ray with a triangle and compute intersection Z.
 */
bool IntersectRayZ(
    const std::vector<float>& _vertices,
    const std::vector<unsigned int>& _facetIndices,
    size_t faceIndex,
    const Vertex& q,
    double& z)
{
    size_t i0 = _facetIndices[faceIndex * 3 + 0];
    size_t i1 = _facetIndices[faceIndex * 3 + 1];
    size_t i2 = _facetIndices[faceIndex * 3 + 2];

    const float* p0 = &_vertices[i0 * 3];
    const float* p1 = &_vertices[i1 * 3];
    const float* p2 = &_vertices[i2 * 3];

    double u, v, w;
    if (PointInTriangle2D(
        q.x, q.y,
        p0[0], p0[1],
        p1[0], p1[1],
        p2[0], p2[1],
        u, v, w))
    {
        z = u * p0[2] + v * p1[2] + w * p2[2];
        return true;
    }
    return false;
}

/**
 * @brief Compute dexel intersections and fill CompressedVolume structure.
 */
void ComputeDexelIntersections(
    const std::vector<float>& _vertices,
    const std::vector<unsigned int>& _facetIndices,
    const AABB& aabb,
    voroffset3d::CompressedVolume& dexels)
{
    auto gridSize = dexels.gridSize();
	auto spacing = dexels.spacing();

    Vertex bbMin, bbMax;
    if (!aabb.Boxes().empty()) {
        bbMin = aabb.Boxes()[1].minCorner;
        bbMax = aabb.Boxes()[1].maxCorner;
    }

    for (int y = 0; y < gridSize[1]; ++y) {
        for (int x = 0; x < gridSize[0]; ++x) {
			
			Eigen::Vector2d center2d = dexels.dexelCenter(x, y);
			Vertex center(center2d.x(), center2d.y(), 0.0f);  // z = 0

            Box queryBox;
            queryBox.minCorner = Vertex(center.x, center.y, bbMin.z - spacing);
            queryBox.maxCorner = Vertex(center.x, center.y, bbMax.z + spacing);

            std::vector<double> intersections;
            aabb.Traverse(queryBox, [&](size_t faceIndex) {
                double z;
                if (IntersectRayZ(_vertices, _facetIndices, faceIndex, center, z)) {
                    intersections.push_back(z / spacing);
                }
            });

            std::sort(intersections.begin(), intersections.end());

			// Write intersections directly into CompressedVolume (same as original compute_sign)
			dexels.at(x, y).resize(intersections.size());
			std::copy(intersections.begin(), intersections.end(), dexels.at(x, y).begin());
        }
    }
}

voroffset3d::CompressedVolume CreateDexelsFromMeshBuffers(
	const std::vector<float>& _vertices,
	const std::vector<unsigned int>& _facetIndices,
	double &_voxelSize,
	int _padding, 
	int _numVoxels
) {

	auto _ComputeBoundingBox = [&](Vertex& bbMin, Vertex& bbMax) {

		if (_vertices.empty()) {
			return;
		}
	
		bbMin = Vertex( std::numeric_limits<float>::max(),
						std::numeric_limits<float>::max(),
						std::numeric_limits<float>::max() );
	
		bbMax = Vertex( std::numeric_limits<float>::lowest(),
						std::numeric_limits<float>::lowest(),
						std::numeric_limits<float>::lowest() );
	
		auto n = _vertices.size() / 3;

		for (auto i = 0; i < n; ++i) {
			
			float x = _vertices[3 * i];
			float y = _vertices[3 * i + 1];
			float z = _vertices[3 * i + 2];
	
			bbMin.x = std::min(bbMin.x, x);
			bbMin.y = std::min(bbMin.y, y);
			bbMin.z = std::min(bbMin.z, z);
	
			bbMax.x = std::max(bbMax.x, x);
			bbMax.y = std::max(bbMax.y, y);
			bbMax.z = std::max(bbMax.z, z);
		}
	};
	
	// Initialize voxel grid and AABB tree
	
	Vertex minCorner, maxCorner;
	_ComputeBoundingBox(minCorner, maxCorner);
	Vertex extent = maxCorner - minCorner;

	if (_numVoxels > 0) {
		// Force number of voxels along longest axis
		double maxExtent = std::max(extent.x, std::max(extent.y, extent.z));
		_voxelSize = maxExtent / _numVoxels;
	}
	
	AABB aabbTree;

	aabbTree.Initialize(
		_facetIndices.size() / 3,  // 一個 triangle 有 3 個 index
		[&](Box& box, size_t faceIndex) {
			size_t i0 = _facetIndices[faceIndex * 3 + 0];
			size_t i1 = _facetIndices[faceIndex * 3 + 1];
			size_t i2 = _facetIndices[faceIndex * 3 + 2];
	
			const float* p0 = &_vertices[i0 * 3];
			const float* p1 = &_vertices[i1 * 3];
			const float* p2 = &_vertices[i2 * 3];
	
			box.minCorner = Vertex(p0[0], p0[1], p0[2]);
			box.maxCorner = box.minCorner;
	
			auto UpdateBounds = [&](const float* p) {
				box.minCorner.x = std::min(box.minCorner.x, p[0]);
				box.minCorner.y = std::min(box.minCorner.y, p[1]);
				box.minCorner.z = std::min(box.minCorner.z, p[2]);
	
				box.maxCorner.x = std::max(box.maxCorner.x, p[0]);
				box.maxCorner.y = std::max(box.maxCorner.y, p[1]);
				box.maxCorner.z = std::max(box.maxCorner.z, p[2]);
			};
	
			UpdateBounds(p1);
			UpdateBounds(p2);
		}
	);
	
	// Dexelize the input mesh
	voroffset3d::CompressedVolume dexels(
		Eigen::Vector3d(minCorner.x, minCorner.y, minCorner.z),
		Eigen::Vector3d(extent.x, extent.y, extent.z),
		_voxelSize, _padding);
	
	ComputeDexelIntersections(_vertices, _facetIndices, aabbTree, dexels);
	return dexels;
}