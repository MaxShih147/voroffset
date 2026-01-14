#include "vor3d/Dexelize.h"
#include "vor3d/AABB.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <iterator>
#include <random>

#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <cstdint>

#ifndef TEST_WID
#define TEST_WID 25
#endif
#ifndef TEST_LEN
#define TEST_LEN 30
#endif
#ifndef TEST_HEI
#define TEST_HEI 40
#endif

void laplacian_smooth(std::vector<float>& _vertices, const std::vector<unsigned int>& _facetIndices) {
    struct v_conncted {
        std::vector<int> facet_id;
        std::vector<unsigned int> id_in_f;
    };

    std::vector<float> orig_vertices = _vertices;

    if (_vertices.size() % 3 != 0 || _facetIndices.size() % 3 != 0) // input error
        return;

    const int vertex_size_orig = static_cast<int>(_vertices.size() / 3);
    const int facet_size = static_cast<int>(_facetIndices.size() / 3);

    if (vertex_size_orig == 0 || facet_size == 0)
        return;

    float weld_eps = 0.01f;

    // smoothing strength: 1.0 = original behavior, >1.0 stronger, <1.0 weaker
    float lambda = 1.0f;
    int iterations = 5;

    float weld_eps2 = weld_eps * weld_eps;

    struct QuantKey {
        int ix, iy, iz;
    };
    struct QuantKeyHash {
        std::size_t operator()(const QuantKey& k) const noexcept {
            std::size_t h = static_cast<std::size_t>(k.ix);
            h ^= static_cast<std::size_t>(k.iy) + 0x9e3779b9u + (h << 6) + (h >> 2);
            h ^= static_cast<std::size_t>(k.iz) + 0x9e3779b9u + (h << 6) + (h >> 2);
            return h;
        }
    };
    struct QuantKeyEq {
        bool operator()(const QuantKey& a, const QuantKey& b) const noexcept {
            return a.ix == b.ix && a.iy == b.iy && a.iz == b.iz;
        }
    };

    std::unordered_map<QuantKey, std::vector<unsigned int>, QuantKeyHash, QuantKeyEq> cells;
    cells.reserve(static_cast<std::size_t>(vertex_size_orig));

    std::vector<unsigned int> rep(static_cast<std::size_t>(vertex_size_orig));

    for (int i = 0; i < vertex_size_orig; ++i) {
        float x = _vertices[i * 3 + 0];
        float y = _vertices[i * 3 + 1];
        float z = _vertices[i * 3 + 2];

        int ix = static_cast<int>(std::floor(x / weld_eps + 0.5f));
        int iy = static_cast<int>(std::floor(y / weld_eps + 0.5f));
        int iz = static_cast<int>(std::floor(z / weld_eps + 0.5f));

        QuantKey key{ ix, iy, iz };
        auto& bucket = cells[key];

        bool found = false;
        for (unsigned int idx : bucket) {
            float x2 = _vertices[idx * 3 + 0];
            float y2 = _vertices[idx * 3 + 1];
            float z2 = _vertices[idx * 3 + 2];
            float ddx = x - x2;
            float ddy = y - y2;
            float ddz = z - z2;
            float dist2 = ddx * ddx + ddy * ddy + ddz * ddz;
            if (dist2 <= weld_eps2) {
                rep[static_cast<std::size_t>(i)] = idx;
                found = true;
                break;
            }
        }

        if (!found) {
            unsigned int iu = static_cast<unsigned int>(i);
            bucket.push_back(iu);
            rep[static_cast<std::size_t>(i)] = iu;
        }
    }

    std::vector<unsigned int> rep_to_new(static_cast<std::size_t>(vertex_size_orig),
        std::numeric_limits<unsigned int>::max());
    std::vector<float> vertices_welded;
    vertices_welded.reserve(static_cast<std::size_t>(vertex_size_orig) * 3);

    unsigned int new_vertex_count = 0;
    for (int i = 0; i < vertex_size_orig; ++i) {
        unsigned int r = rep[static_cast<std::size_t>(i)];
        if (rep_to_new[r] == std::numeric_limits<unsigned int>::max()) {
            rep_to_new[r] = new_vertex_count++;
            vertices_welded.push_back(_vertices[r * 3 + 0]);
            vertices_welded.push_back(_vertices[r * 3 + 1]);
            vertices_welded.push_back(_vertices[r * 3 + 2]);
        }
    }

    // old index → new index
    std::vector<unsigned int> old_to_new(static_cast<std::size_t>(vertex_size_orig));
    for (int i = 0; i < vertex_size_orig; ++i) {
        unsigned int r = rep[static_cast<std::size_t>(i)];
        old_to_new[static_cast<std::size_t>(i)] = rep_to_new[r];
    }

    std::vector<unsigned int> indices_welded(_facetIndices.size());
    for (std::size_t k = 0; k < _facetIndices.size(); ++k) {
        unsigned int old_vid = _facetIndices[k];
        if (old_vid >= static_cast<unsigned int>(vertex_size_orig)) {
            return;
        }
        indices_welded[k] = old_to_new[old_vid];
    }

    int vertex_size = static_cast<int>(vertices_welded.size() / 3);
    std::vector<float> fixed_vertices;
    std::vector<v_conncted> v_connection_table;
    fixed_vertices.resize(static_cast<std::size_t>(vertex_size) * 3u);
    v_connection_table.resize(static_cast<std::size_t>(vertex_size));

    // get connection table
    for (int i = 0; i < facet_size; i++) {
        for (int j = 0; j < 3; j++) {
            unsigned int vid = indices_welded[static_cast<std::size_t>(i) * 3u + static_cast<std::size_t>(j)];
            if (vid < static_cast<unsigned int>(vertex_size)) {
                v_connection_table[vid].facet_id.push_back(i);
                v_connection_table[vid].id_in_f.push_back(static_cast<unsigned int>(j));
            }
            else {
                // wrong index
                return;
            }
        }
    }


    // Boundary / Feature protection
    float feature_angle_deg = 60.0f;     // 夾角超過此值視為銳邊
    float feature_cos_th = std::cos(feature_angle_deg * 3.14159265f / 180.0f);

    std::vector<uint8_t> is_boundary(static_cast<std::size_t>(vertex_size), 0);
    std::vector<uint8_t> is_feature(static_cast<std::size_t>(vertex_size), 0);

    // face normals (normalized). facet_size triangles, indices_welded size = facet_size*3
    struct N3 { float x, y, z; };
    std::vector<N3> face_n(static_cast<std::size_t>(facet_size), { 0,0,0 });
    std::vector<uint8_t> face_valid(static_cast<std::size_t>(facet_size), 0);

    auto cross3 = [](float ax, float ay, float az, float bx, float by, float bz) -> N3 {
        return { ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx };
        };
    auto dot3 = [](const N3& a, const N3& b) -> float { return a.x * b.x + a.y * b.y + a.z * b.z; };

    for (int fi = 0; fi < facet_size; ++fi) {
        unsigned int ia = indices_welded[static_cast<std::size_t>(fi) * 3u + 0u];
        unsigned int ib = indices_welded[static_cast<std::size_t>(fi) * 3u + 1u];
        unsigned int ic = indices_welded[static_cast<std::size_t>(fi) * 3u + 2u];
        if (ia >= (unsigned)vertex_size || ib >= (unsigned)vertex_size || ic >= (unsigned)vertex_size) continue;
        if (ia == ib || ib == ic || ic == ia) continue; // degenerate after weld

        float ax = vertices_welded[static_cast<std::size_t>(ia) * 3u + 0u];
        float ay = vertices_welded[static_cast<std::size_t>(ia) * 3u + 1u];
        float az = vertices_welded[static_cast<std::size_t>(ia) * 3u + 2u];
        float bx = vertices_welded[static_cast<std::size_t>(ib) * 3u + 0u];
        float by = vertices_welded[static_cast<std::size_t>(ib) * 3u + 1u];
        float bz = vertices_welded[static_cast<std::size_t>(ib) * 3u + 2u];
        float cx = vertices_welded[static_cast<std::size_t>(ic) * 3u + 0u];
        float cy = vertices_welded[static_cast<std::size_t>(ic) * 3u + 1u];
        float cz = vertices_welded[static_cast<std::size_t>(ic) * 3u + 2u];

        float e1x = bx - ax, e1y = by - ay, e1z = bz - az;
        float e2x = cx - ax, e2y = cy - ay, e2z = cz - az;

        N3 n = cross3(e1x, e1y, e1z, e2x, e2y, e2z);
        float len2 = n.x * n.x + n.y * n.y + n.z * n.z;
        if (len2 <= 1e-20f) continue;

        float invLen = 1.0f / std::sqrt(len2);
        n.x *= invLen; n.y *= invLen; n.z *= invLen;
        face_n[static_cast<std::size_t>(fi)] = n;
        face_valid[static_cast<std::size_t>(fi)] = 1;
    }

    // Edge map: undirected edge -> up to two adjacent faces
    struct EdgeInfo { int f0 = -1; int f1 = -1; };
    struct U64Hash { std::size_t operator()(std::uint64_t k) const noexcept { return (std::size_t)k; } };

    std::unordered_map<std::uint64_t, EdgeInfo, U64Hash> edge_faces;
    edge_faces.reserve(static_cast<std::size_t>(facet_size) * 3u);

    auto edge_key = [](unsigned int a, unsigned int b) -> std::uint64_t {
        unsigned int lo = (a < b) ? a : b;
        unsigned int hi = (a < b) ? b : a;
        return (std::uint64_t(lo) << 32) | std::uint64_t(hi);
        };

    for (int fi = 0; fi < facet_size; ++fi) {
        unsigned int v[3] = {
            indices_welded[static_cast<std::size_t>(fi) * 3u + 0u],
            indices_welded[static_cast<std::size_t>(fi) * 3u + 1u],
            indices_welded[static_cast<std::size_t>(fi) * 3u + 2u]
        };
        // degenerate triangles can be ignored for edge classification too
        if (v[0] == v[1] || v[1] == v[2] || v[2] == v[0]) continue;

        for (int e = 0; e < 3; ++e) {
            unsigned int a = v[e];
            unsigned int b = v[(e + 1) % 3];
            if (a == b) continue;

            std::uint64_t k = edge_key(a, b);
            auto& info = edge_faces[k];
            if (info.f0 == -1) info.f0 = fi;
            else if (info.f1 == -1) info.f1 = fi;
            else {
                // non-manifold: more than 2 faces share an edge, ignore extras
            }
        }
    }

    // Mark boundary + feature vertices
    for (const auto& kv : edge_faces) {
        std::uint64_t k = kv.first;
        const EdgeInfo& info = kv.second;

        unsigned int a = static_cast<unsigned int>(k >> 32);
        unsigned int b = static_cast<unsigned int>(k & 0xffffffffu);

        if (a >= (unsigned)vertex_size || b >= (unsigned)vertex_size) continue;

        if (info.f1 == -1) {
            // boundary edge
            is_boundary[a] = 1;
            is_boundary[b] = 1;
            continue;
        }

        // feature edge by dihedral angle
        int f0 = info.f0, f1 = info.f1;
        if (f0 < 0 || f1 < 0) continue;
        if (!face_valid[(std::size_t)f0] || !face_valid[(std::size_t)f1]) continue;

        float c = dot3(face_n[(std::size_t)f0], face_n[(std::size_t)f1]); // cos(theta)
        // theta > feature_angle => cos(theta) < cos(threshold)
        if (c < feature_cos_th) {
            is_feature[a] = 1;
            is_feature[b] = 1;
        }
    }


    std::vector<int> neighbors;
    for (int iter = 0; iter < iterations; ++iter) {

        for (int i = 0; i < vertex_size; i++) {
            // protect boundary / feature vertices
            if (is_boundary[static_cast<std::size_t>(i)] || is_feature[static_cast<std::size_t>(i)]) {
                fixed_vertices[static_cast<std::size_t>(i) * 3u + 0u] = vertices_welded[static_cast<std::size_t>(i) * 3u + 0u];
                fixed_vertices[static_cast<std::size_t>(i) * 3u + 1u] = vertices_welded[static_cast<std::size_t>(i) * 3u + 1u];
                fixed_vertices[static_cast<std::size_t>(i) * 3u + 2u] = vertices_welded[static_cast<std::size_t>(i) * 3u + 2u];
                continue;
            }

            neighbors.clear();
            for (std::size_t j = 0; j < v_connection_table[i].id_in_f.size(); j++) {
                unsigned int id_in_f = v_connection_table[i].id_in_f[j];
                int next_v_id = static_cast<int>((id_in_f + 1u) % 3u);
                neighbors.push_back(static_cast<int>(
                    indices_welded[v_connection_table[i].facet_id[j] * 3u + next_v_id]));

                next_v_id = static_cast<int>((id_in_f + 2u) % 3u);
                neighbors.push_back(static_cast<int>(
                    indices_welded[v_connection_table[i].facet_id[j] * 3u + next_v_id]));
            }

            int _s = 0;
            float _x = 0.0f, _y = 0.0f, _z = 0.0f;

            for (std::size_t j = 0; j < neighbors.size(); j++) {
                int nid = neighbors[j];
                if (nid == i) continue;

                bool dup = false;
                for (std::size_t k = 0; k < j; k++) {
                    if (neighbors[k] == nid) { dup = true; break; }
                }
                if (dup) continue;

                _s++;
                _x += vertices_welded[nid * 3u + 0u];
                _y += vertices_welded[nid * 3u + 1u];
                _z += vertices_welded[nid * 3u + 2u];
            }

            if (_s == 0) {
                fixed_vertices[i * 3 + 0] = vertices_welded[i * 3 + 0];
                fixed_vertices[i * 3 + 1] = vertices_welded[i * 3 + 1];
                fixed_vertices[i * 3 + 2] = vertices_welded[i * 3 + 2];
                continue;
            }

            float _g = 0.375f + 0.2f * std::cos(2.0f * 3.14159265f / float(_s));
            float _a = (0.625f - _g * _g) / float(_s);
            float _b = 1.0f - float(_s) * _a;

            float vx = vertices_welded[i * 3 + 0];
            float vy = vertices_welded[i * 3 + 1];
            float vz = vertices_welded[i * 3 + 2];

            float nx = _b * vx + _a * _x;
            float ny = _b * vy + _a * _y;
            float nz = _b * vz + _a * _z;

            fixed_vertices[i * 3 + 0] = vx + lambda * (nx - vx);
            fixed_vertices[i * 3 + 1] = vy + lambda * (ny - vy);
            fixed_vertices[i * 3 + 2] = vz + lambda * (nz - vz);
        }

        // swap back for next iteration
        for (int i = 0; i < vertex_size; ++i) {
            vertices_welded[i * 3 + 0] = fixed_vertices[i * 3 + 0];
            vertices_welded[i * 3 + 1] = fixed_vertices[i * 3 + 1];
            vertices_welded[i * 3 + 2] = fixed_vertices[i * 3 + 2];
        }
    }
    for (int i = 0; i < vertex_size_orig; ++i) {
        unsigned int nid = old_to_new[static_cast<std::size_t>(i)];
        _vertices[static_cast<std::size_t>(i) * 3u + 0u] = vertices_welded[static_cast<std::size_t>(nid) * 3u + 0u];
        _vertices[static_cast<std::size_t>(i) * 3u + 1u] = vertices_welded[static_cast<std::size_t>(nid) * 3u + 1u];
        _vertices[static_cast<std::size_t>(i) * 3u + 2u] = vertices_welded[static_cast<std::size_t>(nid) * 3u + 2u];
    }
}




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

voroffset3d::CompressedVolume voroffset3d::CreateDexelsFromMeshBuffers(
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
        if(_voxelSize <= 0)
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

/**
 * @brief Convert dexel data into mesh buffers.
 *
 * For each cell in the dexel grid, each pair of z-values (span)
 * is converted into a hexahedron (cube) whose base is determined by the
 * cell's x-y position and whose z extents are determined by the span.
 *
 * The hexahedron is then triangulated into 12 triangles.
 *
 * The output vertices are stored in _vertices as a flat float array (x,y,z interleaved),
 * and the triangle facet indices are stored in _facetIndices.
 *
 * @param _dexels [input] The compressed dexel volume.
 * @param _vertices [output] The output vertex array (x, y, z interleaved, float values).
 * @param _facetIndices [output] The output triangle indices (unsigned int), 3 per triangle.
 */
void voroffset3d::DumpDexelsIntoMeshBuffers(
    const vor3d::CompressedVolume &_dexels,
    std::vector<float> &_vertices,
    std::vector<unsigned int> &_facetIndices
) {
    // Get grid dimensions from the dexel volume
    // Assuming gridSize() returns a vector or array with at least two elements: {width, height}
    auto gridSize = _dexels.gridSize();
    int gridW = gridSize[0];
    int gridH = gridSize[1];

    // Get voxel spacing and origin (convert to float)
    float spacing = static_cast<float>(_dexels.spacing());
    auto origin = _dexels.origin(); // Assuming Eigen::Vector3d
    float originArr[3] = { static_cast<float>(origin[0]),
                           static_cast<float>(origin[1]),
                           static_cast<float>(origin[2]) };

    // For each cell in the dexel grid:
    for (int y = 0; y < gridH; ++y) {
        for (int x = 0; x < gridW; ++x) {
            // Get the list of z-span values for cell (x, y)
            const std::vector<double> &spans = _dexels.at(x, y);
            // Process each span (each span is stored as a pair: [z_min, z_max])
            for (size_t i = 0; 2 * i < spans.size(); ++i) {
                // Compute the minimum and maximum 3D coordinates for this cell span
                float xmin = originArr[0] + x * spacing;
                float ymin = originArr[1] + y * spacing;
                float zmin = static_cast<float>(spans[2 * i] * _dexels.spacing());
                float xmax = originArr[0] + (x + 1) * spacing;
                float ymax = originArr[1] + (y + 1) * spacing;
                float zmax = static_cast<float>(spans[2 * i + 1] * _dexels.spacing());

                // Define 8 vertices of the hexahedron (cube)
                // A = (xmin, ymin, zmin)
                // B = (xmax, ymin, zmin)
                // C = (xmin, ymax, zmin)
                // D = (xmax, ymax, zmin)
                // E = (xmin, ymin, zmax)
                // F = (xmax, ymin, zmax)
                // G = (xmin, ymax, zmax)
                // H = (xmax, ymax, zmax)
                float cubeVerts[8][3] = {
                    { xmin, ymin, zmin }, // A, index 0
                    { xmax, ymin, zmin }, // B, index 1
                    { xmin, ymax, zmin }, // C, index 2
                    { xmax, ymax, zmin }, // D, index 3
                    { xmin, ymin, zmax }, // E, index 4
                    { xmax, ymin, zmax }, // F, index 5
                    { xmin, ymax, zmax }, // G, index 6
                    { xmax, ymax, zmax }  // H, index 7
                };

                // Base index for these 8 vertices in the _vertices vector.
                // Each vertex is 3 floats.
                unsigned int baseIndex = static_cast<unsigned int>(_vertices.size() / 3);

                // Append the 8 vertices to _vertices vector.
                for (int vi = 0; vi < 8; ++vi) {
                    _vertices.push_back(cubeVerts[vi][0]);
                    _vertices.push_back(cubeVerts[vi][1]);
                    _vertices.push_back(cubeVerts[vi][2]);
                }

                // Triangulate the hexahedron into 12 triangles.
                // We define the cube faces as follows (using our vertex order above):
                //
                // Back face (z = zmin): A, B, D, C
                // Front face (z = zmax): E, F, H, G
                // Left face (x = xmin): A, C, G, E
                // Right face (x = xmax): B, D, H, F
                // Top face (y = ymax): C, D, H, G
                // Bottom face (y = ymin): A, B, F, E
                //
                // Each face is split into 2 triangles.
                // The triangle vertex indices are offset by baseIndex.
                unsigned int cubeFaces[6][4] = {
                    { 0, 1, 3, 2 }, // Back face: A, B, D, C
                    { 4, 5, 7, 6 }, // Front face: E, F, H, G
                    { 0, 2, 6, 4 }, // Left face: A, C, G, E
                    { 1, 3, 7, 5 }, // Right face: B, D, H, F
                    { 2, 3, 7, 6 }, // Top face: C, D, H, G
                    { 0, 1, 5, 4 }  // Bottom face: A, B, F, E
                };
                // For each face, add 2 triangles
                for (int face = 0; face < 6; ++face) {
                    // First triangle: face[0], face[1], face[2]
                    _facetIndices.push_back(baseIndex + cubeFaces[face][0]);
                    _facetIndices.push_back(baseIndex + cubeFaces[face][1]);
                    _facetIndices.push_back(baseIndex + cubeFaces[face][2]);
                    // Second triangle: face[0], face[2], face[3]
                    _facetIndices.push_back(baseIndex + cubeFaces[face][0]);
                    _facetIndices.push_back(baseIndex + cubeFaces[face][2]);
                    _facetIndices.push_back(baseIndex + cubeFaces[face][3]);
                }
            } // end for each span
        } // end for x
    } // end for y
}

//────────────────────────────────────────────────────────────
// Marching Cubes lookup tables（edgeTable 與 triTable）
static int edgeTable[256]={	// 边表（12位二进制数），对应256种情况；通过8位的顶点状态索引
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0  
};

static int triTable[256][16] ={	// 三角形表，对应256种情况；16个元素的列表，因为最多5个等值面，不为-1的元素表示在该编号的边上有交点
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
    {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
    {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
    {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
    {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
    {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
    {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
    {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
    {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
    {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};


//────────────────────────────────────────────────────────────
// 線性插值，算出等值面上邊線上的交點
Vertex VertexInterp(float iso, const Vertex &p1, const Vertex &p2, float valp1, float valp2) {
    float mu = (iso - valp1) / (valp2 - valp1);
    return {
        p1.x + mu * (p2.x - p1.x),
        p1.y + mu * (p2.y - p1.y),
        p1.z + mu * (p2.z - p1.z)
    };
}

//────────────────────────────────────────────────────────────
// Prototype：將 dexel 轉換為均勻 voxel，並利用 MC 抽取 mesh
void voroffset3d::DumpDexelsToVoxelsMC(
    const vor3d::CompressedVolume &_dexels,
    std::vector<float> &_vertices,          // 依序存放 x,y,z (output)
    std::vector<unsigned int> &_facetIndices  // 每三個為一個三角形 (output)
) {

    bool debugTrigger = false;

    // Step 1: 取得 grid 大小與基本資訊
    auto gridSize = _dexels.gridSize();
    int gridW = gridSize[0];
    int gridH = gridSize[1];
    float spacing = static_cast<float>(_dexels.spacing());
    auto originVec = _dexels.origin();
    float origin[3] = { static_cast<float>(originVec[0]),
                        static_cast<float>(originVec[1]),
                        static_cast<float>(originVec[2]) };

    // Step 2: 計算 z 方向的總 voxel 數
    // 依據所有 dexel span 的最大 z 值決定
    int gridZ = 0;
    for (int y = 0; y < gridH; ++y) {
        for (int x = 0; x < gridW; ++x) {
            const std::vector<double>& spans = _dexels.at(x, y);
            for (size_t i = 0; 2 * i < spans.size(); ++i) {
                float span_zmax = static_cast<float>(spans[2 * i + 1] * _dexels.spacing()) + origin[2];
                int z_index = static_cast<int>(std::ceil((span_zmax - origin[2]) / spacing));
                gridZ = std::max(gridZ, z_index);
            }
        }
    }

    // Step 3: 建立 occupancy grid (尺寸 gridW × gridH × gridZ)
    // 每個 voxel cell 為一 cube，中心判斷是否在 dexel span 內
    std::vector<bool> occupancy(gridW * gridH * gridZ, false);
    auto index3D = [gridW, gridH](int x, int y, int z) -> int {
        return x + gridW * (y + gridH * z);
    };

    for (int y = 0; y < gridH; ++y) {
        for (int x = 0; x < gridW; ++x) {
            const std::vector<double>& spans = _dexels.at(x, y);
            for (size_t i = 0; 2 * i < spans.size(); ++i) {
                float span_zmin = static_cast<float>(spans[2 * i] * _dexels.spacing()) + origin[2];
                float span_zmax = static_cast<float>(spans[2 * i + 1] * _dexels.spacing()) + origin[2];
                // 對該 cell 下所有 z 層，利用 voxel cube 中心判斷是否在 span 內
                for (int z = 0; z < gridZ; ++z) {
                    float voxelCenterZ = origin[2] + (z + 0.5f) * spacing;
                    if (voxelCenterZ >= span_zmin && voxelCenterZ <= span_zmax) {
                        occupancy[index3D(x, y, z)] = true;
                    }
                }
            }
        }
    }

    // debug section
    if (debugTrigger) {
        for (int z = 0; z < gridZ; ++z) {
            for (int y = 0; y < gridH; ++y) {
                for (int x = 0; x < gridW; ++x) {
                    if (!occupancy[index3D(x, y, z)])
                        continue;
    
                    float xmin = origin[0] + x * spacing;
                    float ymin = origin[1] + y * spacing;
                    float zmin = origin[2] + z * spacing;
                    float xmax = origin[0] + (x + 1) * spacing;
                    float ymax = origin[1] + (y + 1) * spacing;
                    float zmax = origin[2] + (z + 1) * spacing;
    
                    float cubeVerts[8][3] = {
                        { xmin, ymin, zmin }, // 0
                        { xmax, ymin, zmin }, // 1
                        { xmin, ymax, zmin }, // 2
                        { xmax, ymax, zmin }, // 3
                        { xmin, ymin, zmax }, // 4
                        { xmax, ymin, zmax }, // 5
                        { xmin, ymax, zmax }, // 6
                        { xmax, ymax, zmax }  // 7
                    };
    
                    unsigned int baseIdx = static_cast<unsigned int>(_vertices.size() / 3);
                    for (int i = 0; i < 8; ++i) {
                        _vertices.push_back(cubeVerts[i][0]);
                        _vertices.push_back(cubeVerts[i][1]);
                        _vertices.push_back(cubeVerts[i][2]);
                    }
    
                    unsigned int faces[6][4] = {
                        { 0, 1, 3, 2 }, // bottom
                        { 4, 5, 7, 6 }, // top
                        { 0, 2, 6, 4 }, // left
                        { 1, 3, 7, 5 }, // right
                        { 2, 3, 7, 6 }, // front
                        { 0, 1, 5, 4 }  // back
                    };
    
                    for (int f = 0; f < 6; ++f) {
                        _facetIndices.push_back(baseIdx + faces[f][0]);
                        _facetIndices.push_back(baseIdx + faces[f][1]);
                        _facetIndices.push_back(baseIdx + faces[f][2]);
    
                        _facetIndices.push_back(baseIdx + faces[f][0]);
                        _facetIndices.push_back(baseIdx + faces[f][2]);
                        _facetIndices.push_back(baseIdx + faces[f][3]);
                    }
                }
            }
        }
    
        return; // 跳過 MC
    }

    int dimX = gridW + 1, dimY = gridH + 1, dimZ = gridZ + 1;
    std::vector<float> scalar(dimX * dimY * dimZ, 1.0f);
    auto idxV = [dimX, dimY](int i, int j, int k) {
        return i + dimX * (j + dimY * k);
    };

    if (debugTrigger) {
        std::cout << ">>> dim(x, y, z) = (" << dimX << ", " << dimY << ", " << dimZ << ")" << std::endl;
    }

    for (int k = 0; k < dimZ; ++k) {
        for (int j = 0; j < dimY; ++j) {
            for (int i = 0; i < dimX; ++i) {
                int count = 0;
                for (int dz = -1; dz <= 0; ++dz)
                    for (int dy = -1; dy <= 0; ++dy)
                        for (int dx = -1; dx <= 0; ++dx) {
                            int x = i + dx, y = j + dy, z = k + dz;
                            if (x >= 0 && x < gridW && y >= 0 && y < gridH && z >= 0 && z < gridZ) {
                                if (occupancy[index3D(x, y, z)]) count++;
                            }
                        }
                scalar[idxV(i,j,k)] = (count == 8) ? 0.0f : 1.0f;
                
                if (debugTrigger) {
                    std::cout 
                    << ">>> iso(" << i << ", " << j << ", " << k << ") = " 
                    << scalar[idxV(i,j,k)] << std::endl;
                }
            }
        }
    }

    const int vertexOffset[8][3] = {
        {0,0,0}, {0,1,0}, {1,1,0}, {1,0,0},
        {0,0,1}, {0,1,1}, {1,1,1}, {1,0,1}
    };

    const float iso = 0.5f;
    for (int k = 0; k < dimZ - 1; ++k) {
        for (int j = 0; j < dimY - 1; ++j) {
            for (int i = 0; i < dimX - 1; ++i) {
                Vertex p[8]; float val[8];
                for (int n = 0; n < 8; ++n) {
                    int dx = vertexOffset[n][0];
                    int dy = vertexOffset[n][1];
                    int dz = vertexOffset[n][2];
                    int ix = i + dx, iy = j + dy, iz = k + dz;
                    p[n] = { origin[0] + ix * spacing,
                             origin[1] + iy * spacing,
                             origin[2] + iz * spacing };
                    val[n] = scalar[idxV(ix, iy, iz)];
                }
                int idx = 0;
                for (int n = 0; n < 8; ++n)
                    if (val[n] < iso) idx |= (1 << n);
                if (edgeTable[idx] == 0) continue;
                Vertex v[12];
                if (edgeTable[idx] & 1) v[0] = VertexInterp(iso, p[0], p[1], val[0], val[1]);
                if (edgeTable[idx] & 2) v[1] = VertexInterp(iso, p[1], p[2], val[1], val[2]);
                if (edgeTable[idx] & 4) v[2] = VertexInterp(iso, p[2], p[3], val[2], val[3]);
                if (edgeTable[idx] & 8) v[3] = VertexInterp(iso, p[3], p[0], val[3], val[0]);
                if (edgeTable[idx] & 16) v[4] = VertexInterp(iso, p[4], p[5], val[4], val[5]);
                if (edgeTable[idx] & 32) v[5] = VertexInterp(iso, p[5], p[6], val[5], val[6]);
                if (edgeTable[idx] & 64) v[6] = VertexInterp(iso, p[6], p[7], val[6], val[7]);
                if (edgeTable[idx] & 128) v[7] = VertexInterp(iso, p[7], p[4], val[7], val[4]);
                if (edgeTable[idx] & 256) v[8] = VertexInterp(iso, p[0], p[4], val[0], val[4]);
                if (edgeTable[idx] & 512) v[9] = VertexInterp(iso, p[1], p[5], val[1], val[5]);
                if (edgeTable[idx] & 1024) v[10] = VertexInterp(iso, p[2], p[6], val[2], val[6]);
                if (edgeTable[idx] & 2048) v[11] = VertexInterp(iso, p[3], p[7], val[3], val[7]);

                for (int t = 0; triTable[idx][t] != -1; t += 3) {
                    unsigned base = static_cast<unsigned>(_vertices.size() / 3);
                    for (int vi = 0; vi < 3; ++vi) {
                        Vertex &pt = v[triTable[idx][t + vi]];
                        _vertices.insert(_vertices.end(), { pt.x, pt.y, pt.z });
                        _facetIndices.push_back(base + vi);
                    }
                }
            }
        }
    }

    laplacian_smooth(_vertices, _facetIndices);

}