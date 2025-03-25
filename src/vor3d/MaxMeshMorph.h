#pragma once

#ifdef MMM_WRAPPER_EXPORTS
#define MMM_WRAPPER_EXPORTS __declspec(dllexport)
#else
#define MMM_WRAPPER_EXPORTS __declspec(dllimport)
#endif

#include <memory>
#include <string>
// #include <thread>

enum class MaxMeshMorphOperator {
	mmm_none = 0,
	mmm_dilation = 1,
	mmm_erosion = 2,
	mmm_closing = 3,
	mmm_opening = 4
};

class /*MM_WRAPPER_EXPORTS*/ MaxMeshMorphArgs {
public:
	double radius = 10;
	double dexels_size = 1;
	int padding = 0;
	int num_dexels = 256;
	// unsigned int num_thread = std::max(1u, std::thread::hardware_concurrency());
	bool force = false;
	bool radius_in_mm = true;
};

class MMM_WRAPPER_EXPORTS MaxMeshMorphWrapper {
public:
	MaxMeshMorphWrapper();
	~MaxMeshMorphWrapper();
	void Execute(
		double* in_vertices, int in_facet_num,
		double*& out_vertices, int& out_facet_num,
		MaxMeshMorphOperator mmm_operator,
		MaxMeshMorphArgs mmm_args
	);
	void Clear();
};