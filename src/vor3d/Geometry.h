#pragma once

#include <array>

class Vertex {
	public:
		Vertex() : x(0.0f), y(0.0f), z(0.0f) {}
		Vertex(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}		
		Vertex operator+(const Vertex& other) const {
			return Vertex(x + other.x, y + other.y, z + other.z);
		}
		Vertex operator-(const Vertex& other) const {
			return Vertex(x - other.x, y - other.y, z - other.z);
		}
		Vertex& operator+=(const Vertex& other) {
			x += other.x;
			y += other.y;
			z += other.z;
			return *this;
		}
		Vertex& operator-=(const Vertex& other) {
			x -= other.x;
			y -= other.y;
			z -= other.z;
			return *this;
		}
	public:
		float x, y, z;
};

class Facet {
	public:
		Facet() : indices{0, 0, 0} {}
		Facet(unsigned int i0, unsigned int i1, unsigned int i2) : indices{{i0, i1, i2}} {}	
		unsigned int& operator[](int i) { return indices[i]; }
		const unsigned int& operator[](int i) const { return indices[i]; }
	public:
		std::array<unsigned int, 3> indices;
};