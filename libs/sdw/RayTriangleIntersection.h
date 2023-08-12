#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include "ModelTriangle.h"

struct RayTriangleIntersection {
	glm::vec3 intersectionPoint;
	double distanceFromCamera;
	ModelTriangle intersectedTriangle;
	size_t triangleIndex;
	bool rayHasHitTriangle;

	RayTriangleIntersection();

	RayTriangleIntersection(const glm::vec3& point, double distance, const ModelTriangle& triangle, size_t index);
	friend std::ostream& operator<<(std::ostream& os, const RayTriangleIntersection& intersection);

	RayTriangleIntersection(const glm::vec3 &point, double distance, const ModelTriangle &triangle, size_t index, bool isNull);
	friend std::ostream &operator<<(std::ostream &os, const RayTriangleIntersection &intersection);
};
