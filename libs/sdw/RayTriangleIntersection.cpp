#include "RayTriangleIntersection.h"

RayTriangleIntersection::RayTriangleIntersection() = default;
RayTriangleIntersection::RayTriangleIntersection(const glm::vec3 &point, double distance, const ModelTriangle &triangle, size_t index) :
	intersectionPoint(point),
	distanceFromCamera(distance),
	intersectedTriangle(triangle),
	triangleIndex(index) {}

RayTriangleIntersection::RayTriangleIntersection(const glm::vec3 &point, double distance, const ModelTriangle &triangle, size_t index, bool isNull) :
		intersectionPoint(point),
		distanceFromCamera(distance),
		intersectedTriangle(triangle),
		rayHasHitTriangle(isNull),
		triangleIndex(index) {}

std::ostream &operator<<(std::ostream &os, const RayTriangleIntersection &intersection) {
	os << "Intersection is at [" << intersection.intersectionPoint[0] << "," << intersection.intersectionPoint[1] << "," <<
	   intersection.intersectionPoint[2] << "] on triangle " << intersection.intersectedTriangle <<
	   " at a distance of " << intersection.distanceFromCamera;
	return os;
}
