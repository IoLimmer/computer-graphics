#include "ModelTriangle.h"
#include <utility>

ModelTriangle::ModelTriangle() = default;

ModelTriangle::ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Colour trigColour) :
		vertices({{v0, v1, v2}}), texturePoints(), vertexNormals(), colour(std::move(trigColour)), normal() {}

std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle) {
	os << "Vertex 0 = (" << triangle.vertices[0].x << ", " << triangle.vertices[0].y << ", " << triangle.vertices[0].z << ")\n";
	os << "Vertex 1 = (" << triangle.vertices[1].x << ", " << triangle.vertices[1].y << ", " << triangle.vertices[1].z << ")\n";
	os << "Vertex 2 = (" << triangle.vertices[2].x << ", " << triangle.vertices[2].y << ", " << triangle.vertices[2].z << ")\n";
	os << "Colour = " << triangle.colour << "\n";
	os << "Face Normal = (" << triangle.normal.x << ", " << triangle.normal.y << ", " << triangle.normal.z << ")\n";
	os << "Vertex Normal 0 = (" << triangle.vertexNormals[0].x << ", " << triangle.vertexNormals[0].y << ", " << triangle.vertexNormals[0].z << ")\n";
	os << "Vertex Normal 1 = (" << triangle.vertexNormals[1].x << ", " << triangle.vertexNormals[1].y << ", " << triangle.vertexNormals[1].z << ")\n";
	os << "Vertex Normal 2 = (" << triangle.vertexNormals[2].x << ", " << triangle.vertexNormals[2].y << ", " << triangle.vertexNormals[2].z << ")\n";
	return os;

}

//std::ostream& operator<<(std::ostream& os, const ModelTriangle& triangle) {
//	os << "Vertex 0 = (" << triangle.vertices[0].x << ", " << triangle.vertices[0].y << ", " << triangle.vertices[0].z << ")\n";
//	os << "Vertex 1 = (" << triangle.vertices[1].x << ", " << triangle.vertices[1].y << ", " << triangle.vertices[1].z << ")\n";
//	os << "Vertex 2 = (" << triangle.vertices[2].x << ", " << triangle.vertices[2].y << ", " << triangle.vertices[2].z << ")\n";
//	os << "Colour = " << triangle.colour << "\n";
//	os << "Face Normal = (" << triangle.normal.x << ", " << triangle.normal.y << ", " << triangle.normal.z << ")\n";
//	os << "Texture Point 0 = (" << triangle.texturePoints[0].x << ", " << triangle.texturePoints[0].y << ")\n";
//	os << "Texture Point 1 = (" << triangle.texturePoints[1].x << ", " << triangle.texturePoints[1].y << ")\n";
//	os << "Texture Point 2 = (" << triangle.texturePoints[2].x << ", " << triangle.texturePoints[2].y << ")\n";
//	os << "Vertex Normal 0 = (" << triangle.vertexNormals[0].x << ", " << triangle.vertexNormals[0].y << ", " << triangle.vertexNormals[0].z << ")\n";
//	os << "Vertex Normal 1 = (" << triangle.vertexNormals[1].x << ", " << triangle.vertexNormals[1].y << ", " << triangle.vertexNormals[1].z << ")\n";
//	os << "Vertex Normal 2 = (" << triangle.vertexNormals[2].x << ", " << triangle.vertexNormals[2].y << ", " << triangle.vertexNormals[2].z << ")\n";
//	return os;
//
//}

//std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle) {
//	os << "(" << triangle.vertices[0].x << ", " << triangle.vertices[0].y << ", " << triangle.vertices[0].z << ")\n";
//	os << "(" << triangle.vertices[1].x << ", " << triangle.vertices[1].y << ", " << triangle.vertices[1].z << ")\n";
//	os << "(" << triangle.vertices[2].x << ", " << triangle.vertices[2].y << ", " << triangle.vertices[2].z << ")\n";
//	return os;
//}
