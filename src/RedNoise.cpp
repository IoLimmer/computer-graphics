#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <TextureMap.h>
#include <CanvasTriangle.h>
#include <ModelTriangle.h>
#include <Colour.h>
#include <RayTriangleIntersection.h>
#include <Utils.h>

#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>
#include <map>
#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <thread>

//using namespace std;
//using namespace glm;

#define WIDTH 640
#define HEIGHT 480
#define noOfCores 8


TextureMap textureMap = TextureMap("../../../src/textures/texture.ppm");

bool lookAtActivated = true;
bool speedyMode = false;
bool incidenceLight = true;
bool softShadows = false;
bool gouraudShading = false;
bool phongShading = false;
float softDist = 0.2f;
float ambient = 0.2f;
glm::vec3 offsets[8] = {
		glm::vec3(0, 0, 0),
		glm::vec3(0, 0, softDist),
		glm::vec3(0, softDist, 0),
		glm::vec3(0, softDist, softDist),
		glm::vec3(0, softDist * -1, softDist * -1),
		glm::vec3(0, 0, softDist * -1),
		glm::vec3(0, softDist * -1, 0),
		glm::vec3(softDist, softDist, softDist),
};
//glm::vec3 offsets[8] = {
//		glm::vec3(0, 0, 0),
//		glm::vec3(0, 0, softDist),
//		glm::vec3(0, softDist, 0),
//		glm::vec3(0, softDist, softDist),
//		glm::vec3(0, softDist * -1, softDist * -1),
//		glm::vec3(0, 0, softDist * -1),
//		glm::vec3(0, softDist * -1, 0),
//		glm::vec3(softDist, softDist, softDist),
//};

float clamp(float val, float min, float max) {
	if (val > max)
		val = max;
	if (val < min)
		val = min;
	return val;
}

// OUTPUT GLM VEC3
void OutputGLMVEC3(glm::vec3 daVec) {
	std::cout << "(" << daVec[0] << ", " << daVec[1] << ", " << daVec[2] << ")\n";
}

// OUTPUT ALLMODELTRIANGLES VECTOR
void OutputALLMODELTRIANGLES(std::vector<std::vector<ModelTriangle>> OBJfile) {
	for (size_t i = 0; i < OBJfile.size(); i++) {
		for (size_t j = 0; j < OBJfile[i].size(); j++) {
			std::cout << OBJfile[i][j] << OBJfile[i][j].colour << "\n\n";
		}
	}
}

// OUTPUT SINGLE VERTEX OF MODEL TRIANGLE
void OutputSingleVertexModelTriangle(ModelTriangle modtri, int i) {
	std::cout << "(" << modtri.vertices[i].x << ", " << modtri.vertices[i].y << ", " << modtri.vertices[i].z << ")\n";
}



// 1D INTERPOLATION
std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	std::vector<float> v;
	float counter = from;
	for (float i = 0; i < numberOfValues; i++) {
		v.push_back(counter);
		counter = counter + ((to - from) / float(numberOfValues));
	}
	return v;
}

//POINT TO POINT INTERPOLATION
std::vector<CanvasPoint> InterpolatePointToPoint(CanvasPoint from, CanvasPoint to) {
	std::vector<CanvasPoint> v;
	CanvasPoint counter = from;

	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	numberOfSteps++;
	float xStepSize = xDiff / numberOfSteps;
	float yStepSize = yDiff / numberOfSteps;

	for (float i = 0.0; i <= numberOfSteps; i++) {
		v.push_back(counter);
		counter.x = int(from.x + (xStepSize * i));
		counter.y = int(from.y + (yStepSize * i));
	}
	return v;
}

//POINT TO POINT INTERPOLATION WITH DEPTH
std::vector<CanvasPoint> InterpolatePointToPointWithDepth(CanvasPoint from, CanvasPoint to) {
	std::vector<CanvasPoint> v;
	CanvasPoint counter = from;

	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float depthDiff = to.depth - from.depth;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	numberOfSteps++;
	float xStepSize = xDiff / numberOfSteps;
	float yStepSize = yDiff / numberOfSteps;
	float depthStepSize = depthDiff / numberOfSteps;

	for (float i = 0.0; i <= numberOfSteps; i++) {
		v.push_back(counter);
		counter.x = int(from.x + (xStepSize * i));
		counter.y = int(from.y + (yStepSize * i));
		counter.depth = int(from.y + (depthStepSize * i));
	}
	return v;
}

//VEC3 TO VEC3 INTERPOLATION
std::vector<glm::vec3> InterpolateVec3ToVec3(glm::vec3 from, glm::vec3 to, int numberOfValues) {
	std::vector<glm::vec3> v;
	glm::vec3 counter = from;
	for (float i = 0; i <= numberOfValues; i++) {
		v.push_back(counter);
		for (float i = 0; i < 3; i++) { counter[i] = counter[i] + ((to[i] - from[i]) / float(numberOfValues)); }
	}
	return v;
}

//DRAW TOP HALF TRIANGLE (TRIANGLE WITH FLAT BASE)
std::vector<CanvasPoint> drawTopHalfTriangle(CanvasPoint topPoint, CanvasPoint point1, CanvasPoint point2, std::vector<CanvasPoint> pointsToFill) {
	std::vector<CanvasPoint> temp;

	float xStep1, xStep2, xCurrent1, xCurrent2;

	xStep1 = (point1.x - topPoint.x) / (point1.y - topPoint.y);
	xStep2 = (point2.x - topPoint.x) / (point2.y - topPoint.y);
	xCurrent1 = topPoint.x;
	xCurrent2 = topPoint.x;

	float depthStep1, depthStep2, depthCurrent1, depthCurrent2;

	depthStep1 = (point1.depth - topPoint.depth) / (point1.y - topPoint.y);
	depthStep2 = (point2.depth - topPoint.depth) / (point2.y - topPoint.y);
	depthCurrent1 = topPoint.depth;
	depthCurrent2 = topPoint.depth;


	for (int currentY = topPoint.y; currentY <= point1.y; currentY++)
	{
		temp = InterpolatePointToPoint(CanvasPoint(int(xCurrent1), currentY, depthCurrent1),
			CanvasPoint(int(xCurrent2), currentY, depthCurrent2));
		temp.push_back(CanvasPoint(int(xCurrent2), currentY, depthCurrent2));
		xCurrent1 += xStep1;
		xCurrent2 += xStep2;
		depthCurrent1 += depthStep1;
		depthCurrent2 += depthStep2;
		for (int j = 0; j < temp.size(); j++) {
			pointsToFill.push_back(temp[j]);
		}
	}

	return pointsToFill;
}

//DRAW BOTTOM HALF TRIANGLE (TRIANGLE WITH FLAT TOP)
std::vector<CanvasPoint> drawBottomHalfTriangle(CanvasPoint bottomPoint, CanvasPoint point1, CanvasPoint point2, std::vector<CanvasPoint> pointsToFill) {
	std::vector<CanvasPoint> temp;

	float xStep1, xStep2, xCurrent1, xCurrent2;
	xStep1 = (bottomPoint.x - point1.x) / (bottomPoint.y - point1.y);
	xStep2 = (bottomPoint.x - point2.x) / (bottomPoint.y - point2.y);
	xCurrent1 = bottomPoint.x;
	xCurrent2 = bottomPoint.x;

	float depthStep1, depthStep2, depthCurrent1, depthCurrent2;

	depthStep1 = (bottomPoint.depth - point1.depth) / (bottomPoint.y - point1.y);
	depthStep2 = (bottomPoint.depth - point2.depth) / (bottomPoint.y - point2.y);
	depthCurrent1 = bottomPoint.depth;
	depthCurrent2 = bottomPoint.depth;

	for (int currentY = bottomPoint.y; currentY > point1.y; currentY--)
	{
		temp = InterpolatePointToPoint(CanvasPoint(int(xCurrent1), currentY, depthCurrent1),
			CanvasPoint(int(xCurrent2), currentY, depthCurrent2));
		temp.push_back(CanvasPoint(int(xCurrent2), currentY));
		xCurrent1 -= xStep1;
		xCurrent2 -= xStep2;
		depthCurrent1 -= depthStep1;
		depthCurrent2 -= depthStep2;
		for (int j = 0; j < temp.size(); j++) {
			pointsToFill.push_back(temp[j]);
		}
	}

	return pointsToFill;
}



// Calculate canvas point from point on 3d object
CanvasPoint GetCanvasIntersectionPoint(glm::vec3 vertexPosition, glm::vec3 cameraPosition, glm::mat3 cameraRotation, float focalLength, float scale) {
	CanvasPoint screenPosition;

	vertexPosition = vertexPosition - cameraPosition;

	vertexPosition = vertexPosition * cameraRotation;

	screenPosition.x = ((focalLength * ((vertexPosition.x * -scale) / vertexPosition.z)) + (WIDTH / 2));
	screenPosition.y = ((focalLength * ((vertexPosition.y * scale) / vertexPosition.z)) + (HEIGHT / 2));

	return screenPosition;
}


float GetPointZDepth(glm::vec3 vertexPosition, glm::vec3 cameraPosition) {
	float depth, z;
	glm::vec3 distanceFromCam;

	distanceFromCam = (vertexPosition - cameraPosition);
	z = sqrt(pow(distanceFromCam.x, 2) + pow(distanceFromCam.y, 2) + pow(distanceFromCam.z, 2));

	//z = abs(vertexPosition.z - cameraPosition.z);
	depth = 1 / (z);
	return depth;
}


// Calculate canvas points for a triangle from face on 3d object
CanvasTriangle ModelTriangleToCanvasTriangle(ModelTriangle modelTriangle, glm::vec3 cameraPosition, glm::mat3 cameraRotation, float focalLength, float scale) {
	CanvasTriangle screenTriangle;
	for (size_t i = 0; i < 3; i++) {
		screenTriangle[i] = GetCanvasIntersectionPoint(modelTriangle.vertices[i], cameraPosition, cameraRotation, focalLength, scale);
		screenTriangle[i].depth = GetPointZDepth(modelTriangle.vertices[i], cameraPosition);
	}
	return screenTriangle;
}


std::vector<ModelTriangle> CalculateVertexNormals(std::vector<ModelTriangle> tris) {

	for (int j = 0; j < tris.size(); j++) {
		for (int n = 0; n < 3; n++) {
			ModelTriangle currTri = tris[j];
			//std::vector<glm::vec3> listOfNormals;

			glm::vec3 vertexNormal = currTri.normal;

			int counter = 1;
			for (int i = 0; i < tris.size(); i++) {

				if ((currTri.vertices[0] == tris[i].vertices[0]) &&
					(currTri.vertices[1] == tris[i].vertices[1]) &&
					(currTri.vertices[2] == tris[i].vertices[2])) {
					// ignore instance of current triangle
					continue;
				}

				if ((currTri.vertices[n] == tris[i].vertices[0]) ||
					(currTri.vertices[n] == tris[i].vertices[1]) ||
					(currTri.vertices[n] == tris[i].vertices[2])) {
					// if the current triangle shares a vertex with a one from the OBJ file

					//listOfNormals.push_back(tris[i].normal);

					vertexNormal += tris[i].normal;
					counter++;
				}
			}
			vertexNormal = glm::normalize(vertexNormal / glm::vec3(counter));
			tris[j].vertexNormals[n] = vertexNormal;


			//std::cout << "OG triangle normal = ";
			//OutputGLMVEC3(currTri.normal);

			//for (int i = 0; i < listOfNormals.size(); i++) {
			//	std::cout << "neighbour " << i << " normal = ";
			//	OutputGLMVEC3(listOfNormals[i]);
			//}

			//std::cout << "\nVertex normal = ";
			//OutputGLMVEC3(vertexNormal);

			//std::cout << "end\n";
		}
		//std::cout << tris[j] << "\n";

	}


	return tris;
}



// LOAD OBJ
std::vector<std::vector<ModelTriangle>> loadOBJFile(std::string filename) {
	std::string line, totalFile;
	std::ifstream fileOBJ, fileMTL;
	std::vector<Colour> mtlColours;


	fileMTL.open("../../../src/models/" + filename + ".mtl");
	if (fileMTL.is_open())
	{
		std::string colourName;
		while (getline(fileMTL, line))
		{
			if (line.substr(0, 7) == "newmtl ") {
				colourName = line.substr(7);
			}
			else if (line.substr(0, 3) == "Kd ") {
				std::istringstream v(line.substr(3));
				double red, green, blue;
				int redToPush, greenToPush, blueToPush;
				v >> red; v >> green; v >> blue;

				redToPush = red * 255.0;
				greenToPush = green * 255.0;
				blueToPush = blue * 255.0;

				Colour colourToPush(colourName, redToPush, greenToPush, blueToPush);
				mtlColours.push_back(colourToPush);
			}
		}
		fileMTL.close();
	}
	else if (fileMTL.fail()) {
		std::cout << "failed to open mtl file\n";
	}


	fileOBJ.open("../../../src/models/" + filename + ".obj");
	std::vector<std::vector<ModelTriangle>> allModelTriangles;
	std::vector<ModelTriangle> currentObject;
	if (fileOBJ.is_open())
	{
		//int objectCounter = 0;

		std::string objectName;
		int red, green, blue;
		std::vector<glm::vec3> vertices;
		std::vector<TexturePoint> texturePoints;
		std::vector<int> faceIndices;
		std::vector<int> texturePointsIndices;
		bool objectIsTextured = false;
		//int vertCounter = 0;
		//glm::vec3 normal;

		while (getline(fileOBJ, line))
		{
			//std::cout << line << "\n";

			if (line.substr(0, 2) == "o ") { //get object name
				objectName = line.substr(2);

				//std::cout << objectName << "\n";

				allModelTriangles.push_back(currentObject); //the biggie
				objectIsTextured = false;
				currentObject.clear();
			}

			else if (line.substr(0, 7) == "usemtl ") { //get object colour
				std::string colourFromFile(line.substr(7));

				for (size_t i = 0; i < mtlColours.size(); i++) {
					if (colourFromFile == mtlColours[i].name) {
						red = mtlColours[i].red;
						green = mtlColours[i].green;
						blue = mtlColours[i].blue;
						break;
					}
				}
				//std::cout << red << ", " << green << ", " << blue << "\n";
			}

			else if (line.substr(0, 2) == "v ") { //get vertices
				std::istringstream v(line.substr(2));
				float x, y, z;
				v >> x; v >> y; v >> z;
				glm::vec3 vert = glm::vec3(x, y, z);
				vertices.push_back(vert);

				//OutputGLMVEC3(vert);
			}

			else if (line.substr(0, 3) == "vt ") { //get vertices
				objectIsTextured = true;
				std::istringstream vt(line.substr(3));
				float x, y;
				vt >> x; vt >> y;
				TexturePoint textPoint = TexturePoint(x, y);
				texturePoints.push_back(textPoint);

				//OutputGLMVEC3(vert);
			}

			else if (line.substr(0, 2) == "f ") { //get faces
				std::istringstream f(line.substr(2));

				faceIndices.clear();
				texturePointsIndices.clear();

				std::string token;
				int faceIndex;
				int texturePointIndex;
				size_t pos = -1;

				if (objectIsTextured) {
					while (f >> token) {
						while ((pos = token.rfind('/'))
							!= std::string::npos) {
							token.replace(pos, 1, " ");
							//token.erase(pos, 1);
						}
						std::istringstream ss(token);

						ss >> faceIndex;
						ss >> texturePointIndex;
						faceIndices.push_back(faceIndex - 1);
						texturePointsIndices.push_back(texturePointIndex - 1);
						//std::cout << faceIndex << "\n";
						//std::cout << texturePointIndex << "\n";
					}
				}
				else {
					while (f >> token) {
						while ((pos = token.rfind('/'))
							!= std::string::npos) {
							token.erase(pos, 1);
						}
						std::istringstream ss(token);

						ss >> faceIndex;
						faceIndices.push_back(faceIndex - 1);
					}
				}

				//std::cout << "\n";

				//std::cout << texturePoints[texturePointsIndices] << "\n";

				// ADD VERTICES TO BIG CURRENTOBJECT
				ModelTriangle faceFromFile;
				glm::vec3 mp0 = vertices[faceIndices[0]];
				glm::vec3 mp1 = vertices[faceIndices[1]];
				glm::vec3 mp2 = vertices[faceIndices[2]];

				Colour colour1(red, green, blue);

				faceFromFile = ModelTriangle(mp0, mp1, mp2, colour1);

				if (objectIsTextured) {
					faceFromFile.texturePoints[0] = texturePoints[texturePointsIndices[0]];
					faceFromFile.texturePoints[1] = texturePoints[texturePointsIndices[1]];
					faceFromFile.texturePoints[2] = texturePoints[texturePointsIndices[2]];
				}

				glm::vec3 normal = glm::normalize(glm::cross((mp1 - mp0), (mp2 - mp0)));
				faceFromFile.normal = normal;

				currentObject.push_back(faceFromFile); //save triangle in current object
			}
		}
		currentObject = CalculateVertexNormals(currentObject);

		allModelTriangles.push_back(currentObject); //the biggie
		fileOBJ.close();
	}
	else if (fileOBJ.fail()) {
		std::cout << "failed to open obj file\n";
	}


	//for (int i = 0; i < allModelTriangles.size(); i++) {
	//	for (int j = 0; j < allModelTriangles[i].size(); j++) {
	//		if (allModelTriangles[i][j].texturePoints[0].x != 0) {
	//			std::cout << allModelTriangles[i][j] << "\n";
	//			std::cout << allModelTriangles[i][j].colour << "\n";
	//			OutputGLMVEC3(allModelTriangles[i][j].normal);
	//			std::cout << allModelTriangles[i][j].texturePoints[0] << ", " << allModelTriangles[i][j].texturePoints[1] << ", " << allModelTriangles[i][j].texturePoints[2] << "\n";
	//			//std::cout << allModelTriangles[i][j].texturePoints << "\n";
	//			std::cout << "\n\n";
	//		}
	//	}
	//}


	return allModelTriangles;
}



glm::vec3 getRaysForScreen(int x, int y, float focalLength) {
	float scale = 40;
	glm::vec3 rayDirection((x - (WIDTH / 2)) / scale, (y - (HEIGHT / 2)) / scale, focalLength);
	//rayDirection.z = rayDirection.z * -1.0f;
	rayDirection = rayDirection * glm::vec3(1.0, -1.0, -1.0);

	return rayDirection;
}

glm::vec3 RayTracerMaths(glm::vec3 startPosition, glm::vec3 rayDirection, ModelTriangle triangle) {
	glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
	glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
	glm::vec3 SPVector = startPosition - triangle.vertices[0];
	glm::mat3 DEMatrix(-rayDirection, e0, e1);
	glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

	return possibleSolution;
}



RayTriangleIntersection getClosestIntersection(glm::vec3 startPoint, glm::vec3 rayDirection, std::vector<std::vector<ModelTriangle>>& OBJfile) {

	glm::vec3 possibleSolution, closestSolution, intersectionPoint1;
	RayTriangleIntersection intersection;
	ModelTriangle triangle;
	float t = 0.0;
	float u = 0.0;
	float v = 0.0;
	float range = 10000000000.0;
	float bestT = 10000000000.0;
	size_t counter = 0;
	int selfHitCounter = 0;
	bool rayHit = false;

	for (size_t i = 0; i < OBJfile.size(); i++) {
		for (size_t j = 0; j < OBJfile[i].size(); j++) {

			triangle = OBJfile[i][j];
			possibleSolution = RayTracerMaths(startPoint, rayDirection, triangle);

			t = possibleSolution.x; // the absolute distance along the ray from the camera to the intersection point
			u = possibleSolution.y; // the proportional distance along the triangle's first edge that the intersection point occurs
			v = possibleSolution.z; // the proportional distance along the triangle's second edge that the intersection point occurs

			intersectionPoint1 = (triangle.vertices[0] +
				(u * (triangle.vertices[1] - triangle.vertices[0])) +
				(v * (triangle.vertices[2] - triangle.vertices[0])));

			if (t > 0 && t < range && (u + v) <= 1 && u >= 0 && v >= 0 && u <= 1 && v <= 1) {

				if (abs(t) < bestT && (glm::length(intersectionPoint1 - startPoint) > 0.0001)) {
					selfHitCounter++;

					bestT = abs(t);
					intersection = RayTriangleIntersection(intersectionPoint1,
						glm::length(startPoint - intersectionPoint1),
						triangle, counter, rayHit);
				}
			}
			if (selfHitCounter > 0) {
				rayHit = true;
			}

			counter++;
			intersection.rayHasHitTriangle = rayHit;
		}
	}
	return intersection;
}


float CalculateLightIntensity(float distanceFromSource) {
	distanceFromSource = abs(distanceFromSource);
	float tempStrength = 100.0;
	double lightIntensity = tempStrength / (4.0f * M_PI * pow(distanceFromSource, 2));

	lightIntensity = clamp(lightIntensity, 0.2, 1.0);

	//lightIntensity = 1.0f;

	return lightIntensity;
}


float CalculateTriArea(std::vector<glm::vec3> points) {
	glm::vec3 A = points[0];
	glm::vec3 B = points[1];
	glm::vec3 C = points[2];

	float AB = sqrt(abs(pow(B.x - A.x, 2) + pow(B.y - A.y, 2)));
	float BC = sqrt(abs(pow(C.x - B.x, 2) + pow(C.y - B.y, 2)));
	float CA = sqrt(abs(pow(A.x - C.x, 2) + pow(A.y - C.y, 2)));

	float semiPerimeter = (AB + BC + CA) / 2.0f;

	float area = sqrt(abs(semiPerimeter * (semiPerimeter - AB) * (semiPerimeter - BC) * (semiPerimeter - CA)));

	return area;
}

//glm::vec3 CalculateBarycentricCoordinates(glm::vec3 pointP, glm::vec3 pointA, glm::vec3 pointB, glm::vec3 pointC) {
glm::vec3 CalculateBarycentricCoordinatesNew(glm::vec3 point, ModelTriangle triangle) {

	glm::vec3 A = glm::vec3(triangle.vertices[0]);
	glm::vec3 B = glm::vec3(triangle.vertices[1]);
	glm::vec3 C = glm::vec3(triangle.vertices[2]);
	glm::vec3 P = glm::vec3(point);


	std::vector<glm::vec3> ABC{ A, B, C };
	std::vector<glm::vec3> CAP{ C, A, P };
	std::vector<glm::vec3> ABP{ A, B, P };
	std::vector<glm::vec3> BCP{ B, C, P };

	float u = CalculateTriArea(CAP) / CalculateTriArea(ABC);
	float v = CalculateTriArea(ABP) / CalculateTriArea(ABC);

	u = clamp(u, 0.0f, 1.0f);
	v = clamp(v, 0.0f, 1.0f);

	float w = 1.0f - u - v;


	return glm::vec3(u, v, w);
}


TexturePoint CanvasPointToTexturePoint(CanvasPoint point) {
	TexturePoint texturePoint;
	texturePoint.x = textureMap.width * (point.x / WIDTH);
	texturePoint.y = textureMap.height * (point.y / HEIGHT);
	return texturePoint;
}


Colour uint32tToColour(uint32_t colour) {
	uint32_t red = (colour & (255 << 16)) >> 16;
	uint32_t green = (colour & (255 << 8)) >> 8;
	uint32_t blue = (colour & 255);
	Colour properColour = Colour(red, green, blue);

	return properColour;
}

//glm::vec3 CalculateBarycentricCoordinates(glm::vec3 pointP, glm::vec3 pointA, glm::vec3 pointB, glm::vec3 pointC) {
glm::vec3 CalculateBarycentricCoordinates(RayTriangleIntersection intersection) {

	glm::vec3 pointP = intersection.intersectionPoint;
	glm::vec3 pointA = intersection.intersectedTriangle.vertices[0];
	glm::vec3 pointB = intersection.intersectedTriangle.vertices[1];
	glm::vec3 pointC = intersection.intersectedTriangle.vertices[2];

	float u, v, w;
	glm::vec3 v0 = pointB - pointA;
	glm::vec3 v1 = pointC - pointA;
	glm::vec3 v2 = pointP - pointA;

	float d00 = glm::dot(v0, v0);
	float d01 = glm::dot(v0, v1);
	float d11 = glm::dot(v1, v1);
	float d20 = glm::dot(v2, v0);
	float d21 = glm::dot(v2, v1);
	float denom = d00 * d11 - d01 * d01;

	u = (d11 * d20 - d01 * d21) / denom;
	v = (d00 * d21 - d01 * d20) / denom;
	w = 1.0f - u - v;

	return glm::vec3(u, v, w);
}

uint32_t TraceLight(glm::vec3 cameraPoint, glm::mat3 cameraRotation, glm::vec3 XYrayDirection, glm::vec3 lightSource, float focalLength, std::vector<std::vector<ModelTriangle>>& OBJfile) {

	Colour outputColour;

	//glm::mat3 tempRotation = cameraRotation * glm::mat3(
	//	1.0, 0.0, 0.0,
	//	0.0, 1.0, 0.0,
	//	0.0, 0.0, 1.0
	//);

	RayTriangleIntersection intersection = getClosestIntersection(cameraPoint, XYrayDirection * cameraRotation, OBJfile);

	if (intersection.rayHasHitTriangle == false) {
		Colour black = Colour(0, 0, 0);
		return (255 << 24) + (black.red << 16) + (black.green << 8) + (black.blue);

	}

	else if (gouraudShading) {

		// get face normals
		ModelTriangle currTri = intersection.intersectedTriangle;
		std::vector<ModelTriangle> neighbourTriangles;

		glm::vec3 currNormal = currTri.normal;

		//Colour red = Colour(255, 0, 0);
		//return ((255 << 24) + (red.red << 16) + (red.green << 8) + (red.blue));

	}

	else if (incidenceLight) {
		glm::vec3 vn0 = intersection.intersectedTriangle.vertexNormals[0];
		glm::vec3 vn1 = intersection.intersectedTriangle.vertexNormals[1];
		glm::vec3 vn2 = intersection.intersectedTriangle.vertexNormals[2];

		// PROXIMITY LIGHTING
		Colour intersectionColour;
		int limit = softShadows ? 8 : 1;
		float intensity = ambient;

		// get barycentric Coordinates
		glm::vec3 baryCoords = CalculateBarycentricCoordinates(intersection);
		float u = baryCoords[0];
		float v = baryCoords[1];
		float w = baryCoords[2];

		w = 1 - u - v;

		glm::vec3 baryPoint = glm::vec3(
			(v * intersection.intersectionPoint.x),
			(u * intersection.intersectionPoint.y),
			(w * intersection.intersectionPoint.z));


		for (int i = 0; i < limit; i++) {
			RayTriangleIntersection closestPointBetweenIntersectionAndLight;

			//if (phongShading) {
			//	closestPointBetweenIntersectionAndLight = getClosestIntersection(
			//		baryPoint,
			//		glm::normalize(glm::vec3(lightSource + offsets[i] - intersection.intersectionPoint)),
			//		OBJfile
			//	);
			//}
			//else {
			//	closestPointBetweenIntersectionAndLight = getClosestIntersection(
			//		intersection.intersectionPoint, ////////////// replace with interpolated value
			//		glm::normalize(glm::vec3(lightSource + offsets[i] - intersection.intersectionPoint)),
			//		OBJfile
			//	);
			//}

			closestPointBetweenIntersectionAndLight = getClosestIntersection(
				intersection.intersectionPoint,
				glm::normalize(glm::vec3(lightSource + offsets[i] - intersection.intersectionPoint)),
				OBJfile
			);
			double shadowRayIntersectionToPlaneIntersection = closestPointBetweenIntersectionAndLight.distanceFromCamera;
			//double lightSourceToPlaneIntersection = glm::length(lightSource + offsets[i] - intersection.intersectionPoint);
			double lightSourceToPlaneIntersection;

			if (phongShading) {
				lightSourceToPlaneIntersection = glm::length(lightSource + offsets[i] - baryPoint);
			}
			else {
				lightSourceToPlaneIntersection = glm::length(lightSource + offsets[i] - intersection.intersectionPoint);
			}

			//intersectionColour = intersection.intersectedTriangle.colour;

			if (closestPointBetweenIntersectionAndLight.rayHasHitTriangle && closestPointBetweenIntersectionAndLight.distanceFromCamera < lightSourceToPlaneIntersection) {
				continue;
			}
			else {
				intensity += CalculateLightIntensity(lightSourceToPlaneIntersection);
			}
		}

		float proxIntensity = intensity / float(limit);
		proxIntensity = clamp(proxIntensity, ambient, 1.0f);


		// ANGLE OF INCIDENCE LIGHTING
		ModelTriangle triangle = intersection.intersectedTriangle;

		//glm::vec3 normal = glm::cross((triangle.vertices[1] - triangle.vertices[0]), (triangle.vertices[2] - triangle.vertices[0]));
		glm::vec3 normal;

		if (phongShading) {
			normal = glm::normalize((u * vn1) + (v * vn2) + (w * vn0));
		}
		else {
			normal = glm::cross((triangle.vertices[1] - triangle.vertices[0]), (triangle.vertices[2] - triangle.vertices[0])); ////////////// replace with interpolated value
		}
		glm::vec3 incidentAngle = glm::normalize(glm::vec3(lightSource - intersection.intersectionPoint));

		float dot = glm::dot(glm::normalize(normal), incidentAngle);
		float incidenceIntensity = clamp(dot, ambient, 1.0f);


		//SPECULAR LIGHTING
		glm::vec3 rPerfectAngle = glm::normalize(glm::reflect(incidentAngle, normal));
		glm::vec3 viewDir = glm::normalize(glm::vec3(cameraPoint - intersection.intersectionPoint));

		float dotOfSpec = glm::dot(rPerfectAngle, incidentAngle);
		glm::vec3 reflectedAngle = glm::normalize(incidentAngle - 2.0f * normal * dotOfSpec);

		int exponent = 256;

		float specIntensity = pow(glm::dot(viewDir, reflectedAngle), exponent) /3;
		//specIntensity = clamp(specIntensity, 0.0f, 1.0f);

		intersectionColour = intersection.intersectedTriangle.colour;


		float relativeIntensity = proxIntensity * incidenceIntensity;
		//float relativeIntensity = proxIntensity * incidenceIntensity + (specIntensity * 0.0f);
		//float relativeIntensity = incidenceIntensity  + (specIntensity * 0.0f);

		//std::cout << "Proximity Intensity = " << proxIntensity << "\n";
		//std::cout << "Incidence Intensity = " << incidenceIntensity << "\n";
		//if (specIntensity > 0.2f) std::cout << "Specular Intensity = " << specIntensity << "\n";
		//std::cout << "Relative Intensity = " << relativeIntensity << "\n\n";

		if (intersection.intersectedTriangle.texturePoints[0].x > 0.0f) {

			//float pee0 = intersection.intersectedTriangle.texturePoints[0].x;
			//float poo0 = intersection.intersectedTriangle.texturePoints[0].y;
			//float pee1 = intersection.intersectedTriangle.texturePoints[1].x;
			//float poo1 = intersection.intersectedTriangle.texturePoints[1].y;
			//float pee2 = intersection.intersectedTriangle.texturePoints[2].x;
			//float poo2 = intersection.intersectedTriangle.texturePoints[2].y;

			glm::vec2 texPoint0 = glm::vec2(
				intersection.intersectedTriangle.texturePoints[0].x * textureMap.width, 
				intersection.intersectedTriangle.texturePoints[0].y * textureMap.height);
			glm::vec2 texPoint1 = glm::vec2(
				intersection.intersectedTriangle.texturePoints[1].x * textureMap.width,
				intersection.intersectedTriangle.texturePoints[1].y * textureMap.height);
			glm::vec2 texPoint2 = glm::vec2(
				intersection.intersectedTriangle.texturePoints[2].x * textureMap.width,
				intersection.intersectedTriangle.texturePoints[2].y * textureMap.height);

			
			//CanvasTriangle canvasTriangle = ModelTriangleToCanvasTriangle(
			//	intersection.intersectedTriangle,
			//	cameraPoint, focalLength, 40.0);
			////CanvasPoint convertedPoint = GetCanvasIntersectionPoint();

			//TexturePoint texPoint0 = CanvasPointToTexturePoint(canvasTriangle[0]);
			//TexturePoint texPoint1 = CanvasPointToTexturePoint(canvasTriangle[1]);
			//TexturePoint texPoint2 = CanvasPointToTexturePoint(canvasTriangle[2]);

			glm::vec3 barry = CalculateBarycentricCoordinatesNew(intersection.intersectionPoint, intersection.intersectedTriangle);

			float u = barry.x;
			float v = barry.y;
			float w = 1.0f - u - v;

			glm::vec2 result = u * texPoint0 + v * texPoint1 + w * texPoint2;

			int textureMapIndex = (std::round(result.x) + std::round(result.y) * textureMap.width);
			textureMapIndex = clamp(textureMapIndex, 0, 189600-1);
			uint32_t colour = textureMap.pixels[textureMapIndex];

			intersectionColour = uint32tToColour(colour);
			//return colour;

		}

		if (phongShading) {
			float relativeIntensity = incidenceIntensity + (specIntensity);
			outputColour = Colour(
				clamp(intersectionColour.red * relativeIntensity + (specIntensity * 255), 0.0f, 255.0f),
				clamp(intersectionColour.green * relativeIntensity + (specIntensity * 255), 0.0f, 255.0f),
				clamp(intersectionColour.blue * relativeIntensity + (specIntensity * 255), 0.0f, 255.0f));
		}
		else {
			float relativeIntensity = proxIntensity * incidenceIntensity + (specIntensity)+ambient;
			outputColour = Colour(
				clamp(intersectionColour.red * relativeIntensity + (specIntensity * 255), ambient, 255.0f),
				clamp(intersectionColour.green * relativeIntensity + (specIntensity * 255), ambient, 255.0f),
				clamp(intersectionColour.blue * relativeIntensity + (specIntensity * 255), ambient, 255.0f));
		}

		

		//float relativeIntensity = 1.0f;
		//outputColour = Colour(
		//	clamp(intersectionColour.red * relativeIntensity, ambient, 255.0f),
		//	clamp(intersectionColour.green * relativeIntensity, ambient, 255.0f),
		//	clamp(intersectionColour.blue * relativeIntensity, ambient, 255.0f));
		//outputColour = Colour(
		//	clamp(intersectionColour.red * relativeIntensity + (specIntensity * 255), ambient, 255.0f),
		//	clamp(intersectionColour.green * relativeIntensity + (specIntensity * 255), ambient, 255.0f),
		//	clamp(intersectionColour.blue * relativeIntensity + (specIntensity *255), ambient, 255.0f));


	}

	return ((255 << 24) + (outputColour.red << 16) + (outputColour.green << 8) + (outputColour.blue));
}



void DrawRandomStrokedTriangle(DrawingWindow& window) {
	CanvasTriangle triangle;
	Colour colour1;
	std::vector<CanvasPoint> v0, v1, v2;

	triangle = CanvasTriangle(CanvasPoint(rand() % WIDTH, rand() % HEIGHT),
		CanvasPoint(rand() % WIDTH, rand() % HEIGHT),
		CanvasPoint(rand() % WIDTH, rand() % HEIGHT));
	colour1 = Colour(rand() % 256, rand() % 256, rand() % 256);

	v0 = InterpolatePointToPoint(triangle[0], triangle[1]); //red
	v1 = InterpolatePointToPoint(triangle[1], triangle[2]); //green
	v2 = InterpolatePointToPoint(triangle[2], triangle[0]); //blue, spans vertical length of triangle

	v0.insert(v0.end(), v1.begin(), v1.end());
	v0.insert(v0.end(), v2.begin(), v2.end());


	// draw triangle ///////////////////////

	for (size_t i = 0; i < v0.size(); i++) {
		uint32_t colour = (255 << 24) + (int(colour1.red) << 16) + (int(colour1.green) << 8) + int(colour1.blue);
		window.setPixelColour(v0[i].x, v0[i].y, colour);
	}
}


void DrawRandomFilledTriangleOld(DrawingWindow& window) {
	CanvasPoint point0, point1, point2;
	CanvasTriangle triangle;
	Colour white, randColour;
	std::vector<CanvasPoint> v0, v1, v2, v2Top, v2Bottom, vOutline, vFill;

	white = Colour(255, 255, 255);
	randColour = Colour(rand() % 255, rand() % 255, rand() % 255);

	triangle = CanvasTriangle(CanvasPoint(rand() % WIDTH, rand() % HEIGHT),
		CanvasPoint(rand() % WIDTH, rand() % HEIGHT),
		CanvasPoint(rand() % WIDTH, rand() % HEIGHT));

	if (triangle[0].y > triangle[2].y) { std::swap(triangle[0], triangle[2]); }
	if (triangle[0].y > triangle[1].y) { std::swap(triangle[0], triangle[1]); }
	if (triangle[1].y > triangle[2].y) { std::swap(triangle[1], triangle[2]); }

	v0 = InterpolatePointToPoint(triangle[0], triangle[1]); //red
	v1 = InterpolatePointToPoint(triangle[1], triangle[2]); //green
	v2 = InterpolatePointToPoint(triangle[2], triangle[0]); //blue, spans vertical length of triangle

	vOutline = v0;
	vOutline.insert(vOutline.end(), v1.begin(), v1.end());
	vOutline.insert(vOutline.end(), v2.begin(), v2.end());

	// we need to find the horizontal line that connects triangle[1] and a point on v2

	for (size_t i = 0; i < v2.size(); i++) {
		if (v2[i].y == triangle[1].y) {
			v2Top = InterpolatePointToPoint(triangle[0], v2[i]);
			v2Bottom = InterpolatePointToPoint(v2[i], triangle[2]);
			break;
		}
	}
	// draw triangles ///////////////////////

	// top half of triangle
	for (size_t i = 0; i < v2Top.size(); i++) {
		vFill.clear();

		for (size_t j = 0; j < v0.size(); j++) {
			if (v0[j].y == v2Top[i].y) {
				vFill = InterpolatePointToPoint(v2Top[i], v0[j]);
				break;
			}
		}
		for (size_t j = 0; j < vFill.size(); j++) {
			uint32_t colour = (255 << 24) + (int(randColour.red) << 16) + (int(randColour.green) << 8) + int(randColour.blue);
			window.setPixelColour(vFill[j].x, vFill[j].y, colour);
		}
	}
	// bottom half of triangle
	for (size_t i = 0; i < v2Bottom.size(); i++) {
		vFill.clear();

		for (size_t j = 0; j < v1.size(); j++) {
			if (v1[j].y == v2Bottom[i].y) {
				vFill = InterpolatePointToPoint(v2Bottom[i], v1[j]);
				break;
			}
		}
		for (size_t j = 0; j < vFill.size(); j++) {
			uint32_t colour = (255 << 24) + (int(randColour.red) << 16) + (int(randColour.green) << 8) + int(randColour.blue);
			window.setPixelColour(vFill[j].x, vFill[j].y, colour);
		}
	}

	// draw triangle outline
	for (size_t i = 0; i < vOutline.size(); i++) {
		uint32_t colour = (255 << 24) + (int(white.red) << 16) + (int(white.green) << 8) + int(white.blue);
		window.setPixelColour(vOutline[i].x, vOutline[i].y, colour);
	}
}


void DrawRandomFilledTriangle(DrawingWindow& window) {
	std::vector<CanvasPoint> v0, v1, v2, topFill, bottomFill;;
	CanvasTriangle triangle;
	CanvasPoint midPoint;
	Colour fillColour;

	//triangle = CanvasTriangle(CanvasPoint(50, 400), CanvasPoint(500, 280), CanvasPoint(200, 20));
	triangle = CanvasTriangle(CanvasPoint(rand() % WIDTH, rand() % HEIGHT),
		CanvasPoint(rand() % WIDTH, rand() % HEIGHT),
		CanvasPoint(rand() % WIDTH, rand() % HEIGHT));

	fillColour = Colour(rand() % 255, rand() % 255, rand() % 255);

	// order triangle points by y axis
	if (triangle[0].y > triangle[2].y) { std::swap(triangle[0], triangle[2]); }
	if (triangle[0].y > triangle[1].y) { std::swap(triangle[0], triangle[1]); }
	if (triangle[1].y > triangle[2].y) { std::swap(triangle[1], triangle[2]); }

	// if triangle is already flat on top
	if (triangle[0].y == triangle[1].y) {
		bottomFill = drawBottomHalfTriangle(triangle[2], triangle[1], triangle[0], bottomFill);
	}

	// if triangle is already flat on bottom
	else if (triangle[1].y == triangle[2].y) {
		topFill = drawTopHalfTriangle(triangle[0], triangle[1], triangle[2], topFill);
	}

	else {
		// find midpoint of triangle opposite triangle[1]
		float gradient, midPointx;
		gradient = (triangle[1].y - triangle[0].y) / (triangle[2].y - triangle[0].y);
		midPointx = (gradient * (triangle[2].x - triangle[0].x)) + triangle[0].x;

		midPoint = CanvasPoint(midPointx, triangle[1].y);

		bottomFill = drawBottomHalfTriangle(triangle[2], triangle[1], midPoint, bottomFill);
		topFill = drawTopHalfTriangle(triangle[0], triangle[1], midPoint, topFill);
	}

	v0.insert(v0.end(), bottomFill.begin(), bottomFill.end());
	v0.insert(v0.end(), topFill.begin(), topFill.end());

	for (size_t i = 0; i < v0.size(); i++) {
		uint32_t colour = (255 << 24) + (int(fillColour.red) << 16) + (int(fillColour.green) << 8) + int(fillColour.blue);
		window.setPixelColour(v0[i].x, v0[i].y, colour);
	}
}


void DrawStrokedTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour lineColour) {
	std::vector<CanvasPoint> v0, v1, v2;

	v0 = InterpolatePointToPoint(triangle[0], triangle[1]);
	v1 = InterpolatePointToPoint(triangle[1], triangle[2]);
	v2 = InterpolatePointToPoint(triangle[2], triangle[0]);

	v0.insert(v0.end(), v1.begin(), v1.end());
	v0.insert(v0.end(), v2.begin(), v2.end());


	// draw triangle ///////////////////////

	for (size_t i = 0; i < v0.size(); i++) {
		if ((v0[i].x > 0 && v0[i].x < WIDTH - 1) && (v0[i].y > 0 && v0[i].y < HEIGHT - 1)) {
			uint32_t colour = (255 << 24) + (int(lineColour.red) << 16) + (int(lineColour.green) << 8) + int(lineColour.blue);
			window.setPixelColour(v0[i].x, v0[i].y, colour);
		}
	}
}


std::vector<std::vector<float>> DrawFilledTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour fillColour, std::vector<std::vector<float>> depthBuffer) {
	std::vector<CanvasPoint> v0, v1, v2, topFill, bottomFill;
	CanvasPoint midPoint;

	// order triangle points by y axis
	if (triangle[0].y > triangle[2].y) { std::swap(triangle[0], triangle[2]); }
	if (triangle[0].y > triangle[1].y) { std::swap(triangle[0], triangle[1]); }
	if (triangle[1].y > triangle[2].y) { std::swap(triangle[1], triangle[2]); }

	// if triangle is already flat on top
	if (std::roundf(triangle[0].y) == std::roundf(triangle[1].y)) {
		bottomFill = drawBottomHalfTriangle(triangle[2], triangle[1], triangle[0], bottomFill);
	}

	// if triangle is already flat on bottom
	else if (std::roundf(triangle[2].y) == std::roundf(triangle[1].y)) {
		topFill = drawTopHalfTriangle(triangle[0], triangle[1], triangle[2], topFill);
	}

	else {
		// find midpoint of triangle opposite triangle[1]
		float gradient = (triangle[1].y - triangle[0].y) / (triangle[2].y - triangle[0].y);
		float midPointx = (gradient * (triangle[2].x - triangle[0].x)) + triangle[0].x;

		midPoint = CanvasPoint(midPointx, triangle[1].y);

		bottomFill = drawBottomHalfTriangle(triangle[2], triangle[1], midPoint, bottomFill);
		topFill = drawTopHalfTriangle(triangle[0], triangle[1], midPoint, topFill);
	}

	v0.insert(v0.end(), bottomFill.begin(), bottomFill.end());
	v0.insert(v0.end(), topFill.begin(), topFill.end());


	for (size_t i = 0; i < v0.size(); i++) {

		if ((v0[i].x > 0 && v0[i].x < WIDTH - 1) && (v0[i].y > 0 && v0[i].y < HEIGHT - 1)) {

			if (v0[i].depth > depthBuffer[v0[i].x][v0[i].y]) {

				uint32_t colour = (255 << 24) + (int(fillColour.red) << 16) + (int(fillColour.green) << 8) + int(fillColour.blue);
				//uint32_t colour = (255 << 24) + (int(fillColour.red / v0[i].depth) << 16) + (int(fillColour.green / v0[i].depth) << 8) + int(fillColour.blue / v0[i].depth);
				//uint32_t colour = (255 << 24) + (int(fillColour.red * v0[i].depth) << 16) + (int(fillColour.green * v0[i].depth) << 8) + int(fillColour.blue * v0[i].depth);
				//uint32_t colour = (255 << 24) + (int(v0[i].depth * 255) << 16) + (int(v0[i].depth * 255) << 8) + int(v0[i].depth * 255);
				window.setPixelColour(v0[i].x, v0[i].y, colour);
				depthBuffer[v0[i].x][v0[i].y] = v0[i].depth;
			}
		}
	}

	return depthBuffer;
}


void drawSection(DrawingWindow& window, glm::vec3 cameraPoint, glm::mat3 cameraRotation, float focalLength, glm::vec3 lightPoint, std::vector<std::vector<ModelTriangle>>& OBJfile, int startY, int endY) {

	Colour intersectionColour;
	uint32_t colour;

	for (int y = startY; y < endY; y++) {
		for (int x = 0; x < WIDTH; x++) {
			if ((((x + y) & 1) > 0) && speedyMode) {
				continue;
			}

			colour = TraceLight(cameraPoint, cameraRotation, getRaysForScreen(x, y, focalLength), lightPoint, focalLength, OBJfile);

			//Colour colourFromTri = intersectionColour;
			//int redTemp = (colourFromTri.red);
			//int greenTemp = (colourFromTri.green);
			//int blueTemp = (colourFromTri.blue);
			//colour = (255 << 24) + (redTemp << 16) + (greenTemp << 8) + (blueTemp);

			//uint32_t colourFromTri = intersectionColour;


			window.setPixelColour(x, y, colour);
		}
	}
}



void drawFilledTriangleScene(DrawingWindow& window, glm::vec3 cameraPoint, glm::mat3 cameraRotation, float focalLength, std::vector<std::vector<ModelTriangle>>& OBJfile) {

	std::vector<std::vector<float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, 0));
	float scale = 40.0;

	// make screen black
	for (int y = 0; y < HEIGHT; y++) {
		for (int x = 0; x < WIDTH; x++) {
			uint32_t colour = (255 << 24) + (int(0) << 16) + (int(0) << 8) + int(0);
			window.setPixelColour(x, y, colour);
		}
	}

	for (size_t i = 0; i < OBJfile.size(); i++) {
		for (size_t j = 0; j < OBJfile[i].size(); j++) {
			CanvasTriangle canvasTriangle1 = ModelTriangleToCanvasTriangle(OBJfile[i][j], cameraPoint, cameraRotation, focalLength, scale);
			depthBuffer = DrawFilledTriangle(window, canvasTriangle1, OBJfile[i][j].colour, depthBuffer);
		}
	}
}


void drawWireframesScene(DrawingWindow& window, glm::vec3 cameraPoint, glm::mat3 cameraRotation, float focalLength, glm::vec3 lightPoint, std::vector<std::vector<ModelTriangle>>& OBJfile) {

	float scale = 40.0;

	// make screen black
	for (int y = 0; y < HEIGHT; y++) {
		for (int x = 0; x < WIDTH; x++) {
			uint32_t colour = (255 << 24) + (int(0) << 16) + (int(0) << 8) + int(0);
			window.setPixelColour(x, y, colour);
		}
	}

	for (size_t i = 0; i < OBJfile.size(); i++) {
		for (size_t j = 0; j < OBJfile[i].size(); j++) {

			CanvasTriangle canvasTriangle1 = ModelTriangleToCanvasTriangle(OBJfile[i][j], cameraPoint, cameraRotation, focalLength, scale);
			DrawStrokedTriangle(window, canvasTriangle1, OBJfile[i][j].colour);
		}
	}

	// ugly bad light ray code

	CanvasPoint lightPointOnScreen = GetCanvasIntersectionPoint(lightPoint, cameraPoint, cameraRotation, focalLength, scale);
	float pointyRaySize = 4.0;
	Colour white = Colour(255, 255, 255);

	ModelTriangle lightTriangleXAxis = ModelTriangle((lightPoint - glm::vec3(pointyRaySize, 0.0, 0.0)),
		(lightPoint - glm::vec3(pointyRaySize, 0.0, 0.0)),
		(lightPoint + glm::vec3(pointyRaySize, 0.0, 0.0)), white);

	ModelTriangle lightTriangleYAxis = ModelTriangle((lightPoint - glm::vec3(0.0, pointyRaySize, 0.0)),
		(lightPoint - glm::vec3(0.0, pointyRaySize, 0.0)),
		(lightPoint + glm::vec3(0.0, pointyRaySize, 0.0)), white);

	ModelTriangle lightTriangleZAxis = ModelTriangle((lightPoint - glm::vec3(0.0, 0.0, pointyRaySize)),
		(lightPoint - glm::vec3(0.0, 0.0, pointyRaySize)),
		(lightPoint + glm::vec3(0.0, 0.0, pointyRaySize)), white);

	CanvasTriangle lightTriangleXAxisCanvas = ModelTriangleToCanvasTriangle(lightTriangleXAxis, cameraPoint, cameraRotation, focalLength, scale);
	CanvasTriangle lightTriangleYAxisCanvas = ModelTriangleToCanvasTriangle(lightTriangleYAxis, cameraPoint, cameraRotation, focalLength, scale);
	CanvasTriangle lightTriangleZAxisCanvas = ModelTriangleToCanvasTriangle(lightTriangleZAxis, cameraPoint, cameraRotation, focalLength, scale);

	DrawStrokedTriangle(window, lightTriangleXAxisCanvas, white);
	DrawStrokedTriangle(window, lightTriangleYAxisCanvas, white);
	DrawStrokedTriangle(window, lightTriangleZAxisCanvas, white);
}



void draw(DrawingWindow& window, glm::vec3 cameraPoint, glm::mat3 cameraRotation, glm::vec3 lightPoint, float focalLength, std::vector<std::vector<ModelTriangle>>& OBJfile, int renderMode) {
	if (renderMode == 1) {
		int step = HEIGHT / noOfCores;
		std::thread* workers[noOfCores];
		for (int i = 0; i < noOfCores; i++)
			workers[i] = new std::thread(drawSection, std::ref(window), cameraPoint, cameraRotation, focalLength, lightPoint, std::ref(OBJfile), step * i, step * (i + 1));
		for (int i = 0; i < noOfCores; i++) {
			workers[i]->join();
			delete workers[i];
		}
		//drawSection(window, cameraPoint, focalLength, lightPoint, OBJfile, 0, HEIGHT);
	}
	else if (renderMode == 2) {
		drawFilledTriangleScene(window, cameraPoint, cameraRotation, focalLength, OBJfile);
	}
	else if (renderMode == 3) {
		drawWireframesScene(window, cameraPoint, cameraRotation, focalLength, lightPoint, OBJfile);
	}
}




void handleEvent(SDL_Event event, DrawingWindow& window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_u) DrawRandomStrokedTriangle(window);
		else if (event.key.keysym.sym == SDLK_f) DrawRandomFilledTriangle(window);
		else if (event.key.keysym.sym == SDLK_0) {
			speedyMode = !speedyMode;
			if (speedyMode) std::cout << "Speedy mode on\n";
			else std::cout << "Speedy mode off\n";
		}
		else if (event.key.keysym.sym == SDLK_p) {
			std::cout << "snap!\n";
			std::string noToFile;

			//if (outputCounter < 10) noToFile = "0000" + std::to_string(outputCounter);
			//if (outputCounter < 100) noToFile = "000" + std::to_string(outputCounter);
			//if (outputCounter < 1000) noToFile = "00" + std::to_string(outputCounter);
			//if (outputCounter < 10000) noToFile = "0" + std::to_string(outputCounter);
			//if (outputCounter < 100000) noToFile = std::to_string(outputCounter);

			//int dd = 1, mm = 9, yy = 1;
			//printf("%02d - %02d - %04d", mm, dd, y6ty);
			//noToFile = ("%06d", outputCounter);


			//std::string filenameppm = "output" + std::to_string(outputCounter) + ".ppm";
			//std::string filenamebmp = "output" + std::to_string(outputCounter) + ".bmp";

			//window.savePPM(filenameppm);
			//window.saveBMP(filenamebmp);

			//window.savePPM("output.ppm");
			//window.saveBMP("output.bmp");

			//outputCounter++;
		}
	}
	//else if (event.type == SDL_MOUSEBUTTONDOWN) {
	//	window.savePPM("output.ppm");
	//	window.saveBMP("output.bmp");
	//}
}


glm::vec3 checkMovementEvent(SDL_Event event, glm::vec3 cameraPoint) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_a) {
			cameraPoint.x = cameraPoint.x - 0.5; // move left
			//std::cout << "cameraPoint = ";
			OutputGLMVEC3(cameraPoint);
		}
		else if (event.key.keysym.sym == SDLK_d) {
			cameraPoint.x = cameraPoint.x + 0.5; // move right
			//std::cout << "cameraPoint = ";
			OutputGLMVEC3(cameraPoint);
		}
		else if (event.key.keysym.sym == SDLK_w) {
			cameraPoint.z = cameraPoint.z - 0.5; // move forwards
			//std::cout << "cameraPoint = ";
			OutputGLMVEC3(cameraPoint);
		}
		else if (event.key.keysym.sym == SDLK_s) {
			cameraPoint.z = cameraPoint.z + 0.5; // move backwards
			//std::cout << "cameraPoint = ";
			OutputGLMVEC3(cameraPoint);
		}
		else if (event.key.keysym.sym == SDLK_e) {
			cameraPoint.y = cameraPoint.y + 0.5; // move up
			//std::cout << "cameraPoint = ";
			OutputGLMVEC3(cameraPoint);
		}
		else if (event.key.keysym.sym == SDLK_q) {
			cameraPoint.y = cameraPoint.y - 0.5; // move down
			//std::cout << "cameraPoint = ";
			OutputGLMVEC3(cameraPoint);
		}
	}
	return cameraPoint;
}

glm::mat3 checkRotationEvent(SDL_Event event, glm::mat3 cameraRotation, int renderMode) {
	float theta = 0.05f;
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_z) { // rotation around x axis
			cameraRotation *= glm::mat3(
				1.0,	0.0,			0.0,
				0.0,	cos(theta),		-sin(theta),
				0.0,	sin(theta),		cos(theta)
			);
		}
		else if (event.key.keysym.sym == SDLK_x) { // rotation around x axis
			cameraRotation *= glm::mat3(
				1.0,	0.0,			0.0,
				0.0,	cos(theta),		sin(theta),
				0.0,	-sin(theta),	cos(theta)
			);
		}
		else if (event.key.keysym.sym == SDLK_v) { // rotation around y axis
			cameraRotation *= glm::mat3(
				cos(theta),		0.0,	sin(theta),
				0.0,			1.0,	0.0,
				-sin(theta),	0.0,	cos(theta)
			);
		}
		else if (event.key.keysym.sym == SDLK_c) { // rotation around y axis
			cameraRotation *= glm::mat3(
				cos(theta),		0.0,	-sin(theta),
				0.0,			1.0,	0.0,
				sin(theta),		0.0,	cos(theta)
			);
		}
		else if (event.key.keysym.sym == SDLK_b) { // rotation around z axis
			cameraRotation *= glm::mat3(
				cos(theta),		sin(theta),		0.0,
				-sin(theta),	cos(theta),		0.0,
				0.0,			0.0,			1.0
			);
		}
		else if (event.key.keysym.sym == SDLK_n) { // rotation around z axis
			cameraRotation *= glm::mat3(
				cos(theta),		-sin(theta),	0.0,
				sin(theta),		cos(theta),		0.0,
				0.0,			0.0,			1.0
			);
		}
	}
	return cameraRotation;
}

glm::vec3 checkLightPosition(SDL_Event event, glm::vec3 lightPoint) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) {
			lightPoint.x = lightPoint.x - 0.5; // move left
			std::cout << "lightPoint = ";
			OutputGLMVEC3(lightPoint);
		}
		else if (event.key.keysym.sym == SDLK_RIGHT) {
			lightPoint.x = lightPoint.x + 0.5; // move right
			std::cout << "lightPoint = ";
			OutputGLMVEC3(lightPoint);
		}
		else if (event.key.keysym.sym == SDLK_UP) {
			lightPoint.z = lightPoint.z - 0.5; // move forwards
			std::cout << "lightPoint = ";
			OutputGLMVEC3(lightPoint);
		}
		else if (event.key.keysym.sym == SDLK_DOWN) {
			lightPoint.z = lightPoint.z + 0.5; // move backwards
			std::cout << "lightPoint = ";
			OutputGLMVEC3(lightPoint);
		}
		else if (event.key.keysym.sym == SDLK_COMMA) {
			lightPoint.y = lightPoint.y + 0.5; // move up
			std::cout << "lightPoint = ";
			OutputGLMVEC3(lightPoint);
		}
		else if (event.key.keysym.sym == SDLK_PERIOD) {
			lightPoint.y = lightPoint.y - 0.5; // move down
			std::cout << "lightPoint = ";
			OutputGLMVEC3(lightPoint);
		}
	}
	return lightPoint;
}

int checkRenderMode(SDL_Event event, int renderMode) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_1) {
			renderMode = 1; // rasterised
			std::cout << "Rasterised mode on\n";
		}
		else if (event.key.keysym.sym == SDLK_2) {
			renderMode = 2; // filled triangles
			std::cout << "Filled triangle mode on\n";
		}
		else if (event.key.keysym.sym == SDLK_3) {
			renderMode = 3; // wireframe triangles
			std::cout << "Wireframe mode on\n";
		}
	}
	return renderMode;
}


glm::mat3 LookAt(glm::mat3 cameraRotation, glm::vec3 cameraPoint, glm::vec3 focus) {

	glm::vec3 forward = glm::normalize(glm::vec3(cameraPoint - focus));
	glm::vec3 right = glm::cross(forward, glm::normalize(glm::vec3(0.0,1.0,0.0)));
	glm::vec3 vertical = glm::cross(forward, right);

	cameraRotation = glm::mat3(
		right, -vertical, forward
	);

	return cameraRotation;
}


int main(int argc, char* argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	glm::vec3 cameraPoint(0.0, 0.0, 8.0);
	//glm::vec3 cameraPoint(-2.0, -2.0, 5.0);
	//glm::vec3 cameraPoint(0.0, 1.2, 5.0);
	
	glm::mat3 cameraRotation = glm::mat3(
		1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0
	);

	//glm::vec3 lightPoint(0.0, 1.5, 2.0);
	glm::vec3 lightPoint(0.0, 2.4, 0.0);
	//glm::vec3 lightPoint(0.2, 2.2, 2.0);

	float focalLength = 10.0;

	int renderMode = 2; // 1 is fully rasterised, 2 is with filled triangles, 3 is with coloured wireframes

	std::vector<std::vector<ModelTriangle>> OBJfile = loadOBJFile("cornell-box");
	//std::vector<std::vector<ModelTriangle>> OBJfile = loadOBJFile("sphere");
	//std::vector<std::vector<ModelTriangle>> OBJfile = loadOBJFile("textured-cornell-box");
	
	//std::cout << textureMap.pixels[50] << "\n";

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		// INPUT
		if (window.pollForInputEvents(event)) {
			handleEvent(event, window);
			cameraPoint = checkMovementEvent(event, cameraPoint);
			if (!lookAtActivated) {
				cameraRotation = checkRotationEvent(event, cameraRotation, renderMode);
			}
			else {
				cameraRotation = LookAt(cameraRotation, cameraPoint, glm::vec3(0.0,0.0,0.0));
			}
			lightPoint = checkLightPosition(event, lightPoint);
			renderMode = checkRenderMode(event, renderMode);
		}

		// PROCESS
		draw(window, cameraPoint, cameraRotation, lightPoint, focalLength, OBJfile, renderMode);


		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		// OUTPUT
		window.renderFrame();
	}

	return 0;
}


int poo(int argc, char* argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	//glm::vec3 lightPoint(0.0, 2.4, 0.0);

	float focalLength = 10.0;
	std::vector<std::vector<ModelTriangle>> OBJfile = loadOBJFile("sphere");
	int renderMode = 1;

	//float diff = 0.5f;
	//std::vector<float> xdiff = {
	//	0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 3.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5,

	//	8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0,
	//	8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0,

	//	7.5, 7.0, 6.5, 6.0, 5.5, 5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0,
	//	-0.5, -1.0, -1.5, -2.0, -2.5, -3.0, -3.5, -3.0, -4.5, -5.0, -5.5, -6.0, -6.5, -7.0, -7.5, -8.0,

	//	-8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0,
	//	-8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0,

	//	-7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0
	//};

	//int n = xdiff.size();

	///

	//std::vector<float> zdiff = {
	//	8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0,

	//	7.5, 7.0, 6.5, 6.0, 5.5, 5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0,
	//	-0.5, -1.0, -1.5, -2.0, -2.5, -3.0, -3.5, -3.0, -4.5, -5.0, -5.5, -6.0, -6.5, -7.0, -7.5, -8.0,

	//	-8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0,
	//	-8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0,

	//	-7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5,
	//	0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 3.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0,

	//	8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0
	//};

	//n = zdiff.size();

	///

	std::vector<glm::vec3> camPositionsFullRotate = {
		glm::vec3(0.5, 0, 8),glm::vec3(1, 0, 8),glm::vec3(1.5, 0, 8),glm::vec3(2, 0, 8),
		glm::vec3(2.5, 0, 8),glm::vec3(3, 0, 8),glm::vec3(3.5, 0, 8),glm::vec3(4, 0, 8),
		glm::vec3(4.5, 0, 8),glm::vec3(5, 0, 8),glm::vec3(5.5, 0, 8),glm::vec3(6, 0, 8),
		glm::vec3(6.5, 0, 8),glm::vec3(7, 0, 8),glm::vec3(7.5, 0, 8),glm::vec3(8, 0, 8),
		glm::vec3(8, 0, 7.5),glm::vec3(8, 0, 7),glm::vec3(8, 0, 6.5),glm::vec3(8, 0, 6),
		glm::vec3(8, 0, 5.5),glm::vec3(8, 0, 5),glm::vec3(8, 0, 4.5),glm::vec3(8, 0, 4),
		glm::vec3(8, 0, 3.5),glm::vec3(8, 0, 3),glm::vec3(8, 0, 2.5),glm::vec3(8, 0, 2),
		glm::vec3(8, 0, 1.5),glm::vec3(8, 0, 1),glm::vec3(8, 0, 0.5),glm::vec3(8, 0, 0),
		glm::vec3(8, 0, -0.5),glm::vec3(8, 0, -1),glm::vec3(8, 0, -1.5),glm::vec3(8, 0, -2),
		glm::vec3(8, 0, -2.5),glm::vec3(8, 0, -3),glm::vec3(8, 0, -3.5),glm::vec3(8, 0, -4),
		glm::vec3(8, 0, -4.5),glm::vec3(8, 0, -5),glm::vec3(8, 0, -5.5),glm::vec3(8, 0, -6),
		glm::vec3(8, 0, -6.5),glm::vec3(8, 0, -7),glm::vec3(8, 0, -7.5),glm::vec3(8, 0, -8),
		glm::vec3(7.5, 0, -8),glm::vec3(7, 0, -8),glm::vec3(6.5, 0, -8),glm::vec3(6, 0, -8),
		glm::vec3(5.5, 0, -8),glm::vec3(5, 0, -8),glm::vec3(4.5, 0, -8),glm::vec3(4, 0, -8),
		glm::vec3(3.5, 0, -8),glm::vec3(3, 0, -8),glm::vec3(2.5, 0, -8),glm::vec3(2, 0, -8),
		glm::vec3(1.5, 0, -8),glm::vec3(1, 0, -8),glm::vec3(0.5, 0, -8),glm::vec3(0, 0, -8),
		glm::vec3(-0.5, 0, -8),glm::vec3(-1, 0, -8),glm::vec3(-1.5, 0, -8),glm::vec3(-2, 0, -8),
		glm::vec3(-2.5, 0, -8),glm::vec3(-3, 0, -8),glm::vec3(-3.5, 0, -8),glm::vec3(-4, 0, -8),
		glm::vec3(-4.5, 0, -8),glm::vec3(-5, 0, -8),glm::vec3(-5.5, 0, -8),glm::vec3(-6, 0, -8),
		glm::vec3(-6.5, 0, -8),glm::vec3(-7, 0, -8),glm::vec3(-7.5, 0, -8),glm::vec3(-8, 0, -8),
		glm::vec3(-8, 0, -7.5),glm::vec3(-8, 0, -7),glm::vec3(-8, 0, -6.5),glm::vec3(-8, 0, -6),
		glm::vec3(-8, 0, -5.5),glm::vec3(-8, 0, -5),glm::vec3(-8, 0, -4.5),glm::vec3(-8, 0, -4),
		glm::vec3(-8, 0, -3.5),glm::vec3(-8, 0, -3),glm::vec3(-8, 0, -2.5),glm::vec3(-8, 0, -2),
		glm::vec3(-8, 0, -1.5),glm::vec3(-8, 0, -1),glm::vec3(-8, 0, -0.5),glm::vec3(-8, 0, 0),
		glm::vec3(-8, 0, 0.5),glm::vec3(-8, 0, 1),glm::vec3(-8, 0, 1.5),glm::vec3(-8, 0, 2),
		glm::vec3(-8, 0, 2.5),glm::vec3(-8, 0, 3),glm::vec3(-8, 0, 3.5),glm::vec3(-8, 0, 4),
		glm::vec3(-8, 0, 4.5),glm::vec3(-8, 0, 5),glm::vec3(-8, 0, 5.5),glm::vec3(-8, 0, 6),
		glm::vec3(-8, 0, 6.5),glm::vec3(-8, 0, 7),glm::vec3(-8, 0, 7.5),glm::vec3(-8, 0, 8),
		glm::vec3(-7.5, 0, 8),glm::vec3(-7, 0, 8),glm::vec3(-6.5, 0, 8),glm::vec3(-6, 0, 8),
		glm::vec3(-5.5, 0, 8),glm::vec3(-5, 0, 8),glm::vec3(-4.5, 0, 8),glm::vec3(-4, 0, 8),
		glm::vec3(-3.5, 0, 8),glm::vec3(-3, 0, 8),glm::vec3(-2.5, 0, 8),glm::vec3(-2, 0, 8),
		glm::vec3(-1.5, 0, 8),glm::vec3(-1, 0, 8),glm::vec3(-0.5, 0, 8),glm::vec3(0, 0, 8)
	};

	std::vector<glm::vec3> camPositionsBitRotate = {
		glm::vec3(0.5, 0, 8),glm::vec3(1, 0, 8),glm::vec3(1.5, 0, 8),glm::vec3(2, 0, 8),
		glm::vec3(2.5, 0, 8),glm::vec3(3, 0, 8),glm::vec3(3.5, 0, 8),glm::vec3(4, 0, 8),
		glm::vec3(3.5, 0, 8),glm::vec3(3, 0, 8),glm::vec3(2.5, 0, 8),glm::vec3(2, 0, 8),
		glm::vec3(1.5, 0, 8),glm::vec3(1, 0, 8),glm::vec3(0.5, 0, 8),glm::vec3(0, 0, 8),

		glm::vec3(-0.5, 0, 8),glm::vec3(-1, 0, 8),glm::vec3(-1.5, 0, 8),glm::vec3(-2, 0, 8),
		glm::vec3(-2.5, 0, 8),glm::vec3(-3, 0, 8),glm::vec3(-3.5, 0, 8),glm::vec3(-4, 0, 8),
		glm::vec3(-3.5, 0, 8),glm::vec3(-3, 0, 8),glm::vec3(-2.5, 0, 8),glm::vec3(-2, 0, 8),
		glm::vec3(-1.5, 0, 8),glm::vec3(-1, 0, 8),glm::vec3(-0.5, 0, 8),glm::vec3(0, 0, 8)
	};

	std::vector<glm::vec3> camPositionsRim = { glm::vec3(0.5, 0, 8)
	,glm::vec3(1, 0, 8)
	,glm::vec3(1.5, 0, 8)
	,glm::vec3(2, 0, 8)
	,glm::vec3(2.5, 0, 8)
	,glm::vec3(3, 0, 8)
	,glm::vec3(3.5, 0, 8)
	,glm::vec3(4, 0, 8)
	,glm::vec3(4, 0.5, 8)
	,glm::vec3(4, 1, 8)
	,glm::vec3(4, 1.5, 8)
	,glm::vec3(4, 2, 8)
	,glm::vec3(4, 2.5, 8)
	,glm::vec3(4, 3, 8)
	,glm::vec3(4, 3.5, 8)
	,glm::vec3(4, 4, 8)
	,glm::vec3(3.5, 4, 8)
	,glm::vec3(3, 4, 8)
	,glm::vec3(2.5, 4, 8)
	,glm::vec3(2, 4, 8)
	,glm::vec3(1.5, 4, 8)
	,glm::vec3(1, 4, 8)
	,glm::vec3(0.5, 4, 8)
	,glm::vec3(0, 4, 8)
	,glm::vec3(-0.5, 4, 8)
	,glm::vec3(-1, 4, 8)
	,glm::vec3(-1.5, 4, 8)
	,glm::vec3(-2, 4, 8)
	,glm::vec3(-2.5, 4, 8)
	,glm::vec3(-3, 4, 8)
	,glm::vec3(-3.5, 4, 8)
	,glm::vec3(-4, 4, 8)
	,glm::vec3(-4, 3.5, 8)
	,glm::vec3(-4, 3, 8)
	,glm::vec3(-4, 2.5, 8)
	,glm::vec3(-4, 2, 8)
	,glm::vec3(-4, 1.5, 8)
	,glm::vec3(-4, 1, 8)
	,glm::vec3(-4, 0.5, 8)
	,glm::vec3(-4, 0, 8)
	,glm::vec3(-4, -0.5, 8)
	,glm::vec3(-4, -1, 8)
	,glm::vec3(-4, -1.5, 8)
	,glm::vec3(-4, -2, 8)
	,glm::vec3(-4, -2.5, 8)
	,glm::vec3(-4, -3, 8)
	,glm::vec3(-4, -3.5, 8)
	,glm::vec3(-4, -4, 8)
	,glm::vec3(-3.5, -4, 8)
	,glm::vec3(-3, -4, 8)
	,glm::vec3(-2.5, -4, 8)
	,glm::vec3(-2, -4, 8)
	,glm::vec3(-1.5, -4, 8)
	,glm::vec3(-1, -4, 8)
	,glm::vec3(-0.5, -4, 8)
	,glm::vec3(0, -4, 8)
	,glm::vec3(0.5, -4, 8)
	,glm::vec3(1, -4, 8)
	,glm::vec3(1.5, -4, 8)
	,glm::vec3(2, -4, 8)
	,glm::vec3(2.5, -4, 8)
	,glm::vec3(3, -4, 8)
	,glm::vec3(3.5, -4, 8)
	,glm::vec3(4, -4, 8)
	,glm::vec3(4, -3.5, 8)
	,glm::vec3(4, -3, 8)
	,glm::vec3(4, -2.5, 8)
	,glm::vec3(4, -2, 8)
	,glm::vec3(4, -1.5, 8)
	,glm::vec3(4, -1, 8)
	,glm::vec3(4, -0.5, 8)
	,glm::vec3(4, 0, 8)
	,glm::vec3(3.5, 0, 8)
	,glm::vec3(3, 0, 8)
	,glm::vec3(2.5, 0, 8)
	,glm::vec3(2, 0, 8)
	,glm::vec3(1.5, 0, 8)
	,glm::vec3(1, 0, 8)
	,glm::vec3(0.5, 0, 8)
	,glm::vec3(0, 0, 8) };

	glm::mat3 cameraRotation = glm::mat3(
		1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0
	);

	int outputCounter = 0;

	std::vector<glm::vec3> ballRotate = {
	glm::vec3(0.0, 1.2, 5.0),
	glm::vec3(0.33, 1.2, 5.0),
	glm::vec3(0.67, 1.2, 5.0),
	glm::vec3(1.0, 1.2, 5.0),
	glm::vec3(1.33, 1.2, 5.0),
	glm::vec3(1.67, 1.2, 5.0),
	glm::vec3(2.0, 1.2, 5.0),
	glm::vec3(-0.33, 1.2, 5.0),
	glm::vec3(-0.67, 1.2, 5.0),
	glm::vec3(-1.0, 1.2, 5.0),
	glm::vec3(-1.33, 1.2, 5.0),
	glm::vec3(-1.67, 1.2, 5.0),
	glm::vec3(-2.0, 1.2, 5.0)
	};

	//glm::vec3 cameraPoint(0.0, 1.2, 5.0);
	//cameraPoint = glm::vec3(-0.33, 1.2, 5.0);j
	//cameraPoint = glm::vec3(-0.67, 1.2, 5.0);
	//cameraPoint = glm::vec3(-1.0, 1.2, 5.0);
	//cameraPoint = glm::vec3(-1.33, 1.2, 5.0);
	//cameraPoint = glm::vec3(-1.67, 1.2, 5.0);
	//cameraPoint = glm::vec3(-2.0, 1.2, 5.0);

	//glm::vec3 cameraPoint(0.0, 0.0, 8.0);
	glm::vec3 cameraPoint(0.0, 1.2, 5.0);

	glm::vec3 lightPoint(0.4, 2.2, 2.0);
	glm::vec3 lightPoint1(-0.4, 2.2, 2.0);


	std::vector<glm::vec3> v0 = InterpolateVec3ToVec3(lightPoint, lightPoint1, 16);


	//glm::vec3 cameraPoint = v0[0];

	cameraRotation = LookAt(cameraRotation, cameraPoint, glm::vec3(0.0, 0.0, 0.0));
	draw(window, cameraPoint, cameraRotation, v0[0], focalLength, OBJfile, renderMode);
	std::string filenamebmp = "output" + std::to_string(outputCounter) + ".bmp";
	window.saveBMP(filenamebmp);
	outputCounter++;


	for (int i = 1; i < v0.size(); i++) {
		//cameraPoint = glm::vec3(v0[i]);
		cameraRotation = LookAt(cameraRotation, cameraPoint, glm::vec3(0.0, 0.0, 0.0));
		draw(window, cameraPoint, cameraRotation, v0[i], focalLength, OBJfile, renderMode);
		filenamebmp = "output" + std::to_string(outputCounter) + ".bmp";
		window.saveBMP(filenamebmp);
		outputCounter++;
	}

	return 0;
}