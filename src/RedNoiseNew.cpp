#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <algorithm>
#include <Colour.h>
#include <typeinfo>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/vec3.hpp> // glm::vec3
#include <stdio.h>


#define WIDTH 640
#define HEIGHT 480


//POINT TO POINT INTERPOLATION
std::vector<CanvasPoint> interpolatePointToPoint(CanvasPoint from, CanvasPoint to) {
    std::vector<CanvasPoint> v;
    CanvasPoint counter = from;

    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
    float xStepSize = xDiff/numberOfSteps;
    float yStepSize = yDiff/numberOfSteps;

    for (float i = 0.0; i < numberOfSteps; i++) {
        v.push_back(counter);
        counter.x = int(from.x + (xStepSize * i));
        counter.y = int(from.y + (yStepSize * i));
    }
    return v;
}


// POINT TO POINT TO POINT INTERPOLATION for lovely triangles
std::vector<CanvasPoint> interpolateTrianglePoints(CanvasTriangle triangle) {
    std::vector<CanvasPoint> v, v1, v2, v3;
    v1 = interpolatePointToPoint(triangle[0], triangle[1]);
    v2 = interpolatePointToPoint(triangle[1], triangle[2]);
    v3 = interpolatePointToPoint(triangle[2], triangle[0]);
    v = v1;
    v.insert(v.end(), v2.begin(), v2.end());
    v.insert(v.end(), v3.begin(), v3.end());
    return v;
}



// OUTLINED TRIANGLE DRAWER
void DrawOutlinedTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour lineColour) {
	std::vector<CanvasPoint> v;
	v = interpolateTrianglePoints(triangle);

//	for (size_t i = 0; i < v.size(); i++) {
//		std::cout << v[i] << "\n";
//	}


	float yMin = std::min(abs(triangle[0].y), abs(triangle[1].y));
	yMin = std::min(abs(yMin), abs(triangle[2].y));
	float yMax = std::max(abs(triangle[0].y), abs(triangle[1].y));
	yMax = std::max(abs(yMax), abs(triangle[2].y));
	float xMin = std::min(abs(triangle[0].x), abs(triangle[1].x));
	xMin = std::min(abs(xMin), abs(triangle[2].x));
	float xMax = std::max(abs(triangle[0].x), abs(triangle[1].x));
	xMax = std::max(abs(xMax), abs(triangle[2].x));

	float red = float(lineColour.red);
	float green = float(lineColour.green);
	float blue = float(lineColour.blue);

	for (size_t y = yMin; y < yMax + 1; y++) { // move along y axis

		for (size_t x = xMin; x < xMax + 1; x++) { // move along x axis

			for (size_t i = 0; i < v.size(); i++) { // scan vector of interpolated values

				if ((x == v[i].x) && (y == v[i].y)) { // if current x and y are in results

					uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
					window.setPixelColour(x, y, colour);
				}
			}
		}
	}
}


// FILLED TRIANGLE DRAWER
void DrawFilledTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour lineColour, Colour fillColour) {
	std::vector<CanvasPoint> v;
	v = interpolateTrianglePoints(triangle);

//	for (size_t i = 0; i < v.size(); i++) {
//		std::cout << v[i] << "\n";
//	}

	float yMin = std::min(abs(triangle[0].y), abs(triangle[1].y));
	yMin = std::min(abs(yMin), abs(triangle[2].y));
	float yMax = std::max(abs(triangle[0].y), abs(triangle[1].y));
	yMax = std::max(abs(yMax), abs(triangle[2].y));
	float xMin = std::min(abs(triangle[0].x), abs(triangle[1].x));
	xMin = std::min(abs(xMin), abs(triangle[2].x));
	float xMax = std::max(abs(triangle[0].x), abs(triangle[1].x));
	xMax = std::max(abs(xMax), abs(triangle[2].x));


	float red = float(fillColour.red);
	float green = float(fillColour.green);
	float blue = float(fillColour.blue);

	for (size_t y = yMin; y < yMax + 1; y++) { // move along y axis

		for (size_t x = xMin; x < xMax + 1; x++) { // move along x axis

			for (size_t i = 0; i < v.size(); i++) { // scan vector of interpolated values

				if ((x == v[i].x) && (y == v[i].y)) { // if current x and y are in results

					uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
					window.setPixelColour(x, y, colour);
				}
			}
		}
	}
}

void CreateOutlinedTriangle(DrawingWindow& window) {
	CanvasPoint point1, point2, point3;
	CanvasTriangle triangle1;
	Colour colour1;

	point1 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
	point2 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
	point3 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
	triangle1 = CanvasTriangle(point1, point2, point3);
	colour1 = Colour(rand() % 256, rand() % 256, rand() % 256);

	DrawOutlinedTriangle(window, triangle1, colour1);
}

void CreateFilledTriangle(DrawingWindow& window) {
	CanvasPoint point1, point2, point3;
	CanvasTriangle triangle1;
	Colour outlineColour, fillColour;

	point1 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
	point2 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
	point3 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
	triangle1 = CanvasTriangle(point1, point2, point3);
	outlineColour = Colour(255,255,255);
	fillColour = Colour(rand() % 256, rand() % 256, rand() % 256);

	DrawFilledTriangle(window, triangle1, outlineColour, fillColour);
}



void draw5(DrawingWindow& window) {
	window.clearPixels(); 
//	CreateTriangle(window);
}




void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;

		else if (event.key.keysym.sym == SDLK_u) CreateOutlinedTriangle(window);
		else if (event.key.keysym.sym == SDLK_f) CreateFilledTriangle(window);

	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {

//    std::vector<glm::vec3> result;
//    glm::vec3 test1(150, 69, 125);
//    glm::vec3 test2(244, 3, 235);
//
//    result = interpolateVec3(test1, test2, 7);
//
//    for(size_t i = 0; i < result.size(); i++) {
//        printf("%f, %f, %f\n", result[i].x, result[i].y, result[i].z);
//    }

	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	draw5(window);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
//		draw2(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
