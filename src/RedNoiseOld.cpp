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


// 1D INTERPOLATION
std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
    std::vector<float> v;
    float counter = from;
    for (float i = 0; i < numberOfValues; i++){
        v.push_back(counter);
        counter = counter + ((to-from)/float(numberOfValues));
    }
    return v;
}

//3D INTERPOLATION
std::vector<glm::vec3> interpolateVec3(glm::vec3 from, glm::vec3 to, int numberOfValues){
    std::vector<glm::vec3> v;
    glm::vec3 counter = from;
    for (float i = 0; i < numberOfValues; i++){
        v.push_back(counter);
		for (float i = 0; i < 3; i++) {counter[i] = counter[i] + ((to[i] - from[i]) / float(numberOfValues));}
    }
    return v;
}

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

//POINT TO POINT INTERPOLATION but you can add in the vector you wanna add numbers to or whatever
std::vector<CanvasPoint> interpolatePointToPointWithVector(std::vector<CanvasPoint> v, CanvasPoint from, CanvasPoint to) {
	/*std::vector<CanvasPoint> v;*/
	CanvasPoint counter = from;

	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff / numberOfSteps;
	float yStepSize = yDiff / numberOfSteps;

	for (float i = 0.0; i < numberOfSteps; i++) {
		v.push_back(counter);
		counter.x = int(from.x + (xStepSize * i));
		counter.y = int(from.y + (yStepSize * i));
	}
	return v;
}

//3 POINT INTERPOLATION
std::vector<CanvasPoint> interpolateTriangle(std::vector<CanvasPoint> v, CanvasTriangle triangle) {
	/*std::vector<CanvasPoint> v;*/

	CanvasPoint counter = triangle[0];
	float xDiff = triangle[0].x - triangle[1].x;
	float yDiff = triangle[0].y - triangle[1].y;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff / numberOfSteps;
	float yStepSize = yDiff / numberOfSteps;

	for (float i = 0.0; i < numberOfSteps; i++) {
		v.push_back(counter);
		counter.x = int(triangle[0].x + (xStepSize * i));
		counter.y = int(triangle[0].y + (yStepSize * i));
	}

	xDiff = triangle[1].x - triangle[2].x;
	yDiff = triangle[1].y - triangle[2].y;
	numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	xStepSize = xDiff / numberOfSteps;
	yStepSize = yDiff / numberOfSteps;

	for (float i = 0.0; i < numberOfSteps; i++) {
		v.push_back(counter);
		counter.x = int(triangle[1].x + (xStepSize * i));
		counter.y = int(triangle[1].y + (yStepSize * i));
	}

	xDiff = triangle[2].x - triangle[0].x;
	yDiff = triangle[2].y - triangle[0].y;
	numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	xStepSize = xDiff / numberOfSteps;
	yStepSize = yDiff / numberOfSteps;

	for (float i = 0.0; i < numberOfSteps; i++) {
		v.push_back(counter);
		counter.x = int(triangle[2].x + (xStepSize * i));
		counter.y = int(triangle[2].y + (yStepSize * i));
	}

	return v;
}


// PROPER LINE DRAWER
void DrawLine(DrawingWindow& window, CanvasPoint from, CanvasPoint to, Colour lineColour) {

	std::vector<CanvasPoint> result;
	result = interpolatePointToPoint(from, to);

	float red = float(lineColour.red);
	float green = float(lineColour.green);
	float blue = float(lineColour.blue);

	//for (size_t i = 0; i < result.size(); i++) {
	//	std::cout << result[i] << "\n";
	//}

	float yMin = std::min(abs(from.y), abs(to.y));
	float yMax = std::max(abs(from.y), abs(to.y));
	float xMin = std::min(abs(from.x), abs(to.x));
	float xMax = std::max(abs(from.x), abs(to.x));

	for (size_t y = yMin; y < yMax+1; y++) { // move along y axis

		for (size_t x = xMin; x < xMax+1; x++) { // move along x axis

			for (size_t i = 0; i < result.size(); i++) { // scan vector of interpolated values

				if ((x == result[i].x) && (y == result[i].y)) { // if current x and y are in results

					uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
					window.setPixelColour(x, y, colour);
				}
			}
		}
	}
}


// TRIANGLE DRAWER
void DrawTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour lineColour) {
	std::vector<CanvasPoint> v;
	v = interpolatePointToPointWithVector(v, triangle[0], triangle[1]);
	v = interpolatePointToPointWithVector(v, triangle[1], triangle[2]);
	v = interpolatePointToPointWithVector(v, triangle[2], triangle[0]);

//	for (size_t i = 0; i < v.size(); i++) {
//		std::cout << v[i] << "\n";
//	}


//	float yMin = std::min(abs(triangle[0].y), abs(triangle[1].y));
//	yMin = std::min(abs(yMin), abs(triangle[2].y));
//	float yMax = std::max(abs(triangle[0].y), abs(triangle[1].y));
//	yMax = std::max(abs(yMax), abs(triangle[2].y));
//	float xMin = std::min(abs(triangle[0].x), abs(triangle[1].x));
//	xMin = std::max(abs(xMin), abs(triangle[2].x));
//	float xMax = std::max(abs(triangle[0].x), abs(triangle[1].x));
//	yMax = std::max(abs(yMax), abs(triangle[2].x));

	float red = float(lineColour.red);
	float green = float(lineColour.green);
	float blue = float(lineColour.blue);

	for (size_t y = 0; y < HEIGHT + 1; y++) { // move along y axis

		for (size_t x = 0; x < WIDTH + 1; x++) { // move along x axis

			for (size_t i = 0; i < v.size(); i++) { // scan vector of interpolated values

				if ((x == v[i].x) && (y == v[i].y)) { // if current x and y are in results

					uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
					window.setPixelColour(x, y, colour);
				}
			}
		}
	}

}


void draw0(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void draw1(DrawingWindow &window) { //black white gradient across window
	window.clearPixels();
	std::vector<float> result;
	result = interpolateSingleFloats(0, 255, WIDTH);

	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = result[x];
			float blue = result[x];
			float green = result[x];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void draw2(DrawingWindow &window) { //green black gradient across window
	window.clearPixels();
	std::vector<float> result;
	result = interpolateSingleFloats(255, 0, WIDTH);

	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = 0.0;
			float blue = 0.0;
			float green = result[x];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void draw3(DrawingWindow &window) { //multicoloured gradient
	window.clearPixels();

    glm::vec3 topLeft(255, 0, 0);        // red
    glm::vec3 topRight(0, 0, 255);       // blue
    glm::vec3 bottomRight(0, 255, 0);    // green
    glm::vec3 bottomLeft(255, 255, 0);   // yellow

    // this is extremely ugly!!!!! i know!!! but something more efficient would take time to write
    // i might do it. i should do it. i WILL do it!!!
	std::vector<float> leftSliceREDCHAN = interpolateSingleFloats(topLeft[0], bottomLeft[0], HEIGHT);
	std::vector<float> leftSliceGREENCHAN = interpolateSingleFloats(topLeft[1], bottomLeft[1], HEIGHT);
	std::vector<float> leftSliceBLUECHAN = interpolateSingleFloats(topLeft[2], bottomLeft[2], HEIGHT);

    std::vector<float> rightSliceREDCHAN = interpolateSingleFloats(topRight[0], bottomRight[0], HEIGHT);
	std::vector<float> rightSliceGREENCHAN = interpolateSingleFloats(topRight[1], bottomRight[1], HEIGHT);
	std::vector<float> rightSliceBLUECHAN = interpolateSingleFloats(topRight[2], bottomRight[2], HEIGHT);

    std::vector<float> topSliceREDCHAN = interpolateSingleFloats(topLeft[0], topRight[0], WIDTH);
	std::vector<float> topSliceGREENCHAN = interpolateSingleFloats(topLeft[1], topRight[1], WIDTH);
	std::vector<float> topSliceBLUECHAN = interpolateSingleFloats(topLeft[2], topRight[2], WIDTH);

    // height wise writer
    for (size_t x = 0; x < window.width; x++) {
		for (size_t y = 0; y < window.height; y++) {
			float red = leftSliceREDCHAN[y];
			float blue = leftSliceBLUECHAN[y];
			float green = leftSliceGREENCHAN[y];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}

    // width wise writer
    /*for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = topSliceREDCHAN[x];
			float blue = topSliceBLUECHAN[x];
			float green = topSliceGREENCHAN[x];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}*/
}

void draw4(DrawingWindow &window) { //lovely new function test
	window.clearPixels();

	glm::vec3 red(255, 0, 0);        // red
	glm::vec3 blue(0, 0, 255);       // blue
	glm::vec3 yellow(255, 255, 0);   // yellow
	glm::vec3 green(0, 255, 0);      // green

	std::vector<glm::vec3> leftSlice = interpolateVec3(red, blue, HEIGHT);
	std::vector<glm::vec3> rightSlice = interpolateVec3(yellow, green, HEIGHT);

    //for(size_t i = 0; i < result.size(); i++) {
    //    printf("%f, %f, %f\n", result[i].x, result[i].y, result[i].z);
    //} //prints vector created by interpolateVec3

// width wise writer
	for (size_t y = 0; y < window.height; y++) {
		std::vector<glm::vec3> horiSlice = interpolateVec3(leftSlice[y], rightSlice[y], WIDTH);
		// for each loop, a new version of this vector is generated
		// which takes the current y values of leftSlice and rightSlice and interpolates them.

		for (size_t x = 0; x < window.width; x++) {
			float red = horiSlice[x].x;
			float green = horiSlice[x].y;
			float blue = horiSlice[x].z;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void LineDraw(DrawingWindow &window) { //lets draw a lovely line
	window.clearPixels();

	CanvasPoint point1;
	CanvasPoint point2;

	point1 = CanvasPoint(0.0, 0.0);
	point2 = CanvasPoint(WIDTH/2, HEIGHT/2);
	//point1 = CanvasPoint(WIDTH / 2, 0.0);
	//point2 = CanvasPoint(WIDTH / 2, HEIGHT);
	//point1 = CanvasPoint(WIDTH/3, HEIGHT/2);
	//point2 = CanvasPoint(WIDTH / 3*2, HEIGHT/2);

	std::vector<CanvasPoint> result;
	result = interpolatePointToPoint(point1, point2);


    for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
		    float red = 0.0;
			float blue = 0.0;
			float green = 255.0;

			for (size_t i = 0; i < result.size(); i++) {
				//std::cout << result[i] << "\n";
				if ((x == result[i].x) && (y == result[i].y)) {
					red = 255.0;
					green = 0.0;
				}
			}

			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}

//    std::cout << lineTest << "\n";

	//for (size_t i = 0; i < result.size(); i++) {
	//	std::cout << result[i] << "\n";
	//}
//        printf("%f, %f, %f\n", result[i].x, result[i].y, result[i].z);
    //prints vector created by interpolateVec3

//    std::cout << specialMan.x << "\n" << specialMan.y << "\n" << typeid(specialMan.x).name() << "\n";
}

void LineDraw2(DrawingWindow &window) { //lets draw a lovely line USING A SEPARATE FUNCTION :O hi drslock if ur reading this uwu
	window.clearPixels();

	CanvasPoint point1, point2, point3, point4, point5, point6;
	Colour white, red, pink;

	point1 = CanvasPoint(0.0, 0.0);
	point2 = CanvasPoint(WIDTH / 2, HEIGHT / 2);
	white = Colour(255, 255, 255);
	DrawLine(window, point1, point2, white);

	point3 = CanvasPoint(WIDTH / 2, 0.0);
	point4 = CanvasPoint(WIDTH / 2, HEIGHT);
	red = Colour(255, 0, 0);
	DrawLine(window, point3, point4, red);

	point5 = CanvasPoint(WIDTH / 3, HEIGHT / 2);
	point6 = CanvasPoint(WIDTH / 3 * 2, HEIGHT / 2);
	pink = Colour(255, 0, 255);
	DrawLine(window, point5, point6, pink);	
}


void CreateTriangle(DrawingWindow& window) {
	CanvasPoint point1, point2, point3;
	CanvasTriangle triangle1;
	Colour colour1;

	point1 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
	point2 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
	point3 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
	triangle1 = CanvasTriangle(point1, point2, point3);
	colour1 = Colour(rand() % 255, rand() % 255, rand() % 255);

	DrawTriangle(window, triangle1, colour1);

	//DrawLine(window, point1, point2, colour1);
	//DrawLine(window, point2, point3, colour1);
	//DrawLine(window, point3, point1, colour1);
}



void draw5(DrawingWindow& window) {
	window.clearPixels(); 
	CreateTriangle(window);
}




void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;

        else if (event.key.keysym.sym == SDLK_0) draw0(window);
		else if (event.key.keysym.sym == SDLK_1) draw1(window);
		else if (event.key.keysym.sym == SDLK_2) draw2(window);
		else if (event.key.keysym.sym == SDLK_3) draw3(window);
		else if (event.key.keysym.sym == SDLK_4) draw4(window);
		else if (event.key.keysym.sym == SDLK_5) LineDraw2(window);
		else if (event.key.keysym.sym == SDLK_6) draw5(window);

		else if (event.key.keysym.sym == SDLK_u) CreateTriangle(window);

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
