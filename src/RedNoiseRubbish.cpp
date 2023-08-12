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

using namespace std;


#define WIDTH 640
#define HEIGHT 480


//VECTOR SELECTION SORT

vector<CanvasPoint> selectionSort(vector<CanvasPoint> v) {
	int i, j, min_idx;
	CanvasPoint temp1, temp2;

	int n = v.size();
	
	for (size_t i = 0; i < n-1; i++) {
		min_idx = i;
		for (size_t j = i+1; j < n; j++) {
			if (v[j].y < v[min_idx].y) {
				min_idx = j;
			}
			if (min_idx != i) {
				temp1 = v[min_idx];
				temp2 = v[i];
				v[min_idx] = temp2;
				v[i] = temp1;
				// swap(&v[min_idx], &v[i])
			}
		}
	}

	for (size_t i = 0; i < v.size(); i++) {
		cout << v[i] << "\n";
	}
	cout << "\n";

	return v;
}


//POINT TO POINT INTERPOLATION
vector<CanvasPoint> interpolatePointToPoint(CanvasPoint from, CanvasPoint to) {
    vector<CanvasPoint> v;
    CanvasPoint counter = from;

    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numberOfSteps = max(abs(xDiff), abs(yDiff));
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
vector<CanvasPoint> interpolateTrianglePoints(CanvasTriangle triangle) {
    vector<CanvasPoint> v, v1, v2, v3;

	if ((triangle[0].y <= triangle[1].y) && (triangle[0].y <= triangle[2].y)) { //this code sucks and i hate it but whatever
		v.push_back(triangle[0]);
		if (triangle[1].y <= triangle[2].y) {
			v.push_back(triangle[1]);
			v.push_back(triangle[2]);
		}
		else {
			v.push_back(triangle[2]);
			v.push_back(triangle[1]);
		}
	}
	else if ((triangle[1].y <= triangle[0].y) && (triangle[1].y <= triangle[2].y)) {
		v.push_back(triangle[1]);
		if (triangle[0].y <= triangle[2].y) {
			v.push_back(triangle[0]);
			v.push_back(triangle[2]);
		}
		else {
			v.push_back(triangle[2]);
			v.push_back(triangle[0]);
		}
	}
	else if ((triangle[2].y <= triangle[0].y) && (triangle[2].y <= triangle[1].y)) {
		v.push_back(triangle[2]);
		if (triangle[0].y <= triangle[1].y) {
			v.push_back(triangle[0]);
			v.push_back(triangle[1]);
		}
		else {
			v.push_back(triangle[1]);
			v.push_back(triangle[0]);
		}
	}

	v1 = interpolatePointToPoint(v[0], v[1]);
	v2 = interpolatePointToPoint(v[1], v[2]);
	v3 = interpolatePointToPoint(v[2], v[0]);
	
	v.clear();
    v = v1;
    v.insert(v.end(), v2.begin(), v2.end());
    v.insert(v.end(), v3.begin(), v3.end());

	// v = selectionSort(v);

    return v;
}



// OUTLINED TRIANGLE DRAWER
void DrawOutlinedTriangle(DrawingWindow& window, CanvasTriangle triangle, Colour lineColour) {
	vector<CanvasPoint> v;
	v = interpolateTrianglePoints(triangle);

//	for (size_t i = 0; i < v.size(); i++) {
//		cout << v[i] << "\n";
//	}

	float yMin = min(abs(triangle[0].y), abs(triangle[1].y));
	yMin = min(abs(yMin), abs(triangle[2].y));
	float yMax = max(abs(triangle[0].y), abs(triangle[1].y));
	yMax = max(abs(yMax), abs(triangle[2].y));
	float xMin = min(abs(triangle[0].x), abs(triangle[1].x));
	xMin = min(abs(xMin), abs(triangle[2].x));
	float xMax = max(abs(triangle[0].x), abs(triangle[1].x));
	xMax = max(abs(xMax), abs(triangle[2].x));

	float red = float(lineColour.red);
	float green = float(lineColour.green);
	float blue = float(lineColour.blue);

	for (size_t y = yMin; y < yMax + 1; y++) { // move along y axis

		for (size_t x = xMin; x < xMax + 1; x++) { // move along x axis

			for (size_t i = 0; i < v.size(); i++) { // scan vector of interpolated values

				if ((x == v[i].x) && (y == v[i].y)) { // if current x and y are in results

					vector<CanvasPoint>::iterator it1 = v.begin() + i;
					v.erase(it1); // remove current point from vector

					uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
					window.setPixelColour(x, y, colour);
				}
			}
		}
	}
}


// FILLED TRIANGLE DRAWER
void DrawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour lineColour, Colour fillColour) {

	vector<CanvasPoint> v, vFill;
	v = interpolateTrianglePoints(triangle);

	float yMin = min(abs(triangle[0].y), abs(triangle[1].y));
	yMin = min(abs(yMin), abs(triangle[2].y));
	float yMax = max(abs(triangle[0].y), abs(triangle[1].y));
	yMax = max(abs(yMax), abs(triangle[2].y));
	float xMin = min(abs(triangle[0].x), abs(triangle[1].x));
	xMin = min(abs(xMin), abs(triangle[2].x));
	float xMax = max(abs(triangle[0].x), abs(triangle[1].x));
	xMax = max(abs(xMax), abs(triangle[2].x));

	float lineRed = float(lineColour.red);
	float lineGreen = float(lineColour.green);
	float lineBlue = float(lineColour.blue);

	float fillRed = float(fillColour.red);
    float fillGreen = float(fillColour.green);
    float fillBlue = float(fillColour.blue);


	for (size_t y = yMin; y < yMax + 1; y++) { // move along y axis

		vFill.clear();
		for (size_t i = 0; i < v.size(); i++) {
			if (v[i].y == y) {
				vFill.push_back(v[i]);
			}
		}

		if (vFill[0].x == vFill[vFill.size() - 1].x) {
			vFill.clear();
			//vFill.push_back(CanvasPoint(0, y));
			//vFill.push_back(CanvasPoint(WIDTH, y));
			//vFill = interpolatePointToPoint(vFill[0], vFill[vFill.size() - 1]);
		} else {
			vFill = interpolatePointToPoint(vFill[0], vFill[vFill.size() - 1]);
		}

		//for (size_t i = 0; i < vFill.size(); i++) {
		//	cout << vFill[i] << "\n";
		//}
		//cout << "\n";


		for (size_t x = xMin; x < xMax + 1; x++) { // move along x axis
			
			for (size_t i = 0; i < vFill.size(); i++) {

				if ((x == vFill[i].x) && (y == vFill[i].y)) { // if current x and y are in results

					uint32_t colour = (255 << 24) + (int(fillRed) << 16) + (int(fillGreen) << 8) + int(fillBlue);
					window.setPixelColour(x, y, colour);
				}
			}

			for (size_t i = 0; i < v.size(); i++) { // scan vector of interpolated values        

				if ((x == v[i].x) && (y == v[i].y)) { // if current x and y are in results

					vector<CanvasPoint>::iterator it1 = v.begin() + i;
					v.erase(it1); // remove current point from vector

					uint32_t colour = (255 << 24) + (int(lineRed) << 16) + (int(lineGreen) << 8) + int(lineBlue);
					window.setPixelColour(x, y, colour);
				}
			}
		}
	}
}

void CreateOutlinedTriangle(DrawingWindow &window) {
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

void CreateFilledTriangle(DrawingWindow &window) {
	CanvasPoint point1, point2, point3;
	CanvasTriangle triangle1;
	Colour outlineColour, fillColour;

	point1 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
	point2 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
	point3 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);

	// point1 = CanvasPoint(rand() % 100, rand() % 100);
	// point2 = CanvasPoint(rand() % 100, rand() % 100);
	// point3 = CanvasPoint(rand() % 100, rand() % 100);

	//(518, 476) (327, 215) (51, 126)

	//point1 = CanvasPoint(518, 476);
	//point2 = CanvasPoint(327, 215);
	//point3 = CanvasPoint(51, 126);

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
		if (event.key.keysym.sym == SDLK_LEFT) cout << "LEFT" << endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) cout << "RIGHT" << endl;
		else if (event.key.keysym.sym == SDLK_UP) cout << "UP" << endl;
		else if (event.key.keysym.sym == SDLK_DOWN) cout << "DOWN" << endl;

		else if (event.key.keysym.sym == SDLK_u) CreateOutlinedTriangle(window);
		else if (event.key.keysym.sym == SDLK_f) CreateFilledTriangle(window);

	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {

//    vector<glm::vec3> result;
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
