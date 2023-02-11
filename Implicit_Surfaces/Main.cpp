#include <stdlib.h>
#include <string>
#include <iostream>
#include <math.h>
#include "IFSclasses.h"
#include "Camera.h"
#include "Vector.h"
#include "VolumeClasses.h"
#include "FieldClasses.h"

#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"
#define M_PI 3.14159265358979323846
using namespace ifs;

#define sq(A) (A*A)

#include <omp.h>



float sphere(float x, float y, float z) {
	return 1.0f - sqrt(x * x + y * y + z * z) ;
}

float steinerPatch(float x, float y, float z) {
	return -(sq(x)*sq(y)+sq(x)*sq(z)+sq(y)*sq(z)-x*y*z);
}

float ellipse(float x, float y, float z0) {
	Vector P(x, y, z0);
	Vector n(0, 0.0f, 1.0f);

	float z = P * n;
	Vector xbar = P - z * n;

	float rmajor = 0.3f;
	float rminor = 0.8f;

	return 1.0f-(sq(z)/sq(rmajor))-(xbar*xbar/sq(rminor));
}

float ellipse2(float x, float y, float z0) {
	Vector P(x, y, z0);
	Vector n(0, 0.0f, 1.0f);

	float z = P * n;
	Vector xbar = P - z * n;

	float rmajor = 0.8f;
	float rminor = 0.3f;

	return 1.0f - (sq(z) / sq(rmajor)) - (xbar * xbar / sq(rminor));
}

float box(float x, float y, float z0) {

	float q = 2;
	float R = 1.0f;

	return pow(R,2*q)- pow(x, 2 * q)- pow(y, 2 * q)- pow(z0, 2 * q);
}

float icosahedron(float x, float y, float z) {
	Vector P(x, y, z);

	float T = 2.0f;
	float res;
	if (P.magnitude() <= 1.8f * M_PI) {
		res = cos(x + T * y) + cos(x - T * y) + cos(y + T * z) + cos(y - T * z) + cos(z - T * x) + cos(z + T * x) - 2.0f;
	}
	else {
		res = -1.8 * M_PI;
	}

	return res;
}

float cone(float x, float y, float z0) {
	Vector P(x, y, z0);
	
	Vector n(0, 0.0f, 1.0f);
	float h = 1.0f;
	float theta = 0.4f;

	float res;

	if (P * n < 0) {
		res = P * n;
	}
	else if (P * n > h) {
		res = h - P * n;
	}
	else {
		res = P * n - P.magnitude() * cos(theta);
	}

	return res;
}

float tours(float x, float y, float z) {
	Vector P(x, y, z);
	Vector n(0.0f, 1.0f, 0.0f);
	Vector xbar = P - (P * n) * n;

	float Rmajor = 1.0f;
	float Rminor = 0.3f;


	return 4* Rmajor* Rmajor*(xbar*xbar)-((P*P+ Rmajor* Rmajor- Rminor* Rminor)* (P * P + Rmajor * Rmajor - Rminor * Rminor));
}

void rayMarch(ifs::ScalarField f, VolumeGrid<float>* DSM, Color lightColor, FieldBase<Color>* colField, Camera c, Vector cameraPosition, float ds, float snear, float sfar, float k, int width, int height, string filename) {

	float dx = (1.0f / (width - 1));
	float dy = (1.0f / (height - 1));
	std::vector<float> buf(width * height * 4);

	int percent = 0;
	
	
	#pragma omp parallel for schedule(dynamic) num_threads(20) collapse(2) shared(percent) 
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			float s = snear;
			//Vector col(0.0f, 0.0f, 0.0f);
			Color col(0.0f, 0.0f, 0.0f, 0.0f);
			Vector ray = c.view(i * dx, j * dy) * s + cameraPosition;
			float T = 1.0f;
			while (s < sfar) {
				float res = f.eval(ray.X(), ray.Y(), ray.Z());
				if (res > 0.0f) {
					//printf("(%f, %f, %f) value: %f\n", ray[0], ray[1], ray[2], res);
					float dT = exp(-k * ds * res);
					float dTL = exp(-k * DSM->eval(ray));
					//Vector colPoint = f.evalColor(ray.X(), ray.Y(), ray.Z());
					Color fieldCol = colField->eval(ray);
					//Vector colPoint = Vector(fieldCol.r, fieldCol.g, fieldCol.b);
					//Vector colPointL = Vector(lightColor.r, lightColor.g, lightColor.b);
					col += fieldCol * T * ((1.0f - dT) / k); // *lightColor * dTL;
					T *= dT;
				}
				ray += c.view(i * dx, j * dy) * ds;
				s += ds;
			}
			buf[4 * i + 4 * j * width] = col.r;
			buf[4 * i + 4 * j * width + 1] = col.g;
			buf[4 * i + 4 * j * width + 2] = col.b;
			buf[4 * i + 4 * j * width + 3] = 1 - T;
		}
		//printf("Scanline %d/%d finished.\n", j, height);
		percent += 1;
		printf("\rCreating frame (%s): %.2f%%", filename.c_str(), (percent * 1.0f / height)*100.0f );
	}

	printf("\n");

	const char* err = NULL;
	int ret = SaveEXR(buf.data(), width, height, 4, 0, filename.c_str(), &err);
	if (ret != TINYEXR_SUCCESS) {
		if (err) {
			printf("Error saving image file: %s\n", err);
			FreeEXRErrorMessage(err);
		}
	}
}

ifs::ScalarField sym(ifs::ScalarField f, ScalarField g, Point P) {
	f = f.max(g.translate(P));
	P.x = -P.x;
	f = f.max(g.translate(P));
	return f;
}

int main() {
	Vector color1(0.3f,0.7f,0.2f);
	Vector color2(0.0f, 0.8f, 0.5f);
	Vector color3(0.6f, 0.0f, 0.9f);


	//// Main Body
	//ifs::ScalarField h(box, color1);
	//h = h.mask();
	//h = h.max(ifs::ScalarField(cone, color3).translate(ifs::Point{ 0.0f,0.0f,3.0f }).mask());
	//h = h.max(ScalarField(sphere, color1).translate(Point{ 0.0f, 0.0f, 1.5f }).mask());
	//

	//// Arms
	//float armscale = 0.7f;
	//h = sym(h, ScalarField(sphere, color2).scale(Point{0.0f,0.0f,0.0f},armscale).mask(), Point{ armscale * 1.0f, 0.0f, 0.7f });
	//h = sym(h, ScalarField(ellipse, color3).scale(Point{ 0.0f,0.0f,0.0f }, armscale).mask(), Point{ armscale * 2.3f, 0.0f, 0.7f });
	//h = sym(h, ScalarField(sphere, color2).scale(Point{ 0.0f,0.0f,0.0f }, armscale).mask(), Point{ armscale * 3.6f, 0.0f, 0.7f });
	//h = sym(h, ScalarField(ellipse, color3).scale(Point{ 0.0f,0.0f,0.0f }, armscale).mask(), Point{ armscale * 4.6f, 0.0f, 0.7f });
	//h = sym(h, ScalarField(icosahedron, color1).scale(Point{ 0.0f,0.0f,0.0f }, armscale*0.1f).mask(), Point{ armscale * 5.6f, 0.0f, 0.7f });

	//// Legs
	//float T = -0.5f;
	//float legscale = 0.5f;
	//h = sym(h, ScalarField(steinerPatch, color2).scale(Point{ 0.0f,0.0f,0.0f }, legscale).mask(), Point{  0.8f, 0.0f, legscale * -1.0f + T });
	//h = sym(h, ScalarField(tours, color3).scale(Point{ 0.0f,0.0f,0.0f }, legscale*0.5f).mask(), Point{ 0.8f, 0.0f, legscale * -1.5f + T });
	//h = sym(h, ScalarField(ellipse2, color1).scale(Point{ 0.0f,0.0f,0.0f }, legscale).mask(), Point{ 0.8f, 0.0f, legscale * -2.5f + T });
	//h = sym(h, ScalarField(sphere, color2).scale(Point{ 0.0f,0.0f,0.0f }, legscale).mask(), Point{ 0.8f, 0.0f, legscale * -3.5f + T });
	//h = sym(h, ScalarField(ellipse2, color1).scale(Point{ 0.0f,0.0f,0.0f }, legscale).mask(), Point{ 0.8f, 0.0f, legscale * -4.5f + T });
	//h = sym(h, ScalarField(cone, color3).scale(Point{ 0.0f,0.0f,0.0f }, legscale).mask(), Point{ 0.8f, 0.0f, legscale * -5.3f + T });
	//-------------------------------------------------------------------->

	ScalarField f = ScalarField(sphere);
	ScalarField g = ScalarField(sphere).scale(Point{ 0.0f,0.0f,0.0f }, 0.5f).translate(Point{ 1.0f,0.0f,0.0f });
	VolumeGrid<float> t = VolumeGrid<float>(f.max(g), 500, 500, 500, 0.02f, 0.02f, 0.02f, Vector(-5.0f, -5.0f, -5.0f));

	VolumeGrid<float> dsm;
	//VolumeGrid<float> L = createGrid("bunny.obj", 300, 300, 300, 0.033f, 0.033f, 0.033f, Vector(-4.0f, -4.0f, -4.0f));


	ScalarField h(&t, color3);
	ColorFieldMask colField = ColorFieldMask(h, Color(0.2f, 0.1f, 0.2f, 1.0f));


	// Input parameters
	int width = 1920;
	int height = 1080;
	float ds = 0.01f;
	float sfar = 8.0f;
	float snear = 1.0f;
	float k = 0.4f;

	Camera c;
	float r = 5.0f;
	float theta = M_PI / 2.0f;
	int N = 120;
	float dt = 2.0f * M_PI / N;


	for (int i = 0; i < 120; i++) {
		Vector cameraCenter = Vector(cos(theta + i * dt) * r, sin(theta + i * dt) * r, 0);
		c.setEyeViewUp(cameraCenter, -1.0f * cameraCenter, Vector(0, 0, 1));


		char name[100];
		sprintf_s(name, 100, "../../exrFiles/out%03d.exr", i + 1);
		rayMarch(h, &dsm, Color(1.0f, 0.5f, 0.0f, 1.0f), &colField, c, cameraCenter, ds, snear, sfar, k, width, height, name);
	}
	

}