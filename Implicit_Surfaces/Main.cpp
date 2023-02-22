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
	return 0.5f - sqrt(x * x + y * y + z * z) ;
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

void rayMarch(scalarFieldT f, std::vector<VolumeGrid<float>* > DSM, std::vector<Color> lightColor, ColorField colField, Camera c, Vector cameraPosition, float ds, float snear, float sfar, float k, int width, int height, string filename) {

	float dx = (1.0f / (width - 1));
	float dy = (1.0f / (height - 1));
	std::vector<float> buf(width * height * 4);

	if (DSM.size() != lightColor.size()) {
		std::cout << "Please supply the same number of colors as there are deep shadow maps.\n";
		return;
	}

	int percent = 0;
	
	#pragma omp parallel for schedule(dynamic) num_threads(30) collapse(2) shared(percent, f, DSM) 
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			float s = snear;
			//Vector col(0.0f, 0.0f, 0.0f);
			Color col(0.0f, 0.0f, 0.0f, 0.0f);
			Vector ray = c.view(i * dx, j * dy) * s + cameraPosition;
			float T = 1.0f;
			while (s < sfar) {
				//printf("The ray: %s\n", ray.__str__());
				float res = f->eval(ray);
				if (res > 0.0f) {
					float dT = exp(-k * ds * res);
					Color dtls= lightColor[0] * exp(-k * DSM[0]->eval(ray));
					//cout << dtls.print() << endl;
					/*for (int q = 0; q < DSM.size(); q++) {
						dtls += lightColor[q] * exp(-k * DSM[q]->eval(ray));
					}*/
					Color fieldCol = colField->eval(ray);
					col += fieldCol * 1.0e-6f + fieldCol * T * ((1.0f - dT)/k) * dtls;
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

void sym(scalarFieldT& f, scalarFieldT& g, Vector P, ColorField& curField, Color col) {
	scalarFieldT f2 = f.max((g).translate(P));
	curField = curField + colorMaskField(f2, col);
	P.set(-P.X(), P.Y(), P.Z());
	f = f2.max(g.translate(P).rotate(Vector(0.0f, 0.0f, 1.0f), -M_PI / 2.0f).translate(Vector(1.0f, 0.0f, 0.0f)));
	curField = curField + colorMaskField(f2.max(g.translate(P).rotate(Vector(0.0f, 0.0f, 1.0f), -M_PI / 2.0f)).translate(Vector(1.0f, 0.0f, 0.0f)), col);
}

int main() {
	Color color1(0.0f, 1.0f, 0.0f, 1.0f);
	Color color2(0.0f, 0.8f, 0.5f, 1.0f);
	Color color3(0.6f, 0.0f, 0.9f, 1.0f);


	//// Main Body
	//
	//scalarFieldT h = funcField(box).scale(0.6f);
	//ColorField colField = colorMaskField(h, color1);
	////h = h.mask();
	//scalarFieldT temp = funcField(cone).scale(0.6f).translate(Vector(0.0f, 0.0f, -1.5f));
	//h = h.max(temp);
	//colField = colField + colorMaskField(temp, color3);

	//temp = funcField(sphere).scale(0.5f).translate(Vector(0.0f, 0.0f, -2.0f));
	//h = h.max(temp);
	//colField = colField + colorMaskField(temp, color1);

	//temp = funcField(box).scale(0.6f).translate(Vector(0.0f, 0.0f, -1.5f));
	//h = h.max(temp);
	//colField = colField + colorMaskField(temp, color2);

	////// Arms
	//float armscale = 0.7f;
	//temp = funcField(sphere).scale(armscale);
	//sym(h, temp, Vector(-armscale * 1.0f, 0.0f, -0.7f), colField, color1);
	//temp = funcField(ellipse).scale(armscale);
	//sym(h, temp, Vector(-armscale * 2.3f, 0.0f, -0.7f), colField, color3);
	//temp = funcField(sphere).scale(armscale);
	//sym(h, temp, Vector(-armscale * 3.6f, 0.0f, -0.7f), colField, color2); //h = sym(h, ScalarField(sphere, color2).scale(Point{ 0.0f,0.0f,0.0f }, armscale).mask(), Point{ armscale * 3.6f, 0.0f, 0.7f });
	//temp = funcField(ellipse).scale(armscale);
	//sym(h, temp, Vector(-armscale * 4.6f, 0.0f, -0.7f), colField, color3); //h = sym(h, ScalarField(ellipse, color3).scale(Point{ 0.0f,0.0f,0.0f }, armscale).mask(), Point{ armscale * 4.6f, 0.0f, 0.7f });
	//temp = funcField(icosahedron).scale(armscale * 0.1f);
	//sym(h, temp, Vector(-armscale * 5.6f, 0.0f, -0.7f), colField, color1); //h = sym(h, ScalarField(icosahedron, color1).scale(Point{ 0.0f,0.0f,0.0f }, armscale*0.1f).mask(), Point{ armscale * 5.6f, 0.0f, 0.7f });

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

	/*scalarFieldT f = funcField(sphere);
	scalarFieldT g = funcField(sphere).translate(Vector(-1.0f, 0.f, 0.f));
	scalarFieldT s = f+g;

	ColorField f1 = colorMaskField(f, Color(0.2f, 0.8f, 0.2f, 1.0f));
	ColorField g1 = colorMaskField(g, Color(0.7f, 0.3f, 0.8f, 1.0f));
	ColorField colField = f1+g1; */


	scalarFieldT h = gridField("teapot.obj", 200, 200, 200, 0.05f, 0.05f, 0.05f, Vector(-5.0f, -5.0f, -5.0f)).rotate(Vector(1.0f, 0.0f, 0.0f), M_PI / 2.0f).translate(Vector(0.0f, 0.0f, 1.0f));
	ColorField colField = colorMaskField(h, color1);

	/*h = h.max(funcField(sphere).translate(Vector(-3.0, 0.0, 0.0)));
	colField = colField + colorMaskField(funcField(sphere).translate(Vector(-3.0, 0.0, 0.0)), color2);*/
	VolumeGrid<float> dsm(h, Vector(8.0, 8.0, 0.0), 300, 300, 300, 0.033f, 0.033f, 0.033f, Vector(-5.0f, -5.0f, -5.0f));


	std::vector< VolumeGrid<float>* > dsmMap = { &dsm };
	std::vector< Color > lightColorMap = { Color(1.0f, 1.0f, 1.0f, 1.0f) };

	// Input parameters
	int width = 1920;
	int height = 1080;
	float ds = 0.01f;
	float sfar = 8.0f;
	float snear = 1.0f;
	float k = 0.5f;

	Camera c;
	float r = 7.0f;
	float theta = M_PI / 2.0f;
	int N = 120;
	float dt = 2.0f * M_PI / N;


	for (int i = 0; i < 120; i++) {
		Vector cameraCenter = Vector(cos(theta + i * dt) * r, sin(theta + i * dt) * r, -2.0f);
		c.setEyeViewUp(cameraCenter, -1.0f * cameraCenter, Vector(0, 0, 1));


		char name[100];
		sprintf_s(name, 100, "../../exrFiles/out%03d.exr", i + 1);
		//f = (&f)->rotate(Vector(0.0f, 0.0f, 1.0f), (M_PI/120));
		//f1 = new ColorFieldMask(f, color3);
		rayMarch(h, dsmMap, lightColorMap, colField, c, cameraCenter, ds, snear, sfar, k, width, height, name);
	}
	

}