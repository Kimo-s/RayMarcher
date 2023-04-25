#include <stdlib.h>
#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
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


float cylinder(float x, float y, float z) {
	float R = 0.6f;

	Vector normal(0.0, 0.0, 1.0);

	float height = 1.0f;

	Vector pos(x, y, z);

	if (fabs(z) > height) {
		return -1.0f;
	}
	else {
		return R - (pos - (pos * normal) * normal).magnitude();
	}
}

Vector tornadoFunc(float x, float y, float z){
	Vector pos(-y,x,0);
	pos.normalize();
	return pos;
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
	
	#pragma omp parallel for schedule(dynamic) shared(percent, f, DSM) 
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			float s = snear;
			Color col(0.0f, 0.0f, 0.0f, 0.0f);
			Vector ray = c.view(i * dx, j * dy) * s + cameraPosition;
			float T = 1.0f;

			
			while (s < sfar) {
				float res = f->eval(ray);
				if (res > 0.0f) {
					float dT = exp(-k * ds * res);
					Color dtls(0.0f, 0.0f, 0.0f, 0.0f);
					for (int q = 0; q < DSM.size(); q++) {
						dtls += lightColor[q] * exp(-k * DSM[q]->evalFrust(ray));
						//cout << DSM[q]->evalFrust(ray) << endl;
					}
					Color fieldCol = colField->eval(ray);
					col += fieldCol * T * ((1.0f - dT)/k) * dtls;
					T *= dT;
				}
				ray += c.view(i * dx, j * dy) * ds;
				s += ds;
			}
			/*if(T != 1.0f)
				printf("T = %f\n", T);*/
			buf[4 * i + 4 * j * width] = col.r;
			buf[4 * i + 4 * j * width + 1] = col.g;
			buf[4 * i + 4 * j * width + 2] = col.b;
			buf[4 * i + 4 * j * width + 3] = 1.0f - T;
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


int main() {
	Color color1(1.0f, 1.0f, 0.0f, 1.0f);
	Color color2(1.0f, 1.0f, 1.0f, 1.0f);
	Color color3(1.0f, 1.0f, 5.0f, 1.0f);


	// Input parameters
	int width = 1920;
	int height = 1080;
	float ds = 0.0005f;
	float sfar = 6.0f;
	float snear = 2.0f;
	float k = 1.0f;

	Camera c;
	float r = 5.0f;
	float theta = M_PI / 2.0f;
	c.setAspectRatio(1.0 * width / height);


	float dx = 0.12f;
	VolumeParms volparms = {
		50,
		50,
		50,
		dx,
		dx,
		dx,
		Vector(-3.0,-3.0,-3.0)
	};

	scalarFieldT inith = funcField(sphere).scale(0.4).translate(Vector(0,0,0.4));
	VectorField V = constantVectorField(Vector(0.0,0.0,-0.3));

	VectorField xmap = identityVectorField();

	int N = 100;
	float dt = 0.04;
	int i = 0;

	{
		VolumeGrid<float> dsmKey(inith, Vector(0.0f, 0.0f, -5.0f), 300, 300, 300, 2.0f, 8.0f, 45.0*(M_PI/180.0));

		std::vector< VolumeGrid<float>* > dsmMap = {&dsmKey};
		std::vector< Color > lightColorMap = {Color(1.0,1.0,1.0,1.0)};


		ColorField colField = colorMaskField(inith, color1);

		char name[100];

		Vector cameraCenter = Vector(0, r, -1.0f);
		c.setEyeViewUp(cameraCenter, -1.0f * cameraCenter, Vector(0, 0, 1));
		sprintf(name, "orignial.exr");
		rayMarch(inith, dsmMap, lightColorMap, colField, c, cameraCenter, ds, snear, sfar, k, width, height, name);
	}

	for (int frame = 0; frame < 120; frame++) {

		// with characteristic map 
		xmap = advect(xmap, V, dt);
		xmap = xmap - identityVectorField();
		xmap = gridField(xmap, 50, 50, 50, dx, dx, dx, Vector(-3.0, -3.0, -3.0));
		xmap = xmap + identityVectorField();
		scalarFieldT h = warp(inith, xmap);
		
		// h = advect(h, V, dt);
		h = gridField(h, N, N, N, 6.0/N, 6.0/N, 6.0/N, Vector(-3.0, -3.0, -3.0));

		V = advect(V, V, dt) + (funcField(tornadoFunc)*h)*dt;
		// incompress(V, &volparms);

		VolumeGrid<float> dsmKey(h, Vector(0.0f, 0.0f, -5.0f), 400, 400, 400, 2.0f, 8.0f, 45.0*(M_PI/180.0));

		std::vector< VolumeGrid<float>* > dsmMap = {&dsmKey};
		std::vector< Color > lightColorMap = {Color(1.0,1.0,1.0,1.0)};


		ColorField colField = colorMaskField(h, color1);

		char name[100];

		Vector cameraCenter = Vector(0, r, -1.0f);
		c.setEyeViewUp(cameraCenter, -1.0f * cameraCenter, Vector(0, 0, 1));
		sprintf(name, "out%03d.exr", i);
		rayMarch(h, dsmMap, lightColorMap, colField, c, cameraCenter, ds, snear, sfar, k, width, height, name);
		i+=1;
	}

}