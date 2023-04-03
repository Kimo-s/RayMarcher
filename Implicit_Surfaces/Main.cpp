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

bool interesect(Vector pos, Vector d, Vector minCorner, Vector maxCorner) {
	double tmin = -INFINITY, tmax = INFINITY;

	if (d.X() != 0.0) {
		double tx1 = (minCorner.X() - pos.X()) / d.X();
		double tx2 = (maxCorner.X() - pos.X()) / d.X();

		tmin = max(tmin, min(tx1, tx2));
		tmax = min(tmax, max(tx1, tx2));
	}

	if (d.Y() != 0.0) {
		double ty1 = (minCorner.Y() - pos.Y()) / d.Y();
		double ty2 = (maxCorner.Y() - pos.Y()) / d.Y();

		tmin = max(tmin, min(ty1, ty2));
		tmax = min(tmax, max(ty1, ty2));
	}

	return tmax >= tmin;
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
	
	#pragma omp parallel for schedule(dynamic) num_threads(20) collapse(2) shared(percent, f, DSM) 
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			float s = snear;
			Color col(0.0f, 0.0f, 0.0f, 0.0f);
			Vector ray = c.view(i * dx, j * dy) * s + cameraPosition;
			float T = 1.0f;

			
			if (interesect(cameraPosition, c.view(i * dx, j * dy), Vector(-5.0f, -5.0f, -5.0f), Vector(5.0f, 5.0f, 5.0f))) {
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

void sym(scalarFieldT& f, scalarFieldT& g, Vector P, ColorField& curField, Color col) {
	scalarFieldT temp = g.translate(P);
	f = f.max(temp);
	curField = curField + colorMaskField(temp, col);

	P.set(-P.X(), P.Y(), P.Z());

	temp = g.translate(P).rotate(Vector(0.0f, 0.0f, 1.0f), -M_PI / 2.0f).translate(Vector(1.0f, 0.0f, 0.0f));
	f = f.max(temp);
	curField = curField + colorMaskField(temp, col);
}

void sym2(scalarFieldT& f, scalarFieldT& g, Vector P, ColorField& curField, Color col) {
	scalarFieldT temp = g.translate(P);
	f = f.max(temp);
	curField = curField + colorMaskField(temp, col);

	P.set(-P.X(), P.Y(), P.Z());

	temp = g.translate(P);
	f = f.max(temp);
	curField = curField + colorMaskField(temp, col);
}

int main() {
	Color color1(1.0f, 1.0f, 0.0f, 1.0f);
	Color color2(5.0f, 1.0f, 1.0f, 1.0f);
	Color color3(1.0f, 1.0f, 5.0f, 1.0f);

	//scalarFieldT h = funcField(sphere);

	/*h = h.addWispParticale(VolumeParms{ 500, 500, 500, 0.01f, 0.01f, 0.01f, Vector(0.0f, 0.0f, 0.0f) },
		FSPNParms{ 1.0f, 5, 2.5f, 3.2f, 1.7f, Vector(0.0f, 0.0f, 0.0f) }, 
		FSPNParms{ 1.0f, 4, 1.5f, 3.2f, 1.2f, Vector(1.0f, 0.0f, 0.0f) },
		Vector(0.0f, 0.0, 0.0),
		0.3f,
		1.0f,
		2.0f,
		1.0 / 2.0,
		10000000
	);*/
	//VolumeGrid<Vector> a(gradField(h), 500, 500, 500, 0.01f, 0.01f, 0.01f, Vector(0.0f, 0.0f, 0.0f));

	//h = gridField(h, 800, 800, 800, 0.0125f, 0.0125f, 0.0125f, Vector(-5.0f, -5.0f, -5.0f));
	
	// Input parameters
	int width = 1920;
	int height = 1080;
	float ds = 0.01f;
	float sfar = 10.0f;
	float snear = 3.0f;
	float k = 4.0f;

	Camera c;
	float r = 6.0f;
	float theta = M_PI / 2.0f;
	int N = 120;
	float dt = 2.0f * M_PI / N;
	c.setAspectRatio(1.0 * width / height);


	/*for (int i = 0; i < 120; i++) {
		Vector cameraCenter = Vector(cos(theta + i * dt) * r, sin(theta + i * dt) * r, -3.0f);
		c.setEyeViewUp(cameraCenter, -1.0f * cameraCenter, Vector(0, 0, 1));


		char name[100];
		sprintf_s(name, 100, "../../exrFiles/out%03d.exr", i + 1);
		rayMarch(h, dsmMap, lightColorMap, colField, c, cameraCenter, ds, snear, sfar, k, width, height, name);
	}*/


	//each parameter change should take 125 frames
	int frames = 125;

	Vector cameraCenter = Vector(0.0f, r, -3.0f);
	c.setEyeViewUp(cameraCenter, -1.0f * cameraCenter, Vector(0, 0, 1));

	FSPNParms params = {
		3.0f,
		3,
		1.9f,
		1.9f,
		0.6f,
		Vector(0.0, 0.0, 0.0),
		1.5f
	};

	scalarFieldT h = plane(Vector(0.0,0.0,1.0), Vector(0,0,0), params);
	//h = gridField(h, 500, 500, 500, 0.01f, 0.01f, 0.01f, Vector(-2.5f, -2.5f, -2.5f));

	VolumeGrid<float> dsmKey(h, Vector(2.5f, 0.0f, -4.0f), 500, 500, 500, 0.3f, 15.0f, M_PI/2.0);

	std::vector< VolumeGrid<float>* > dsmMap = {&dsmKey};
	std::vector< Color > lightColorMap = {Color(1.0,0.2,0.6,1.0)};


	ColorField colField = colorMaskField(h, color1);
	char name[100];
	sprintf_s(name, 100, "../../exrFiles/out%03d.exr", 0);
	rayMarch(h, dsmMap, lightColorMap, colField, c, cameraCenter, ds, snear, sfar, k, width, height, name);

}