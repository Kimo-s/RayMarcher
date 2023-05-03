#include <stdlib.h>
#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include "IFSclasses.h"
#include "Camera.h"
#include "Vector.h"
#include "VolumeClasses.h"
#include "ScalarFieldFuncs.h"
#include "FieldClasses.h"

#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"
#define M_PI 3.14159265358979323846
using namespace ifs;

#define DB_PERLIN_IMPL
#include "db_perlin.hpp"

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

float flatPlane(float x, float y, float z){
	Vector pos(x,y,z);

	if (fabs(pos.X()) > 1.50) {
		return 0.0f;
	}
	else if (fabs(pos.Y()) > 1.50) {
		return 0.0f;
	}

	float val = pos * Vector(0.0,0.0,1.0);

	if(z > 0.5){
		val = -val;
	}

	return val;
}

Vector tornadoFunc(float x, float y, float z){
	Vector pos(-y,x,0);
	pos.normalize();
	//pos = pos * Vector(x,y,0).magnitude();
	Vector n = Vector(x,y,0) ;
	//float mag = n.magnitude();
	
	FSPNParms parms = {
		1.0,
		2,
		0.6f,
		5.0f,
		0.7f,
		Vector(0,0,0),
		1.5f
	};

	//float sign = -n.magnitude() + z;
	//if(sign > 0.0){
	//	sign = 30.0f;
	//}
	//n.normalize();

	float noise = evalFSPN(parms, Vector(x,y,z));
	//n = n * -fabs(evalFSPN(parms, Vector(x,y,z)));
	pos = pos * 2.0 - n * 9.0 * fabs(noise);
	pos[2] = noise;
	
	if(x == 0 && y == 0){
		return Vector(0,0,0.1);
	}

	if(z > 0.0){
		pos[2] += -4;
	} else if (z < 0.0) {
		pos[2] += 4;
	}

	return pos * 9.0f;
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
	
	#pragma omp parallel for schedule(static) shared(percent, f, DSM) 
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
	Color color1(0.75f, 0.75f, 0.75f, 1.0f);
	Color color2(1.0f, 0.0f, 1.0f, 1.0f);
	Color color3(1.0f, 1.0f, 5.0f, 1.0f);


	// Input parameters
	int width = 1920;
	int height = 1080;
	float ds = 0.002f;
	float sfar = 9.0f;
	float snear = 3.0f;
	float k = 0.06f;

	Camera c;
	float r = 6.0f;
	c.setAspectRatio(1.0 * width / height);

	int N = 400;
	float dx = 6.0/N;
	int Nvec = 50;
	float dxvec = 6.0/Nvec;
	VolumeParms volparms = {
		Nvec,
		Nvec,
		Nvec,
		dxvec,
		dxvec,
		dxvec,
		Vector(-3.0,-3.0,-3.0)
	};

	scalarFieldT h = funcField(sphere).scale(0.4).translate(Vector(0.3,0,0.6)).mask() * constantField(0.15);
	//scalarFieldT obs = funcField(sphere).scale(0.5).translate(Vector(0,0,-0.9));
	// scalarFieldT obs = funcField(flatPlane).translate(Vector(0,0,-1.7));
	scalarFieldT obs = constantField(0.0f);
	scalarFieldT source = funcField(sphere).scale(0.1).translate(Vector(0,0,1.2)).mask();

	VectorField V = constantVectorField(Vector(0.0,0.0,-0.1));


	VectorField xmap = identityVectorField();

	//float dt = 0.01;
	float dt = 0.1;
	int i = 0;

	// {
	// 	VolumeGrid<float> dsmKey(inith, Vector(0.0f, 0.0f, -5.0f), 300, 300, 300, 2.0f, 8.0f, 45.0*(M_PI/180.0));

	// 	std::vector< VolumeGrid<float>* > dsmMap = {&dsmKey};
	// 	std::vector< Color > lightColorMap = {Color(1.0,1.0,1.0,1.0)};


	// 	ColorField colField = colorMaskField(inith, color1);

	// 	char name[100];

	// 	Vector cameraCenter = Vector(0, r, -1.0f);
	// 	c.setEyeViewUp(cameraCenter, -1.0f * cameraCenter, Vector(0, 0, 1));
	// 	sprintf(name, "orignial.exr");
	// 	rayMarch(inith, dsmMap, lightColorMap, colField, c, cameraCenter, ds, snear, sfar, k, width, height, name);
	// }

	int timesteps = 1;
	scalarFieldT h2 = constantField(0.0);
	for (int frame = 0; frame < 360; frame++) {

		// for(int step = 0; step < timesteps; step++){
		// 	xmap = advect(xmap, V, dt);
		// 	xmap = xmap - identityVectorField();
		// 	xmap = gridField(xmap, Nvec, Nvec, Nvec, dxvec, dxvec, dxvec, Vector(-3.0,-3.0,-3.0));
		// 	xmap = xmap + identityVectorField();
			
		// 	//h = gridField(h, N, N, N, dx, dx, dx, Vector(-5.0, -5.0, -5.0));
		// 	h2 = warp(h, xmap);

		// 	VectorField forcingField = funcField(tornadoFunc);

		// 	V = advect(V, V, dt) + (forcingField*h2)*dt;
		// 	V = incompress(V, obs, &volparms);
		// }



		for(int step = 0; step < timesteps; step++){
			h = advect(h, V, dt);
			h = gridField(h, N, N, N, dx, dx, dx, Vector(-3.0, -3.0, -3.0));

			VectorField forcingField = funcField(tornadoFunc);


			V = advect(V, V, dt) + (forcingField*h)*dt;
			V = incompress(V, obs, &volparms);
		} 


		scalarFieldT finalField = h;
		VolumeGrid<float> dsmKey(finalField, Vector(0.0f, 0.0f, -7.0f), 700, 700, 700, 0.0f, 7.0f+3.0f, 70.0*(M_PI/180.0));
		//VolumeGrid<float> dsmKey2(finalField, Vector(0.0f, 7.0f, -3.0f), 700, 700, 700, 0.0f, 7.0f+3.0f, 70.0*(M_PI/180.0));

		std::vector< VolumeGrid<float>* > dsmMap = {&dsmKey};
		std::vector< Color > lightColorMap = {Color(1.0,1.0,1.0,1.0)};


		ColorField colField = colorMaskField(h, color1);
		//colField = colField + colorMaskField(obs, color2);

		char name[100];

		Vector cameraCenter = Vector(0, r, -2.0f);
		c.setEyeViewUp(cameraCenter, -1.0f * cameraCenter, Vector(0, 0, 1));
		sprintf(name, "out%03d.exr", i);
		rayMarch(finalField, dsmMap, lightColorMap, colField, c, cameraCenter, ds, snear, sfar, k, width, height, name);
		i+=1;
	}

}
