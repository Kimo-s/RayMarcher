#include "VolumeClasses.h"
#include "IFSclasses.h"
#include "Camera.h"
#include "Color.h"
#include "OBJ_Loader.h"
#include "FieldClasses.h"
#include <string>
#include <memory>
#include <omp.h>
#include <limits>

template<>
ifs::VolumeGrid<float>::VolumeGrid<float>(int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos) {
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;
	this->deltax = deltax;
	this->deltay = deltay;
	this->deltaz = deltaz;
	this->startPos = startPos;
}

// Volume Grid
ifs::VolumeGrid<float>::VolumeGrid<float>(scalarFieldT vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos) {
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;
	this->deltax = deltax;
	this->deltay = deltay;
	this->deltaz = deltaz;
	this->startPos = startPos;

	data = new float*[(Nx * Ny * Nz) / 64];
	//data = new float* [(Nx * Ny * Nz) / 64];
	for (int p = 0; p < (Nx * Ny * Nz) / 64; p++) {
		data[p] = NULL;
	}

	int finished = 0;

	#pragma omp parallel for schedule(dynamic) num_threads(20) shared(finished) 
	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
			for (int i = 0; i < Nx; i++) {
				Vector pos = Vector(i * deltax + startPos[0], j * deltay + startPos[1], k * deltaz + startPos[2]);
				set(i, j, k, vol->eval(pos));
			}
		}
		finished += 1;
		printf("\rCreating scalar grid from scalarfield: %0.2f%%", (finished + 1) * 1.0f / Nz);
	}
	printf("\n");
}

//template<>
//ifs::VolumeGrid<Color>::VolumeGrid<Color>(ColorField vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
//{
//	this->Nx = Nx;
//	this->Ny = Ny;
//	this->Nz = Nz;
//	this->deltax = deltax;
//	this->deltay = deltay;
//	this->deltaz = deltaz;
//	this->startPos = startPos;
//
//	data = make_unique<Color* []>((Nx * Ny * Nz) / 64);
//	//data = new float* [(Nx * Ny * Nz) / 64];
//	for (int p = 0; p < (Nx * Ny * Nz) / 64; p++) {
//		data[p] = NULL;
//	}
//
//	int finished = 0;
//
//	#pragma omp parallel for schedule(dynamic) num_threads(20) shared(finished) 
//	for (int k = 0; k < Nz; k++) {
//		for (int j = 0; j < Ny; j++) {
//			for (int i = 0; i < Nx; i++) {
//				Vector pos = Vector(i * deltax + startPos[0], j * deltay + startPos[1], k * deltaz + startPos[2]);
//				set(i, j, k, vol->eval(pos));
//			}
//		}
//		finished += 1;
//		printf("\rCreating scalar grid from scalarfield: %0.2f%%", (finished + 1) * 1.0f / Nz);
//	}
//	printf("\n");
//}


float rayMarchLight(scalarFieldT f, Vector point, Vector lightPos, float ds) {

	Vector direction = (point - lightPos).unitvector();
	float distance = (point - lightPos).magnitude();

	float s = 0.0f;
	Vector ray = lightPos + s * direction;
	float T = 0.0f;
	while ((point-ray).magnitude() > ds) {
		float res = f->eval(ray);
		if (res > 0.0f) {
			T += res;
		}
		ray += direction * ds;
		s += ds;
	}

	return T;
}


// Deep shadow maps
VolumeGrid<float>::VolumeGrid<float>(scalarFieldT vol, Vector lightPos, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos) {
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;
	this->deltax = deltax;
	this->deltay = deltay;
	this->deltaz = deltaz;
	this->startPos = startPos;

	data = new float* [(Nx * Ny * Nz) / 64];
	//data = new float* [(Nx * Ny * Nz) / 64];
	for (int p = 0; p < (Nx * Ny * Nz) / 64; p++) {
		data[p] = NULL;
	}

	int finished = 0;

	#pragma omp parallel for schedule(dynamic) num_threads(20) shared(finished) 
	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
			for (int i = 0; i < Nx; i++) {
				Vector pos = Vector(i * deltax + startPos[0], j * deltay + startPos[1], k * deltaz + startPos[2]);
				set(i, j, k, rayMarchLight(vol, pos, lightPos, 0.01f));
			}
		}
		finished += 1;
		printf("\rCreating light DSM: %0.2f%%", (finished + 1) * 1.0f / Nz);
	}
	printf("\n");
}


float distanceToTriangle(Vector p, Vector a, Vector b, Vector c) {
	const Vector ab = b - a;
	const Vector ac = c - a;
	const Vector ap = p - a;

	const float d1 = ab * ap;
	const float d2 = ac * ap;
	if (d1 <= 0.f && d2 <= 0.f) return (a - p).magnitude(); //#1

	const Vector bp = p - b;
	const float d3 = ab * bp;
	const float d4 = ac * bp;
	if (d3 >= 0.f && d4 <= d3) return (p - b).magnitude(); //#2

	const Vector cp = p - c;
	const float d5 = ab * cp;
	const float d6 = ac * cp;
	if (d6 >= 0.f && d5 <= d6) return (c - p).magnitude(); //#3

	const float vc = d1 * d4 - d3 * d2;
	if (vc <= 0.f && d1 >= 0.f && d3 <= 0.f)
	{
		const float v = d1 / (d1 - d3);
		return ((a + v * ab) - p).magnitude(); //#4
	}

	const float vb = d5 * d2 - d1 * d6;
	if (vb <= 0.f && d2 >= 0.f && d6 <= 0.f)
	{
		const float v = d2 / (d2 - d6);
		return ((a + v * ac) - p).magnitude(); //#5
	}

	const float va = d3 * d6 - d5 * d4;
	if (va <= 0.f && (d4 - d3) >= 0.f && (d5 - d6) >= 0.f)
	{
		const float v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		return ((b + v * (c - b)) - p).magnitude(); //#6
	}

	const float denom = 1.f / (va + vb + vc);
	const float v = vb * denom;
	const float w = vc * denom;
	return ((a + v * ab + w * ac) - p).magnitude(); //#0
}



bool rayTriangleIntersect(
	const Vector& orig, const Vector& dir,
	const Vector& v0, const Vector& v1, const Vector& v2,
	double& t)
{
	float kEpsilon = 1.0e-5f;

	// compute the plane's normal
	Vector v0v1 = v1 - v0;
	Vector v0v2 = v2 - v0;
	// no need to normalize
	Vector N = v0v1^v0v2; // N
	float area2 = N.magnitude();

	// Step 1: finding P

	// check if the ray and plane are parallel.
	float NdotRayDirection = N*dir;
	if (fabs(NdotRayDirection) < kEpsilon) // almost 0
		return false; // they are parallel, so they don't intersect! 

	// compute d parameter using equation 2
	double d = -N*v0;

	// compute t (equation 3)
	t = -(N*orig + d) / NdotRayDirection;

	// check if the triangle is behind the ray
	if (t < 0) return false; // the triangle is behind

	// compute the intersection point using equation 1
	Vector P = orig + t * dir;

	// Step 2: inside-outside test
	Vector C; // vector perpendicular to triangle's plane

	// edge 0
	Vector edge0 = v1 - v0;
	Vector vp0 = P - v0;
	C = edge0^vp0;
	if (N*C < 0) return false; // P is on the right side

	// edge 1
	Vector edge1 = v2 - v1;
	Vector vp1 = P - v1;
	C = edge1^vp1;
	if (N*C < 0)  return false; // P is on the right side

	// edge 2
	Vector edge2 = v0 - v2;
	Vector vp2 = P - v2;
	C = edge2^vp2;
	if (N*C < 0) return false; // P is on the right side;

	return true; // this ray hits the triangle
}

//template<>
//Color VolumeGrid<Color>::eval<Color>(float x, float y, float z) {
//	Vector p(1.0f * (x - startPos[0]) / deltax, 1.0f * (y - startPos[1]) / deltay, 1.0f * (z - startPos[2]) / deltaz);
//
//	Vector lowerLeft(static_cast<int>(floor(p.X())), static_cast<int>(floor(p.Y())), static_cast<int>(floor(p.Z())));
//
//	float weight[3];
//	Color value = 0.0f;
//
//	for (int i = 0; i < 2; ++i) {
//		int cur_x = lowerLeft[0] + i;
//		weight[0] = 1.0f - abs(p[0] - cur_x);
//		for (int j = 0; j < 2; ++j)
//		{
//			int cur_y = lowerLeft[1] + j;
//			weight[1] = 1.0 - std::abs(p[1] - cur_y);
//			for (int k = 0; k <= 1; ++k)
//			{
//				int cur_z = lowerLeft[2] + k;
//				weight[2] = 1.0 - std::abs(p[2] - cur_z);
//				if (index2(cur_x, cur_y, cur_z) < Nx * Ny * Nz && index2(cur_x, cur_y, cur_z) > 0) {
//					value += weight[0] * weight[1] * weight[2] * get(cur_x, cur_y, cur_z);
//				}
//			}
//		}
//	}
//
//	return value;
//}


ifs::VolumeGrid<float>::VolumeGrid<float>(const char* filename, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos) {

	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;
	this->deltax = deltax;
	this->deltay = deltay;
	this->deltaz = deltaz;
	this->startPos = startPos;
	this->data = new float* [(Nx * Ny * Nz) / 64];
	this->defaultValue = -1000.0f;
	for (int p = 0; p < (Nx * Ny * Nz) / 64; p++) {
		this->data[p] = NULL;
	}

	objl::Loader loader;
	bool data = loader.LoadFile(filename);

	if (!data) {
		std::cout << "Error loading the obj file. Please make sure the name of the file is correct." << endl;
		return;
	}

	int finished = 0;

	int bandwidth = 30;
	for (int ind = 0; ind < loader.LoadedIndices.size(); ind += 3) {

		Vector p1 = Vector(loader.LoadedVertices[loader.LoadedIndices[ind]].Position.X, loader.LoadedVertices[loader.LoadedIndices[ind]].Position.Y, loader.LoadedVertices[loader.LoadedIndices[ind]].Position.Z);
		Vector p2 = Vector(loader.LoadedVertices[loader.LoadedIndices[ind + 1]].Position.X, loader.LoadedVertices[loader.LoadedIndices[ind + 1]].Position.Y, loader.LoadedVertices[loader.LoadedIndices[ind + 1]].Position.Z);
		Vector p3 = Vector(loader.LoadedVertices[loader.LoadedIndices[ind + 2]].Position.X, loader.LoadedVertices[loader.LoadedIndices[ind + 2]].Position.Y, loader.LoadedVertices[loader.LoadedIndices[ind + 2]].Position.Z);

		//file << "T" << j / 3 << ": " << loader.LoadedVertices.Indices[j] << ", " << loader.LoadedVertices.Indices[j + 1] << ", " << loader.LoadedVertices.Indices[j + 2] << "\n";


		Vector LLC(std::min(std::min(p1.X(), p2.X()), p3.X())
			, std::min(std::min(p1.Y(), p2.Y()), p3.Y())
			, std::min(std::min(p1.Z(), p2.Z()), p3.Z()));
		Vector URC(std::max(std::max(p1.X(), p2.X()), p3.X())
			, std::max(std::max(p1.Y(), p2.Y()), p3.Y())
			, std::max(std::max(p1.Z(), p2.Z()), p3.Z()));

		int li, lj, lk;
		int ui, uj, uk;

		li = static_cast<int>(floor((LLC.X() - startPos[0]) / deltax)) - bandwidth;
		lj = static_cast<int>(floor((LLC.Y() - startPos[1]) / deltay)) - bandwidth;
		lk = static_cast<int>(floor((LLC.Z() - startPos[2]) / deltaz)) - bandwidth;

		ui = static_cast<int>(floor((URC.X() - startPos[0]) / deltax)) + bandwidth;
		uj = static_cast<int>(floor((URC.Y() - startPos[1]) / deltay)) + bandwidth;
		uk = static_cast<int>(floor((URC.Z() - startPos[2]) / deltaz)) + bandwidth;

		#pragma omp parallel for schedule(dynamic) num_threads(20)
		for (int k = lk; k < uk; k++) {
			for (int j = lj; j < uj; j++) {
				for (int i = li; i < ui; i++) {
					Vector pos = Vector(i * deltax + startPos[0], j * deltay + startPos[1], k * deltaz + startPos[2]);

					float d = distanceToTriangle(pos, p1, p2, p3);
					
					//printf("Getting (%d,%d,%d)\n", i, j, k);
					if (fabs(this->get(i, j, k)) > d) {
						this->set(i, j, k, d);
						//printf("(%d,%d,%d) = %f, d = %f\n", i, j, k, grid.get(i, j, k), d);
					}
				}
			}
		}
		finished += 1;
		printf("\rCreating sgd from (%s): %0.2f%%", filename, (finished + 1) * 1.0f / (loader.LoadedIndices.size()/3));
	}
	printf("\n");


	finished = 0;

	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
			for (int i = 0; i < Nx; i++) {
				this->set(i, j, k, -fabs(this->get(i, j, k)));
			}
		}
	}

	#pragma omp parallel for schedule(dynamic) num_threads(20) shared(finished)
	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {

			bool exitwhile = false;
			int numOfIntersactions = 0;
			int lasti = 0;

			Vector nextpos = Vector(deltax + startPos[0], j * deltay + startPos[1], k * deltaz + startPos[2]);
			Vector pos = Vector(startPos[0], j * deltay + startPos[1], k * deltaz + startPos[2]);
			Vector d = (nextpos - pos).unitvector();
			std::vector<int> ranges = {0};


			while(!exitwhile){
				bool intersacted = false;
				double t = numeric_limits<double>::max();
				double oldt;


				for (int ind2 = 0; ind2 < loader.LoadedIndices.size(); ind2 += 3) {
					Vector A = Vector(loader.LoadedVertices[loader.LoadedIndices[ind2]].Position.X, loader.LoadedVertices[loader.LoadedIndices[ind2]].Position.Y, loader.LoadedVertices[loader.LoadedIndices[ind2]].Position.Z);
					Vector B = Vector(loader.LoadedVertices[loader.LoadedIndices[ind2 + 1]].Position.X, loader.LoadedVertices[loader.LoadedIndices[ind2 + 1]].Position.Y, loader.LoadedVertices[loader.LoadedIndices[ind2 + 1]].Position.Z);
					Vector C = Vector(loader.LoadedVertices[loader.LoadedIndices[ind2 + 2]].Position.X, loader.LoadedVertices[loader.LoadedIndices[ind2 + 2]].Position.Y, loader.LoadedVertices[loader.LoadedIndices[ind2 + 2]].Position.Z);
	
					double tn;
					bool test;
					test = rayTriangleIntersect(pos, d, A, B, C, tn);
					if (test && tn > 0.0f && tn < t) {
						intersacted = true;
						oldt = t;
						t = tn;
					}
				}
				//std::cout << "value of intersacted " << intersacted  << " Value of t = " << t << endl;
				if (!intersacted) {
					exitwhile = true;
				}
				else if (intersacted && t != numeric_limits<float>::max()) {
					Vector pointOnTriang = pos + t * d;
					int hiti = static_cast<int>(ceil((pointOnTriang.X() - startPos[0]) / deltax));
					//Vector hit2detect = pos + oldt * d;
					//int hit2 = static_cast<int>(ceil((hit2detect.X() - startPos[0]) / deltax));
					//if (hit2 == hiti) {
					//	//printf("lol %d %d\n", hit2, lasti);
					//	set(hiti, j, k, fabs(get(hiti, j, k)));
					//	set(hit2, j, k, fabs(get(hit2, j, k)));
					//} 
					//ranges.push_back(hiti+1);
					if (numOfIntersactions % 2 == 1) {
						for (int q = lasti; q < hiti; q++) {
							float val = get(q, j, k);
							set(q, j, k, fabs(val));
						}
					}
					pos = pointOnTriang + 1.0e-6 * d;
					//pos.set((hiti+1)* deltax + startPos[0], pointOnTriang.Y(), pointOnTriang.Z());
					lasti = hiti;
					numOfIntersactions += 1;
				}
			}

			/*if (ranges.size() != 0) {
				for (int q = 0; q < ranges.size(); q++) {
					if (q % 2 == 1) {
						for (int r = ranges[q]; r < ranges[q + 1]; r++) {
							set(r, j, k, fabs(get(r, j, k)));
						}
					}
				}
			}*/
			
		}
		finished += 1;
		printf("\rChecking intersections (%s): %0.2f%%", filename, (finished + 1) * 1.0f / Nz);
	}
	printf("\n");
}

ifs::VolumeColorGrid::VolumeColorGrid(ColorField vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;
	this->deltax = deltax;
	this->deltay = deltay;
	this->deltaz = deltaz;
	this->startPos = startPos;

	data = make_unique<Color* []>((Nx * Ny * Nz) / 64);
	//data = new float* [(Nx * Ny * Nz) / 64];
	for (int p = 0; p < (Nx * Ny * Nz) / 64; p++) {
		data[p] = NULL;
	}

	int finished = 0;

	#pragma omp parallel for schedule(dynamic) num_threads(20) shared(finished) 
	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
			for (int i = 0; i < Nx; i++) {
				Vector pos = Vector(i * deltax + startPos[0], j * deltay + startPos[1], k * deltaz + startPos[2]);
				set(i, j, k, vol->eval(pos));
			}
		}
		finished += 1;
		printf("\rCreating Color grid field from scalarfield: %0.2f%%", (finished + 1) * 1.0f / Nz);
	}
	printf("\n");
}
