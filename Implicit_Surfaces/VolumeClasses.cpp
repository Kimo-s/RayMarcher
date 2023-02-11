#include "VolumeClasses.h"
#include "IFSclasses.h"
#include <string>
#include "Camera.h"
#include "FieldBase.h"
#include "Color.h"
#include <omp.h>
#include "OBJ_Loader.h"

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
ifs::VolumeGrid<float>::VolumeGrid<float>(ScalarField vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos) {
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;
	this->deltax = deltax;
	this->deltay = deltay;
	this->deltaz = deltaz;
	this->startPos = startPos;

	data = new float* [(Nx * Ny * Nz) / 64];
	for (int p = 0; p < (Nx * Ny * Nz) / 64; p++) {
		data[p] = NULL;
	}

	int finished = 0;

	#pragma omp parallel for schedule(dynamic) num_threads(20) shared(finished) 
	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
			for (int i = 0; i < Nx; i++) {
				Vector pos = Vector(i * deltax + startPos[0], j * deltay + startPos[1], k * deltaz + startPos[2]);
				set(i, j, k, vol.eval(pos));
			}
		}
		finished += 1;
		printf("\rCreating scalar grid from scalarfield: %0.2f%%", (finished + 1) * 1.0f / Nz);
	}
	printf("\n");
}


float rayMarchLight(ifs::ScalarField f, Vector point, Vector lightPos, float ds) {

	Vector direction = (point - lightPos).unitvector();
	float distance = (point - lightPos).magnitude();
	Color returnColor;

	float s = 1.0f;
	Vector ray = lightPos + s * direction;
	float T = 0.0f;
	while ((point-ray).magnitude() < distance) {
		float res = f.eval(ray.X(), ray.Y(), ray.Z());
		if (res > 0.0f) {
			T += res;
		}
		ray += direction * ds;
		s += ds;
	}

	return T;
}


// Deep shadow maps
ifs::VolumeGrid<float>::VolumeGrid<float>(ScalarField vol, Vector lightPos, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos) {
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;
	this->deltax = deltax;
	this->deltay = deltay;
	this->deltaz = deltaz;
	this->startPos = startPos;

	data = new float* [(Nx * Ny * Nz) / 64];
	for (int p = 0; p < (Nx * Ny * Nz) / 64; p++) {
		data[p] = NULL;
	}

	int finished = 0;

	#pragma omp parallel for schedule(dynamic) num_threads(20) shared(finished) 
	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
			for (int i = 0; i < Nx; i++) {
				Vector pos = Vector(i * deltax + startPos[0], j * deltay + startPos[1], k * deltaz + startPos[2]);
				set(i, j, k, rayMarchLight(vol, pos, lightPos, 0.05f));
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
	float& t, float& u, float& v)
{
	float kEpsilon = 1.0e-20f;

	Vector v0v1 = v1 - v0;
	Vector v0v2 = v2 - v0;
	Vector pvec = dir^v0v2;
	float det = v0v1*pvec;
	// if the determinant is negative, the triangle is 'back facing'
	// if the determinant is close to 0, the ray misses the triangle
	//if (det < kEpsilon) return false;
	// ray and triangle are parallel if det is close to 0
	if (fabs(det) < kEpsilon) return false;
	float invDet = 1 / det;

	Vector tvec = orig - v0;
	u = tvec * pvec * invDet;
	if (u < 0 || u > 1) return false;

	Vector qvec = tvec^v0v1;
	v = dir * qvec * invDet;
	if (v < 0 || u + v > 1) return false;

	t = v0v2 * qvec * invDet;

	return true;
}



VolumeGrid<float> ifs::createGrid(string filename, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos) {

	VolumeGrid<float> grid(Nx, Ny, Nz, deltax, deltay, deltaz, startPos);
	grid.Nx = Nx;
	grid.Ny = Ny;
	grid.Nz = Nz;
	grid.deltax = deltax;
	grid.deltay = deltay;
	grid.deltaz = deltaz;
	grid.startPos = startPos;
	grid.data = new float* [(Nx * Ny * Nz) / 64];
	grid.defaultValue = -1.0e4f;
	for (int p = 0; p < (Nx * Ny * Nz) / 64; p++) {
		grid.data[p] = NULL;
	}

	objl::Loader loader;
	bool data = loader.LoadFile(filename);

	int finished = 0;

	int bandwidth = 20;
	for (int ind = 0; ind < loader.LoadedIndices.size(); ind += 3) {

		Vector p1 = Vector(loader.LoadedVertices[ind].Position.X, loader.LoadedVertices[ind].Position.Y, loader.LoadedVertices[ind].Position.Z);
		Vector p2 = Vector(loader.LoadedVertices[ind + 1].Position.X, loader.LoadedVertices[ind + 1].Position.Y, loader.LoadedVertices[ind + 1].Position.Z);
		Vector p3 = Vector(loader.LoadedVertices[ind + 2].Position.X, loader.LoadedVertices[ind + 2].Position.Y, loader.LoadedVertices[ind + 2].Position.Z);

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
					if (fabs(grid.get(i, j, k)) > d) {
						grid.set(i, j, k, d);
						//printf("(%d,%d,%d) = %f, d = %f\n", i, j, k, grid.get(i, j, k), d);
					}

				}
			}
		}
		finished += 1;
		printf("\rCreating sgd from (%s): %0.2f%%", filename.c_str(), (finished + 1) * 1.0f / (loader.LoadedIndices.size()/3));
	}
	printf("\n");


	finished = 0;

	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
			for (int i = 0; i < Nx; i++) {
				float val = grid.get(i, j, k);
				if (val > 0.0f) {
					grid.set(i, j, k, -val);
				}
			}
		}
	}

	#pragma omp parallel for schedule(dynamic) num_threads(20) shared(finished)
	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {

			int intersactions = 0;
			bool exitwhile = false;
			int numOfIntersactions = 0;
			Vector nextpos = Vector(deltax + startPos[0], j * deltay + startPos[1], k * deltaz + startPos[2]);
			int lasti = 0;
			Vector pos = Vector(startPos[0], j * deltay + startPos[1], k * deltaz + startPos[2]);
			Vector d = (nextpos - pos).unitvector();
			while(!exitwhile){
				bool intersacted = false;
				for (int ind = 0; ind < loader.LoadedIndices.size(); ind += 3) {
					Vector A = Vector(loader.LoadedVertices[ind].Position.X, loader.LoadedVertices[ind].Position.Y, loader.LoadedVertices[ind].Position.Z);
					Vector B = Vector(loader.LoadedVertices[ind + 1].Position.X, loader.LoadedVertices[ind + 1].Position.Y, loader.LoadedVertices[ind + 1].Position.Z);
					Vector C = Vector(loader.LoadedVertices[ind + 2].Position.X, loader.LoadedVertices[ind + 2].Position.Y, loader.LoadedVertices[ind + 2].Position.Z);
					
					float t, u, v;
					intersacted = rayTriangleIntersect(pos, d, A, B, C, t, u, v);

					if (intersacted && t > 0.0f) {
						Vector pointOnTriang = pos + t * d;
						
						int hiti = static_cast<int>(ceil((pointOnTriang.X() - startPos[0]) / deltax));
						//if (k == 100 && j == 100) {
							//printf("(%d,%d) intersected new hiti %d and lasti %d, point on triangle pos:(%f,%f,%f) TriangleInteresection:(%f,%f,%f) value of t %f\n", j, k, hiti, lasti, pos[0], pos[1], pos[2], pointOnTriang[0], pointOnTriang[1], pointOnTriang[2], t);
						//}
						if (numOfIntersactions % 2 == 1) {
							for (int q = lasti; q < hiti; q++) {
								float val = grid.get(q, j, k);
								if (val < 0.0f) {
									grid.set(q, j, k, fabs(val));
								}
							}
						}
						lasti = hiti;
						pos = pointOnTriang;
						pos[0] = pointOnTriang[0] + 0.0001f;
						numOfIntersactions += 1;
						break;
					}
				}
				if (!intersacted) {
					exitwhile = true;
				}

			}

		}
		finished += 1;
		printf("\rChecking intersections (%s): %0.2f%%", filename.c_str(), (finished + 1) * 1.0f / Nz);
	}
	printf("\n");

	return grid;
}