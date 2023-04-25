#pragma once
#include <vector>
#include "Vector.h"
#include "VolumeBase.h"
#include "Color.h"
#include <string>
#include <memory>
#include <math.h>
#include "Camera.h"

using namespace ifs;
using namespace std;

namespace ifs 
{

	struct VolumeParms {
		int Nx;
		int Ny;
		int Nz;
		float dx;
		float dy;
		float dz;
		Vector startPos;
	};


	// Wrapping a function volume
	class VolumeFunc : public VolumeBase<float> {
	public:
		float (*func)(float x, float y, float z);


		VolumeFunc(float (*implicitFunction)(float x, float y, float z)) {
			func = implicitFunction;
		}

		float eval(float x, float y, float z) {
			return func(x, y, z);
		}
	};

	class ScalarField;
	class scalarFieldT;
	class VectorField;
	class ColorField;
	

	template <class T>
	class VolumeGrid : public VolumeBase<T> {
	public:

		float deltax = 0.1f;
		float deltay = 0.1f;
		float deltaz = 0.1f;
		T defaultValue{};
		int Nx, Ny, Nz;
		Vector startPos;
		int blocksize;
		std::unique_ptr<T* []> data;
		Camera c;

		VolumeGrid() {
			Nx = 1;
			Ny = 1;
			Nz = 1;
			startPos = Vector(0.0f, 0.0f, 0.0f);
		}

		VolumeGrid(const VolumeGrid<T> &c) {
			this->Nx = c.Nx;
			this->Ny = c.Ny;
			this->Nz = c.Nz;
			this->deltax = c.deltax;
			this->deltay = c.deltay;
			this->deltaz = c.deltaz;
			this->startPos = c.startPos;
			this->defaultValue = c.defaultValue;
			this->data = c.data.release();
		}

		~VolumeGrid() {
			for (int p = 0; p < (Nx * Ny * Nz) / 64; p++) {
				//cout << "Deleting array at " << p << endl;
				if (data[p] != NULL) {
					delete [] data[p];
				}
			}
			data.reset();
		};
		VolumeGrid(int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);

		VolumeGrid(scalarFieldT vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		VolumeGrid(VectorField vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);

		VolumeGrid(scalarFieldT vol, Vector lightPos, int Nx, int Ny, int Nz, float nearPlane, float farPlane, float fov); // Frustum shaped light deep shadow map
		VolumeGrid(scalarFieldT vol, Vector lightPos, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos); // Light shadow map
		VolumeGrid(const char* filename, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos); // Triangle Mesh to SDF

		T eval(Vector pos) {
			return eval(pos[0], pos[1], pos[2]);
		}

		T eval(float x, float y, float z) {
			Vector p(1.0f * (x - startPos[0]) / deltax, 1.0f * (y - startPos[1]) / deltay, 1.0f * (z - startPos[2]) / deltaz);

			/*int p1 = static_cast<int>(floor(p.X()));
			int p2 = static_cast<int>(floor(p.Y()));
			int p3 = static_cast<int>(floor(p.Z()));*/
			Vector lowerLeft(static_cast<int>(floor(p.X())), static_cast<int>(floor(p.Y())), static_cast<int>(floor(p.Z())));

			/*if (index2(p1, p2, p3) > Nx * Ny * Nz || index2(p1, p2, p3) < 0) {
				return defaultValue;
			}*/

			float weight[3];
			//T value = 0.0f;
			T value{};
			value = value - value;

			for (int i = 0; i < 2; ++i) {
				int cur_x = lowerLeft[0] + i;
				weight[0] = 1.0f - abs(p[0] - cur_x);
				for (int j = 0; j < 2; ++j)
				{
					int cur_y = lowerLeft[1] + j;
					weight[1] = 1.0 - std::abs(p[1] - cur_y);
					for (int k = 0; k <= 1; ++k)
					{
						int cur_z = lowerLeft[2] + k;
						weight[2] = 1.0 - std::abs(p[2] - cur_z);
						//if (index2(cur_x, cur_y, cur_z) < Nx * Ny * Nz && index2(cur_x, cur_y, cur_z) > 0) {
							value += weight[0] * weight[1] * weight[2] * get(cur_x, cur_y, cur_z);
						//}
					}
				}
			}

			return value;
		}

		void splat(Vector pos, T val) {
			Vector p(1.0f * (pos.X() - startPos[0]) / deltax, 1.0f * (pos.Y() - startPos[1]) / deltay, 1.0f * (pos.Z() - startPos[2]) / deltaz);


			Vector lowerLeft(static_cast<int>(floor(p.X())), static_cast<int>(floor(p.Y())), static_cast<int>(floor(p.Z())));


			Vector weight = p - lowerLeft;

			set(lowerLeft[0] , lowerLeft[1], lowerLeft[2], get(lowerLeft[0], lowerLeft[1], lowerLeft[2])
				+ val * (1.0 - weight[0]) * (1.0 - weight[1]) * (1.0 - weight[2]));
			set(lowerLeft[0] + 1, lowerLeft[1], lowerLeft[2], get(lowerLeft[0] + 1, lowerLeft[1], lowerLeft[2])
				+ val * ( weight[0]) * (1.0 - weight[1]) * (1.0 - weight[2]));
			set(lowerLeft[0], lowerLeft[1] + 1, lowerLeft[2], get(lowerLeft[0], lowerLeft[1] + 1, lowerLeft[2])
				+ val * (1.0 - weight[0]) * (weight[1]) * (1.0 - weight[2]));
			set(lowerLeft[0], lowerLeft[1], lowerLeft[2] + 1, get(lowerLeft[0], lowerLeft[1], lowerLeft[2] + 1)
				+ val * (1.0 - weight[0]) * (1.0 - weight[1]) * (weight[2]));
			set(lowerLeft[0] + 1, lowerLeft[1] + 1, lowerLeft[2], get(lowerLeft[0] + 1, lowerLeft[1] + 1, lowerLeft[2])
				+ val * (weight[0]) * (weight[1]) * (1.0 - weight[2]));
			set(lowerLeft[0] + 1, lowerLeft[1], lowerLeft[2] + 1, get(lowerLeft[0] + 1, lowerLeft[1], lowerLeft[2] + 1)
				+ val * (weight[0]) * (1.0 - weight[1]) * (weight[2]));
			set(lowerLeft[0], lowerLeft[1] + 1, lowerLeft[2] + 1, get(lowerLeft[0], lowerLeft[1] + 1, lowerLeft[2] + 1)
				+ val * (1.0 - weight[0]) * (weight[1]) * (weight[2]));
			set(lowerLeft[0] + 1, lowerLeft[1] + 1, lowerLeft[2] + 1, get(lowerLeft[0] + 1, lowerLeft[1] + 1, lowerLeft[2] + 1)
				+ val * (weight[0]) * (weight[1]) * (weight[2]));

		}


		int getBlock(int i, int j, int k) const {
			if (i < 0 || j < 0 || k < 0) {
				return -1;
			}

			int blocki, blockj, blockk;
			blocki = i/4;
			blockj = j/4;
			blockk = k/4;
			int toRet = blocki + blockj * Nx / 4 + blockk * Nx / 4 * Ny / 4;
			if (blocki < 0 || blockj < 0 || blockk < 0 || blocki >= Nx / 4 || blockj >= Ny / 4 || blockk >= Nz / 4 || toRet < 0 || toRet >= blocksize) {
				return -1;
			}
			return toRet;
		}

		int index2(int i, int j, int k) {
			return i + j * Nx + k * Nx * Ny;
		}

		void set(int i, int j, int k, T val) {
			if (index2(i, j, k) > Nx * Ny * Nz || index2(i, j, k) < 0 || val <= defaultValue) {
				return;
			}

			if (getBlock(i, j, k) == -1) {
				return;
			}

			if (data[getBlock(i, j, k)] == NULL) {
				T* dataCube = new T[4 * 4 * 4];
				for (int i = 0; i < 4 * 4 * 4; i++) {
					dataCube[i] = defaultValue;
				}
				data[getBlock(i, j, k)] = dataCube;
				data[getBlock(i, j, k)][(i % 4) + (j % 4) * 4 + (k % 4) * 4 * 4] = val;
			}
			else {
				data[getBlock(i, j, k)][(i % 4) + (j % 4) * 4 + (k % 4) * 4 * 4] = val;
			}

		}

		T get(int i, int j, int k) const{
			//printf("\r%d %d %d", i, j, k);
			int block = getBlock(i, j, k);
			if (block == -1 /*|| index2(i, j, k) > Nx * Ny * Nz || index2(i, j, k) < 0 || i < 0 || j < 0 || k < 0 || i >= Nx || j >= Ny || k >= Nz*/) {
				return defaultValue;
			}
			else if (data[block] == NULL) {
				return defaultValue;
			}
			else {
				return data[block][(i % 4) + (j % 4) * 4 + (k % 4) * 4 * 4];
			}
		}

		T evalFrust(Vector pos) {
			float x, y, z;
			double u, v;
			c.XY(pos, u, v);

			x = (0.5) * ((u / tan(c.fov() / 2)) + 1);
			y = (0.5) * ((v * c.aspectRatio() / tan(c.fov() / 2)) + 1);
			z = ((pos - c.eye()).magnitude() - c.nearPlane()) / (c.farPlane() - c.nearPlane());

			Vector p(x * Nx, y * Ny, z * Nz);
			Vector lowerLeft(static_cast<int>(floor(x*Nx)), static_cast<int>(floor(y*Ny)), static_cast<int>(floor(z*Nz)));

			float weight[3];
			T value{};
			value = value - value;

			for (int i = 0; i < 2; ++i) {
				int cur_x = lowerLeft[0] + i;
				weight[0] = 1.0f - abs(p[0] - cur_x);
				for (int j = 0; j < 2; ++j)
				{
					int cur_y = lowerLeft[1] + j;
					weight[1] = 1.0 - std::abs(p[1] - cur_y);
					for (int k = 0; k <= 1; ++k)
					{
						int cur_z = lowerLeft[2] + k;
						weight[2] = 1.0 - std::abs(p[2] - cur_z);
						value += weight[0] * weight[1] * weight[2] * get(cur_x, cur_y, cur_z);
					}
				}
			}

			return value;
		}

		void setUpCamera(int width, int height, float dnear, float dfar,
			float f, Vector cameraCenter) {

			float a = 1.0f * width / height;

			Camera c;
			c.setAspectRatio(a);
			c.setFov(f);
			c.setEyeViewUp(cameraCenter, -1.0f * cameraCenter, Vector(1.0, 0.0, 0.0));
			c.setFarPlane(dfar);
			c.setNearPlane(dnear);
		}

		void setUpCamera(Camera t) {
			c = t;
		}

		Vector getFrustPos(int i, int j, int k) const {

			float u, v;
			u = (2 * (1.0f * i / Nx) - 1) * tan(c.fov() / 2);
			v = (2 * (1.0f * j / Nx) - 1) * tan(c.fov() / 2)/c.aspectRatio();
			float z = c.nearPlane() + (1.0 * k / Nz) * (c.farPlane() - c.nearPlane());

			Vector q = u * c.right() + v * c.up();

			return c.eye() + (c.view() + q).unitvector() * z;
		}

		friend ostream& operator<<(ostream& os, const VolumeGrid& grid){
			os << grid.deltax;
			os << "\n";
			os << grid.deltay;
			os << "\n";
			os << grid.deltaz;
			os << "\n";
			os << grid.defaultValue;
			os << "\n";
			os << grid.blocksize;
			os << "\n";
			os << grid.Nx;
			os << "\n";
			os << grid.Ny;
			os << "\n";
			os << grid.Nz;
			os << "\n";
			for(int i = 0; i < grid.blocksize; i++){
				if(grid.data[i] == NULL){
					os << 0;
				} else {
					for(int q = 0; q < 4*4*4; q++){
						os << grid.data[i][q];
						os << " ";
					}
				}
				os << "\n";
			}


			return os;
		}


		// friend ifstream &operator>>(ifstream& is, const VolumeGrid<float>& grid){
		// 	string line;
		// 	is.getline(line, 50);
		// 	grid.deltax = stof(line);
		// 	is.getline(line, 50);
		// 	grid.deltay  = stof(line);
		// 	is.getline(line, 50);
		// 	grid.deltaz  = stof(line);
		// 	is.getline(line, 50);
		// 	grid.defaultValue = stof(line);
		// 	is.getline(line, 50);
		// 	grid.blocksize = stoi(line);
		// 	is.getline(line, 50);
		// 	grid.Nx = stoi(line);
		// 	is.getline(line, 50);
		// 	grid.Ny = stoi(line);
		// 	is.getline(line, 50);
		// 	grid.Nz = stoi(line);
			
		// 	// for(int i = 0; i < grid.blocksize; i++){
		// 	// 	char buf[30];
		// 	// 	is >> buf;
		// 	// 	if((buf)
		// 	// 	is >> grid.data[i][q];
		// 	// }

		// 	cout << "Testing the read: " << grid.deltax << grid.Nx << grid.blocksize << endl;


		// 	return is;
		// }
	};


	class VolumeColorGrid {
	public:

		float deltax = 0.1f;
		float deltay = 0.1f;
		float deltaz = 0.1f;
		Color defaultValue = Color(0.0f, 0.0f, 0.0f, 0.0f);
		int Nx, Ny, Nz;
		Vector startPos;
		std::unique_ptr<Color* []> data;

		VolumeColorGrid() {
			Nx = 1;
			Ny = 1;
			Nz = 1;
			startPos = Vector(0.0f, 0.0f, 0.0f);
		}

		VolumeColorGrid(const VolumeColorGrid& c) {
			this->Nx = c.Nx;
			this->Ny = c.Ny;
			this->Nz = c.Nz;
			this->deltax = c.deltax;
			this->deltay = c.deltay;
			this->deltaz = c.deltaz;
			this->startPos = c.startPos;
			this->defaultValue = c.defaultValue;
			//data(new T*[](c.data));
		}

		VolumeColorGrid(ColorField vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);


		Color eval(Vector pos) {
			return eval(pos[0], pos[1], pos[2]);
		}

		Color eval(float x, float y, float z) {
			Vector p(1.0f * (x - startPos[0]) / deltax, 1.0f * (y - startPos[1]) / deltay, 1.0f * (z - startPos[2]) / deltaz);

			int p1 = static_cast<int>(floor(p.X()));
			int p2 = static_cast<int>(floor(p.Y()));
			int p3 = static_cast<int>(floor(p.Z()));
			Vector lowerLeft(static_cast<int>(floor(p.X())), static_cast<int>(floor(p.Y())), static_cast<int>(floor(p.Z())));

			/*if (index2(p1, p2, p3) > Nx * Ny * Nz || index2(p1, p2, p3) < 0) {
				return defaultValue;
			}*/

			float weight[3];
			Color value = Color(0.0f,0.0f,0.0f,0.0f);

			for (int i = 0; i < 2; ++i) {
				int cur_x = lowerLeft[0] + i;
				weight[0] = 1.0f - abs(p[0] - cur_x);
				for (int j = 0; j < 2; ++j)
				{
					int cur_y = lowerLeft[1] + j;
					weight[1] = 1.0 - std::abs(p[1] - cur_y);
					for (int k = 0; k <= 1; ++k)
					{
						int cur_z = lowerLeft[2] + k;
						weight[2] = 1.0 - std::abs(p[2] - cur_z);
						//if (index2(cur_x, cur_y, cur_z) < Nx * Ny * Nz && index2(cur_x, cur_y, cur_z) > 0) {
							value += weight[0] * weight[1] * weight[2] * get(cur_x, cur_y, cur_z);
						//}
					}
				}
			}

			return value;
		}


		int getBlock(int i, int j, int k) const {
			if (i < 0 || j < 0 || k < 0) {
				return -1;
			}

			int blocki, blockj, blockk;
			blocki = i/4;
			blockj = j/4;
			blockk = k/4;
			int toRet = blocki + blockj * Nx / 4 + blockk * Nx / 4 * Ny / 4;
			if (blocki < 0 || blockj < 0 || blockk < 0 || blocki >= Nx / 4 || blockj >= Ny / 4 || blockk >= Nz / 4 || toRet < 0 || toRet >= (Nx*Ny*Nz)/64) {
				return -1;
			}
			return toRet;
		}

		int index2(int i, int j, int k) {
			return i + j * Nx + k * Nx * Ny;
		}

		void set(int i, int j, int k, Color val) {
			if (index2(i, j, k) > Nx * Ny * Nz || index2(i, j, k) < 0 || val == defaultValue) {
				return;
			}

			if (getBlock(i, j, k) == -1) {
				return;
			}

			if (data[getBlock(i, j, k)] == NULL) {
				Color* dataCube = new Color[4 * 4 * 4];
				for (int i = 0; i < 4 * 4 * 4; i++) {
					dataCube[i] = defaultValue;
				}
				data[getBlock(i, j, k)] = dataCube;
				data[getBlock(i, j, k)][(i % 4) + (j % 4) * 4 + (k % 4) * 4 * 4] = val;
			}
			else {
				data[getBlock(i, j, k)][(i % 4) + (j % 4) * 4 + (k % 4) * 4 * 4] = val;
			}

		}

		Color get(int i, int j, int k) {
			int block = getBlock(i, j, k);
			if (block == -1 /*|| index2(i, j, k) > Nx * Ny * Nz || index2(i, j, k) < 0 || i < 0 || j < 0 || k < 0 || i >= Nx || j >= Ny || k >= Nz*/) {
				return defaultValue;
			}
			else if (data[block] == NULL) {
				return defaultValue;
			}
			else {
				return data[block][(i % 4) + (j % 4) * 4 + (k % 4) * 4 * 4];
			}
		}
	};
}
