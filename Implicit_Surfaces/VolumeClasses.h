#pragma once
#include <vector>
#include "Vector.h"
#include "VolumeBase.h"
#include "Color.h"
#include <string>
#include <memory>

using namespace ifs;
using namespace std;

namespace ifs 
{
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
	class ColorField;
	

	template <class T>
	class VolumeGrid : public VolumeBase<T> {
	public:

		float deltax = 0.1f;
		float deltay = 0.1f;
		float deltaz = 0.1f;
		T defaultValue = 0.0f;
		int Nx, Ny, Nz;
		Vector startPos;
		T** data;

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
			//data(new T*[](c.data));
		}


		VolumeGrid(int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		VolumeGrid(scalarFieldT vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		VolumeGrid(scalarFieldT vol, Vector lightPos, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		VolumeGrid(const char* filename, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);

		T eval(Vector pos) {
			return eval(pos[0], pos[1], pos[2]);
		}

		T eval(float x, float y, float z) {
			Vector p(1.0f * (x - startPos[0]) / deltax, 1.0f * (y - startPos[1]) / deltay, 1.0f * (z - startPos[2]) / deltaz);

			int p1 = static_cast<int>(floor(p.X()));
			int p2 = static_cast<int>(floor(p.Y()));
			int p3 = static_cast<int>(floor(p.Z()));
			Vector lowerLeft(static_cast<int>(floor(p.X())), static_cast<int>(floor(p.Y())), static_cast<int>(floor(p.Z())));

			if (index2(p1, p2, p3) > Nx * Ny * Nz || index2(p1, p2, p3) < 0) {
				return defaultValue;
			}

			float weight[3];
			T value = 0.0f;

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
						if (index2(cur_x, cur_y, cur_z) < Nx * Ny * Nz && index2(cur_x, cur_y, cur_z) > 0) {
							value += weight[0] * weight[1] * weight[2] * get(cur_x, cur_y, cur_z);
						}
					}
				}
			}

			return value;
		}


		int getBlock(int i, int j, int k) {
			int toRet = i + j * Nx / 4 + k * Nx / 4 * Ny / 4;
			if (i < 0 || j < 0 || k < 0 || i > Nx / 4 || j > Ny / 4 || k > Nz / 4 || toRet < 0 || toRet > ((Nx * Ny * Nz) / 64)) {
				return -1;
			}
			return i + j * Nx / 4 + k * Nx / 4 * Ny / 4;
		}

		int index2(int i, int j, int k) {
			return i + j * Nx + k * Nx * Ny;
		}

		void set(int i, int j, int k, T val) {
			if (index2(i, j, k) > Nx * Ny * Nz || index2(i, j, k) < 0 || val <= defaultValue) {
				return;
			}

			if (getBlock(i / 4, j / 4, k / 4) == -1) {
				return;
			}

			if (data[getBlock(i / 4, j / 4, k / 4)] == NULL) {
				T* dataCube = new T[4 * 4 * 4];
				for (int i = 0; i < 4 * 4 * 4; i++) {
					dataCube[i] = defaultValue;
				}
				data[getBlock(i / 4, j / 4, k / 4)] = dataCube;
				data[getBlock(i / 4, j / 4, k / 4)][(i % 4) + (j % 4) * 4 + (k % 4) * 4 * 4] = val;
			}
			else {
				data[getBlock(i / 4, j / 4, k / 4)][(i % 4) + (j % 4) * 4 + (k % 4) * 4 * 4] = val;
			}

		}

		T get(int i, int j, int k) {
			if (getBlock(i / 4, j / 4, k / 4) == -1) {
				return defaultValue;
			}
			else if (data[getBlock(i / 4, j / 4, k / 4)] == NULL) {
				return defaultValue;
			}
			else {
				return data[getBlock(i / 4, j / 4, k / 4)][(i % 4) + (j % 4) * 4 + (k % 4) * 4 * 4];
			}
		}
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

			if (index2(p1, p2, p3) > Nx * Ny * Nz || index2(p1, p2, p3) < 0) {
				return defaultValue;
			}

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
						if (index2(cur_x, cur_y, cur_z) < Nx * Ny * Nz && index2(cur_x, cur_y, cur_z) > 0) {
							value += weight[0] * weight[1] * weight[2] * get(cur_x, cur_y, cur_z);
						}
					}
				}
			}

			return value;
		}


		int getBlock(int i, int j, int k) {
			if (i < 0 || j < 0 || k < 0 || i > Nx / 4 || j > Ny / 4 || k > Nz / 4) {
				return -1;
			}
			//int ind = i + j * Nx / 4 + k * Nx / 4 * Ny / 4;
			//if (ind > Nx * Ny * Nz / 64) {
			//	return -1;
			//}
			return i + j * Nx / 4 + k * Nx / 4 * Ny / 4;
		}

		int index2(int i, int j, int k) {
			return i + j * Nx + k * Nx * Ny;
		}

		void set(int i, int j, int k, Color val) {
			if (index2(i, j, k) > Nx * Ny * Nz || index2(i, j, k) < 0 || val == defaultValue) {
				return;
			}

			if (data[getBlock(i / 4, j / 4, k / 4)] == NULL) {
				Color* dataCube = new Color[4 * 4 * 4];
				for (int i = 0; i < 4 * 4 * 4; i++) {
					dataCube[i] = defaultValue;
				}
				data[getBlock(i / 4, j / 4, k / 4)] = dataCube;
				data[getBlock(i / 4, j / 4, k / 4)][(i % 4) + (j % 4) * 4 + (k % 4) * 4 * 4] = val;
			}
			else {
				data[getBlock(i / 4, j / 4, k / 4)][(i % 4) + (j % 4) * 4 + (k % 4) * 4 * 4] = val;
			}

		}

		Color get(int i, int j, int k) {
			/*if (index2(i, j, k) > Nx * Ny * Nz || index2(i, j, k) < 0) {
				return defaultValue;
			}*/
			if (getBlock(i / 4, j / 4, k / 4) == -1) {
				return defaultValue;
			}
			else if (data[getBlock(i / 4, j / 4, k / 4)] == NULL) {
				return defaultValue;
			}
			else {
				return data[getBlock(i / 4, j / 4, k / 4)][(i % 4) + (j % 4) * 4 + (k % 4) * 4 * 4];
			}
		}
	};
}
