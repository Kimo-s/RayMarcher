#pragma once
#include "IFSclasses.h"
#include "Color.h"
#include "VolumeBase.h"
#include <iostream>
#define sq(A) (A*A)
using namespace std;


namespace ifs {

	template <typename T>
	class FieldBase {
	public:
		typedef T volumeDataType;

		volumeDataType eval(float x, float y, float z) const {
			return eval(Vector(x, y, z));
		}

		virtual const volumeDataType eval(const Vector& pos) const = 0; // {
		//	volumeDataType base; 
		//	//base = 0; 
		//	return base;
		//}

	};

	// Operations on color fields
	//template <typename T>
	//class AddFields : public FieldBase<T> {
	//public:
	//	FieldBase<T>* a;
	//	FieldBase<T>* b;

	//	AddFields(FieldBase<T>* a, FieldBase<T>* b) {
	//		this->a = a;
	//		this->b = b;
	//		/*a.eval(Vector(0.0f, 0.0f, 0.0f));
	//		b.eval(Vector(0.0f, 0.0f, 0.0f));*/
	//	}

	//	T eval(const Vector pos) const {
	//		//printf("Testing A.eval:%s B.eval:%s\n", a->eval(pos).print(), b->eval(pos).print());
	//		return a->eval(pos) + b->eval(pos);
	//	}

	//};

	//template <typename T>
	//class SubBase : public FieldBase<T> {
	//public:
	//	FieldBase<T>* a;
	//	FieldBase<T>* b;

	//	SubBase(FieldBase<T>* a, FieldBase<T>* b) {
	//		this->a = a;
	//		this->b = b;
	//	}

	//	T eval(const Vector pos) const {
	//		return a->eval(pos) - b->eval(pos);
	//	}

	//};

	//template <typename T>
	//class TranslateField : public FieldBase<T> {
	//public:
	//	FieldBase<T>* a;
	//	Vector translateVec;

	//	TranslateField(FieldBase<T>* a, Vector translateVec) {
	//		this->a = a;
	//		this->translateVec = translateVec;
	//	}

	//	T eval(const Vector pos) const {
	//		return a->eval(pos - translateVec);
	//	}

	//};

	//template <typename T>
	//class RotateField : public FieldBase<T> {
	//public:
	//	FieldBase<T>* a;
	//	//float al, b, g;
	//	Vector u;
	//	float theta;
	//	float m11, m12, m13, m21, m22, m23, m31, m32, m33;

	//	RotateField(FieldBase<T>* a, Vector u, float theta) {
	//		this->a = a;
	//		this->u = u.unitvector();
	//		this->theta = theta;
	//		/*this->al = al;
	//		this->b = b;
	//		this->gamma = g;*/
	//		m11 = (cos(theta) + sq(u.X()) * (1.0f - cos(theta)));
	//		m12 = ((u.X() * u.Y() * (1.0f - cos(theta)) - u.Z() * sin(theta)));
	//		m13 = ((u.X() * u.Z() * (1.0f - cos(theta)) + u.Y() * sin(theta)));

	//		m21 = ((u.X() * u.Y() * (1.0f - cos(theta)) + u.Z() * sin(theta)));
	//		m22 = (cos(theta) + sq(u.Y()) * (1.0f - cos(theta)));
	//		m23 = ((u.Y() * u.Z() * (1.0f - cos(theta)) - u.X() * sin(theta)));

	//		m31 = ((u.Z() * u.X() * (1.0f - cos(theta)) - u.Y() * sin(theta)));
	//		m32 = ((u.Z() * u.Y() * (1.0f - cos(theta)) + u.X() * sin(theta)));
	//		m33 = (cos(theta) + sq(u.Z()) * (1.0f - cos(theta)));
	//	}

	//	T eval(const Vector pos) const{
	//		/*pos.xyz[0] = cos(al) * cos(b) * pos.X() + (cos(al) * sin(b) * sin(g) - sin(al) * cos(g)) * pos.Y() + (cos(al) * sin(b) * cos(g) + sin(al) * sin(g)) * pos.Z();
	//		pos.xyz[1] = sin(al) * cos(b) * pos.X() + (sin(al) * sin(b) * sin(g) + cos(al) * cos(g)) * pos.Y() + (sin(al) * sin(b) * cos(g) - cos(al) * sin(g)) * pos.Z();
	//		pos.xyz[2] = -sin(b) * pos.X() + cos(b) * sin(g) * pos.Y() + cos(b) * cos(g) * posZ();*/
	//		float xc = m11 * pos.X() + m12 * pos.Y() + m13 * pos.Z();
	//		float yc = m21 * pos.X() + m22 * pos.Y() + m23 * pos.Z();
	//		float zc = m31 * pos.X() + m32 * pos.Y() + m33 * pos.Z();
	//		return a->eval(Vector(xc, yc, zc));
	//	}

	//};

	//template <typename T>
	//class ScaleField : public FieldBase<T> {
	//public:
	//	FieldBase<T>* a;
	//	float scaleFactor;

	//	ScaleField(FieldBase<T>* a, float scaleFactor) {
	//		this->a = a;
	//		this->scaleFactor = scaleFactor;
	//	}

	//	T eval(const Vector pos) const {
	//		return a->eval(pos/scaleFactor);
	//	}

	//};




	class scalarFieldT : public std::shared_ptr<FieldBase<float> > {
	public:

		//const VolumeBase<float> volume;


		scalarFieldT(): std::shared_ptr<FieldBase<float> >() {
		}
		scalarFieldT(FieldBase<float>* f): std::shared_ptr<FieldBase<float> >(f) {
		};
		~scalarFieldT() {};

		/*const float eval(const Vector& pos) const {
			return volume.eval(pos);
		}*/

		static scalarFieldT add(scalarFieldT& v1, const scalarFieldT& v2);\
		scalarFieldT max(const scalarFieldT& v2);
		scalarFieldT translate(const Vector transVec);
		scalarFieldT scale(float scaleFactor);
		scalarFieldT mask();
		scalarFieldT rotate(Vector u, float theta);
		scalarFieldT operator+(const scalarFieldT& second);

	};

	class ColorField : public std::shared_ptr<FieldBase<Color> > {
	public:

		~ColorField() {
		}

		ColorField() {
		}

		ColorField(FieldBase<Color>* f) : std::shared_ptr<FieldBase<Color> >(f) {
			//(ColorFieldMask*)f.f
		}

		/*Color eval(const Vector& pos) const {
			if (f->eval(pos) > 0.0f) {
				return col;
			}
			else {
				return Color(0.0f, 0.0f, 0.0f, 0.0f);
			}
		}*/

		ColorField add(ColorField& v1, const ColorField& v2);
		ColorField translate(const Vector transVec);
		ColorField scale(float scaleFactor);
		ColorField rotate(Vector u, float theta);
		ColorField operator+(const ColorField& second);
	};

	scalarFieldT funcField(float(*func)(float x, float y, float z));
	ColorField gridColorField(ColorField& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos); //Will be implmented in the future
	ColorField colorMaskField(const scalarFieldT& f, Color col);
	scalarFieldT gridField(int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
	scalarFieldT gridField(scalarFieldT& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
	scalarFieldT gridField(scalarFieldT& vol, Vector lightPos, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
	scalarFieldT gridField(const char* filename, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
}