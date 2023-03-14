#pragma once
#define sq(A) (A*A)
#include "FieldClasses.h"
#include "VolumeClasses.h"

namespace ifs {

	class IdentityVectorField : public FieldBase<Vector> {
	public:

		IdentityVectorField() {};

		const Vector eval(const Vector& pos) const {
			return pos;
		};
	};

	class AddVectorFields : public FieldBase<Vector> {
	public:
		AddVectorFields(const VectorField& a, const VectorField& b);

		const Vector eval(const Vector& pos) const;
	private:
		const VectorField e1;
		const VectorField e2;
	};

	class gradVectorField : public FieldBase<Vector> {
	public:
		const scalarFieldT a;

		gradVectorField(const scalarFieldT& a): a(a) {};

		const Vector eval(const Vector& pos) const {
			return a->grad(pos);
		};
	};
	

	class SubVectorFields : public FieldBase<Vector> {
	public:
		SubVectorFields(const VectorField& a, const VectorField& b);

		const Vector eval(const Vector& pos) const;
	private:
		const VectorField e1;
		const VectorField e2;
	};

	class ScaleVectorField : public FieldBase<Vector> {
	public:
		float scalefactor;
		const VectorField e2;

		ScaleVectorField(const VectorField& a, const float& scalefactor);

		const Vector eval(const Vector& pos) const;

	};

	class warpVectorField : public FieldBase<Vector> {
	public:
		const VectorField F;
		const VectorField V;

		warpVectorField(const VectorField& a, const VectorField& b) : F(a), V(b) {};

		const Vector eval(const Vector& pos) const {
			return F->eval(V->eval(pos));
		};

	};

	class TranslateVectorField : public FieldBase<Vector> {
	public:
		Vector v;
		const VectorField e2;

		TranslateVectorField(const VectorField& a, const Vector& v);

		const Vector eval(const Vector& pos) const;

	};

	class RotateVectorField : public FieldBase<Vector> {
	public:
		const VectorField a;
		//float al, b, g;
		const Vector u;
		float theta;
		float m11, m12, m13, m21, m22, m23, m31, m32, m33;

		RotateVectorField(const VectorField& a, const Vector& u, float theta);

		const Vector eval(const Vector& pos) const;

	};

	class MaxVectorField : public FieldBase<Vector> {
	public:
		const VectorField a;
		const VectorField b;

		MaxVectorField(const VectorField& a, const VectorField& b);

		const Vector eval(const Vector& pos) const;

	};

	class CutVectorField : public FieldBase<Vector> {
	public:
		const VectorField a;
		const VectorField b;

		CutVectorField(const VectorField& a, const VectorField& b);

		const Vector eval(const Vector& pos) const;

	};

	class MaskVectorField : public FieldBase<Vector> {
	public:
		const VectorField a;

		MaskVectorField(const VectorField& a);

		const Vector eval(const Vector& pos) const;

	};

	//class FuncScalarField : public FieldBase<float> {
	//public:
	//	float (*func)(float x, float y, float z);


	//	FuncScalarField(float (*implicitFunction)(float x, float y, float z)) {
	//		func = implicitFunction;
	//	}

	//	float eval(float x, float y, float z) {
	//		return eval(Vector(x, y, z));
	//	}

	//	const float eval(const Vector& pos) const {
	//		return func(pos.X(), pos.Y(), pos.Z());
	//	}
	//};

	/*class GridScalarField : public FieldBase<float> {
	public:
		VolumeGrid<float>* grid;

		GridScalarField(int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		GridScalarField(scalarFieldT& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		GridScalarField(scalarFieldT& vol, Vector lightPos, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		GridScalarField(const char* filename, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);

		const float eval(const Vector& pos) const;
	};*/

}