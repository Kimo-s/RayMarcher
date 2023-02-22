#pragma once
#define sq(A) (A*A)
#include "FieldClasses.h"
#include "VolumeClasses.h"

namespace ifs {

	class AddScalarFields : public FieldBase<float> {
	public:
		AddScalarFields(const scalarFieldT& a, const scalarFieldT& b);

		const float eval(const Vector& pos) const;
	private:
		const scalarFieldT e1;
		const scalarFieldT e2;
	};

	class ScaleScalarField : public FieldBase<float> {
	public:
		float scalefactor;
		const scalarFieldT e2;

		ScaleScalarField(const scalarFieldT& a, const float& scalefactor);

		const float eval(const Vector& pos) const;

	};

	class TranslateScalarField : public FieldBase<float> {
	public:
		Vector v;
		const scalarFieldT e2;

		TranslateScalarField(const scalarFieldT& a, const Vector& v);

		const float eval(const Vector& pos) const;

	};

	class RotateScalarField : public FieldBase<float> {
	public:
		const scalarFieldT a;
		//float al, b, g;
		const Vector u;
		float theta;
		float m11, m12, m13, m21, m22, m23, m31, m32, m33;

		RotateScalarField(const scalarFieldT& a, const Vector& u, float theta);

		const float eval(const Vector& pos) const;

	};

	class MaxScalarField : public FieldBase<float> {
	public:
		const scalarFieldT a;
		const scalarFieldT b;

		MaxScalarField(const scalarFieldT& a, const scalarFieldT& b);

		const float eval(const Vector& pos) const;

	};

	class MaskScalarField : public FieldBase<float> {
	public:
		const scalarFieldT a;

		MaskScalarField(const scalarFieldT& a);

		const float eval(const Vector& pos) const;

	};

	class FuncScalarField : public FieldBase<float> {
	public:
		float (*func)(float x, float y, float z);


		FuncScalarField(float (*implicitFunction)(float x, float y, float z)) {
			func = implicitFunction;
		}

		float eval(float x, float y, float z) {
			return eval(Vector(x,y,z));
		}

		const float eval(const Vector& pos) const{
			return func(pos.X(), pos.Y(), pos.Z());
		}
	};

	class GridScalarField : public FieldBase<float> {
	public:
		VolumeGrid<float>* grid;

		GridScalarField(int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		GridScalarField(scalarFieldT& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		GridScalarField(scalarFieldT& vol, Vector lightPos, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		GridScalarField(const char* filename, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);

		const float eval(const Vector& pos) const;
	};

}