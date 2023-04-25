#pragma once
#define sq(A) (A*A)
#include "FieldClasses.h"
#include "VolumeClasses.h"

namespace ifs {

	float evalFSPN(FSPNParms parm, Vector pos);

	class constScalarField : public FieldBase<float> {
	public:
		float C;

		constScalarField(float c) : C(c) {};

		const float eval(const Vector& pos) const {
			return C;
		};
	};

	class planeScalarField : public FieldBase<float> {
	public:
		Vector normal, center;
		//float thickness, width, height;
		FSPNParms params;
		FSPNParms paramsNegtive;

		planeScalarField(Vector normal, Vector center, FSPNParms params) :
			normal(normal.unitvector()),
			center(center),
			params(params)
		{

			paramsNegtive = {
				2.3f,
				params.N,
				1.4f,
				7.8f,
				0.7f,
				params.xt,
				0.4f
			};
		
		};

		const float eval(const Vector& pos) const;
	};


	class AddScalarFields : public FieldBase<float> {
	public:
		AddScalarFields(const scalarFieldT& a, const scalarFieldT& b);

		const float eval(const Vector& pos) const;
	private:
		const scalarFieldT e1;
		const scalarFieldT e2;
	};

	class MultiplyScalarFields : public FieldBase<float> {
	public:
		MultiplyScalarFields(const scalarFieldT& a, const scalarFieldT& b);

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

	class warpScalarField : public FieldBase<float> {
	public:
		const scalarFieldT F;
		const VectorField V;

		warpScalarField(const scalarFieldT& a, const VectorField& b) : F(a), V(b) {};

		const float eval(const Vector& pos) const {
			return F->eval(V->eval(pos));
		};

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

	class CutScalarField : public FieldBase<float> {
	public:
		const scalarFieldT a;
		const scalarFieldT b;

		CutScalarField(const scalarFieldT& a, const scalarFieldT& b);

		const float eval(const Vector& pos) const;

	};

	class MaskScalarField : public FieldBase<float> {
	public:
		const scalarFieldT a;

		MaskScalarField(const scalarFieldT& a);

		const float eval(const Vector& pos) const;

	};

	class wispScalarField : public FieldBase<float> {
	public:
		VolumeGrid<float>* grid;
		FSPNParms fspn1;
		FSPNParms fspn2;

		~wispScalarField() {
			delete grid;
		}
		wispScalarField(VolumeParms parms, FSPNParms fspn1, FSPNParms fspn2, Vector guidPos, float density, float pscale, float wisp_displacement, float clump, int wispCount);

		const float eval(const Vector& pos) const;

	};

	class pyroclasticScalarField : public FieldBase<float> {
	public:
		const scalarFieldT a;
		FSPNParms params;

		pyroclasticScalarField(const scalarFieldT& a, FSPNParms params);

		const float eval(const Vector& pos) const;

	};

	class addGuideParticaleScalarField : public FieldBase<float> {
	public:
		VolumeGrid<float>* grid;

		~addGuideParticaleScalarField() {
			delete grid;
		}
		addGuideParticaleScalarField(Vector u, FSPNParms params, float fade, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz);

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

		~GridScalarField() {
			delete grid;
		}
		GridScalarField(int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		GridScalarField(VolumeGrid<float>* thegrid) {
			grid = new VolumeGrid<float>(thegrid->Nx, thegrid->Ny, thegrid->Nz, thegrid->deltax, thegrid->deltay, thegrid->deltaz, thegrid->startPos);
			grid->defaultValue = thegrid->defaultValue;

			for (int k = 0; k < thegrid->Nz; k++) {
				for (int j = 0; j < thegrid->Ny; j++) {
					for (int i = 0; i < thegrid->Nx; i++) {
						grid->set(i, j, k, thegrid->get(i, j, k));
					}
				}
			}
		};
		GridScalarField(scalarFieldT& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		GridScalarField(scalarFieldT& vol, Vector lightPos, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		GridScalarField(const char* filename, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);

		const float eval(const Vector& pos) const;
	};

}