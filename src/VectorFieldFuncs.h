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

	class ConstVectorField : public FieldBase<Vector> {
	public:
		Vector C;

		ConstVectorField(Vector c) : C(c) {};

		const Vector eval(const Vector& pos) const {
			return C;
		};
	};

	class VectorTimesScalarField : public FieldBase<Vector> {
	public:
		const VectorField e1;
		const scalarFieldT e2;

		VectorTimesScalarField(const VectorField& V, const scalarFieldT& F) : e1(V), e2(F) {};

		const Vector eval(const Vector& pos) const {
			return e2->eval(pos) * e1->eval(pos);
		};
	};

	class FuncVectorField : public FieldBase<Vector> {
	public:
		Vector (*func)(float x, float y, float z);


		FuncVectorField(Vector (*implicitFunction)(float x, float y, float z)) {
			func = implicitFunction;
		}

		const Vector eval(const Vector& pos) const{
			return func(pos.X(), pos.Y(), pos.Z());
		}
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

	class MultiVectorField : public FieldBase<Vector> {
	public:
		float constant;
		const VectorField a;

		MultiVectorField(const VectorField& a, const float& constant);

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

	class NoiseVectorField : public FieldBase<Vector> {
	public:
		FSPNParms parms;

		NoiseVectorField(FSPNParms parms) : parms(parms) {};

		const Vector eval(const Vector& pos) const;

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

	class GridVectorField : public FieldBase<Vector> {
	public:
		VolumeGrid<Vector>* grid;

		//~GridVectorField() {
		//	delete grid;
		//}

		GridVectorField(VectorField& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
		GridVectorField(VolumeGrid<Vector>* thegrid) {
			grid = new VolumeGrid<Vector>(thegrid->Nx, thegrid->Ny, thegrid->Nz, thegrid->deltax, thegrid->deltay, thegrid->deltaz, thegrid->startPos);
			grid->defaultValue = thegrid->defaultValue;

			for (int k = 0; k < thegrid->Nz; k++) {
				for (int j = 0; j < thegrid->Ny; j++) {
					for (int i = 0; i < thegrid->Nx; i++) {
						grid->set(i, j, k, thegrid->get(i, j, k));
					} 
				}
			}
		};

		const Vector eval(const Vector& pos) const;
	};

	//class incpompressedVectorField : public FieldBase<Vector> {
	//public:
	//	VolumeGrid<Vector>* grid;

	//	~incpompressedVectorField() {
	//		delete grid;
	//	}

	//	incpompressedVectorField(VolumeGrid<Vector>& a){
	//		
	//	
	//	};

	//	const Vector eval(const Vector& pos) const {
	//		return (grid)->grad(pos);
	//	};
	//};

}