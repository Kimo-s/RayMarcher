#pragma once
#include "FieldClasses.h"
#define sq(A) (A*A)

namespace ifs {


	class AddColorFields : public FieldBase<Color> {
	public:
		const ColorField e1;
		const ColorField e2;

		AddColorFields(const ColorField& a, const ColorField& b);

		const Color eval(const Vector& pos) const;

	};

	class ScaleColorField : public FieldBase<Color> {
	public:
		float scalefactor;
		ColorField e2;

		ScaleColorField(const ColorField& a, const float& scalefactor);

		const Color eval(const Vector& pos) const;

	};

	class TranslateColorField : public FieldBase<Color> {
	public:
		Vector v;
		ColorField e2;

		TranslateColorField(const ColorField& a, const Vector& v);

		const Color eval(const Vector& pos) const;

	};

	class RotateColorField : public FieldBase<Color> {
	public:
		ColorField a;
		//float al, b, g;
		Vector u;
		float theta;
		float m11, m12, m13, m21, m22, m23, m31, m32, m33;

		RotateColorField(const ColorField& a, const Vector& u, float theta);

		const Color eval(const Vector& pos) const;

	};

	class ColorMaskField : public FieldBase<Color> {
	public:
		const scalarFieldT f;
		Color col;

		ColorMaskField(const scalarFieldT& f, Color col);

		const Color eval(const Vector& pos) const;
	};

	class GridColorField : public FieldBase<Color> {
	public:
		VolumeColorGrid* grid;

		GridColorField(ColorField& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);


		const Color eval(const Vector& pos) const;
	};
}