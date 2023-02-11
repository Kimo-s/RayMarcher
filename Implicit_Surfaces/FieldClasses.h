#pragma once
#include "IFSclasses.h"
#include "Color.h"
#include "FieldBase.h"
#include <iostream>

using namespace std;


namespace ifs {

	// Operations on color fields
	template <class T>
	class AddColor : public FieldBase<T> {
	public:
		FieldBase<T>* a;
		FieldBase<T>* b;

		AddColor(FieldBase<T>* a, FieldBase<T>* b) {
			this->a = a;
			this->b = b;
			/*a.eval(Vector(0.0f, 0.0f, 0.0f));
			b.eval(Vector(0.0f, 0.0f, 0.0f));*/
		}

		T eval(Vector pos) {
			//printf("Testing A.eval:%s B.eval:%s\n", a->eval(pos).print(), b->eval(pos).print());
			return a->eval(pos) + b->eval(pos);
		}

	};

	template <class T>
	class SubBase : public FieldBase<T> {
	public:
		FieldBase<T>* a;
		FieldBase<T>* b;

		SubBase(FieldBase<T>* a, FieldBase<T>* b) {
			this->a = a;
			this->b = b;
		}

		T eval(Vector pos) {
			return a->eval(pos) - b->eval(pos);
		}

	};

	class ColorFieldMask : public FieldBase<Color> {
	public:

		ScalarField f;
		Color col;

		ColorFieldMask() {
		}

		ColorFieldMask(const ColorFieldMask& field) {
			this->f = field.f;
			this->col = field.col;
		}

		ColorFieldMask(ScalarField f, Color col) {
			this->f = f;
			this->col = col;
		}

		Color eval(Vector pos) {
			if (f.eval(pos) > 0.0f) {
				return col;
			}
			else {
				return Color(0.0f, 0.0f, 0.0f, 0.0f);
			}
		}	

	};


}