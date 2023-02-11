#pragma once
#include <vector>
#include "Vector.h"
#include <algorithm> 
#include "VolumeClasses.h"
#include <omp.h>

using namespace ifs;
using namespace std;

namespace ifs 
{

	struct Point {
		float x;
		float y;
		float z;
	};

	enum class CSGop
	{
		MIN,
		MAX,
		CUT,
		TRANSLATE,
		SCALE,
		MASK,
		ADD,
		None,
	};




	class ScalarField
	{
	public:
		ScalarField(VolumeBase<float>* vol, Vector col);
		ScalarField(float (*func)(float x, float y, float z));
		ScalarField(float (*func)(float x, float y, float z), Vector col);
		ScalarField() {
			volume = NULL;
			g = NULL;
			f = NULL;
			operation = CSGop::None;
			translateVector = { 0.0f,0.0f,0.0f };
		}
		ScalarField(const ScalarField &field);
		~ScalarField();


		float eval(float x, float y, float z) const {
			float result;

			if (operation == CSGop::None) {
				result = volume->eval(x, y, z);
				//printf("Value at (%f,%f,%f): %f\n", x, y, z, result);
			}
			else if (operation == CSGop::MIN) {
				result = std::min(f->eval(x, y, z), g->eval(x, y, z));
				//printf("(min)Value at (%f,%f,%f): %f\n", x, y, z, result);
			}
			else if (operation == CSGop::MAX) {
				result = std::max(f->eval(x, y, z), g->eval(x, y, z));
				//printf("(max)Value at (%f,%f,%f): %f\n", x, y, z, result);

			}
			else if (operation == CSGop::CUT) {
				result = std::min(f->eval(x, y, z), -g->eval(x, y, z));
				//printf("(cut)Value at (%f,%f,%f): %f\n", x, y, z, result);
			}
			else if (operation == CSGop::TRANSLATE) {
				result = f->eval(x + translateVector.x, y + translateVector.y, z + translateVector.z);
				//printf("(translate)Value at (%f,%f,%f): %f\n", x, y, z, result);
			}
			else if (operation == CSGop::SCALE) {
				result = f->eval(x / scaleFactor, y / scaleFactor, z / scaleFactor);
				//printf("(scale)Value at (%f,%f,%f): %f\n", x, y, z, result);
			}
			else if (operation == CSGop::MASK) {
				result = f->eval(x, y, z) < 0.0f ? 0.0f : 1.0f;
			}
			else if (operation == CSGop::ADD) {
				result = f->eval(x, y, z) + g->eval(x, y, z);
			}
			else {
				return volume->eval(x, y, z);
			}

			return result;
		}

		float operator()(float x, float y, float z) {
			return eval(x, y, z);
		}

		float operator()(ifs::Point pointToEval) {
			return eval(pointToEval.x, pointToEval.y, pointToEval.z);
		}

		float eval(ifs::Point pointToEval) {
			return eval(pointToEval.x, pointToEval.y, pointToEval.z);
		}

		float eval(Vector v) {
			return eval(v.X(), v.Y(), v.Z());
		}
		
		Vector evalColor(float x, float y, float z) {
			if (operation == CSGop::None) {
				return color;
			}
			else if (operation == CSGop::MASK){
				return f->evalColor(x, y, z);
			}
			else if (operation == CSGop::MIN) {
				if (f->eval(x, y, z) < g->eval(x, y, z)) {
					return f->evalColor(x, y, z);
				}
				else {
					return g->evalColor(x, y, z);
				}
			}
			else if (operation == CSGop::MAX) {
				if (f->eval(x, y, z) > g->eval(x, y, z)) {
					return f->evalColor(x, y, z);
				}
				else {
					return g->evalColor(x, y, z);
				}
			}
			else if (operation == CSGop::CUT) {
				if (f->eval(x, y, z) > -g->eval(x, y, z)) {
					return f->evalColor(x, y, z);
				}
				else {
					return g->evalColor(x, y, z);
				}
			}
			else {
				return color;
			}
		}

		Vector evalColor(Point p) {
			return evalColor(p.x, p.y, p.z);
		}

		ScalarField min(ScalarField g) {
			ScalarField h = init_(&g);
			h.operation = CSGop::MIN;
			return h;
		}

		ScalarField max(ScalarField g) {
			ScalarField h = init_(&g);
			h.operation = CSGop::MAX;
			return h;
		}

		ScalarField add(ScalarField g) {
			ScalarField h = init_(&g);
			h.operation = CSGop::ADD;
			return h;
		}

		ScalarField operator+(ScalarField g) {
			return this->add(g);
		}

		ScalarField cut(ScalarField g) {
			ScalarField h = init_(&g);
			h.operation = CSGop::CUT;
			return h;
		}

		ScalarField mask() {
			ScalarField h;
			h.f = new ScalarField(*this);
			h.operation = CSGop::MASK;
			h.color = this->color;
			return h;
		}

		ScalarField translate(Point p) {
			ScalarField h;
			h.f = new ScalarField(*this);
			h.operation = CSGop::TRANSLATE;
			h.translateVector = {p.x, p.y, p.z};
			h.color = h.f->color;
			return h;
		}

		ScalarField scale(Point p, float factor) {
			ScalarField h;

			h.operation = CSGop::SCALE;
			h.scaleFactor = factor;
			h.f = new ScalarField(this->translate(p));

			h = h.translate(Point { -p.x, -p.y, -p.z });
			h.color = this->color;

			return h;
		}

		void setColor(Vector p) {
			color = p;
		}

		Vector color = Vector(0.0f, 0.0f, 0.0f);

	private:
		//float (*func)(float x, float y, float z);
		VolumeBase<float>* volume;
		ScalarField* f;
		ScalarField* g;
		CSGop operation;

		Point translateVector;
		float scaleFactor = 1.0f;


		ScalarField init_(ScalarField* g) {
			ScalarField h;
			h.f = new ScalarField(*this);
			h.g = new ScalarField(*g);
			return h;
		}

	};


}
