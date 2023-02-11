#pragma once

namespace ifs {

	template <class T>
	class FieldBase {
	public:

		T eval(float x, float y, float z) {
			return eval(Vector(x, y, z));
		}

		virtual T eval(Vector pos) {
			printf("Calling the virtual function, error.\n");
			return T(0.0f, 0.0f, 0.0f, 0.0f);
		}

		FieldBase* add(FieldBase* second);
		FieldBase* operator-(FieldBase* second);
		FieldBase* operator+(FieldBase * second);



	};


}