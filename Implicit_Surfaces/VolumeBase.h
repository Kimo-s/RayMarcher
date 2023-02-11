#pragma once
#include "Vector.h"

using namespace ifs;

namespace ifs {

	// The base abstract volume class
	template <class T>
	class VolumeBase {
	public:

		virtual T eval(float x, float y, float z) {
			printf("The volumebase eval method called!\n");
			return 0.0f;
		};

		T operator()(float x, float y, float z) {
			return eval(x, y, z);
		}

		T eval(Vector x) {
			return eval(x.X(), x.Y(), x.Z());
		}
	};

}