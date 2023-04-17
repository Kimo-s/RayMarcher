#pragma once
#include "Vector.h"
#include "Matrix.h"

using namespace ifs;

namespace ifs {


	// The base abstract volume class
	template <class T>
	class VolumeBase {
	public:

		virtual ~VolumeBase() {}

		virtual T eval(float x, float y, float z) {
			printf("The volumebase eval method called!\n");
			T value{};
			value = value - value;
			return value;
		};

		T operator()(float x, float y, float z) {
			return eval(x, y, z);
		}

		T eval(Vector x) {
			return eval(x.X(), x.Y(), x.Z());
		}


	};

}