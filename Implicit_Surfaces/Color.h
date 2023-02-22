#pragma once
#include <iostream>

namespace ifs {

	class Color {
	public:
		float r;
		float g;
		float b;
		float a;

		Color() {
			this->r = 0.0f;
			this->g = 0.0f;
			this->b = 0.0f;
			this->a = 0.0f;
		}

		Color(float r, float g, float b, float a) {
			this->r = r;
			this->g = g;
			this->b = b;
			this->a = a;
		}

		const Color operator+ (const Color& col) const {
			return Color((this->r + col.r), (this->g + col.g), (this->b + col.b), (this->a + col.a));
		}

		const bool operator< (const Color& col) const {
			float magn1 = sqrt(r * r + g * g + b * b + a * a);
			float magn2 = sqrt(col.r * col.r + col.g * col.g + col.b * col.b + col.a * col.a);
			return magn1 < magn2;
		}


		Color operator- (const Color& col) {
			return Color((this->r - col.r), (this->g - col.g), (this->b - col.b),  (this->a - col.a));
		}

		Color operator* (float scalar) {
			return Color(scalar * this->r, scalar * this->g, scalar * this->b, scalar * this->a);
		}

		Color operator* (Color col) {
			return Color(col.r * this->r, col.g * this->g, col.b * this->b, col.a * this->a);
		}

		Color& operator+= (const Color& v)
		{
			r += v.r; b += v.b; g += v.g; a += v.a; return *this;
		}

		const bool operator== (const Color& v) const
		{
			if (v.r == r && v.b == b && v.g == g && v.a == a) {
				return true;
			}
			else {
				return false;
			}
		}

		friend const Color operator* (const double scalar, const Color& v) {
			return v * scalar;
		}

		const Color operator*  (const double scalar) const {
			return Color(scalar * r, scalar * g, scalar * b, scalar * a);
		}

		const Color operator/  (const double scalar) const {
			return Color(r / scalar, g / scalar, b / scalar, a / scalar);
		}

		friend const Color operator/ (const double scalar, const Color& v) {
			return v / scalar;
		}

		char* print() {
			char* buf = new char[100];
			sprintf_s(buf, 100, "(R:%.2f,G:%.2f,B:%.2f,A:%.2f)", r, g, b, a);
			return buf;
		}
	};

}