//#include <math.h>
//#include <string>
//
//#define sq(A) (A*A)
//
//float sphere(float x, float y, float z) {
//	return 0.5f - sqrt(x * x + y * y + z * z);
//}
//
//float steinerPatch(float x, float y, float z) {
//	return -(sq(x) * sq(y) + sq(x) * sq(z) + sq(y) * sq(z) - x * y * z);
//}
//
//float ellipse(float x, float y, float z0) {
//	Vector P(x, y, z0);
//	Vector n(0, 0.0f, 1.0f);
//
//	float z = P * n;
//	Vector xbar = P - z * n;
//
//	float rmajor = 0.3f;
//	float rminor = 0.8f;
//
//	return 1.0f - (sq(z) / sq(rmajor)) - (xbar * xbar / sq(rminor));
//}
//
//float ellipse2(float x, float y, float z0) {
//	Vector P(x, y, z0);
//	Vector n(0, 0.0f, 1.0f);
//
//	float z = P * n;
//	Vector xbar = P - z * n;
//
//	float rmajor = 0.8f;
//	float rminor = 0.3f;
//
//	return 1.0f - (sq(z) / sq(rmajor)) - (xbar * xbar / sq(rminor));
//}
//
//float box(float x, float y, float z0) {
//
//	float q = 2;
//	float R = 1.0f;
//
//	return pow(R, 2 * q) - pow(x, 2 * q) - pow(y, 2 * q) - pow(z0, 2 * q);
//}
//
//float icosahedron(float x, float y, float z) {
//	Vector P(x, y, z);
//
//	float T = 2.0f;
//	float res;
//	if (P.magnitude() <= 1.8f * M_PI) {
//		res = cos(x + T * y) + cos(x - T * y) + cos(y + T * z) + cos(y - T * z) + cos(z - T * x) + cos(z + T * x) - 2.0f;
//	}
//	else {
//		res = -1.8 * M_PI;
//	}
//
//	return res;
//}
//
//float cone(float x, float y, float z0) {
//	Vector P(x, y, z0);
//
//	Vector n(0, 0.0f, 1.0f);
//	float h = 1.0f;
//	float theta = 0.4f;
//
//	float res;
//
//	if (P * n < 0) {
//		res = P * n;
//	}
//	else if (P * n > h) {
//		res = h - P * n;
//	}
//	else {
//		res = P * n - P.magnitude() * cos(theta);
//	}
//
//	return res;
//}
//
//float tours(float x, float y, float z) {
//	Vector P(x, y, z);
//	Vector n(0.0f, 1.0f, 0.0f);
//	Vector xbar = P - (P * n) * n;
//
//	float Rmajor = 1.0f;
//	float Rminor = 0.3f;
//
//
//	return 4 * Rmajor * Rmajor * (xbar * xbar) - ((P * P + Rmajor * Rmajor - Rminor * Rminor) * (P * P + Rmajor * Rmajor - Rminor * Rminor));
//}