#include "ScalarFieldFuncs.h"
#include "FieldClasses.h"
#include "VolumeClasses.h"
#include "Vector.h"
#include <cstdlib>
#include <cassert>
#include <chrono>
#include <random>
#include <omp.h>

#define DB_PERLIN_IMPL
#include "db_perlin.hpp"


using namespace ifs;

float evalFSPN(FSPNParms parm, Vector pos) {
	Vector temp = (pos - parm.xt) * parm.f;
	float toAdd = db::perlin(temp.X(), temp.Y(), temp.Z());

	for (int q = 1; q <= parm.N; q++) {
		temp = (pos - parm.xt) * parm.f * pow(parm.fj, q);
		toAdd += pow(parm.roughness, q) * db::perlin(temp.X(), temp.Y(), temp.Z());
	}
	toAdd *= (1 - parm.roughness) / (1 - pow(parm.roughness, parm.N));

	return pow(toAdd, parm.gamma) * parm.A;
}

AddScalarFields::AddScalarFields(const scalarFieldT& a, const scalarFieldT& b): e1(a), e2(b) {

}

const float AddScalarFields::eval(const Vector& pos) const
{
	return e1->eval(pos) + e2->eval(pos);
}

ScaleScalarField::ScaleScalarField(const scalarFieldT& a1, const float& scalefactor1): scalefactor(scalefactor1), e2(a1){
}

const float ScaleScalarField::eval(const Vector& pos) const {
	return e2->eval(pos / this->scalefactor);
}

TranslateScalarField::TranslateScalarField(const scalarFieldT& a, const Vector& v) : v(v), e2(a) {
}

const float TranslateScalarField::eval(const Vector& pos) const {
	return e2->eval(pos - v);
}

RotateScalarField::RotateScalarField(const scalarFieldT& a, const Vector& u, float theta) : a(a), u(u.unitvector()), theta(theta) {
	m11 = (cos(theta) + sq(u.X()) * (1.0f - cos(theta)));
	m12 = ((u.X() * u.Y() * (1.0f - cos(theta)) - u.Z() * sin(theta)));
	m13 = ((u.X() * u.Z() * (1.0f - cos(theta)) + u.Y() * sin(theta)));

	m21 = ((u.X() * u.Y() * (1.0f - cos(theta)) + u.Z() * sin(theta)));
	m22 = (cos(theta) + sq(u.Y()) * (1.0f - cos(theta)));
	m23 = ((u.Y() * u.Z() * (1.0f - cos(theta)) - u.X() * sin(theta)));

	m31 = ((u.Z() * u.X() * (1.0f - cos(theta)) - u.Y() * sin(theta)));
	m32 = ((u.Z() * u.Y() * (1.0f - cos(theta)) + u.X() * sin(theta)));
	m33 = (cos(theta) + sq(u.Z()) * (1.0f - cos(theta)));
}

const float RotateScalarField::eval(const Vector& pos) const {
	float xc = m11 * pos.X() + m12 * pos.Y() + m13 * pos.Z();
	float yc = m21 * pos.X() + m22 * pos.Y() + m23 * pos.Z();
	float zc = m31 * pos.X() + m32 * pos.Y() + m33 * pos.Z();
	return a->eval(Vector(xc, yc, zc));
}

ifs::MaxScalarField::MaxScalarField(const scalarFieldT& a, const scalarFieldT& b): a(a), b(b)
{
}

const float ifs::MaxScalarField::eval(const Vector& pos) const
{
	float v1, v2;
	v1 = a->eval(pos);
	v2 = b->eval(pos);
	if (v1 > v2) {
		return v1;
	}
	else {
		return v2;
	}
}

ifs::GridScalarField::GridScalarField(int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	grid = new VolumeGrid<float>(Nx, Ny, Nz, deltax, deltay, deltaz, startPos);
}

GridScalarField::GridScalarField(scalarFieldT& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos) {
	grid = new VolumeGrid<float>(vol, Nx, Ny, Nz, deltax, deltay, deltaz, startPos);
}

ifs::GridScalarField::GridScalarField(scalarFieldT& vol, Vector lightPos, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	grid = new VolumeGrid<float>(vol, lightPos, Nx, Ny, Nz, deltax, deltay, deltaz, startPos);
}

ifs::GridScalarField::GridScalarField(const char* filename, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	grid = new VolumeGrid<float>(filename, Nx, Ny, Nz, deltax, deltay, deltaz, startPos);
}

const float GridScalarField::eval(const Vector& pos) const {
	return grid->eval(pos);
}

ifs::MaskScalarField::MaskScalarField(const scalarFieldT& a): a(a)
{
}


const float MaskScalarField::eval(const Vector& pos) const {
	if (a->eval(pos) > 0.0f) {
		return 1.0f;
	}
	else {
		return 0.0f;
	}
}

ifs::CutScalarField::CutScalarField(const scalarFieldT& a, const scalarFieldT& b) : a(a), b(b)
{
}

const float ifs::CutScalarField::eval(const Vector& pos) const
{
	float v1, v2;
	v1 = a->eval(pos);
	v2 = -b->eval(pos);
	if (v1 > v2) {
		return v1;
	}
	else {
		return v2;
	}
}

ifs::pyroclasticScalarField::pyroclasticScalarField(const scalarFieldT& a, FSPNParms params) :
	params(params),
	a(a)
{
}

const float ifs::pyroclasticScalarField::eval(const Vector& pos) const
{
	float val = a->eval(pos);
	Vector cpt = pos - val * a->grad(pos);

	float toAdd = 0.0f;
	for (int i = 0; i < params.N; i++) {
		Vector temp = (cpt - params.xt) * params.f * pow(params.fj, i);
		toAdd += pow(params.roughness, i) * db::perlin(temp.X(), temp.Y(), temp.Z());
	}
	toAdd *= (1 - params.roughness) / (1 - pow(params.roughness, params.N));

	return val + pow(fabs(toAdd), params.gamma) * params.A;
}

ifs::addGuideParticaleScalarField::addGuideParticaleScalarField(Vector u, FSPNParms params, float fade, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz)
{

	Vector startpos = u - Vector(deltax * Nx / 2.0f, deltay * Ny / 2.0f, deltaz * Nz / 2.0f);
	grid = new VolumeGrid<float>(Nx, Ny, Nz, deltax, deltay, deltaz, startpos);
	float pscale = 1.0f;

	#pragma omp parallel for schedule(dynamic) num_threads(20) collapse(2)
	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
			for (int i = 0; i < Nx; i++) {
				Vector pos(i * deltax + startpos[0], j * deltay + startpos[1], k * deltaz + startpos[2]);
				float falloff = pow((pos - u).magnitude(), fade) / pscale;

				if (1.0f - falloff > 1.0f) {
					falloff = 1.0f;
				} else if (1.0f - falloff < 0.0f) {
					falloff = 0.0f;
				}

				float toAdd = evalFSPN(params, pos);

				grid->set(i, j, k, toAdd * falloff);
			}
		}
	}

}

const float ifs::addGuideParticaleScalarField::eval(const Vector& pos) const
{
	return grid->eval(pos);
}

ifs::wispScalarField::wispScalarField(VolumeParms parms, FSPNParms fspn1, FSPNParms fspn2, Vector guidPos, float density, float pscale, float wisp_displacement, float clump, int wispCount)
{
	Vector startpos = guidPos - Vector(parms.dx * parms.Nx / 2.0f, parms.dy * parms.Ny / 2.0f, parms.dz * parms.Nz / 2.0f);
	grid = new VolumeGrid<float>(parms.Nx, parms.Ny, parms.Nz, parms.dx, parms.dy, parms.dz, startpos);

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(-1.0, 1.0);

	float x = distribution(generator);
	float y = distribution(generator);
	float z = distribution(generator);

	for (int w = 0; w < wispCount; w++) {

		float vx = distribution(generator);
		float vy = distribution(generator);
		float vz = distribution(generator);

		float corr = 0.5f;
		x = corr * x + (1 - corr) * vx;
		y = corr * y + (1 - corr) * vy;
		z = corr * z + (1 - corr) * vz;

		float radius = sqrt(x * x + y * y + z * z);
		float xsphere = x / radius;
		float ysphere = y / radius;
		float zsphere = z / radius;

		float radial_displacement = pow(fabs(evalFSPN(fspn1,Vector(x,y,z))), clump);
		xsphere *= radial_displacement;
		ysphere *= radial_displacement;
		zsphere *= radial_displacement;

		float xdot = guidPos.X() + xsphere * pscale;
		float ydot = guidPos.Y() + ysphere * pscale;
		float zdot = guidPos.Z() + zsphere * pscale;

		float shift = 0.1f;

		float xfsn = evalFSPN(fspn2, Vector(xsphere, ysphere, zsphere)) * wisp_displacement;
		float yfsn = evalFSPN(fspn2, Vector(xsphere + shift, ysphere + shift, zsphere + shift)) * wisp_displacement;
		float zfsn = evalFSPN(fspn2, Vector(xsphere - shift, ysphere - shift, zsphere - shift)) * wisp_displacement;
		xdot += xfsn;
		ydot += yfsn;
		zdot += zfsn;

		//cout << "Value of wisp " << w << ": " << floor((xdot - startpos[0]) / parms.dx) << "," << floor((ydot - startpos[1]) / parms.dy) << "," << floor((zdot - startpos[2]) / parms.dz) << endl;
//		grid->set(static_cast <int>(floor((xdot - startpos[0]) / parms.dx)),
//			static_cast <int>(floor((ydot - startpos[1]) / parms.dy)),
	//		static_cast <int>(floor((zdot - startpos[2]) / parms.dz)), density);

		grid->splat(Vector(xdot, ydot, zdot), density);
	}
}

const float ifs::wispScalarField::eval(const Vector& pos) const
{
	return grid->eval(pos);
}

const float ifs::planeScalarField::eval(const Vector& pos) const
{

	if (fabs(pos.X()) > 2.0) {
		return 0.0f;
	}
	else if (fabs(pos.Y()) > 2.0) {
		return 0.0f;
	}


	float val = (pos - center) * normal;
	Vector cpt = pos - val * normal;

	float noise = evalFSPN(params, cpt);
	if (noise < 0.0f) {
		noise *= 0.2f;
	}

	return val + noise;
}
