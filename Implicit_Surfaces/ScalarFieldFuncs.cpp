#include "ScalarFieldFuncs.h"
#include "FieldClasses.h"
#include "VolumeClasses.h"
#include "Vector.h"

using namespace ifs;

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