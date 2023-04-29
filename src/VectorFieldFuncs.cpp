#include "VectorFieldFuncs.h"
#include "FieldClasses.h"
#include "VolumeClasses.h"
#include "Vector.h"
#include "ScalarFieldFuncs.h"

using namespace ifs;

AddVectorFields::AddVectorFields(const VectorField& a, const VectorField& b) : e1(a), e2(b) {

}

const Vector AddVectorFields::eval(const Vector& pos) const
{
	return e1->eval(pos) + e2->eval(pos);
}

ScaleVectorField::ScaleVectorField(const VectorField& a1, const float& scalefactor1) : scalefactor(scalefactor1), e2(a1) {
}

const Vector ScaleVectorField::eval(const Vector& pos) const {
	return e2->eval(pos / this->scalefactor);
}

TranslateVectorField::TranslateVectorField(const VectorField& a, const Vector& v) : v(v), e2(a) {
}

const Vector TranslateVectorField::eval(const Vector& pos) const {
	return e2->eval(pos - v);
}

RotateVectorField::RotateVectorField(const VectorField& a, const Vector& u, float theta) : a(a), u(u.unitvector()), theta(theta) {
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

const Vector RotateVectorField::eval(const Vector& pos) const {
	float xc = m11 * pos.X() + m12 * pos.Y() + m13 * pos.Z();
	float yc = m21 * pos.X() + m22 * pos.Y() + m23 * pos.Z();
	float zc = m31 * pos.X() + m32 * pos.Y() + m33 * pos.Z();
	return a->eval(Vector(xc, yc, zc));
}

ifs::MaxVectorField::MaxVectorField(const VectorField& a, const VectorField& b) : a(a), b(b)
{
}

const Vector ifs::MaxVectorField::eval(const Vector& pos) const
{
	Vector v1, v2;
	v1 = a->eval(pos);
	v2 = b->eval(pos);
	if (v1 > v2) {
		return v1;
	}
	else {
		return v2;
	}
}

ifs::MaskVectorField::MaskVectorField(const VectorField& a) : a(a)
{
}


const Vector MaskVectorField::eval(const Vector& pos) const {
	if (a->eval(pos).magnitude() > 1.0f) {
		return a->eval(pos).unitvector();
	}
	else {
		return a->eval(pos);
	}
}

ifs::CutVectorField::CutVectorField(const VectorField& a, const VectorField& b) : a(a), b(b)
{
}

const Vector ifs::CutVectorField::eval(const Vector& pos) const
{
	Vector v1, v2;
	v1 = a->eval(pos);
	v2 = -b->eval(pos);
	if (v1 > v2) {
		return v1;
	}
	else {
		return v2;
	}
}

ifs::SubVectorFields::SubVectorFields(const VectorField& a, const VectorField& b) : e1(a), e2(b)
{
}

const Vector ifs::SubVectorFields::eval(const Vector& pos) const
{
	return e1->eval(pos)-e2->eval(pos);
}

const Vector ifs::NoiseVectorField::eval(const Vector& pos) const
{
	float dt = 0.01f;
	float x = evalFSPN(parms, pos + dt * Vector(1.0, 0.0, 0.0));
	float y = evalFSPN(parms, pos + dt * Vector(0.0, 1.0, 0.0));
	float z = evalFSPN(parms, pos + dt * Vector(0.0, 0.0, 1.0));
	return Vector((y - z) / dt, (z - x) / dt, (x - y) / dt);
}

ifs::MultiVectorField::MultiVectorField(const VectorField& a, const float constant) : a(a), constant(constant)
{
	this->constant = constant;
}

const Vector ifs::MultiVectorField::eval(const Vector& pos) const
{
	return a->eval(pos)*constant;
}

ifs::GridVectorField::GridVectorField(VectorField& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	grid = new VolumeGrid<Vector>(vol, Nx, Ny, Nz, deltax, deltay, deltaz, startPos);
}

const Vector ifs::GridVectorField::eval(const Vector& pos) const
{
	return grid->eval(pos);
}
