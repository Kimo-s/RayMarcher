#pragma once
#include "ColorFieldFuncs.h"
#include "FieldClasses.h"

AddColorFields::AddColorFields(const ColorField& a, const ColorField& b) : e1(a), e2(b) {

}

const Color AddColorFields::eval(const Vector& pos) const {
	return e1->eval(pos) + e2->eval(pos);
}

ScaleColorField::ScaleColorField(const ColorField& a, const float& scalefactor) : scalefactor(scalefactor), e2(a) {
}

const Color ScaleColorField::eval(const Vector& pos) const {
	return e2->eval(pos / scalefactor);
}

TranslateColorField::TranslateColorField(const ColorField& a, const Vector& v) : v(v), e2(a) {
}

const Color TranslateColorField::eval(const Vector& pos) const {
	return e2->eval(pos - v);
}

RotateColorField::RotateColorField(const ColorField& a, const Vector& u, float theta): a(a){
	this->u = u.unitvector();
	this->theta = theta;
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

const Color RotateColorField::eval(const Vector& pos) const {
	float xc = m11 * pos.X() + m12 * pos.Y() + m13 * pos.Z();
	float yc = m21 * pos.X() + m22 * pos.Y() + m23 * pos.Z();
	float zc = m31 * pos.X() + m32 * pos.Y() + m33 * pos.Z();
	return a->eval(Vector(xc, yc, zc));
}

ifs::ColorMaskField::ColorMaskField(const scalarFieldT& f, Color col): f(f), col(col)
{
}

const Color ifs::ColorMaskField::eval(const Vector& pos) const
{
	if (f->eval(pos) > 0.0f) {
		return col;
	}
	else {
		return Color(0.0f, 0.0f, 0.0f, 0.0f);
	}
}

ifs::GridColorField::GridColorField(ColorField& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	grid = new VolumeColorGrid(vol, Nx, Ny, Nz, deltax, deltay, deltaz, startPos);
}

const Color ifs::GridColorField::eval(const Vector& pos) const
{
	return grid->eval(pos);
}
