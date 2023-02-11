#include "FieldClasses.h"

using namespace std;
using namespace ifs;

FieldBase<Color>* FieldBase<Color>::add(FieldBase* second) {
	AddColor<Color> toRet = AddColor<Color>(this, second);
	return &toRet;
}

FieldBase<Color>* FieldBase<Color>::operator+(FieldBase* second) {
	AddColor<Color> toRet = AddColor<Color>(this, second);
	return &toRet;
}

FieldBase<Color>* FieldBase<Color>::operator-(FieldBase* second) {
	SubBase<Color> toRet = SubBase<Color>(this, second);
	return &toRet;
}

// Constructers
//ScalarField::ScalarField()
//{
//	volume = NULL;
//	g = NULL;
//	f = NULL;
//	operation = CSGop::None;
//	translateVector = { 0.0f,0.0f,0.0f };
//}

ScalarField::ScalarField(VolumeBase<float>* vol, Vector col)
{
	volume = vol;
	g = NULL;
	f = NULL;
	operation = CSGop::None;
	color = col;
	translateVector = { 0.0f,0.0f,0.0f };
}

ScalarField::ScalarField(float(*func)(float x, float y, float z))
{
	volume = new VolumeFunc(func);
	g = NULL;
	f = NULL;
	operation = CSGop::None;
	translateVector = { 0.0f,0.0f,0.0f };
}

ScalarField::ScalarField(float(*func)(float x, float y, float z), Vector col)
{
	volume = new VolumeFunc(func);
	g = NULL;
	f = NULL;
	operation = CSGop::None;
	color = col;
	translateVector = { 0.0f,0.0f,0.0f };
}

ScalarField::ScalarField(const ScalarField& field)
{
	volume = field.volume;
	g = field.g;
	f = field.f;
	operation = field.operation;
	translateVector = field.translateVector;
	scaleFactor = field.scaleFactor;
	color = field.color;
}

ScalarField::~ScalarField()
{
}