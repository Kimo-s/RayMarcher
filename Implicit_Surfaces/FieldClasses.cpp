#include "FieldClasses.h"
#include "ScalarFieldFuncs.h"
#include "ColorFieldFuncs.h"
#include <memory>

using namespace std;
using namespace ifs;


//template<class T>
//FieldBase<T>* FieldBase<T>::add(FieldBase<T>& second) {
//	AddFields<T> toRet = AddFields<T>(this, second);
//	return &toRet;
//}

//FieldBase<float>* FieldBase<float>::add(FieldBase<float>& second) {
//	return new AddFields<float>(this, &second);
//}
//
//FieldBase<Color>* FieldBase<Color>::add(FieldBase<Color>& second) {
//	return new AddFields<Color>(this, &second);
//}

//template<class T>
//FieldBase<T>* FieldBase<T>::translate(Vector transVec) {
//	return new TranslateField<T>(this, transVec);
//}

//FieldBase<float>* FieldBase<float>::translate(Vector transVec)
//{
//	return new TranslateField<float>(this, transVec);
//}
//
//FieldBase<Color>* FieldBase<Color>::translate(Vector transVec)
//{
//	return new TranslateField<Color>(this, transVec);
//}
//
//FieldBase<float>* FieldBase<float>::rotate(Vector u, float theta)
//{
//	return new RotateField<float>(this, u, theta);
//}
//
//FieldBase<Color>* FieldBase<Color>::rotate(Vector u, float theta)
//{
//	return new RotateField<Color>(this, u, theta);
//}
//
//FieldBase<float>* FieldBase<float>::scale(float scalefactor)
//{
//	return new ScaleField<float>(this, scalefactor);
//}
//
//FieldBase<Color>* FieldBase<Color>::scale(float scalefactor)
//{
//	return new ScaleField<Color>(this, scalefactor);
//}
//
//template<class T>
//FieldBase<T>* FieldBase<T>::scale(float scaleFactor)
//{
//	ScaleField<T> toRet = ScaleField<T>(this, scaleFactor);
//	return &toRet;
//}
//
//template<class T>
//FieldBase<T>* FieldBase<T>::operator+(FieldBase* second) {
//	AddFields<T> toRet = AddFields<T>(this, second);
//	return &toRet;
//}
//
//
//template<class T>
//FieldBase<T>* FieldBase<T>::operator-(FieldBase* second) {
//	SubBase<T> toRet = SubBase<T>(this, second);
//	return &toRet;
//}


//scalarFieldT::scalarFieldT(VolumeBase<float>& vol)
//{
//	volume = &vol;
//}
//
scalarFieldT ifs::funcField(float(*func)(float x, float y, float z))
{
	return scalarFieldT(new FuncScalarField(func));
}

ColorField ifs::gridColorField(ColorField& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	return ColorField(new GridColorField(vol, Nx, Ny, Nz, deltax, deltay, deltaz, startPos));
}

ColorField ifs::colorMaskField(const scalarFieldT& f, Color col)
{
	return ColorField(new ColorMaskField(f, col));
}

scalarFieldT ifs::gridField(int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	return scalarFieldT(new GridScalarField(Nx, Ny, Nz, deltax, deltay, deltaz, startPos));
}

scalarFieldT ifs::gridField(scalarFieldT& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	return scalarFieldT(new GridScalarField(vol, Nx, Ny, Nz, deltax, deltay, deltaz, startPos));
}

scalarFieldT ifs::gridField(scalarFieldT& vol, Vector lightPos, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	return scalarFieldT(new GridScalarField(vol, lightPos, Nx, Ny, Nz, deltax, deltay, deltaz, startPos));
}

scalarFieldT ifs::gridField(const char* filename, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	return scalarFieldT(new GridScalarField(filename, Nx, Ny, Nz, deltax, deltay, deltaz, startPos));
}

scalarFieldT scalarFieldT::add(scalarFieldT& v1, const scalarFieldT& v2) { return scalarFieldT(new AddScalarFields(v1, v2)); }
scalarFieldT ifs::scalarFieldT::max(const scalarFieldT& v2) { return scalarFieldT(new MaxScalarField(*this, v2)); };
scalarFieldT scalarFieldT::translate(const Vector transVec) { return scalarFieldT(new TranslateScalarField(*this, transVec)); };
scalarFieldT scalarFieldT::scale(float scaleFactor) { return scalarFieldT(new ScaleScalarField(*this, scaleFactor)); };
scalarFieldT scalarFieldT::mask() { return scalarFieldT(new MaskScalarField(*this)); };
scalarFieldT scalarFieldT::rotate(Vector u, float theta) { return scalarFieldT(new RotateScalarField(*this, u, theta)); };
scalarFieldT scalarFieldT::operator+(const scalarFieldT& second) { return add(*this, second); };

ColorField ColorField::add(ColorField& v1, const ColorField& v2) { return ColorField(new  AddColorFields(v1, v2)); };
ColorField ColorField::translate(const Vector transVec) { return ColorField(new TranslateColorField(*this, transVec)); };
ColorField ColorField::scale(float scaleFactor) { return ColorField(new ScaleColorField(*this, scaleFactor)); };
ColorField ColorField::rotate(Vector u, float theta) { return ColorField(new RotateColorField(*this, u, theta)); };
ColorField ColorField::operator+(const ColorField& second) { return add(*this, second); };

// Trash delete soon
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

