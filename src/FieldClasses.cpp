#include "FieldClasses.h"
#include "ScalarFieldFuncs.h"
#include "ColorFieldFuncs.h"
#include "VectorFieldFuncs.h"
#include "VolumeClasses.h"
#include <memory>

using namespace std;
using namespace ifs;

const std::vector< std::vector<double> > FDGradCoefficientGenerator()
{
	std::vector<double> coeffs;
	std::vector< std::vector<double> > grad_coefficients;
	// n = 1
	coeffs.push_back(0.5);
	grad_coefficients.push_back(coeffs);
	coeffs.clear();

	// n = 2
	coeffs.push_back(6.666666666666666296592e-01);
	coeffs.push_back(-8.333333333333332870740e-02);
	grad_coefficients.push_back(coeffs);
	coeffs.clear();

	// n = 3
	coeffs.push_back(7.500000000000002220446e-01);
	coeffs.push_back(-1.500000000000001332268e-01);
	coeffs.push_back(1.666666666666666296592e-02);
	grad_coefficients.push_back(coeffs);
	coeffs.clear();

	// n = 4
	coeffs.push_back(7.999999999999992672528e-01);
	coeffs.push_back(-2.000000000000004829470e-01);
	coeffs.push_back(3.809523809523807785782e-02);
	coeffs.push_back(-3.571428571428579123309e-03);
	grad_coefficients.push_back(coeffs);
	coeffs.clear();

	// n = 5
	coeffs.push_back(8.333333333333301506940e-01);
	coeffs.push_back(-2.380952380952371660872e-01);
	coeffs.push_back(5.952380952380949968861e-02);
	coeffs.push_back(-9.920634920635062331540e-03);
	coeffs.push_back(7.936507936507907227594e-04);
	grad_coefficients.push_back(coeffs);
	coeffs.clear();

	// n = 6
	coeffs.push_back(8.571428571428428844214e-01);
	coeffs.push_back(-2.678571428571295820475e-01);
	coeffs.push_back(7.936507936506739802063e-02);
	coeffs.push_back(-1.785714285714170776465e-02);
	coeffs.push_back(2.597402597401952516892e-03);
	coeffs.push_back(-1.803751803751641390964e-04);
	grad_coefficients.push_back(coeffs);
	coeffs.clear();

	// n = 7
	coeffs.push_back(8.749999999997666311202e-01);
	coeffs.push_back(-2.916666666664030627132e-01);
	coeffs.push_back(9.722222222201161445643e-02);
	coeffs.push_back(-2.651515151521397634093e-02);
	coeffs.push_back(5.303030303006155132817e-03);
	coeffs.push_back(-6.798756798778969601127e-04);
	coeffs.push_back(4.162504162505081361893e-05);
	grad_coefficients.push_back(coeffs);
	coeffs.clear();

	// n = 8
	coeffs.push_back(8.888888888890622563821e-01);
	coeffs.push_back(-3.111111111113007976492e-01);
	coeffs.push_back(1.131313131315020842349e-01);
	coeffs.push_back(-3.535353535356308696258e-02);
	coeffs.push_back(8.702408702345020702351e-03);
	coeffs.push_back(-1.554001554001562023302e-03);
	coeffs.push_back(1.776001776006136932424e-04);
	coeffs.push_back(-9.712509712526897105722e-06);
	grad_coefficients.push_back(coeffs);
	coeffs.clear();

	// n = 9
	coeffs.push_back(8.999999999973695707922e-01);
	coeffs.push_back(-3.272727272682948718163e-01);
	coeffs.push_back(1.272727272692814326494e-01);
	coeffs.push_back(-4.405594405380468259192e-02);
	coeffs.push_back(1.258741258753500423528e-02);
	coeffs.push_back(-2.797202796916069093835e-03);
	coeffs.push_back(4.495504494691525475096e-04);
	coeffs.push_back(-4.627725215224719976801e-05);
	coeffs.push_back(2.285296402912368231850e-06);
	grad_coefficients.push_back(coeffs);
	coeffs.clear();

	// n = 10
	coeffs.push_back(9.090909090725343144612e-01);
	coeffs.push_back(-3.409090908829451871398e-01);
	coeffs.push_back(1.398601398464504319552e-01);
	coeffs.push_back(-5.244755243406638844927e-02);
	coeffs.push_back(1.678321677751935456224e-02);
	coeffs.push_back(-4.370629367345699872738e-03);
	coeffs.push_back(8.814714692904230272999e-04);
	coeffs.push_back(-1.285479226007764982226e-04);
	coeffs.push_back(1.202787579588283707100e-05);
	coeffs.push_back(-5.412544109527078402898e-07);
	grad_coefficients.push_back(coeffs);
	coeffs.clear();

	// n = 11
	coeffs.push_back(9.166666665040884565130e-01);
	coeffs.push_back(-3.525641023000493645689e-01);
	coeffs.push_back(1.510989006669366252478e-01);
	coeffs.push_back(-6.043956026669643211147e-02);
	coeffs.push_back(2.115384606316183732644e-02);
	coeffs.push_back(-6.221719517138330109163e-03);
	coeffs.push_back(1.481361764454423206663e-03);
	coeffs.push_back(-2.728824305737073381908e-04);
	coeffs.push_back(3.638432414571888845901e-05);
	coeffs.push_back(-3.118656339922096169245e-06);
	coeffs.push_back(1.288700981451494452140e-07);
	grad_coefficients.push_back(coeffs);
	coeffs.clear();

	// n = 12
	coeffs.push_back(9.230769225029789026848e-01);
	coeffs.push_back(-3.626373617033492591233e-01);
	coeffs.push_back(1.611721583806974000819e-01);
	coeffs.push_back(-6.799450484120614368599e-02);
	coeffs.push_back(2.559793110528202006448e-02);
	coeffs.push_back(-8.295625777082871188384e-03);
	coeffs.push_back(2.245432612438852566783e-03);
	coeffs.push_back(-4.911883622080408370522e-04);
	coeffs.push_back(8.316416682614131267431e-05);
	coeffs.push_back(-1.020651206170876503291e-05);
	coeffs.push_back(8.068388089184792311508e-07);
	coeffs.push_back(-3.081676248363273868301e-08);
	grad_coefficients.push_back(coeffs);
	coeffs.clear();

	// n = 13
	coeffs.push_back(9.285714267400464461133e-01);
	coeffs.push_back(-3.714285685139248061049e-01);
	coeffs.push_back(1.702380751958524618406e-01);
	coeffs.push_back(-7.510503985705131724249e-02);
	coeffs.push_back(3.004201550432860148843e-02);
	coeffs.push_back(-1.054105789791678036982e-02);
	coeffs.push_back(3.162317968689574727154e-03);
	coeffs.push_back(-7.905793074257984895739e-04);
	coeffs.push_back(1.597129844416282369157e-04);
	coeffs.push_back(-2.499855371995233626345e-05);
	coeffs.push_back(2.840742476654620619051e-06);
	coeffs.push_back(-2.083212542482251737132e-07);
	coeffs.push_back(7.396022670499950410872e-09);
	grad_coefficients.push_back(coeffs);
	coeffs.clear();

	return grad_coefficients;
}

FDGradHandler::FDGradHandler() :
	_nb_grad(1),
	_dx(0.01),
	_dy(0.01),
	_dz(0.01)
{
	_grad_coefficients = FDGradCoefficientGenerator();
}

template<>
const Vector ifs::FDGradient(const float& x, const float& y, const float& z)
{
	float xx = x;
	float yy = y;
	float zz = z;
	if (isnan(x)) { cout << "Gradient NAN x\n"; xx = 0.0; }
	if (isnan(y)) { cout << "Gradient NAN y\n"; yy = 0.0; }
	if (isnan(z)) { cout << "Gradient NAN z\n"; zz = 0.0; }
	return Vector(xx, yy, zz);
}

template<>
const Matrix ifs::FDGradient(const Vector& x, const Vector& y, const Vector& z)
{
	Matrix m(x, y, z);
	return m.transpose();
}

scalarFieldT ifs::constantScalarField(float c)
{
	return scalarFieldT(new constScalarField(c));
}

VectorField ifs::constantVectorField(Vector c)
{
	return VectorField(new ConstVectorField(c));
}

scalarFieldT ifs::funcField(float(*func)(float x, float y, float z))
{
	return scalarFieldT(new FuncScalarField(func));
}

VectorField ifs::funcField(Vector(*func)(float x, float y, float z))
{
	return VectorField(new FuncVectorField(func));
}

ColorField ifs::gridColorField(ColorField& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	return ColorField(new GridColorField(vol, Nx, Ny, Nz, deltax, deltay, deltaz, startPos));
}

ColorField ifs::colorMaskField(const scalarFieldT& f, Color col)
{
	return ColorField(new ColorMaskField(f, col));
}

scalarFieldT ifs::gridField(VolumeGrid<float> &grid){
	return scalarFieldT(new GridScalarField(&grid));
}

scalarFieldT ifs::gridField(int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	return scalarFieldT(new GridScalarField(Nx, Ny, Nz, deltax, deltay, deltaz, startPos));
}

scalarFieldT ifs::constantField(float constant)
{
	return scalarFieldT(new constScalarField(constant));
}

scalarFieldT ifs::plane(Vector normal, Vector center, FSPNParms params)
{
	return scalarFieldT(new planeScalarField(normal, center, params));
}

scalarFieldT ifs::gridField(scalarFieldT& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	return scalarFieldT(new GridScalarField(vol, Nx, Ny, Nz, deltax, deltay, deltaz, startPos));
}

VectorField ifs::gridField(VectorField& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	return VectorField(new GridVectorField(vol, Nx, Ny, Nz, deltax, deltay, deltaz, startPos));
}

scalarFieldT ifs::gridField(scalarFieldT& vol, Vector lightPos, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	return scalarFieldT(new GridScalarField(vol, lightPos, Nx, Ny, Nz, deltax, deltay, deltaz, startPos));
}

scalarFieldT ifs::gridField(const char* filename, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos)
{
	return scalarFieldT(new GridScalarField(filename, Nx, Ny, Nz, deltax, deltay, deltaz, startPos));
}

scalarFieldT ifs::GuideParticale(Vector u, FSPNParms params, float fade, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz)
{
	return scalarFieldT(new addGuideParticaleScalarField(u, params, fade, Nx, Ny, Nz, deltax, deltay, deltaz));
}

VectorField ifs::warp(VectorField soruce, VectorField wraper)
{
	return VectorField(new warpVectorField(soruce, wraper));
}

VectorField ifs::identityVectorField()
{
	return VectorField(new IdentityVectorField());
}

scalarFieldT ifs::warp(scalarFieldT soruce, VectorField wraper)
{
	return scalarFieldT(new warpScalarField(soruce, wraper));
}

float ifs::divergence(VectorField V, int i, int j, int k, float dt, Vector startPos) {

	Vector pos1 = Vector((i - 1) * dt + startPos[0], j * dt + startPos[1], k * dt + startPos[2]);
	Vector pos2 = Vector((i + 1) * dt + startPos[0], j * dt + startPos[1], k * dt + startPos[2]);
	float x = (V->eval(pos2).X() - V->eval(pos1).X()) / (2 * dt);

	pos1 = Vector(i * dt + startPos[0], (j-1) * dt + startPos[1], k * dt + startPos[2]);
	pos2 = Vector(i * dt + startPos[0], (j+1) * dt + startPos[1], k * dt + startPos[2]);
	float y = (V->eval(pos2).Y() - V->eval(pos1).Y()) / (2 * dt);

	pos1 = Vector(i * dt + startPos[0], j * dt + startPos[1], (k-1) * dt + startPos.Z());
	pos2 = Vector(i * dt + startPos[0], j * dt + startPos[1], (k+1) * dt + startPos.Z());
	float z = (V->eval(pos2).Z() - V->eval(pos1).Z()) / (2 * dt);
	return x + y + z;
}

VectorField ifs::incompress(VectorField V, VolumeParms* parms)
{
	VolumeGrid<Vector>* newV = new VolumeGrid<Vector>(parms->Nx, parms->Ny, parms->Nz, parms->dx, parms->dx, parms->dx, parms->startPos);
	VolumeGrid<float> p(parms->Nx, parms->Ny, parms->Nz, parms->dx, parms->dx, parms->dx, parms->startPos);
	VolumeGrid<float> oldp(parms->Nx, parms->Ny, parms->Nz, parms->dx, parms->dx, parms->dx, parms->startPos);
	p.defaultValue = 0.0f;
	p.minval = -100.0f;
	oldp.defaultValue = 0.0f;
	oldp.minval= -100000.0f;
	float dt = parms->dx;

	for (int k = 0; k < p.Nz; k++) {
		for (int j = 0; j < p.Ny; j++) {
			for (int i = 0; i < p.Nx; i++) {
				p.set(i, j, k, 0.0f);
				oldp.set(i, j, k, divergence(V, i, j, k, dt, parms->startPos));
			}
		}
	}

	//cout << "Value at the mid point " << divergence(V, 20, 20, 20, dt, parms->startPos) << endl;
	int numIters = 20;
	for (int iter = 0; iter < numIters; iter++) {
		//#pragma omp parallel for schedule(dynamic)
		for (int k = 0; k < p.Nz; k++) {
			for (int j = 0; j < p.Ny; j++) {
				for (int i = 0; i < p.Nx; i++) {

					// float div = oldp->get(i,j,k);
					float div = oldp.get(i,j,k);

					float avg = (p.get(i + 1, j, k) + p.get(i - 1, j, k)
						+ p.get(i, j + 1, k) + p.get(i, j - 1, k)
						+ p.get(i, j, k + 1) + p.get(i, j, k - 1)) / 6.0;
					float newp = avg - ((dt * dt) / 6.0) * div;
					//cout << "lol: " << newp << "\n";
					p.set(i, j, k, newp);
					//cout << "test p  " << p.get(i, j, k) << "\n";
				}
			}
		}
	}

	scalarFieldT pfield = scalarFieldT(new GridScalarField(&p));


	for (int k = 0; k < newV->Nz; k++) {
		for (int j = 0; j < newV->Ny; j++) {
			for (int i = 0; i < newV->Nx; i++) {
				Vector pos = Vector(i * dt + parms->startPos[0], j * dt + parms->startPos[1], k * dt + parms->startPos[2]);
				newV->set(i, j, k, V->eval(pos) - pfield->grad(pos));
			}
		}
	}

	//cout << "Test vector field: " << res->eval(Vector(0, 0, 0)).__str__() << endl;
	//
	//cout << "Value at the mid point before" << divergence(V, 20, 20, 20, dt, parms->startPos) << endl;
	//V = res;
	//cout << "Value at the mid point after " << divergence(V, 20, 20, 20, dt, parms->startPos) << endl;
	//cout << "Value at the mid point after (V) " << divergence(V, 20, 20, 20, dt, parms->startPos) << endl;

	return VectorField(new GridVectorField(newV));
}

scalarFieldT ifs::advect(scalarFieldT density, VectorField vel, float dt)
{
	scalarFieldT p1 = warp(density, identityVectorField() - vel * dt);
	scalarFieldT p2 = warp(p1, identityVectorField() - vel * (-dt));

	return p1 + (density-p2)*0.5;
}

VectorField ifs::advect(VectorField s, VectorField vel, float dt)
{
	return warp(s, identityVectorField() - vel * dt);
}

VectorField ifs::gradField(scalarFieldT a) {
	return VectorField(new gradVectorField(a));
}

VectorField ifs::noiseVectorField(FSPNParms parms)
{

	return VectorField(new NoiseVectorField(parms));
}

scalarFieldT scalarFieldT::add(scalarFieldT& v1, const scalarFieldT& v2) { return scalarFieldT(new AddScalarFields(v1, v2)); }
scalarFieldT ifs::scalarFieldT::max(const scalarFieldT& v2) { return scalarFieldT(new MaxScalarField(*this, v2)); }
scalarFieldT ifs::scalarFieldT::cut(const scalarFieldT& v2) { return scalarFieldT(new CutScalarField(*this, v2)); }
scalarFieldT scalarFieldT::translate(const Vector transVec) { return scalarFieldT(new TranslateScalarField(*this, transVec)); };
scalarFieldT scalarFieldT::scale(float scaleFactor) { return scalarFieldT(new ScaleScalarField(*this, scaleFactor)); }
scalarFieldT ifs::scalarFieldT::pyroclasticNoise(FSPNParms params) { return scalarFieldT( new pyroclasticScalarField(*this, params)); }
scalarFieldT scalarFieldT::mask() { return scalarFieldT(new MaskScalarField(*this)); };
scalarFieldT scalarFieldT::rotate(Vector u, float theta) { return scalarFieldT(new RotateScalarField(*this, u, theta)); }
scalarFieldT ifs::scalarFieldT::addGuideParticale(Vector u, FSPNParms params, float fade, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz)
{
	return this->max(scalarFieldT(new addGuideParticaleScalarField(u, params, fade, Nx, Ny, Nz, deltax, deltay, deltaz)));
}
scalarFieldT ifs::scalarFieldT::addWispParticale(VolumeParms parms, FSPNParms fspn1, FSPNParms fspn2, Vector guidPos, float density, float pscale, float wisp_displacement, float clump, int wispCount)
{
	return this->max(scalarFieldT(new wispScalarField(parms, fspn1, fspn2, guidPos, density, pscale, wisp_displacement, clump, wispCount)));
}
scalarFieldT scalarFieldT::operator+(const scalarFieldT& second) { return add(*this, second); }
scalarFieldT ifs::scalarFieldT::operator-(const scalarFieldT &second)
{
    return add(*this, constantField(-1)*second);
}
scalarFieldT ifs::scalarFieldT::operator*(const scalarFieldT& second)
{
	return scalarFieldT(new MultiplyScalarFields(*this, second));
}
scalarFieldT ifs::scalarFieldT::operator*(const float num)
{
    return scalarFieldT(new MultiplyScalarFields(*this, constantField(num)));
}

VectorField VectorField::add(VectorField& v1, const VectorField& v2) { return VectorField(new AddVectorFields(v1, v2)); }
VectorField ifs::VectorField::max(const VectorField& v2) { return VectorField(new MaxVectorField(*this, v2)); }
VectorField ifs::VectorField::cut(const VectorField& v2) { return VectorField(new CutVectorField(*this, v2)); }
VectorField VectorField::translate(const Vector transVec) { return VectorField(new TranslateVectorField(*this, transVec)); };
VectorField VectorField::scale(float scaleFactor) { return VectorField(new ScaleVectorField(*this, scaleFactor)); };
VectorField VectorField::mask() { return VectorField(new MaskVectorField(*this)); };
VectorField VectorField::rotate(Vector u, float theta) { return VectorField(new RotateVectorField(*this, u, theta)); };
VectorField VectorField::operator+(const VectorField& second) { return add(*this, second); }
VectorField VectorField::operator-(const VectorField& second) { return VectorField(new SubVectorFields(*this, second)); }
VectorField ifs::VectorField::operator*(const float constant) { return VectorField(new MultiVectorField(*this, constant)); }
VectorField ifs::VectorField::operator*(const scalarFieldT& constant) { return VectorField(new VectorTimesScalarField(*this, constant)); }

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

