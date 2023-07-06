#pragma once
#include "IFSclasses.h"
#include "Color.h"
#include "VolumeBase.h"
// #include "VectorFieldFuncs.h"
#include <iostream>
#define sq(A) (A*A)
using namespace std;


namespace ifs {

	struct FSPNParms {
		float gamma = 1.0f;
		int N = 3;
		float roughness = 1.1f;
		float f = 2.0f;
		float fj = 1.0f;
		Vector xt = Vector(0.0, 0.0, 0.0);
		float A = 1.0f;
	};

	class FDGradHandler
	{
	public:

		FDGradHandler();

		~FDGradHandler() {}

		void setNbTerms(const int n)
		{
			if (n <= (int)_grad_coefficients.size() && n > 0)
			{
				_nb_grad = n;
			}
		}

		const int nbTerms() const { return _nb_grad; }

		void setStep(const double dx)
		{
			if (dx > 0) { _dx = dx; _dy = dx; _dz = dx; }
		}

		void setStep(const double dx, const double dy, const double dz)
		{
			if (dx > 0 && dy > 0 && dz > 0) { _dx = dx; _dy = dy; _dz = dz; }
		}

		const double step() const { return _dx; }
		const double step_x() const { return _dx; }
		const double step_y() const { return _dy; }
		const double step_z() const { return _dz; }

		const double coefficient(const int n) const
		{
			if (n > 0 && n <= _nb_grad) { return _grad_coefficients[_nb_grad - 1][n - 1]; }
			return 0.0;
		}


	private:
		int _nb_grad;
		double _dx, _dy, _dz;
		std::vector< std::vector<double> > _grad_coefficients;
	};

	//-----------------------------------------------------------------------------
	// Setting up logic to be able to determine the data type of the gradient 
	template <typename U>
	struct GradType
	{
		typedef int GType;
	};

	template<>
	struct GradType<float>
	{
		typedef Vector GType;
	};

	template<>
	struct GradType<Vector>
	{
		typedef Matrix GType;
	};

	template<typename U, typename G>
	const G FDGradient(const U& x, const U& y, const U& z)
	{
		return 0;
	};

	template<>
	const Vector FDGradient(const float& x, const float& y, const float& z);

	template<>
	const Matrix FDGradient(const Vector& x, const Vector& y, const Vector& z);


	template <typename T>
	class FieldBase {
	public:

		typedef T volumeDataType;
		typedef typename GradType<T>::GType volumeGradType;
		FDGradHandler gradParams;

		FieldBase() {}
		virtual ~FieldBase() {}

		virtual const volumeDataType eval(const Vector& pos) const = 0;

		virtual const volumeGradType grad(const Vector& P) const
		{
			volumeDataType valueX{}, valueY{}, valueZ{};
			valueX = valueX - valueX;
			valueY = valueY - valueY;
			valueZ = valueZ - valueZ;
			Vector dx(gradParams.step_x(), 0, 0), dy(0, gradParams.step_y(), 0), dz(0, 0, gradParams.step_z());
			for (size_t i = 1; i <= (size_t)gradParams.nbTerms(); i++)
			{
				double coeff = gradParams.coefficient((int)i);
				if (isnan(coeff)) { std::cout << "Volume grad NAN " << i << "   dx " << dx.X() << " " << dx.Y() << " " << dx.Z() << std::endl; }
				valueX += (eval(P + i * dx) - eval(P - i * dx)) * coeff / gradParams.step_x();
				valueY += (eval(P + i * dy) - eval(P - i * dy)) * coeff / gradParams.step_y();
				valueZ += (eval(P + i * dz) - eval(P - i * dz)) * coeff / gradParams.step_z();
			}
			return FDGradient<volumeDataType, volumeGradType>(valueX, valueY, valueZ);
		}


		volumeDataType eval(float x, float y, float z) const {
			return eval(Vector(x, y, z));
		}


		void setFDSize(int nb) { gradParams.setNbTerms(nb); }
		void setFDStep(double dx) { gradParams.setStep(dx); }
		void setFDStep(double dx, double dy, double dz) { gradParams.setStep(dx, dy, dz); }

	};

	class VectorField : public std::shared_ptr<FieldBase<Vector> > {
	public:

		//const VolumeBase<float> volume;


		VectorField() : std::shared_ptr<FieldBase<Vector> >() {
		}
		VectorField(FieldBase<Vector>* f) : std::shared_ptr<FieldBase<Vector> >(f) {
		};
		~VectorField() {};


		static VectorField add(VectorField& v1, const VectorField& v2);
		VectorField max(const VectorField& v2);
		VectorField cut(const VectorField& v2);
		VectorField translate(const Vector transVec);
		VectorField scale(float scaleFactor);
		VectorField mask();
		VectorField rotate(Vector u, float theta);
		VectorField operator+(const VectorField& second);
		VectorField operator-(const VectorField& second);
		// friend const VectorField operator*(const float constant, const VectorField& v){
		// 	return VectorField(new MultiVectorField(v, constant));
		// };
		// friend const VectorField operator*(const VectorField& f, const scalarFieldT& v){
		// 	return VectorField(new VectorTimesScalarField(f, v));
		// };
		VectorField operator*(const scalarFieldT& constant);
		VectorField operator*(const float constant);
	};

	class scalarFieldT : public std::shared_ptr<FieldBase<float> > {
	public:

		//const VolumeBase<float> volume;


		scalarFieldT(): std::shared_ptr<FieldBase<float> >() {
		}
		scalarFieldT(FieldBase<float>* f): std::shared_ptr<FieldBase<float> >(f) {
		};
		~scalarFieldT() {};

		static scalarFieldT add(scalarFieldT& v1, const scalarFieldT& v2);
		scalarFieldT max(const scalarFieldT& v2);
		scalarFieldT cut(const scalarFieldT& v2);
		scalarFieldT translate(const Vector transVec);
		scalarFieldT scale(float scaleFactor);
		scalarFieldT pyroclasticNoise(FSPNParms params);
		scalarFieldT mask();
		scalarFieldT rotate(Vector u, float theta);
		scalarFieldT addGuideParticale(Vector u, FSPNParms params, float fade, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz);
		scalarFieldT operator+(const scalarFieldT& second);
		scalarFieldT operator-(const scalarFieldT& second);
		scalarFieldT operator*(const scalarFieldT& second);
		scalarFieldT operator*(const float num);
		scalarFieldT addWispParticale(VolumeParms parms, FSPNParms fspn1, FSPNParms fspn2, Vector guidPos, float density, float pscale, float wisp_displacement, float clump, int wispCount);
	};

	class ColorField : public std::shared_ptr<FieldBase<Color> > {
	public:

		~ColorField() {
		}

		ColorField() {
		}

		ColorField(FieldBase<Color>* f) : std::shared_ptr<FieldBase<Color> >(f) {
		}


		ColorField add(ColorField& v1, const ColorField& v2);
		ColorField translate(const Vector transVec);
		//ColorField max(const ColorField& second);
		ColorField scale(float scaleFactor);
		ColorField rotate(Vector u, float theta);
		ColorField operator+(const ColorField& second);
	};



	scalarFieldT constantScalarField(float c);
	VectorField constantVectorField(Vector c);
	scalarFieldT funcField(float(*func)(float x, float y, float z));
	VectorField funcField(Vector(*func)(float x, float y, float z));

	ColorField gridColorField(ColorField& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos); 
	ColorField colorMaskField(const scalarFieldT& f, Color col);
	scalarFieldT gridField(int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);

	scalarFieldT plane(Vector normal, Vector center, FSPNParms params);
	scalarFieldT gridField(scalarFieldT& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
	scalarFieldT gridField(VolumeGrid<float> &grid);
	VectorField gridField(VectorField& vol, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
	scalarFieldT gridField(scalarFieldT& vol, Vector lightPos, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);
	scalarFieldT gridField(const char* filename, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz, Vector startPos);

	scalarFieldT GuideParticale(Vector u, FSPNParms params, float fade, int Nx, int Ny, int Nz, float deltax, float deltay, float deltaz);

	scalarFieldT constantField(float constant); //duplicated???
	VectorField warp(VectorField soruce, VectorField wraper);
	VectorField identityVectorField();
	VectorField gradField(scalarFieldT f);
	VectorField noiseVectorField(FSPNParms parms);
	scalarFieldT warp(scalarFieldT soruce, VectorField wraper);

	VectorField incompress(VectorField V, scalarFieldT obs, VolumeParms* parms);
	scalarFieldT advect(scalarFieldT density, VectorField vel, float dt);
	VectorField advect(VectorField s, VectorField vel, float dt);
	float divergence(VectorField V, int i, int j, int k, float dt, Vector startPos);
}