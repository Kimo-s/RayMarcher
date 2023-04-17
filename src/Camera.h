
#ifndef __CAMERA_H__
#define __CAMERA_H__


#include "Vector.h"

namespace ifs
{
class Camera
{

  public:

    Camera();
   ~Camera(){}

    void setEyeViewUp( const ifs::Vector& eye, const ifs::Vector& view, const ifs::Vector& up );
    const ifs::Vector& eye() const  { return position; }
    const ifs::Vector& view() const { return axis_view; }
    const ifs::Vector& up() const   { return axis_up; }
    const ifs::Vector& right() const { return axis_right; }

    // view direction of a pixel at the fractional position x,y.
    // Nominally 0 <= x <= 1 and 0 <= y <= 1 for the primary fov,
    // but the values can extend beyond that
    virtual const ifs::Vector view( const double x, const double y ) const;
    virtual void XY( const ifs::Vector& P, double& x, double& y ) const;
    virtual void XYZ( const ifs::Vector& P, double& x, double& y, double& z ) const;

    void setFov( const double fov );
    const double& fov() const { return FOV; }

    void setAspectRatio( const double ar );
    const double& aspectRatio() const { return aspect_ratio; }

    void setNearPlane( const double n ){ near = n; }
    const double& nearPlane() const { return near; }

    void setFarPlane( const double n ){ far = n; }
    const double& farPlane() const { return far; }

    bool isVisible( const ifs::Vector& P ) const;

    char *__str__(); 

  protected:


    
    double FOV, aspect_ratio;
    double htanfov, vtanfov;
    double near, far;

    ifs::Vector position;
    ifs::Vector axis_right, axis_up, axis_view;



};
}
#endif


