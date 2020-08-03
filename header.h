#include <math.h>   
#include <stdlib.h> 
#include <stdio.h> 

double Pi=3.14159;
double Pi_inv=1.0/Pi;

struct RGB {
	double r;
	double g;
	double b;
};

struct Vec 
{        
  double x, y, z;                 
  Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; }
	  Vec operator+(const Vec &b) const 
	  {
		  return Vec(x+b.x,y+b.y,z+b.z); 
	  }
	  Vec operator-(const Vec &b) const
	  { 
		  return Vec(x-b.x,y-b.y,z-b.z);
	  }
	  Vec operator*(double b) const 
	  {
		  return Vec(x*b,y*b,z*b); 
	  }
	  Vec mult(const Vec &b) const 
	  { 
		  return Vec(x*b.x,y*b.y,z*b.z);
	  }
	  Vec& norm()
	  { 
		  return *this = *this * (1/sqrt(x*x+y*y+z*z));
	  }
	  double dot(const Vec &b) const
	  { 
		  return x*b.x+y*b.y+z*b.z;
	  }
	  Vec operator%(Vec&b)
	  {
		  return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);
	  }
};

struct Ray 
{ 
	Vec origin, direction;
	Ray(Vec o_, Vec d_) : origin(o_), direction(d_) {}
};

enum reflectance
{ 
	diff, spec
};  

struct Sphere
{
	double rad;       
	Vec pos;
	Vec emt;
	Vec clr;      
	reflectance refl;      
  
	Sphere(double radius, Vec position, Vec emittance, Vec color, reflectance reflected):
	rad(radius), pos(position), emt(emittance), clr(color), refl(reflected) {}
  
	double intersect(const Ray &r) const 
	{	
		// return distance or 0 if no hit
		// formula= t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		Vec op;
		op = pos-r.origin; 
		double t,b;
		double dcrm;
		double eps;
	
		eps=1e-4;
		b =op.dot(r.direction);
		dcrm=b*b-op.dot(op)+rad*rad;
		if (dcrm<0)
			return 0;
		else
			dcrm=sqrt(dcrm);
		t=b-dcrm;
		if(t>eps)
			return t;
		else 
		{
			t=b+dcrm;
			if(t>eps)
				return t;
			else
				return 0;
		}
	}

};

Vec calculateRadiance(const Ray &r, int length, unsigned short *Random,int E=1);
