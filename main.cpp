#include "header.h"

double randomNum( unsigned short xsubi[3])
{
	return (double)rand()/(double)RAND_MAX;
}

double makeRandom( unsigned short Random[3])
{
	double r1=2*randomNum(Random);
	double dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
	return dx;
}

void savebmp (const char *filename, int w, int h, int dpi, RGB *data) {
	FILE *f;
	int k = w*h;
	int s = 4*k;
	int filesize = 54 + s;
	
	double factor = 39.375;
	int m = static_cast<int>(factor);
	
	int ppm = dpi*m;
	
	unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0,0,0, 54,0,0,0};
	unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
	
	bmpfileheader[ 2] = (unsigned char)(filesize);
	bmpfileheader[ 3] = (unsigned char)(filesize>>8);
	bmpfileheader[ 4] = (unsigned char)(filesize>>16);
	bmpfileheader[ 5] = (unsigned char)(filesize>>24);
	
	bmpinfoheader[ 4] = (unsigned char)(w);
	bmpinfoheader[ 5] = (unsigned char)(w>>8);
	bmpinfoheader[ 6] = (unsigned char)(w>>16);
	bmpinfoheader[ 7] = (unsigned char)(w>>24);
	
	bmpinfoheader[ 8] = (unsigned char)(h);
	bmpinfoheader[ 9] = (unsigned char)(h>>8);
	bmpinfoheader[10] = (unsigned char)(h>>16);
	bmpinfoheader[11] = (unsigned char)(h>>24);
	
	bmpinfoheader[21] = (unsigned char)(s);
	bmpinfoheader[22] = (unsigned char)(s>>8);
	bmpinfoheader[23] = (unsigned char)(s>>16);
	bmpinfoheader[24] = (unsigned char)(s>>24);
	
	bmpinfoheader[25] = (unsigned char)(ppm);
	bmpinfoheader[26] = (unsigned char)(ppm>>8);
	bmpinfoheader[27] = (unsigned char)(ppm>>16);
	bmpinfoheader[28] = (unsigned char)(ppm>>24);
	
	bmpinfoheader[29] = (unsigned char)(ppm);
	bmpinfoheader[30] = (unsigned char)(ppm>>8);
	bmpinfoheader[31] = (unsigned char)(ppm>>16);
	bmpinfoheader[32] = (unsigned char)(ppm>>24);
	
	f = fopen(filename,"wb");
	
	fwrite(bmpfileheader,1,14,f);
	fwrite(bmpinfoheader,1,40,f);
	
	for (int i = 0; i < k; i++) {
		RGB rgb = data[i];
		
		double red = (data[i].r)*255;
		double green = (data[i].g)*255;
		double blue = (data[i].b)*255;
		
		unsigned char color[3] = {(int)floor(blue),(int)floor(green),(int)floor(red)};
		
		fwrite(color,1,3,f);
	}
	
	fclose(f);
}

Sphere spheres[] = 
{	
	Sphere(180,
		Vec(50,-160,47),
		Vec(0,0,0),
		Vec(0.05,0.74,0.25)*.999,
		diff),
	Sphere(16.5,
		Vec(75,71,47),
		Vec(0,0,0),
		Vec(1,1,1)*.999,
		diff),
	Sphere(16.5,
		Vec(27,34,47),
		Vec(0,0,0),
		Vec(1,1,1)*.999,
		diff),
	Sphere(1.0,
		Vec(50,100,81.6),
		Vec(30,30,30)*100,
		Vec(0,0,0),
		diff),
};

int numSpheres = sizeof(spheres)/sizeof(Sphere);

inline double clamp(double x)
{ 
	if(x<0)
		x=0;
	else if (x>1)
		x=1;
	else 
		x=x;
	return  x;
}

inline int makeInt(double x)
{
	double y = clamp(x);
	double y1 = pow(y,1/2.2);
	int y2 = int(y1*255+.5);
	return y2; 
}

inline bool intersectS(const Ray &r, double &t, int &id)
{
  double n;
  n= sizeof(spheres)/sizeof(Sphere);
  double d;
  double inf;
  inf=t=1e20;
  for(int i=int(n);i--;)
	  if((d=spheres[i].intersect(r)) && d<t )
	  {
		  t=d;
		  id=i;
	  }
  return t<inf;
}

Vec diffReflection(const Sphere obj, const Ray &r, int length,double t, int id, unsigned short *Random,int E=1)
{
	Vec hitPoint, Snormal;
	Vec Snl, Scolor;
	double r1, r2, r2s;
	Vec w, u, v, d;

	hitPoint=r.origin+r.direction*t;
	Snormal=(hitPoint-obj.pos).norm();
	if(Snormal.dot(r.direction)<0)
		Snl=Snormal;
	else
		Snl=Snormal*-1;
	Scolor=obj.clr;
		
	r1=2*Pi*randomNum(Random);
	r2=randomNum(Random);
	r2s=sqrt(r2);
		
	w=Snl;
	if(fabs(w.x)>0.1)
		u=Vec(0,1,0)%w;
	else
		u=Vec(1,0,0)%w;
	u=u.norm();
	v=w%u;
	d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2));
	d=d.norm();
    Vec e;
    for (int i=0; i<numSpheres; i++)
	{
      const Sphere &s = spheres[i];
      if (s.emt.x<=0 && s.emt.y<=0 && s.emt.z<=0) 
		  continue; 
      
      Vec sw, su, sv;

	  sw=s.pos-hitPoint;
	  if(fabs(sw.x)>.1)
		  su=Vec(0,1,0)%sw;
	  else
		  su=Vec(1,0,0)%sw;
	  su=su.norm();
	  sv=sw%su;

      double cos_a_max, cos_a, sin_a;
	  double eps1, eps2;
	  double srad2,hspos2;
	  double phi;

	  srad2=s.rad*s.rad;
	  hspos2=(hitPoint-s.pos).dot(hitPoint-s.pos);
	  cos_a_max = sqrt(1-srad2/hspos2);
      eps1 = randomNum(Random); 
	  eps2 = randomNum(Random);
      cos_a = 1-eps1+eps1*cos_a_max;
      sin_a = sqrt(1-cos_a*cos_a);
      phi = 2*Pi*eps2;

      Vec l = su*cos(phi)*sin_a + sv*sin(phi)*sin_a + sw*cos_a;
      l=l.norm();
      
	  if (intersectS(Ray(hitPoint,l), t, id) && id==i)
	  {
        double omega;
		omega = 2*Pi*(1-cos_a_max);
		Vec semto=s.emt*l.dot(Snl)*omega;
        e = e + Scolor.mult(semto)*Pi_inv; 
      }
    }
	Vec nradiance= calculateRadiance(Ray(hitPoint,d),length,Random,0);
    Vec snradiance= Scolor.mult(nradiance);
	Vec res=obj.emt*E+e+snradiance;
    return res;
}

Vec calculateRadiance(const Ray &r, int length, unsigned short *Random,int E)
{
  double t;                               
  int id=0;                               
  if (!intersectS(r, t, id)) 
	  return Vec();

  const Sphere &obj = spheres[id];

  Vec hitPoint, Snormal;
  Vec Snl, Scolor;

  hitPoint=r.origin+r.direction*t;
  Snormal=(hitPoint-obj.pos).norm();
  if(Snormal.dot(r.direction)<0)
	  Snl=Snormal;
  else
	  Snl=Snormal*-1;
  Scolor=obj.clr;

  double p;
  if(Scolor.x>Scolor.y && Scolor.x>Scolor.z)
	  p= Scolor.x;
  else if(Scolor.y>Scolor.z)
	  p= Scolor.y;
  else
	  p= Scolor.z;

  if (++length>5||!p)
  {
	  if (randomNum(Random)<p)
		  Scolor=Scolor*(1/p); 
	  else 
		  return obj.emt*E;
  }

  if (obj.refl == diff)
  {   
	  Vec res= diffReflection(obj,r,length,t,id,Random,1);
	  return res;
  } 
  else 
  {
	  Vec sNdotR= Snormal*2*Snormal.dot(r.direction);
	  Vec sDir= r.direction-sNdotR;
	  Ray sRay= Ray(hitPoint,sDir);
	  Vec sRadiance= calculateRadiance(sRay,length,Random);
	  Vec sC= Scolor.mult(sRadiance);
	  Vec res=obj.emt + sC; 
	  return res;
  }
}

int main(int argc, char *argv[])
{
	int w=640;
	int h=480;
	int samps = 10;

	Vec CamOrigin = Vec(50,52,295.6);
	Vec CamDirection = Vec(0,-0.042612,-1);
	CamDirection = CamDirection.norm();
	Ray cam(CamOrigin,CamDirection);

	double wh=w*0.5135/h;
	Vec crx = Vec(wh,0,0);
	Vec cry = crx%cam.direction;
	cry = cry.norm()*0.5135;

	Vec rd;
	Vec *cres=new Vec[w*h];
	printf("Ray Tracing Project:\n\n");

#pragma omp parallel for schedule(dynamic, 1) private(r)  

  for (int y=0; y<h; y++)
  {                       
    fprintf(stderr,"\rRendering ... %3.2f%%",100.*y/(h-1));
	unsigned short Random[3]={0,0,y*y*y};
    for (unsigned short x=0; x<w; x++)
	{
		int  i=(h-y-1)*w+x;
		for (int u=0; u<2; u++)
		{     
			for (int v=0; v<2; v++)
			{   
				rd=Vec();
				for (int s=0; s<samps; s++)
				{
					double dx, dy;
					Vec dir;
					dx=makeRandom(Random);
					dy=makeRandom(Random);            
					Vec vcx = crx*(((v+0.5+dx)/2+x)/w-0.5);
					Vec vcy = cry*(((u+0.5+dy)/2+y)/h-0.5);
					dir =  vcx + vcy + cam.direction;
					Ray ranRay=Ray(cam.origin+dir*140,dir.norm());            
					rd = rd + calculateRadiance(ranRay,0,Random)*(1./samps);          
				} 
				double rdx=clamp(rd.x);
				double rdy=clamp(rd.y);
				double rdz=clamp(rd.z);
				cres[i] = cres[i] + Vec(rdx,rdy,rdz)*0.25;
			}
		}
	}
  }
  RGB *pixels = new RGB[w*h];
	for (int i=0; i<w*h; i++)
	{
		pixels[i].r=cres[(w*h)-i].x;
		pixels[i].g=cres[(w*h)-i].y;
		pixels[i].b=cres[(w*h)-i].z;
	}
	savebmp("scene.bmp",w,h,72,pixels);
	delete cres;
	delete pixels;
}
