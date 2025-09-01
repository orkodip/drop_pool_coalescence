/*UNIVERSAL HEADER FILE FOR DIFFERENT CFD APPLICATIONS

CREATED BY - ORKODIP MOOKHERJEE FOR TESTING PURPOSE

CONTENTS

1) SIGN FUNCTION
2) THOMAS' TDMA ALGORITHM WITH BACKWARD SUBSTITUTION
3) MAXIMUM OF TWO NUMBERS
4) MINIMUM OF TWO NUMBERS
5) CHECK WHETHER A NUMBER IS GREATER THAN ZERO OR NOT
6) MAXIMUM OF THREE NUMBERS
7) MINIMUM OF THREE NUMBERS
8) VECTOR CLASS (USES THE CONCEPT OF OPERATOR OVERLOADING)
9) MINMOD FUNCTION OF TWO NUMBERS
10) REAL ROOTS OF A CUBIC EQUATION (CARDANO'S METHOD OR TRIGONOMETRIC SUBSTITUTION)
11) REAL ROOTS OF A QUADRATIC EQUATION

*/

#define SGN(A) ((A)>=0.0 ? 1.0 : -1.0)

void tdma_bs(int n, double *a, double *b, double *c, double *d)
{
	//TDMA ALGORITHM STARTS
	for (int i=0;i<n;i++)
	{
		if (i==0)
		{
			c[0]/=b[0];
			d[0]/=b[0];
		}
		else
		{
			c[i]/=(b[i]-a[i]*c[i-1]);
			d[i]=(d[i]-a[i]*d[i-1])/(b[i]-a[i]*c[i-1]);
		}
	}
	//TDMA ALGORITHM ENDS

	//BACKWARD SUBSTITUTION BEGINS

	a[n-1]=d[n-1];
	for (int i=(n-2);i>=0;i--)
		a[i]=(d[i]-c[i]*a[i+1]);

	//BACKWARD SUBSTITUTION ENDS
}

double MAX2 (double a, double b)
{
	if (a>=b)
		return a;
	else
		return b;
}

double MIN2 (double a, double b)
{
	if (a<=b)
		return a;
	else
		return b;
}

double CHK(double a)
{
	if (a>=0.0)
		return 1.0;
	else
		return 0.0;
}

double MAX3(double a,double b,double c)
{
	if((a>=b)&&(a>=c)) return a;
	else if((b>=a)&&(b>=c)) return b;
	else return c;
}

double MIN3(double a,double b,double c)
{
	if((a<=b)&&(a<=c)) return a;
	else if((b<=a)&&(b<=c)) return b;
	else return c;
}

class VECTOR
{
	public:
			double x,y,z;	//x,y,z components of the vector
			VECTOR() {x=0.0; y=0.0; z=0.0;}	//null vector initialization
			VECTOR(double a,double b,double c=0.0) {x=a; y=b; z=c;}	//initialize the vector with prescribed values
			double mag() {return sqrt(x*x+y*y+z*z);}	//return magnitude of the vector
			VECTOR unit();	//return unit vector
			VECTOR operator+(VECTOR V);	//vector addition
			VECTOR operator-();	//change the direction of the vector
			VECTOR operator-(VECTOR V);	//vector substraction
			friend double dot(VECTOR V1,VECTOR V2);	//dot product
			friend VECTOR operator*(double m,VECTOR V);	//scalar multiplication (type 1)
			friend VECTOR operator*(VECTOR V,double m);	//scalar multiplication (type 2)
			friend VECTOR operator*(VECTOR V1,VECTOR V2);	//cross product
};
VECTOR VECTOR::unit()
{
	double T; VECTOR uV;
	if((abs(x)<=EPS)&&(abs(y)<=EPS)) { x=1.0; y=0.0; }
	if(abs(x)>abs(y))	//safeguard against overflow or underflow
	{
		T=y/x; uV.x=SGN(x)*pow((1.0+T*T),-0.5); uV.y=T*uV.x;
	}
	else
	{
		T=x/y; uV.y=SGN(y)*pow((1.0+T*T),-0.5); uV.x=T*uV.y;
	}
	return uV;
}
VECTOR VECTOR::operator+(VECTOR V)
{
	VECTOR add;
	add.x=x+V.x;
	add.y=y+V.y;
	add.z=z+V.z;
	return add;
}
VECTOR VECTOR::operator-()
{
	VECTOR s(-x,-y,-z);
	return s;
}
VECTOR VECTOR::operator-(VECTOR V)
{
	VECTOR sub;
	sub.x=x-V.x;
	sub.y=y-V.y;
	sub.z=z-V.z;
	return sub;
}
double dot(VECTOR V1,VECTOR V2) {return (V1.x*V2.x+V1.y*V2.y+V1.z*V2.z);}
VECTOR operator*(double m,VECTOR V)
{
	VECTOR mul(m*V.x,m*V.y,m*V.z);
	return mul;
}
VECTOR operator*(VECTOR V,double m)
{
	VECTOR mul(m*V.x,m*V.y,m*V.z);
	return mul;
}
VECTOR operator*(VECTOR V1,VECTOR V2)
{
	VECTOR cp;
	cp.x=V1.y*V2.z-V2.y*V1.z;	//calculate the cross product
	cp.y=-(V1.x*V2.z-V2.x*V1.z);
	cp.z=V1.x*V2.y-V2.x*V1.y;
	return cp;
}

double MINMOD(double a,double b)
{
	if(abs(a)<=abs(b)) return a;
	else return b;
}

int r_cubic(double a3,double a2,double a1,double a0,double *x)
{
	double NR,NQ;
	NR=(2.0*pow(a2,3.0)-9.0*a1*a2*a3+27.0*a0*pow(a3,2.0));
	NQ=(pow(a2,2.0)-3.0*a1*a3);
	if((0.25*pow(NR,2.0))<=pow(NQ,3.0))	//use trigonometric solution to obtain three real roots
	{
		a0=(abs(NQ)<=EPS)?0:(acos(0.5*NR/pow(NQ,1.5)));	//variable is reused
		for(int k=-1;k<=1;k++)
			x[k+1]=0.333333333*(-2.0*sqrt(NQ)*cos(0.333333333*(a0+2.0*k*M_PI))-a2)/a3;
		return 3;
	}
	else	//use Cardano's method to obtain the real root
	{
		a0=pow((abs(NR)+sqrt(pow(NR,2.0)-4.0*pow(NQ,3.0))),0.333333333);	//variables are reused
		a1=-0.264566842*SGN(NR)*a0/a3;
		NR=(abs(a0)<=EPS)?0:(-0.419973683*NQ/(SGN(NR)*a0*a3));
		x[0]=a1+NR-0.333333333*a2/a3;
		return 1;
	}
}

void r_quad(double a2,double a1,double a0,double *x)
{
	a1=-0.5*(a1+SGN(a1)*sqrt(pow(a1,2.0)-4.0*a2*a0));
	x[0]=a1/a2; x[1]=a0/a1;
}
