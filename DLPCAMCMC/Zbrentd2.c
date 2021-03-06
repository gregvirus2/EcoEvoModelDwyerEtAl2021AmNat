#include <math.h>
#define NRANSI
//#include "nrutil.h"
#define ITMAX 100
#define EPS 3.0e-8

double zbrentd2(double (*func)(double, double,double,double), double k, double N, double Z, double x1, double x2, double tol)
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(k,N,Z,a),fb=(*func)(k,N,Z,b),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
	/*
	 printf("root not bracketed...\n");
	*/
	 return -1.0;
	}
	 /*
		nrerror("Root must be bracketed in zbrent");
	 */
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=(*func)(k,N,Z,b);
	}
//	nrerror("Maximum number of iterations exceeded in zbrent");
	return 0.0;
}

double nuEqmfx(double slope,double bherit, double V, double rbase, double i){

	double noobar = (1 - rbase*(1-i))/(rbase*slope*pow(1-i,V+1));
	double bb = bherit;
	
	double numerator = pow(1-i,bb*V)*(1 + slope*noobar*(bb*V+1)*pow(1-i,bb*V));
	double denominator = 1 + slope*noobar*pow(1-i,bb*V);

	return(1 - (numerator/denominator));

}

//double zbrentd2(double (*func)(double, double,double,double), double k, double N, double Z, double x1, double x2, double tol)


double zbrentd3(double (*func)(double, double,double,double, double), double slope, double bherit, double V,double rbase, double x1, double x2, double tol)
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(slope,bherit,V,rbase,a),fb=(*func)(slope,bherit,V,rbase,b),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
	/*
	 printf("root not bracketed...\n");
	*/
	 return -1.0;
	}
	 /*
		nrerror("Root must be bracketed in zbrent");
	 */
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		//fb=(*func)(k,N,Z,b);
		fb=(*func)(slope,bherit,V,rbase,b);		
	}
//	nrerror("Maximum number of iterations exceeded in zbrent");
	return 0.0;
}


#undef ITMAX
#undef EPS
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 'j'$79L. */
