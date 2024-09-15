#include <stdio.h> 
#include <sys/time.h>
#include <math.h>
#include "mycom.h"
double interm(double x); 
double interm(double x) { 
    return (1 - x)/(1 + x);
}
double plus(double x); 
double plus(double x) { 
    return (((1 - x)/(1 + x)) + x);
}
double minus(double x); 
double minus(double x) {
    return (((1 - x)/(1 + x)) - x);
}
double divid(double x); 
double divid(double x) { 
    return (((1 - x)/(1 + x)) / x);
}
double mult(double x); 
double mult(double x) { 
    return ((1 - x)/(1 + x) * x);
}
double poww(double x); 
double poww(double x) { 
    return pow((1 - x)/(1 + x), x);
}
double expp(double x); 
double expp(double x) { 
    return exp((1 - x)/(1 + x));
}
double logg(double x); 
double logg(double x) { 
    return log((1 - x)/(1 + x));
}
double sinn(double x); 
double sinn(double x) { 
    return sin((1 - x)/(1 + x));
}
int main(int argc, char *argv[]) {
 int nc=1000000000;
 double t1, t2, dt, sum;
 t1 = mytime(0); sum = integrate(interm,0.0,1.0,nc); t1 = mytime(1);
 t2 = mytime(0); sum = integrate(plus,0.0,1.0,nc); t2 = mytime(1);
 dt = 1.0/dabs(t2-t1);
 printf("Time: %lf %lf sec plus perf.: %le GFlops\n",t1,t2,dt);
 t2 = mytime(0); sum = integrate(minus,0.0,1.0,nc); t2 = mytime(1);
 dt = 2.0/dabs(t2-t1);
 printf("Time: %lf %lf sec minus perf.: %le GFlops\n",t1,t2,dt);
 t2 = mytime(0); sum = integrate(divid,0.0,1.0,nc); t2 = mytime(1);
 dt = 2.0/dabs(t2-t1);
 printf("Time: %lf %lf sec div perf.: %le GFlops\n",t1,t2,dt);
 t2 = mytime(0); sum = integrate(mult,0.0,1.0,nc); t2 = mytime(1);
 dt = 2.0/dabs(t2-t1);
 printf("Time: %lf %lf sec mult perf.: %le GFlops\n",t1,t2,dt);
 t2 = mytime(0); sum = integrate(poww,0.0,1.0,nc); t2 = mytime(1);
 dt = 2.0/dabs(t2-t1);
 printf("Time: %lf %lf sec pow perf.: %le GFlops\n",t1,t2,dt);
 t2 = mytime(0); sum = integrate(expp,0.0,1.0,nc); t2 = mytime(1);
 dt = 2.0/dabs(t2-t1);
 printf("Time: %lf %lf sec exp perf.: %le GFlops\n",t1,t2,dt);
 t2 = mytime(0); sum = integrate(logg,0.0,1.0,nc); t2 = mytime(1);
 dt = 2.0/dabs(t2-t1);
 printf("Time: %lf %lf sec log perf.: %le GFlops\n",t1,t2,dt);
 t2 = mytime(0); sum = integrate(sinn,0.0,1.0,nc); t2 = mytime(1);
 dt = 2.0/dabs(t2-t1);
 printf("Time: %lf %lf sec sin perf.: %le GFlops\n",t1,t2,dt);
 return 0;
}