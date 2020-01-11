#include <cstdio>
#include <cstdlib>
#include <cmath>

const double x0 = 0.01;
const double v0 = 0.0;
const double dt0 = 1.0;
const double S = 0.75;
const double p = 2.0;

struct xv
{
  double x;
  double v;
};

double a11(void)
{
  return 1.0;
}

double a12(double dt)
{
  return (-1.0)*dt/2.0;
}

double a21(double dt, double alfa, double xn1k, double vn1k)
{
  return (-1.0)*(dt/2.0)*(-2.0*alfa*xn1k*vn1k - 1.0);
}

double a22(double dt, double alfa, double xn1k)
{
  return 1.0 - (dt/2.0)*alfa*(1.0 - xn1k*xn1k);
}

double f(double v, double x = 0.0, double t = 0.0)
{
  return v;
}

double g(double alfa, double x, double v, double t = 0.0)
{
  return alfa*(1.0 - x*x)*v - x;
}

double F(double xn1, double xn, double vn1, double vn, double dt){

  return xn1 - xn - (dt/2)*(f(vn)+f(vn1));
}

double G(double xn1, double xn, double vn1, double vn, double dt, double alfa){

  return vn1 - vn - (dt/2)*(g(alfa,xn,vn)+g(alfa,xn1,vn1));
}


xv metoda_RK2(double xn, double vn, double dt, double alfa){
    double k1x = f(vn);
    double k1v = g(alfa, xn, vn);

    double k2x = f(vn+dt*k1v);
    double k2v = g(alfa, xn + dt*k1x, vn + dt*k1v);

    xv wynik;
    wynik.x = xn + (k1x + k2x)*dt/2.0;
    wynik.v = vn + (k1v + k2v)*dt/2.0;
    return wynik;

}


xv metoda_trapezow(double xn, double vn, double dt, double alfa){

  double xn1 = xn;
  double vn1 = vn;
  double dx, dv;

  double delta = pow(10.0,-10.0);

  do
  {
      dx = ((-1.0)*F(xn1,xn, vn1, vn, dt)*a22(dt,alfa,xn1)-(-1.0)*G(xn1,xn,vn1,vn,dt, alfa)*a12(dt)) / (a11()*a22(dt,alfa,xn1) - a12(dt)*a21(dt,alfa,xn1,vn1));

      dv = (a11()*(-1.0)*G(xn1,xn,vn1,vn,dt, alfa) - a21(dt,alfa,xn1,vn1)*(-1.0)*F(xn1,xn,vn1,vn,dt)) / (a11()*a22(dt,alfa,xn1) - a12(dt)*a21(dt,alfa,xn1,vn1));

      xn1 = xn1 + dx;
      vn1 = vn1 + dv;
      

  }while((fabs(dx) > delta) || (fabs(dv) > delta));


  xv wynik;

  wynik.x = xn1; //+ (f(vn) + f(vn1))*dt/2.0;
  wynik.v = vn1;// + (g(alfa,xn,vn) + g(alfa,xn1,vn1))*dt/2.0;

  return wynik;
}

void algorytmKontroliKroku(FILE* plik, double TOL, xv (*schemat_numeryczny)(double xn, double vn, double dt, double alfa))
{
  double t = 0.0;
  double dt = dt0;
  double xn = x0;
  double vn = v0;
  double tmax = 40.0;
  double alfa = 5.0;

  fprintf(plik, "%f %f %f %f\n", t, dt, xn, vn);

// do policzenia
  double Ex = 0.0; 
  double Ev = 0.0;
  xv tmp1, tmp2;

  do
  {
    tmp2 = schemat_numeryczny(xn,vn,dt,alfa);
    tmp2 = schemat_numeryczny(tmp2.x,tmp2.v,dt,alfa);

    //stawiamy jeden krok 2dt
    tmp1 = schemat_numeryczny(xn,vn,2*dt,alfa);

    //liczymy Ex Ev
        Ex = (tmp2.x - tmp1.x)/(pow(2.0,p) - 1.0);
        Ev = (tmp2.v - tmp1.v)/(pow(2.0,p) - 1.0);

    if(fmax(fabs(Ex),fabs(Ev)) < TOL )
    {
        t = t + 2.0*dt;
        xn = tmp2.x;
        vn = tmp2.v;
        fprintf(plik, "%f %f %f %f \n", t, dt, xn, vn);
    }
    dt = (pow(S*TOL/(fmax(fabs(Ex),fabs(Ev))), (1.0/(p+1.0))) * dt);
  }while(t < (tmax - dt));
}

int main(void)
{
  FILE* plik;
  double TOL = pow(10.0, -2.0);

  // zadanie 1
  plik = fopen("metoda_RK2_TOL_-2.txt", "w");
  algorytmKontroliKroku(plik, TOL, metoda_RK2);
  fclose(plik);

  plik = fopen("metoda_RK2_TOL_-5.txt", "w");
  TOL = pow(10.0, -5.0);
  algorytmKontroliKroku(plik, TOL, metoda_RK2);
  fclose(plik);

  // zadanie 2
  TOL = pow(10.0, -2.0);
  plik = fopen("metoda_trapezow_TOL-2.txt", "w");
  algorytmKontroliKroku(plik, TOL, metoda_trapezow);
  fclose(plik);

  plik = fopen("metoda_trapezow_TOL-5.txt", "w");
  TOL = pow(10.0, -5.0);
  algorytmKontroliKroku(plik, TOL, metoda_trapezow);
  fclose(plik);
  
  return 0;
}
