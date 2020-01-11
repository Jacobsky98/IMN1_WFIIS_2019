#include <cstdlib>
#include <cstdio>
#include <cmath>


const double t_min = 0.0;
const double t_max = 5.0;
const double lambda = -1.0;

void metodaJawnaEulera(FILE* bledy, FILE* rozwiazanie, double dt)
{

  double y0 = 1.0;
  double y = 1.0;
  double yd = 1.0;

  for(double t = t_min; t <= t_max; t+=dt)
  {
    yd = exp(lambda*t	);
    fprintf(rozwiazanie, "%f %f %f\n", t, y, yd);
    fprintf(bledy, "%f %f\n", t, y-yd);
    y = y + dt*lambda*y;

    y0 = y;
  }

}

void metodaJawnaRK2(FILE* bledy, FILE* rozwiazanie, double dt)
{

  double y0 = 1.0;
  double y = 1.0;
  double yd = 1.0;
  for(double t = t_min; t <= t_max; t+=dt)
  {
    
    
    yd = exp(lambda*t);
    fprintf(rozwiazanie, "%f %f %f\n", t, y, yd);
    fprintf(bledy, "%f %f\n", t, y-yd);

    double k1 = lambda*y;
    double k2 = lambda*(y+(dt*k1));
    y = y + (dt/2.0)*(k1+k2);
    y0 = y;
  }

}

void metodaJawnaRK4(FILE* bledy, FILE* rozwiazanie, double dt)
{

  double y0 = 1.0;
  double y = 1.0;
  double yd = 1.0;
  for(double t = t_min; t <= t_max; t+=dt)
  {
    
    yd = exp(lambda*t);
    fprintf(rozwiazanie, "%f %f %f\n", t, y, yd);
    fprintf(bledy, "%f %f\n", t, y-yd);
    double k1 = lambda*y;
    double k2 = lambda*(y+(dt*k1/2.0));
    double k3 = lambda*(y+(dt*k2/2.0));
    double k4 = lambda*(y+(dt*k3));
    y = y + (dt/6.0)*(k1 + (2.0*k2) + (2.0*k3) + k4);
    y0 = y;
  }

}

double obliczV(double omega, double t)
{
  return 10.0*sin(omega*t);
}

void RRZ2rzedu(FILE * plikQ, FILE* plikI, double mnoznik_omega)
{
  double R = 100.0;
  double L = 0.1;
  double C = 0.001;
  double omega0 = 1/sqrt(L*C);
  double omega = mnoznik_omega*omega0;
  double T0 = 2.0*M_PI/omega0;
  double t0 = 0.0, tk = 4*T0, dt = 0.001;

  double Q = 0.0;
  double I = 0.0;


  for(double t = t0; t <= tk; t+=dt)
  {
    double V = obliczV(omega, t);

    double k1Q = I;
    double k1I = V/L - Q/(L*C) - R*I/L;

    double k2Q = I + k1I*dt/2.0;
    double k2I = obliczV(omega, t+(dt/2.0))/L + (Q + k1Q*dt/2.0)/(L*C) - (R/L)*(I + k1I*dt/2.0);

    double k3Q = I + k2I*dt/2.0;
    double k3I = obliczV(omega, t+(dt/2.0))/L + (Q + k2Q*dt/2.0)/(L*C) - (R/L)*(I + k2I*dt/2.0);

    double k4Q = I + dt*k3I;
    double k4I = obliczV(omega, t+dt)/L + (Q + k3Q*dt)/(L*C) - (R/L)*(I + k3I*dt);




    fprintf(plikQ, "%f %f\n", t, Q);
    fprintf(plikI, "%f %f\n", t, I);

    Q = Q + (dt/6.0)*(k1Q + 2.0*k2Q + 2.0*k3Q + k4Q);
    I = I + (dt/6.0)*(k1I + 2.0*k2I + 2.0*k3I + k4I);
  }

  

}

int main(void)
{

// zadanie 1
  FILE* bledy, *rozwiazanie;
  bledy = fopen("zad_1_bledy_0.01.dat", "w");
  rozwiazanie = fopen("zad_1_rozwiazanie_0.01.dat", "w");
  metodaJawnaEulera(bledy, rozwiazanie, 0.01);
  fclose(bledy);
  fclose(rozwiazanie);

  bledy = fopen("zad_1_bledy_0.1.dat", "w");
  rozwiazanie = fopen("zad_1_rozwiazanie_0.1.dat", "w");
  metodaJawnaEulera(bledy, rozwiazanie, 0.1);
  fclose(bledy);
  fclose(rozwiazanie);
  
  bledy = fopen("zad_1_bledy_1.dat", "w");
  rozwiazanie = fopen("zad_1_rozwiazanie_1.dat", "w");
  metodaJawnaEulera(bledy, rozwiazanie, 1);
  fclose(bledy);
  fclose(rozwiazanie);
  
// zadanie 2

  bledy = fopen("zad_2_bledy_0.01.dat", "w");
  rozwiazanie = fopen("zad_2_rozwiazanie_0.01.dat", "w");
  metodaJawnaRK2(bledy, rozwiazanie, 0.01);
  fclose(bledy);
  fclose(rozwiazanie);

  bledy = fopen("zad_2_bledy_0.1.dat", "w");
  rozwiazanie = fopen("zad_2_rozwiazanie_0.1.dat", "w");
  metodaJawnaRK2(bledy, rozwiazanie, 0.1);
  fclose(bledy);
  fclose(rozwiazanie);
  
  bledy = fopen("zad_2_bledy_1.dat", "w");
  rozwiazanie = fopen("zad_2_rozwiazanie_1.dat", "w");
  metodaJawnaRK2(bledy, rozwiazanie, 1);
  fclose(bledy);
  fclose(rozwiazanie);


// zadanie 3

  bledy = fopen("zad_3_bledy_0.01.dat", "w");
  rozwiazanie = fopen("zad_3_rozwiazanie_0.01.dat", "w");
  metodaJawnaRK4(bledy, rozwiazanie, 0.01);
  fclose(bledy);
  fclose(rozwiazanie);

  bledy = fopen("zad_3_bledy_0.1.dat", "w");
  rozwiazanie = fopen("zad_3_rozwiazanie_0.1.dat", "w");
  metodaJawnaRK4(bledy, rozwiazanie, 0.1);
  fclose(bledy);
  fclose(rozwiazanie);
  
  bledy = fopen("zad_3_bledy_1.dat", "w");
  rozwiazanie = fopen("zad_3_rozwiazanie_1.dat", "w");
  metodaJawnaRK4(bledy, rozwiazanie, 1);
  fclose(bledy);
  fclose(rozwiazanie);


// zadanie 4

  FILE *plikQ, *plikI;

  plikQ = fopen("zad_4_Q_0.5.dat", "w");
  plikI = fopen("zad_4_I_0.5.dat", "w");
  RRZ2rzedu(plikQ, plikI, 0.5);
  fclose(plikQ);
  fclose(plikI);

  plikQ = fopen("zad_4_Q_0.8.dat", "w");
  plikI = fopen("zad_4_I_0.8.dat", "w");
  RRZ2rzedu(plikQ, plikI, 0.8);
  fclose(plikQ);
  fclose(plikI);
  
  plikQ = fopen("zad_4_Q_1.dat", "w");
  plikI = fopen("zad_4_I_1.dat", "w");
  RRZ2rzedu(plikQ, plikI, 1);
  fclose(plikQ);
  fclose(plikI);

  plikQ = fopen("zad_4_Q_1.2.dat", "w");
  plikI = fopen("zad_4_I_1.2.dat", "w");
  RRZ2rzedu(plikQ, plikI, 1.2);
  fclose(plikQ);
  fclose(plikI);


  return 0;
  
}
