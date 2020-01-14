#include <iostream>
#include <string.h>
#include <cstdio>
#include <stdio.h>
#include <fstream>
#include <string>
#include <cmath>

const int it_max = 10000; // do zmiany
const int nx = 400;
const int ny = 90;
const int i1 = 200;
const int i2 = 210;
const int j_1 = 50;
const double delta = 0.01;
const double sigma = 10.0 * delta; // rozmycie
const double xa = 0.45; // poczatek srodka gestosci
const double ya = 0.45; // poczatek srodka gestosci

double v[401][91] = {0.0};
double u0[nx+1][ny+1] = {0.0};
double u1[nx+1][ny+1] = {0.0};
double vx[nx+1][ny+1] = {0.0};
double vy[nx+1][ny+1] = {0.0};
double x[nx+1] = {0.0};
double y[ny+1] = {0.0};
double deltaT = 0.0;
double D = 0.0;

void readFromFile()
{
  std::ifstream file("psi.dat");
  std::string str; 
  while (std::getline(file, str))
  {
    int i = std::stoi(str.substr(3,5));
    int j = std::stoi(str.substr(10,11));
    v[i][j] = std::stod(str.substr(21,28));
  }
}

void polePredkosci() // podpunkt 3
{
  for(int i = 1; i <= nx-1; i++)
  {
    for(int j = 1; j <= ny-1; j++)
    {
      if(i >= i1 && i <= i2 && j >= 0 && j <= j_1) // na zakladce
      {
        vx[i][j] = 0.0;
        vy[i][j] = 0.0;
      }
      else
      {
        vx[i][j] = (v[i][j+1] - v[i][j-1])/(2.0*delta);
        vy[i][j] = -1.0*(v[i+1][j] - v[i-1][j])/(2.0*delta);
      }
    }
  }
  for(int i = 1; i <= nx-1; i++) // gorny i dolny brzeg
  {
    vx[i][0] = 0.0;
    vy[i][ny] = 0.0;
  }
  for(int j = 0; j <= ny; j++) // lewy i prawy brzeg
  {
    vx[0][j] = vx[1][j];
    vx[nx][j] = vx[nx-1][j];
  }

  // szukanie vmax;
  double vmax = 0.0;
  for(int i = 0; i <= nx; i++)
  {
    for(int j = 0; j <= ny; j++)
    {
      double tmp = abs(sqrt( pow(vx[i][j],2.0) + pow(vy[i][j],2.0) ));
      if(vmax < tmp)
        vmax = tmp;
    }
  }
  deltaT = delta/(4.0*vmax);
    deltaT /= 10.0;
std::cout << deltaT;
}

double oblicz_u(int i, int j, double mi)
{
  double wynik = 0.0;
  if(i == 0)
  {
    wynik = u0[i][j] - (deltaT/2.0)*vx[i][j]*((u0[i+1][j] - u0[nx][j])/(2.0*delta) + (u1[i+1][j] - u1[nx][j])/(2.0*delta));
    wynik -= (deltaT/2.0)*vy[i][j]*((u0[i][j+1] - u0[i][j-1])/(2.0*delta) + (u1[i][j+1] - u1[i][j-1])/(2.0*delta));
    wynik += (deltaT*D/2.0) * ((u0[i+1][j] + u0[nx][j] + u0[i][j+1] + u0[i][j-1] - 4.0*u0[i][j])/pow(delta,2.0)
                              + ( u1[i+1][j] + u1[nx][j] + u1[i][j+1] + u1[i][j-1] )/pow(delta,2.0) 
                              );
    wynik *= (1.0/(1.0 + 2.0*D*deltaT/(delta*delta)));
  }
  else if(i == nx)
  {
    wynik = u0[i][j] - (deltaT/2.0)*vx[i][j]*((u0[0][j] - u0[i-1][j])/(2.0*delta) + (u1[0][j] - u1[i-1][j])/(2.0*delta));
    wynik -= (deltaT/2.0)*vy[i][j]*((u0[i][j+1] - u0[i][j-1])/(2.0*delta) + (u1[i][j+1] - u1[i][j-1])/(2.0*delta));
    wynik += (deltaT*D/2.0) * ((u0[0][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4.0*u0[i][j])/pow(delta,2.0)
                              + ( u1[0][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1] )/pow(delta,2.0)
                              );
    wynik *= (1.0/(1.0 + 2.0*D*deltaT/(delta*delta)));
  }
  else
  {
    wynik = u0[i][j] - (deltaT/2.0)*vx[i][j]*((u0[i+1][j] - u0[i-1][j])/(2.0*delta) + (u1[i+1][j] - u1[i-1][j])/(2.0*delta));
    wynik -= (deltaT/2.0)*vy[i][j]*((u0[i][j+1] - u0[i][j-1])/(2.0*delta) + (u1[i][j+1] - u1[i][j-1])/(2.0*delta));
    wynik += (deltaT*D/2.0) * ((u0[i+1][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4.0*u0[i][j])/pow(delta,2.0)
                              + ( u1[i+1][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1] )/pow(delta,2.0)
                              );
    wynik *= (1.0/(1.0 + 2.0*D*deltaT/(delta*delta)));
  }
  return wynik;
}


void zapisMap()
{
  FILE* mapa_vx[5];
  FILE* mapa_vy[5];
  if(D <= 0.05)
  {
    mapa_vx[0] = fopen("wyniki/mapa_vx1_0.0.txt", "w");
    mapa_vx[1] = fopen("wyniki/mapa_vx2_0.0.txt", "w");
    mapa_vx[2] = fopen("wyniki/mapa_vx3_0.0.txt", "w");
    mapa_vx[3] = fopen("wyniki/mapa_vx4_0.0.txt", "w");
    mapa_vx[4] = fopen("wyniki/mapa_vx5_0.0.txt", "w");
    mapa_vy[0] = fopen("wyniki/mapa_vy1_0.0.txt", "w");
    mapa_vy[1] = fopen("wyniki/mapa_vy2_0.0.txt", "w");
    mapa_vy[2] = fopen("wyniki/mapa_vy3_0.0.txt", "w");
    mapa_vy[3] = fopen("wyniki/mapa_vy4_0.0.txt", "w");
    mapa_vy[4] = fopen("wyniki/mapa_vy5_0.0.txt", "w");
  }
  else
  {
    mapa_vx[0] = fopen("wyniki/mapa_vx1_0.1.txt", "w");
    mapa_vx[1] = fopen("wyniki/mapa_vx2_0.1.txt", "w");
    mapa_vx[2] = fopen("wyniki/mapa_vx3_0.1.txt", "w");
    mapa_vx[3] = fopen("wyniki/mapa_vx4_0.1.txt", "w");
    mapa_vx[4] = fopen("wyniki/mapa_vx5_0.1.txt", "w");
    mapa_vy[0] = fopen("wyniki/mapa_vy1_0.1.txt", "w");
    mapa_vy[1] = fopen("wyniki/mapa_vy2_0.1.txt", "w");
    mapa_vy[2] = fopen("wyniki/mapa_vy3_0.1.txt", "w");
    mapa_vy[3] = fopen("wyniki/mapa_vy4_0.1.txt", "w");
    mapa_vy[4] = fopen("wyniki/mapa_vy5_0.1.txt", "w");
  }

  for(int i = 0; i <= nx; i++)
    {
        for(int j = 0; j <= ny; j++)
        {
          for(int k = 1; k <= 1; k++)
          {
            fprintf(mapa_vx[k-1], "%.2f %.2f %.3f\n", x[i], y[j], vx[i][j]);
            fprintf(mapa_vy[k-1], "%.2f %.2f %.3f\n", x[i], y[j], vy[i][j]);
         }

        }

        fprintf(mapa_vx[0], "\n");
        fprintf(mapa_vy[0], "\n");
  }

  for(int i = 0; i < 5; i++)
  {
    fclose(mapa_vx[i]);
    fclose(mapa_vy[i]);
  }
}
void algorytmAdwekcjiDyfuzji()
{
  FILE* calka;
  FILE* sredniaPolozenia;
  if(D <= 0.05)
  {
    calka = fopen("wyniki/calka_0.0.txt", "w");
    sredniaPolozenia = fopen("wyniki/sredniaPolozenia_0.0.txt", "w");
  }
  else
  {
    calka = fopen("wyniki/calka_0.1.txt", "w");
    sredniaPolozenia = fopen("wyniki/sredniaPolozenia_0.1.txt", "w");
  }
  double t_max = it_max/delta;
  
  for(int i = 0; i <= nx; i++)
  {
    x[i] = i*delta;
  }
  for(int j = 0; j <= ny; j++)
  {
    y[j] = j*delta;
  }

  // poczatkowe u0 i u1
  for(int i = 0; i <= nx; i++)
  {
    for(int j = 0; j <= ny; j++)
    {
      u0[i][j] = exp(-1.0*(pow(x[i]-xa,2.0) + pow(y[j]-ya,2.0))/(2.0*sigma*sigma)) / (2.0*M_PI*sigma*sigma);
      u1[i][j] = 0.0;
    }
  }
  polePredkosci();

  for(int it = 1; it <= it_max; it++)
  {
    for(int i = 0; i <= nx; i++)
    {
      for(int j = 1; j <= ny; j++)
      {
        u1[i][j] = u0[i][j];
      }
    }
    for(int K = 1; K <= 20; K++)
    {
      for(int i = 0; i <= nx; i++)
      {
        for(int j = 1; j <= ny-1; j++)
        {
          if(i >= i1 && i <= i2 && j >= 0 && j <= j_1)
          {

          }
          else if(i == 0 || i == nx)
          {
            u1[i][j] = oblicz_u(i, j, 1);
          }
          else
          {
            u1[i][j] = oblicz_u(i, j, 1);
          }
          
        }
      }

      double wartoscCalki = 0.0; 
      double xsr = 0.0;
      double tn = deltaT * it;
      for(int i = 0; i <= nx; i++)
      {
        for(int j = 0; j <= ny; j++)
        {
          u0[i][j] = u1[i][j];
          wartoscCalki += u0[i][j]*delta*delta;
          xsr += x[i]*u0[i][j]*delta*delta;
        }
      }

      fprintf(calka, "%f, %f\n", tn, wartoscCalki);
      fprintf(sredniaPolozenia, "%f, %f\n", tn, xsr);
        //std::cout << it << "  " << wartoscCalki << "  " << xsr << std::endl;
    }

    
  }

  fclose(calka);
  fclose(sredniaPolozenia);
  
}

int main(void)
{
  readFromFile();
  D = 0.0;
  algorytmAdwekcjiDyfuzji();
  zapisMap();
  D = 0.1;
  algorytmAdwekcjiDyfuzji();
  zapisMap();
  return 0;
}
