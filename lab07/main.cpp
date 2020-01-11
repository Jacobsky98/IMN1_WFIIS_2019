#include <cstdio>
#include <cmath>
#include <stdio.h>

const int IT_MAX = 20000;

const double delta = 0.01;
const double ro = 1.0;
const double u = 1.0;
const int nx = 200;
const int ny = 90;
const int i_1 = 50;
const int j_1 = 55;

double psi_wzor8(double psi[nx+1][ny+1], double ksi[nx+1][ny+1], int i, int j)
{
  return 0.25*(psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - delta*delta*ksi[i][j]);
}

double ksi_wzor9(double psi[nx+1][ny+1], double ksi[nx+1][ny+1], int i, int j, double omega)
{
  return 0.25*(ksi[i+1][j] + ksi[i-1][j] + ksi[i][j+1] + ksi[i][j-1])
          - omega*(ro/(16.0*u))*( (psi[i][j+1] - psi[i][j-1])*(ksi[i+1][j] - ksi[i-1][j]) - (psi[i+1][j] - psi[i-1][j])*(ksi[i][j+1] - ksi[i][j-1]));
}

void WB_ksi(double psi[nx+1][ny+1], double ksi[nx+1][ny+1], double* x, double* y, double Qwe, double Qwy)
{
  // brzeg A (wejscie)
  for(int j = j_1; j <= ny; j++)
    ksi[0][j] = (Qwe/(2.0*u))*(2.0*y[j] - y[j_1] - y[ny]);

  // brzeg C (wyjscie)
  for(int j = 0; j <= ny; j++)
    ksi[nx][j] = (Qwy/(2.0*u))*(2.0*y[j] - y[ny]);

  // brzeg B
  for(int i = 1; i <= nx-1; i++)
    ksi[i][ny] = (psi[i][ny-1] - psi[i][ny])*2.0/pow(delta, 2.0);

  // brzeg D
  for(int i = i_1+1; i <= nx-1; i++)
    ksi[i][0] = (psi[i][1] - psi[i][0])*2.0/pow(delta, 2.0);

  // brzeg E
  for(int j = 1; j <= j_1; j++)
    ksi[i_1][j] = (psi[i_1+1][j] - psi[i_1][j])*2.0/pow(delta, 2.0);

  // brzeg F
  for(int i = 1; i <= i_1; i++)
    ksi[i][j_1] = (psi[i][j_1+1] - psi[i][j_1])*2.0/pow(delta, 2.0);

  // wierzcholek E/F
  ksi[i_1][j_1] = 0.5*(ksi[i_1-1][j_1] + ksi[i_1][j_1-1]);
}

void WB_psi(double psi[nx+1][ny+1], double ksi[nx+1][ny+1], double* x, double* y, double Qwe, double Qwy)
{
  // brzeg A (wejscie)
  for(int j = j_1; j <= ny; j++)
    psi[0][j] = (Qwe/(2.0*u))*( pow(y[j],3.0)/3.0 - pow(y[j],2.0)*(y[j_1]+y[ny])/2.0 + y[j]*y[j_1]*y[ny] );

  // brzeg C (wyjscie)
  for(int j = 0; j <= ny; j++)
    psi[nx][j] = (Qwy/(2.0*u))*(pow(y[j],3.0)/3.0 - pow(y[j],2.0)*y[ny]/2.0) + Qwe*pow(y[j_1],2.0)*(3*y[ny]-y[j_1])/(12.0*u);

  // brzeg B
  for(int i = 1; i <= nx-1; i++)
    psi[i][ny] = psi[0][ny];

  // brzeg D
  for(int i = i_1; i <= nx-1; i++)
    psi[i][0] = psi[0][j_1];

  // brzeg E
  for(int j = 1; j <= j_1; j++)
    psi[i_1][j] = psi[0][j_1];

  // brzeg F
  for(int i = 1; i <= i_1; i++)
    psi[i][j_1] = psi[0][j_1];
}

void relaksacja(double psi[nx+1][ny+1], double ksi[nx+1][ny+1], double* x, double* y, double Qwe, double Qwy)
{
  double omega = 0;
  for(int it=0; it < IT_MAX; it++)
  {
    if(it < 2000)
      omega= 0.0;
    else
      omega = 1.0;

    for(int i = 1; i <= nx-1; i++)
    {
      for(int j = 1; j <= ny-1; j++)
      {
        bool srodek = (i <= i_1 && j <= j_1);//(i>i_1 && i < nx && j < ny && j > 0) || (j > j_1 && j < ny && i > 0 && i <= i_1);
        if(!srodek)
        {
          psi[i][j] = psi_wzor8(psi,ksi,i,j);
          ksi[i][j] = ksi_wzor9(psi,ksi,i,j,omega);
        }
      }
    }

    WB_ksi(psi, ksi, x, y, Qwe, Qwy);

  }
}


void relaksacjaNS(double Qwe, int numerZadania)
{
  FILE* konturPsi; 
  FILE* konturKsi;
  FILE* mapaPredkoscPoziomaU;
  FILE* mapaPredkoscPionowaV;
  if(numerZadania == 4)
  {
    konturPsi = fopen("wyniki/wykresKonturowyPsi-1000.txt", "w");
    konturKsi = fopen("wyniki/wykresKonturowyKsi-1000.txt", "w");
    mapaPredkoscPoziomaU = fopen("wyniki/mapaPozioma-1000.txt", "w");
    mapaPredkoscPionowaV = fopen("wyniki/mapaPionowa-1000.txt", "w");
  }
  else if(numerZadania == 5)
  {
    konturPsi = fopen("wyniki/wykresKonturowyPsi-4000.txt", "w");
    konturKsi = fopen("wyniki/wykresKonturowyKsi-4000.txt", "w");
    mapaPredkoscPoziomaU = fopen("wyniki/mapaPozioma-4000.txt", "w");
    mapaPredkoscPionowaV = fopen("wyniki/mapaPionowa-4000.txt", "w");
  }
  else if(numerZadania == 6)
  {
    konturPsi = fopen("wyniki/wykresKonturowyPsi+4000.txt", "w");
    konturKsi = fopen("wyniki/wykresKonturowyKsi+4000.txt", "w");
    mapaPredkoscPoziomaU = fopen("wyniki/mapaPozioma+4000.txt", "w");
    mapaPredkoscPionowaV = fopen("wyniki/mapaPionowa+4000.txt", "w");
  }
  else
  {
    printf("Zly numer zadania w wywolaniu funkcji\n");
    return;
  }
  double x[nx+1];
  double y[ny+1];
  double ksi[nx+1][ny+1];
  double psi[nx+1][ny+1];

  for(int i = 0; i <= nx; i++)
  {
    for(int j = 0; j <= ny; j++)
    {
      ksi[i][j] = 0.0;
      psi[i][j] = 0.0;
    }
  }

  for(int i=0; i <= nx; i++)
    x[i] = delta*i;

  for(int j=0; j <= ny; j++)
    y[j] = delta*j;

  double Qwy = Qwe*(pow(y[ny],3.0) - pow(y[j_1],3.0) - 3.0*y[j_1]*pow(y[ny],2.0) + 3.0*pow(y[j_1],2.0)*y[ny])/pow(y[ny],3.0);


  
  WB_psi(psi, ksi, x, y, Qwe, Qwy);
  // WB_ksi(psi, ksi, x, y, Qwe, Qwy);
  relaksacja(psi, ksi, x, y, Qwe, Qwy);

  
  for(int i = 0; i < nx+1; i += 1)
  {
    for(int j = 0; j < ny+1; j += 1)
    {


    double c=psi[i][j];
      if(i <= i_1 && j <= j_1)c=psi[0][j_1];
      fprintf(konturPsi, "%d  %d  %f \n ",i,j, c);
      fprintf(konturKsi, "%d  %d  %f \n ",i,j, ksi[i][j]);
      fprintf(mapaPredkoscPoziomaU, "%f %f %f \n", (double)i*delta, (double)j*delta, psi[i][j]);
     fprintf(mapaPredkoscPionowaV, "%f %f %f \n", (double)i*delta, (double)j*delta, ksi[i][j]);
    }
    fprintf(konturPsi, "\n");
    fprintf(konturKsi, "\n");
    fprintf(mapaPredkoscPoziomaU, "\n");
    fprintf(mapaPredkoscPionowaV, "\n");
  }

  fclose(konturPsi); 
  fclose(konturKsi);
  fclose(mapaPredkoscPoziomaU);
  fclose(mapaPredkoscPionowaV);
}

int main(void)
{
  relaksacjaNS(-1000, 4);
  relaksacjaNS(-4000, 5);
  relaksacjaNS(4000, 6);
  return 0;
}