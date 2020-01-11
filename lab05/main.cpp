#include <cstdio>
#include <cstdlib>
#include <cmath>

const double delta = 0.2;
const int nx = 128;
const int ny = 128;
const double xMax = delta*(double)nx;
const double yMax = delta*(double)ny;
const double TOL = pow(10.0, -8.0);



double calka(double V[][ny+1], int k)
{
  double s = 0.0;
  for(int i = 0; i < nx-k+1; i += k)
  {
    for(int j = 0; j < ny-k+1; j += k)
    {
      s += 0.5*pow(k*delta, 2.0)*(pow((V[i+k][j] - V[i][j])/(2.0*k*delta) + (V[i+k][j+k] - V[i][j+k])/(2.0*k*delta), 2.0) + pow((V[i][j+k] - V[i][j])/(2.0*k*delta) + (V[i+k][j+k] - V[i+k][j])/(2.0*k*delta), 2.0));
    }
  }
  return s;
}

void rownaniePoissona(int k)
{
  double V[nx+1][ny+1];
  double s = 0.0, s1 = 0.0;
  int iteracja = 0;
  FILE* plikCalka; FILE* plikMapa;

  for(int i = 0; i < nx+1; i++)
  {
    for(int j = 0; j < ny+1; j++)
    {
      V[i][j] = 0.0;
    }
  }

  for(int i = 0; i < nx+1; i++)
  {
    V[i][0] = sin(2.0*M_PI*delta*i/xMax);
    V[i][ny] = -1.0*sin(2.0*M_PI*delta*i/xMax);
  }

  for(int j = 0; j < ny+1; j++)
  {
    V[0][j] = sin(M_PI*delta*j/yMax);
    V[nx][j] = sin(M_PI*delta*j/yMax);
  }

  while(k >= 1)
  {
    if(k == 16)
    {
      plikCalka = fopen("wyniki/calka_16.txt", "w");
      plikMapa = fopen("wyniki/mapa_16.txt", "w");
    }
    else if(k == 8)
    {
      plikCalka = fopen("wyniki/calka_8.txt", "w");
      plikMapa = fopen("wyniki/mapa_8.txt", "w");
    }
    else if(k == 4)
    {
      plikCalka = fopen("wyniki/calka_4.txt", "w");
      plikMapa = fopen("wyniki/mapa_4.txt", "w");
    }
    else if(k == 2)
    {
      plikCalka = fopen("wyniki/calka_2.txt", "w");
      plikMapa = fopen("wyniki/mapa_2.txt", "w");
    }
    else if(k == 1)
    {
      plikCalka = fopen("wyniki/calka_1.txt", "w");
      plikMapa = fopen("wyniki/mapa_1.txt", "w");
    }
    
    s1 = calka(V, k);
    do{
      for(int i = k; i <= nx-k; i += k)
      {
        for(int j = k; j <= ny-k; j += k)
        {
          V[i][j] = 0.25*(V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);
        }
      }

      s = s1;
      s1 = calka(V, k);
      fprintf(plikCalka, "%d %f \n", iteracja, s);
      iteracja++;
    }while(fabs((s1-s)/s) >= TOL);



    for(int i = 0; i < nx+1; i += k)
    {
      for(int j = 0; j < ny+1; j += k)
      {
        fprintf(plikMapa, "%f %f %f \n", delta*i, delta*j, V[i][j]);
      }
      fprintf(plikMapa, "\n");
    }


    for(int i = 0; i <= nx-k; i += k)
    {
      for(int j = 0; j <= ny-k; j += k)
      {
        V[i + k/2][j + k/2] = 0.25*(V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);
        if(i != 0)
          V[i][j + k/2] = 0.5*(V[i][j] + V[i][j+k]);
        else if(i != (nx-k))
          V[i + k][j + k/2] = 0.5*(V[i+k][j] + V[i+k][j+k]);
        
        if(j != 0)
          V[i + k/2][j] = 0.5*(V[i][j] + V[i+k][j]);
        else if(j != (ny-k))
          V[i + k/2][j + k] = 0.5*(V[i][j+k] + V[i+k][j+k]);
        
      }
    }
  k = k/2;
  fclose(plikCalka);
  fclose(plikMapa);
  }
}

int main(void)
{
  int k = 16;
  rownaniePoissona(k);
  return 0;
}