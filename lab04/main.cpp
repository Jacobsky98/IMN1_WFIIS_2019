#include <cstdio>
#include <cstdlib>
#include <cmath>

const int nx = 150.0;
const int ny = 100.0;
const int maxIt = 1000;

const double epsilon = 1.0;
const double delta = 0.1;
const double V1 = 10.0;
const double V2 = 0.0;
const double xmax = delta*nx;
const double ymax = delta*ny;
const double sigmaX = 0.1*xmax;
const double sigmaY = 0.1*ymax;
const double TOL = pow(10.0, -8.0);


double ro1(double x, double y)
{
  return exp(-1.0*( pow(x-0.35*xmax, 2.0)/pow(sigmaX, 2.0) ) - ( pow(y-0.5*ymax, 2.0)/pow(sigmaY, 2.0) ));
}

double ro2(double x, double y)
{
  return -1.0*exp(-1.0*( pow(x-0.65*xmax, 2.0)/pow(sigmaX, 2.0) ) - ( pow(y-0.5*ymax, 2.0)/pow(sigmaY, 2.0) ));
}

double ro(double x, double y)
{
  return ro1(x,y) + ro2(x,y);
}

void relaksacjaGlobalna(double wG, FILE* plikS, FILE* plikMapa, FILE* plikBlad)
{
  double V[nx+1][ny+1];
  double Vs[nx+1][ny+1];

  for(int i = 0 ; i <= nx; i++)
  {
    for(int j = 0; j <= ny; j++)
    {
      V[i][j] = 0.0;
      Vs[i][j] = 0.0;
    }
  }

  for(int i = 0; i <= ny; i++)
  {
    V[0][i] = V2;
    V[nx][i] = V1;
    Vs[0][i] = V2;
    Vs[nx][i] = V1;
  }



  double sit = 100.0, sit_1 = 0.0;
	int it = 1;
	do {

		for(int i = 1; i < nx; i++) 
    {
			for(int j = 1; j < ny; j++)
      {
				V[i][j] = (1.0/4.0)*(Vs[i+1][j] + Vs[i-1][j] + Vs[i][j+1] + Vs[i][j-1] +
					(pow(delta, 2.0)*ro(i, j)/epsilon));
			}
		}

		for(int j = 1; j < ny; j++) 
    {
			V[0][j] = V[1][j];
			V[nx][j] = V[nx-1][j];
		}

		for(int i = 0; i <= nx; i++) 
    {
			for(int j = 0; j< ny; j++)
      {
				Vs[i][j] = (1.0 - wG)*Vs[i][j] + wG*V[i][j];
			}
		}

		sit_1 = sit;
		sit = 0.0;
		for(int i = 0; i < nx; i++) 
    {
			for(int j = 0; j < ny; j++) 
      {
				sit += pow(delta, 2.0)*(pow((V[i+1][j] - V[i][j])/delta, 2.0)/2.0	+ pow((V[i][j+1] - V[i][j])/delta, 2.0)/2.0 - ro(i, j)*V[i][j]);
			}
		}
    fprintf(plikS, "%d, %f\n", it, sit);
		it++;
	} while((fabs((sit - sit_1)/sit_1) >= TOL) && (it < maxIt));

    for (int i = 0; i < nx+1; i++) 
    {
      for (int j = 0; j < ny+1; j++) 
      {
        // if (i == -1 && j == -1) { fprintf(plikMapa, "0 "); fprintf(plikBlad, "0 "); }
        // else {
          fprintf(plikMapa, "%f ", V[i][j]);
          fprintf(plikBlad, "%f ", pow(delta, 2.0) * V[i][j] + ro(i,j));
        //}
		}
		fprintf(plikMapa, "\n");
		fprintf(plikBlad, "\n");


}
}

void relaksacjaLokalna(double wL, FILE* plik)
{
  double V[nx+1][ny+1];

  for(int i = 0 ; i <= nx; i++)
  {
    for(int j = 0; j <= ny; j++)
    {
      V[i][j] = 0.;
    }
  }

  for(int i = 0; i <= ny; i++)
  {
    V[0][i] = V2;
    V[nx][i] = V1;
  }

  	double sit = 1, sit_1;
	int it = 1;
	do {

		for(int i = 1; i < nx-1; i++) {
			for(int j = 1; j < ny-1; j++){
				V[i][j] = (1.0 - wL)*V[i][j] + wL*(V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + (pow(delta, 2.0)*ro(delta*i, delta*j))/epsilon)/4.0;
			}
		}

		for(int j = 1; j < ny; j++) {
			V[0][j] = V[1][j];
			V[nx-1][j] = V[nx-2][j];
		}

		sit_1 = sit;
		sit = 0.0;
		for(int i = 0; i < nx-1; i++) {
			for(int j = 0; j < ny-1; j++) {
				sit += pow(delta, 2.0)*(pow((V[i+1][j] - V[i][j])/delta, 2.0)/2.0 + pow((V[i][j+1] - V[i][j])/delta, 2.0)/2.0 - ro(delta*i, delta*j)*V[i][j]);
			}
		}

    fprintf(plik, "%d, %f\n", it, sit);

		it++;
	} while((fabs((sit - sit_1)/sit_1) > TOL) &&( it < maxIt));
}



int main(void)
{
  FILE* plik; FILE* plikMapa; FILE* plikBlad;

  plik = fopen("wyniki/relaksacja_globalna_0.6.txt", "w");
  plikMapa = fopen("wyniki/mapa_globalna_0.6.txt", "w");
  plikBlad = fopen("wyniki/mapa_globalna_bledu_0.6.txt", "w");
  relaksacjaGlobalna(0.6, plik, plikMapa, plikBlad);
  fclose(plik);
  fclose(plikBlad);
  fclose(plikMapa);

  plik = fopen("wyniki/relaksacja_globalna_1.0.txt", "w");
  plikMapa = fopen("wyniki/mapa_globalna_1.0.txt", "w");
  plikBlad = fopen("wyniki/mapa_globalna_bledu_1.0.txt", "w");
  relaksacjaGlobalna(1.0, plik, plikMapa, plikBlad);
  fclose(plik);
  fclose(plikBlad);
  fclose(plikMapa);

        plik = fopen("wyniki/relaksacja_lokalna_1.0.txt", "w");
        relaksacjaLokalna(1.0, plik);
        fclose(plik);

        plik = fopen("wyniki/relaksacja_lokalna_1.4.txt", "w");
        relaksacjaLokalna(1.4, plik);
        fclose(plik);

        plik = fopen("wyniki/relaksacja_lokalna_1.8.txt", "w");
        relaksacjaLokalna(1.8, plik);
        fclose(plik);
        
        plik = fopen("wyniki/relaksacja_lokalna_1.9.txt", "w");
        relaksacjaLokalna(1.9, plik);
        fclose(plik);
}
