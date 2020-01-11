#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include "mgmres.h"

const double delta = 0.1;

int j(int nx, int L){
  return (L/(nx+1));
}

int i(int nx, int L){
  return L-j(nx,L)*(nx+1);
}


double El(int nx, int L, int epsilon1, int epsilon2){
  if(i(nx,L) <= nx/2)
    return epsilon1;
  else
    return epsilon2;
}

double ro1(double x, double y, double xmax, double ymax, double sigma)
{
  return exp( -1.0*pow(x-0.25*xmax, 2.0)/pow(sigma,2.0) - pow(y-0.5*ymax, 2.0)/pow(sigma, 2.0) );
}

double ro2(double x, double y, double xmax, double ymax, double sigma)
{
  return -1.0*exp( -1.0*pow(x-0.75*xmax, 2.0)/pow(sigma,2.0) - pow(y-0.5*ymax, 2.0)/pow(sigma, 2.0) );
}



void metodaPoissonaAlgebraiczna(const char* zad = "")
{
  bool saveTestMatrix = false; // testowanie generowania macierzy
  bool saveVector = false;
  bool saveMap = true;
  bool debug = true;
  
  FILE* mapa, *wektor;

  printf("Metoda Poissona %g", zad);

  int V1 = 10;
  int V3 = 10;
  int V2 = -10;
  int V4 = -10;
  double sigma = 0.0;
  double xmax = 0.0;
  double ymax = 0.0;
  int nx = 4;
  int ny = 4;
  int epsilon1 = 1;
  int epsilon2 = 1;

  if(strcmp("macierzTestowa", zad) == 0)
  {
    if(debug) printf("\nZadanie macierz testowa\n");
    saveTestMatrix = true; // testowanie generowania macierzy
    saveVector = true;
    saveMap = false;
    wektor = fopen("wektor.txt", "w");
  }
  if(strcmp("zad1a", zad) == 0)
  {
    if(debug) printf("\nZadanie 1a\n");
    nx = ny= 50;
    mapa = fopen("mapa_zad1a.txt", "w");
  }
  if(strcmp("zad1b", zad) == 0)
  {
    if(debug) printf("\nZadanie 1b\n");
    nx = ny= 100;
    mapa = fopen("mapa_zad1b.txt", "w");
  }
  if(strcmp("zad1c", zad) == 0)
  {
    if(debug) printf("\nZadanie 1c\n");
    nx = ny= 200; 
    mapa = fopen("mapa_zad1c.txt", "w");
  }
  if(strcmp("zad2a", zad) == 0)
  {
    if(debug) printf("\nZadanie 2a\n");
    V1 = V2 = V3 = V4 = 0;
    nx = ny = 100;
    xmax = delta*nx;
    ymax = delta*ny;
    sigma = xmax/10.0;
    epsilon1 = 1;
    epsilon2 = 1;
    mapa = fopen("mapa_zad2a.txt", "w");
  }
  if(strcmp("zad2b", zad) == 0)
  {
    if(debug) printf("\nZadanie 2b\n");
    V1 = V2 = V3 = V4 = 0;
    nx = ny = 100;
    xmax = delta*nx;
    ymax = delta*ny;
    sigma = xmax/10.0;
    epsilon1 = 1;
    epsilon2 = 2;
    mapa = fopen("mapa_zad2b.txt", "w");
  }
  if(strcmp("zad2c", zad) == 0)
  {
    if(debug) printf("\nZadanie 2c\n");
    V1 = V2 = V3 = V4 = 0;
    nx = ny = 100;
    xmax = delta*nx;
    ymax = delta*ny;
    sigma = xmax/10.0;
    epsilon1 = 1;
    epsilon2 = 10;
    mapa = fopen("mapa_zad2c.txt", "w");
  }
  
  if(debug) printf("Parametry: \nV1 %d\tV2 %d\tV3 %d\tV4 %d\nnx %d\tny %d\nxmax %f\tymax %f\tsigma %f\nepsilon1 %d\tepsilon2 %d\n\n", V1,V2,V3,V4,nx,ny,xmax,ymax,sigma,epsilon1,epsilon2);

  int N = (nx+1)*(ny+1); 
  
  double a[5*N];
  int ja[5*N];
  int ia[N+1];
  double b[N];
  double V[N];

  for(int i = 0; i < N+1; i++)
  {
    ia[i] = -1;
  }

  if(debug) printf("Dirichlet start\n");
// wypelnianie Dirichleta
  int k = -1; // numeruje niezerowe elementy A
  int nz_num = 0; // liczba niezerowych elementow
  for(int L = 0; L < N; L++)
  {
    int brzeg = 0;    // wskaznik polozenia 0-srodek obszaru, 1-brzeg
    double vb = 0.0; // potencjal na brzegu

    if(i(nx, L) == 0) // lewy brzeg
    {
      brzeg = 1;
      vb = V1;
    } 
    else if(i(nx, L) == nx) // gorny brzeg
    {
      brzeg = 1;
      vb = V3;
    }
		else if(j(nx, L) == ny) // prawy brzeg
    {
			brzeg = 1;
			vb = V2;
		}
		else if(j(nx, L) == 0) // dolny brzeg
    {
			brzeg = 1;
			vb = V4;
		}
    

    // wypelniamy od razu wektor wyrazow wolnych
    b[L] = -1.0*(ro1(delta*i(nx,L),delta*j(nx,L),xmax,ymax,sigma) + ro2(delta*i(nx,L),delta*j(nx,L),xmax,ymax,sigma));

    if(brzeg == 1)
    {
      b[L] = vb; // wymuszamy wartosc potencjalu na brzegu
    }

    // wypelniamy elementy macierzy a
    ia[L] = -1; // wskaznik dla pierwszego elementu w wierszu

    // lewa skrajna przekatna
    if(L - nx - 1 > 0 && brzeg == 0)
    {
      k++;
      if(ia[L] < 0)
        ia[L] = k;
      
      a[k] = El(nx, L, epsilon1, epsilon2)/pow(delta,2.0);
      ja[k] = L-nx-1;
    }

    // poddiagonala
    if(L - 1 > 0 && brzeg == 0)
    {
      k++;
      if(ia[L] < 0)
        ia[L] = k;
      
      a[k] = El(nx, L, epsilon1, epsilon2)/pow(delta,2.0);
      ja[k] = L-1;
    }

    // diagonala
    k++;
    if(ia[L] < 0)
      ia[L] = k;

    if(brzeg == 0)
      a[k] = -1.0*(2*El(nx, L, epsilon1, epsilon2) + El(nx, L+1, epsilon1, epsilon2) + El(nx, L+nx+1, epsilon1, epsilon2))/pow(delta,2.0);
		else
			a[k] = 1;

    ja[k] = L;

    // naddiagonala
    if(L < N && brzeg == 0)
    {
      k++;
      a[k] = El(nx, L+1, epsilon1, epsilon2)/pow(delta, 2.0);
      ja[k] = L + 1;
    }

    //prawa skrajna przekatna
		if(L < N-nx-1 && brzeg == 0){
			k++;
			a[k] = El(nx,L+nx+1,epsilon1,epsilon2)/pow(delta, 2.0);
			ja[k] = L + nx + 1;
		}
    if(saveVector)
    {
        if(L%5 == 0 && L != 0)
          fprintf(wektor, "\n");
        fprintf(wektor,"%d %d %d %f \n", L, i(nx,L), j(nx,L), b[L]);
    }
  }
  nz_num = k+1;
	ia[N] = nz_num;
  
  if(saveVector)
    fclose(wektor);
  if(debug) printf("Dirichlet stop\n");

  if(saveTestMatrix)
  {
    FILE* testMacierz = fopen("testowa_macierz.txt", "w");
    for(int it = 0; it < 5*N; it++)
      fprintf(testMacierz, "%d %d %d %f\n", it, i(nx,it), j(nx,it), a[it]);
    fclose(testMacierz);
    return;
  }

  int itr_max = 500;
  int mr = 500;
  double tol_abs = pow(10.0, -8.0);
  double tol_rel = pow(10.0, -8.0);
  if(debug) printf("pmgmres_ilu_cr start\n");
  pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);
  if(debug) printf("pmgmres_ilu_cr stop\n");

  //zapis mapy
	if(saveMap){
		double tmp = 0.0;
		for(int it = 0; it < N; it++){
			if(delta*j(nx,it) > tmp)
				fprintf(mapa,"\n");
			fprintf(mapa,"%f %f %f \n", delta*j(nx,it), delta*i(nx,it), V[it]);
			tmp = delta*j(nx,it);
		}
    fclose(mapa);
	}

}

int main(void)
{
    
  metodaPoissonaAlgebraiczna("macierzTestowa");

  metodaPoissonaAlgebraiczna("zad1a");
  metodaPoissonaAlgebraiczna("zad1b");
  metodaPoissonaAlgebraiczna("zad1c");
  
  metodaPoissonaAlgebraiczna("zad2a");
  metodaPoissonaAlgebraiczna("zad2b");
  metodaPoissonaAlgebraiczna("zad2c");
    
  return 0;
}
