#include <cstdio>
#include <cstdlib>
#include <cmath>

const double beta = 0.001;
const int N = 500;
const double my_gamma = 0.1;
const double dt = 0.1;
const double tMax = 100.0;
const double miMax = 20;

// zadanie 1
double alfa(void)
{
    return beta*N - my_gamma;
}
double helpfulFunction(double u)
{
	return alfa()*u - beta*u*u;
}

void metodaPicarda(FILE* plik, double u0, double TOL)
{
	double un = u0;
	double un1 = un;
	double umi = 0.0;
	int mi_iteracja = 0;

	fprintf(plik, "%.1f %f %f \n" , 0.0, un1, N-un1);

	for(double t = 0.1; t <= tMax; t += dt)
	{
		un = un1;
		mi_iteracja = 0;
		while(mi_iteracja < miMax)
		{
			mi_iteracja++;
			umi = un1;
			un1 = un + (dt/2.0)*(helpfulFunction(un) + helpfulFunction(umi));
		}
		fprintf(plik, "%.1f %f %f \n", t, un1, N-un1);

	}
}


// zadanie 2
void metodaIteracjaNewtona(FILE* plik, double u0, double TOL)
{
	double un = u0;
	double un1 = un;
	double umi = 0.0;
	int mi_iteracja = 0;

	fprintf(plik, "%.1f %f %f \n" , 0.0, un1, N-un1);

	for(double t = 0.1; t <= tMax; t += dt)
	{
		un = un1;
		mi_iteracja = 0;
		while(fabs(umi - un) >= TOL && mi_iteracja < miMax)
		{
			mi_iteracja++;
			umi = un1;
			double licznik = umi - un - (dt/2.0)*(helpfulFunction(un) + helpfulFunction(umi));
			double mianownik = 1.0 - (dt/2.0)*(alfa() - 2.0*beta*umi);
			un1 = umi - licznik/mianownik;
		
		}
		fprintf(plik, "%.1f %f %f \n", t, un1, N-un1);

	}
}



// zadanie 3

double m11(double a11, double U1)
{
	return 1.0 - dt*a11*(alfa() - 2*beta*U1);
}

double m12(double a12, double U2)
{
	return -1.0*dt*a12*(alfa() - 2*beta*U2);
}

double m21(double a21, double U1)
{
	return -1.0*dt*a21*(alfa() - 2*beta*U1);
}

double m22(double a22, double U2)
{
	return 1.0 - dt*a22*(alfa() - 2*beta*U2);
}

double F1(double un, double U1, double U2, double a11, double a12)
{
    return U1 - un - dt*(a11*helpfulFunction(U1) + a12*helpfulFunction(U2));
}

double F2(double un, double U1, double U2, double a21, double a22)
{
    return U2 - un - dt*(a21*helpfulFunction(U1) + a22*helpfulFunction(U2));
}

void metodaRK2(FILE* plik, double u0, double TOL)
{
	double a11 = 0.25;
	double a12 = 0.25 - sqrt(3.0)/6.0;
	double a21 = 0.25 + sqrt(3.0)/6.0;
	double a22 = 0.25;
	double b, b1, b2; b = b1 = b2 = 0.5; // b1 = b2 wiec oznaczamy to jako b
    double c1 = 0.5 - sqrt(3.0)/6.0;
    double c2 = 0.5 + sqrt(3.0)/6.0;
	
	double un = u0;
	double un1 = un;
	double umi = 0.0;
	int mi_iteracja = 0;

	double U1 = 0.0;
	double U2 = 0.0;
	double U1mi = 0.0;
	double U2mi = 0.0;

	fprintf(plik, "%.1f %f %f \n" , 0.0, un1, N-un1);

	for(double t = 0.1; t <= tMax; t += dt)
	{
		mi_iteracja = 0;
		un = un1;
		U1mi = 0.;
		U2mi = 0.;
		U1 = un;
		U2 = un;

		while((fabs(U1 - U1mi) >= TOL || fabs(U2 - U2mi) >= TOL) && mi_iteracja < miMax  )
        {
            mi_iteracja++;			
            U1mi = U1;
			U2mi = U2;

			double dU1 = (F2(un,U1,U2,a21,a22) * m12(a12,U2) - F1(un,U1,U2,a11,a12) * m22(a22,U2))/
																					(m11(a11, U1)*m22(a22,U2) - m12(a12,U2) * m21(a21,U1));

			double dU2 = (F1(un,U1,U2,a11,a12) * m21(a21, U1) - F2(un, U1, U2, a21, a22) * m11(a11,U1)) / (m11(a11,U1)*m22(a22,U2) - m12(a12,U2)*m21(a21,U1));

			U1 = U1mi+dU1;
			U2 = U2mi + dU2;

		}
		un1 = un + dt*(b*(alfa()*U1 - beta*U1*U1) + b*(alfa()*U2 - beta*U2*U2));
		fprintf(plik, "%.1f %f %f \n", t, un1, N-un1);

	}
}


int main(void)
{
	FILE* plik;
	double u0 = 1.0;
	double TOL = pow(10.0, -6.0);

	// zadanie 1
	plik = fopen("wynik_zad_1.txt", "w");
	metodaPicarda(plik, u0, TOL);
	fclose(plik);


	// zadanie 2
	plik = fopen("wynik_zad_2.txt", "w");
	metodaIteracjaNewtona(plik, u0, TOL);
	fclose(plik);


	// zadanie 3
	plik = fopen("wynik_zad_3.txt", "w");
	metodaRK2(plik, u0, TOL);
	fclose(plik);
}
