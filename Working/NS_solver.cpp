// A very simple Navier-Stokes solver for a drop falling in a rectangular domain
// The viscosity is taken to be a constant
// A forward in time, centered in space discretization is used
// The density is advected by a simple upwind scheme

#include <iostream> 
using namespace std;

// Domain size and physical variables
float Lx = 1;
float Ly = 1;
int gx = 0;
int gy = -100;
int rho1 = 1;
int rho2 = 2;
double m0 = 0.01;
int rro = rho1;
int unorth = 0;
int usouth = 0;
int veast = 0;
int vwest = 0;
int time = 0;

// Initial drop size and location
double rad = 0.15;
double xc = 0.5;
double yc = 0.7;

// Numerical variables
static const int nx = 32;
static const int ny = 32;
double dt = 0.00125;
int nstep = 100;
int maxIt = 200;
double maxError = 0.001;
double beta = 1.2;

// Initialise arrays 
double u[nx + 1][ny + 2] = {};
double ut[nx + 1][ny + 2] = {};
double uu[nx + 1][ny + 1] = {};
double v[nx + 2][ny + 1] = {};
double vt[nx + 2][ny + 1] = {};
double vv[nx + 1][ny + 1] = {};
double p[nx + 2][ny + 2] = {};
double tmp1[nx + 2][ny + 2] = {};

// Set the grid
double dx = (Lx / nx);
double dy = (Ly / ny);
double x[nx + 2] = {};
double y[ny + 2] = {};
int main() {

	for (int i = 0; i < nx + 2; i++) {
		x[i] = double(-0.5+i)*dx;
 	}

	for (int i = 0; i < ny + 2; i++) {
		y[i] = double(-0.5+i)*dy;
	}

	// Set the density 
	double r[nx + 2][ny + 2] = {}; 

	for (int i = 0; i < nx + 2; i++) {
		for (int j = 0; j < ny + 2; j++) {
			if ((pow((x[i] - xc),2) + pow((y[j] - yc),2)) < pow(rad,2)) {
				r[i][j] = rho2;
			}
			else {
				r[i][j] = rho1;
			}
		}
	}

	// Print
	for (int j = 0; j < ny + 2; j++) {
		for (int i = 0; i < nx + 2; i++) {
			cout << r[i][j];
		}
		cout << endl;
	}

	for (int i = 0; i < nx + 2; i++) {
		cout << x[i];
		cout << endl;
	}

	return 0;
}