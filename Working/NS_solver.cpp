// A very simple Navier-Stokes solver for a drop falling in a rectangular domain
// The viscosity is taken to be a constant
// A forward in time, centered in space discretization is used
// The density is advected by a simple upwind scheme

#include <iostream> 
#include <algorithm>
#include <iomanip>
using namespace std;

// Domain size and physical variables
float Lx = 1;
float Ly = 1;
double gx = 0;
double gy = -100;
double rho1 = 1;
double rho2 = 2;
double m0 = 0.01;
double rro = rho1;
double unorth = 0;
double usouth = 0;
double veast = 0;
double vwest = 0;
double timep = 0;

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
double maxAllowedError = 0.001;
double beta = 1.2;

// Initialise arrays 
double u[nx + 1][ny + 2] = {0};
double ut[nx + 1][ny + 2] = {0};
double uu[nx + 1][ny + 1] = {0};
double v[nx + 2][ny + 1] = {0};
double vt[nx + 2][ny + 1] = {0};
double vv[nx + 1][ny + 1] = {0};
double p[nx + 2][ny + 2] = {0};
double tmp1[nx + 2][ny + 2] = {0};
double tmp2[nx + 2][ny + 2] = {0};

// Set the grid
double dx = (Lx / nx);
double dy = (Ly / ny);
double x[nx + 2] = {};
double y[ny + 2] = {};

double findMaxError(double old[][ny + 2], double p[][ny + 2]) {
	double maxError = 0;
	double error = 0;

	for (int i = 0; i < nx + 2; i++) {
		for (int j = 0; j < ny + 2; j++) {
			error = abs(old[i][j] - p[i][j]);
			if (error > maxError) {
				maxError = error;
			}
		}
	}

	return maxError;
}

int main() {

	for (int i = 0; i < nx + 2; i++) {
		x[i] = double(-0.5+i)*dx;
 	}

	for (int i = 0; i < ny + 2; i++) {
		y[i] = double(-0.5+i)*dy;
	}

	// Set the density 
	double r[nx + 2][ny + 2] = {}; 
	double rt[nx + 2][ny + 2] = {};

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
	 
	// Print initial state 
	for (int j = 0; j < ny + 2; j++) {
		for (int i = 0; i < nx + 2; i++) {
			cout << fixed << setprecision(1) << r[i][j];
		}
		cout << endl;
	}

	// Start the time-loop
	for (int t = 0; t < nstep; t++) {

		// Calculate velocity in ghost cells
		for (int j = 0; j < nx + 1; j++) {
			u[j][0] = 2 * usouth - u[j][1];
			u[j][ny + 1] = 2 * unorth - u[j][ny];
		}

		for (int j = 0; j < ny + 1; j++) {
			v[0][j] = 2 * vwest - v[1][j];
			v[nx + 1][j] = 2 * veast - v[nx][j];
		}

		// Calculate temporary velocity field
		for (int j = 1; j < nx; j++) {
			for (int k = 1; k < ny + 1; k++) {
				ut[j][k] = u[j][k] + dt * (0.25 * (
					(1 / dx) * (pow(u[j + 1][k] + u[j][k], 2) - pow(u[j][k] + u[j - 1][k], 2))
					+ (1 / dy) * ((u[j][k + 1] + u[j][k]) * (v[j + 1][k] + v[j][k])
						- (u[j][k] + u[j][k - 1]) * (v[j + 1][k - 1] + v[j][k - 1])))
					+ (m0 / (0.5 * (r[j + 1][k] + r[j][k]))) * (
						(u[j + 1][k] - 2 * u[j][k] + u[j - 1][k]) / pow(dx, 2)
						+ (u[j][k + 1] - 2 * u[j][k] + u[j][k - 1]) / pow(dy, 2)) + gx);
			}
		}

		for (int j = 1; j < nx + 1; j++) {
			for (int k = 1; k < ny; k++) {
				vt[j][k] = v[j][k] + dt * (0.25 * (
					(1 / dx) * ((u[j][k + 1] + u[j][k]) * (v[j + 1][k] + v[j][k])
						- (u[j + 1][k] + u[j + 1][k - 1]) * (v[j][k] + v[j - 1][k]))
					+ (1 / dy) * (pow(v[j][k + 1] + v[j][k], 2) - pow(v[j][k] + v[j][k - 1], 2)))
					+ (m0 / (0.5 * (r[j][k + 1] + r[j][k]))) * (
						(v[j + 1][k] - 2 * v[j][k] + v[j - 1][k]) / pow(dx, 2)
						+ (v[j][k + 1] - 2 * v[j][k] + v[j][k - 1]) / pow(dy, 2)) + gy);
			}
		}

		// Source term
		double lrg = 1000;
		for (int j = 0; j < nx + 2; j++) {
			for (int k = 0; k < ny + 2; k++) {
				if (j == 0 || j == nx + 1 || k == 0 || k == ny + 1) {
					rt[j][k] = lrg;
				}
				else {
					rt[j][k] = r[j][k];
					tmp1[j][k] = (0.5 / dt) * ((ut[j][k] - ut[j - 1][k]) / dx + (vt[j][k] - vt[j][k - 1]) / dy);
					tmp2[j][k] = 1 / ((1 / dx) * (1 / (dx * (rt[j + 1][k] + rt[j][k]))
						+ 1 / (dx * (rt[j - 1][k] + rt[j][k])))
						+ (1 / dy) * (1 / (dy * (rt[j][k + 1] + rt[j][k]))
							+ 1 / (dy * (rt[j][k - 1] + rt[j][k]))));
				}
			}
		}

		double oldArray[nx + 2][ny + 2] = {};

		// Solve for pressure
		for (int j = 0; j < maxIt; j++) {

			for (int k = 0; k < nx + 2; k++) {
				for (int g = 0; g < ny + 2; g++) {
					oldArray[k][g] = p[k][g];
				}
			}
			
			for (int k = 1; k < nx + 1; k++) {
				for (int g = 1; g < ny + 1; g++) {
					p[k][g] = (1 - beta) * p[k][g] + beta * tmp2[k][g] * (
						(1 / dx) * (p[k + 1][g] / (dx * (rt[k + 1][g] + rt[k][g])) +
							p[k - 1][g] / (dx * (rt[k - 1][g] + rt[k][g]))) +
						(1 / dy) * (p[k][g + 1] / (dy * (rt[k][g + 1] + rt[k][g])) +
							p[k][g - 1] / (dy * (rt[k][g - 1] + rt[k][g]))) - tmp1[k][g]
						);
				}
			}

			if (findMaxError(oldArray, p) < maxAllowedError) {
				break;
			}
		}

		// Correct the u-velocity
		for (int i = 1; i < nx; i++) {
			for (int j = 1; j < ny + 1; j++) {
				u[i][j] = ut[i][j] - dt * (2 / dx) * (p[i + 1][j] - p[i][j]) / (r[i + 1][j] + r[i][j]);
			}
		}

		// Correct the v-velocity 
		for (int i = 1; i < nx + 1; i++) {
			for (int j = 1; j < ny; j++) {
				v[i][j] = vt[i][j] - dt * (2 / dy) * (p[i][j + 1] - p[i][j]) / (r[i][j + 1] + r[i][j]);
			}
		}

		// Advect density using centered difference plus diffusion
		double ro[nx + 2][ny + 2] = {};
		for (int i = 0; i < nx + 2; i++) {
			for (int j = 0; j < ny + 2; j++) {
				ro[i][j] = r[i][j];
			}
		}

		for (int i = 1; i < nx + 1; i++) {
			for (int j = 1; j < ny + 1; j++) {
				r[i][j] = ro[i][j] - (0.5 * dt / dx) * (u[i][j] * (ro[i + 1][j] + ro[i][j])
					- u[i - 1][j] * (ro[i - 1][j] + ro[i][j]))
					- (0.5 * dt / dy) * (v[i][j] * (ro[i][j + 1] + ro[i][j])
						- v[i][j - 1] * (ro[i][j - 1] + ro[i][j]))
					+ (m0 * dt / dx / dx) * (ro[i + 1][j] - 2 * ro[i][j] + ro[i - 1][j])
					+ (m0 * dt / dy / dy) * (ro[i][j + 1] - 2 * ro[i][j] + ro[i][j - 1]);
			}
		}
	}

	// Print end state 
	for (int j = 0; j < ny + 2; j++) {
		for (int i = 0; i < nx + 2; i++) {
			cout << fixed << setprecision(1) << r[i][j];
			cout << " ";
		}
		cout << endl;
	}

	return 0;

}