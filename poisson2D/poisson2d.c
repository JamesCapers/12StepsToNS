#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static void *xcalloc(size_t count, size_t size);
static void print2d(double **array, int nx, int ny);
static void save2d(double **array, int nx, int ny, char *fname);
static double **new2Darray(int nx, int ny);
static void free2dArray(double **array, int nx, int ny);
static double* linspace(double start, double stop, int N);
static void setEqual2D(int nx, int ny, double **a1, double **a2);
static void saveTwoArrays(int N, double *array1, double *array2, char *fname);

int main(int argc, char **argv)
{
    int nx      = 50;
    int ny      = 50;
    double *x   = linspace(0.0, 5.0, nx);
    double *y   = linspace(0.0, 5.0, ny);
    double dx   = fabs(x[1] - x[0]); // Delta x
    double dy   = fabs(y[1] - y[0]); // Delta y
    double endT = 10.0;
    double dt   = 0.01;

    // Allocate the arrays
    double **p           = new2Darray(nx, ny);
    double **p_prev      = new2Darray(nx, ny);
    double **source      = new2Darray(nx, ny);

    // Setup the source
    source[nx/4][ny/4]      = 100.0;
    source[3*nx/4][3*ny/4]  = -100.0;

    setEqual2D(nx, ny, p_prev, source);

    saveTwoArrays(nx, x, y, "coords.txt");
    save2d(p_prev, nx, ny, "initial.txt");

    for(double t = 0; t <= endT; t += dt){
        for(int i = 0; i < nx; i++){
            for(int j = 0; j < ny; j++){

                if(i == 0 || i == nx-1 || j == 0 || j == ny-1 ){
                    p[i][j] = 0.0;
                }else{
                    p[i][j] = (dy*dy)*( p_prev[i+1][j] + p_prev[i-1][j] ) + (dx*dx)*(p_prev[i][j+1] + p_prev[i][j-1]) - source[i][j] * (dx*dx*dy*dy);
                    p[i][j] *= (1.0)/(2.0*(dx*dx + dy*dy));
                }
            }
        }

        setEqual2D(nx, ny, p_prev, p);
    }

    save2d(p, nx, ny, "final.txt");

    free(x);
    free(y);
    free2dArray(p, nx, ny);
    free2dArray(p_prev, nx, ny);
    free2dArray(source, nx, ny);

    return 0;
}

// Set a1 = a2
static void setEqual2D(int nx, int ny, double **a1, double **a2){
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            a1[i][j] = a2[i][j];
        }
    }
}

// Safe way to allocate memory and set all blocks to 0
static void *xcalloc(size_t count, size_t size)
{
	void *ptr;
	ptr = calloc(count, size);

	if(ptr == NULL){
		fprintf(stderr, "Could not allocate memory\n");
		exit(1);
	}else{
		return ptr;
	}
}

/* FUNCTIONS TO HELP WITH 2D/3D ARRAYS */

// Print 2d array to terminal
static void print2d(double **array, int nx, int ny){
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            printf("%lf\t", array[i][j]);
        }
        printf("\n");
    }
}

static void save2d(double **array, int nx, int ny, char *fname){
    FILE *fp = fopen(fname, "w");

    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            fprintf(fp, "%lf\t", array[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

// Allocate new 2d array, setting all elements to zero
// Access elements using array[i][j]
static double **new2Darray(int nx, int ny){
    double **array;
    array = (double**)xcalloc(nx, sizeof(double*));
    for(int i = 0; i < nx; i++){
        array[i] = (double*)xcalloc(ny, sizeof(double));
    }
    return array;
}

// Free a 2d array
static void free2dArray(double **array, int nx, int ny){
    for(int i = 0; i < nx; i++){
        free(array[i]);
    }
    free(array);
}

static double* linspace(double start, double stop, int N)
{
    double *array = (double*)malloc(sizeof(double) * N);
    if(array != NULL){
        double step = (stop - start) / (N - 1.0);
        for(int i = 0; i < N; i++)
        {
            array[i] = start + step*i;
        }

        return array;
    }else{
        fprintf(stderr, "Could not allocate memory\n");
		exit(1);
    }
}

// Arrays must be the same size - we do not check for this
static void saveTwoArrays(int N, double *array1, double *array2, char *fname){
    FILE *fp = fopen(fname, "w");
    for(int i = 0; i < N; i++){
        fprintf(fp, "%lf\t%lf\n", array1[i], array2[i]);
    }
    fclose(fp);
}
