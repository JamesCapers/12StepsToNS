#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static double* linspace(double start, double stop, int N);
static void setAll(int N, double *array, double value);
static void printArray(int N, double *array);
static void *xcalloc(size_t count, size_t size);
static void makeTopHat(int N, double *x, double *u);
static void saveArray(int N, double *array, char *fname);
static void saveTwoArrays(int N, double *array1, double *array2, char *fname);
static void setEqual(int N, double *array1, double *array2);
static double norm(int N, double *array);

int main(int argc, char const *argv[])
{
    int nx      = 100; // Number of x points
    // double *x   = (double*)xcalloc(nx, sizeof(double));
    double *x   = linspace(0.0, 3.0, nx);
    double dx   = fabs(x[1] - x[0]); // Delta x
    int nt      = 50; // Number of time steps
    double dt   = 0.01; // Time increment
    double c    = 2.0; // Wavespeed set to 1.0

    printf("nx = %i\n", nx);
    printf("dx = %lf\n", dx);
    printf("nt = %i\n", nt);
    printf("dt = %lf\n", dt);
    printf("c = %lf\n", c);

    double u[nx], u_prev[nx];
    setAll(nx, u, 1.0);
    makeTopHat(nx, x, u);
    setEqual(nx, u_prev, u);

    saveTwoArrays(nx, x, u, "initial.txt");

    // Step through time
    for(int k = 0; k < nt; k++){
        printf("t = %lf, |u| = %lf\n", k*dt, norm(nx, u));
        for(int i = 1; i < nx; i++){
            u[i] = u_prev[i] - c * (dt/dx) * (u_prev[i] - u_prev[i-1]);
        }
        setEqual(nx, u_prev, u);
    }

    saveTwoArrays(nx, x, u, "final.txt");

    free(x);

    return 0;
}
static double norm(int N, double *array){
    double n = 0.0;
    for(int i = 0; i < N; i++){
        n += array[i]*array[i];
    }
    return sqrt(n);
}
// Sets array1 = array2
static void setEqual(int N, double *array1, double *array2){
    for(int i = 0; i < N; i++){
        array1[i] = array2[i];
    }
}
// Make a top hat function from an array
static void makeTopHat(int N, double *x, double *u){
    for(int i = 0; i < N; i++){
        if(x[i] >= 0.75 && x[i] <= 1.25){
            u[i] = 2.0;
        }
    }
}
static void saveArray(int N, double *array, char *fname){
    FILE *fp = fopen(fname, "w");
    for(int i = 0; i < N; i++){
        fprintf(fp, "%lf\n", array[i]);
    }
    fclose(fp);
}
// Arrays must be the same size - we do not check for this
static void saveTwoArrays(int N, double *array1, double *array2, char *fname){
    FILE *fp = fopen(fname, "w");
    for(int i = 0; i < N; i++){
        fprintf(fp, "%lf\t%lf\n", array1[i], array2[i]);
    }
    fclose(fp);
}
static void printArray(int N, double *array){
    for(int i = 0; i < N; i++){
        printf("%lf\n", array[i]);
    }
}
static void setAll(int N, double *array, double value)
{
    for(int i = 0; i < N; i++){
        array[i] = value;
    }
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
