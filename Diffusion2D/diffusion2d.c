#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define nx 100
#define ny 100

static double* linspace(double start, double stop, int N);
static void setAll(int Nx, int Ny, double *array, double value);
static void printArray(int N, double *array);
static void *xcalloc(size_t count, size_t size);
static void makeTopHat(int N, double *x, double *u);
static void saveArray(int N, double *array, char *fname);
static void saveTwoArrays(int N, double *array1, double *array2, char *fname);
static void setEqual(int Nx, int Ny, double *array1, double *array2);
static double norm(int Nx, int Ny, double *array);

static void setElem(double *array, int i, int j, double val);
static double getElem(double *array, int i, int j);
static void print2dArray(int Nx, int Ny, double *array);
static void save2dArray(int Nx, int Ny, double *array, char *fname);
static void make2dTopHat(int Nx, int Ny, double *x, double *y, double *u);

int main(int argc, char const *argv[])
{
    double *x   = linspace(0.0, 5.0, nx);
    double *y   = linspace(0.0, 5.0, ny);
    double dx   = fabs(x[1] - x[0]); // Delta x
    double dy   = fabs(y[1] - y[0]); // Delta y
    double endT = 10.0;
    double sigma = 0.25;
    double nu   = 0.05; // Diffusion coefficient
    double dt   = sigma*dx*dy/nu; // Time increment
    double c_x  = 1.0; // Wavespeed set to 1.0
    double c_y  = 1.0; // Wavespeed set to 1.0

    printf("nx = %i\n", nx);
    printf("dx = %lf\n", dx);
    printf("endT = %lf\n", endT);
    printf("dt = %lf\n", dt);
    printf("cx = %lf\n", c_x);
    printf("cy = %lf\n", c_y);

    double *u       = (double*)xcalloc(nx * ny, sizeof(double));
    double *u_prev  = (double*)xcalloc(nx * ny, sizeof(double));
    setAll(nx, ny, u, 0.0);
    make2dTopHat(nx, ny, x, y, u);
    setEqual(nx, ny, u_prev, u);

    saveTwoArrays(nx, x, y, "coords.txt");
    save2dArray(nx, ny, u, "initial.txt");

    // Time loop
    double tmp = 0.0;

	for(double time = 0; time <= endT; time += dt)
    {
        printf("t = %lf\t|u| = %lf\n", time, norm(nx, ny, u));
        for(int i = 0; i < nx; i++)
        {
            for(int j = 0; j < ny; j++)
            {
                if(i == 0 || j == 0 || i == nx-1 || j == ny-1 )
                {
                    tmp = 0.0;
                }else{
                    tmp = getElem(u_prev, i, j) + ( (nu*dt/(dx*dx)) * ( getElem(u, i+1, j) - 2.0*getElem(u, i, j) + getElem(u, i-1, j) ) )
                                        + ( (nu*dt/(dy*dy)) * ( getElem(u, i, j+1) - 2.0*getElem(u, i, j) + getElem(u, i, j-1) ) );
                }
                setElem(u, i, j, tmp);
            }
        }
        setEqual(nx, ny, u_prev, u);
    }

    save2dArray(nx, ny, u_prev, "final.txt");


    free(x);
    free(y);
    free(u);
    free(u_prev);

    return 0;
}
static void setElem(double *array, int i, int j, double val){
    array[nx * i + j] = val;
}
static double getElem(double *array, int i, int j){
    return array[nx * i + j];
}
static double norm(int Nx, int Ny, double *array){
    double n = 0.0;
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            n += getElem(array, i, j) * getElem(array, i, j);
        }

    }
    return sqrt(n);
}
// Sets array1 = array2
static void setEqual(int Nx, int Ny, double *array1, double *array2){
    double tmp;
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            tmp = getElem(array2, i, j);
            setElem(array1, i, j, tmp);
        }
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
// Make a top hat function from an array, in 2D
static void make2dTopHat(int Nx, int Ny, double *x, double *y, double *u){
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            if(x[i] > 2.0 && x[i] < 3.0 && y[j] > 2.0 && y[j] < 3.0){
                setElem(u, i, j, 2.0);
            }
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
static void save2dArray(int Nx, int Ny, double *array, char *fname){
    FILE *fp = fopen(fname, "w");
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            fprintf(fp, "%lf\t", getElem(array, i, j));
        }
        fprintf(fp, "\n");
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
static void print2dArray(int Nx, int Ny, double *array){
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            printf("%lf\t", getElem(array, i, j));
        }
        printf("\n");
    }
}
static void setAll(int Nx, int Ny, double *array, double value)
{
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            setElem(array, i, j, value);
        }
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
