#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define nx 31
#define ny 31
#define tolerence 1e-4

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
static double L1norm(int Nx, int Ny, double *array1, double *array2);

int main(int argc, char const *argv[])
{
    double *x   = linspace(0.0, 2.0, nx);
    double *y   = linspace(0.0, 1.0, ny);
    double dx   = fabs(x[1] - x[0]); // Delta x
    double dy   = fabs(y[1] - y[0]); // Delta y

    double c_x  = 1.0; // Wavespeed set to 1.0
    double c_y  = 1.0; // Wavespeed set to 1.0

    printf("nx = %i\n", nx);
    printf("dx = %lf\n", dx);
    printf("cx = %lf\n", c_x);
    printf("cy = %lf\n", c_y);

    double *p       = (double*)xcalloc(nx * ny, sizeof(double));
    double *p_prev       = (double*)xcalloc(nx * ny, sizeof(double));
    setAll(nx, ny, p, 0.0);

    for(int j = 0; j < ny; j++){
        double bVal = y[j];
        setElem(p, nx-1, j, bVal);
    }

    setEqual(nx, ny, p_prev, p);
    saveTwoArrays(nx, x, y, "coords.txt");
    save2dArray(nx, ny, p, "initial.txt");

    // Time loop
    double tmp = 0.0;
    int iter = 0;
    double error = L1norm(nx, ny, p, p_prev);
	while(1)
    {
        for(int i = 0; i < nx; i++)
        {
            for(int j = 0; j < ny; j++)
            {
                if(i == 0){
                    tmp = 0.0;
                }else if(i == nx-1){
                    tmp = y[j];
                }else if(j == 0){
                    tmp = getElem(p_prev, i, 1);
                }else if(j == ny-1){
                    tmp = getElem(p_prev, i, ny-2);
                }else{
                    tmp = (dy*dy)*( getElem(p_prev, i, j+1) + getElem(p_prev, i, j-1) ) + (dx*dx)*( getElem(p_prev, i+1, j) + getElem(p_prev, i-1, j) );
                    tmp /= 2.0*(dx*dx + dy*dy);
                }
                setElem(p, i, j, tmp);
            }
        }
        error = L1norm(nx, ny, p, p_prev);
        printf("Iter = %i, L1norm(p-p_prev) = %lf\n", iter, error);
        if(error < tolerence){ break; }
        setEqual(nx, ny, p_prev, p);
        iter++;
    }

    save2dArray(nx, ny, p, "final.txt");

    free(x);
    free(y);
    free(p);

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
static double L1norm(int Nx, int Ny, double *array1, double *array2){
    double n = 0.0;
    double sum = 0.0;
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            n += fabs(getElem(array1, i, j) - getElem(array2, i, j));
            sum += fabs(getElem(array1, i, j));
        }
    }
    return n/sum;
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
