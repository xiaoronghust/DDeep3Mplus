#include "mex.h"
#include "math.h"
#ifdef MAX
#undef MAX
#endif
#define MAX(a, b) ((a)>(b)?(a):(b))
#define n 3

/* This function
 *
 *
 *
 */

/* Eigen decomposition code for symmetric 3x3 matrices, copied from the public
 * 本征分解代码对称3x3矩阵，从public复制
 * domain Java Matrix library JAMA. */
static double hypot2(double x, double y) { return sqrt(x*x+y*y); } 
/* 返回平方和 */

__inline double absd(double val){ if(val>0){ return val;} else { return -val;} };
/* 返回绝对值 */

/* Symmetric Householder reduction to tridiagonal form.
   对称矩阵减少到三对角线的形式。 */
static void tred2(double V[n][n], double d[n], double e[n]) {
    
/*  This is derived from the Algol procedures tred2 by */
/*  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for */
/*  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding */
/*  Fortran subroutine in EISPACK. */
    int i, j, k;
    double scale;
    double f, g, h;
    double hh;
    for (j = 0; j < n; j++) {d[j] = V[n-1][j]; }
    
    /* Householder reduction to tridiagonal form. */
    
    for (i = n-1; i > 0; i--) {
        /* Scale to avoid under/overflow. */
        scale = 0.0;
        h = 0.0;
        for (k = 0; k < i; k++) { scale = scale + fabs(d[k]); }
        if (scale == 0.0) {
            e[i] = d[i-1];
            for (j = 0; j < i; j++) { d[j] = V[i-1][j]; V[i][j] = 0.0;  V[j][i] = 0.0; }
        } else {
            
            /* Generate Householder vector. */
            
            for (k = 0; k < i; k++) { d[k] /= scale; h += d[k] * d[k]; }
            f = d[i-1];
            g = sqrt(h);
            if (f > 0) { g = -g; }
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            for (j = 0; j < i; j++) { e[j] = 0.0; }
            
            /* Apply similarity transformation to remaining columns. */
            
            for (j = 0; j < i; j++) {
                f = d[j];
                V[j][i] = f;
                g = e[j] + V[j][j] * f;
                for (k = j+1; k <= i-1; k++) { g += V[k][j] * d[k]; e[k] += V[k][j] * f; }
                e[j] = g;
            }
            f = 0.0;
            for (j = 0; j < i; j++) { e[j] /= h; f += e[j] * d[j]; }
            hh = f / (h + h);
            for (j = 0; j < i; j++) { e[j] -= hh * d[j]; }
            for (j = 0; j < i; j++) {
                f = d[j]; g = e[j];
                for (k = j; k <= i-1; k++) { V[k][j] -= (f * e[k] + g * d[k]); }
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
            }
        }
        d[i] = h;
    }
    
    /* Accumulate transformations. */
    
    for (i = 0; i < n-1; i++) {
        V[n-1][i] = V[i][i];
        V[i][i] = 1.0;
        h = d[i+1];
        if (h != 0.0) {
            for (k = 0; k <= i; k++) { d[k] = V[k][i+1] / h;}
            for (j = 0; j <= i; j++) {
                g = 0.0;
                for (k = 0; k <= i; k++) { g += V[k][i+1] * V[k][j]; }
                for (k = 0; k <= i; k++) { V[k][j] -= g * d[k]; }
            }
        }
        for (k = 0; k <= i; k++) { V[k][i+1] = 0.0;}
    }
    for (j = 0; j < n; j++) { d[j] = V[n-1][j]; V[n-1][j] = 0.0; }
    V[n-1][n-1] = 1.0;
    e[0] = 0.0;
}

/* Symmetric tridiagonal QL algorithm.  对称三对角QL算法*/
 
void eigen_decomposition(double A[n][n], double V[n][n], double d[n]) {
	/* 
	A[n][n]：输入n*n矩阵
	V[n][n]：特征向量n*n矩阵
	d[n]：特征值大小为n的数组
	*/
    double e[n];
    double da[3];
    double dt, dat;
    double vet[3];
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            V[i][j] = A[i][j];
        }
    }
    tred2(V, d, e);
    tql2(V, d, e);
    
    /* Sort the eigen values and vectors by abs eigen value */
    da[0]=absd(d[0]); da[1]=absd(d[1]); da[2]=absd(d[2]);
    if((da[0]>=da[1])&&(da[0]>da[2]))
    {
        dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
        d[2]=d[0]; da[2]=da[0];  V[0][2] = V[0][0]; V[1][2] = V[1][0]; V[2][2] = V[2][0];
        d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2]; 
    }
    else if((da[1]>=da[0])&&(da[1]>da[2]))  
    {
        dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
        d[2]=d[1]; da[2]=da[1];  V[0][2] = V[0][1]; V[1][2] = V[1][1]; V[2][2] = V[2][1];
        d[1]=dt;   da[1]=dat;    V[0][1] = vet[0];  V[1][1] = vet[1];  V[2][1] = vet[2]; 
    }
    if(da[0]>da[1])
    {
        dt=d[1];   dat=da[1];    vet[0]=V[0][1];    vet[1]=V[1][1];    vet[2]=V[2][1];
        d[1]=d[0]; da[1]=da[0];  V[0][1] = V[0][0]; V[1][1] = V[1][0]; V[2][1] = V[2][0];
        d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2]; 
    }

    
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
	/* nlhs：输入变量个数，mxArray *plhs[]：输入指针
	nrhs：输出变量个数，const mxArray *prhs[]：输出指针 */

	double *Dxx, *Dxy, *Dxz, *Dyy, *Dyz, *Dzz;
    double *Dvecx, *Dvecy, *Dvecz, *Deiga, *Deigb, *Deigc;

    float *Dxx_f, *Dxy_f, *Dxz_f, *Dyy_f, *Dyz_f, *Dzz_f;
    float *Dvecx_f, *Dvecy_f, *Dvecz_f, *Deiga_f, *Deigb_f, *Deigc_f;

    mwSize output_dims[2]={1, 3};
    // mwSize是Matlab中一种类型，它实际上是int的别名，即int
　　typedef int mwSize;
    double Ma[3][3];
    double Davec[3][3];
    double Daeig[3];
    
    /* Loop variable */
    int i;
    
    /* Size of input */
    const mwSize *idims;
    int nsubs=0;
    
    /* Number of pixels */
    int npixels=1;
    
    /* Check for proper number of arguments. */
    if(nrhs!=6) {
        mexErrMsgTxt("Six inputs are required.");
    } else if(nlhs<3) {
        mexErrMsgTxt("Three or Six outputs are required");
    }
    
   
    /*  Get the number of dimensions */
    // mxGetNumberOfDimensions：获得阵列的维数。
??// int mxGetNumberOfDimensions(const mxArray *array_ptr);
    nsubs = mxGetNumberOfDimensions(prhs[0]);
    // 此时为nsubs =  3
    /* Get the sizes of the inputs */
    // mxGetDimensions：获得一个指向阵列维数大小的数组指针。
??// const int *mxGetDimensions(const mxArray *array_ptr);
    idims = mxGetDimensions(prhs[0]);
    // 每一维的大小 如果是512*512*512 那idims=[512,512,512]
    for (i=0; i<nsubs; i++) { npixels=npixels*idims[i]; }
    
    if(mxGetClassID(prhs[0])==mxDOUBLE_CLASS) {
       /* Assign pointers to each input.输入的是Hessian矩阵的Dxx，Dyy，Dzz... */
        Dxx = (double *)mxGetPr(prhs[0]);
        Dxy = (double *)mxGetPr(prhs[1]);
        Dxz = (double *)mxGetPr(prhs[2]);
        Dyy = (double *)mxGetPr(prhs[3]);
        Dyz = (double *)mxGetPr(prhs[4]);
        Dzz = (double *)mxGetPr(prhs[5]);

        /* Assign pointers to each output. */
        plhs[0] = mxCreateNumericArray(nsubs, idims, mxDOUBLE_CLASS, mxREAL);
        plhs[1] = mxCreateNumericArray(nsubs, idims, mxDOUBLE_CLASS, mxREAL);
        plhs[2] = mxCreateNumericArray(nsubs, idims, mxDOUBLE_CLASS, mxREAL);
        Deiga = mxGetPr(plhs[0]); Deigb = mxGetPr(plhs[1]); Deigc = mxGetPr(plhs[2]);
        if(nlhs==6) {
            /* Main direction (larged eigenvector) */
            plhs[3] = mxCreateNumericArray(nsubs, idims, mxDOUBLE_CLASS, mxREAL);
            plhs[4] = mxCreateNumericArray(nsubs, idims, mxDOUBLE_CLASS, mxREAL);
            plhs[5] = mxCreateNumericArray(nsubs, idims, mxDOUBLE_CLASS, mxREAL);
            Dvecx = mxGetPr(plhs[3]); Dvecy = mxGetPr(plhs[4]); Dvecz = mxGetPr(plhs[5]);
        }
        
		/* 对每一个像素点都构造一个3*3的矩阵Ma */
        for(i=0; i<npixels; i++) {
            Ma[0][0]=Dxx[i]; Ma[0][1]=Dxy[i]; Ma[0][2]=Dxz[i];
            Ma[1][0]=Dxy[i]; Ma[1][1]=Dyy[i]; Ma[1][2]=Dyz[i];
            Ma[2][0]=Dxz[i]; Ma[2][1]=Dyz[i]; Ma[2][2]=Dzz[i];
            eigen_decomposition(Ma, Davec, Daeig);
            Deiga[i]=Daeig[0]; Deigb[i]=Daeig[1]; Deigc[i]=Daeig[2];
            if(nlhs==6) {
                /* Main direction (smallest eigenvector) */
                Dvecx[i]=Davec[0][0];
                Dvecy[i]=Davec[1][0];
                Dvecz[i]=Davec[2][0];
            }
        }
    }
    else if(mxGetClassID(prhs[0])==mxSINGLE_CLASS) {
       /* Assign pointers to each input. */
        Dxx_f = (float *)mxGetPr(prhs[0]);
        Dxy_f = (float *)mxGetPr(prhs[1]);
        Dxz_f = (float *)mxGetPr(prhs[2]);
        Dyy_f = (float *)mxGetPr(prhs[3]);
        Dyz_f = (float *)mxGetPr(prhs[4]);
        Dzz_f = (float *)mxGetPr(prhs[5]);

        /* Assign pointers to each output. */
        plhs[0] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
        plhs[1] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
        plhs[2] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
        Deiga_f  = (float *)mxGetPr(plhs[0]); 
        Deigb_f  = (float *)mxGetPr(plhs[1]); 
        Deigc_f  = (float *)mxGetPr(plhs[2]);
        if(nlhs==6) {
            /* Main direction (smallest eigenvector) */
            plhs[3] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
            plhs[4] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
            plhs[5] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
            Dvecx_f  = (float *)mxGetPr(plhs[3]); 
            Dvecy_f  = (float *)mxGetPr(plhs[4]); 
            Dvecz_f  = (float *)mxGetPr(plhs[5]);
        }
        
        
        for(i=0; i<npixels; i++) {
            Ma[0][0]=(double)Dxx_f[i]; Ma[0][1]=(double)Dxy_f[i]; Ma[0][2]=(double)Dxz_f[i];
            Ma[1][0]=(double)Dxy_f[i]; Ma[1][1]=(double)Dyy_f[i]; Ma[1][2]=(double)Dyz_f[i];
            Ma[2][0]=(double)Dxz_f[i]; Ma[2][1]=(double)Dyz_f[i]; Ma[2][2]=(double)Dzz_f[i];
            eigen_decomposition(Ma, Davec, Daeig);
            Deiga_f[i]=(float)Daeig[0]; 
            Deigb_f[i]=(float)Daeig[1]; 
            Deigc_f[i]=(float)Daeig[2];
            if(nlhs==6) {
                /* Main direction (smallest eigenvector) */
                Dvecx_f[i]=(float)Davec[0][0]; 
                Dvecy_f[i]=(float)Davec[1][0]; 
                Dvecz_f[i]=(float)Davec[2][0];
            }
        }

    }
}
