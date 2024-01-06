/***********************************************/
/* This is the source code for PSF calculation *
 * using AMD/NVDA GPU based on OpenCL frame.   *
 * 											   *
 * Maintained by Xuanwen Hua				   *
 * (xwghua@gmail.com)						   */
/***********************************************/

#define PYOPENCL_DEFINE_CDOUBLE
//#include <pyopencl-complex.h>
//#include <math.h>
#include <pyopencl-bessel-j.cl>
//#include <E:\Google Drive\thesisproj\PostOBJgpu\integral\qags_re.cl>
//#include <E:\Google Drive\thesisproj\PostOBJgpu\integral\qags_im.cl>
//#include <qags_re.cl>
//#include <qags_im.cl>
// Identify the sub-functions that will be used in the following codes.
static inline
double Integrate_re(double theta, double alpha, double u, double v);
static inline
double Integrate_im(double theta, double alpha, double u, double v);
static inline
double gk_re(__global double *int_gkx, __global double *int_gkw1, double int_alpha, double int_u, double int_v);
static inline
double gk_im(__global double *int_gkx, __global double *int_gkw1, double int_alpha, double int_u, double int_v);

static inline
double Integrate_re(double theta, double alpha, double u, double v){
    return (sqrt(cos(theta)))*(1+cos(theta))*
        (cos((u/2)*(pow(sin(theta/2),2))/pow(sin(alpha/2),2)))*
        (bessel_j0((bessel_j_scalar_type)(sin(theta)/sin(alpha)*v)))*
        (sin(theta));
}

static inline
double Integrate_im(double theta, double alpha, double u, double v){
    return (sqrt(cos(theta)))*(1+cos(theta))*
        (sin((u/2)*(pow(sin(theta/2),2))/pow(sin(alpha/2),2)))*
        (bessel_j0((bessel_j_scalar_type)(sin(theta)/sin(alpha)*v)))*
        (sin(theta));
}

static inline
double gk_re(__global double *int_gkx, __global double *int_gkw1, double int_alpha, double int_u, double int_v){
    double U0_re = 0; int i;
    for (i=0; i<201; i++){
        U0_re = U0_re + Integrate_re(int_gkx[i],int_alpha,int_u,int_v) * int_gkw1[i];
        //printf("U0_re = %lf\n", U0_re);
    }
    return U0_re;
}

static inline
double gk_im(__global double *int_gkx, __global double *int_gkw1, double int_alpha, double int_u, double int_v){
    double U0_im = 0; int i;
    for (i=0; i<201; i++){
        U0_im = U0_im + Integrate_im(int_gkx[i],int_alpha,int_u,int_v) * int_gkw1[i];
    }
    return U0_im;
}

/////////////////////////////////////////////////////////
// Here is the kernel function denoted as "__kernel"
// which means it can only run on the OpenCL platform.
/////////////////////////////////////////////////////////
__kernel void calcPSFgpu(double dev_p1,
                        double dev_p2,
                        double dev_p3,
                        double dev_fobj,
                        double dev_k,
                        double dev_alpha, //*1e9
                        double dev_M,
                        double dev_wavelen,
                        int dev_boundary,
                        int dev_centerPT,
                        //__global float* dev_zeroline,
                        __global double* dev_xspace, //*1e6
                        __global double* dev_yspace, //*1e6
                        __global double* dev_gkx,
                        __global double* dev_gkw1,
                        __global double* dev_pattern_quart_re,
                        __global double* dev_pattern_quart_im) {
/* 
Note that the pattern_quart_xx is from the boundary(or boundary-1 in C
indeces) to the centerPT (or the centerPT-1 in C indeces), and calculation
only works when the row & col indeces meet the 1/8 limitation.
This scheme actually minimize the calculation consumption in case of large
matrices being introduced.
*/
    //printf("########## Kernel calculation started ##########");
	int gid = get_global_id(0);
    dev_pattern_quart_re[gid] = 0.0;
    dev_pattern_quart_im[gid] = 0.0;
    int col_pt_total = dev_centerPT - dev_boundary;
    int rowpt = gid/col_pt_total; // row a
    int colpt = gid%col_pt_total; // col b

    if (colpt>=rowpt)
    { // identify the 1/8 boundary
        double x1 = dev_xspace[rowpt + dev_boundary ]; //*1e6
        double x2 = dev_yspace[colpt + dev_boundary ]; //*1e6
        //printf("%lf,%lf\n",x1,x2);
        double xL2normsq = (sqrt(x1*x1 + x2*x2))/(dev_M*1e6);
        //printf("%lf\n", sqrt(x1*x1 + x2*x2));
        double v = dev_k*xL2normsq*sin(dev_alpha);
        double u = 4*dev_k*1e-6*dev_p3*pow(sin(dev_alpha/2),2);
        //printf("%lf, %lf\n", u,v);
        //double Koi_re = dev_M/(pow(dev_fobj*dev_wavelen,2))*cos(-u/(4*pow(sin(dev_alpha/2),2)));
        //double Koi_im = dev_M/(pow(dev_fobj*dev_wavelen,2))*sin(-u/(4*pow(sin(dev_alpha/2),2)));
        //printf("%lf, %lf\n", Koi_re,Koi_im);
        //printf("===== ");
        double U0_re = gk_re(dev_gkx, dev_gkw1, dev_alpha, u, v);
        double U0_im = gk_im(dev_gkx, dev_gkw1, dev_alpha, u, v);

        // printf("%lf, %lf \n", v,U0_re);
        //printf("start running qags............|..|..............\n");
        
        //dev_pattern_quart_re[gid] = U0_re*Koi_re - U0_im*Koi_im; //*1e18
        //dev_pattern_quart_im[gid] = U0_re*Koi_im + U0_im*Koi_re; //*1e18
        dev_pattern_quart_re[gid] = U0_re; //*1e18
        dev_pattern_quart_im[gid] = U0_im; //*1e18
        //printf("Koi_re, Koi_im = (%f, %f)\n",Koi_re,Koi_im);
        //printf("dev_pattern_quart_re, dev_pattern_quart_im = (%f, %f)\n",U0_re*Koi_re - U0_im*Koi_im, U0_re*Koi_im + U0_im*Koi_re);
        //printf("u,v: %lf, %lf\n",u,v);
        //printf("U0_re, U0_im = (%f, %f)\n",U0_re,U0_im);
        
        
    }
    else
        ;
    
    //printf("########## kernel completed! ##########\n");
}



/************************
// ------------------------------
// This part is for reference as it's a CPU scheme
// ------------------------------
    for (unsigned int a = dev_boundary,a<=dev_centerPT,a++){
        int x1 = dev_xspace[a-1];
        float* patternLine = dev_zeroline;
        for (unsigned int b = a,b<=dev_centerPT,b++){
            int x2 = dev_yspace[b-1];
            float xL2normsq = (sqrt(pow(x1+dev_M*dev_p1,2) + pow(x2+dev_M*dev_p2,2)))/dev_M;
            float v = dev_k*xL2normsq*sin(dev_alpha);
            float u = 4*dev_k*dev_p3*pow(sin(dev_alpha/2),2);
            float Koi_re = dev_M/(pow(dev_fobj*dev_wavelen,2))*cos(u/(4*pow(sin(dev_alpha/2),2)));
            float Koi_im = dev_M/(pow(dev_fobj*dev_wavelen,2))*cos(u/(4*pow(sin(dev_alpha/2),2)));
            float U_re = NumInt_re(dev_alpha,u,v,0,dev_alpha,1e-4);
            float U_im = NumInt_im(dev_alpha,u,v,0,dev_alpha,1e-4);
            float U0_re = U_re*Koi_re - U_im*Koi_im;
            float U0_im = U_re*Koi_im + U_im*Koi_re;

        }
    }
************************/