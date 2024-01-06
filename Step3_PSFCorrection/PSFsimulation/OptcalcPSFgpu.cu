/**********************************************
 * This is the source code for PSF calculation *
 * using AMD/NVDA GPU based on OpenCL frame.   *
 * 											   *
 * Maintained by Xuanwen Hua				   *
 * (xwghua@gmail.com)						   *
***********************************************/


#include <math.h>


/* Identify the sub-functions that will be used in the following codes.
/////////////////////////////////////////////////////////
// Here is the kernel function denoted as "__kernel"
// which means it can only run on the OpenCL platform.
**********************************************************/

__global__ void calcPSFgpu(double *pdev_p1,
                        double *pdev_p2,
                        double *pdev_p3,
                        double *pdev_fobj,
                        double *pdev_k,
                        double *pdev_alpha, //*1e9
                        double *pdev_M,
                        double *pdev_wavelen,
                        int *pdev_boundary,
                        int *pdev_centerPT,
                        double *pn_obj,
                        double *pn_sample,
                        //__global float* dev_zeroline,
                        double *dev_xspace, //*1e6
                        double *dev_yspace, //*1e6
                        double *dev_gkx,
                        double *dev_gkw1,
                        double *dev_pattern_quart_re,
                        double *dev_pattern_quart_im) {
/* 
Note that the pattern_quart_xx is from the boundary(or boundary-1 in C
indeces) to the centerPT (or the centerPT-1 in C indeces), and calculation
only works when the row & col indeces meet the 1/8 limitation.
This scheme actually minimize the calculation consumption in case of large
matrices being introduced.
*/
    // double dev_p1 = *pdev_p1;
    // double dev_p2 = *pdev_p2;
    double dev_p3 = *pdev_p3;
    // double dev_fobj = *pdev_fobj;
    double dev_k = *pdev_k;
    double dev_alpha = *pdev_alpha;
    double dev_M = *pdev_M;
    // double dev_wavelen = *pdev_wavelen;
    int dev_boundary = *pdev_boundary;
    int dev_centerPT = *pdev_centerPT;
    double n_obj = *pn_obj;
    double n_sample = *pn_sample;

    __device__ double cos(double x);
    __device__ double sin(double x);
    __device__ double pow(double x, double y);
    __device__ double sqrt(double x);
    __device__ __CUDA_MATH_CRTIMP double j0(double x);
    __device__ __CUDA_MATH_CRTIMP double j1(double x);
    __device__ __CUDA_MATH_CRTIMP double jn(int bessel_n, double x);
    //printf("########## Kernel calculation started ##########");

    int col_pt_total = dev_centerPT - dev_boundary;

    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    long int gid = blockId * (blockDim.x * blockDim.y)+ (threadIdx.y * blockDim.x) + threadIdx.x;

    
    
    int rowpt = gid/col_pt_total; // row a
    int colpt = gid%col_pt_total; // col b
    
    double dev_gkx2i,tau_si,tau_pi,kphid;
    double rim_d = 1e-7; //100 nm
    // double rim_d = -10e-6; //10 um
    // double rim_d = 30e-6; //30 um
    // double rim_d = 50e-6; //50 um
    // double rim_d = 70e-6; //70 um
    // double rim_d = 90e-6; //90 um
    // double rim_d = 110e-6; //110 um
    // double rim_d = 130e-6; //130 um
    // double rim_d = 150e-6; //150 um
    // double rim_d = 170e-6; //170 um
    // printf("col_pt_total = %d\n", col_pt_total);
    // printf("rowpt, colpt, gid = %d, %d, %d\n", rowpt, colpt, gid);
    if ((rowpt<dev_centerPT) && (colpt<dev_centerPT))
    {
        // printf("%d | %d,%d | %lf,%lf,\n",gid,
        //     rowpt, colpt,
        //     dev_xspace[rowpt],dev_yspace[colpt]);

        dev_pattern_quart_re[gid] = 0.0;
        dev_pattern_quart_im[gid] = 0.0;

        if (colpt>=rowpt)
        { // identify the 1/8 boundary
            double x1 = dev_xspace[rowpt + dev_boundary ]; //*1e6
            double x2 = dev_yspace[colpt + dev_boundary ]; //*1e6
            //printf("%lf,%lf\n",x1,x2);
            double xL2normsq = (sqrt(x1*x1 + x2*x2))/(dev_M*1e6);
            //printf("%lf\n", sqrt(x1*x1 + x2*x2));
            double v = dev_k*xL2normsq*sin(dev_alpha);
            double u = 4*dev_k*1e-6 / n_obj * n_sample *dev_p3*pow(sin(dev_alpha/2),2);
            // double phi_p = atan(x2/x1);
            //printf("%lf, %lf\n", u,v);
            //printf("===== ");

            double U0_re = 0; int i;
            for (i=0; i<201; i++){
                // printf("i = %d\n", i);
                dev_gkx2i = asin(sin(dev_gkx[i]) * n_obj / n_sample);
                tau_si = 2 * sin(dev_gkx2i) * cos(dev_gkx[i]) / sin(dev_gkx[i] + dev_gkx2i);
                tau_pi = tau_si / cos (dev_gkx[i] - dev_gkx2i);
                kphid = (0.0-dev_k) /n_obj * rim_d * (n_obj * cos(dev_gkx[i]) - n_sample * cos(dev_gkx2i));
                // printf("dev_gkx2i,tau_si,tau_pi,kphid = (%f,%f, %f,%f)\n",dev_gkx2i,tau_si,tau_pi,kphid);
                //////////////////////////////////////////////
                U0_re = U0_re + (sqrt(cos(dev_gkx[i])))*
                (tau_si+tau_pi * cos(dev_gkx2i))*
                (cos(kphid + (u/2)*(pow(sin(dev_gkx2i/2),2))/pow(sin(dev_alpha/2),2)))*
                (j0(sin(dev_gkx[i])/sin(dev_alpha)*v))*(sin(dev_gkx[i])) *
                dev_gkw1[i];
                // if (isnan(U0_re)==1){
                //     printf("i, dev_gkx[i],dev_gkx2i,tau_si,tau_pi,kphid = (%d,%f,%f,%f, %f,%f)\n",i,dev_gkx[i],dev_gkx2i,tau_si,tau_pi,kphid);
                // }
                // else
                //     ;
            }
            

            double U0_im = 0; int j;
            for (j=0; j<201; j++){
                // printf("j = %d\n", j);
                dev_gkx2i = asin(sin(dev_gkx[j]) * n_obj / n_sample);
                tau_si = 2 * sin(dev_gkx2i) * cos(dev_gkx[j]) / sin(dev_gkx[j] + dev_gkx2i);
                tau_pi = tau_si / cos (dev_gkx[j] - dev_gkx2i);
                kphid = (0.0-dev_k) /n_obj * rim_d * (n_obj * cos(dev_gkx[j]) - n_sample * cos(dev_gkx2i));
                //////////////////////////////////////////////
                U0_im = U0_im + (sqrt(cos(dev_gkx[j])))*
                (tau_si+tau_pi * cos(dev_gkx2i))*
                (sin(kphid + (u/2)*(pow(sin(dev_gkx2i/2),2))/pow(sin(dev_alpha/2),2)))*
                (j0(sin(dev_gkx[j])/sin(dev_alpha)*v))*(sin(dev_gkx[j])) *
                dev_gkw1[j];
            }

            // double U0_re2 = 0; int ii;
            // for (ii=0; ii<201; ii++){
            //     // printf("i = %d\n", i);
            //     dev_gkx2i = asin(sin(dev_gkx[ii]) * n_obj / n_sample);
            //     tau_si = 2 * sin(dev_gkx2i) * cos(dev_gkx[ii]) / sin(dev_gkx[ii] + dev_gkx2i);
            //     tau_pi = tau_si / cos (dev_gkx[ii] - dev_gkx2i);
            //     kphid = -dev_k /n_obj * rim_d * (n_obj * cos(dev_gkx[ii]) - n_sample * cos(dev_gkx2i));
            //     //////////////////////////////////////////////
            //     U0_re2 = U0_re2 + (sqrt(cos(dev_gkx[ii])))*
            //     (tau_pi * sin(dev_gkx2i))*
            //     (cos(kphid + (u/2)*(pow(sin(dev_gkx2i/2),2))/pow(sin(dev_alpha/2),2)))*
            //     (j0(sin(dev_gkx[ii])/sin(dev_alpha)*v))*(sin(dev_gkx[ii])) *
            //     dev_gkw1[ii];
            // //     //printf("U0_re = %lf\n", U0_re);
            // }

            // double U0_im2 = 0; int jj;
            // for (jj=0; jj<201; jj++){
            //     // printf("j = %d\n", j);
            //     dev_gkx2i = asin(sin(dev_gkx[jj]) * n_obj / n_sample);
            //     tau_si = 2 * sin(dev_gkx2i) * cos(dev_gkx[jj]) / sin(dev_gkx[jj] + dev_gkx2i);
            //     tau_pi = tau_si / cos (dev_gkx[jj] - dev_gkx2i);
            //     kphid = -dev_k /n_obj * rim_d * (n_obj * cos(dev_gkx[jj]) - n_sample * cos(dev_gkx2i));
            //     //////////////////////////////////////////////
            //     U0_im2 = U0_im2 + (sqrt(cos(dev_gkx[jj])))*
            //     (tau_pi * sin(dev_gkx2i))*
            //     (sin(kphid + (u/2)*(pow(sin(dev_gkx2i/2),2))/pow(sin(dev_alpha/2),2)))*
            //     (j0(sin(dev_gkx[jj])/sin(dev_alpha)*v))*(sin(dev_gkx[jj])) *
            //     dev_gkw1[jj];
            // }

            double U0_re3 = 0; int iii;
            for (iii=0; iii<201; iii++){
                // printf("i = %d\n", i);
                dev_gkx2i = asin(sin(dev_gkx[iii]) * n_obj / n_sample);
                tau_si = 2 * sin(dev_gkx2i) * cos(dev_gkx[iii]) / sin(dev_gkx[iii] + dev_gkx2i);
                tau_pi = tau_si / cos (dev_gkx[iii] - dev_gkx2i);
                kphid = (0.0-dev_k) /n_obj * rim_d * (n_obj * cos(dev_gkx[iii]) - n_sample * cos(dev_gkx2i));
                //////////////////////////////////////////////
                U0_re3 = U0_re3 + (sqrt(cos(dev_gkx[iii])))*
                (tau_si-tau_pi * cos(dev_gkx2i))*
                (cos(kphid + (u/2)*(pow(sin(dev_gkx2i/2),2))/pow(sin(dev_alpha/2),2)))*
                (jn(2,sin(dev_gkx[iii])/sin(dev_alpha)*v))*(sin(dev_gkx[iii])) *
                // (j1(sin(dev_gkx[iii])/sin(dev_alpha)*v))*(sin(dev_gkx[iii])) *
                dev_gkw1[iii];
            //     //printf("U0_re = %lf\n", U0_re);
            }

            double U0_im3 = 0; int jjj;
            for (jjj=0; jjj<201; jjj++){
                // printf("j = %d\n", j);
                dev_gkx2i = asin(sin(dev_gkx[jjj]) * n_obj / n_sample);
                tau_si = 2 * sin(dev_gkx2i) * cos(dev_gkx[jjj]) / sin(dev_gkx[jjj] + dev_gkx2i);
                tau_pi = tau_si / cos (dev_gkx[jjj] - dev_gkx2i);
                kphid = (0.0-dev_k) /n_obj * rim_d * (n_obj * cos(dev_gkx[jjj]) - n_sample * cos(dev_gkx2i));
                //////////////////////////////////////////////
                U0_im3 = U0_im3 + (sqrt(cos(dev_gkx[jjj])))*
                (tau_si-tau_pi * cos(dev_gkx2i))*
                (sin(kphid + (u/2)*(pow(sin(dev_gkx2i/2),2))/pow(sin(dev_alpha/2),2)))*
                (jn(2,sin(dev_gkx[jjj])/sin(dev_alpha)*v))*(sin(dev_gkx[jjj])) *
                dev_gkw1[jjj];
            }
            /*************************************/
            // double U0_re = gk_re(dev_gkx, dev_gkw1, dev_alpha, u, v);
            // double U0_im = gk_im(dev_gkx, dev_gkw1, dev_alpha, u, v);

            /***********************************/
            dev_pattern_quart_re[gid] = U0_re-U0_re3; //*1e18
            dev_pattern_quart_im[gid] = U0_im-U0_im3; //*1e18
            //printf("Koi_re, Koi_im = (%f, %f)\n",Koi_re,Koi_im);
            //printf("dev_pattern_quart_re, dev_pattern_quart_im = (%f, %f)\n",U0_re*Koi_re - U0_im*Koi_im, U0_re*Koi_im + U0_im*Koi_re);
            //printf("u,v: %lf, %lf\n",u,v);
            // printf("U0_re,U0_re3, U0_im,U0_im3 = (%f,%f, %f,%f)\n",U0_re,U0_re3,U0_im,U0_im3);
            
            
        }
        else
            ;
    }
    
    
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
***********************/