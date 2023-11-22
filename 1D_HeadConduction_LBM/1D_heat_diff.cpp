//D1Q3 Lattice
//BGK Model

#include <iostream>
#include <math.h>
#include <omp.h>
#include <chrono>
#include<fstream>

using namespace std;
int LX=101;
double alpha = 0.5/3.0;
const int zeta=3;
double tauinv = 1;
int n_iter = 10000000;
int save_every=10000;

double Wlist[zeta] = {2.0/3.0, 1.0/6.0, 1.0/6.0};
int C_x[zeta] = {0,1,-1};

double* temp =new double[LX];
double* temp_old =new double[LX];
double** f;
double** f_eq;
double** ftemp; //intermediate discrete velcotiy distr func 

FILE* tempfile;
char tempdat[256];

void allocatememory()
{
    f = new double*[LX];
    f_eq = new double*[LX];
    ftemp = new double*[LX];
    for(int i=0; i<LX; i++)
    {
        f[i] = new double[zeta];
        f_eq[i] = new double[zeta];
        ftemp[i] = new double[zeta];
    }
}


int main()
{

    system("rm -rf output_folders");
    system("mkdir output_folders");
    allocatememory(); //memory allocation

    for(int i=0; i<LX; i++)
    {
        temp[i]=0.0;
        temp_old[i]=0.0;
    }

    temp[(LX-1)/(2)]=1000.0;
    temp[0]=100.0;
    temp[LX-1]=100.0;

    //inititalzing fi to feq
    
    // for(int i=0; i<LX; i++)
    // {
    //     for(int k=0; k<zeta; k++)
    //     {
    //         f_eq[i][k]=Wlist[k]*temp[i];
    //     }
    // }

    for(int i=1; i<LX-1; i++)
    {
        for(int k=0; k<zeta; k++)
        {
            f_eq[i][k]=Wlist[k]*temp[i];
        }
    }  

    for(int k=0; k<zeta; k++)
    {
       f_eq[0][k]=Wlist[k]*temp[0];
       f_eq[LX-1][k]=Wlist[k]*temp[LX-1]; 
    }

    for(int i=0; i<LX; i++)
    {
        for(int k=0; k<zeta; k++)
        {
            f[i][k]=f_eq[i][k];
        }
    }

    //iteration
    for(int t=0; t<n_iter; t++)
    {

         for(int i=1; i<LX-1; i++)
        {
            for(int k=0; k<zeta; k++)
            {
                f_eq[i][k]=Wlist[k]*temp[i];
            }
        }

        for(int k=0; k<zeta; k++)
        {
            f_eq[0][k]=Wlist[k]*temp[0];//100.0;
            f_eq[LX-1][k]=Wlist[k]*temp[LX-1];//100.0; 
        }

        //collision
        for(int i=0; i<LX; i++)
        {
            for(int k=0; k<zeta; k++)
            {
                ftemp[i][k] = f[i][k]*(1-tauinv) + f_eq[i][k]*tauinv;
            }
        }        

        //[propagation + Periodic boundary connditions] (calling elements to node)
        for(int i=0; i<LX; i++) 
        {
            for(int k=0; k<zeta; k++)
            {
                int destX = i + C_x[k];
                int x = (destX+LX)%LX;
                f[i][k]=ftemp[x][k];
            }
        }

        for(int i=1; i<LX-1; i++)
        {
            double temp_temp=0.0;
            for(int k=0; k<zeta; k++)
            {
                temp_temp+=f[i][k];
            }
            temp[i]=temp_temp;
        }

        float num_sum=0.0;
        float deno_sum=0.0;
        float error = 0.0;

        for(int j=1; j<LX-1; j++)
		{ 
			num_sum += pow((temp[j]-temp_old[j]),2);
			deno_sum += pow((temp_old[j]),2);
		}
		error = sqrt(num_sum/deno_sum);

        std::cout<<"Iteration ="<<" "<<t<<"   |   "<<"Tolerance ="<<" "<<error<<endl;

		if (error<1e-15)
		break;
		
		for(int k=1; k<LX-1; k++)
		{
			temp_old[k]=temp[k];
		}
        
        if(t%save_every==0)
        {
            sprintf(tempdat, "output_folders/temp_output_t%d", t);
            tempfile = fopen(tempdat, "w");
            for(int i=0; i<LX; i++)
            {
                fprintf(tempfile, "%d %f\n", i, temp[i]);
            }
            fclose(tempfile);
        }

        for(int i=0; i<LX; i++)
        {
            temp_old[i]=temp[i];
        }

        
    }

    ofstream f;
    f.open("output.dat", ios::out);
    {
        for(int i=0; i<LX; i++)
        {
            f<<i<<" "<<temp[i]<<endl;
        }
    }
    f.close();

}
