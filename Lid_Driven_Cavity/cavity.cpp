#include<iostream>
#include<math.h>
#include<cmath>
#include<chrono>
#include<omp.h>
#include<cfloat>
#include<iomanip>
#include<fstream>
using namespace std;

int LX=101;
int LY=101;
int zeta=9;
int last_itr=100000;
int itr_save=10000;

double kinVisc = 0.1;
double U_wall=0.1;
double tauinv = (1.0)/(3.0*kinVisc + 0.5);

double ds=(1.0/(LX-1));
const int n=201;
double a[101], b[101];
double x,y;
double RMSerror;

FILE *outputFile;
char vel[200];

int Cx[9] = {0, 1, 0,-1, 0, 1, -1, -1, 1};
int Cy[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
double Wlist[9] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
double invList[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

double** density;
double** u_x;
double** u_y;
double** density_old;
double** u_x_old;
double** u_y_old;

double*** f;
double*** f_eq;
double*** ftemp;

void allocatememory()
{
    density = new double*[LX];
    u_x = new double*[LX];
    u_y = new double*[LX];
    density_old = new double*[LX];
    u_x_old = new double*[LX];
    u_y_old = new double*[LX];


    f = new double**[LX];
    f_eq = new double**[LX];
    ftemp = new double**[LX];

    for(int i=0; i<LX; i++)
    {
        density[i] = new double[LY];
        u_x[i]= new double[LY];
        u_y[i]= new double[LY];
        density_old[i] = new double[LX];
        u_x_old[i] = new double[LX];
        u_y_old[i] = new double[LX];

        f[i]= new double*[LY];
        f_eq[i] = new double*[LY];
        ftemp[i] = new double*[LY];
    }

    for(int i=0; i<LX; i++)
    {
        for(int j=0; j<LY; j++)
        {
            f[i][j]= new double[zeta];
            f_eq[i][j] = new double[zeta];
            ftemp[i][j] = new double[zeta];
        }
    }
}

void boundaryConditions()
{
    for(int i=1; i<LX-1; i++) //upper wall
    {
        int j= (LY-2);
        for(int k=0; k<zeta; k++)
        {   
            if(k==5 || k==6 || k==2)
            {
            int ktemp=invList[k];
            f[i][j][ktemp] = ftemp[i][j][k] - 6*Wlist[k]*density[i][j]*(Cx[k]*U_wall);
            }
        }    
    }

    for(int i=1; i<LX-1; i++) //bottom wall
    {
        int j=1;
        for(int k=0; k<zeta; k++)
        {
        if(k==7 || k==4 || k==8)
        {
        int ktemp=invList[k];
        f[i][j][ktemp] = ftemp[i][j][k];
        }
        }
    }

    for(int j=1; j<LY-1; j++) //right wall
    {
        int i=(LX-2);
        for(int k=0; k<zeta; k++)
        {
        if(k==5 || k==1 || k==8)
        {
        int ktemp=invList[k];
        f[i][j][ktemp] = ftemp[i][j][k];
        }
        }
    }

    for(int j=1; j<LY-1; j++) //left wall
    {
        int i=1;
        for(int k=0; k<zeta; k++)
        {
        if(k==7 || k==6 || k==3)
        {
        int ktemp=invList[k];
        f[i][j][ktemp] = ftemp[i][j][k];
        }
        }
    }

}

int main()
{

    system("rm -rf output");
    system("mkdir output");

    allocatememory();

    for(int i=0; i<LX; i++) //initilization
    {
        for(int j=0; j<LY; j++)
        {
            u_x[i][j]=0.0;
            u_y[i][j]=0.0;
            density[i][j]=1.0;
            if(j == LY - 1){
                u_x[i][j]=U_wall;
            }
        }
    }

    for(int i=1; i<LX-1; i++)
    {
        for(int j=1; j<LY-1; j++)
        {
            double udotu = u_x[i][j]*u_x[i][j] + u_y[i][j]*u_y[i][j];
            for(int k=0; k<zeta; k++)
            {
                double udotc = u_x[i][j]*Cx[k] + u_y[i][j]*Cy[k];
                f_eq[i][j][k] = Wlist[k]*density[i][j]*(1 + 3*(udotc) + 4.5*(udotc*udotc) - 1.5*(udotu));
            }
        }
    }
    
    for(int i=1; i<LX-1; i++)
    {
        for(int j=1; j<LY-1; j++)
        {
            for(int k=0; k<zeta; k++)
            {
                f[i][j][k]=f_eq[i][j][k];
            }
        }
    }

    for(int itr=0; itr<last_itr; itr++) //iteration starts
    {

        for(int i=1; i<LX-1; i++)
        {
            for(int j=1; j<LY-1; j++)
            {
                double udotu = u_x[i][j]*u_x[i][j] + u_y[i][j]*u_y[i][j];
                for(int k=0; k<zeta; k++)
                {
                    double udotc = u_x[i][j]*Cx[k] + u_y[i][j]*Cy[k];
                    f_eq[i][j][k] = Wlist[k]*density[i][j]*(1 + 3*(udotc) + 4.5*(udotc*udotc) - 1.5*(udotu));
                }
            }
        }

        // for(int i=1; i<LX-1; i++)
        // {
        //     for(int j=1; j<LY-1; j++)
        //     {
        //         for(int k=0; k<zeta; k++)
        //         {
        //             f[i][j][k]=f_eq[i][j][k];
        //         }
        //     }
        // }

        for(int i=1; i<LX-1; i++) //collision
        {
            for(int j=1; j<LY-1; j++)
            {
                for(int k=0; k<zeta; k++)
                {
                    ftemp[i][j][k]=f[i][j][k]*(1-tauinv) + f_eq[i][j][k]*tauinv;
                }
            }
        }

        for(int i=1; i<LX-1; i++) //streaming
        {
            for(int j=1; j<LY-1; j++)
            {
                for(int k=0; k<zeta; k++)
                {
                    int destX = i-Cx[k];
                    int destY = j-Cy[k];
                    f[i][j][k]=ftemp[destX][destY][k];
                }
            }
        }

        boundaryConditions();

        for(int i=1; i<LX-1; i++)
        {
            for(int j=1; j<LY-1; j++)
            {
               density_old[i][j]=density[i][j];
               u_x_old[i][j]=u_x[i][j];
               u_y_old[i][j]=u_y[i][j];
            }
        }

        for(int i=1; i<LX-1; i++) //updating density
        {
            for(int j=1; j<LY-1; j++)
            {
                double densitytemp = 0.0;
                double uXtemp = 0.0;
                double uYtemp = 0.0;
                for(int k=0; k<zeta; k++)
                {
                    densitytemp += f[i][j][k];
                    uXtemp += Cx[k]*f[i][j][k];
                    uYtemp += Cy[k]*f[i][j][k];
                }
                density[i][j]=densitytemp;
                u_x[i][j] = uXtemp/(density[i][j]+FLT_MIN);
                u_y[i][j] = uYtemp/(density[i][j]+FLT_MIN);
            }
        }

        double Denstolerance_temp = 0.0;
        double uxtolerance_temp = 0.0;
        double uytolerance_temp = 0.0;
        double Denstolerance_tempdeno = 0.0;
        double uxtolerance_tempdeno = 0.0;
        double uytolerance_tempdeno = 0.0;
        for(int i=1; i<LX-1; i++)
        {
            for(int j=1; j<LY-1; j++)
            {
                Denstolerance_temp += (density_old[i][j]-density[i][j])*(density_old[i][j]-density[i][j]);
                Denstolerance_tempdeno += (density_old[i][j]*density_old[i][j]);

                uxtolerance_temp += (u_x_old[i][j]-u_x[i][j])*(u_x_old[i][j]-u_x[i][j]);
                uxtolerance_tempdeno += (u_x_old[i][j]*u_x_old[i][j]);

                uytolerance_temp += (u_y_old[i][j]-u_y[i][j])*(u_y_old[i][j]-u_y[i][j]);
                uytolerance_tempdeno += (u_y_old[i][j]*u_y_old[i][j]);
            }
        }
        double Denstolerance = sqrt(Denstolerance_temp/Denstolerance_tempdeno);
        double uxtolerance = sqrt(uxtolerance_temp/uxtolerance_tempdeno);
        double uytolerance = sqrt(uytolerance_temp/uytolerance_tempdeno);

        if(uxtolerance<=1e-7 && uytolerance<=1e-5)
        {
            break;
        }

        if(itr%itr_save == 0)
        {
            sprintf(vel, "output/output_t%d.dat", itr);
            outputFile = fopen(vel, "w");
            for(int i=0; i<LX; i++)
            {
                for(int j=0; j<LY; j++)
                {
                    fprintf(outputFile, "%d %d %f %f %f %f\n", i,j,u_x[i][j],u_y[i][j],sqrt(pow(u_x[i][j],2) + pow(u_y[i][j],2)),density[i][j]);
                }
            }
            fclose(outputFile);
        }

        int space = 20;
        cout<<"iteration= "<<itr;//<<right<< setfill(' ') << setw(space)<<" | ";
        cout.precision(12);
        cout.width(16);
        cout<<"resRho= "<<round(Denstolerance * 1e10)/1e10;//<<right << setfill(' ') << setw(space)<<" | ";
        cout.precision(12);
        cout.width(16);
        cout<<"resU_x= "<<round(uxtolerance * 1e10)/1e10;//<<right << setfill(' ') << setw(space)<<" | ";
        cout.precision(12);
        cout.width(16);
        cout<<"resU_y= "<<round(uytolerance * 1e10)/1e10<<endl;//<<right << setfill(' ') << setw(space)<<endl;
    }

    ofstream f;
    f.open("U_Re100.dat", ios::out);
    {
        for(int j=1; j<LY-1; j++)
        {
            int i = (LX-1)/(2.0);
            f<<j<<" "<<u_x[i][j]<<endl;
        }
    }
    f.close();

    f.open("V_Re100.dat", ios::out);
    {
        for(int i=1; i<LY-1; i++)
        {
            int j = (LY-1)/(2.0);
            f<<i<<" "<<u_y[i][j]<<endl;
        }
    }
    f.close();
    
    cout<<"\n"<<endl;
    cout<<"Iterations complete"<<endl;

}
