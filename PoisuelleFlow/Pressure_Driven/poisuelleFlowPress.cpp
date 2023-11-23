//Re=100;
//D2Q9 Lattice- BGK model

#include<iostream>
#include<math.h>
#include<cmath>
#include<chrono>
#include<omp.h>
#include<cfloat>
#include<iomanip>
#include<fstream>

using namespace std;

int LX=501;
int LY=51;
int zeta=9;
int last_itr=100000;
int save = 10000;
double kinVisc=0.05;
double ds=0.5;
double delRho = 0.03;
double tauinv = (1.0)/(3.0*kinVisc + 0.5);

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
    for(int j=1; j<LY-1; j++) //left wall
    {
        int i=1;
        double u_wX = u_x[i][j] + (1.0/2.0)*(u_x[i][j] - u_x[i+1][j]);
        double u_wY = u_y[i][j] + (1.0/2.0)*(u_y[i][j] - u_y[i+1][j]);
        // double u_wY = u_y[i][j] + (1.0/2.0)*(u_y[i][j] - u_y[i][j+1]);
        double uDotu = u_wX*u_wX + u_wY*u_wY;
        for(int k=0; k<zeta; k++)
        {
            double cDotu_w = (Cx[k])*(u_wX) + (Cy[k])*(u_wY);
            int tempK = invList[k];
            if(k==7 || k==3 || k==6)
            {
            f[i][j][tempK] = -ftemp[i][j][k] + 2.0*Wlist[k]*density[i-1][j]*(1 + 4.5*(cDotu_w*cDotu_w)- 1.5*(uDotu));
            }
        }
    }

    for(int j=1; j<LY-1; j++) //right wall
    {
        int i=(LX-2);
        double u_wX = u_x[i][j] + (1.0/2.0)*(u_x[i][j] - u_x[i-1][j]);
        double u_wY = u_y[i][j] + (1.0/2.0)*(u_y[i][j] - u_y[i-1][j]);
        // double u_wY = u_y[i][j] + (1.0/2.0)*(u_y[i][j] - u_y[i][j+1]);
        double uDotu = u_wX*u_wX + u_wY*u_wY;
        for(int k=0; k<zeta; k++)
        {
            double cDotu_w = (Cx[k])*(u_wX) + (Cy[k])*(u_wY);
            int tempK = invList[k];
            if(k==5 || k==1 || k==8)
            {
            f[i][j][tempK] = -ftemp[i][j][k] + 2.0*Wlist[k]*density[i+1][j]*(1 + 4.5*(cDotu_w*cDotu_w)- 1.5*(uDotu));
            }
        }
    }

    for(int i=0; i<LX; i++) //bottomWall
    {
        int j=1;
        for(int k=0; k<zeta; k++)
        {
            if(k==7 || k==8|| k==4)
            {
                int tempK = invList[k];
                f[i][j][tempK]=ftemp[i][j][k];
            }
        }
    }


    for(int i=0; i<LX; i++) //Top Wall
    {
        int j=LY-2;
        for(int k=0; k<zeta; k++)
        {
            if(k==5 || k==2|| k==6)
            {
                int tempK = invList[k];
                f[i][j][tempK]=ftemp[i][j][k];
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
        for(int j=1; j<LY-1; j++)
        {
            {

            u_x[i][j]=0.0;
            u_y[i][j]=0.0;
            density[i][j]=1.0;

            if(i== 0)
            {
                density[i][j] = 1.0 + delRho*0.5;
            }

            if(i == (LX-1))
            {
                density[i][j] = 1.0 - delRho*0.5;
            }

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
                f_eq[i][j][k] = Wlist[k]*density[i][j]*(1 + 3.0*(udotc) + 4.5*(udotc*udotc) - 1.5*(udotu));
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
                    f_eq[i][j][k] = Wlist[k]*density[i][j]*(1 + 3.0*(udotc) + 4.5*(udotc*udotc) - 1.5*(udotu));
                }
            }
        }

    for(int i=1; i<LX-1; i++) //collison
    {
        for(int j=1; j<LY-1; j++)
        {
            for(int k=0; k<zeta; k++)
            {
                ftemp[i][j][k] = f[i][j][k]*(1-tauinv) + f_eq[i][j][k]*(tauinv);
            }
        }
    }

    for(int i=1; i<LX-1; i++) //streaming
    {
        for(int j=1; j<LY-1; j++)
        {
            for(int k=0; k<zeta; k++)
            {
                int destX = i - Cx[k];
                int destY = j - Cy[k];
                int x = ((destX+LX)%LX);
                int y = ((destY+LY)%LY);
                f[i][j][k] = ftemp[x][y][k];
            }
        }
    }

    boundaryConditions();

    // for(int i=1; i<LX-1; i++)
    // {
    //     for(int j=1; j<LY-1; j++)
    //     {
    //         double tempDensity = 0.0;
    //         for(int k=0; k<zeta; k++)
    //         {
    //             tempDensity += f[i][j][k];
    //         }
    //         density[i][j] = tempDensity;

    //     }
    // }

    for(int i=0; i<LX; i++) //residue calculation 
    {
        for(int j=0; j<LY; j++)
        {
            u_x_old[i][j]=u_x[i][j];
            u_y_old[i][j]=u_y[i][j];
            density_old[i][j]=density[i][j];
        }
    }

    for(int i=1; i<LX-1; i++)
    {
        for(int j=1; j<LY-1; j++)
        {
            double tempDensity = 0.0;
            double tempUx = 0.0;
            double tempUy = 0.0;
            for(int k=0; k<zeta; k++)
            {
                tempDensity += f[i][j][k];
                tempUx += (Cx[k]*f[i][j][k]);
                tempUy += (Cy[k]*f[i][j][k]);
            }
            density[i][j] = tempDensity;
            u_x[i][j] = tempUx/density[i][j];
            u_y[i][j] = tempUy/density[i][j];

        }
    }

    double tempToluxnum = 0.0;
    double tempToluynum = 0.0;
    double tempTolrhonum = 0.0;

    double tempToluxden = 0.0;
    double tempToluyden = 0.0;
    double tempTolrhoden = 0.0;

    for(int i=0; i<LX; i++)
    {
        for(int j=0; j<LY; j++)
        {
            tempToluxnum += (u_x_old[i][j]-u_x[i][j])*(u_x_old[i][j]-u_x[i][j]);
            tempToluynum += (u_y_old[i][j]-u_y[i][j])*(u_y_old[i][j]-u_y[i][j]);
            tempTolrhonum += (density_old[i][j]-density[i][j])*(density_old[i][j]-density[i][j]);

            tempToluxden += (u_x_old[i][j])*(u_x_old[i][j]);
            tempToluyden += (u_y_old[i][j])*(u_y_old[i][j]);
            tempTolrhoden += (density_old[i][j])*(density_old[i][j]);
        }
    }

    double UxTolerance = sqrt(tempToluxnum/tempToluxden);
    double UyTolerance = sqrt(tempToluynum/tempToluyden);
    double RhoTolerance = sqrt(tempTolrhonum/tempTolrhoden);

    cout<<"Density Tolerance ="<<" "<<RhoTolerance<<"  |  "<<"Ux Tolerance ="<<" "<<UxTolerance<<"  |  "<<"Uy Tolerance ="<<UyTolerance<<endl;


    if(itr%save == 0)
    {
        sprintf(vel, "output/output_t%d.dat", itr);
        outputFile = fopen(vel, "w");
            for(int i=1; i<LX-1; i++)
            {
                for(int j=1; j<LY-1; j++)
                {
                    fprintf(outputFile, "%d %d %f %f %f %f\n", i,j,u_x[i][j],u_y[i][j],sqrt(pow(u_x[i][j],2) + pow(u_y[i][j],2)),density[i][j]);
                }
            }
        fclose(outputFile);
    }

    cout<<"iteration ="<<itr<<endl;

     if(UxTolerance<1e-7 && UyTolerance<1e-7 && RhoTolerance<1e-5)
    {
        break;
    }

    }


    double tempj[LY];
    tempj[0]=0;
    for(int j=1; j<LY-1; j++)
    {
        tempj[j]=tempj[j-1]+ds;
    }
    for(int j=2; j<LY-2; j++)
    {
        tempj[j]=tempj[j-1]+1;
    }

    ofstream f;
    f.open("test", ios::out);
    {
        for(int j=0; j<LY; j++)
        {
            int i=(LX-2);
            f<<j<<" "<<u_x[i][j]<<endl;
        }
    }
    f.close();



}



