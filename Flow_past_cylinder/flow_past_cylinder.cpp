//Re=100;
//D2Q9 - BGK Model

#include<iostream>
#include<math.h>
#include<cmath>
#include<chrono>
#include<omp.h>
#include<cfloat>
#include<iomanip>
#include<fstream>

using namespace std;

int LX=1001;
int LY=101;
int zeta=9;
double tauinv = (1.0/0.65);
int last_itr=100000;
int save = 1000;
double kinVisc=0.05;
double radius = 25.0;
double centX = (LX-1)/2;
double centY = (LY-1)/2;
double u_inlet = 0.1;
double forceX;
double forceY;

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
int** obstacle;

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
    obstacle = new int*[LX];


    f = new double**[LX];
    f_eq = new double**[LX];
    ftemp = new double**[LX];

    for(int i=0; i<LX; i++)
    {
        density[i] = new double[LY];
        u_x[i]= new double[LY];
        u_y[i]= new double[LY];
        density_old[i] = new double[LY];
        u_x_old[i] = new double[LY];
        u_y_old[i] = new double[LY];
        obstacle[i] = new int[LY];

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

void cylinder()
{
    for(int i=0; i<LX; i++)
    {
        for(int j=0; j<LY; j++)
        {
            double distfrmcntr = sqrt((i-((LX-1)/5))*(i-((LX-1)/5)) + (j-((LY-1)/2))*(j-((LY-1)/2)));
            if(distfrmcntr<=radius)
            {
                obstacle[i][j]=1;
            }
            else
            {
                obstacle[i][j]=0;
            }
        }
    }
}

void boundaryConditions()
{
   for(int j=0; j<LY; j++) //left wall - velocity inlet
    {
        int i=1;
        // double u_wY = u_y[i][j] + (1.0/2.0)*(u_y[i][j] - u_y[i][j+1]);
        // double uDotu = u_wX*u_wX + u_wY*u_wY;
        for(int k=0; k<zeta; k++)
        {
            double CdotU = u_inlet*Cx[k];
            int tempK = invList[k];
            if(k==7 || k==3 || k==6)
            {
            f[i][j][tempK] = ftemp[i][j][k] - 6.0*Wlist[k]*density[i-1][j]*(CdotU);
            }
        }
    }

    for(int j=0; j<LY; j++) //right wall- zeroGradient
    {
        int i=(LX-1);
        // double u_wX = u_x[i][j] + (1.0/2.0)*(u_x[i][j] - u_x[i-1][j]);
        // double u_wY = u_y[i][j] + (1.0/2.0)*(u_y[i][j] - u_y[i-1][j]);
        // double u_wY = u_y[i][j] + (1.0/2.0)*(u_y[i][j] - u_y[i][j+1]);
        // double uDotu = u_wX*u_wX + u_wY*u_wY;
        for(int k=0; k<zeta; k++)
        {
            f[i][j][k] = f[i-1][j][k];
        }
    }

    // for(int j=0; j<LY; j++) //right wall
    // {
    //     int i=(LX-2);
    //     double u_wX = u_x[i][j] + (1.0/2.0)*(u_x[i][j] - u_x[i-1][j]);
    //     double u_wY = u_y[i][j] + (1.0/2.0)*(u_y[i][j] - u_y[i-1][j]);
    //     // double u_wY = u_y[i][j] + (1.0/2.0)*(u_y[i][j] - u_y[i][j+1]);
    //     double uDotu = u_wX*u_wX + u_wY*u_wY;
    //     for(int k=0; k<zeta; k++)
    //     {
    //         double cDotu_w = (Cx[k])*(u_wX) + (Cy[k])*(u_wY);
    //         int tempK = invList[k];
    //         if(k==5 || k==1 || k==8)
    //         {
    //         f[i][j][tempK] = -ftemp[i][j][k] + 2.0*Wlist[k]*density[i+1][j]*(1 + 4.5*(cDotu_w*cDotu_w)- 1.5*(uDotu));
    //         }
    //     }
    // }


    // for(int i=1; i<LX-1; i++) //bottom wall (Bounceback)
    // {
    //     int j=1;
    //     for(int k=0; k<zeta; k++)
    //     {
    //         if(k==7 || k==7 || k==8)
    //         {
    //             int tempK = invList[k];
    //             f[i][j][k] = ftemp[i][j][tempK];
    //         }
    //     }
    // }

    //  for(int i=1; i<LX-1; i++) //Top wall (Bounceback)
    // {
    //     int j=(LY-2);
    //     for(int k=0; k<zeta; k++)
    //     {
    //         if(k==2 || k==5 || k==6)
    //         {
    //             int tempK = invList[k];
    //             f[i][j][k] = ftemp[i][j][tempK];
    //         }
    //     }
    // }
}   


void computeForce()
{
    double delmomentumX = 0.0;
    double delmomentumY = 0.0;
    for(int i=0; i<LX; i++)
    {
        for(int j=0; j<LY; j++)
        {
            if(obstacle[i][j]!=1)
            {
                for(int k=0; k<zeta; k++)
                {
                    int destX = i - Cx[k];
                    int destY = j - Cy[k];
                    int x = ((destX+LX)%LX);
                    int y = ((destY+LY)%LY);
                    int tempK = invList[k];
                    if(obstacle[x][y]==1)
                    {
                        delmomentumX += (f[i][j][k]+f[i][j][tempK])*Cx[k];
                        delmomentumY += (f[i][j][k]+f[i][j][tempK])*Cy[k];
                    }
                }
            }
        }
    }
    forceX = (delmomentumX)/(1);
    forceY = (delmomentumY)/(1); 
}


int main()
{
    system("rm -rf output");
    system("mkdir output");
    allocatememory();
    cylinder();

    for(int i=0; i<LX; i++) //initilization
    {
        for(int j=0; j<LY; j++)
        {
            if(obstacle[i][j]!=1)
            {
            u_x[i][j]=0.0;
            u_y[i][j]=0.0;
            density[i][j]=1.0;
            }
            // if(i==0)
            // {
            //     density[i][j]=1.0;
            // }
            // if(i==(LX-1))
            // {
            //     density[i][j]=1.0;
            // }
        }
    }

    for(int j=0; j<LY; j++)
    {
        int i=0;
        u_x[i][j]=u_inlet;
    }

    for(int i=1; i<LX; i++)
    {
        for(int j=0; j<LY; j++)
        {
            if(obstacle[i][j]!=1)
            {
            double udotu = u_x[i][j]*u_x[i][j] + u_y[i][j]*u_y[i][j];
            for(int k=0; k<zeta; k++)
            {
                double udotc = u_x[i][j]*Cx[k] + u_y[i][j]*Cy[k];
                f_eq[i][j][k] = Wlist[k]*density[i][j]*(1 + 3.0*(udotc) + 4.5*(udotc*udotc) - 1.5*(udotu));
            }
            }
        }
    }
    
    for(int i=1; i<LX; i++)
    {
        for(int j=0; j<LY; j++)
        {
            if(obstacle[i][j]!=1)
            {
            for(int k=0; k<zeta; k++)
            {
                f[i][j][k]=f_eq[i][j][k];
            }
            }
        }
    }

    for(int itr=0; itr<last_itr; itr++) //iteration starts
    {   
        for(int i=1; i<LX; i++)
        {
            for(int j=0; j<LY; j++)
            {
                if(obstacle[i][j]!=1)
                {
                double udotu = u_x[i][j]*u_x[i][j] + u_y[i][j]*u_y[i][j];
                for(int k=0; k<zeta; k++)
                {
                    double udotc = u_x[i][j]*Cx[k] + u_y[i][j]*Cy[k];
                    f_eq[i][j][k] = Wlist[k]*density[i][j]*(1 + 3.0*(udotc) + 4.5*(udotc*udotc) - 1.5*(udotu));
                }
                }
            }
        }

    #pragma omp parallel for collapse(2)
    for(int i=1; i<LX; i++) //collison
    {
        for(int j=0; j<LY; j++)
        {
            if(obstacle[i][j]!=1)
            {
            for(int k=0; k<zeta; k++)
            {
                ftemp[i][j][k] = f[i][j][k]*(1-tauinv) + f_eq[i][j][k]*(tauinv);
            }
            }
        }
    }

    #pragma omp parallel for collapse(2)
    for(int i=1; i<LX; i++) //streaming
    {
        for(int j=0; j<LY; j++)
        {
            if(obstacle[i][j]!=1)
            {
            for(int k=0; k<zeta; k++)
            {
                int destX = i - Cx[k];
                int destY = j - Cy[k];
                int x = ((destX+LX)%LX);
                int y = ((destY+LY)%LY);
                if(obstacle[x][y]!=1)
                {
                f[i][j][k] = ftemp[x][y][k];
                }
                else if(obstacle[x][y]==1)
                {
                    int tempK = invList[k];
                    f[i][j][k]=ftemp[i][j][tempK];
                }
            }
            }
        }
    }

    boundaryConditions();

    #pragma omp parallel for collapse(2)
    for(int i=1; i<LX; i++)
    {
        for(int j=0; j<LY; j++)
        {
            double tempDensity = 0.0;
            double tempUx = 0.0;
            double tempUy = 0.0;
            if(obstacle[i][j]!=1)
            {
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
    }



    computeForce();

    if(itr%save == 0)
    {
        sprintf(vel, "output/output_t%d.dat", itr);
        outputFile = fopen(vel, "w");
            for(int i=0; i<LX; i++)
            {
                for(int j=0; j<LY; j++)
                {
                    fprintf(outputFile, "%d %d %f %f %f %f %f %f\n", i,j,u_x[i][j],u_y[i][j],sqrt(pow(u_x[i][j],2) + pow(u_y[i][j],2)),density[i][j],forceX,forceY);
                }
            }
        fclose(outputFile);
    }

    cout<<"iteration ="<<itr<<endl;

    }




}


