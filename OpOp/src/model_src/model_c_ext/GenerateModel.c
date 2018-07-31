#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MT_random.h"
//#include "omp.h"


/*
Generate velocity using trial-error algorithm.
It generate a sphrical symmetric velocity distribution.
The random numbers are provided by the MT method.
*/



void v_gen(double const *particle_pot, int const *particle_pot_idx,  double const *particle_maxdf, double const *grid_pot, double const *df_grid, int const Npart, int const Ngrid , double *vx, double *vy, double *vz, double *v)
{
  double v1,v2,v3,v_tmp,u,norm;
  double df_max,pot,df_e,e,check_out;
  int check,max_indx,k;

  init_genrand(time(NULL));
  int i;


  //#pragma omp parallel  for private(i,k,check,v_tmp,u,v1,v2,v3,max_indx,df_max,pot,df_e,e,norm)
  for(i=0; i<Npart; i++)
  {
    v_tmp=2;
    check=-1;

    max_indx=particle_pot_idx[i]; //Index on the pot grid of the value of potmax=pot(R) of the particle
                                  //it is calculated in the python wrapper as searchsort, is the index k
                                  //where to insert the potmax in the grid pot to mantain the order.
                                  //e.g. potgrid=[1,3,4,7,8] potmax=6 m   max_indx=k=3 to have [1,3,4,6,7,8]
    df_max=particle_maxdf[i];
    pot=particle_pot[i];


    while( v_tmp>=1 || check<0 )
    {
        v1=2*genrand_real3()-1;
        v2=2*genrand_real3()-1;
        v3=2*genrand_real3()-1;

        v_tmp=sqrt(v1*v1+v2*v2+v3*v3);
        e=pot*(1-v_tmp*v_tmp);


    if(v_tmp<1)
    {

        //Now we want to find the df(e), first
        //we search the low and up values on the pot grid that bound the value e
        //Since e is always lower than potmax that in the grid is lower than the value potgrid[k] (see before)
        //we start the search from k. When we found a k where e<=pot_grid[k], we find the two values on the grid
        //that bound the value e, i.e.: pot_grid[k-1]<e<pot_grid[k].
        //Now we can interpolate in this interval to find f(e).
        k=max_indx;
        check_out=0;
        while(e<=grid_pot[k])
        {
            k=k+1;


			//If the cyclec reach k=nrt, we are out the grid, so we assume to be in the last point of the grid
            if ( k==(Ngrid) )
            {
            	break;
            	check_out=1;
            }



        }


        //Simply a 1D interpolation: (y-y1)/(y2-y1)=(x-x1)/(x2-x1)
        if(check_out==0)
        {
        	df_e= df_grid[k-1] + (e-grid_pot[k-1])*(df_grid[k] - df_grid[k-1])/(grid_pot[k] - grid_pot[k-1]);
        }
        else
        {
        //if outgrid, dfe always 0, i.e. always reject the value
        	df_e= 0;
        	//df_e= df_grid[k-1] + (e-grid_pot[k-1])*(df_grid[k-1])/(grid_pot[k-1]);
        	//df_e= df_grid[k-1] + (e-grid_pot[k-1])*(df_grid[k] - df_grid[k-1])/(grid_pot[k] - grid_pot[k-1]);
        }

        u=genrand_real3();



        if (df_e>=u*df_max)
        {
            check=1;
        }

    }


    }
    norm=sqrt(2*pot);

    vx[i]=v1*norm;
    vy[i]=v2*norm;
    vz[i]=v3*norm;
    v[i]=v_tmp*norm;

   }


}




/*
//QUesto qui sotto non funziona jamais, MISTERO
void v_gen_alt(double const *r_grid,double const *pot_grid, double const *m_grid, double const *df_grid, int const nrt, int const npart,double *r,double *vx, double *vy, double *vz, double *pot, double *outk)
{
   double u,mnorm,xr,du,psit,v1,v2,v3,pp,v,ff,dpsi,qv,fc;
   int k,i,c,j,ii,check,w,iter;



   for (w=0; w<npart; w++)
   {


       //Assegna raggio
       k=0;

           u=genrand_real3();
           mnorm=m_grid[nrt-1]; //maxmass


           //assign mass grid point
           for(i=1;i<nrt;i++)
           {
                if(u>m_grid[i]/mnorm)
                {
                    k=i;
                }
           }
           xr=r_grid[k] + (u*mnorm-m_grid[k])*(r_grid[k+1]-r_grid[k])/(m_grid[k+1]-m_grid[k]);


       r[w]=xr;
       du=(xr-r_grid[k])/(r_grid[k+1]-r_grid[k]);






      //Assegna potenziale
      if(k==nrt-1)
      {
        psit=pot_grid[k];
      }
      else
      {
        psit=pot_grid[k]+(pot_grid[k+1]-pot_grid[k])*du;
      }

      pot[w]=psit;


      v=2;
      check=-1;
      iter=0;
      while( (v>1) || (check<0))
      {
        v1=2*genrand_real3()-1;
        v2=2*genrand_real3()-1;
        v3=2*genrand_real3()-1;

        v=v1*v1+v2*v2+v3*v3;

      if (v<=1)
      {


              pp=genrand_real3();

              //dpsi=(psit-pot_grid[k])/(pot_grid[k+1] - pot_grid[k]);

              //pp=df_grid[k]*pp;
              if(k==nrt-1)
              {
               fc=df_grid[k];
              }
              else
              {
                if(k>0 && k<nrt-1)
                {
                    fc=( df_grid[k] + (df_grid[k+1] - df_grid[k])*du);
                }
                else
                {
                  fc=0.5*(df_grid[k+1] - (df_grid[k+2] - df_grid[k+1]) ) * (1-du);
                }
              }

              pp=pp*fc;

              qv=psit*(1-v);
              j=0;
              for (ii=k; i<nrt; i++)
              {
                if(pot_grid[ii]>=qv)
                {
                    j=ii;
                }
              }

              if (j!=(nrt-1))
              {
                dpsi=(qv - pot_grid[j])/( pot_grid[j+1] - pot_grid[j] );
                ff=df_grid[j] + (df_grid[j+1] -df_grid[j])*dpsi;
              }
              else
              {
                ff=df_grid[j];
              }

              //ff=df_grid[j];
              //ff= df_grid[j-1] + (qv-pot_grid[j-1])*(df_grid[j] - df_grid[j-1])/(pot_grid[j] - pot_grid[j-1]);
              //dpsi=(qv-pot_grid[j])/(pot_grid[j+1] - pot_grid[j]);
              //ff=df_grid[j] + (df_grid[j+1] - df_grid[j])*dpsi;


              if(pp>ff)
              {
                check=-1;
              }
              else
              {
                check=1;
              }
      }

    iter=iter+1;
    }
    vx[w]=v1;
    vy[w]=v2;
    vz[w]=v3;
    outk[w]=iter;



}
}
*/