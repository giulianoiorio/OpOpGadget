
void evalmass(double const *R, double const *dens, int const len, double *mass)
{
    for (int i=0; i<len; i++)
    {
        if (i==0)
            {
                mass[i]=R[i]*R[i]*R[i]*dens[i]/3.;
            }

        else
            {

                mass[i]=mass[i-1]+0.5*(R[i]-R[i-1])*(dens[i]*R[i]*R[i] + dens[i-1]*R[i-1]*R[i-1] );
            }

    }

}


void evalpot(double const *R, double const *dens, double const *mass  ,int const len, double *pot)
{
    int k;
    int i;

    for(k=0; k<len; k++)
    {
        pot[k]=0;
        for(i=k; i<len; i++)
        {
            if(i==0)
            {
                pot[k]=pot[k] + 0.5*(R[i+1]-R[i])*(dens[i+1]*R[i+1]+dens[i]*R[i]);
                //pot[k]=pot[k] + (R[i]*R[i]) * (dens[k]+(dens[k] - dens[k+1])*(0.5*R[k]/(R[k+1]-R[k])));
            }
            else
            {
                pot[k]=pot[k] + 0.5*(R[i]-R[i-1])*(dens[i]*R[i]+dens[i-1]*R[i-1]);
            }
        }

        pot[k]=pot[k]+mass[k]/R[k];
    }
}


