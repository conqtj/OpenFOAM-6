#include"spectralMaths.H"

#define PI 3.14159265359

void Foam::quadgen(int Nq, scalarList& w, scalarList& g)
{
    scalarList gg(2*Nq);
    scalarList ww(2*Nq);

    gausscheb2(gg,ww,2*Nq);

    for (int i=0;i<Nq;i++)
    {
        (g)[i]=-gg[i+Nq];
        (w)[i]=ww[i+Nq];
    }
    scalar sum=0.;

    for(int i=0;i<Nq;i++)
    {
        sum+=(w)[i]; 
    }
   
    for(int i=0;i<Nq;i++) (w)[i]/=sum;

    return ;
}

void Foam::gausscheb2(scalarList& x, scalarList& wq, int n)
{
    scalar sum,theta;

    for(int k= 1;k<=(n+1)/2;k++)
    {
        theta=scalar(k)*PI/(n+1.0);
        x[k-1]= cos(theta);
        sum= 0.;
        for(int m=1;m<=(n+1)/2;m++)
        {
            sum= sum+sin((2*m-1)*theta)/(2.0*m-1.0);
        }
        wq[k-1]= 4.0*sin(theta)*sum/(n+1.0);
        x[n-k]=-x[k-1];
        wq[n-k]= wq[k-1];
    }

    return ;
}


Foam::scalarList Foam::kPowerLaw(scalar kmin,scalar kmax, scalar pwr,int n)
{
    /*! function generate a list of k values between kmin and kmax
    ! according to power law with power "pwr"
    */

    scalarList k(n);
    scalar pwrk_min, pwrk_max, pwrk_step;

    pwrk_min = pow(kmin,pwr);
    pwrk_max = pow(kmax,pwr);
    pwrk_step = (pwrk_max-pwrk_min)/(n-1);
    for (int i=0;i<n;i++)
    {
        k[i] = pwrk_min+scalar(i-1)*pwrk_step;
        k[i] = pow(k[i],1.0/pwr);
    }

    return k;
}


Foam::scalarList Foam::linearInterpMono
(int nop,const scalarList & xx,const scalarList & yy, int ni,const scalarList& xi)
{
    /* simple linear interpolation using closest points
    constant value for extrapolation
    assume xx and xi are monotonically increasing
    (1) It is assumed that xx is in ascendant order, which is typical for
    linear interpolation
    (2) xi is also assumed to be in ascendant order, which improves search
    speed
    (3) Constant extrapolation using two boundary values.
    */


    scalarList yi(ni);
    int n;
    n = 0;
    for (int i=0;i<ni;i++)
    {  
        while(xi[i]>=xx[n])
        {
            // xi(i) is not in interval xx(n-1) to xx(n)
            n = n + 1 ;//! then move to next interval

            if(n > nop-1) 
            {    
                //then ! xi(i) is larger than the largest xx
                // out of bound and use constant value
                for(int j=i;j<ni;j++)yi[j]= yy[nop-1]; 
                return yi; //exit loopi ! no need to perform further calculation
            }
        }  
        // now xi[i] is between xx[n-1] and xx[n], except xi[i] is smaller than xx[1]
        if(n==0) 
        { 
            // xi[i] is smaller than the smallest of xx
            yi[i] = yy[0] ;
        }
        else 
        {
            yi[i] = yy[n]+(yy[n-1]-yy[n])*(xi[i]-xx[n])/(xx[n-1]-xx[n]);
        } 
    }
   
    return yi; 
}
