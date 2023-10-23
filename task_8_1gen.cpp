#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double m_1=0; //choose variant of approximation 
double m_2=1;


double real=1; //evident variant 
double image=0; //unevident vatiant


double u; //current value of u for method_1
double v; //current value of v for method_1

double old_u[]={0,0,0,0,0}; // store previous iteration of u
double old_v[]={0,0,0,0,0}; // store previous iteration of v

const double l_edge=0; // left edge for t
const double r_edge=5; // right edge for t
const int N=100000;      // total amount of samples
const double step=(r_edge-l_edge)/N; // iteration step for t
bool st_flag=1; // flag meaning that calculating can be started
const double u_0=1.0;     // initial condition for u: u(t=l_edge)=u_0
const double v_0=1.0;     // initial condition for v: v(t=l_edge)=v_0
double t=l_edge; // current value of t


double a=998;  //u'=a*u+b*v 
double b=1998; //v'=c*u+d*v
double c=-999;
double d=-1999;

double true_u(double s)
{
    static double C_0=-u_0-v_0;
    static double C_1=2*v_0+u_0;

    return -2*C_0*exp(-s)-C_1*exp(-1000*s);


}
double true_v(double s)
{
    static double C_0=-u_0-v_0;
    static double C_1=2*v_0+u_0;

    return C_0*exp(-s)+C_1*exp(-1000*s);
}




double string_1(double u, double v) // u'= 998*u+1998*v
{
    return 998*u+1998*v;
}

double string_2(double u, double v) // v'=-999*u-1999*v
{
    return -999*u-1999*v;
}

double initialziation(int id) //use the simplest explicit method to determine the initial values for the advanced methods
{
    old_u[id-1]=old_u[id]+step*string_1(old_u[id],old_v[id]);
    old_v[id-1]=old_v[id]+step*string_2(old_u[id],old_v[id]);

    return 0;
}



int method_1(FILE* out)
{




if(m_2)
{

int max_id=sizeof(old_u)/sizeof(double)-1;
if (st_flag)
{
    old_u[max_id-1]=u_0;
    old_v[max_id-1]=v_0;
    for(int j=(max_id-1);j>0;j--) 
    {
        initialziation(j);
        fprintf(out,"%lf %lf %lf %lf %lf\n",t,old_u[j],old_v[j], true_u(t), true_v(t));
        t+=step;
    }
}
else
{



for (int i=(sizeof(old_u)/sizeof(double)-1);i>0;i--)
{
old_u[i]=old_u[i-1];
old_v[i]=old_v[i-1];
}


if (real)
{
    u=old_u[1]+step/24*(55*string_1(old_u[1],old_v[1])-59*string_1(old_u[2],old_v[2])+37*string_1(old_u[3],old_v[3])-9*string_1(old_u[4],old_v[4]));
    v=old_v[1]+step/24*(55*string_2(old_u[1],old_v[1])-59*string_2(old_u[2],old_v[2])+37*string_2(old_u[3],old_v[3])-9*string_2(old_u[4],old_v[4]));
}


if (image)
{
    double str_1=old_u[1]+(19*(a*old_u[1]+b*old_v[1])-5*(a*old_u[2]+b*old_v[2])+(a*old_u[3]+b*old_v[3]))*step/24;
    double str_2=old_v[1]+(19*(c*old_u[1]+d*old_v[1])-5*(c*old_u[2]+d*old_v[2])+(c*old_u[3]+d*old_v[3]))*step/24;
    double k=1-9/24*step*a;
    double l=step*b*9/24;
    double m=1-9/24*step*d;
    double n=step*c*9/24;

    u=(n/m*str_1/k+str_2/m)/(1-n*l/(m*k));
    v=(l/k*str_2/m+str_1/k)/(1-n*l/(m*k));
    //printf("%lf %lf\n",u,v);
}

old_u[0]=u;
old_v[0]=v;
}
}

if (m_1)
{
if (st_flag) 
    
    {
        u=u_0;
        v=v_0;
        fprintf(out,"%lf %lf %lf %lf %lf %lf %lf\n",t,u,v, true_u(t), true_v(t), log(abs(true_u(t)-u)), log(abs(true_v(t))-v));
        t+=step;
    }

    else
    {

double old_u=u;
double old_v=v;

if (real)
{
    u=old_u+step*string_1(old_u,old_v);
    //printf("%lf\n",u);
    v=old_v+step*string_2(old_u,old_v);
}

if (image)
{
    u=old_u+string_1(old_u,old_v)*step/2+(b*step/2)/(1-d*step/2)*(old_v+string_2(old_u,old_v)*step/2);
    u/=1-a*step/2-b*step/2*c*step/2/(1-d*step/2);
    v=old_v+string_2(old_u,old_v)*step/2+(c*step/2)/(1-a*step/2)*(old_u+string_1(old_u,old_v)*step/2);
    v/=1-d*step/2-c*step/2*b*step/2/(1-a*step/2);
}
    
    
    }
}


if (!st_flag)
{
fprintf(out,"%lf %lf %lf %lf %lf %lf %lf\n",t,u,v, true_u(t), true_v(t), log(abs(true_u(t)-u)), log(abs(true_v(t))-v));
t+=step;    
}

else  st_flag=0; 
    return 0;
}






int main()
{
    double b=1.0;
    int g=sizeof(b);
    FILE* out_1=fopen("m_1_r.txt","w");
    FILE* out_2=fopen("m_1_i.txt","w");
    FILE* out_3=fopen("m_2_r.txt","w");
    FILE* out_4=fopen("m_2_i.txt","w");

    if (out_1==NULL || out_2==NULL )
        {
            printf("cannot open or create output file");
            return 1;
        }

    double flag_m1[]={1,1,0,0};
    double flag_m2[]={0,0,1,1};
    double flag_r[]={1,0,1,0};
    double flag_i[]={0,1,0,1};
    FILE* out[]={out_1,out_2,out_3,out_4};

    for (int i=0;i<4;i++)
    {
        m_1=flag_m1[i];
        m_2=flag_m2[i];
        real=flag_r[i];
        image=flag_i[i];

        for (int j=0; j<N; j++) method_1(out[i]);
        t=l_edge;
        st_flag=1;
    }

    fclose(out_1);
    fclose(out_2);
    fclose(out_3);
    fclose(out_4);
    system("pause");
    return 0;
}

//plot "m_1_r.txt" u 1:4 w l, "" u 1:5 w l,"m_2_r.txt" u 1:4 w l, "" u 1:5 w l,"m_1_i.txt" u 1:4 w l, "" u 1:5 w l