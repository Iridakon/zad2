#include <iostream>
#include "windows.h"
#include "math.h"

double f(double x){
    return ((x*sin(x))*(x*sin(x)));
}
double Ans(double a, double b)
{
    return abs(a-b);
}
double LPramoug(double left, double right, int N)//Метод левых прямоугольников
{
    double sum=0;
    double shag=(right-left)/N;
    double pos=left;
    for (int i = 0; i < N; ++i) {
        sum+=f(pos)*shag;

        pos+=shag;

    }
    return sum;
}
double CPramoug(double left, double right, int N)//метод центральных пр-в
{
    double sum=0;
    double shag=(right-left)/N;
    double pos=left+shag/2;
    for (int i = 0; i < N; ++i) {
        sum+=f(pos)*shag;
        //printf("%f %f \n",pos,sum);
        pos+=shag;
    }
    return sum;
}
double Trap(double left, double right, int N)//метод трапеций
{
    double sum=0;
    double shag=(right-left)/N;
    double pos=left;
    for (int i = 0; i < N; ++i) {
        sum+=(f(pos)+f(pos+shag))/2*shag;
        //printf("%f %f \n",pos,pos+shag);
        pos+=shag;
    }
    return sum;
}
double Simp(double left, double right, int N)//метод Симnсона
{
    double sum=0;
    double shag=(right-left)/N;
    double pos=left;
    for (int i = 0; i < N; ++i) {
        sum+=shag/6*(f(pos)+4*f((pos+pos+shag)/2)+f(pos+shag));
        pos+=shag;
    }
    return sum;
}
double Gauss2(double left, double right, int N)//метод Симnсона
{
    double sum=0;
    double shag=(right-left)/N;
    double pos=left;
    double kor,ab;
    for (int i = 0; i < N; ++i) {
        kor=shag/(2* sqrt(3));
        ab=pos+shag/2;
        sum+=(f(ab+kor)+f(ab-kor))*shag/2;
        pos+=shag;
        //printf("%f %f\n",pos,sum);
    }
    return sum;
}
double Gauss3(double left, double right, int N)//метод Симnсона
{
    double sum=0;
    double shag=(right-left)/N;
    double pos=left;
    double kor,ab;
    for (int i = 0; i < N; ++i) {
        kor=shag/2* sqrt(0.6);
        ab=pos+shag/2;
        sum+=((f(ab-kor)+f(ab+kor))*((float)5/9.0)+f(ab)*((float)8/9))*shag/2;
        pos+=shag;
        //printf("%f %f\n",f(ab+kor),f(ab+kor)*((float)5/9));
    }
    return sum;
}
double Better(double In,double In2,int p)
{
    double alpha=1/(1-pow(2,p));
    return abs(alpha*In+(1-alpha)*In2);
}
int main() {
    SetConsoleOutputCP(CP_UTF8);
    double a=0;
    double b=M_PI/2;
    double ans=pow(M_PI,3)/48+M_PI/8;
    int N=10;
    //Trap(a,b,N);

    printf("Кол-во узлов   \t%i\t\t%i \t\tИзменение \tУлучш\n",N,N*2);
    double y1,y2;
    /*
    y1=Ans(LPramoug(a,b,N),ans);
    y2=Ans(LPramoug(a,b,2*N),ans);
    printf("Левые пр        %8.12f   %8.12f  %5.2f   %8.12f\n",y1,y2,y1/y2, Better(y1,y2,1));
    y1=Ans(CPramoug(a,b,N),ans);
    y2=Ans(CPramoug(a,b,N*2),ans);
    printf("Центральные пр  %8.12f   %8.12f  %5.2f   %8.12f\n",y1,y2,y1/y2,Better(y1,y2,2));
    y1=Ans(Trap(a,b,N),ans);
    y2=Ans(Trap(a,b,N*2),ans);
    printf("Трапеции        %8.12f   %8.12f  %5.2f   %8.12f\n",y1,y2,y1/y2,Better(y1,y2,2));
    y1=Ans(Simp(a,b,N),ans);
    y2=Ans(Simp(a,b,N*2),ans);
    printf("Сипсон          %8.12f   %8.12f  %5.2f   %8.12f\n",y1,y2,y1/y2,Better(y1,y2,4));
     */
    y1=Ans(Gauss2(a,b,N),ans);
    y2=Ans(Gauss2(a,b,N*2),ans);
    printf("Гаусс по 2      %8.12f   %8.12f  %5.2f   %8.12f\n",y1,y2,y1/y2,Better(y1,y2,4));
    y1=Ans(Gauss3(a,b,N),ans);
    y2=Ans(Gauss3(a,b,N*2),ans);
    printf("Гаусс по 3      %8.12f   %8.12f  %5.2f   %8.12f\n",y1,y2,y1/y2,Better(y1,y2,1));

    return 0;
}
