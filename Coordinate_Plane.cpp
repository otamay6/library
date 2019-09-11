#include<iostream>
#include<math.h>
constexpr double eps=1e-9;
const double PI=acos(-1);

struct GP{
    double x,y;
    GP(double X=0,double Y=0):x(X),y(Y){}
    GP(const GP &cp){
        x=cp.x;y=cp.y;
    }
    GP& set(double X,double Y){
        x=X;y=Y;
        return *this;
    }
    GP& operator=(const GP &p){return set(p.x , p.y);}
    GP operator+(const GP &p)const{return GP(x+p.x , y+p.y);}
    GP operator-(const GP &p)const{return GP(x-p.x , y-p.y);}
    GP operator*(double s)const{return GP(x*s,y*s);}
    GP operator/(double s)const{return GP(x/s,y/s);}
    friend GP operator*(double s,const GP &r){return GP(r.x*s,r.y*s);}
    friend GP operator/(double s,const GP &r){return GP(r.x/s,r.y/s);}
    double operator*(GP p)const{return dot(p);}
    double operator^(GP p)const{return cross(p);}
    GP& operator+=(const GP &p){return set(x+p.x , y+p.y);}
    GP& operator-=(const GP &p){return set(x-p.x , y-p.y);}
    GP& operator*=(double s){return set(x*s,y*s);}
    GP& operator/=(double s){return set(x/s,y/s);}
    double dot(const GP &p)const{return x*p.x + y*p.y;} 
    double cross(const GP &p)const{return x*p.y-p.x*y;}
    int ort()const{
        if(fabs(x)<eps && fabs(y)<eps) return 0;
        if(y>0) return (x>0)?1:2;
        else return (x<=0)?3:4;
    }
    bool operator<(const GP &p)const{
        int o=ort(),po=p.ort();
        if(o!=po) return o<po;
        else{
            double cr=cross(p);
            if(cr==0){
                double d=norm(),pd=p.norm();
                return d<pd;
            }
            return cr>0;
        }
    }
    friend std::istream& operator>>(std::istream &is, GP &x){double valX,valY; is>>valX>>valY; x.set(valX,valY); return is;}
    friend std::ostream& operator<<(std::ostream &os, const GP &v){ os << v.x<<" "<<v.y<<endl; return os; }
    double norm2()const{return x*x+y*y;}
    double norm()const{return sqrt(norm2());}
    double rad()const{
        if(x==0){
            if(y>0) return PI/2.0;
            if(y==0) return -1.0;
            if(y<0) return -PI/2.0;
        }else{
            if(y==0){
                if(x<0) return PI;
                if(x>0) return 0.0;
            }
            else return atan2(y,x);
        }
    }
    void rotate(double theta){
        double X=x,Y=y;
        x=cos(theta)*X-sin(theta)*Y;
        y=sin(theta)*X+cos(theta)*Y;
    }
};
struct Line{
    GP A,B;
    Line(GP p1=GP(0,0),GP p2=GP(1,1)):A(p1),B(p2){}
    double len()const{return (B-A).norm();}
    GP shortest_point(const GP &x)const{
        double a=B.x-A.x,b=B.y-A.y;
        double t=-(a*(A.x-x.x)+b*(A.y-x.y))/(a*a+b*b);
        return GP(a*t+A.x,b*t+A.y);
    }
    double dist(const GP &x)const{
        return (x-shortest_point(x)).norm();
    }
    bool online(const GP &x)const{
        return ((B-A)^(x-A))==0&&(x-A).norm()<=len()&&(x-B).norm()<=len();
    }
    
    friend std::istream& operator>>(std::istream &is, Line &x){GP X,Y; is>>X>>Y;x.A=X;x.B=Y;return is;}
    friend std::ostream& operator<<(std::ostream &os, const Line &x){ os << x.A << x.B; return os; }
};
struct Triangle{
    GP a,b,c;
    Triangle(GP p1=GP(0,0),GP p2=GP(1,0),GP p3=(0,1)):a(p1),b(p2),c(p3){}

    double A()const{double res=fabs((b-a).rad()-(c-a).rad());if(res>=PI) res-=PI;return res;}
    double B()const{double res=fabs((a-b).rad()-(c-b).rad());if(res>=PI) res-=PI;return res;}
    double C()const{double res=fabs((a-c).rad()-(b-c).rad());if(res>=PI) res-=PI;return res;}
    Line AB()const{return Line(a,b);}
    Line BC()const{return Line(b,c);}
    Line CA()const{return Line(c,a);}

    GP G()const{return (a+b+c)/3;}
    GP O()const{return (sin(2*A())*a+sin(2*B())*b+sin(2*C())*c)/(sin(2*A())+sin(2*B())+sin(2*C()));}
    double R()const{return a.norm()/(2*sin(A()));}
    GP I()const{return (a.norm()*a+b.norm()*b+c.norm()*c)/(a.norm()+b.norm()+c.norm());}
    double r()const{return 2*area()/((a-b).norm()+(b-c).norm()+(c-a).norm());}
    GP H()const{return (tan(A())*a+tan(B())*b+tan(C())*c)/(tan(A())+tan(B())+tan(C()));}
    double AH()const{return 2*R()*cos(A());}
    double BH()const{return 2*R()*cos(B());}
    double CH()const{return 2*R()*cos(C());}
    double area()const{return ((b-a)^(c-a))/2;}

    bool inside(const GP &x)const{
        return ((b-a)^(x-a))*((c-a)^(x-a))<0&&
               ((a-b)^(x-b))*((c-b)^(x-b))<0&&
               ((a-c)^(x-c))*((b-c)^(x-c))<0;
    }
    bool online(const GP &x)const{return AB().online(x)||BC().online(x)||CA().online(x);}
    bool outside(const GP &x)const{return (!inside(x)&&!online(x));}

    friend std::istream& operator>>(std::istream &is, Triangle &x){GP X,Y,Z; is>>X>>Y>>Z; x.a=X;x.b=Y;x.c=Z;return is;}
    friend std::ostream& operator<<(std::ostream &os, const Triangle &v){ os << v.a << v.b; return os; }
};