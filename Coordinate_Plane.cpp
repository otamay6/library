#include<iostream>
#include<math.h>
constexpr double eps=1e-9;
const double PI=acos(-1);

struct Coordinate_Plane{
    double x,y;
    Coordinate_Plane(double X=0,double Y=0):x(X),y(Y){}
    Coordinate_Plane(const Coordinate_Plane &cp){
        x=cp.x;y=cp.y;
    }
    Coordinate_Plane& set(double X,double Y){
        x=X;y=Y;
        return *this;
    }
    Coordinate_Plane& operator=(const Coordinate_Plane &p){return set(p.x , p.y);}
    Coordinate_Plane operator+(const Coordinate_Plane &p)const{return Coordinate_Plane(x+p.x , y+p.y);}
    Coordinate_Plane operator-(const Coordinate_Plane &p)const{return Coordinate_Plane(x-p.x , y-p.y);}
    Coordinate_Plane operator*(double s)const{return Coordinate_Plane(x*s,y*s);}
    Coordinate_Plane operator/(double s)const{return Coordinate_Plane(x/s,y/s);}
    friend Coordinate_Plane operator*(double s,const Coordinate_Plane &r){return Coordinate_Plane(r.x*s,r.y*s);}
    friend Coordinate_Plane operator/(double s,const Coordinate_Plane &r){return Coordinate_Plane(r.x/s,r.y/s);}
    double operator*(Coordinate_Plane p)const{return dot(p);}
    double operator^(Coordinate_Plane p)const{return cross(p);}
    Coordinate_Plane& operator+=(const Coordinate_Plane &p){return set(x+p.x , y+p.y);}
    Coordinate_Plane& operator-=(const Coordinate_Plane &p){return set(x-p.x , y-p.y);}
    Coordinate_Plane& operator*=(double s){return set(x*s,y*s);}
    Coordinate_Plane& operator/=(double s){return set(x/s,y/s);}
    double dot(const Coordinate_Plane &p)const{return x*p.x + y*p.y;} 
    double cross(const Coordinate_Plane &p)const{return x*p.y-p.x*y;}
    int ort()const{
        if(fabs(x)<eps && fabs(y)<eps) return 0;
        if(y>0) return (x>0)?1:2;
        else return (x<=0)?3:4;
    }
    bool operator<(const Coordinate_Plane &p)const{
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
    friend std::istream& operator>>(std::istream &is, Coordinate_Plane &x){double valX,valY; is>>valX>>valY; x.set(valX,valY); return is;}
    friend std::ostream& operator<<(std::ostream &os, const Coordinate_Plane &v){ printf("%.6f %.6f",v.x,v.y); return os; }
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
        // no reach
        return 0.0;
    }
    void rotate(double theta){
        double X=x,Y=y;
        x=cos(theta)*X-sin(theta)*Y;
        y=sin(theta)*X+cos(theta)*Y;
    }
};
struct Line{
    Coordinate_Plane A,B;
    Line(Coordinate_Plane p1=Coordinate_Plane(0,0),Coordinate_Plane p2=Coordinate_Plane(1,1)):A(p1),B(p2){}
    double len()const{return (B-A).norm();}
    Coordinate_Plane shortest_point(const Coordinate_Plane &x)const{
        double a=B.x-A.x,b=B.y-A.y;
        double t=-(a*(A.x-x.x)+b*(A.y-x.y))/(a*a+b*b);
        return Coordinate_Plane(a*t+A.x,b*t+A.y);
    }
    double dist(const Coordinate_Plane &x)const{
        return (x-shortest_point(x)).norm();
    }
    bool online(const Coordinate_Plane &x)const{
        return fabs((B-A)^(x-A))<eps&&(x-A).norm()<len()+eps&&(x-B).norm()<len()+eps;
    }
    bool intersect(const Coordinate_Plane &x,const Coordinate_Plane &y){
        return ((x.x-y.x)*(A.y-x.y) + (x.y-y.y)*(x.x-A.x))*
                ((x.x-y.x)*(B.y-x.y) + (x.y-y.y)*(x.x-B.x)) < 0;
    }

    bool Lintersect(const Line &L)const{
        return  ((L.A.x-L.B.x)*(A.y-L.A.y) + (L.A.y-L.B.y)*(L.A.x-A.x))*
                ((L.A.x-L.B.x)*(B.y-L.A.y) + (L.A.y-L.B.y)*(L.A.x-B.x))<0&&
                ((A.x-B.x)*(L.A.y-A.y) + (A.y-B.y)*(A.x-L.A.x))*
                ((A.x-B.x)*(L.B.y-A.y) + (A.y-B.y)*(A.x-L.B.x))<0;
    }
    friend std::istream& operator>>(std::istream &is, Line &x){Coordinate_Plane X,Y; is>>X>>Y;x.A=X;x.B=Y;return is;}
    friend std::ostream& operator<<(std::ostream &os, const Line &x){ os << x.A <<" "<< x.B; return os; }
};
struct Triangle{
    Coordinate_Plane a,b,c;
    Triangle(Coordinate_Plane p1=Coordinate_Plane(0,0),Coordinate_Plane p2=Coordinate_Plane(1,0), Coordinate_Plane p3=Coordinate_Plane(0,1)):a(p1),b(p2),c(p3){}

    double A()const{double res=fabs((b-a).rad()-(c-a).rad());if(res>=PI) res-=PI;return res;}
    double B()const{double res=fabs((a-b).rad()-(c-b).rad());if(res>=PI) res-=PI;return res;}
    double C()const{double res=fabs((a-c).rad()-(b-c).rad());if(res>=PI) res-=PI;return res;}
    Line AB()const{return Line(a,b);}
    Line BC()const{return Line(b,c);}
    Line CA()const{return Line(c,a);}

    Coordinate_Plane G()const{return (a+b+c)/3;}
    Coordinate_Plane O()const{return (sin(2*A())*a+sin(2*B())*b+sin(2*C())*c)/(sin(2*A())+sin(2*B())+sin(2*C()));}
    double R()const{return a.norm()/(2*sin(A()));}
    Coordinate_Plane I()const{return (a.norm()*a+b.norm()*b+c.norm()*c)/(a.norm()+b.norm()+c.norm());}
    double r()const{return 2*area()/((a-b).norm()+(b-c).norm()+(c-a).norm());}
    Coordinate_Plane H()const{return (tan(A())*a+tan(B())*b+tan(C())*c)/(tan(A())+tan(B())+tan(C()));}
    double AH()const{return 2*R()*cos(A());}
    double BH()const{return 2*R()*cos(B());}
    double CH()const{return 2*R()*cos(C());}
    double area()const{return ((b-a)^(c-a))/2;}

    bool inside(const Coordinate_Plane &x)const{
        return ((b-a)^(x-a))*((c-a)^(x-a))<0&&
               ((a-b)^(x-b))*((c-b)^(x-b))<0&&
               ((a-c)^(x-c))*((b-c)^(x-c))<0;
    }
    bool online(const Coordinate_Plane &x)const{return AB().online(x)||BC().online(x)||CA().online(x);}
    bool outside(const Coordinate_Plane &x)const{return (!inside(x)&&!online(x));}

    friend std::istream& operator>>(std::istream &is, Triangle &x){Coordinate_Plane X,Y,Z; is>>X>>Y>>Z; x.a=X;x.b=Y;x.c=Z;return is;}
    friend std::ostream& operator<<(std::ostream &os, const Triangle &v){ os << v.a << v.b; return os; }
};
struct Circle{
    Coordinate_Plane center;
    double radius;
    Circle(double r=0,Coordinate_Plane c={0,0}):radius(r),center(c){}
    bool inside(const Coordinate_Plane &C)const{
        return (C-center).norm2()<radius*radius+eps;
    }
    bool intersect(const Circle &C)const{
        return (C.center-center).norm2()<(C.radius+radius)*(C.radius+radius)+eps;
    }
    std::pair<Coordinate_Plane,Coordinate_Plane> interP(const Circle &C)const{
        Coordinate_Plane Q=center-C.center;
        double x1=Q.x,y1=Q.y,r1=radius,r2=C.radius;
        double a=(x1*x1+y1*y1+r1*r1-r2*r2)/2.0,t=Q.norm2();
        Coordinate_Plane A((a*x1+y1*sqrt(t*r1*r1-a*a))/t,(a*y1-x1*sqrt(t*r1*r1-a*a))/t);
        Coordinate_Plane B((a*x1-y1*sqrt(t*r1*r1-a*a))/t,(a*y1+x1*sqrt(t*r1*r1-a*a))/t);
        A+=C.center;
        B+=C.center;
        return std::make_pair(A,B);
    }
    friend std::ostream& operator<<(std::ostream &os, const Circle &v){ os << v.center<<" ";printf("%.06f",v.radius); return os; }
};