#include<iostream>
#include<math.h>

#include "common.h"

/// @brief 2次元座標平面上の点、ベクトル、直線
/// @details CoordinatePonitは点、ベクトルを表し、Lineは直線
/// @note CoordinatePonitの内積は*、外積は^で表現される
/// @note CoordinatePonitのort()はベクトルの向きを1~4で表現する
struct CoordinatePonit{
    double x,y;
    CoordinatePonit(double X=0,double Y=0):x(X),y(Y){}
    CoordinatePonit(const CoordinatePonit &cp){
        x=cp.x;y=cp.y;
    }
    CoordinatePonit& set(double X,double Y){
        x=X;y=Y;
        return *this;
    }
    CoordinatePonit& operator=(const CoordinatePonit &p){return set(p.x , p.y);}
    /// @brief ベクトルとベクトルの和、差
    CoordinatePonit operator+(const CoordinatePonit &p)const{return CoordinatePonit(x+p.x , y+p.y);}
    CoordinatePonit operator-(const CoordinatePonit &p)const{return CoordinatePonit(x-p.x , y-p.y);}
    /// @brief ベクトルとスカラーの積、商
    CoordinatePonit operator*(double s)const{return CoordinatePonit(x*s,y*s);}
    CoordinatePonit operator/(double s)const{return CoordinatePonit(x/s,y/s);}
    /// @brief スカラーとベクトルの積
    friend CoordinatePonit operator*(double s,const CoordinatePonit &r){return CoordinatePonit(r.x*s,r.y*s);}
    /// @brief ベクトルとベクトルの内積、外積
    double operator*(CoordinatePonit p)const{return dot(p);}
    double operator^(CoordinatePonit p)const{return cross(p);}
    /// @brief ベクトルとベクトルの和、差、スカラー倍、商の代入演算子
    CoordinatePonit& operator+=(const CoordinatePonit &p){return set(x+p.x , y+p.y);}
    CoordinatePonit& operator-=(const CoordinatePonit &p){return set(x-p.x , y-p.y);}
    CoordinatePonit& operator*=(double s){return set(x*s,y*s);}
    CoordinatePonit& operator/=(double s){return set(x/s,y/s);}
    double dot(const CoordinatePonit &p)const{return x*p.x + y*p.y;} 
    double cross(const CoordinatePonit &p)const{return x*p.y-p.x*y;}
    /// @brief 点の位置を1~4で表現する
    /// @return 1:第一象限 2:第二象限 3:第三象限 4:第四象限 0:原点
    int ort()const{
        if(fabs(x)<mathlib::eps && fabs(y)<mathlib::eps) return 0;
        if(y>0) return (x>0)?1:2;
        else return (x<=0)?3:4;
    }
    /// @brief ベクトルの向きを比較する
    /// @param p 比較対象のベクトル
    /// @return thisの方がpより反時計回りにあるならtrue
    bool operator<(const CoordinatePonit &p)const{
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
    friend std::istream& operator>>(std::istream &is, CoordinatePonit &x){double valX,valY; is>>valX>>valY; x.set(valX,valY); return is;}
    friend std::ostream& operator<<(std::ostream &os, const CoordinatePonit &v){ printf("%.6f %.6f",v.x,v.y); return os; }
    /// @brief ベクトルの大きさ、ノルム
    /// @return ノルムの2乗
    double norm2()const{return x*x+y*y;}
    /// @brief ベクトルの大きさ、ノルム
    /// @return ノルム
    double norm()const{return sqrt(norm2());}
    /// @brief ベクトルの角度を求める
    /// @return x軸からの反時計回りの角度、ラジアン
    double rad()const{
        if(x==0){
            if(y>0) return mathlib::PI/2.0;
            if(y==0) return 0.0;
            if(y<0) return -mathlib::PI/2.0;
        }else{
            if(y==0){
                if(x<0) return mathlib::PI;
                if(x>0) return 0.0;
            }
            return atan2(y,x);
        }
        // no reach
        return 0.0;
    }
    /// @brief 点を原点中心にtheta回転させる
    /// @param theta 回転角、ラジアン
    void rotate(double theta){
        double X=x,Y=y;
        x=cos(theta)*X-sin(theta)*Y;
        y=sin(theta)*X+cos(theta)*Y;
    }
};

/// @brief 直線
/// @details Lineは2つの点を持つ
struct Line{
    CoordinatePonit A,B;
    Line(CoordinatePonit p1=CoordinatePonit(0,0),CoordinatePonit p2=CoordinatePonit(1,1)):A(p1),B(p2){}
    double len()const{return (B-A).norm();}
    CoordinatePonit shortest_point(const CoordinatePonit &x)const{
        double a=B.x-A.x,b=B.y-A.y;
        double t=-(a*(A.x-x.x)+b*(A.y-x.y))/(a*a+b*b);
        return CoordinatePonit(a*t+A.x,b*t+A.y);
    }
    /// @brief 点と直線の距離
    /// @param x 距離を求める点
    /// @return 点xと直線ABの距離
    double dist(const CoordinatePonit &x)const{
        return (x-shortest_point(x)).norm();
    }
    /// @brief 点が直線上にあるかを判定する
    /// @param x 判定対象の点
    bool online(const CoordinatePonit &x)const{
        return fabs((B-A)^(x-A))<mathlib::eps&&(x-A).norm()<len()+mathlib::eps&&(x-B).norm()<len()+mathlib::eps;
    }
    /// @brief 線分と直線の交差判定
    /// @param x 直線の一端点
    /// @param y 直線のもう一端点
    /// @return 交差しているならtrue
    bool intersect(const CoordinatePonit &x,const CoordinatePonit &y){
        return ((x.x-y.x)*(A.y-x.y) + (x.y-y.y)*(x.x-A.x))*
                ((x.x-y.x)*(B.y-x.y) + (x.y-y.y)*(x.x-B.x)) < 0;
    }

    /// @brief 線分と線分の交差判定
    /// @param L 判定対象の線分
    /// @return 交差しているならtrue
    bool Lintersect(const Line &L)const{
        return  ((L.A.x-L.B.x)*(A.y-L.A.y) + (L.A.y-L.B.y)*(L.A.x-A.x))*
                ((L.A.x-L.B.x)*(B.y-L.A.y) + (L.A.y-L.B.y)*(L.A.x-B.x))<0&&
                ((A.x-B.x)*(L.A.y-A.y) + (A.y-B.y)*(A.x-L.A.x))*
                ((A.x-B.x)*(L.B.y-A.y) + (A.y-B.y)*(A.x-L.B.x))<0;
    }
    friend std::istream& operator>>(std::istream &is, Line &x){CoordinatePonit X,Y; is>>X>>Y;x.A=X;x.B=Y;return is;}
    friend std::ostream& operator<<(std::ostream &os, const Line &x){ os << x.A <<" "<< x.B; return os; }
};

/// @brief 三角形
/// @details Triangleは3つの頂点を持つ
/// @note A(),B(),C()はそれぞれの頂点の角度
struct Triangle{
    CoordinatePonit a,b,c;
    Triangle(CoordinatePonit p1=CoordinatePonit(0,0),CoordinatePonit p2=CoordinatePonit(1,0), CoordinatePonit p3=CoordinatePonit(0,1)):a(p1),b(p2),c(p3){}

    /// @brief 各頂点の角度
    /// @return それぞれの頂点の角度、ラジアン
    double A()const{double res=fabs((b-a).rad()-(c-a).rad());if(res>=mathlib::PI) res-=mathlib::PI;return res;}
    double B()const{double res=fabs((a-b).rad()-(c-b).rad());if(res>=mathlib::PI) res-=mathlib::PI;return res;}
    double C()const{double res=fabs((a-c).rad()-(b-c).rad());if(res>=mathlib::PI) res-=mathlib::PI;return res;}
    /// @brief 各辺の直線
    /// @return 辺AB,BC,CAに対応する直線
    Line AB()const{return Line(a,b);}
    Line BC()const{return Line(b,c);}
    Line CA()const{return Line(c,a);}

    /// @brief 重心、外心、内心、高心、傍心、外接円の半径、内接円の半径、各辺と対応する高さ
    CoordinatePonit G()const{return (a+b+c)/3;}
    /// @brief 外心
    CoordinatePonit O()const{return (sin(2*A())*a+sin(2*B())*b+sin(2*C())*c)/(sin(2*A())+sin(2*B())+sin(2*C()));}
    /// @brief 外接円の半径
    double R()const{return a.norm()/(2*sin(A()));}
    /// @brief 内心
    CoordinatePonit I()const{return (a.norm()*a+b.norm()*b+c.norm()*c)/(a.norm()+b.norm()+c.norm());}
    /// @brief 内接円の半径
    double r()const{return 2*area()/((a-b).norm()+(b-c).norm()+(c-a).norm());}
    /// @brief 傍心
    /// @details 各頂点から対辺に下ろした垂線の交点
    CoordinatePonit H()const{return (tan(A())*a+tan(B())*b+tan(C())*c)/(tan(A())+tan(B())+tan(C()));}
    /// @brief 各辺と対応する高さ
    /// @return 辺aに対応する高さAH、辺bに対応する高さBH、辺cに対応する高さCH
    /// @note AHは頂点Aから辺BCに下ろした垂線の長さ
    double AH()const{return 2*R()*cos(A());}
    double BH()const{return 2*R()*cos(B());}
    double CH()const{return 2*R()*cos(C());}
    /// @brief 三角形の面積
    double area()const{return ((b-a)^(c-a))/2;}

    /// @brief 点が三角形の内側、辺上、外側のどこにあるかを判定する
    /// @param x 判定対象の点
    /// @return inside()が内側、online()が辺上、outside()が外側を判定する
    bool inside(const CoordinatePonit &x)const{
        return ((b-a)^(x-a))*((c-a)^(x-a))<0&&
               ((a-b)^(x-b))*((c-b)^(x-b))<0&&
               ((a-c)^(x-c))*((b-c)^(x-c))<0;
    }
    /// @brief 点が三角形の辺上にあるかを判定する
    /// @param x 判定対象の点
    bool online(const CoordinatePonit &x)const{return AB().online(x)||BC().online(x)||CA().online(x);}
    bool outside(const CoordinatePonit &x)const{return (!inside(x)&&!online(x));}

    friend std::istream& operator>>(std::istream &is, Triangle &x){CoordinatePonit X,Y,Z; is>>X>>Y>>Z; x.a=X;x.b=Y;x.c=Z;return is;}
    friend std::ostream& operator<<(std::ostream &os, const Triangle &v){ os << v.a << v.b; return os; }
};

/// @brief 円
/// @details Circleは中心と半径を持つ
/// @note inside()は点が円の内側にあるかどうかを判定する
struct Circle{
    CoordinatePonit center;
    double radius;
    Circle(double r=0,CoordinatePonit c={0,0}):radius(r),center(c){}
    bool inside(const CoordinatePonit &C)const{
        return (C-center).norm2()<radius*radius + mathlib::eps;
    }
    bool intersect(const Circle &C)const{
        return (C.center-center).norm2()<(C.radius+radius)*(C.radius+radius)+mathlib::eps;
    }
    std::pair<CoordinatePonit,CoordinatePonit> interP(const Circle &C)const{
        CoordinatePonit Q=center-C.center;
        double x1=Q.x,y1=Q.y,r1=radius,r2=C.radius;
        double a=(x1*x1+y1*y1+r1*r1-r2*r2)/2.0,t=Q.norm2();
        CoordinatePonit A((a*x1+y1*sqrt(t*r1*r1-a*a))/t,(a*y1-x1*sqrt(t*r1*r1-a*a))/t);
        CoordinatePonit B((a*x1-y1*sqrt(t*r1*r1-a*a))/t,(a*y1+x1*sqrt(t*r1*r1-a*a))/t);
        A+=C.center;
        B+=C.center;
        return std::make_pair(A,B);
    }
    friend std::ostream& operator<<(std::ostream &os, const Circle &v){ os << v.center<<" ";printf("%.06f",v.radius); return os; }
};