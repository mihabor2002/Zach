#pragma once
#include <vector>
#include <fstream>
//211-Borovikov-Mikhail-2021
bool eq(float a, float b);
class Point
{
public:
    float x;//x koordinata tochki
    float y;//y koordinata tochki
    int a;//moment vremeni, v kotorii eta tochka "gorit" v volnovom algoritme 
    int b;//moment vremeni +1, kogda prediduschaya tochka gorit 
    char DB;//class tochki v algorotme DBSCAN, 0 - vneschnya, 1 - vnutrennaya, 2 - pereferiynaya
    int ID;//poryadkovii nomer tochki v polye
    int for_k_sr;//peremennaya dlya algoritma k-srednikh
    float P_c;
    int Cl_forel;
    Point()
    {
        x = 0;
        y = 0;
        for_k_sr = -1;
        DB = 0;
        a = 0;
        b = 0;
        ID = 0;
        P_c = 0;
    }
    Point(double x, double y);
    Point(const Point& v);
    std::vector <float> c_po_x;
    std::vector <float> x_po_c;
    void printToFile(std::ofstream& out, int k)
    {
        if (k == 0)
        {
            out << x << " " << y << std::endl;
        }
    }
    const Point& operator=(const Point& v);
    friend Point operator+(const Point v, const Point u);
    friend Point operator-(const Point v, const Point u);
    friend Point operator*(double lambda, const Point& v);
    friend Point operator*(const Point& v, double lambda);
    double length();
    double get_first_coordinate();
    double get_second_coordinate();
    double get_function_value();
    int get_Point_id();
    void assign_Point_id(int id);
    void print_Point(std::ofstream& file);
};
class Line
{
public:
    float A;
    float B;
    float C;
};
class Circle
{
public:
    Point O;
    float R;
};
float dist(Point A, Point B);
bool P_eq(Point A, Point B);
bool po_chas(Point A, Point B, Point C);
Line points_to_line(Point P, Point Q);
bool prinadlezhit_pryamoy(Point A, Line l);
bool prinadlezhit_otrezku(Point A, Point P, Point Q);
int kolvo_peresech(Line M, Line N);
Point peresech(Line M, Line N);
Point podelit_otrezok_v_otnoshenii(Point X, Point Y, double a, double b);
Point middle_otrezok(Point X, Point Y);
Line perpendikul(Point P, Line l);
Line serper(Point A, Point B);
bool po_odnu_storonu(Point P, Point Q, Line l);

