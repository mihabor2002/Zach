#include <vector>
#include <cmath>
#include "point.h"
//211-Borovikov-Mikhail-2021
#define epsi 1e-5
bool eq(float a, float b)
{
    float eps = 1e-6;
    return fabs(a - b) < eps;
}
float dist(Point A, Point B)//finkciya prinimaet na vhod 2 tochki A i B i schitaet rasstoyanie mezhdu nimi
{
    return std::sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
}
bool P_eq(Point A, Point B)
{
    if (std::fabs(A.x - B.x) < epsi)
    {
        if (std::fabs(A.y - B.y) < epsi)
        {
            return 1;
        }
        return 0;
    }
    return 0;
}
bool po_chas(Point A, Point B, Point C)/*по часовой ли стрелке лежит точка С по отношения к вектору АВ*/
{
    return B.x * A.y - A.x * B.y + C.x * B.y - C.y * B.x + A.x * C.y - C.x * A.y >= 0;
}
Line points_to_line(Point P, Point Q)/*прямая по точкам P, Q*/
{
    Line F;
    F.A = P.y - Q.y;
    F.B = Q.x - P.x;
    F.C = -F.A * P.x - F.B * P.y;
    return F;
}
bool prinadlezhit_pryamoy(Point A, Line l)/*принадлежит ли точка А прямой l*/
{
    return eq(l.A * A.x + l.B * A.y + l.C, 0);
}
bool prinadlezhit_otrezku(Point A, Point P, Point Q)/*принадлежит ли точка А отрезку PQ*/
{
    Line PQ = points_to_line(P, Q);
    double x_min = fmin(P.x, Q.x), x_max = fmax(P.x, Q.x), y_min = fmin(P.y, Q.y), y_max = fmax(P.y, Q.y), ax = A.x, ay = A.y;
    return prinadlezhit_pryamoy(A, PQ) && x_min <= ax && ax <= x_max && y_min <= ay && ay <= y_max;
}
int kolvo_peresech(Line M, Line N)
{
    if (eq(M.A * N.B, N.A * M.B))
    {
        if (eq(N.A * M.C, N.C * M.A) && eq(N.B * M.C, N.C * M.B))
        {
            return 2;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return 1;
    }
}
Point peresech(Line M, Line N)
{
    Point P;
    if (kolvo_peresech(M, N) == 1)
    {
        P.x = (N.C * M.B - N.B * M.C) / (M.A * N.B - M.B * N.A);
        P.y = -(N.C * M.A - N.A * M.C) / (M.A * N.B - M.B * N.A);
        return P;
    }
    else
    {
        return P;
    }
}
Point podelit_otrezok_v_otnoshenii(Point X, Point Y, double a, double b)
{
    Point P;
    P.x = (b * X.x + a * Y.x) / (a + b);
    P.y = (b * X.y + a * Y.y) / (a + b);
    return P;
}
Point middle_otrezok(Point X, Point Y)
{
    return podelit_otrezok_v_otnoshenii(X, Y, 1, 1);
}
Line perpendikul(Point P, Line l)
{
    Line M;
    M.A = -l.B;
    M.B = l.A;
    M.C = -M.A * P.x - M.B * P.y;
    return M;
}
Line serper(Point A, Point B)
{
    return perpendikul(middle_otrezok(A, B), points_to_line(A, B));
}
bool po_odnu_storonu(Point P, Point Q, Line l)
{
    Line PQ = points_to_line(Q, P);
    Point A = peresech(PQ, l);
    if (kolvo_peresech(PQ, l) == 1 && fmin(Q.x, P.x) <= A.x && A.x <= fmax(P.x, Q.x) && fmin(Q.y, P.y) <= A.y && A.y <= fmax(Q.y, P.y))
    {
        return false;
    }
    else
    {
        return true;
    }
}

Point::Point(double x, double y)
{
    x = x;
    y = y;
    ID = 0;
}

Point::Point(const Point& v)
{
    x = v.x;
    y = v.y;
    ID = v.ID;
}
std::vector <Point> obol(std::vector <Point> A, int n)/*функция вычисляет оболочку множества n множества точек A, и записывает в нужном порядке номера точек,образующих выпуклую оболочку, из массива р в массив nom_obol, функция возвращает количество точек выпуклой оболочки*/
{
    int first, q, next, i, j;/*first - переменная для поиска первой точки, q - переменная для номера текущей точки,next - переменная для поиска следующей точки, i - переменная для прохода по циклу,j-количество точек выпуклой оболочки */
    j = 0;
    first = 0;
    for (i = 1; i < n; ++i)/*находим самую нижнюю из самых левых точек*/
    {
        if (A[i].x < A[first].x || (A[i].x == A[first].x && A[i].y < A[first].y))
        {
            first = i;
        }
    }
    std::vector <int> nom_obol;
    q = first; /* текущая точка*/
    do
    {
        nom_obol.push_back(q);/*добавляем текущую точку в */
        next = q;/*ищем следующую точку*/
        for (i = n - 1; i >= 0; --i)
        {
            if (A[i].x != A[q].x || A[i].y != A[q].y)/*если это не текущая точка*/
            {
                if (next == q || po_chas(A[q], A[i], A[next]) || (prinadlezhit_otrezku(A[next], A[q], A[i])))
                {
                    next = i;/*вычисляем следующую точку, которая должна быть самой левой,берём первую, если ещё не брали в этом цикле, если не первая, проверяем, левее(относительно текущей точки) ли она самой левой из прошлых, или она лежит на отрезке из прошлой самой левой и текущей*/
                }
            }
        }
        q = next;
    }     while (q != first);
    std::vector <Point> otv;
    for (i = 0; i < nom_obol.size(); i++)
    {
        otv.push_back(A[nom_obol[i]]);
    }
    return otv;
}
const Point& Point::operator=(const Point& v)
{
    x = v.x;
    y = v.y;
    ID = v.ID;
    return *this;
}

Point operator+(const Point v, const Point u)
{
    Point P;
    P.x = v.x + u.x;
    P.y = v.y + u.y;
    return P;
}

Point operator-(const Point v, const Point u)
{
    Point P;
    P.x = v.x-u.x;
    P.y = v.y - u.y;
    return P;
}

Point operator*(double lambda, const Point& v)
{
    return Point(lambda * v.x, lambda * v.y);
}

Point operator*(const Point& v, double lambda)
{
    return Point(lambda * v.x, lambda * v.y);
}
double Point::length() { return sqrt(x * x + y * y); }
double Point::get_function_value()
{
    return sin(x * x + y * y);
}