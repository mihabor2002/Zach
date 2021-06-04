#include "field.h"
#include <iostream>
#include <math.h>
//211-Borovikov-Mikhail-2021
int vspom(int k, std::vector <int>& v)
{
	int n;
	for (int i = 0; i < v.size(); i++)
	{
		if (v[i] == k)
		{
			n = i;
			std::cout << "pep";
			break;
		}
	}
	return n;
}
void Field::find_min_y(std::vector <Point>& points)
{
    int k = 0;
    for (int i = 0; i < points.size(); i++)
    {
        if ((i != k && points[i].y < points[k].y) ||
            (i != k && points[i].y == points[k].y &&
                points[i].x < points[k].x))
        {
            k = i;
        }
    }
    rotate(points.begin(), points.begin() + k, points.end());
}

void Field::polar_angle_sort(vector <Point>& points)
{
    vector <double> angles;
    Point v;
    double a;

    for (int i = 0; i < points.size(); i++)
    {
        if (i != 0)
        {
            if ((points[i] - points[0]).length() < 0.0001) angles.push_back(0);
            else angles.push_back(acos((points[i] - points[0]).x / (points[i] - points[0]).length()));
        }
        else angles.push_back(0);
    }
    for (int i = 0; i < points.size() - 1; i++)
    {
        for (int j = 0; j < points.size() - 1 - i; j++)
        {
            if (angles[j] > angles[j + 1])
            {
                a = angles[j];
                v = points[j];
                angles[j] = angles[j + 1];
                points[j] = points[j + 1];
                angles[j + 1] = a;
                points[j + 1] = v;
            }
        }
    }
}

float Field::left_rotate(Point p, Point q, Point r)
{
    return (q.x-p.x) * (r.y - q.y) - (q.y-p.y) * (r.x-q.x);
}

vector <Point> Field::convex_hull(vector <Point> points)
{
    vector <Point> q;

    if (points.size() < 3) return points;
    find_min_y(points);
    q.push_back(points[0]);
    polar_angle_sort(points);
    q.push_back(points[1]);
    for (int i = 2; i < points.size(); i++)
    {
        while (left_rotate(q[q.size() - 2], q[q.size() - 1], points[i]) < 0 && q.size() > 1) q.pop_back();
        q.push_back(points[i]);
    }
    q.push_back(points[0]);
    return q;
}

Delaunay_triangulation Field::generate_delaunay_trinagulation(vector <Point>& points)
{
    Delaunay_triangulation triangulation;
    vector <Point> vertices;
    int i, j, k=0;
    Point a, b, c;
    vector <int> ind(points.size(), 0);

    for (i = 0; i < points.size() - 1; i++)
    {
        for (j = 0; j < points.size() - 1 - i; j++)
        {
            if (points[j].x > points[j + 1].x)
            {
                a = points[j];
                points[j] = points[j + 1];
                points[j + 1] = a;
            }
        }
    }
    a = points[0];
    b = points[1];
    c = points[2];
    triangulation.add_triangle(Triangle(a, b, c));
    ind[0] = ind[1] = ind[2] = 1;
    c = points[2];
    for (i = 3; i < points.size(); i++)
    {
        vertices.clear();
        for (j = 0; j < points.size(); j++)
        {
            if (ind[j] == 1) vertices.push_back(points[j]);
        }
        vertices = convex_hull(vertices);
        for (j = 0; j < vertices.size(); j++)
        {
            if (vertices[j].x == c.x &&
                vertices[j].y == c.y) k = j;
        }
        rotate(vertices.begin(), vertices.begin() + k, vertices.end());
        vertices.push_back(vertices[0]);
        a = c - points[i];
        a.ID=c.ID;
        b = vertices[1] - points[i];
        b.ID=vertices[1].ID;
        for (j = 2; j < vertices.size(); j++)
        {
            if (a.x * b.y - a.y * b.x < 0)
            {
                Point p = a + points[i];
                Point q = b + points[i];
                p.ID=a.ID;
                q.ID=b.ID;
                Triangle T = Triangle(points[i], p, q);
                T.assing_neighbouring_triangle_id(triangulation.find_triangle(p, q));
                triangulation.add_triangle(T);
                ind[i] = 1;
            }
            if (a.x * b.y - a.y * b.x > 0) break;
            a = b;
            b = vertices[j] - points[i];
            b.ID=vertices[j].ID;
        }
        reverse(vertices.begin(), vertices.end());
        a = c - points[i];
        a.ID=c.ID;
        b = vertices[1] - points[i];
        b.ID=vertices[1].ID;
        for (j = 2; j < vertices.size(); j++)
        {
            if (a.x * b.y - a.y * b.x > 0)
            {
                Point p = a + points[i];
                Point q = b + points[i];
                p.ID=a.ID;
                q.ID=b.ID;
                Triangle T = Triangle(points[i], p, q);
                T.assing_neighbouring_triangle_id(triangulation.find_triangle(p, q));
                triangulation.add_triangle(T);
                ind[i] = 1;
            }
            if (a.x * b.y - a.y * b.x < 0) break;
            a = b;
            b = vertices[j] - points[i];
            b.ID=vertices[j].ID;
        }
        c = points[i];
    }
    return triangulation;
}