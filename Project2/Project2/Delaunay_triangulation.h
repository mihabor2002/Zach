#pragma once
#include "Triangle.h"
#include <vector>
#include <fstream>
using namespace std;
//211-Borovikov-Mikhail-2021


class Delaunay_triangulation
{
private:
	vector <Triangle> triangles;
	int number_triangles, number_points;
	vector < vector < int > > triangle_indicators;
public:
	Delaunay_triangulation();
	Delaunay_triangulation(const Delaunay_triangulation& Del_tr);
	const Delaunay_triangulation& operator=(const Delaunay_triangulation& Del_tr);
	~Delaunay_triangulation() = default;
	void add_triangle(Triangle T);
	Triangle get_triangle(int id);
	int get_number_triangles();
	int find_triangle(Point p, Point q);
	void create_triangle_indicators();
	vector <Point> find_neighbouring_points(Point p);
	void print_Delaunay_triangulation();
	void print(int n);
	void print_gif(int n);
};