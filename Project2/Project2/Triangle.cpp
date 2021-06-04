#include "Triangle.h"
//211-Borovikov-Mikhail-2021

Triangle::Triangle() { A = B = C = Point(0, 0); id = 0; }

Triangle::Triangle(Point a, Point b, Point c) { A = a; B = b; C = c; id = 0; }

Triangle::Triangle(const Triangle& T) { A = T.A; B = T.B; C = T.C; id = T.id; neighbouring_triangle_id = T.neighbouring_triangle_id; }

const Triangle& Triangle::operator=(const Triangle& T)
{
	A = T.A; B = T.B; C = T.C; id = T.id; neighbouring_triangle_id = T.neighbouring_triangle_id;
	return *this;
}

Point Triangle::get_A() { return A; }

Point Triangle::get_B() { return B; }

Point Triangle::get_C() { return C; }

int Triangle::get_id() { return id; }

int Triangle::get_neighbouring_triangle_id() { return neighbouring_triangle_id; }

void Triangle::assign_A(Point a) { A = a; }

void Triangle::assign_B(Point b) { B = b; }

void Triangle::assign_C(Point c) { C = c; }

void Triangle::assign_id(int k) { id = k; }

void Triangle::assing_neighbouring_triangle_id(int k) { neighbouring_triangle_id = k; }
void Triangle::print_triangle(std::ofstream& file)
{
    file << A.x << " " << A.y << std::endl;
	file << B.x << " " << B.y << std::endl;
	file << C.x << " " << C.y << std::endl;
}