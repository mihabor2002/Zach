#pragma once
#include "point.h"
#include <fstream>
//211-Borovikov-Mikhail-2021

class Triangle
{
private:
	Point A, B, C;
	int id;
	int neighbouring_triangle_id;
public:
	Triangle();
	Triangle(Point a, Point b, Point c);
	Triangle(const Triangle& T);
	~Triangle() = default;
	const Triangle& operator=(const Triangle& T);
	Point get_A();
	Point get_B();
	Point get_C();
	int get_id();
	int get_neighbouring_triangle_id();
	void assign_A(Point a);
	void assign_B(Point b);
	void assign_C(Point c);
	void assign_id(int k);
	void assing_neighbouring_triangle_id(int k);
	void print_triangle(std::ofstream& file);
	float Square(Point aA, Point aB, Point aC) { return abs(aB.x * aC.y - aC.x * aB.y - aA.x * aC.y + aC.x * aA.y + aA.x * aB.y - aB.x * aA.y); }
	int Belongs_P(Point P) 
	{ 
		return abs(Square(P,A,B)+Square(P,B,C)+Square(P,A,C)-Square(A,B,C))<0.01; 
	};
};