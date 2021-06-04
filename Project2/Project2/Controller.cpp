#include <cmath>]
//211-Borovikov-Mikhail-2021
float sum_200(void)//funkcia vidaet summu 200 sluchainikh chisel
{
	float s = 0, l;
	int n;
	for (int i = 0; i < 200; i++)
	{
		n = rand() & 20000;
		l = float(n);
		s += l / 100 - 100;
	}
	return s;
}
