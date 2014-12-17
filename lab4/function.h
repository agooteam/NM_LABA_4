#include <stdio.h>
#include <conio.h>
#include "Math.h"
#include "iostream"
using namespace std;

typedef double TYPE;
struct Vector{
	TYPE *mas;
};

struct param{
	int N, M; //число переменных, число уравнений
	TYPE eps1, eps2;//ограничения на невязку, на бетта
	int maxiter;//максимальное число итераций
};

class SLAU {
	protected:
		TYPE nevyazka;
		TYPE Norm;
		TYPE betta;
		TYPE **Jacoby;
		Vector x_0;
		Vector x;
		Vector delta_x;
		Vector F;

	public:
		bool read_param(const char *file,param &M);
		bool read_x0(const char *file, param M);
		TYPE calc_norm(Vector *v, int N);
		void calc_F(param M);
		void allocation_memory(param M);
		void clear_memory();
		void calc_Jacoby(param M,int method);
		void Newton(param M, int k);
		void Gaus(param M,FILE *fp);
		void sort(param M);
		void svertka(param M);
		int FindMax_F(Vector x, int N, int str0);
		void ChangeStr(int N, int str1, int str2);
		int FindMax_Stolb(param M,int stolb,int str0);
		TYPE calc_betta(param M);
		void calc_x(param T);
};

 