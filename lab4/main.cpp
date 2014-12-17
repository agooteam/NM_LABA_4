#include "function.h"

SLAU A;
param p;
char *file = "param.txt";
int method = 2;
void main(){
	setlocale(LC_CTYPE,"rus");
	if(A.read_param(file,p)){
		A.allocation_memory(p);
		A.read_x0("x0.txt",p);
		A.Newton(p,method);
		A.clear_memory();
	}
	system("pause");
};

