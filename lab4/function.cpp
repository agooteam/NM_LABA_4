#include "function.h"
TYPE h = 0.00001;
//массив функций
/*TYPE f1(Vector x){return  ((x.mas[0]-1.0)*(x.mas[0]-1.0)+(x.mas[1]-2.0)*(x.mas[1]-2.0)-1.0);}
TYPE f2(Vector x){return  ((x.mas[0]-4.0)*(x.mas[0]-4.0)+(x.mas[1]-2.0)*(x.mas[1]-2.0)-1.0);}

TYPE func1(int i, Vector x){
	if(i==0) return (2.0*x.mas[0]-2.0);
	if(i==1) return (2.0*x.mas[1]-4.0);
}
TYPE func2(int i, Vector x){
	if(i==0) return (2.0*x.mas[0]-8.0);
	if(i==1) return (2.0*x.mas[1]-4.0);
}

TYPE (*func[])(int, Vector)={func1, func2};//производные
TYPE (*f[])(Vector)={f1, f2};//функции
*/
TYPE c = 100000;
TYPE f1(Vector x,TYPE h,int i){
	if(i==0) return  c*((x.mas[0]+h-2.0)*(x.mas[0]+h-2.0)+(x.mas[1]-2.0)*(x.mas[1]-2.0)-1.0 );
	if(i==1) return  c*((x.mas[0]-2.0)*(x.mas[0]-2.0)+(x.mas[1]+h-2.0)*(x.mas[1]+h-2.0)-1.0 );
	if(i==2) return  c*((x.mas[0]-2.0)*(x.mas[0]-2.0)+(x.mas[1]-2.0)*(x.mas[1]-2.0)-1.0 );
}
TYPE f2(Vector x,TYPE h,int i){
	if(i==0)return  c*((x.mas[0]+h-4.0)*(x.mas[0]+h-4.0)+(x.mas[1]-2.0)*(x.mas[1]-2.0)-1.0 );
	if(i==1)return  c*((x.mas[0]-4.0)*(x.mas[0]-4.0)+(x.mas[1]+h-2.0)*(x.mas[1]+h-2.0)-1.0 );
	if(i==2)return  c*((x.mas[0]-4.0)*(x.mas[0]-4.0)+(x.mas[1]-2.0)*(x.mas[1]-2.0)-1.0 );
}
TYPE f3(Vector x,TYPE h,int i){
	if(i==0)return  c*(x.mas[1] - 2*(x.mas[0] + h) + 4 );
	if(i==1)return  c*(x.mas[1] + h - 2*x.mas[0] + 4);
	if(i==2)return  c*(x.mas[1] - 2*x.mas[0] + 4);
}

TYPE func1(int i, Vector x){
	if(i==0) return c*(2.0*x.mas[0]-4.0);
	if(i==1) return c*(2.0*x.mas[1]-4.0);
	if(i==2) return 0;
}
TYPE func2(int i, Vector x){
	if(i==0) return c*(2.0*x.mas[0]-8.0);
	if(i==1) return c*(2.0*x.mas[1]-4.0);
	if(i==2) return 0;
}
TYPE func3(int i, Vector x){
	if(i==0) return c*-2;
	if(i==1) return c*1;
	if(i==2) return 0;
}

/*TYPE f1(Vector x,TYPE h,int i){
	if(i==0) return  (x.mas[1] - sin(x.mas[0] + h));
	if(i==1) return  (x.mas[1] + h - sin(x.mas[0]));
	if(i==2) return  ((x.mas[1] - sin(x.mas[0])));
}
TYPE f2(Vector x,TYPE h,int i){
	if(i==0)return  (x.mas[1] - (x.mas[0]+h)/2);
	if(i==1)return  (x.mas[1] + h - x.mas[0]/2);
	if(i==2)return  (x.mas[1] - x.mas[0]/2);
}

TYPE func1(int i, Vector x){
	if(i==0) return cos(x.mas[0]);
	if(i==1) return 1;
	if(i==2) return 0;
}
TYPE func2(int i, Vector x){
	if(i==0) return -1;
	if(i==1) return 1;
	if(i==2) return 0;
}
*/
TYPE (*func[])(int, Vector)={func1, func2,func3};//производные аналитика
TYPE (*f[])(Vector,TYPE, int)={f1, f2,f3};//функции


TYPE Norma_f(Vector x, int M){
	TYPE res = 0.0;
	for(int i = 0; i < M; i++)
		res += f[i](x,0,i)*f[i](x,0,i);
	res = sqrt(res);
	return res;
}

bool SLAU::read_param(const char *file, param &T){
	bool result = true;
	FILE *fp = fopen(file,"r");
	if(fp){
		fscanf(fp,"%d",&T.N);
		fscanf(fp,"%d",&T.M);
		fscanf(fp,"%d",&T.maxiter);
		fscanf(fp,"%lf",&T.eps1);
		fscanf(fp,"%lf",&T.eps2);
	}
	else result = false;
	return result;
};

void SLAU::allocation_memory(param T){
	Jacoby = new TYPE *[T.M];
	for(int i = 0; i < T.M; i++) Jacoby[i] = new TYPE[T.N];
	x_0.mas = new TYPE[T.N];
	x.mas = new TYPE[T.N];
	delta_x.mas = new TYPE[T.N];
	F.mas = new TYPE[T.M];
};

bool SLAU::read_x0(const char *file, param T){
	bool result = true;
	FILE *fp = fopen(file,"r");
	if(fp) for(int i = 0 ; i < T.N ; i++) fscanf(fp,"%lf",&x_0.mas[i]);
	else result = false;
	return result;
};
void SLAU::calc_Jacoby(param T,int method){
	for(int i = 0; i < T.M; i++) {
		if(method == 1)for(int j = 0; j < T.N; j++) Jacoby[i][j] = func[i](j, x);
		else if(method == 2)for(int j = 0; j < T.N; j++) Jacoby[i][j] = (f[i](x,h,j)-f[i](x,0,j))/h;
		
		
	}
};

void SLAU::calc_F(param T){
	for(int i = 0 ; i < T.M ;i++) F.mas[i] = -f[i](x,0,i);
};

void SLAU::sort(param T){
	int str1, str2;
	for(str1 = 0; str1 < T.M; str1++) {
		str2 = FindMax_F(F, T.M, str1);
		ChangeStr(T.M, str1, str2);
	}
};

void SLAU::svertka(param T){
	for(int i = T.N; i < T.M; i++) { //возвели строки и пр.часть в квадрат
		for(int j = 0; j < T.N; j++)
			Jacoby[i][j] = Jacoby[i][j]*2.0*f[i](x,0,j);
		F.mas[i] = -f[i](x,0,i)*f[i](x,0,i);
	}
	for(int i = T.N; i < T.M; i++) {//прибавили к N-ному
		for(int j = 0; j < T.N; j++)
			Jacoby[T.N-1][j] += Jacoby[i][j];
			F.mas[T.N-1] += F.mas[i];
	}

};

void SLAU::ChangeStr(int N, int str1, int str2){
	int j;
	TYPE ch;
	for(j = 0; j < N; j++) {
		ch = Jacoby[str1][j];
		Jacoby[str1][j] = Jacoby[str2][j];
		Jacoby[str2][j] = ch;
	}
	ch = F.mas[str1];
	F.mas[str1] = F.mas[str2];
	F.mas[str2] = ch;
};

int SLAU::FindMax_F(Vector x, int N, int str0){
	int numStr = str0;
	TYPE maxEl = x.mas[str0];
	for(int i = str0; i < N; i++){
		if(fabs(x.mas[i]) > fabs(maxEl)){
			maxEl = x.mas[i];
			numStr = i;
		}
	}
	return numStr;
};
int SLAU::FindMax_Stolb(param T,int stolb, int str0)	{
	int numStr = str0;
	TYPE maxEl = Jacoby[str0][stolb];
	for(int i = str0; i < T.N; i++) {
		if(fabs(Jacoby[i][stolb]) > fabs(maxEl)) {
			maxEl = Jacoby[i][stolb];
			numStr = i;
		}
	}
	if(maxEl = 0) 
		numStr = -1;
	return numStr;
}
void SLAU::Gaus(param T,FILE *fp){
	int i, j, k;
	TYPE t, sum;

	for (k = 0; k < T.N; k++) { //номер этапа
		i = FindMax_Stolb(T, k, k);
		if(i == -1) {
			fprintf(fp,"\nError. det(A) = 0\n");
			system("pause");
			exit(0);
		}
		if(i != k) ChangeStr(T.N, i, k);
		for(i = k+1; i < T.N; i++) { //номер строки
			t = Jacoby[i][k]/Jacoby[k][k];
			F.mas[i] -= t*F.mas[k];
			for(j = k+1; j < T.N; j++) {
				Jacoby[i][j] -= t*Jacoby[k][j];
			}
		}
	}
	delta_x.mas[T.N-1] = F.mas[T.N-1]/Jacoby[T.N-1][T.N-1];
	for(k = T.N-2; k >= 0; k--){
		sum = 0;
		for(j = k+1; j < T.N; j++) sum += Jacoby[k][j]*delta_x.mas[j];
		delta_x.mas[k]=(F.mas[k]-sum)/Jacoby[k][k];
	}
};

void SLAU::calc_x(param T){
	for(int i = 0; i < T.N; i++)
		x.mas[i] = x_0.mas[i]+betta*delta_x.mas[i];
}

TYPE SLAU::calc_betta(param T){
	TYPE norm1, norm2;
	while(betta > T.eps2) {
		calc_x(T);
		norm1 = Norma_f(x,T.M);
		norm2 = Norma_f(x_0,T.M);
		if (norm1 > norm2) 
			betta = betta/2.0;
		else return betta;
	}
	return betta;
}

void SLAU::Newton(param T,int method){
	nevyazka = 1;
	int i,flag = 0;
	betta = 1;
	FILE *fp = fopen("MRES.txt","w");
	if(fp){
		Norm = Norma_f(x_0, T.M);
		for(i = 0 ; i < T.N ; i++) x.mas[i] = x_0.mas[i];
		for (i = 0; i < T.maxiter && flag == 0; i++) {
			fprintf(fp, "%d\t%.16lf\t%.16lf\t%le\t%lf\n",i, x_0.mas[0], x_0.mas[1],nevyazka,betta);
			calc_Jacoby(T,method);//+
			calc_F(T);//+
			if(T.M > T.N){
				if(method == 1) sort(T);
				else if(method == 2) svertka(T);
			}
			Gaus(T,fp);//delta_x
			betta = 1;
			betta = calc_betta(T);
			nevyazka = Norma_f(x, T.M)/Norm;
			for(int h = 0 ; h < T.N; h++) x_0.mas[h] = x.mas[h];//косяк"
			if(fabs(nevyazka) < T.eps1) flag = 1;
			if(fabs(betta) < T.eps2) flag = 2;
		}
		fprintf(fp, "%d\t%.16lf\t%.16lf\t%le\t%lf\n",i, x.mas[0], x.mas[1],nevyazka,betta);
		if(fabs(nevyazka) > T.eps2) fprintf(fp, "решение не найдено\n");	
		if(flag == 1) fprintf(fp, "вышли по невязке\n");
		if(flag == 2) fprintf(fp, "вышли по бетта\n");
		if(flag == 0) fprintf(fp, "аварийный выход\n");
		fclose(fp);
	}
};

void SLAU::clear_memory(){
	delete []Jacoby;
	delete []x_0.mas;
	delete []x.mas;
	delete []delta_x.mas;
	delete []F.mas;
};