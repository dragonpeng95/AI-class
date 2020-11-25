// 遗传算法.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <ctime>
#include <math.h>
using namespace std;

const double pi = 3.14159265;
const double pcross = 0.7;//交叉概率
const double pmutate = 0.001;//变异概率
const int len = 22;//22位的染色体
const int daishu = 500;//更迭代数
const int Size = 500;//种群规模
double bestval;//适应值最大值
typedef struct node {//染色体结构体
	bool chromo[len];
}node;
node bestchromo;//记录最优个体
node group[Size];//记录种群中的个体的数组
node temp[Size];//记录种群中的个体的临时数组
void gouzao(node& c) {//对单个染色体随机赋值
	for (int i = 0; i < len; i++) {
		c.chromo[i] = rand() % 2;
	}
}
void decode(node& c, double& x) {//二进制解码操作
	double num = 4194394;//即2的22次方
	double tem = 0;
	for (int i = 0; i < len;i++) {
		tem += c.chromo[i] * pow(2, i);		
	}
	x = (3 / num * tem)-1;
}
double f(double x) {//目标函数
	return x * sin(10 * pi * x) + 2.0;
}
double fitness(node& c) {//适应度函数
	double x;
	decode(c, x);
	return f(x);
}
void cross(node& c1, node& c2, int point) {//交叉操作
	node c3 = c1;
	for (int i = 0; i < len - point; i++) {
		c1.chromo[point + i ] = c2.chromo[point + i ];
	}
	for (int j = 0; j < len - point; j++) {
		c2.chromo[point + j ] = c3.chromo[point + j ];
	}
}
void mutate(node& c) {//变异操作
	int i = rand() % len;
	c.chromo[i] = !c.chromo[i];
}
double inline rand0() {//产生0到1的随机小数
	return rand() % 10000 / 10000.0;
}
void select(node group[Size]) {//选择操作
	double fitnessval[Size];
	double sum = 0;
	double avgfitness[Size];
	int id[Size];
	for (int i = 0; i < Size; i++) {
		fitnessval[i] = fitness(group[i]);
	}
	for (i = 0; i < Size; i++) {//适应度总和
		sum += fitnessval[i];
	}
	for (i = 0; i < Size; i++) {
		avgfitness[i] = fitnessval[i] / sum;
	}
	for (i = 1; i < Size; i++) {//适应度累加
		avgfitness[i] += avgfitness[i - 1];
	}
	for (i = 0; i < Size; i++) {//轮盘赌选择法
		double rannum = rand0();//产生0到1随机数
		int j;
		for (j = 0; j < Size - 1; j++) {
			if (rannum < avgfitness[j]) {
				id[i] = j;
				break;
			}
		}
		if (j == Size - 1) {
			id[i] = j;
		}
	}
	for (i = 0; i < Size; i++) {//将新个体替换旧个体
		temp[i] = group[i];
	}
	for (i = 0; i < Size; i++) {
		group[i] = temp[id[i]];
	}
}
int getBest(node group[Size], double& x, double& number) {//取得最优个体对应的位置
	double fitnessval[Size];
	for (int i = 0; i < Size; i++) {
		fitnessval[i] = fitness(group[i]);
	}
	int id = 0;
	for (i = 1; i < Size; i++) {
		if (fitnessval[i] > fitnessval[id]) {
			id = i;
		}
	}
	decode(group[id], x);
	number = f(x);
	return id;
}
void GA(double& x, double& number) {//遗传算法流程
	for (int i = 0; i < Size; i++) {
		gouzao(group[i]);
	}
	bestchromo = group[getBest(group, x, bestval)];
	for (i = 0; i < daishu; i++) {
		select(group);//选择操作
		int p = rand() % len;
		for (int j = 0, pre = -1; j < Size; j++) {//根据概率交叉		
			if (rand0() < pcross) {
				if (pre == -1)
					pre = j;
				else {
					cross(group[pre], group[j], p);
					pre = -1;
				}
			}
		}
		for (int k = 0; k < Size; k++) {//根据概率进行变异
			if ((rand0() < pmutate)) {
				mutate(group[k]);
			}
		}
		getBest(group, x, number);
		cout << "第" << i << "代" << "最优x值为:" << x << "函数值为" << f(x) << endl;//结果的输出
	}
}
int main() {
	srand((unsigned)time(0));//产生随机数种子
	double x;
	double max;
	GA(x, max);
	return 0;
}



