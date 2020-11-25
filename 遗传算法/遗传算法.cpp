// �Ŵ��㷨.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <ctime>
#include <math.h>
using namespace std;

const double pi = 3.14159265;
const double pcross = 0.7;//�������
const double pmutate = 0.001;//�������
const int len = 22;//22λ��Ⱦɫ��
const int daishu = 500;//��������
const int Size = 500;//��Ⱥ��ģ
double bestval;//��Ӧֵ���ֵ
typedef struct node {//Ⱦɫ��ṹ��
	bool chromo[len];
}node;
node bestchromo;//��¼���Ÿ���
node group[Size];//��¼��Ⱥ�еĸ��������
node temp[Size];//��¼��Ⱥ�еĸ������ʱ����
void gouzao(node& c) {//�Ե���Ⱦɫ�������ֵ
	for (int i = 0; i < len; i++) {
		c.chromo[i] = rand() % 2;
	}
}
void decode(node& c, double& x) {//�����ƽ������
	double num = 4194394;//��2��22�η�
	double tem = 0;
	for (int i = 0; i < len;i++) {
		tem += c.chromo[i] * pow(2, i);		
	}
	x = (3 / num * tem)-1;
}
double f(double x) {//Ŀ�꺯��
	return x * sin(10 * pi * x) + 2.0;
}
double fitness(node& c) {//��Ӧ�Ⱥ���
	double x;
	decode(c, x);
	return f(x);
}
void cross(node& c1, node& c2, int point) {//�������
	node c3 = c1;
	for (int i = 0; i < len - point; i++) {
		c1.chromo[point + i ] = c2.chromo[point + i ];
	}
	for (int j = 0; j < len - point; j++) {
		c2.chromo[point + j ] = c3.chromo[point + j ];
	}
}
void mutate(node& c) {//�������
	int i = rand() % len;
	c.chromo[i] = !c.chromo[i];
}
double inline rand0() {//����0��1�����С��
	return rand() % 10000 / 10000.0;
}
void select(node group[Size]) {//ѡ�����
	double fitnessval[Size];
	double sum = 0;
	double avgfitness[Size];
	int id[Size];
	for (int i = 0; i < Size; i++) {
		fitnessval[i] = fitness(group[i]);
	}
	for (i = 0; i < Size; i++) {//��Ӧ���ܺ�
		sum += fitnessval[i];
	}
	for (i = 0; i < Size; i++) {
		avgfitness[i] = fitnessval[i] / sum;
	}
	for (i = 1; i < Size; i++) {//��Ӧ���ۼ�
		avgfitness[i] += avgfitness[i - 1];
	}
	for (i = 0; i < Size; i++) {//���̶�ѡ��
		double rannum = rand0();//����0��1�����
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
	for (i = 0; i < Size; i++) {//���¸����滻�ɸ���
		temp[i] = group[i];
	}
	for (i = 0; i < Size; i++) {
		group[i] = temp[id[i]];
	}
}
int getBest(node group[Size], double& x, double& number) {//ȡ�����Ÿ����Ӧ��λ��
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
void GA(double& x, double& number) {//�Ŵ��㷨����
	for (int i = 0; i < Size; i++) {
		gouzao(group[i]);
	}
	bestchromo = group[getBest(group, x, bestval)];
	for (i = 0; i < daishu; i++) {
		select(group);//ѡ�����
		int p = rand() % len;
		for (int j = 0, pre = -1; j < Size; j++) {//���ݸ��ʽ���		
			if (rand0() < pcross) {
				if (pre == -1)
					pre = j;
				else {
					cross(group[pre], group[j], p);
					pre = -1;
				}
			}
		}
		for (int k = 0; k < Size; k++) {//���ݸ��ʽ��б���
			if ((rand0() < pmutate)) {
				mutate(group[k]);
			}
		}
		getBest(group, x, number);
		cout << "��" << i << "��" << "����xֵΪ:" << x << "����ֵΪ" << f(x) << endl;//��������
	}
}
int main() {
	srand((unsigned)time(0));//�������������
	double x;
	double max;
	GA(x, max);
	return 0;
}



