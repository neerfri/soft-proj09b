/*
//============================================================================
// Name        : soft-proj09b.cpp
// Author      : Neer Friedman & Lea Stolowicz
// Version     :
// Copyright   : Your copyright notice
// Description : Gives an approximation to pi, according to the exercise 
directions
//============================================================================
*/
#include <stdio.h>
#include <stdlib.h>

struct MyPoint {
	double x;
	double y;
};

double approximate_pi(int number_of_runs);
extern double pow(double, double);

void randomize_point(struct MyPoint *p) {
	(*p).x = (double)rand()/RAND_MAX;
	(*p).y = (double)rand()/RAND_MAX;
	/*printf("[%f, %f]\n",(*p).x ,(*p).y);*/	
}

double approximate_pi(int number_of_runs) {
	struct MyPoint p;
	int i, count = 0;
	double pi;
	for (i=0; i<number_of_runs; i++) {
		randomize_point(&p);
		if ((pow(p.x, 2.0) + pow(p.y, 2.0)) <= 1) {
			count++;
		}
	}
	pi = (double)count/number_of_runs;
	pi = pi * 4;
	return(pi);
}

int main(void) {
	int i, n;
	srand(1234);
	for(i=0; i<8; i++) {
		n = pow(10,i+1);
		/*printf("%d:", n)*/
		printf("%f\n", approximate_pi(n));
	}
	return(0);
}
