/*
//============================================================================
// Name        : soft-proj09b.cpp
// Author      : Neer Friedman
// Version     :
// Copyright   : Your copyright notice
// Description : Checks division by 3 of a number
//============================================================================
*/
#include <stdio.h>
#include <stdlib.h>

int digits_sum(int n) {
	int sum = 0;
	while(n > 0) {
		sum = sum + (n % 10);
		n = n/10;
	}
	return(sum);
}

int main(void) {
	int n;
	scanf("%d", &n);
	while(n > 9) {
		printf("%d\n", n);
		n = digits_sum(n);
	}
	printf("%d\n", n);
	if (n == 9 || n==6 || n==3)
		printf("DIVIDABLE");
	else
		printf("NOT DIVIDABLE");
	return EXIT_SUCCESS;
}
