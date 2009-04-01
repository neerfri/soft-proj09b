/*
//============================================================================
// Name        : soft-proj09b.cpp
// Author      : Neer Friedman & Lea Stolowicz
// Version     :
// Copyright   : Your copyright notice
// Description : Checks password according to policy given in the exercise
//============================================================================
*/
#include <stdio.h>
#include <stdlib.h>

int main(void) {
	
	int ch;
	int is_length = 0; /* will be true if length od password is between 8 and 14 */
	int is_upper = 0; /* will be true if the password starts with an upper letter */
	int is_digit = 0;  /* true if contains at least one digit */
	int is_lower_case = 0; /* true if contains at least one lower case letter */
	int is_character = 0; /* true if contains at least one char which is not a digit or a letter */
	int length = 0;
	
	ch = getchar();
	if (('A' <= ch) && (ch <= 'Z'))
	{
		is_upper = 1;
	} 
	else 
	{
		printf("BAD PASSWORD");
		return 0;
	}
		
	do
	{
		ch = getchar();
		if (ch != '\n')
		{
			if (('0' <= ch) && (ch <= '9'))
				is_digit = 1;
			if (('a' <= ch) && (ch <= 'z'))
				is_lower_case = 1;
			if ((!(('0' <= ch) && (ch <= '9'))) 
				&& (!(('a' <= ch ) && (ch  <= 'z')))
				&&(!(('A' <= ch) && (ch <= 'Z'))))
				is_character = 1;
		}
		++length;
	}
	while (ch != '\n');
	if ((length >= 8) && (length <=14))
		is_length= 1;
	
	if (is_length && is_upper && is_digit && is_lower_case && is_character)
	{
		printf("GOOD PASSWORD");
	}
	else
	{
		printf("BAD PASSWORD");
	}	
	return 0;
}
