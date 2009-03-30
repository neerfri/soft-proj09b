/*
//============================================================================
// Name        : soft-proj09b.cpp
// Author      : Neer Friedman
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================
*/
#include <stdio.h>
#include <stdlib.h>

int main(void) {
	
	int ch;
	int is_length = 0; /* will be true if length od password is between 8 and 14 */
	int is_upper = 0; /* will be true if the password starts with an upper letter 
							65-90 in ascii code */
	int is_digit = 0;  /* true if contains at least one digit 
					 48-57 in ascii code */
	int is_lower_case = 0; /* true if contains at least one lower case letter
							 97-122 in ascii code */
	int is_character = 0; /* true if contains at least one char which is not a digit or a letter */
	int length = 0;
	
	ch = getchar();
	if (('A' <= ch) && (ch <= 'Z'))
	{
		is_upper = 1;
	} else 
	{
		printf("BAD PASSWORD\n");
		printf("First letter in not upper case\n");
		return 0;
	}
		
	do
	{
		ch = getchar();
		
		if ('0' <= ch) && (ch <= '9'))
			is_digit = 1;
		if (('a' <= ch) && (ch <= 122))
			is_lower_case = 1;
		if (!('0' <= ch) && (ch <= '9')) && (!97 <= ch <= 122))
			is_character = 1;
		++length;
	
		printf ("** %d **", ch);
	}
	while (ch != '\n');
	
	printf ("length:%d\n", length);
	
	if (is_length && is_upper && is_digit && is_lower_case && is_character)
	{
		printf("GOOD PASSWORD\n");
	}
	else
	{
		printf("BAD PASSWORD\n");
		if (is_length == 0)
			printf("Password len not between 8 and 14\n");
		if (is_lower_case == 0)
			printf("No lower case letters in password\n");
		if (is_character == 0)
		 printf("Password len not between 8 and 14\n");
	}
	
	return 0;
}
