/*
//============================================================================
// Name        : string_editor.c
// Author      : Neer Friedman
// Version     :
// Copyright   : Your copyright notice
// Description : A String Editing Utility
//============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include "my_list.h"

void print_element_data(LINK element) {
	printf("%c", element->d);
}

int char_finder(LINK element, void *user_data) {
	char *d = (char *)user_data;
	return((element->d == *d ? 1 : 0));
}

int main(void) {
	LINK head, found;
	char chr = 'c';
	printf("Hello string editor !\n");
	head = linked_list_create_with_data('a');
	linked_list_insert_list_in_index(head, linked_list_create_with_data('c'), 0);
	linked_list_insert_list_in_index(head, linked_list_create_with_data('d'), 1);
	linked_list_insert_list_in_index(head, linked_list_create_with_data('b'), 0);
	linked_list_for_each(head, print_element_data);
	printf("\n");
	found = linked_list_find(head, char_finder, &chr);
	linked_list_for_each(found, print_element_data);
	printf("\n");
	return EXIT_SUCCESS;
}
