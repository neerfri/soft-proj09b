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

void main_loop() {
	char chr;
	char c1, c2;
	int index;
	LINK *head;
	
	while((chr = tolower(getchar()))) {
		switch(chr)
		 {
		  case 'q':
		    return;
		  case 'i':
		    c1 = getchar();
		    c2 = getchar();
		    InsertAfterToList(head, c1, c2);
		    break;
		  case 'a':
		    c1 = getchar();
		    AppendToList(head, c1);
		  case 'd':
		    scanf("%d", &index);
		    DeleteCharFromList(head, index);
		  case 'p':
		    PrintList(*head);
		    break;
		  case '\n':
		    break;
		  default:
		    printf("unknown command.\n");
		 }
		 /*getchar();  to waste the new line char */
	}
}

int main(void) {
	/* 
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
	*/
	LINK *head;
	char c1='a',c2='b',c3='c',c4='d',c5='e',c6='f';
	*head = linked_list_create_empty();
	/*DeleteCharFromList(*head, 0);*/
	AppendToList(head, c1);
	InsertAfterToList(*head, c1, c3);
	InsertAfterToList(*head, c1, c2);
	AppendToList(head, c4);
	/*DeleteCharFromList(head, 1);*/
	DeleteList(head);
	AppendToList(head, c5);
	/*
	InsertAfterToList(*head, c1, c2);
	ReplaceInList(*head, 'd', 'b');
	ReplaceInList(*head, 'b', 'f');
	*/
	  
	/*printf("linked_list_count: %d\n", linked_list_count(head));*/
	/*head = linked_list_create_with_data(c1);*/
	/*AppendToList(head, c2);*/
	/*EraseFromList(head, 'f');
	ReverseList(head);
	*/
	PrintList(*head);
	printf("\n");
	/* DeleteList(*head);*/
	main_loop();
	return EXIT_SUCCESS;
}
