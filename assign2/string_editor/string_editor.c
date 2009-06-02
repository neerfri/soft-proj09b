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
#include <ctype.h>
#include "my_list.h"

#define CHECK_ALPHA_NUMERIC_OR_EXIT(chr) if (!isalnum(chr)) printf("Illegal char enterd: '%c'\nAborting.", chr); if (!isalnum(chr)) exit(EXIT_FAILURE); 



void main_loop() {
	LINK *head;
	char chr;
	char c1, c2;
	int index;
	
	head = malloc(sizeof(LINK *));
	*head = linked_list_create_empty();
	
	while(1) {
		while((chr = tolower(getchar()))) {
			if(isalnum(chr)) {
				break;
			}
		}
		switch(chr) {
		case 'q':
			return;
		case 'i':
			c1 = getchar();
			c2 = getchar();
			CHECK_ALPHA_NUMERIC_OR_EXIT(c1)
			CHECK_ALPHA_NUMERIC_OR_EXIT(c2)
			InsertAfterToList(head, c1, c2);
			break;
		case 'a':
			c1 = getchar();
			CHECK_ALPHA_NUMERIC_OR_EXIT(c1)
			AppendToList(head, c1);
			break;
		case 'd':
			scanf("%d", &index);
			DeleteCharFromList(head, index);
			break;
		case 'r':
			c1 = getchar();
			c2 = getchar();
			CHECK_ALPHA_NUMERIC_OR_EXIT(c1)
			CHECK_ALPHA_NUMERIC_OR_EXIT(c2)
			ReplaceInList(*head, c1, c2);
			break;
		case 'e':
			c1 = getchar();
			CHECK_ALPHA_NUMERIC_OR_EXIT(c1)
			EraseFromList(head, c1);
			break;
		case 'm':
			ReverseList(head);
			break;
		case 'f':
			DeleteList(head);
			break;
		case 'p':
			printf("[");
			PrintList(*head);
			printf("]\n");
			break;
		case '\n':
			break;
		default:
			printf("unknown command '%c'.\n",chr);
		}
	}
}

int main(void) {
	main_loop();
	return EXIT_SUCCESS;
}
