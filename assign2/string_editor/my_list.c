/*
//============================================================================
// Name        : my_list.c
// Author      : Neer Friedman
// Version     :
// Copyright   : Your copyright notice
//============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "my_list.h"



LINK linked_list_create_empty() {
	return NULL;
}

LINK linked_list_create_with_data(DATA data) {
	LINK empty_list = (ELEMENT *)malloc(sizeof(ELEMENT));
	if (empty_list == NULL) {
		printf("Error allocating memory for empty list");
		exit(EXIT_FAILURE);
	}
	empty_list->d = data;
	empty_list->next = NULL;
	return empty_list;
}

LINK *linked_list_pointer_to_index(LINK *head, int count) {
	LINK *ptr = head;
	assert(count < linked_list_count(*ptr));
	for(;count>0;count--) {
		ptr = &((*ptr)->next);
	}
	return ptr;
}

int linked_list_count(LINK head) {
	LINK ptr = head;
	int i;
	for(i = 0; ptr; ptr = ptr->next) {
		i++;
	}
	return i;
}

LINK *linked_list_pointer_to_end(LINK *head) {
	return linked_list_pointer_to_index(head, linked_list_count(*head)-1);
}

void linked_list_insert_list(LINK *element_before, LINK list_to_insert) {
	LINK element_after;
	if (list_to_insert == NULL)
		return;
	if ((*element_before) == NULL) {
		*element_before = list_to_insert;
	} else {
		element_after = (*element_before)->next;
		(*element_before)->next = list_to_insert;
		(*linked_list_pointer_to_end(&list_to_insert))->next = element_after;
	}
}

void linked_list_for_each(LINK head, void (*f)(LINK element)) {
	LINK ptr = head;
		int i;
		for(i = 0;ptr!=NULL; ptr = ptr->next) {
			f(ptr);
		}
}

LINK linked_list_find(LINK head, int (*f)(LINK element, void *user_data),  void *user_data_to_func) {
	LINK ptr = head;
	int i;
	for(i = 0;ptr!=NULL; ptr = ptr->next) {
		if (f(ptr, user_data_to_func)) {
			return ptr;
		}
	}
	
	return NULL;
}

LINK linked_list_find_first_data(LINK head, DATA data) {
	return(linked_list_find(head, prvt_linked_list_data_finder, &data));
}

LINK linked_list_find_last_data_or_tail(LINK head, DATA data) {
	LINK ptr, ret = NULL;
	for(ptr=head;(ptr = linked_list_find_first_data(ptr, data)) != NULL; ptr = ptr->next) {
		ret = ptr;
	}
	return ret;
}

LINK linked_list_remove(LINK *head, LINK to_remove) {
	LINK ptr;
	if ((*head) == NULL) {
		return NULL;
	}
	if ((*head) == to_remove) {
		ptr = *head;
		*head = (*head)->next;
		free(ptr);
		return(*head);
	} else {
		return(linked_list_remove(&((*head)->next), to_remove));
	}
}

/* helper methods */
int prvt_linked_list_data_finder(LINK element, void *user_data) {
	DATA *d = (DATA *)user_data;
	return((element->d == *d ? 1 : 0));
}

void prvt_linked_list_data_printer(LINK element) {
	printf("%c", element->d);
}

void prvt_linked_list_data_shifter(LINK element) {
	if (element->next) {
		element->d = element->next->d;
		if (element->next->next == NULL) {
			/*
			 *   Ex:  a -> b -> c -> ||
			 *   If we want to remove 'b' (a.k.a element) we shift the 'c' char to it,
			 *   free 'c' and set element->next to NULL
			 *   This will always happen since we always run till the end of the list 
			 */
			free(element->next);
			element->next = NULL;
		}
	}
	
}


/*
 * REQUIRED METHODS FOR THE ASSIGNMENT
 */

void InsertAfterToList(LINK head, DATA after, DATA to_insert) {
	LINK element_before = linked_list_find_last_data_or_tail(head, after);
	linked_list_insert_list(&element_before, linked_list_create_with_data(to_insert));
}

void AppendToList(LINK *head, DATA data) {
	LINK new_list = linked_list_create_with_data(data);
	LINK *tail = linked_list_pointer_to_end(head);
	linked_list_insert_list(tail, new_list);
}

void DeleteCharFromList(LINK head, int index) {
	LINK element_to_remove = *linked_list_pointer_to_index(&head, index);
	linked_list_for_each(element_to_remove, prvt_linked_list_data_shifter);
}

void ReplaceInList(LINK head, DATA search, DATA replace) {
	LINK ptr = linked_list_find_first_data(head, search);
	while(ptr) {
		ptr->d = replace;
		ptr = linked_list_find_first_data(ptr->next, search);
	}
}

void EraseFromList(LINK *head, DATA to_erase) {
	LINK ptr;
	for(ptr = linked_list_find_first_data(*head, to_erase); ptr; ptr = linked_list_find_first_data(ptr, to_erase)) {
		ptr = linked_list_remove(head, ptr); 
	}
}

void ReverseList(LINK *head) {
	LINK ptr, new_list, temp_ptr;
	new_list = linked_list_create_empty();
	for(ptr = *head; ptr; ptr=ptr->next) {
		temp_ptr = linked_list_create_with_data(ptr->d);
		linked_list_insert_list(&temp_ptr, new_list);
		new_list = temp_ptr;
	}
	ptr = *head;
	*head = new_list;
	DeleteList(ptr);
}

void PrintList(LINK head) {
	linked_list_for_each(head, prvt_linked_list_data_printer);
}

void DeleteList(LINK head) {
	if (head != NULL) {
		DeleteList(head->next);
		free(head);
	}
}