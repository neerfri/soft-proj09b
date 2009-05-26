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
	return empty_list;
}

LINK linked_list_jump_to_index(LINK head, int count) {
	LINK ptr = head;
	for(;count>0;count--) {
		if (head->next == NULL) {
			LIST_OUT_OF_RANGE_ERROR("linked_list_jump_forward")
		}
		ptr = ptr->next;
	}
	return ptr;
}

int linked_list_count(LINK head) {
	LINK ptr = head;
	int i;
	for(i = 0;ptr!=NULL; ptr = ptr->next) {	
	}
	return i;
}

LINK linked_list_jump_to_end(LINK head) {
	return linked_list_jump_to_index(head, linked_list_count(head));
}

void linked_list_insert_list(LINK element_before, LINK list_to_insert) {
	LINK element_after;
	assert(element_before!=NULL);
	element_after = element_before->next;
	element_before->next = list_to_insert;
	linked_list_jump_to_end(list_to_insert)->next = element_after;
}

void linked_list_insert_list_in_index(LINK head, LINK list_to_insert, int index) {
	LINK ptr = linked_list_jump_to_index(head, index);
	linked_list_insert_list(ptr, list_to_insert);
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
