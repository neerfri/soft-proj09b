/*
//============================================================================
// Name        : my_list.h
// Author      : Neer Friedman
// Version     :
// Copyright   : Your copyright notice
//============================================================================
*/

typedef   char   DATA;           /* will use char in examples */
extern void *malloc(size_t size);

struct linked_list {
   DATA                 d;
   struct linked_list   *next;
};

typedef   struct linked_list   ELEMENT;
typedef   ELEMENT              *LINK;

/* MACROs */
#define LIST_OUT_OF_RANGE_ERROR(method_name) printf("Error allocating memory for empty list in %s\n", method_name); exit(EXIT_FAILURE);

/* function prototype */

/* This function creates an empty list, currently this just means returning NULL
 * seems pretty unuseful :-)
 **/
LINK linked_list_create_empty();

/* returns a new created list with one element containing the supplied data:
 * ex: LINK head = linked_list_create_with_data('a');
 */
LINK linked_list_create_with_data(DATA data);

/* returns a pointer to an element in the list that is found after going 'index'
 * times to the next element
 */
LINK linked_list_jump_to_index(LINK head, int count);

/* returns the number of elements until the end of the list
 */
int linked_list_count(LINK head);

/* returns a pointer to the last element of the list
 */
LINK linked_list_jump_to_end(LINK head);

/* inserets a list after 'element_before'
 */
void linked_list_insert_list(LINK element_before, LINK list_to_insert);

/* inserets a list 'index' elements after 'head'
 */
void linked_list_insert_list_in_index(LINK head, LINK list_to_insert, int index);

/* do this method for every element in the list
 */
void linked_list_for_each(LINK head, void (*f)(LINK element));

/* finds the first element for which 'f' returns a non-zero value
 */ 
LINK linked_list_find(LINK head, int (*f)(LINK element, void *user_data),  void *user_data_to_func);
