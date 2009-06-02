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
#define ILLEGAL_LIST_OPERATION(description) printf(stderr, "ILLEGAL_LIST_OPERATION: %s", description); exit(EXIT_FAILURE);

/* function prototype */

/* This function creates an empty list, currently this just means returning NULL
 * seems pretty unuseful :-)
 **/
LINK linked_list_create_empty();

/* returns a new created list with one element containing the supplied data:
 * ex: LINK head = linked_list_create_with_data('a');
 */
LINK linked_list_create_with_data(DATA data);

/* returns a pointer to an element in the list that is found after
 * going 'index' times to the next element
 */
LINK *linked_list_pointer_to_index(LINK *head, int count);

/* returns the number of elements until the end of the list
 */
int linked_list_count(LINK head);

/* returns a pointer to the last element of the list
 */
LINK *linked_list_pointer_to_end(LINK *head);

/* inserets a list after 'element_before'
 */
void linked_list_insert_list(LINK *element_before, LINK list_to_insert);

/* do this method for every element in the list
 */
void linked_list_for_each(LINK head, void (*f)(LINK element));

/* finds the first element for which 'f' returns a non-zero value
 */ 
LINK linked_list_find(LINK head, int (*f)(LINK element, void *user_data),  void *user_data_to_func);

/* find the first element where 'ELEMENT->d == data'  
 */
LINK linked_list_find_first_data(LINK head, DATA data);

/* find the last occurrence of data or the last element of the list
 */
LINK linked_list_find_last_data_or_tail(LINK head, DATA data);

/* remove element LINK from the list returns the LINK to the following ptr (the one who took his place)
 */
LINK linked_list_remove(LINK *head, LINK to_remove);

/*
 * Methods for private use
 */
int prvt_linked_list_data_finder(LINK element, void *user_data);
void prvt_linked_list_data_printer(LINK element);
void prvt_linked_list_data_shifter(LINK element);


/*
 * REQUIRED METHODS FOR THE ASSIGNMENT
 */

/*insert 'to_insert' after the last occurrence of 'after' in the list
 * if after is not in the list, insert 'to_insert' at the end
 */
void InsertAfterToList(LINK *head, DATA after, DATA to_insert);

/* insert 'data' at the end of the list 
 */
void AppendToList(LINK *head, DATA data);

/* delete element with index 'index' from list 'head'
 */
void DeleteCharFromList(LINK *head, int index);

/* replace each occurrence of 'search' with 'replace'
 */
void ReplaceInList(LINK head, DATA search, DATA replace);

/* erase all occurrences of 'to_erase'
 */
void EraseFromList(LINK *head, DATA to_erase);

/* revereses the list
 */
void ReverseList(LINK *head);

/* prints the list
 */
void PrintList(LINK head);

/* delete the whole list
 */
void DeleteList(LINK *head);