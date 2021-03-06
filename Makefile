#########################################################
#    Main Makefile for eclipse project
#        Software Project Assigns
#
#         Author: Neer Friedman
########################################################

APPROXIMATE_PI = assign1/approximate_pi
CHECK_POLICY = assign1/check_policy
IS_DIVIDABLE = assign1/is_dividable
STRING_EDITOR = assign2/string_editor
FINAL_PROJECT = final_project


SUBDIRS = $(APPROXIMATE_PI) $(CHECK_POLICY) $(IS_DIVIDABLE) $(STRING_EDITOR) $(FINAL_PROJECT)

# Default Target
all:
	@echo 
	@echo "#######################################"
	@echo "### BUILDING ALL TARGETS ###"
	@echo "#######################################"
	@echo 
	
	for i in $(SUBDIRS) ; do ( cd $$i ; make ) ; done
	
clean:
	@echo 
	@echo "#######################################"
	@echo "### CLEANING ALL TARGETS ###"
	@echo "#######################################"
	@echo 
	
	for i in $(SUBDIRS) ; do ( cd $$i ; make clean ) ; done