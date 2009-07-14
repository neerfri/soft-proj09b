#########################################################
#    Main Makefile for final project
#            Software Project
#
#         Author: Neer Friedman
########################################################

SUBDIRS = "shared" "modularity_mat" "sparse_mlpl"

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