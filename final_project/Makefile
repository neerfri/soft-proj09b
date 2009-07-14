#########################################################
#    Main Makefile for final project
#            Software Project
#
#         Author: Neer Friedman
########################################################

SUBDIRS = "shared" "modularity_mat" "sparse_mlpl" "eigen_pair"

# Default Target
all:
	@echo 
	@echo "########################################"
	@echo "### BUILDING ALL TARGETS FOR FINAL PROJ#"
	@echo "########################################"
	@echo 
	
	for i in $(SUBDIRS) ; do ( cd $$i ; make ) ; done
	
clean:
	@echo 
	@echo "#######################################"
	@echo "### CLEANING ALL FOR FINAL PROJ     ###"
	@echo "#######################################"
	@echo 
	
	for i in $(SUBDIRS) ; do ( cd $$i ; make clean ) ; done