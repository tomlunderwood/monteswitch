# Makefile for the 'monteswitch' package
# Author: Tom L Underwood


# README: * Adjust the variables 'FC', 'FC_MPI', 'FC_FLAGS', and 'FC_MPI_FLAGS' to suit your
#           requirements before 'make' is invoked.
#         * The 'BINS' variable is the list of binaries to be compiled. Delete 'monteswitch_mpi'
#           if you do not have MPI, and hence will be unable to compile 'monteswitch_mpi'
#         * Type 'make help' for a list of callable targets.




# A note on .mod file rules:
#
# The rules pertaining to .mod files 'touch' the .mod file once it has been
# created or updated. This is to address the fact that gfortran - which I primarily used -
# only updates the timestamp on .mod files if they have been modified. However, gfortran started
# doing this to address a bug, but for me it makes sense to update the timestamp,
# otherwise make will remake things which it doesn't need to.


#
# Variables
#

# The fortran compiler (consider gfortran)
FC = gfortran
# The mpi fortran compiler
FC_MPI = mpif90
# Extra flags to be used in the compiling using the mpi compiler (consider -O3, -static, -pg, -g -fcheck=all -fbounds-check)
FC_FLAGS =   -Wall -fbounds-check -O3
# Extra flags to be used in the compiling using the mpi fortran compiler
FC_MPI_FLAGS =  -Wall -fbounds-check




#
# Macros
#

# List of binaries to be compiled
BINS = monteswitch monteswitch_post monteswitch_mpi lattices_in_hcp_fcc
# List of html documents to be compiled
DOCS = kinds_mod_docs.html metropolis_mod_docs.html rng_mod_docs.html \
monteswitch_mod_docs.html interactions_docs.html monteswitch_docs.html \
monteswitch_post_docs.html monteswitch_mpi_docs.html \
lattices_in_hcp_fcc_docs.html




#
# Utility targets
#

.PHONY: all help clean

# Default target
# target: all - Makes all binaries and docs (default target).
all: $(BINS) $(DOCS)

# target: help - Display callable targets (obtained from the Malefile itself)
help:
	@echo "Here is a list of callable targets: "; egrep "^# target:" Makefile | sed "s/# target://"

# target: clean - Removes all binaries, object files, module files, and html documents
clean:
	rm -f *.o *.mod $(BINS) $(DOCS)




#
# Real targets
#

# lattices_in_hcp_fcc

# target: lattices_in_hcp_fcc_docs.html
lattices_in_hcp_fcc_docs.html: lattices_in_hcp_fcc.f95
	./tomdocs.sh '\!\! ?' $< > $@

# target: lattices_in_hcp_fcc
lattices_in_hcp_fcc: lattices_in_hcp_fcc.o kinds_mod.o
	$(FC) $(FC_FLAGS) -o $@ $^

# target: lattices_in_hcp_fcc.o
lattices_in_hcp_fcc.o: lattices_in_hcp_fcc.f95 kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<


# monteswitch_mpi

# target: monteswitch_mpi_docs.html
monteswitch_mpi_docs.html: monteswitch_mpi.f95
	./tomdocs.sh '\!\! ?' $< > $@

# target: monteswitch_mpi
monteswitch_mpi: monteswitch_mpi.o monteswitch_mod.o kinds_mod.o rng_mod.o mersenne_twister.o \
metropolis_mod.o
	$(FC_MPI) $(FC_MPI_FLAGS) -o $@ $^

# target: monteswitch_mpi.o
monteswitch_mpi.o: monteswitch_mpi.f95 monteswitch_mod.mod kinds_mod.mod rng_mod.mod \
mersenne_twister.mod metropolis_mod.mod
	$(FC_MPI) $(FC_MPI_FLAGS) -c $<


# monteswitch

# target: monteswitch_docs.html
monteswitch_docs.html: monteswitch.f95
	./tomdocs.sh '\!\! ?' $< > $@

# target: monteswitch
monteswitch: monteswitch.o monteswitch_mod.o kinds_mod.o rng_mod.o mersenne_twister.o \
metropolis_mod.o
	$(FC) $(FC_FLAGS) -o $@ $^

# target: monteswitch.o
monteswitch.o: monteswitch.f95 monteswitch_mod.mod kinds_mod.mod rng_mod.mod \
mersenne_twister.mod
	$(FC) $(FC_FLAGS) -c $<


# monteswitch_post

# target: monteswitch_post_docs.html
monteswitch_post_docs.html: monteswitch_post.f95
	./tomdocs.sh '\!\! ?' $< > $@

# target: monteswitch_post
monteswitch_post: monteswitch_post.o monteswitch_mod.o kinds_mod.o rng_mod.o mersenne_twister.o \
metropolis_mod.o
	$(FC) $(FC_FLAGS) -o $@ $^

# target: monteswitch_post.o
monteswitch_post.o: monteswitch_post.f95 monteswitch_mod.mod kinds_mod.mod rng_mod.mod \
mersenne_twister.mod
	$(FC) $(FC_FLAGS) -c $<


# monteswitch module and interactions

# target: monteswitch_mod_docs.html
monteswitch_mod_docs.html: monteswitch_mod.f95
	./tomdocs.sh '\!\! ?' $< > $@

# target: interactions_docs.html
interactions_docs.html: interactions.f95
	./tomdocs.sh '\!\! ?' $< > $@

# target: monteswitch_mod.o
monteswitch_mod.o: monteswitch_mod.f95 interactions.f95 kinds_mod.mod rng_mod.mod mersenne_twister.mod \
metropolis_mod.mod
	$(FC) $(FC_FLAGS) -c $<

# target: monteswitch_mod.mod
monteswitch_mod.mod: monteswitch_mod.f95 interactions.f95 kinds_mod.mod rng_mod.mod mersenne_twister.mod \
metropolis_mod.mod
	$(FC) $(FC_FLAGS) -c $<
	touch $@


# metropolis module

# target: metropolis_mod_docs.html
metropolis_mod_docs.html: metropolis_mod.f95
	./tomdocs.sh '\!\! ?' $< > $@

# target: metropolis_mod.o
metropolis_mod.o: metropolis_mod.f95 kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<

# target: metropolis_mod.mod
metropolis_mod.mod: metropolis_mod.f95 kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<
	touch $@


# rng module, mersenne_twister

# target: rng_mod_docs.html
rng_mod_docs.html: rng_mod.f95
	./tomdocs.sh '\!\! ?' $< > $@

# target: rng_mod.o
rng_mod.o: rng_mod.f95 mersenne_twister.mod kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<

# target: rng_mod.mod
rng_mod.mod: rng_mod.f95 mersenne_twister.mod kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<
	touch $@

# target: mersenne_twister.o
mersenne_twister.o: mersenne_twister.f90 kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<

# target: mersenne_twister.mod
mersenne_twister.mod: mersenne_twister.f90 kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<
	touch $@


# kinds module

# target: kinds_mod_docs.html
kinds_mod_docs.html: kinds_mod.f95
	./tomdocs.sh '\!\! ?' $< > $@

# target: kinds_mod.o
kinds_mod.o: kinds_mod.f95
	$(FC) $(FC_FLAGS) -c $<

# target: kinds_mod.mod
kinds_mod.mod: kinds_mod.f95
	$(FC) $(FC_FLAGS) -c $<
	touch $@

