# Makefile for the 'monteswitch' package
# Author: Tom L Underwood


# *************** README ***************
#
# Installation instructions:
# 1) Adjust the variables 'FC', 'FC_MPI', 'FC_FLAGS', and 'FC_MPI_FLAGS' (described below)
#    to suit your system before 'make' is invoked. If MPI is not available then don't worry,
#    an option is available to compile only the serial binaries.
# 2) Type 'make help' for a list and description of callable targets, and invoke the most
#    appropriate one.

# The fortran compiler
FC = gfortran
# The mpi fortran compiler
FC_MPI = mpif90
# Extra flags to be used in the compiling using the mpi compiler
FC_FLAGS = -O3
# Extra flags to be used in the compiling using the mpi fortran compiler
FC_MPI_FLAGS = -O3

# ************ END OF README ***********


# A note on .mod file rules:
# The rules pertaining to .mod files 'touch' the .mod file once it has been
# created or updated. This is to address the fact that gfortran - which I primarily used -
# only updates the timestamp on .mod files if they have been modified. However, gfortran started
# doing this to address a bug, but for me it makes sense to update the timestamp,
# otherwise make will remake things which it doesn't need to.


# MACROS

# List of binaries to be compiled with the 'serial' target
BINS_SERIAL = monteswitch monteswitch_post lattices_in_hcp_fcc lattices_in_bcc_fcc lattices_in_bcc_hcp
# List of binaries to be compiled with the 'mpi' target
BINS_MPI = monteswitch_mpi monteswitch monteswitch_post \
lattices_in_hcp_fcc lattices_in_bcc_fcc lattices_in_bcc_hcp
# List of html documents to be compiled ('srcdocs' target)
DOCS = kinds_mod_docs.html metropolis_mod_docs.html rng_mod_docs.html \
monteswitch_mod_docs.html monteswitch_docs.html monteswitch_post_docs.html monteswitch_mpi_docs.html \
lattices_in_hcp_fcc_docs.html lattices_in_bcc_fcc_docs.html lattices_in_bcc_hcp_docs.html \
interactions_EAM_docs.html interactions_EC_docs.html interactions_EC_NPT_docs.html \
interactions_LJ_docs.html interactions_LJ_hcp_fcc_docs.html interactions_pair_table_docs.html


# UTILITY TARGETS

.PHONY: all help clean

# target: help - List targets (default target).
help:
	@echo "Here is a list of targets: "; egrep "^# target:" Makefile | sed "s/# target://"

# target: serial - Make serial (non-MPI) binaries only
serial: $(BINS_SERIAL)
	$(MAKE) clean

# target: mpi - Make serial and MPI binaries
mpi: $(BINS_MPI)
	$(MAKE) clean

# target: srcdocs - Make html source code documentation
srcdocs: $(DOCS)

# target: uninstall - Remove all binaries and html source code documentation
uninstall:
	rm -f *.o *.mod $(BINS_SERIAL) $(BINS_MPI) $(DOCS)

# Remove object files and module files
clean:
	rm -f *.o *.mod


# REAL TARGETS

# interactions files

interactions_EAM_docs.html: interactions_EAM.f95
	sh docmaker.sh '\!\! ?' $< > $@

interactions_EC_docs.html: interactions_EC.f95
	sh docmaker.sh '\!\! ?' $< > $@

interactions_EC_NPT_docs.html: interactions_EC_NPT.f95
	sh docmaker.sh '\!\! ?' $< > $@

interactions_LJ_docs.html: interactions_LJ.f95
	sh docmaker.sh '\!\! ?' $< > $@

interactions_LJ_hcp_fcc_docs.html: interactions_LJ_hcp_fcc.f95
	sh docmaker.sh '\!\! ?' $< > $@

interactions_pair_table_docs.html: interactions_pair_table.f95
	sh docmaker.sh '\!\! ?' $< > $@


# lattices_in_bcc_hcp

lattices_in_bcc_hcp_docs.html: lattices_in_bcc_hcp.f95
	sh docmaker.sh '\!\! ?' $< > $@

lattices_in_bcc_hcp: lattices_in_bcc_hcp.o kinds_mod.o
	$(FC) $(FC_FLAGS) -o $@ $^

lattices_in_bcc_hcp.o: lattices_in_bcc_hcp.f95 kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<


# lattices_in_bcc_fcc

lattices_in_bcc_fcc_docs.html: lattices_in_bcc_fcc.f95
	sh docmaker.sh '\!\! ?' $< > $@

lattices_in_bcc_fcc: lattices_in_bcc_fcc.o kinds_mod.o
	$(FC) $(FC_FLAGS) -o $@ $^

lattices_in_bcc_fcc.o: lattices_in_bcc_fcc.f95 kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<


# lattices_in_hcp_fcc

lattices_in_hcp_fcc_docs.html: lattices_in_hcp_fcc.f95
	sh docmaker.sh '\!\! ?' $< > $@

lattices_in_hcp_fcc: lattices_in_hcp_fcc.o kinds_mod.o
	$(FC) $(FC_FLAGS) -o $@ $^

lattices_in_hcp_fcc.o: lattices_in_hcp_fcc.f95 kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<


# monteswitch_mpi

monteswitch_mpi_docs.html: monteswitch_mpi.f95
	sh docmaker.sh '\!\! ?' $< > $@

monteswitch_mpi: monteswitch_mpi.o monteswitch_mod.o kinds_mod.o rng_mod.o mersenne_twister.o \
metropolis_mod.o
	$(FC_MPI) $(FC_MPI_FLAGS) -o $@ $^

monteswitch_mpi.o: monteswitch_mpi.f95 monteswitch_mod.mod kinds_mod.mod rng_mod.mod \
mersenne_twister.mod metropolis_mod.mod
	$(FC_MPI) $(FC_MPI_FLAGS) -c $<


# monteswitch

monteswitch_docs.html: monteswitch.f95
	sh docmaker.sh '\!\! ?' $< > $@

monteswitch: monteswitch.o monteswitch_mod.o kinds_mod.o rng_mod.o mersenne_twister.o \
metropolis_mod.o
	$(FC) $(FC_FLAGS) -o $@ $^

monteswitch.o: monteswitch.f95 monteswitch_mod.mod kinds_mod.mod rng_mod.mod \
mersenne_twister.mod
	$(FC) $(FC_FLAGS) -c $<


# monteswitch_post

monteswitch_post_docs.html: monteswitch_post.f95
	sh docmaker.sh '\!\! ?' $< > $@

monteswitch_post: monteswitch_post.o monteswitch_mod.o kinds_mod.o rng_mod.o mersenne_twister.o \
metropolis_mod.o
	$(FC) $(FC_FLAGS) -o $@ $^

monteswitch_post.o: monteswitch_post.f95 monteswitch_mod.mod kinds_mod.mod rng_mod.mod \
mersenne_twister.mod
	$(FC) $(FC_FLAGS) -c $<


# monteswitch module and interactions (no html docs for interactions)

monteswitch_mod_docs.html: monteswitch_mod.f95
	sh docmaker.sh '\!\! ?' $< > $@

monteswitch_mod.o: monteswitch_mod.f95 interactions.f95 kinds_mod.mod rng_mod.mod mersenne_twister.mod \
metropolis_mod.mod
	$(FC) $(FC_FLAGS) -c $<

monteswitch_mod.mod: monteswitch_mod.f95 interactions.f95 kinds_mod.mod rng_mod.mod mersenne_twister.mod \
metropolis_mod.mod
	$(FC) $(FC_FLAGS) -c $<
	touch $@


# metropolis module

metropolis_mod_docs.html: metropolis_mod.f95
	sh docmaker.sh '\!\! ?' $< > $@

metropolis_mod.o: metropolis_mod.f95 kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<

metropolis_mod.mod: metropolis_mod.f95 kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<
	touch $@


# rng module, mersenne_twister

rng_mod_docs.html: rng_mod.f95
	sh docmaker.sh '\!\! ?' $< > $@

rng_mod.o: rng_mod.f95 mersenne_twister.mod kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<

rng_mod.mod: rng_mod.f95 mersenne_twister.mod kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<
	touch $@

mersenne_twister.o: mersenne_twister.f90 kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<

mersenne_twister.mod: mersenne_twister.f90 kinds_mod.mod
	$(FC) $(FC_FLAGS) -c $<
	touch $@


# kinds module

kinds_mod_docs.html: kinds_mod.f95
	sh docmaker.sh '\!\! ?' $< > $@

kinds_mod.o: kinds_mod.f95
	$(FC) $(FC_FLAGS) -c $<

kinds_mod.mod: kinds_mod.f95
	$(FC) $(FC_FLAGS) -c $<
	touch $@

