#!/bin/bash
# compiling idea from:
# https://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules

# alias
    #compile='gfortran -Og -Wall -Wextra -fcheck=all -fbacktrace'
    compile='gfortran -Ofast'

# directories
    SRC=src
    OBJ=obj
    BIN=bin

# reads
    $compile -o $BIN/read_zero $SRC/read_zero.f90
    $compile -o $BIN/read_one $SRC/read_one.f90

# total_infected
    $compile -o $BIN/total_infected_sir $SRC/total_infected_sir.f90

# mt generator
    $compile -J$OBJ -c $SRC/generator.f90 -o $OBJ/generator.o

# subroutines
    $compile -J$OBJ -c $SRC/subroutines.f90 -o $OBJ/subroutines.o

# sis
    $compile -o $BIN/sis $SRC/sis.f90

# sis_1event
    $compile -I$OBJ -o $BIN/sis_1event $OBJ/subroutines.o $SRC/sis_1event.f90

# sis endemic
    $compile -I$OBJ -o $BIN/sis_endemic $OBJ/subroutines.o $SRC/sis_endemic.f90

# sir
    $compile -o $BIN/sir $SRC/sir.f90

# sirs
    $compile -I$OBJ -o $BIN/sirs $OBJ/subroutines.o $SRC/sirs.f90

# sir dynamic
    $compile -I$OBJ -o $BIN/sir_dynamic $OBJ/subroutines.o $SRC/sir_dynamic.f90

# sir_mut
    $compile -I$OBJ -o $BIN/sir_mut $OBJ/subroutines.o $SRC/sir_mut.f90

# sir_mut_prev
    $compile -I$OBJ -o $BIN/sir_mut_prev $OBJ/subroutines.o $SRC/sir_mut_prev.f90

# sir_mut_periodic
    $compile -I$OBJ -o $BIN/sir_mut_periodic $OBJ/subroutines.o $SRC/sir_mut_periodic.f90

# sir_mut_joined
    $compile -I$OBJ -o $BIN/sir_mut_joined $OBJ/subroutines.o $SRC/sir_mut_joined.f90

# sir_mut_i
    $compile -I$OBJ -o $BIN/sir_mut_i $OBJ/subroutines.o $SRC/sir_mut_i.f90

# prevalence
    $compile -o $BIN/prevalence $SRC/prevalence.f90

# sir_mut_i_inf_prev
    $compile -I$OBJ -o $BIN/sir_mut_i_inf_prev $OBJ/subroutines.o $SRC/sir_mut_i_inf_prev.f90

# histo_inf_prev
    $compile -o $BIN/histo_inf_prev $SRC/histo_inf_prev.f90
