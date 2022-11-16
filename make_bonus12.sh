#!/bin/sh
gfortran -c cl_option.f
gfortran -fno-align-commons -fopenmp  -o bonus12 cl_option.f bonus12.f epgenb.f rc_modb.f sampproton.f dn.f wfn.f pinterp.f
#!
