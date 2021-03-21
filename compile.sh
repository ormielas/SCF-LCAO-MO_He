#!/bin/bash

for file in *.out; do
	if [ -f $file ]; then
		rm -r $file
		echo " $file removed"
	fi
done

gfortran main.f90 tred2.f90 tqli.f90 pqrsint.f90 frobenius.f90 -o main.out

