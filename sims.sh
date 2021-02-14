#!/bin/bash
# A sample Bash script
for i in {51..99..1}; 
do 
	python3 matrix_generation.py $i &
done 
for i in {51..99..1}; 
do 
	python3 metapop_strats.py $i c &
done 
for i in {51..99..1}; 
do 
	python3 metapop_strats.py $i duration &
done 



