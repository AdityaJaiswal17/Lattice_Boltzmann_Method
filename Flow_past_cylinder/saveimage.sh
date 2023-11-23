#!/bin/bash

rm -rf images
mkdir images
for i in {1000..99000..1000}
do 

gnuplot -e "set terminal png size 1000,200; set output 'images/$i.png'; set pm3d map; set size ratio 1/10; sp 'output/output_t$i.dat' u 1:2:5 w image"

done
