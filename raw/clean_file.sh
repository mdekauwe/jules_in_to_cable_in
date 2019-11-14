#!/bin/bash

#gawk '$1=$1' Loobos_1997.dat > tmp; mv tmp Loobos_1997.dat
sed 's/  \+/ /g' Loobos_1997.dat > tmp; mv tmp Loobos_1997.dat
