#!/bin/bash
shopt -s globstar
for file in **/*\ *
do 
    mv "$file" "${file// /_}"       
done