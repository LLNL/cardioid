#!/bin/bash 

echo cells
grep cells $1 | awk '{ print $5 }'

echo time
grep 'us sec' $1 | awk '{ print $3 }'

