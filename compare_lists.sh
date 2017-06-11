#!/bin/bash

if [ $# == 0 ]; then
    echo "Your command line contains no arguments"
    echo "Usage: compare content of two lists; What is in file2 but not in file1"
fi

awk 'FNR == NR {data[ $1 ];next;}{if ( ! ($1 in data) ) {print $0;}}' $1 $2




