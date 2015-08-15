#!/bin/bash
#
# Author: Tom Underwood
#
# This script is designed for generating documentation from source 
# code or something similar. Lines beginning with matches to the
# regular expression of the first argument are extracted from the file
# specified in the second argument, and printed without the
# corresponding match or any leading whitespace.
#

# Check script arguments
if [ $# != 2 ]; then
    echo "$0: This script requires 2 arguments - a regular expression and a file."
    exit 1
fi
if ! [ -e $2 ]; then
    echo "$0: File specified in argument does not exist."
    exit 1
fi

# This command selects lines beginning with the regular expression $1
# (ignoring leading whitespace), and prints them - after removing the
# match corresponding to the regular expression.
egrep "^[[:space:]]*$1" $2 | sed -E "s/^[[:space:]]*$1//"

