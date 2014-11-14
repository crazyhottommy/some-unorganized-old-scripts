#!/bin/bash

files=`cat file_list` #Note the use of backtick to execute
		      #a command and put the result into a
		      #variable--` is key to left of 1

for i in $files
do
	cp ~/pcfb/examples/$i ~/pcfb/sandbox/
done
