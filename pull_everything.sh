#!/bin/bash -x
#find . -type d -depth 1 -exec git --git-dir={}/.git --work-tree=$PWD/{} pull origin main \;
for DN in $(find . -type d -depth 1)
do
	git --git-dir=${DN}/.git --work-tree=${DN} pull origin main
done
