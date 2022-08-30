#!/usr/bin/env bash
# Zip up local files.
DTE=`date +"%Y-%m-%d"`
FNAME=random_${DTE}.tar.gz
echo "Zip to ${FNAME} ..."
tar --exclude=init_*.sh --exclude=*.gp --exclude=random_stream --exclude=*.tar.gz --exclude=*.o --exclude=*.d  --exclude=._* --exclude=.DS_Store -czvf ${FNAME} *
