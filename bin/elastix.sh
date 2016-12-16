#!/bin/bash

#Usage:elastix.sh <elastix arguments>  [log-file-path]

if [ -z "$1"  ] 
then
  echo "Error: No input parameters to run elastix "  1>&2
  exit 1
fi

if ! [ -z "$2" ]
then 
  OUT_LOG_FILE="$2"/log-out.txt 
  ERR_LOG_FILE="$2"/log-err.txt 
else
  OUT_LOG_FILE=log-out.txt 
  ERR_LOG_FILE=log-err.txt 
fi

ELASTIX_TEST=$(which elastix)
if [ "$?" -eq "1" ]
then
  if [ -z "$ELASTIX_PATH" ]
  then
    echo "Error: The environment variable ELASTIX_PATH is not set" 1>&2
    exit 1
  fi    
else
  ELASTIX_PATH=$(dirname $ELASTIX_TEST | sed 's/bin/lib/') 
fi

export DYLD_LIBRARY_PATH="$ELASTIX_PATH:$DYLD_LIBRARY_PATH"

elastix $1   >$OUT_LOG_FILE 2>$ERR_LOG_FILE
