#!/bin/bash

export DYLD_LIBRARY_PATH="$ELASTIX_PATH:$DYLD_LIBRARY_PATH"
transformix $1
