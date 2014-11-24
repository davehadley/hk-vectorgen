#!/bin/bash

scriptdir=$(dirname $(readlink -f $0));
srcdir=$(readlink -f $scriptdir/src);

export PYTHONPATH=${srcdir}:${PYTHONPATH}

#echo ${PYTHONPATH}

