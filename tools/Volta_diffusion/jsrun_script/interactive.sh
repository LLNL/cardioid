#!/bin/bash
#
# sample lsf bsub to run an interactive job, optionally on a selected host.
#
# pick a host to land on.
#host=${1:-fstgb001}

#
# the -Is says you want an interactive session
#               the s says you want a terminal session.
#
# shared_int is the "shared interactive queue"
if [ -z $LSB_BATCH_JID ]; then
    set -x
    bsub  \
        -R "select[ngpus=4] rusage[ngpus_shared=20]" \
        -env "LSB_START_JOB_MPS=N" \
        -Is \
        -n 1 \
        -q excl_int \
        /bin/bash
fi
