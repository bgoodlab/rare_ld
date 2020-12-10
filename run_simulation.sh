#!/bin/bash

export type=$1
for idx in $(python parameters.py idxs ${type}); do
    export params=$(python parameters.py get_params ${type} ${idx})
    echo ${params}
    ./simulate_twolocus ${params} | gzip -c > output/output_${type}_${idx}.txt.gz
done
