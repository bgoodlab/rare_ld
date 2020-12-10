#!/bin/bash

export type="selA"
for idx in $(python parameters.py idxs ${type}); do
    export params=$(python parameters.py get_params ${type} ${idx})
    echo ${params}
    mv output/output_${type}_${idx}.txt.gz output/output_r_${type}_${idx}.txt.gz
done
