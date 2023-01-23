#!/bin/bash
set -e

image="singularity/mktest.sif"
def="singularity/Singularity"

singularity build -s ${image} ${def}
