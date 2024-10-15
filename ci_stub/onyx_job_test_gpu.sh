#!/bin/bash

env | grep JOB_COUNT

JOB_COUNT=1

cd ./KORC/build_gpu && ctest -j ${JOB_COUNT:-1} --output-on-failure

