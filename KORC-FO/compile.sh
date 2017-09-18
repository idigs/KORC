#!/bin/bash

clear

make clean

make info PRG_ENV=$1
make -B  PRG_ENV=$1
