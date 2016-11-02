#!/bin/bash


echo `date +%Y-%m-%d` >> log.txt
for model in ../tests/models/*

do
    ./qh.out $model >> log.txt
done
