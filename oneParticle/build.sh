#!/bin/bash
rm output.txt simu
g++ simu.cpp -o simu
./simu
python3 plot.py
