#!/usr/bin/env python
from os import system,popen
import string
import argparse


parser = argparse.ArgumentParser(description="Takes in a file generated by Adam's symmetric matcher mover and creates a constraints file for HID")
parser.add_argument('PDB', type=str, nargs=1, help="score file to grab score type headers")
args = parser.parse_args()
input_file = args.PDB[0]


scores_array = map(string.split,open(input_file,'r').readlines())
output_string = ""
for score in scores_array[1]:
    output_string += "'" + score + "':[], "

print(output_string)

