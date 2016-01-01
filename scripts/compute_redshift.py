#!/usr/bin/env python

import numpy as np
import argparse
import sys


def find_nearest_value(array,value):
        idx = (np.abs(array-value)).argmin()
        return array[idx]

def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return idx

def format(value):
        return "%.6e" % value


# arguments handling

parser = argparse.ArgumentParser(
    prog='COMPUTE REDSHIFT',
    description='computes the redhift from values of "a" and "H" within "DATA_FILE"',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('data_file', metavar='DATA_FILE', type=str,
    help='the data file')
args = parser.parse_args()

with open(args.data_file,'r') as data_file:
        H0_line = data_file.readline()
        H0_line = H0_line.strip()
        H0_str = H0_line[H0_line.find('=')+1:]
        H0 = float(H0_str)
        header = data_file.readline()
        names = header.split()
        data = np.loadtxt(args.data_file, skiprows=2)

variables = {}
for j in range(len(names)):
        variables[names[j]] = data[:,j]

if 'H' in list(variables.keys()):
        index_a0 = find_nearest(variables['H'],H0)
        H0_nearest = find_nearest_value(variables['H'],H0)
        print ' nearest value to H0 =', H0,'is H =', H0_nearest
        print ' computing the redshift using a0 corresponding to this value'
        if 'a' in list(variables.keys()):
                a0 = variables['a'][index_a0]
                z = a0/variables['a']-1
                print ' this corresponds to a =', a0
        else:
                sys.exit('The input file does not contain any instance of "a"')
        if 't_sync' in list(variables.keys()):
                t_sync_a0 = variables['t_sync'][index_a0]
                print ' this corresponds to a value of the synchronous time of', t_sync_a0, '(Gyr)'
else:
        sys.exit('The input file does not contain any instance of "H"')


variables['z'] = z

with open(args.data_file,'w') as output_file:
        output_file.write(H0_line+'\n')
        for item in variables:
                output_file.write('%15s' % item)
        output_file.write('\n')
        for j in range(len(variables[names[1]])):
                for item in variables:
                        formatted = format(variables[item][j]) 
                        output_file.write('%15s' % str(formatted))
                output_file.write('\n')
                        
with open('rectangle_data_t.gnu','w') as rect_file:
        rect_file.write('set style rect fc lt -1 fs solid 0.1 border 2\n')
        rect_file.write('set obj rect from '+str(t_sync_a0)+' , graph 0 to graph 1, graph 1\n')
with open('rectangle_data_a.gnu','w') as rect_file:
        rect_file.write('set style rect fc lt -1 fs solid 0.1 border 2\n')
        rect_file.write('set obj rect from '+str(a0)+' , graph 0 to graph 1, graph 1\n')


