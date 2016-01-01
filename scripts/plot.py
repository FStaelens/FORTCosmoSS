#!/usr/bin/env python

import numpy
import matplotlib.pyplot as plt
import argparse

# arguments handling

parser = argparse.ArgumentParser(
    prog='PLOT',
    description='plot [y] vs x from file',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('data_file', metavar='data_file', type=str,
    help='the data file')
parser.add_argument('x_array', metavar='x', type=str,
    help='the name of the x-axis array')
parser.add_argument('y_arrays', metavar='y', type=str, nargs='+',
    help='the set of arrays to plot vs x')
parser.add_argument('--scale', metavar='SCALE', type=str,
    help='choose semilogx, semilogy, loglog, default is linear')
parser.add_argument('--x_range', metavar=('X_MIN','X_MAX'), type=float , nargs=2,
    help='the x-range')
parser.add_argument('--y_range', metavar=('Y_MIN','Y_MAX'), type=float , nargs=2,
    help='the y-range')
args = parser.parse_args()

fig = plt.figure()
ax = plt.axes()

ax.set_xlabel(args.x_array)

if args.x_range != None:
    ax.set_xlim(args.x_range[0],args.x_range[1])

if args.y_range != None:
    ax.set_ylim(args.y_range[0],args.y_range[1])

with open(args.data_file,'r') as data_file:
    time = data_file.readline()
    header = data_file.readline()
    names = header.split()
    data = numpy.loadtxt(args.data_file, skiprows=2)

variables = {}
lines = {}
for j in range(len(names)):
    variables[names[j]] = data[:,j]

for arg in args.y_arrays:
    if arg in list(variables.keys()):
	if args.scale == None:
            lines[arg], = ax.plot([], [], label=arg, lw=2, ls='-')
        elif args.scale == 'semilogx':
            lines[arg], = ax.semilogx([], [], label=arg, lw=2, ls='-')
        elif args.scale == 'semilogy':
            lines[arg], = ax.semilogy([], [], label=arg, lw=2, ls='-')
        elif args.scale == 'loglog':
            lines[arg], = ax.loglog([], [], label=arg, lw=2, ls='-')
    else:
        print 'ignoring argument:', arg,'Check presence in', args.data_file
    
    x = variables[args.x_array]
    for key in list(lines.keys()):
        y = variables[key]
        lines[key].set_data(x,y)

    ax.relim()
    ax.autoscale_view(True,True,True)

plt.title(time)
plt.grid()
plt.legend()
plt.show()
