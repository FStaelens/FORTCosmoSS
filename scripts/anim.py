#!/usr/bin/env python

import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import argparse
import os, os.path

# arguments handling

parser = argparse.ArgumentParser(
	prog='ANIM',
	description='animate [y] vs x from ./data/',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('y_arrays', metavar='y', type=str, nargs='+',
                    help='the set of arrays to plot vs x')
parser.add_argument('--x_range', metavar=('X_MIN','X_MAX'), type=long , nargs=2,
                    help='the x-range')
parser.add_argument('--y_range', metavar=('Y_MIN','Y_MAX'), type=float , nargs=2,
    		    help='the y-range')
args = parser.parse_args()

DIR='./data'

fig = plt.figure()
ax = plt.axes()

if args.x_range != None:
    ax.set_xlim(args.x_range[0],args.x_range[1])

if args.y_range != None:
    ax.set_ylim(args.y_range[0],args.y_range[1])

with open(DIR+'/data0000.dat','r') as data_file:
    data_file.readline()
    header = data_file.readline()
    names = header.split()
    data = numpy.loadtxt(DIR+'/data0000.dat', skiprows=2)

variables = {}
lines = {}
for j in range(len(names)):
    variables[names[j]] = data[:,j]

for arg in args.y_arrays:
    if arg in list(variables.keys()):
        lines[arg], = ax.plot([], [], label=arg, lw=1, ls='-', marker='+')


def init():
    with open(DIR+'/data0000.dat','r') as data_file:
        time = data_file.readline()
        header = data_file.readline()
    
    time = time.replace('#','')
    names = header.split()
    data = numpy.loadtxt(DIR+'/data0000.dat', skiprows=2)

    variables = {}
    for j in range(len(names)):
        variables[names[j]] = data[:,j]

    plt.title(time)
    
    x = variables['x']
    for key in list(lines.keys()):
        y = variables[key]
        lines[key].set_data(x,y)

    ax.relim()
    ax.autoscale_view(True,True,True)

    return lines,

def animate(i):
    string_i = str(i).zfill(4)
    with open(DIR+'/data'+string_i+'.dat','r') as data_file:
        time = data_file.readline()
        header = data_file.readline()

    time = time.replace('#','')
    names = header.split()
    data = numpy.loadtxt(DIR+'/data'+string_i+'.dat', skiprows=2)

    variables = {}
    for j in range(len(names)):
        variables[names[j]] = data[:,j]
 
    plt.title(time)
    plt.legend()

    x = variables['x']
    for key in list(lines.keys()):
        y = variables[key]
        lines[key].set_data(x,y)

    ax.relim()
    ax.autoscale_view(True,True,True)
    ax.grid(True)

    return lines,

data_number = len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])-1    # number of data files
anim = animation.FuncAnimation(fig,animate,init_func=init,frames=data_number,interval=20,blit=False)

anim.save('movie.mp4',writer='ffmpeg',dpi=200,bitrate=-1)
print ' --> output produced in "./movie.mp4"'

#plt.show()




