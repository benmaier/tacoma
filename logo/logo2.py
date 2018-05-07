"""
============================
Draw flat objects in 3D plot
============================

Demonstrate using pathpatch_2d_to_3d to 'draw' shapes and text on a 3D plot.
"""

import matplotlib as mpl
import matplotlib.pyplot as pl
from matplotlib.patches import Circle, PathPatch, Rectangle
# register Axes3D class with matplotlib by importing Axes3D
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D

import numpy as np

mpl.rcParams['font.family'] = 'monospace'
#mpl.rcParams['text.latex.preamble'] = r'\usepackage{pslatex}'


def text3d(ax, xyz, s, zdir="z", size=None, angle=0, usetex=False, **kwargs):

    x, y, z = xyz
    if zdir == "y":
        xy1, z1 = (x, z), y
    elif zdir == "y":
        xy1, z1 = (y, z), x
    else:
        xy1, z1 = (x, y), z

    text_path = TextPath((0, 0), s, size=size, usetex=usetex)
    trans = Affine2D().rotate(angle).translate(xy1[0], xy1[1])

    p1 = PathPatch(trans.transform_path(text_path), **kwargs)
    ax.add_patch(p1)
    art3d.pathpatch_2d_to_3d(p1, z=z1, zdir=zdir)


fig = pl.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_proj_type('ortho')

text = "TACOMA"
x = 2*np.arange(10)

for i in range(len(text)):

    text3d(ax, (3+i, 1, x[i]),
           text[i],
           zdir="x", size=2, usetex=False,
           ec="None", fc="k")

    #if i == 0 or i == 5:
    if True:
        text3d(ax, (1.5+i, 0.7, x[i]),            
               '$t_{%d}$' % (i),
               zdir="x", size=1.1, usetex=False,
               ec="None", fc="k")    
    if i != 0 and i != 5:
        ax.plot([x[i],x[i]],[2.4+i,4.1+i],[0.7,0.7],':',c=[0.7]*3,lw=1)

p = Rectangle((2.7,0.7),6,8,capstyle='round',fill=False)
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=0, zdir="x")

p = Rectangle((2.7,0.7),6,8,capstyle='round',fill=False)
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=i*2, zdir="x")

ax.plot([0,0],[6,7],[6,7],'-k')
ax.plot([10,10], [5.8,6.9],[4.5,3.7],'-k')

ax.plot([10,10],[6,7.3],[6,5.4],'-k')
ax.plot([10,10],[7,7.3],[7,5.4],'-k')
ax.plot([10,10],[6,7],[6,7],'-k')


p = Circle((7,7),0.2,color='w',ec='k')
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=0, zdir="x")

p = Circle((6,6),0.2,color='w',ec='k')
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=0, zdir="x")

p = Circle((7.3,5.4),0.2,color='w',ec='k')
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=0, zdir="x")

p = Circle((5.8, 4.5),0.2,color='w',ec='k')
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=0, zdir="x")

p = Circle((6.9,3.7),0.2,color='w',ec='k')
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=0, zdir="x")

" ================= "

p = Circle((7,7),0.2,color='w',ec='k')
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=10, zdir="x")

p = Circle((6,6),0.2,color='w',ec='k')
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=10, zdir="x")

p = Circle((7.3,5.4),0.2,color='w',ec='k')
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=10, zdir="x")

p = Circle((5.8, 4.5),0.2,color='w',ec='k')
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=10, zdir="x")

p = Circle((6.9,3.7),0.2,color='w',ec='k')
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=10, zdir="x")

#ax.plot([0,0], [5.8,6.9],[4.5,3.7],'-k')

ax.axis('off')

ax.set_xlim3d(0, 10)
ax.set_ylim3d(0, 10)
ax.set_zlim3d(0, 10)

ax.view_init(35,-47)

fig.tight_layout()
fig.savefig('logo.pdf')
fig.savefig('logo.png')

pl.show()
