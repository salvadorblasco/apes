#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

# This script reads a data file containing the labels and the data of equilibria
# calculated whith HySS. What does is:
#  1. Open the file and read the data
#  2. Transform data
#  3. Calculate label position
#      The label x-position is right over the maximum. The y-position is calculated
#     for the labels not to collide.


import sys
#from scipy import optimize
import numpy
from matplotlib import pyplot
import matplotlib 
import matplotlib.colors, matplotlib.collections, matplotlib.cm
import matplotlib.text, matplotlib.transforms
import matplotlib.afm, matplotlib.axis
import re

class DragableLabel:
    def __init__(self, label):
        self.label = label
        self.press = None

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.label.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.label.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.label.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.label.axes: return

        contains, attrd = self.label.contains(event)
        if not contains: return
        #print 'event contains', self.label.get_position()
        x0, y0 = self.label.get_position()
        self.press = x0, y0, event.xdata, event.ydata

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.press is None: return
        if event.inaxes != self.label.axes: return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        #print 'x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f'%(x0, xpress, event.xdata, dx, x0+dx)
        self.label.set_x(x0+dx)
        self.label.set_y(y0+dy)

        self.label.figure.canvas.draw()

    def on_release(self, event):
        'on release we reset the press data'
        self.press = None
        self.label.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.label.figure.canvas.mpl_disconnect(self.cidpress)
        self.label.figure.canvas.mpl_disconnect(self.cidrelease)
        self.label.figure.canvas.mpl_disconnect(self.cidmotion)

class Iterate_Over_File:
    def __init__(self, f):
        self.f = f
        self.line_counter = 0

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        t = self.f.next()
        self.line_counter += 1
        
        while t[0] == '#':          # Ignore comments is the only
            t = self.f.next()       # purpose of this class
            self.line_counter += 1

        return t
        
    def close(self):
        self.f.close()

    def current_line(self):
        return self.line_counter 

def read_range(r):
    '''From a range format, returns a list of numbers.
    e.g.  
        2,5-9,11 -> [2, 5, 6, 7, 8, 9, 11] '''

    assert isinstance(r, str)
    list_of_ranges = r.split(',')
    rlist = list()

    for c in list_of_ranges:
        d = c.strip().split('-')
        if len(d) == 1:
            assert d[0].isdigit()  
            rlist.append(int(d[0]))
        elif len(d) == 2:
            assert d[0].isdigit() and d[1].isdigit()  
            rlist.extend(range(int(d[0]), int(d[1])+1))
        else:
            print "Parse error at line %d" % line_counter
            sys.exit(1)

    return rlist


# sumary of variables
#
# int  line_counter  current line number of input file
# int  special_n     number of specials
# dict stoichio      array of ints. stoichiometric coefficientes

# +---------------------------+
# |         Read data         |
# +---------------------------+

# Datafile format
# ------------------------
# lines meaning              format
# 1     heading: contains a list of %variable = value. End in blank line.
# 2     y-label              "
# 3     x data column        x=<int>
# 4     y data column        y=(<int>-<int>),...
# 5     operation on data    python code to execute over data
# 6     x-shift              label x shift correction
# 7     y-shift              label y shift correction
# 8     data labels          short space-separated strings
# 9-end data                 space-separated float values

#float = "\d\+\.(\d\+)?"

# Open file 
#f = Iterate_Over_File(open(sys.argv[1], "r"))
f = Iterate_Over_File(sys.stdin)

line_counter = 0

# Read headings
keywords=["title","ylabel", "datacolor", "datastyle", "unlabel_lower_than", "update_labels", "ignore_lower_than",
    "special_sumof","arrow"]
array_keywords=["override_datacolor","override_datastyle", "legend_xy", "special_title"]

override_datastyle = dict()
override_datacolor = dict()

re_kw  = re.compile("(" + '|'.join(keywords) + ")=(.*)")
re_kwa = re.compile("(" + '|'.join(array_keywords) + ")\[(-?\d+)\]=(.*)")

title = ""
y_title = ""

for l in f:
    line_counter += 1
    #print str(line_counter) + ": " + l[0:-1]
    if l.strip() == '':
        break

    kw_match=re_kw.match(l)
    if kw_match:
        #print kw_match.groups()
        #print kw_match.groups()[0] == 'unlabel_lower_than'
        if kw_match.groups()[0] == keywords[0]:        # title
            title=kw_match.group()[1]
            #print "found title"
        if kw_match.groups()[0] == keywords[1]:        # ylabel 
            y_title=kw_match.groups()[1]
            #print "found y_label"
        if kw_match.groups()[0] == keywords[4]:        # unlabel_lower_than
            unlabel_lower_than=float(kw_match.groups()[1])
            #print "found unlabel_lower_than"
        if kw_match.groups()[0] == keywords[6]:        # ignore_lower_than 
            ignore_lower_than=float(kw_match.groups()[1])
            #print "ignore_lower_than found"
        if kw_match.groups()[0] == keywords[7]:        # special_sumof
            if 'special_sumof' not in dir():
                special_sumof = list()
            assert isinstance(special_sumof, list)
            special_sumof.append(kw_match.groups()[1])
        if kw_match.groups()[0] == keywords[8]:        # arrow
            if 'arrows' not in dir():
                arrows = list()
            assert isinstance(arrows, list)
            arrows.append(kw_match.groups()[1])
        continue

    kwa_match=re_kwa.match(l)
    if kwa_match:
        w = kwa_match.groups()[0]
        k = kwa_match.groups()[1]
        v = kwa_match.groups()[2]

        assert isinstance(int(k), int)

        if w == array_keywords[0]:      # override_datacolor
            #print "override_datacolor not implemented yet @line %d" % line_counter
            override_datacolor[int(k)] = v

        if w == array_keywords[1]:      # override_datastyle 
            override_datastyle[int(k)] = v

        if w == array_keywords[2]:      # legend_xy
            print "legend_xy not implemented yet @line %d" % line_counter

        if w == array_keywords[3]:      # special_title
            if 'special_title' not in dir():
                special_title = dict()

            special_title[k] = v

        continue

    print "Syntax error at line " + str(line_counter) + " ('" + l.strip() + "')"
    f.close()
    sys.exit(1)

x_line = f.next().strip()
line_counter += 1
re_x = re.compile(r"x=(\d+)(\s+from\s+(?P<from>\d+\.\d+)\s+to\s(?P<to>\d+\.\d+))?")
#print r"x=(\d+)\s+from\s+(?P<from>\d+\.\d+)\s+to\s(?P<to>\d+\.\d+)"
#print x_line
x_match = re_x.match(x_line)

if not x_match:
    print "Syntax error at line " + str(line_counter) + ": '" + x_line + "'"
    f.close()
    sys.exit(1)

x_col = int(x_match.group(1)) - 1
#print "x_col is " + str(x_col)

y_toplot = re.match(r"y=(.+)", f.next().strip()).group(1)
transform_operator = f.next()
line_counter += 2

stoichio = list()
for l in f:
    line_counter += 1
    if l.strip() == '':
        break
    else:
        stoichio.append(l.split())

#print stoichio

# Read eigth line -> species names
species_names = f.next().split()
line_counter += 1
x_axis_title = species_names[x_col]

# Read other lines -> data matrix
raw_list = list()
for l in f:
    if l.strip() == '':
        break
    raw_list.append([float(i) for i in l.split()])


#raw_data = numpy.fromstring(txt, dtype=float,sep='\t') 
raw_data = numpy.array(raw_list)
#print raw_data
f.close()


# +---------------------------+
# |     Transform data        |
# +---------------------------+

# Get columns to plot
y_ranges = y_toplot.split(',')
y_cols = list()
for c in y_ranges:
    d = c.split('-')
    if len(d) == 1:
        y_cols.append(int(d[0]))
    elif len(d) == 2:
        y_cols.extend(range(int(d[0]), int(d[1])+1))
    else:
        print "Parse error in ranges"
        sys.exit(1)

special_n = 0
if 'special_sumof' in dir():
    special_n += len(special_sumof)      # TODO add new specials if any

print raw_data.shape
x_data = raw_data[:,x_col]
y_data = numpy.zeros((x_data.shape[0], len(y_cols) + special_n))

if transform_operator.strip() != 'y':
    re_n = re.compile("n\[(\d+)\]").search(transform_operator)
    re_S = re.compile("\$\[(\d+)\]").search(transform_operator)
    if not re_S:
        print "Must specify data column as $[number]"
        sys.exit(1)

    for i, c in enumerate(y_cols):
        rule = transform_operator.replace('$[' + re_S.groups(1)[0] + ']', 'raw_data[:,' + re_S.groups(1)[0] + '-1]' )
        if re_n:
            rule = rule.replace('n['+ re_n.groups(1)[0] + ']', 'int(stoichio['+ re_n.groups(1)[0] + '][' + str(c-1) + '])')
        rule = rule.replace('y', 'raw_data[:,' + str(c-1) + ']')
        rule = 'y_data[:,' + str(i) + '] = ' + rule
        #print rule.strip()
        #exec("y_data[:,0] = raw_data[:,5]/raw_data[:,3-1]*(stoichio[1][5])*100")
        exec rule
else:
    y_data=raw_data[:,:]

# add specials

if 'special_sumof' in dir():
    for n, s in enumerate(special_sumof):
        r = read_range(s)
        t = [y_cols.index(i) for i in r]
        tmp = numpy.sum(y_data[:,t],1)
        y_data[:,len(y_cols)+n] = tmp

labels_text = [ species_names[i-1] for i in y_cols]
y_cols.extend(range(-1,-special_n-1,-1))  # y_cols is extended with specials (-1, -2, ...)

if  x_match.group('from') is not None:
    x_min = float(x_match.group('from'))
else:
    x_min = numpy.min(x_data)   

if  x_match.group('to') is not None:
    x_max = float(x_match.group('to'))
else:
    x_max = numpy.max(x_data)   

print "plotting from %f to %f" % (x_min, x_max)

# Get maxima -> fit x
labels_y = numpy.max(y_data, 0)
for i, l in enumerate(labels_y):
    if numpy.isnan(l) or numpy.isinf(l):
        labels_y[i] = 0.0

labels_x = list()
plot_size = len(y_cols)
for r, m in enumerate(labels_y):
    ff = numpy.where(y_data[:,r] == m)
    if len(ff[0]) == 0:
        a = x_min
    else:
        a = ff[0][0]
    labels_x.append(
        x_data[a]
    )

#print labels_x
#print labels_y

# specials names
for i in range(1, special_n + 1):
    if str(-i) not in special_title.keys():
        print "warning: special -%d has no title" % i
        t = None
    else:
        t = special_title[str(-i)]

    labels_text.append(t)

#x_label_shift = [ float(full_x_label_shift[i-1]) for i in y_cols ] 
#y_label_shift = [ float(full_y_label_shift[i-1]) for i in y_cols ] 
x_label_shift = [ 0 for i in y_cols ] 
y_label_shift = [ 0 for i in y_cols ] 

#if 'unlabel_lower_than' in dir():
absolute_maxima = numpy.amax(y_data,0)
#    print absolute_maxima 


# initialise matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', size=24.0)
matplotlib.rc('axes', linewidth=2)
matplotlib.rc('figure.subplot', bottom=0.12)
matplotlib.rc('figure.subplot', right=0.95)
matplotlib.rc('figure.subplot', top=0.95)

#colormap=matplotlib.cm.spectral(i*256/plot_size),
#colormap=matplotlib.cm.spectral(i*256/plot_size)

default_linestyle="-"
default_datacolor="b"

#pyplot.bone()

#pylab.title(title)

labels = list()
drlbls = list()

#plot_size = len(y_cols)

#print override_datastyle

for i, ny in enumerate(y_cols):

    #if ('ignore_lower_than' in dir()):
    #    print "ignoring lower than ", ignore_lower_than

    if ('ignore_lower_than' in dir()) and \
          (absolute_maxima[i] < ignore_lower_than):
        continue

    if ny in override_datastyle.keys():
        current_datastyle = override_datastyle[ny]
    else:
        current_datastyle = default_datacolor + default_linestyle

    #if ny in override_datastyle.keys():
    #    current_linestyle = override_datastyle[ny]
    #    #print "current_linestyle = " + override_datastyle[ny]
    #else:
    #    current_linestyle = default_linestyle

    #if ny in override_datacolor.keys():
    #    current_color = override_datacolor[ny]
    #else:
    #    current_color = default_datacolor

    # plot curve
    pyplot.plot(
        x_data, 
        y_data[:,i], 
        current_datastyle,
        linewidth=2,
        solid_joinstyle='bevel')

    # fit label position
    if i == 0: 
        align = 'right'
    elif i == plot_size-1:
        align = 'left'
    else:
        align = 'center'

    # put label

    if (labels_text[i] is not None) and ((not ('unlabel_lower_than' in dir())) or \
          (absolute_maxima[i] >= unlabel_lower_than)):

        labels.append(
            pyplot.text(
                labels_x[i] + x_label_shift[i], 
                labels_y[i] + y_label_shift[i], 
                labels_text[i], #r'$\mathrm{' + labels_text[i] + r'}$', 
                verticalalignment='bottom',
                horizontalalignment=align,
                color='k',
                fontname="Times"
            )
        )

        drlbl = DragableLabel(labels[-1])
        drlbl.connect()
        drlbls.append(drlbl)

        x_aux, y_aux = labels[-1].get_position()

# draw arrows
if 'arrows' in dir():
    for a in arrows:
        x1,y1,x2,y2 = map(float,a.split(','))
        pyplot.annotate('',xy=(x2,y2),xytext=(x1,y1),xycoords='data',textcoords='data',color='black', arrowprops=dict(arrowstyle="-|>",lw=2.0,color='black'))

pyplot.xlim(x_min, x_max)
pyplot.ylim(0,100)
#print "setting limits from " + str(x_min) + " and " + str(x_max)

pyplot.xlabel(x_axis_title)
pyplot.ylabel(y_title)

pyplot.show()




