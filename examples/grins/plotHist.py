import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plot
from numpy import loadtxt                 
import numpy as np
import glob
import fileinput
import pylab as P

# rc('text',usetex=True)
# font={'family' : 'normal',
#     'weight' : 'normal',
#     'size' : 12}
# matplotlib.rc('font',**font)

# cycle marker
markers = []
markers.append('^')
markers.append('s')
markers.append('<')
markers.append('p')
markers.append('>')
markers.append('o')
markers.append('v')
markers.append('D')
markers.append('h')

#############################################################
# adjoint solutions
# plot.figure()
#############################################################
plot.xlabel(r"$x$")
plot.title("Histogram of Q")
dataFile='samples.dat'
data = loadtxt(dataFile,comments="%")

n,bins,patches = P.hist(data[:],50,normed=1,histtype='stepfilled')
P.setp(patches,'facecolor','g','alpha',0.75)

# P.ylim((0.0,1.0))



P.gcf().subplots_adjust(bottom=0.15)
P.gcf().subplots_adjust(right=0.75)

P.legend(loc='upper right')
# P.grid(True)

P.savefig('histOfQ.pdf')
P.show()
#############################################################


# #############################################################
# # convergence plot
# plot.figure()
# #############################################################
# plot.title("errors ")
# plot.xlabel("n\_active\_dofs()")
# plot.ylabel("QoI error")

# uniform_file ='uniform_error.txt'
# adapt_file ='adaptive_error.txt'
# adapt_nonlinear_file ='adaptive_error_nonlinear.txt'

# uniformData = loadtxt(uniform_file)
# adaptData = loadtxt(adapt_file)
# adaptNonlinearData = loadtxt(adapt_nonlinear_file)

# # plot.loglog(uniformData[:,0],uniformData[:,1], 'k--s', ms=10, label='uniform-exact')
# plot.loglog(uniformData[:,0],uniformData[:,7], 'k-o', ms=10, label='uniform-estimate')

# # plot.loglog(adaptData[:,0],adaptData[:,1], 'b--^', ms=10, label='adapt-exact')
# plot.loglog(adaptData[:,0],adaptData[:,7], 'b-s', ms=10, label='adapt-estimate')

# # plot.loglog(adaptNonlinearData[:,0],adaptNonlinearData[:,1], 'r--o', ms=10, label='adapt-nonlinear-exact')
# plot.loglog(adaptNonlinearData[:,0],adaptNonlinearData[:,7], 'r--o', ms=10, label='adapt-nonlinear-estimate')

# plot.gcf().subplots_adjust(bottom=0.15)
# plot.gcf().subplots_adjust(left=0.15)

# plot.legend(loc="upper right", numpoints=1)
# plot.grid(True)

# plot.savefig("error.pdf")
# plot.show()
# #############################################################
