#!/usr/bin/gnuplot

set samples 10000

###############################################################################
### First plot
###############################################################################

set term pngcairo size 400, 400
set output "images/sin_first_plot.png"
plot sin(x)
set output
set terminal x11


###############################################################################
### Ranges
###############################################################################

set term pngcairo size 400, 400
set output "images/sin_ranges_1.png"
plot[0:2] sin(x)
set output
set terminal x11

# -----------------------------------------------------------------------------

set term pngcairo size 400, 400
set output "images/sin_ranges_2.png"
plot[0:2 * pi][-2:3] sin(x)
set output
set terminal x11

# -----------------------------------------------------------------------------

set term pngcairo size 400, 400
set output "images/sin_ranges_3.png"
plot[0:] sin(x)
set output
set terminal x11


###############################################################################
### Functions
###############################################################################

set term pngcairo size 500, 300
set output "images/functions_1.png"
plot[0:4] cos(2 * x) * cos(50 * x) * exp(-x), exp(-x), -exp(-x)
set output
set terminal x11

# -----------------------------------------------------------------------------

set term pngcairo size 400, 400
set output "images/functions_2.png"
foo(x) = (x - 2) * (x - 1) * x * (x + 1) * (x + 2)
bar(x) = 1 / (x + 0.5) - 1
plot[-3:3][-3:3] foo(x), bar(x)
set output
set terminal x11


###############################################################################
### Title
###############################################################################

set term pngcairo size 400, 400
set output "images/titles_1.png"
plot[0:pi][-1:2] sin(x) title "Hello world!", 2 * x title "", \
                 1 / sqrt(x) notitle
set output
set terminal x11

# -----------------------------------------------------------------------------

set term pngcairo size 400, 400
set output "images/titles_2.png"
set key right bottom
plot x + sin(x) title "set key right bottom"
set key right top
set output
set terminal x11


###############################################################################
### Data files
###############################################################################

set term pngcairo enhanced size 400, 400
set output "images/data_files_1.png"
! python3 -c \
'for x in range(10): \
    print("{}\t{}".format(x, 2**-x)) \
' > data/data_files_1.dat
set grid
plot[-1:10][-0.2:1.2] "data/data_files_1.dat" title "Calculated 2^-^x"
set output
set terminal x11

# -----------------------------------------------------------------------------

set term pngcairo enhanced size 400, 400
set output "images/data_files_2.png"
! python3 -c \
'for x in range(1,10): \
    print("{:.2e}\t{:.2e}\t{:.2e}".format(x, 2**-x, 1./x)) \
' > data/data_files_2.dat
set grid
plot[-1:10][-0.2:1.2] "data/data_files_2.dat" using 1:2 title "2^-^x", \
                      "" u 1:3 title "1/x", \
                      "" u 1:(1. / ($1 - 2)) tit "1/(x-2)"
set output
set terminal x11

# -----------------------------------------------------------------------------

set term pngcairo enhanced size 400, 400
set output "images/data_files_3.png"
set key right bottom
plot "data/data_files_2.dat" using 0:($0 + 1) title "x"
set key right top
set output
set terminal x11


###############################################################################
### Style
###############################################################################

set term pngcairo enhanced size 400, 400
set output "images/styles_1.png"
! python3 -c \
'from math import *; \
 list(map(lambda x: print("{}\t{}".format(x, sin(0.5 * x))), range(20))) \
' > data/styles_1.dat
set grid
plot[-1:20][-8:8] "data/styles_1.dat" with boxes title "boxes", \
                  "" u 1:(-$2) with impulses title "impulses", \
                  "" u 1:($2 - 2) with lines title "lines", \
                  "" u 1:($2 - 4) with linespoints title "linespoints", \
                  "" u 1:($2 - 6) with points title "boxes"
set output
set terminal x11

# -----------------------------------------------------------------------------

set term pngcairo enhanced size 800, 600
set output "images/test.png"
# lt ---> linetype
# lc ---> linecolor
# lw ---> linewidth
# ps ---> pointsize
test
set output
set terminal x11

# -----------------------------------------------------------------------------

set term pngcairo enhanced size 400, 400
set output "images/styles_2.png"
! python3 -c \
'import numpy as np; \
 mu, sigma = 0, 4; \
 N = 20; \
 errs = np.random.normal(mu, sigma, N); \
 t = np.linspace(10, 30, N); \
 x = 3 * t + 4 + errs; \
 data = np.dstack((t,x,errs))[0]; \
 np.savetxt("data/styles_2.dat", data, delimiter="\t") \
'
set grid
set key right bottom
# p  ---> points
# lt ---> linetype
# lc ---> linecolor
# ps ---> pointsize
set xlabel "Time, t [s]"
set ylabel "Particle position, X [cm]"
plot[8:32] "data/styles_2.dat" with p lt 7 ps 1 lc rgb "forest-green" \
                               title "Measured data"
set output
set terminal x11

# -----------------------------------------------------------------------------

set term pngcairo enhanced size 400, 400
set output "images/styles_3.png"
plot[8:32] "data/styles_2.dat" u 1:2:3 w yerrorbars \
                               lt 7 ps 1 lc rgb "forest-green" \
                               title "Measured data"
set output
set terminal x11


###############################################################################
### Useful Stuff
###############################################################################

set term pngcairo enhanced size 400, 400
set output "images/useful_1.png"
position(t) = V * t + X0
set fit quiet
set fit errorvariables
fit position(x) "data/styles_2.dat" u 1:2:3 via V, X0
plot[8:32] "data/styles_2.dat" u 1:2:3 w yerrorbars \
                               lt 7 ps 1 lc rgb "forest-green" \
                               title "Measured data", \
           position(x) lt 1 lw 2 lc rgb "dark-blue" \
                       title "Approximation"
print "V = ", V, " +/- ", V_err
set output
set terminal x11

# -----------------------------------------------------------------------------

set term pngcairo enhanced size 400, 400
set output "images/useful_2.png"
! python3 -c \
'import numpy as np; \
 \
 mu, sigma = 0, 4; \
 N = 400; \
 errs = np.random.normal(mu, sigma, N); \
 t = np.linspace(10, 30, N); \
 x = 3 * t + 4 + errs; \
 data = np.dstack((t,x,errs))[0]; \
 np.savetxt("data/useful_2.dat", data, delimiter="\t") \
'
position(t) = V * t + X0
set fit quiet
set fit errorvariables
fit position(x) "data/useful_2.dat" u 1:2:3 via V, X0
plot[8:32] "data/useful_2.dat" u 1:2 w p \
                               lt 7 ps 0.4 lc rgb "forest-green" \
                               title "Measured data", \
           position(x) lt 1 lw 2 lc rgb "dark-blue" \
                       title "Approximation"
set output
set terminal x11

# -----------------------------------------------------------------------------

set term pngcairo enhanced size 400, 400
set output "images/useful_3.png"
set key right top
set xlabel "Time, t [s]"
set ylabel "Deviation, dX [cm]"
deviation(t, x) = x - position(t)
plot[8:32][-15:15] "data/useful_2.dat" u 1:(deviation($1, $2)) w impulses \
                                       lt 7 lw 0.4 lc rgb "forest-green" \
                                       title "Deviation from theory"
set output
set terminal x11

# -----------------------------------------------------------------------------

set term pngcairo enhanced size 400, 400
set output "images/useful_4.png"
# This genius way of making histograms in Gnuplot is
# taken from http://stackoverflow.com/a/2538846.
#
set xlabel "Deviation, dX [cm]"
set ylabel "Number of points, N [1]"
# We want our boxes to be filled. Here, we make them 50%
# transparent with solid black border. For more, see
# `help set style fill`
set style fill transparent solid 0.5 border rgb "black" 
# The following function rounds `x` down to the 
# nearest multiple of `width`
bin(x, width) = width * floor(x / width)
# Get some statistics about the second column of our data
stats "data/useful_2.dat" u 2 name "S" nooutput
binwidth = 5 * (S_max - S_min) / S_records 
# Magic...
plot[-15:15] "data/useful_2.dat" u (bin(deviation($1,$2), binwidth)):(1.0) \
                                 smooth frequency w boxes \
                                 lt 1 lc rgb "forest-green" notitle
set output
set terminal x11

# -----------------------------------------------------------------------------

set term pngcairo enhanced size 400, 400
set output "images/useful_5.png"
gauss(x) = A * exp(- (x - mu)**2 / (2 * sigma**2))
A = 50.0
mu = 0.1
sigma = 5.0
set table ".temp.dat"
plot "data/useful_2.dat" u (bin(deviation($1,$2), binwidth)):(1.0) \
                         smooth frequency
unset table
fit gauss(x) ".temp.dat" via A, sigma, mu
! rm .temp.dat
plot[-15:15][0:50] "data/useful_2.dat" u (bin(deviation($1,$2), binwidth)):(1.0) \
                                       smooth frequency w boxes \
                                       lt 1 lc rgb "forest-green" \
                                       title "Measured distribution", \
                   gauss(x) lt 1 lw 2 lc rgb "dark-blue" \
                            title "Gaussian distribution"
print "A = ", A, " +/- ", A_err
print "mu = ", mu, " +/- ", mu_err
print "sigma = ", sigma, " +/- ", sigma_err
set output
set terminal x11


###############################################################################
### Difficult Stuff
###############################################################################

set term epslatex size 10cm,10cm
set output "my-plot.tex"
set xlabel "Deviation, $\\Delta x$ [cm]"
set ylabel "Number of points, $N$ [1]"
plot[-15:15][0:55] \
    "data/useful_2.dat" u (bin(deviation($1,$2), binwidth)):(1.0) \
                        smooth frequency w boxes \
                        lt 1 lc rgb "forest-green" \
                        title "Measured", \
    gauss(x) lt 1 lw 2 lc rgb "dark-blue" \
             title "$A\\cdot\\exp\\left(\\frac{(x-\\mu)^2}{2\\sigma^2}\\right)$ fit"
set output
set terminal x11
! pdflatex document.tex && pdfcrop document.pdf
! gs -sDEVICE=png16m \
     -dTextAlphaBits=4 \
     -r150 \
     -o "images/difficult_1.png" document-crop.pdf
! rm -f document.aux document.log document.pdf \
        document-crop.pdf \
	my-plot.tex my-plot.eps my-plot-eps-converted-to.pdf

# -----------------------------------------------------------------------------

set term epslatex size 10cm,10cm
set output "my-plot.tex"
set xlabel "Deviation, $\\Delta x$ [cm]"
set ylabel "Number of points, $N$ [1]"
set key width -7.5 height 0.5 spacing 1.5 box # Adjust the spacing between 
                                              # lines and create a frame
                                              # around the titles.
set object 1 rectangle from 4,37.5 to 14.1,44 front
set label 1   "$\\mu = " . sprintf("%.1f", mu) \
            . "\\pm " . sprintf("%.1f", ceil(10 * mu_err) / 10.) . "$" \
	    at 4.2,42 front
set label 2   "$\\sigma = " . sprintf("%.1f", sigma) \
            . "\\pm " . sprintf("%.1f", ceil(10 * sigma_err) / 10.) . "$" \
	    at 4.2,39 front
plot[-15:15][0:55] \
    "data/useful_2.dat" u (bin(deviation($1,$2), binwidth)):(1.0) \
                        smooth frequency w boxes \
                        lt 1 lc rgb "forest-green" \
                        title "Measured", \
    gauss(x) lt 1 lw 2 lc rgb "dark-blue" \
             title "$A\\cdot\\exp\\left(\\frac{(x-\\mu)^2}{2\\sigma^2}\\right)$ fit"
unset object 1
unset label 1
unset label 2
set output
set terminal x11
! pdflatex document.tex && pdfcrop document.pdf
! gs -sDEVICE=png16m \
     -dTextAlphaBits=4 \
     -r150 \
     -o "images/difficult_2.png" document-crop.pdf
! rm -f document.aux document.log document.pdf \
        document-crop.pdf \
        my-plot.tex my-plot.eps my-plot-eps-converted-to.pdf



###############################################################################
### Zoom
###############################################################################

exit 0

set multiplot layout 1, 2
set xlabel "Time, t [s]"
set ylabel "Particle position, X [cm]"
set key right bottom
plot "datafile.dat" u 1:2 w p \
     lt 7 ps 0.5 lc rgb "forest-green" title "Measured data", \
     position(x) lt 1 lw 2 lc rgb "dark-blue" title "Approximation"
set xlabel "Deviation from the theoretical value"
set ylabel "Number of points"
set key right top
set offsets 1, 1, graph 0.1, 0
plot "datafile.dat" u (bin($2 - position($1), binwidth)):(1.0) \
     smooth freq with boxes \
     fs solid 0.50 border lt 7 lt 1 lc rgb "forest-green" title "Histogram"
unset multiplot

set output
set terminal x11
