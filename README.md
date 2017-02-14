# Intro to Gnuplot

![Gnuplot Homepage Screenshot](images/gnuplot_homepage.png)
From the homepage, `Gnuplot` is a
> portable command-line driven graphing utility for Linux, OS/2, MS Windows, 
> OSX, VMS, and many other platforms.


* Download (http://gnuplot.info/download.html)
* Install
* Start!


This is how `Gnuplot` start screen looks on my machine:
```
tom@wheelhorse:~$ gnuplot

	G N U P L O T
	Version 4.6 patchlevel 6    last modified September 2014
	Build System: Linux x86_64

	Copyright (C) 1986-1993, 1998, 2004, 2007-2014
	Thomas Williams, Colin Kelley and many others

	gnuplot home:     http://www.gnuplot.info
	faq, bugs, etc:   type "help FAQ"
	immediate help:   type "help"  (plot window: hit 'h')

Terminal type set to 'x11'
gnuplot> 
```

It's just like `Python` shell: type in your command, press `ENTER`, and
your command gets executed!

First of all, here are some useful commands:
* `exit` or `quit` to exit the application;
* `save <filename>` to save all your current setup to file;
* `load <filename>` -- the inverse -- to load settings/commands from file;
* `help <topic>` to get help!



### Simple plotting
From `help plot`:
> `plot` is the primary command for drawing plots with `gnuplot`.  It creates
> plots of  functions and data in many, many ways.  `plot` is used to draw 2D
> functions and data;

#### Syntax
```
plot  [<ranges>]
      [<iteration>]
      <function> | <datafile> [datafile-modifiers]
      [axes <axes>] [<title-spec>] [with <style>]
      [, ...]
```

-------------------------------------------------------------------------------

```gp
gnuplot> plot sin(x)
```
![](images/sin_first_plot.png)



##### Ranges
From `help ranges`:
> The optional axis ranges at the start of a plot command specify the region of
> the graph that will be displayed.  These override any ranges established by a
> previous `set range` statement.

```gp
gnuplot> plot[0:2] sin(x)
```
![](images/sin_ranges_1.png)

-------------------------------------------------------------------------------

```gp
plot[0:2 * pi][-2:3] sin(x)
```
![](images/sin_ranges_2.png)

-------------------------------------------------------------------------------

```gp
gnuplot> plot[0:] sin(x)
```
![](images/sin_ranges_3.png)



##### Function
From `help function`:
> Built-in or user-defined functions can be displayed by the `plot` and `splot`
> commands in addition to, or instead of, data read from a file. 

```gp
gnuplot> plot[0:4] cos(2 * x) * cos(50 * x) * exp(-x), exp(-x), -exp(-x)
```
![](images/functions_1.png)

-------------------------------------------------------------------------------

```gp
gnuplot> foo(x) = (x - 2) * (x - 1) * x * (x + 1) * (x + 2)
gnuplot> bar(x) = 1 / (x + 0.5) - 1
gnuplot> plot[-3:3][-3:3] foo(x), bar(x)
```
![](images/functions_2.png)


List of all `Gnuplot` built-in functions: 
http://gnuplot.sourceforge.net/docs_4.2/node53.html


##### Title
From `help plot title`:
> By default each plot is listed in the key by the corresponding function or 
> file name. You can give an explicit plot title instead using the `title` 
> option.

```gp
gnuplot> plot[0:pi][-1:2] sin(x) title "Hello world!", \
                          2 * x title "", \
                          1 / sqrt(x) notitle
```
![](images/titles_1.png)

Notice the `pi` in xrange -- `Gnuplot` has it predefined.
```gp
gnuplot> print pi
3.14159265358979
```
To see all defined variables use `show variables all`; `show variables pi` 
prints variables starting with `pi`.

Sometimes it is desired to change the key position. This can be done with the
`set key` command. For example,  `set key left center` places the  key on the
left horizontally and in the center vertically; `set key center bottom` --
in the center horizontally and at the bottom vertically. The default is
`set key right top`. For more, see `help set key`.

```gp
gnuplot> set key right bottom
gnuplot> plot x + sin(x) title "set key right bottom"
```
![](images/titles_2.png)


From `help plot title`:
> There is also an option that will interpret the first entry in a column of
> input data (i.e. the column header) as a text field, and use it as the key
> title.  See `datastrings`.




##### Data files
From `help datafile`:
> Discrete data contained in a file can be displayed by specifying the name of
> the data file (enclosed in single or double quotes) on the `plot` command line.

`Gnuplot` supports both text and binary files. Text files must consist of one or
more columns separated by `TAB`s. Just like in `gp` a hashtag `#` starts a
comment. It is common practice to give files of such format the `.dat` extension.

```gp
# data_files_1.dat
#
###############################################
# DO use comments! Otherwise you end up       #
# having multitude of files containing 'some' #
# useful data and no way of figuring out what #
# this data means.                            #
###############################################
#
# X	2**(-X)

0	1
1	0.5
2	0.25
3	0.125
4	0.0625
5	0.03125
6	0.015625
7	0.0078125
8	0.00390625
9	0.001953125
```

```gp
! python3 -c \
'for x in range(10): \
    print("{}\t{}".format(x, 2**-x)) \
' > data/data_files_1.dat
set grid
plot[-1:10][-0.2:1.2] "data/data_files_1.dat" title "Calculated 2^-^x"
```
![](images/data_files_1.png)

-------------------------------------------------------------------------------

```gp
! python3 -c \
'for x in range(1,10): \
    print("{:.2e}\t{:.2e}\t{:.2e}".format(x, 2**-x, 1./x)) \
' > data/data_files_2.dat
set grid
plot[-1:10][-0.2:1.2] "data/data_files_2.dat" using 1:2 title "2^-^x", \
                      "" u 1:3 title "1/x", \
                      "" u 1:(1. / ($1 - 2)) tit "1/(x-2)"
```
![](images/data_files_2.png)

Use `<N1>:<N2>` to specify columns you want to plot. Here `<N1>` is the "x-column"
and `<N2>` is the "y-column". You can also perform mathematical operations on
columns -- use `$<N>` to get the numerical value of column `N`.

Column `0` gives you the number of the point.
```gp
set key right bottom
plot "data/data_files_2.dat" using 0:($0 + 1) title "x"
```
![](images/data_files_3.png)


##### Style
```gp
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
```
![](images/styles_1.png)

-------------------------------------------------------------------------------

To test which "style"s are available, type `test`:
![](images/test.png)

-------------------------------------------------------------------------------


Let's emulate a real world example. Suppose, we are measuring position of a 
particle as a function of time. We thus have three columns in our file:
* time point `t` in seconds, 
* `X`-coordinate in cm, and
* error in `X`.
```gp
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
```
![](images/styles_2.png)


But we need errorbars!

From `help errorbars`:
> Error bars are supported for 2D data file plots by reading one to four
> additional columns (or `using` entries); these additional values are used in
> different ways by the various errorbar styles.

Possible data layouts:
```
	(x, y, ydelta)
	(x, y, ylow, yhigh)
	(x, y, xdelta)
	(x, y, xlow, xhigh)
	(x, y, xdelta, ydelta)
	(x, y, xlow, xhigh, ylow, yhigh)
```

```gp
plot[8:32] "data/styles_2.dat" u 1:2:3 w yerrorbars \
                               lt 7 ps 1 lc rgb "forest-green" \
                               title "Measured data"
```
![](images/styles_3.png)



### The useful stuff

We didn't do the measurements for nothing, right? Let's determine particle's speed! 
Guess: `position(t) = V * t + X0`. We want to find `V`. This is called a `fit`.
```gp
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
```
![](images/useful_1.png)

Sample output: `V = 2.96750427263468 +/- 0.051196205997724`. The actual value
of `V` is `3.0`.

-------------------------------------------------------------------------------

We all learned about the "Gaussian" distribution etc. Can we visualize it?
First of all, we need more points.
```gp
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
```
![](images/useful_2.png)

Let's now determine the deviations from the theoretical values.
```gp
set xlabel "Time, t [s]"
set ylabel "Deviation, dX [cm]"
deviation(t, x) = x - position(t)
plot[8:32][-15:15] "data/useful_2.dat" u 1:(deviation($1, $2)) w impulses \
                                       lt 7 lw 0.4 lc rgb "forest-green" \
                                       title "Deviation from theory"
```
![](images/useful_3.png)

I wouldn't call this "clear". A histogram?
```gp
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
```
![](images/useful_4.png)

From `help smooth frequency`:
> The `frequency` option makes the data monotonic in x; points with the same
> x-value are replaced by a single point having the summed y-values.  The
> resulting points are then connected by straight line segments.

So we did the following:
* Make points that fall into the same "bin" have the same values. `bin()`
  does that.
* Make all points have the same y-value -- `1.0`.
* When `frequency` options sums y-value, the result turns out to be the
  number of points inside a bin.

```gp
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
```
![](images/useful_5.png)



### The difficult stuff


#### LaTeX

Suppose now that we want to include the generated plot in our `LaTeX` document.
The `epslatex` terminal lets us create `EPS` images that can be included in the
`PDF`. If you're thinking "I can include a `PNG` file..." don't! `PNG` images 
don't scale, your fonts won't match etc.

Here are the steps:
* Use `epslatex` terminal;
* Set output to `*.tex` (substiture `*` by your desired file name);
* Copy both the `*.tex` and `*.eps` files to where your LaTeX document relies;
* `\include` the `*.tex` file into your LaTeX document.

`LaTeX` document:
```tex
% This is a minimal LaTeX file that shows how to include
% Gnuplot generated plots in your document.

\documentclass{article}

\usepackage{graphics}        % This package is required.
\usepackage{nopageno}        % Allows to get rid of the page numbers,
                             % we don't need them here.

\begin{document}
\begin{center}
\input{my-plot.tex}          % NB: you include the `.tex` file,
                             %     but `.eps` file is needed, too!
\end{center}
\end{document}
```

`Gnuplot` script:
```gp
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
```
![](images/difficult_1.png)

There are a couple of issues with this plot:
* Our keys overlap;
* There is no grid on the left of the formulas;
* We didn't specify the calculated values of `mu` and `sigma`.

Here is a possible fix:
```gp
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
```
![](images/difficult_2.png)


### Good websites
* http://www.gnuplotting.org/
* http://lowrank.net/gnuplot/index-e.html 


### Exercises

Guess what? `Gnuplot` is installed on all university computers except for those at the "Practicum Afdeling"... Fortunately, this problem is easily solved:
* Go <a href=https://sourceforge.net/projects/gnuplot/files/gnuplot/5.0.5/>here</a>
* Download `gp505-win64-mingw.zip` or `gp505-win32-mingw.zip`.
* Unpack the archive.
* `Gnuplot` binary is `gnuplot/bin/gnuplot.exe`.

__NB__: When you've started `Gnuplot`, type `pwd` to see the present working directory. You can change it using the `cd` command. Example: `cd "C:\Documents\Gnuplot-tutorial"`.


##### 0. [__Very easy__]
Try plotting the following functions:

* `f(x) = x**2`
* `f(x) = atan(x)`
* `f(x) = sin(2 * x)**2 * exp(-x)`
* `f(x, y) = exp(-(x**2 + y**2)) * sin(x**2 + y**2)`. __Hint__: use `splot` in place of `plot`. You may also want to set isosamples to a somewhat higher value (`set isosamples 40,40`) and use a smaller range (`splot[-3:3][-3:3] ...`).


##### 1. [__Easy__]
Create a plot of the following function:
```
       / 1 / x**2,   for x >= 1
f(x) = |
       \ x,          for 0 <= x < 1
```
for `xrange = [-0.1:3]` and `yrange = [-0.1:1]`. __Mind you__: this is a single function, thus there should only be one value for each `x`. __Hint__: Read about the <a href=https://en.wikipedia.org/wiki/%3F:#Conditional_assignment>ternary operator</a> and use the fact that `1/0` means _undefined_ in `Gnuplot`, i.e. `plot 1/0` will plot nothing.


##### 2. [__Easy__]
Solve the following equation numerically: `x = 5 * (1 - exp(-x))` (this equation comes from WC1 in QM1a). __Hint__: `x11` terminal (if you're running Linux) and `wxt` (if you're running Windows) can be used interactively:

* mouse scrolling to move up / down,
* Shift + mouse scrolling to move left / right,
* Ctrl + mouse scrolling to zoom in / out,
* ...


##### 3. [__Hard__]
One of the first exercises in the "Statistical Mechanics" course at Radboud is to check how good the <a href=https://en.wikipedia.org/wiki/Stirling's_approximation>Stirling's approximation</a> (`stirling(n) = n*log(n) - n`) of `log(n!)` is. Please, create a nice plot that illustrates it. Cover all `n`s from 0 to 100. You still need to make sure that your plot shows how bad the Stirling's approximation is for small `n`s. For inspiration, have a look at the following picture:

![](images/exercise_3.png)

__Hint__: read about `multiplot` functionality of `Gnuplot`: `help multiplot`.


##### 4. [__Hard__]
During the "KE17: Bepaling van elementaire lading volgens Millikan" experiment you have to measure the velocity of oil particles with the help of a camera. After using `Tracker` you most likely have something like this:

![](images/exercise_4.png)

Next step is to measure velocities `V1` and `V2` of going up and down respectively. Please, create a script to automate this process. Define a variable `file_name` to make it easier to run this script for multiple files. Make no assumptions about the `x` and `y` ranges or the position of the maximum. Print the velocities and errors to standard output. Also create a `.png` file so that you can check how good the fits are. __Hint__: read about `stats` and `every` commands.

__Bonus__: Generalize your script to account for multiple ups and downs. Average over obtained results. I'd be really impressed if you can figure this one out. __Hint__: read about `do` blocks and `if/else` statements.
