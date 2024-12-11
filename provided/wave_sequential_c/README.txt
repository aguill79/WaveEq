To compile serial wave file type
gcc wave.c -o wave_exc -O3

To create an anomated gif file for
a 64x64 playfield for 500 iterations
gen_gif.sh 64 500
the gif file that is created is named
wave.gif which can be observed with
a standard web browser.

In general you would type 
gen_gif.sh <Size of square wave field> <Number of interations>

Note if in the middle of the script it stops with a message
similar to that shown below:
X11 connection rejected because of wrong authentication.

gnuplot: unable to open display 'localhost:12.0'
gnuplot: X11 aborted.

Then just hit return and it should proceed onward. This message
is that gnuplot does not have a valid X11 connection which
is not needed for gif file creation.

Also Note: The prespective that is present artificially 
renders a line or two on the plot which makes it appear
that there is a wave ridge. It becomes clear soon after
the animation begins that this is not a real wavefront.
The idea is that you should get very similar plots if
you switch out the wave_exc with your parallel version.
I would use your visual inspection of similar plots 
for both the serial and parallel versions to be used
for verification of correctness but have your timing
come from your parallel version executed in a stand
alone manner.

BEW

