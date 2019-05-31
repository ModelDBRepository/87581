This is the readme for the models associated with the paper

     Edwards FR, Hirst GD (2006) An electrical analysis of slow wave
     propagation in the guinea-pig gastric antrum. 
     J Physiol 571:179-89

which is the companion paper to another which describes the
physiologically acquired data:

     Hirst GD, Garcia-Londono AP, Edwards FR (2006) Propagation of
     slow waves in the guinea-pig gastric antrum.
     J Physiol 571:165-77

Usage:

These Matlab files generate parts of figure 4 - a
longitudinal string of 20 ICC compartments with attached longitudinal
muscle compartments.

This file is called Figure4.m and requires ICs180.dat.  If
you run it as it is, it will produce figure 4A as Matlab figure 50
(green line) after about 10 seconds on my Dell laptop (3/29/07).

If you modify line 26 to read PASSIVE=0, then it will produce a
version of figure 4C as Matlab figure 51 after about 37 minutes.

If you subsequently modify line 45 to read UNITS=0, and rerun it then
it will produce a smooth looking version of figure 4C after about 1
minute.  This option uses simple analytic functions to approximate
conductance modulations, rather than random number generated
stochastic discharges of unitary conductance modulations.  It runs a
lot faster and was useful for shaping things up.  Also useful for
quick verification.

Also attached is Figure6.m which requires FCs134.dat.  This models a
circumferential string of 16 compartments of antral circular muscle,
with the first 6 connected to ICC compartments each with an attached
longitudinal muscle compartment.  Driving the first ICC compartment is
a corporal circular muscle compartment.

Run as it is, it will produce figure 6A as Matlab figure 20 (green
line) after about 15 seconds.

If you modify line 31 to read PASSIVE=0, then it will produce a
version of figure 6C as Matlab figure 21 after about 37 minutes.

If you subsequently modify line 51 to read UNITS=0, and rerun it then
it will produce a smooth looking version of figure 6C after about 45
seconds.

The only flagrant contradiction I have noticed in the paper is on line
13 of p184.  The value in brackets should read 10 mV, not 1 mV.  This
value appears as NUDGE in both Matlab files.

I apologize in advance for the lack of elegance of the runnable files.
They are replete with unused options and commented out tracts to draw
intermediate graphs and stuff.  Many parameters are dependent or
various multiplications, and so appear cryptic.... my target was
development speed rather than Matlab's programming advantages.

PlatCond and CMBCond are the most interesting (more or less identical)
subroutines.  They define the time locations of stochastic unitary
conductance modulations whose mean frequency varies with membrane
potential.

Dr Frank Edwards
3/29/07
