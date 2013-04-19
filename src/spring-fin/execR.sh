#!/bin/bash
# Lines starting with #PBS are treated by bash as comments, 
# but interpreted by qsub as arguments.  For more details 
# about usage of these arguments see "man qsub"

# Name the job.

#PBS -N DIME-compare

# Set a walltime for the job. The time format is HH:MM:SS

# Run for 1 hour:
#PBS -l walltime=4:00:00

# Select one node, and only 1 processors per node
#you can see what properties are available
# with the "pbsnodes -a" command which will list 
# all nodes and their properties

#PBS -l nodes=1:ppn=1

# # Join the Output and Errors in one file. you can also 
# # set the # path for this output.
# (see "man qsub" for details.)

#PBS -j oe

#
# Source the Dotkit init so you can pull in extra software 
# via the "reuse" command:

. /curc/tools/utils/dkinit

#
# cd to the jobs working directory, which you can set above 
# with a #PBS # directive,
# (see "man qsub" for details)

# cd $PBS_O_WORKDIR


use Torque
use R-2.13

#
# cd to the jobs working directory, which you can set above with a #PBS # directive,
# (see "man qsub" for details)

cd /projects/diesh

#
# Execute the program.
R --save < DIME-compare.R

# This script needs to be submitted via qsub to run on the cluster.