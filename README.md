This code is the author's implementation of the algorithm presented in:

Carl Doersch, Abhinav Gupta, and Alexei A. Efros, "Mid Level Visual Element 
Discovery as Discriminative Mode Seeking" in NIPS 2013.

The majority of the code was written by Carl Doersch (cdoersch at cs dot cmu dot edu),
although there are major contributions from Saurabh Singh and some from others who
are noted in the code.

This is officially unsupported research code, but our goal is that this will be
as useful as possible.  You are encouraged to ask questions via e-mail,
and strongly encouraged to give feedback if you find the code 
counter-intuitive.  I plan to update this code as issues are discovered.

General setup:

0) This code is only designed to run on linux.  It might run on a mac or other
   unix system since it only  uses unix commands, but it has not been tested.
   make sure you clone it with 'git clone --recursive' to pull down dswork.

1) Install libsvm and edit the file 'myaddpath.m' so that it adds libsvm to the path.

2) Check that you have the statistics toolbox.  If you plan to use the inter-element
   communication and are not running 64-bit linux, you need to copy the file
   [MATLAB ROOT]/toolbox/stats/stats/private/linkagemex.(your extension)
   onto the matlab path (e.g. in the directory containing this file).  This is to
   work around a bug I've experienced in Matlab's linkage implementation; if 
   you don't want to use the workaround, you can switch the call in findOverlapping3.m
   from linkage2 to linkage.

3) Edit the indoor67_main file with the path of your download of the indoor67 dataset
   (there's directions in the file).  If you want to use your own data, setdataset.m gives 
   instructions on how to do this.

4) If you aren't running on 64-bit linux, go to the hog/ directory
   and run 'mex features.cc'.

5) [optional] If you want to measure coverage pixel-wise (the code doesn't by default), compile
   the code in the clipper directory--i.e. cd into the directory and run
   "mex greedySelectDetrsCoveragemex.cpp clipper.hpp clipper.cpp"

6) Disable the matlab toolbox cache, either in settings->General->uncheck "Enable toolbox path cache"--
   and exiting matlab to save your changes--,
   or by modifying .matlab/VERSION/matlab.prf and adding the line "GeneralUseToolboxCache=Bfalse"
   (or changing the line for GeneralUseToolboxCache if it exists).
   The code may work without this change, but by default matlab writes a toolbox cache file
   to the home directory when it exits.  If many matlabs exit at the same time,
   they can corrupt this file and cause some matlab processes to error out.  The dswork
   code restarts the workers occasionally, which means the execution will get stuck if
   matlab can no longer start properly.

7) Run the clustering code in Matlab.  The main script for the indoor67 experiment
   is indoor67_main.m.  

Running the code on indoor67 should require about 8GB of RAM per machine.  indoor67_main
will attempt to estimate the number of jobs to run on each machine based on the RAM
each machine has.  The program also needs about 30GB of free disk space and about 300GB
in the temporary local directories (that's total; it can be distributed across different
machines).  If you're not using local directories, then that 300GB will need to be on
the shared filesystem.


Parallelizing the code with dswork:

This codebase uses heavily the dswork framework.  The README in the dswork
directory gives full documentation, but here's a tl;dr summary.  

dswork has two main features.  (1) it establishes a mapping between
some directory on the filesystem and the variable 'ds' in your
workspace.  Hence, you can call 

dssetout('/tmp');
ds.mydirectory.myvariable=rand(100);
dssave;

This causes the variable ds.myvariable to be saved to '/tmp/ds/mydirectory/myvariable.mat'.
dswork supports filesystem command analogous to unix, including dsmv, dsdelete,
dssymlink (though this implementation is incomplete), and dscd.  To make the syntax
as concise as possible, the format that things are saved in depends on the variable
suffix--thus far, the suffixes img and html and txt have special meanings.

On top of this, dswork supports some basic distributed processing features, including
multiple MATLAB's on one machine, and multiple matlabs on different machines.  For these
to work, the directory where dswork saves its files needs to be shared among all machines
you are using.

At a high level, dsmapredopen() sets up a pool of workers that are essentially stateless.  
Using dsrundistributed() or dsmapreduce() will assign work to each worker, allows the
workers to load data from the shared storage, and tracks the variables that
get written.  Additionally, dsmapreduce() supports mapreduce-like communication between
the workers that does not go through the shared storage.  This instead uses ssh to
share data.  This is important because, if you share all your data using an nfs-share 
of a disk on a single machine (as I did), I/O will likely become a bottleneck for the 
main element mining loop.  Sharing data directly between machines alleviates this
issue.  To use this, the only setup required should be to give dswork a local directory
where it can save temporary files on each machine, via the dssetlocaldir() funciton. 


Running multiple matlab's locally should not require any additional setup, but running
distributed will require a machine that supports qsub.  All of the experiments
for this project were performed using Starcluster on EC2, which sets up an
OGS cluster with data shared over nfs.  See dsmapredopen for instructions on
starting the distributed session.


So you want to understand the code...

My coding style is developed around rapid prototyping, and is probably different from
what you're used to.  Here's a few patterns that I tend to use.

1) I generally use parallel arrays where other programmers would use arrays of structs or
   arrays of objects.  This is the case because I often need quick access to all values of a single
   field.  Matlab's struct arrays support this, but it is extremely inefficient.
   The distributeby/invertdistributeby have become a sort of swiss army knife for handling
   parallel arrays in my code.  You should memorize what distributeby does.

2) To ease dealing with parallel arrays, the effstr... commands are designed to deal
   with a struct holding multiple parallel arrays (effstr means 'efficient replacement
   for matlab struct arrays').  The motivation is that I can add temporary data to an
   object and keep track of it alongside those objects, all with minimal modification
   of the code.  

3) If I have a collection of n bounding boxes, they will be stored in an n-by-8 array with
   the following column order: [x1 y1 x2 y2 detection_score detector_id image_id flip].
   x- and y- coordinates are in terms of pixels *in the space returned by getimg*.  flipped
   detections have flip=1, but are still in terms of the coordinates of the un-flipped image.
   These are used so frequently in the code that this format is used without comment; you
   should memorize the order.

