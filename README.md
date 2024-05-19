Hello, reader!  These are Becca's notes for compiling and running
G4d2o (on Ubuntu Linux).  I've listed four (and a half) steps, but if
you already have the CRY module installed, skip ahead to step 3.

I should note here: I use cmake and makefiles, but using Xcode users
should follow the same procedure -- simply replace the "makefile" in
Step 3's command line arguments with "Xcode"!

*********************************************************************************
Step 1: Install CRY module!

I downloaded the "Install from source" package from
[https://nuclear.llnl.gov/simulation/main.html].  Use the makefile to
install the module.  Note that git is currently NOT tracking the CRY
module, so I recommend installing it outside the G4d2o directory so
that our terminals won't be grumpy geese if we try to change branches
and forget that we installed the module.

*********************************************************************************
Step 2: Set all environmental variables

In the location where you built the CRY module, there will be a setup
file.  (Note that you might have to run ./setup.create if you changed
the directory).  In setupBuild.sh, change the path on line 3 to the
location of this setup file.  While you have setupBuild.sh open,
adjust the paths for the environmental variables one lines 4 and 5 to
find the two root files -- everything after G4d2o should be the same.

NOTE: You may have to re-run the setupBuild.sh when opening a new
terminal -- if you are annoyed by this, feel free to copy Lines 3, 4,
and 5 into your terminal setup scripts!

*********************************************************************************
Step 3: Build the G4d2o application

There exist two bash scripts to construct the executable, used with
the following command executions:
. ./setupBuild.sh makefile
./compileApp.sh makefile

*NOTE*: The executable does NOT run from the G4d2o-build
   directory.  You should run the executable from the main project
   directory.  The executable is copied to the project directory
   within compileApp.sh!!
   
-----------------------------------------------------------------------------------------------------
Step 3.5: If you get linking errors when "make"-ing your
application...

Edit CMakeLists.txt in the G4d2o main directory.  I needed to make 3
changes:

1. Anywhere ${CRYHOME} appears, change it to $ENV{CRYHOME}
2. Find where we are GLOBing the sources and headers, and add the following line:
"file(GLOB crylibs $ENV{CRYHOME}/lib/libCRY.a)"
3. Find where target_link_libraries occurs.  Comment out the existing, and replace it with:
"target_link_libraries G4d2o ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} simEvent ${crylibs}"

*********************************************************************************
Step 4: Run the program!!

The biggest thing to recognize is that beamOn.dat controls the sim.
Make sure the PhysicsOn variable is set to true!

"./G4d2o CommandLine" will tell you the current setting of the
variables that will be implemented when using ./G4d2o
