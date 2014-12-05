hk-vectorgen
============

Script for generating vectors for possible HK near detectors. 

There are two main scripts: run_neut.py and run_genie.py.

NEUT
====

There are two steps to generating NEUT MC.
  1. Compute maximum interaction probabilities for the given geometry and interaction model.
  2. Generate events.

The script "run_neut.py" generates a TITUS geometry file and runs both steps. 
Each step is evaluated lazily, if the target files already exist the step is skipped.

Prerequisities
--------------
  * You must have compiled and setup CERNLIB.
  * You must have compiled and setup NEUT (both neutgeom and neutsmpl sub-packages).
  * The environment variables:
    * NEUT_ROOT points to NEUT_ROOT directory containing "src/neutgeom" and "src/neutsmpl".
    * NEUTGEOM (points to NEUT_ROOT/src/neutgeom).
    * RANFILE (points to ${NEUT_ROOT}/neutsmpl/random.tbl).

Setup
-----
  * Source setup.sh in the root directory of this package. This will set the PYTHONPATH to include the vectorgen package.
  * You must have a set of JNUBEAM flux files. At present these are hard-coded in run_neut.py. You will need to edit this script to point at the correct files.

Run
---

The script "run_neut.py" will run all steps required to generate NEUT files in the TITUS near detector. You can run it with the command:

    python -m vectorgen.run_neut "polarity" "radius" "length"
  
  where
  
    polarity = -1 or 1 for RHC and FHC modes.
    radius = radius of detector in metres (this is interpreted as side length for cuboid detector geometry).
    length = length of detector in metres.
    
  In addition you can supply the optional arguments:
  
    --pdg= -12 -14 12 or 14 which will generate only neutrinos with that PDG code.
    --geometry= cylinder or cuboid
    --copyflux will copy flux files to /tmp which may improve IO performance on some systems.

for example the following command will generate anti-nue events in reverse-horn-current mode with a cylindrical geometry with radius=5.51m length=22.0m:

    python -m vectorgen.run_neut 1 5.51 22.0 --geometry=cylinder --pdg=-12

If you have an LSF batch system the script will automatically submit jobs when provided with the following flag:

    --batch

Batch mode has only been tested at Warwick. The scripts may need some alteration to run elsewhere.

GENIE
=====

GENIE works in a similar way to NEUT. The script "run_genie.py" takes most of the same arguments as described above for "run_neut.py"

Prerequisities
--------------
  * The GENIE environment variable must be set.
  * You must have compiled and setup GENIE.
  * The T2K support binaries must have been included in compilation. You can check this by looking for $GENIE/bin/gevgen_t2k. If this does not exist you may need to re-run configure with the enable-t2k option switched on and then recompile.
  * The GENIE_SPLINES enivronment variable must set and must point to a valid splines file.
  * The script uses an alternate event generator list called "DefaultWithMEC" which is the same as the default model, but has MEC switched on. This is enabled by adding the following to $GENIE/config/EventGeneratorListAssembler.xml:
~~~
  <param_set name="DefaultWithMEC"> 
     <param type="int" name="NGenerators">   14                                 </param>
     <param type="alg" name="Generator-0">   genie::EventGenerator/QEL-CC       </param>
     <param type="alg" name="Generator-1">   genie::EventGenerator/QEL-NC       </param>
     <param type="alg" name="Generator-2">   genie::EventGenerator/RES-CC       </param>
     <param type="alg" name="Generator-3">   genie::EventGenerator/RES-NC       </param>
     <param type="alg" name="Generator-4">   genie::EventGenerator/DIS-CC       </param>
     <param type="alg" name="Generator-5">   genie::EventGenerator/DIS-NC       </param>
     <param type="alg" name="Generator-6">   genie::EventGenerator/COH-CC       </param>
     <param type="alg" name="Generator-7">   genie::EventGenerator/COH-NC       </param>
     <param type="alg" name="Generator-8">   genie::EventGenerator/DIS-CC-CHARM </param>
     <param type="alg" name="Generator-9">   genie::EventGenerator/QEL-CC-CHARM </param>
     <param type="alg" name="Generator-10">  genie::EventGenerator/NUE-EL       </param>
     <param type="alg" name="Generator-11">  genie::EventGenerator/IMD          </param>
     <param type="alg" name="Generator-12">  genie::EventGenerator/IMD-ANH      </param>
     <param type="alg" name="Generator-13">  genie::EventGenerator/MEC-CC       </param>
     <!-- <param type="alg" name="Generator-14">  genie::EventGenerator/MEC-NC       </param> -->
  </param_set>
~~~ 


WCSIM
=====
