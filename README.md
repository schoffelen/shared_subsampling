# Towards a more robust non-invasive assessment of functional connectivity

Britta U. Westner, Jan Kujala, Joachim Gross, Jan-Mathijs Schoffelen

## Repository
This repository contains the simulation code and data for the project _Towards a more robust non-invasive assessment of functional connectivity_.

This paper will soon be available at _Imaging Neuroscience_.

## Usage

This software is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version. See the file COPYING for more details.

The code and associated data are referenced here for the interested user. The
code has not been curated to work out-of-the-box. For instance, there are hard
coded file paths in some of the files that are pointing to the filesystem on
which the original computations were performed. This needs to be adjusted to
the local situation. Also, some code was written to run efficiently a large
batch of jobs on a torque compute cluster, which might or might not work on
other compute environments. Also, some of the code relies on compiled mex-files
for efficiency, these are platform (and version of the operating system)
dependent and might require local compilation. 
