This folder contains experimental scripts that must not be included in any release. The whole sandbox folder is actually listed in the .Rbuildignore file to be ignored when building the bundled package. These scripts are intended to the development of new features that eventually will make his way to the package structure. 

The R standard "tests" folder is reserved for its standard usage: functions for testing the package. So, the prototyping of new features must life in the "sandbox" folder.

Naming conventions: in order to easily track the author of files/folders in the "sandbox" folder each file/folder must start with the initials of the author followed by an underscore (e.g. prs_testing_new_things.R). 
