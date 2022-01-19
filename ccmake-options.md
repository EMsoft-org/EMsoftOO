## Compiling sections of EMsoftOO

EMsoftOO contains a large number of programs for a number of electron scattering modalities as well as other modalities (optical microscopy, XRD, ...). Compilation of different modalities can be controlled by means of *CMake* switches that can be accessed readily using the *ccmake* program. When the following command is executed inside the *EMsoftOOBuild/Release* folder in a command window:

```fortran
  ccmake ../../EMsoftOO
```
a simple interface showing name-value pairs appears; a screen shot of the interface on Mac OS X is shown here: 

![Screen shot of ccmake, showing the EMsoft configuration switches](resources/ccmake-screenshot.jpeg)

The parameters are listed on the left, their values on the right inside the dark background box. Use the arrow keys to navigate the cursor to a particular value field; hit *enter* to activate editing of the field if it is a text field, or to toggle its value between all possible values. In an editable field, press *enter* again to leave editing mode.  The highlighted entries represent a number of switches that can be used to control compilation; for instance, to turn compilation of all utility programs off, navigate to the *EMsoft\_enable\_Utilities* switch and toggle its value to *OFF*.  When all parameters have been set correctly, hit the *c*-key to configure, then the *g*-key to generate the new make files; if any asterisks appear in front of value fields, press the *c*-key again until no more asterisks are present. The program will terminate when the generation process ends; recompile the EMsoftOO package using *make -j* or *nmake*, depending on your platform. If only EBSD dictionary indexing is needed, it is sufficient to turn all switches to *OFF*, except for *EMsoft\_enable\_DictionaryIndexing*, *EMsoft\_enable\_SEM*, *EMsoft\_enable\_HDF5\_SUPPORT*, *EMsoft\_enable\_OpenCL\_SUPPORT*, and *EMsoft\_enable\_Utilities*.

Note that the *EMsoftWorkbench* is not compatible (yet) with EMsoftOO; you will need to continue to use EMsoft version 5.0 to build the workbench.
