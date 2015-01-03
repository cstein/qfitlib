QFITLIB
===

This is an API intended to be used to fit partial atomic charges to a molecular electrostatic potential (MEP or molecular ESP) that is evaluated on a grid. The charges can be determined using constraints on the total charge and the total permanent dipole moment.

This library includes code to generate a simple [Connolly surface](http://en.wikipedia.org/wiki/Accessible_surface_area).

## Basic usage

Basic usage is carried out through the `qfit` module which you can include in your own project using invoking

    use qfit

Then, given a set of atomic coordinates `R` and atomic nuclear charges `Z` the workflow to determine the set of potential derived charges is as follows

    call qfit_initialize( R, Z )

which will initialize the charge fitting library with the default settings. When one has the electron density matrix in AO basis `DAO`, fitting the charges to the ESP is done through

    call qfit_fit( DAO )

which will generate the grid, evaluate the ESP and solve the set of linear equations. To extract the results one simply has to invoke

    call qfit_get_results( Q )

where `Q` is a vector with the length of the number of atoms in the molecule. This will hold the potential derived charges fitted to the ESP. Finally, to free the memory allocated by the `qfit`-library, it is important to invoke

    call qfit_finalize
