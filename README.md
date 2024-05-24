# Trition

This is a finite volume numerical solver employed to solve a set of hyperbolic / parabolic differential equations using both staggered and collocated grid. Currently it supports a 2nd order central difference scheme for incompressible flow and is run in serial. A variety of scenerios has been tested successfully including the lid driven cavity.

## Features

- Toggle for viscosity
- Toggle for staggered grid
- Dimension specification
- Automatic CFL determination
- Reynold,s number/ viscous coefficient prescription

## Instructions

- All the features can be set the `main.h` header file
- Commpile all the files using gcc
- Run the binary
- Output for grid and data will be produced in `.txt` files

## Upcoming update:

- Support for parallelization
- makefile support
- .vtk ASCII output
