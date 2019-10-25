# HuckelTheory
A package for simplifying [Huckel theory](https://en.wikipedia.org/wiki/H%C3%BCckel_method) calculations and visualizations using the Mathematica 12 [Molecule](https://reference.wolfram.com/language/ref/Molecule.html) functionality

## Background 
[Hückel molecular orbital theory](https://en.wikipedia.org/wiki/H%C3%BCckel_method) describes the pi electrons in planar molecules.  For background theory, see Schrier, [Introduction to Computational Physical Chemistry](https://amzn.to/2Jj4Yp4) (University Science Books, 2017), Chapter 6.  A [preview version of this chapter is available at the publisher's website]( http://www.uscibooks.com/schrier.htm)

The code in this package largely follows the approach described in that chapter.  However, as the book is meant to be pedagogical for students who are new to programming and easily implemented in different programming languages, it does not take full advantage of the functional programing paradigm in Mathematica.  In contrast, this package does.  This package also makes extensive use of the [Molecule](https://reference.wolfram.com/language/ref/Molecule.html) functionality that was released in Mathematica 12 (2019).

## Live Demo

![Live demo of a sample calculation](HuckelTheory_livedemo.gif )

## Installation

Note:  To install Packages in Mathematica, go to File >> Install... 

The relevant package file is `HuckelTheory.wl`

## Usage

All functions take a [Molecule](https://reference.wolfram.com/language/ref/Molecule.html) as an input.  Outputs are in units of the C-C coupling index, t.

The `HuckelTheory` package is general enough to treat most heteroaromatic systems (B, N, O, F, Cl, Br), using the Streitwieser parameters. (See [Introduction to Computational Physical Chemistry](https://amzn.to/2Jj4Yp4), Table 6.1, p. 154).

Atoms considered part of the pi system are:
- Any sp2 atom
- oxygens, nitrogens, and halogens connected to the pi system

*Potential limitation:* Charge-density based calculations (bond orders, net charges, etc.) depend on having an even number of pi-electrons.

Functions provided:

- `HuckelHamiltonianMatrix`
- `HuckelMO`
- `HuckelChargeDensityMatrix`

- `HuckelMOPlot`
- `HuckelBondOrderPlot`
- `HuckelTotalElectronsPlot`
- `HuckelNetChargePlot`

The file `HuckelTheory_demo.nb` in this repository demonstrates each of the functions.

## Why not put this on the Wolfram Function Repository?

It's not just one public function, and they are all closely dependent on one another.  Maybe in the future.
