#!/bin/bash

# Visualization window geometry
window 0 0 600 400

# Initial solution
mesh Q.mesh Perturbed_0000.gf

# Setup the GLVis scene. Executed after pressing the space bar.
{
  keys R
  zoom 1.5
}
