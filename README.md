# SibCIOM-coupler
Coupler of SibCIOM model

This program alows to couple 4 climate components: ice, ocean, atmosphere and land.
It is possible to use your own models. To connect your model you should to do next two steps:
  
  1) Initialise your model using mpi procedures in MPI directory;
  
  2) Set up the exchange between coupler and model using mpi procedures from MPI directory.
     (send and recieve functions)
     
Coupler has interpolation procedures that can interpolate data from one frid to another.
