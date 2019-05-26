The plankton-simulator is a set of matlab scripts to simulate he dynamics of structured microbial communities. 
The simulator is based on a multiscale biophysical model that incorporates species traits, consumer-resource dynamics and stochastic ecosystem assembly. 
The model hierarchically simulates features of ecosystem at different scales by defining three classes in the ecosystem using object-oriented programming in MATLAB: Species, spatially-isolated subcompartments (called Demes from now on) and Plankton (collection of Demes).
Classes are defined by their properties and methods that can be applied to them, this makes it modular and easily extendable.
