# Plankton Simulator

The plankton-simulator is a set of matlab scripts to simulate he dynamics of structured microbial communities. 
The simulator is based on a multiscale biophysical model that incorporates species traits, consumer-resource dynamics and stochastic ecosystem assembly. 
The model hierarchically simulates features of ecosystem at different scales by defining three classes in the ecosystem using object-oriented programming in MATLAB: Species, spatially-isolated subcompartments (called Demes from now on) and Plankton (collection of Demes).
Classes are defined by their properties and methods that can be applied to them, this makes it modular and easily extendable.


# Overview of the files in the plankton-simulator folder:

* Species.m: defines the Species class and the functions that can be applied its objects
* Deme.m: defines the Deme class and its functions
* Plankton.m: defines the Plankton class and its functions
* dynamics.m: contains the ordinary differential equations that describe the population and concentration dynamics within ecosystem
* config.m: script that generates configurations files (congig.mat) with speciefied simulation options and parameters
* simulator.m: script that reads config files and simulate ecosystem and save data as config_data.mat files
* data_reader.m: script that reads config_data.mat files and produces plots




