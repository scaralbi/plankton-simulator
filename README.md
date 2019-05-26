# Plankton Simulator

The plankton-simulator is a set of matlab scripts to simulate he dynamics of structured microbial communities. 
The simulator is based on a multiscale biophysical model that incorporates species traits, consumer-resource dynamics and stochastic ecosystem assembly. 
The model hierarchically simulates features of ecosystem at different scales by defining three classes in the ecosystem using object-oriented programming in MATLAB: Species, spatially-isolated subcompartments (called Demes from now on) and Plankton (collection of Demes).
Classes are defined by their properties and methods that can be applied to them, this makes it modular and easily extendable.

# How to run the Simulator

1. Download and extract the repository; alternatively, clone the repository using Git with the command
```git clone https://github.com/scaralbi/plankton-simulator.git```

2. Open MATLAB 

3. Open ```config.m``` 

4. Specify desired parameters (```params```), as described in file ```params_guide.txt```   

5. Run ```config.m``` 

6. Open and execute ```simulator.m```

7. To visualise simulation plots and statistics execute ```data_reader.m``` 



# Overview of the files in the plankton-simulator folder:

* [`Species.m`](https://github.com/scaralbi/plankton-simulator/blob/master/Species.m )
: defines the Species class and the functions that can be applied its objects
* [`Deme.m`](https://github.com/scaralbi/plankton-simulator/blob/master/Deme.m )
: defines the Deme class and its functions
* [`Plankton.m`](https://github.com/scaralbi/plankton-simulator/blob/master/Plankton.m ): defines the Plankton class and its functions
* [`dynamics.m`](https://github.com/scaralbi/plankton-simulator/blob/master/dynamics.m ): contains the ordinary differential equations that describe the population and concentration dynamics within ecosystem
* [`config.m`](https://github.com/scaralbi/plankton-simulator/blob/master/config.m ): script that generates configurations files (config.mat) with speciefied simulation options and parameters
* [`simulator.m`](https://github.com/scaralbi/plankton-simulator/blob/master/simulator.m ): script that reads config files and simulate ecosystem and save data as config_data.mat files
* [`data_reader.m`](https://github.com/scaralbi/plankton-simulator/blob/master/data_reader.m ): script that reads config_data.mat files and produces plots

# simulator.m 
## Stochastic Ecosystem Assembly Algorithm
![alt text](https://github.com/scaralbi/plankton-simulator/blob/master/flowchart.png)

# Species.m 
Every species is defined by metabolic strategies 

```Matlab
classdef Species    
    properties
        name = 0;
        traits = [];
        population = 0;
        budget = 0;
        mutation_rate = 0;
```


# Deme.m 
```Matlab
classdef Deme    
    properties
        nutrient_concentrations = [];
        nutrient_supply = [];
        species composition = [];
```


# Plankton.m 
```Matlab
classdef Plankton    
    properties
        demes = [];
        size = [];
    end
```



