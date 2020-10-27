# biobuilders-hyphae-model

Repository for Mycemulator by DTU BioBuilders 2020.
====================

Current version: Mycemulator 1.0
--------------------


Welcome to the home of Mycemulator - the fungal mycelium growth simulator made by DTU Biobuilders as part of the 2020 iGEM competition.


Description
-----------

Mycemulator is a tool for making growth simulations of mycelial networks in filamentous fungi. It is a stochastic model designed to incorporate a range of experimentally determined parameter. The default values come from experimental data collected by the 2020 DTU Biobuilders.

Getting started
---------------

### Dependencies

Mycemulator is written in Python 3.6 and uses common Python 3 packages, all of which can be found in the requirements.txt file in this repository.

### Installation

To install, clone the repository or simply download the newest version of the simulation.

```
mycemulator_1.py
```
Make the file executable using

```
chmod +x mycemulator_1.py
```

### Executing the program

The program can be run as is or can be modified by a number of additional arguments.

#### Optional arguments
```
-ani
```
Name of folder of time-dependent png images and gif animation.
Type: string. Default: None
                                                                                 
```
-pdf
```
Name of pdf of time-dependent subplots.
Type: string. Default: None

```
-img
```
Use parameters from csv-output of image analysis tool.
Type: string.

```
--hours
```
Number of hours for which to simulate exponential hyphal growth.
Type: float. Default: 12

```
-tstep
```
Hours per simulation round.
Type: float. Default: 1/60 (1 minute per round)

```
-field
```
Substrate field type. Choose whether to use a uniform field (u) or a radial gradient field (g).
Type: string. Default: g (gradient field)

```
-source
```
Determines whether the simulation area will have an infinite substrate source at its edge.
Type: boolean. Default: False (no source)

```
-mu_max
```
Maximum growth rate for the given strain in exponential phase.
Type: Float. Default: 0.2757333 (biolector measurement for ATCC 1015)

```
-q
```
Branching frequency per hyphal element per hour (if not taken from image analysis output).
Type: float. Default: 0.016925731922538392 (taken from image analysis of ATCC 1015")

```
-lat_sub_min
```
Minimum substrate concentration for lateral branching to occur.
Type: float. Default: 10

```
-ap_sub_min
```
Minimum substrate concentration for apical branching to occur.
Type: float. Default: 10

```
-ext_sub_min
```
Minimum substrate concentration for tip extension to occur
Type: float. Default: 5

```
-S0
```
Initial uniform substrate concentration. Also boundary condition if '-source' is applied.
Type: float. Default: 100

```
-D
```
Diffusion constant.
Type: float. Default: 0.99

```-S_tip
```
Substrate consumption per hyphal tip per simulation round.
Type: float. Default: 10

```
-S_nontip
```
Substrate consumption of a non-tip hyphal element
Type: float. Default: 1

#### Running the program

Run the program from the terminal.

```
./mycemulator_1.py <additional arguments> 
```

Acknowledgements
---------------
Original R model:
2018 DTU Biobuilders
Chris Workman

