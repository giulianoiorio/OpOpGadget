# OpOpGadget
The aim of this module is to provide a number of tools  for the generation of initial conditions and analysis of N-body simulations. It can manage files in the format used by the N-body simulations for the codes FVFPS  and Gadget-2. 


**WARNING:**  OpOpGadget and its documentation are still under development. Some examples can be found in the folder 'tutorial'.
Comments and suggestions are welcome! If you are interested in using it and/or in introducing new functionalities you can send me an email (giuliano.iorio89@gmail.com).

## Install

### Requirements 
This module is for Python3 only. It also uses C to build Python extensions, so you need  a C compiler too. 

Before to install OpOpGadget, you have to install the following packages:

- numpy ([https://numpy.org](https://numpy.org))
- scipy ([https://scipy.org](https://scipy.org))
- matplotlib ([https://matplotlib.org](https://matplotlib.org))
- colossus (use pip install colossus, or [https://bitbucket.org/bdiemer/colossus/src/master/](https://bitbucket.org/bdiemer/colossus/src/master/))
- roteasy ([https://github.com/giulianoiorio/roteasy.git](https://github.com/giulianoiorio/roteasy.git))
- FERMI (pip install fermi)
- Cython (pip install cython)
- CythonGSL ([https://github.com/twiecki/CythonGSL.git](https://github.com/twiecki/CythonGSL.git))

Finally, install OpOpGadget running pip install . or python setup.py install
