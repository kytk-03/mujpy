*****
MuJPy
*****

A Python MuSR data analysis graphical interface, based on classes, designed for jupyter.

Released under the MIT licence.

Use virtualenv (see virtualenv.pypa.io).
Make sure you have python(3) and jupyter, or install them (see jupyter.readthedoc.io).
Install mujpy.

Start jupyter:
>>> jupyter notebook
In the first cell type
>>>%matplotlib notebook 
Add a cell and write
>>>from mujpy.mugui import mugui as MG
>>>MuJPy = MG()
>>>MuJPy.start()
Run the cells
