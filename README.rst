*****
MuJPy
*****

A Python MuSR data analysis graphical interface, based on classes, designed for jupyter.

Released under the MIT licence.

Linux instructions for the time being. VAlid on WIN10, that has a linux shell
* Use virtualenv (see virtualenv.pypa.io).
* Make sure you have python(3) and jupyter, or install them (see jupyter.readthedoc.io).
* Install mujpy. Download from https://github.com/RDeRenzi/mujpy, unzip into the directory of your choice,
   cd mujpy/mujpy/musr2py
   make
   sudo make install
   cd ../..
   # now launch your python3 virtualenv   
   python setup.py install
* Start jupyter:
   jupyter notebook
* In the first cell type::
  >>>%matplotlib
* Add a cell and write
   from mujpy.mugui import mugui
   MuJPy = mugui()
* Run the cells
