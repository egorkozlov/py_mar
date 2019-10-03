#!/usr/bin/env python3
# -*- coding: utf-8 -*-

How to wrap fortran routines?

This note is based on https://numpy.org/devdocs/f2py/getting-started.html

Here is two example with f_test.f90 and f_modules.f90


1. You need gfrotran and stuff to be working. If you are not sure it does, try
to run execute.command file using bash in compiler-test folder. It contains
some modification of the code I use in this test
2. You need to create a signature file to wrap the f_test.f90. You do this by
typing " f2py f_test.f90 -m Ftest -h f_test.pyf " in the terminal.
3. After you created a signature file you need to compile everything. You do 
this by typing " f2py -c f_test.pyf f_test.f90 "
4. After this, in python you can do " import Ftest ", that imports the module
ftest, that contains the subroutine vfi that is written in f_test.f90 file. 
You can call it as r = Ftest.vfi(sigma), and it returns r.
5. If this worked you can experiment more.
6. Including toolbox.f90 file did not really work well, as comppiler did not 
like some parts of it. There is a chance that by feeding options to f2py we can
overcome this.