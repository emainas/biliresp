# This module performs iterative restraint fitting of the unrestrained ESP charges using a hyperbolic restraint
# This module needs to be written with flexibility in mind. The restraint fitting can be done in many different ways
# and it is not clear which way is the best. Therefore, the module should be able to adapt to different ways of restraint fitting.

# We will need the following functions for the most basic implementation:
# 1. A function to read the unrestrained ESP charges from a file
# 2. A restraint potential function 
# the functional form of this is given under root/extra/restraint_potential.py -- It takes two arguments a and b. This restraint it only applied to 
# non-hydrogen atoms so we will need to load the xyz file to get the atom elements, we have capability for that so far.
# 3. loss function to minimize
# 4. The algorithm that performs the iterative fitting - gradient descent, conjugate gradient, etc.
# for this we will need choice of algorithm, epochs, learning rate, tolerance, etc.
# 5. A function to load the answer from Terachem's resp.out in the raw data and compare with it to see if we are close. 
# 6. Appropriate logging and saving of intermediate results (for example loss at every step) and reporting of them to make sure we are doing it right.
# 7. The hardest part of this code will be to integrate it with the symmetry module. After permorming simple RESP fitting, we will need to
# identify the symmetry groups and then apply the restraints accordingly. This will require some thought and planning. 