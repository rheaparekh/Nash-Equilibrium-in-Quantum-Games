# Nash Equilibrium in Quantum Games

This is a brute force method similar to gradient descent method inorder to find Nash Equilibrium points in random quantum games.

## Prerequisite

This is a MATLAB program. This program requires QETLAB (Quantum Entanglement Theory LABoratory) which is a MATLAB toolbox for exploring quantum entanglement theory.

To install QETLAB, vist this [page](http://www.qetlab.com/Installation).

## Running the code

* First run the file `PartialTraceModified.m`: This is a modified version of the inbuilt `PartialTrace` function included in QETLAB. The modification allows us to calculate the partial trace of symbolic matrices.

* Next run the file `generate_random_game.m`: This file is use to generate a random quantum game. This file will have two inputs:
  * Number of strategies available for player A
  * Number of strategies available for player B

* Run the file `find_equilibrium.m`: This file will run the brute force algorithm to find the equilibrium for the random quantum game generated in the previous step. The important parameters in this file are:
  * Set `linear_update_method = true` to use the Linear update method and set `linear_update_method = false` to use the matrix exponential update method
  * Set `total_iterations` to a desired value. Current value is `total_iterations = 1500`
  * Set `weight` to desired value. The recommended `weight` in linear update method is in the range 0.1 to 10. The recommended `weight` in matrix exponential update method is in the range 1 to 15.
  * Set the error tolerance by changing parameter `epsilon`.