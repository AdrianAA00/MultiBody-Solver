# MultiBody Solver

Multi Body Solver. 

It solves N body problem with given constrains. Integration done with Newmark integration scheme.

Documentation needed to start running the code:
- main_general_general: binary file (compiled Fortran code).
- initial_input: document containing initial conditions of the N body problem. Position, Euler Angles, Velocities and Euler Angle Velocities of each body have to be specified at t=0. Moreover,the constrains of each body are defined.
- iter_input: in each time interation certain information has to be provided to the integrator. Forces and moments applied to each body and forced accelerations of each coordinate (in case there are forced accelerations).

Output files:
- data_numerical_general.txt: output file containing the specified coordined of each body in the code. This is constantly modified as I am currenty developing the code.


Test_Cases_General: folder where test cases are available to test the code. 


Under development. For more questions contact me.
