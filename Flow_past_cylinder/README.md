Flow past an obstruction (cylinder).

D2Q9 Lattice which follows BGK model, default code Reynolds number is set at 100 (U_inlet = 0.1, viscosity=0.05 and radius of cylinder = 25).

Default boundary conditions; (can be later changed accordingly)
Left wall - const velocity inlet.
Right wall - zero gradient.
Top & Bottom wall - periodic boundary condition.

A folder is created which is named "output" which stores all the simulation information.
1st column & 2nd column - x and y
3rd & 4th column - x velocity and y velocity
5th column - magnitude of resultant velocity
6th column - density in the domain
7th & 8th column - Resultant Force experienced by the obstacle (cylinder) in x and y directions respectively

Another bash script named "saveimage.sh" can be run to save all the images after the simulation has been stopped which will produce all images upto the last saved iteration values.

Output files are created and:
1st & 2nd column - x and y values
3rd & 4th column - x & y velocities
5th column - resultant velocity magnitude
6th column - density
7th & 8th column - resultant force experienced by the obstacle (cylinder) in x & y directions respectively
