A low residual criteria (1e-2) is used for demonstrating purpose.
For a "real simulation", it is recommend to at least use a criteria of 1e-3.
This can be set the system/fvSolution-file.

Time on my laptop (using 4 cores):
1e-2: 9 seconds
1e-3: 96 seconds

In a study of the single V80 turbine, it was found that the power changed by 
around 1% by going from 1e-2 to 1e-4. Going from 1e-3 to 1e-4 only changed the 
power by around 0.1%, hence the recommendation of 1e-3. Also, the wake profiles
were found to be "residual-independent" at around 1e-3 to 1e-4.
However, for full wind farms a more strict residual criteria might be needed
due the more complex flow interaction (to be investigated quantitatively).
