This test case simulates an idealized left ventricle with a Holzapfel-Ogden material model
contracting due to time-dependent active stress, and subject to a time-dependent
pressure load on the endocardial surface. The full problem is described in
case 1A of the cardiac elastodynamcis benchmark paper by Aróstica et al. (2025)[1]. A comparison of the displacement of two points throughout the cardiac cycle as computed by multiple solvers including svMultiphysics is shown below:

![Displacement Benchmark](comparison_plots_p0_p1_step_1_nonblinded.png)

Aditionally, we present a pressure-volume loop for the idealized left ventricle. Note that the time-dependent pressure load in this problem is not intended to reflect physiological conditions.

![P-V loop](p-v_loop.png)

These plots can be generated in figures.py by dowloading the benchmark dataset from https://zenodo.org/records/14260459 and running solver.xml with 1000 time steps to obtain the svMultiphysics results for the entire cycle.

[1]: Reidmen Aróstica, David Nolte, Aaron Brown, Amadeus Gebauer, Elias Karabelas, Javiera Jilberto, Matteo Salvador, Michele Bucelli, Roberto Piersanti, Kasra Osouli, Christoph Augustin, Henrik Finsberg, Lei Shi, Marc Hirschvogel, Martin Pfaller, Pasquale Claudio Africa, Matthias Gsell, Alison Marsden, David Nordsletten, Francesco Regazzoni, Gernot Plank, Joakim Sundnes, Luca Dede’, Mathias Peirlinck, Vijay Vedula, Wolfgang Wall, Cristóbal Bertoglio,
A software benchmark for cardiac elastodynamics,
Computer Methods in Applied Mechanics and Engineering,
Volume 435,
2025,
117485,
ISSN 0045-7825,
https://doi.org/10.1016/j.cma.2024.117485.