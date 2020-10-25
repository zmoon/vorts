# vorts

Integrate a system of *N* point vortices.

Created for the PSU class METEO 523 – Modeling the climate system.

## Equations

<div align="center">
<img src="./examples/img/N-vortex_evolution_equations.png"
    alt="N-vortex system of equations in 2 dimensions."
    width=300>
</div>

## Example visualizations

Releasing tracers into the rotating *N*-vortex system:
<!-- <div align="center"><img src="examples/tracer_art_1.jpg" width=300 alt="Tracer art example 1"></div> -->
![tracer art example 1](./examples/img/tracer_art_1.jpg)
<!-- <div align="center"><img src="examples/tracer_art_2.jpg" width=300 alt="Tracer art example 2"></div> -->
![tracer art example 2](./examples/img/tracer_art_2.png)

Also includes tool for making Poincare plots, given a sufficiently long run with a large number (e.g. 100) of tracers. You can make pretty pictures like this:
![example Poincare section plot](./examples/img/ps_theta60deg.png)

## TODO:

* [x] Add the Poincare output as an option to the Fortran code
* [ ] Convert the two test scripts to a Jupyter Notebook (Jupytext) for example
* [ ] Write Python interface to the Fortran version as practice?
