# Testing mixing layer problem with time-periodic mean velocity profile

This case solves the forced Navier-Stokes equations for the velocity fluctuations that arise for an imposed mean mixing layer flow profile $$\bar{u}(y,t)=\frac{U}{2}\left(1+\frac{1}{2}\sin(F t 2\pi)\right)\tanh\left(\frac{y}{2\delta}\right)$$

The Reynolds number is defined as $$\text{Re} = \dfrac{U \delta}{\nu}$$

The forcing scheme is based on [Dhandapani, C. & Blanquart, G. 2020 From isotropic turbulence in triply periodic cubic domains to sheared turbulence with inflow/outflow. PRF 5, 124605](https://doi.org/10.1103/PhysRevFluids.5.124605).

The case is loosely based on [VanDine, A. & Pham, H. T. & Sarkar, S. 2021 Turbulent shear layers in a uniformly stratified background: DNS at high Reynolds number. JFM 916, A42](https://doi.org/10.1017/jfm.2021.212).