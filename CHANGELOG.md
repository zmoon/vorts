# Changelog

## Pre-GitHub change notes

History
* 16-Mar-18: created vorton class. intitial version v0.1, testing the i/o
* 17-Mar-18: v0.2, v0.3  
  ... in between, working on the Python version and wrapper
* 21-Mar-18: v0.4
* 23-Mar-18: v0.4.1

Versions
* 0.2: FT only
* 0.3: added RK4, moved main to separate file and created Makefile
* 0.4: discovered problem with RK4 implementation: this is now fixed following the Python implementation
  - now all dxdt's and dydt's calculated at once  
    and these calculations no longer part of vorton class, following the Python version
  - vorton class used only for history
  - switched from allocatable arrays in the calculations to static, passing in number of vortons and using size declarations
* 0.4.1: attempts to speed up write-out, which is a bottleneck for runs with large number (50+) of vortons/tracers
