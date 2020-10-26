URLS=[
"index.html",
"model_f.html",
"model_py.html",
"plot.html",
"py/index.html",
"py/integ.html",
"vortons.html"
];
INDEX=[
{
"ref":"vorts",
"url":0,
"doc":"vorts  point vortex models  vorts Integrate a system of  N point vortices.  [pdoc3](https: pdoc3.github.io/pdoc/) documentation: [zmoon.github.io/vorts](https: zmoon.github.io/vorts)  GitHub: [zmoon/vorts](https: github.com/zmoon/vorts)  Equations     Example visualizations Releasing tracers into the rotating  N -vortex system:     > ![tracer art example 1](https: raw.githubusercontent.com/zmoon/vorts/master/examples/img/tracer_art_1.jpg)     > ![tracer art example 2](https: raw.githubusercontent.com/zmoon/vorts/master/examples/img/tracer_art_2.png) Poincar\u00e9 plots. Given a sufficiently long run with a large number (e.g. 100) of tracers, you can make pretty pictures like this: ![example Poincar\u00e9 section plot](https: raw.githubusercontent.com/zmoon/vorts/master/examples/img/ps_theta60deg.png)  Notes Originally created for the PSU class METEO 523 \u2013 Modeling the climate system."
},
{
"ref":"vorts.model_f",
"url":1,
"doc":"Python driver for the Fortran version."
},
{
"ref":"vorts.model_f.fort_bool",
"url":1,
"doc":"Convert Python boolean to string of the Fortran form.",
"func":1
},
{
"ref":"vorts.model_f.Model_f",
"url":1,
"doc":"Thin wrapper for functionality of the Fortran model in  src/ . In this implementation we communicate with the Fortran program via text files. Initialize and create inputs for the Fortran model. Parameters      vortons : Vortons default: equilateral triangle with all G=1 tracers : Tracers (optional) default: no tracers dt : float time step for the output for the integrators,  dt is used as the constant or maximum integration time step depending on the integration scheme nt : int number of time steps to run (not including t=0) int_scheme_name : str default: 'RK4_3' (handwritten basic RK4 stepper)"
},
{
"ref":"vorts.model_f.Model_f.create_inputs",
"url":1,
"doc":"Create input files for the Fotran model describing the initial condition and the simulation settings.",
"func":1
},
{
"ref":"vorts.model_py",
"url":2,
"doc":"Driver for the Python version."
},
{
"ref":"vorts.model_py.Model_py",
"url":2,
"doc":"Model in Python. Create model. Parameters      vortons : Vortons default: equilateral triangle with all G=1 tracers : Tracers (optional) default: no tracers dt : float time step for the output for the integrators,  dt is used as the constant or maximum integration time step depending on the integration scheme nt : int number of time steps to run (not including t=0) int_scheme_name : str default: 'RK4_3' (handwritten basic RK4 stepper)  int_scheme_kwargs passed on to  integrate_manual or  integrate_scipy see their signatures"
},
{
"ref":"vorts.plot",
"url":3,
"doc":"Plotting routines"
},
{
"ref":"vorts.plot.plot_vorton_trajectories",
"url":3,
"doc":"Plot lines: one for each vorton's trajectory.  kwargs are passed on to  plt.subplots() ",
"func":1
},
{
"ref":"vorts.plot.plot_tracer_trajectories",
"url":3,
"doc":"Plot tracer trajectories.",
"func":1
},
{
"ref":"vorts.plot.ps_data",
"url":3,
"doc":"From full set of data, extract data corresponding to times for the Poincare Section. We find the times when the reference vorton is in a certain place or passing through a certain plane (TODO). Parameters      iv_ref : int index of the vorton to use for reference",
"func":1
},
{
"ref":"vorts.plot.plot_ps",
"url":3,
"doc":"Poincare section plot. Here using the data set of all data.  kwargs are passed on to either  ps_data() (if applicable) or  plt.subplots() ",
"func":1
},
{
"ref":"vorts.plot.frame_only",
"url":3,
"doc":"Remove ticks and tick labels from  ax .",
"func":1
},
{
"ref":"vorts.plot.remove_frame",
"url":3,
"doc":"Remove ticks, tick labels, and frame (spines).",
"func":1
},
{
"ref":"vorts.py",
"url":4,
"doc":""
},
{
"ref":"vorts.py.integ",
"url":5,
"doc":"Integration routines and time steppers in Python. In addition to the sub-par handwritten RK and FT schemes, SciPy RK routines (which are written in Python as well) are also available."
},
{
"ref":"vorts.py.integ.integrate_scipy",
"url":5,
"doc":"Integrate using  scipy.integrate.solve_ivp .  options passed through to  solve_ivp ",
"func":1
},
{
"ref":"vorts.py.integ.integrate_manual",
"url":5,
"doc":"Integration routine for use with my handwritten FT/RK4 steppers. Optional naive adaptive time-stepping. t_eval : array_like times to store not including 0, which has already been stored",
"func":1
},
{
"ref":"vorts.py.integ.calc_lsqd_xy",
"url":5,
"doc":"Calculate intervortical distance $l^2$.",
"func":1
},
{
"ref":"vorts.py.integ.calc_lsqd_diff",
"url":5,
"doc":"Calculate intervortical distance $l^2$.",
"func":1
},
{
"ref":"vorts.py.integ.calc_C",
"url":5,
"doc":"Calculate $C$ at time $l$. $C$ is supposed to be a conserved quantity in this system. - Chamecki (2005) eq. 15 We use deparature from $C_0$ (initial value of $C$) in the adaptive stepping to see if we need to go back and step with smaller dt.",
"func":1
},
{
"ref":"vorts.py.integ.calc_tend_vec",
"url":5,
"doc":"Calculate both x- and y-tend written in a vectorized way.",
"func":1
},
{
"ref":"vorts.py.integ.calc_tend_vec_premesh",
"url":5,
"doc":"Calculate both x- and y-tend vectorized, but must pass in x and y in meshgrid form. Numba doesn't support  np.meshgrid .",
"func":1
},
{
"ref":"vorts.py.integ.calc_tend",
"url":5,
"doc":"Calculate tendencies for each vorton (or tracer). Using explicit loops intead of vector(ized) operations.",
"func":1
},
{
"ref":"vorts.py.integ.calc_tend_one",
"url":5,
"doc":"Calculate tendencies for one position based on others. Using explicit loops intead of vector(ized) operations.",
"func":1
},
{
"ref":"vorts.py.integ.FT_step_1b1",
"url":5,
"doc":"Step using 1st-O forward-in-time. Calculate tendencies / integrate one vorton by one. This doesn't affect the result for FT, but for RK4 it does (since sub-steps are used) .",
"func":1
},
{
"ref":"vorts.py.integ.FT_step",
"url":5,
"doc":"Step using 1st-O forward-in-time. Calculate tendencies / integrate all vortons at once.",
"func":1
},
{
"ref":"vorts.py.integ.RK4_step_1b1",
"url":5,
"doc":"Step using RK4. One vorton at a time.  This doesn't work for RK4 since we have sub-steps where the tendencies due to all other vortons need to be up to date or we introduce error. ",
"func":1
},
{
"ref":"vorts.py.integ.RK4_step",
"url":5,
"doc":"Step using RK4. Whole system at once  matrix math version.",
"func":1
},
{
"ref":"vorts.vortons",
"url":6,
"doc":" Vorton class and  Vortons container class."
},
{
"ref":"vorts.vortons.Vorton",
"url":6,
"doc":"Vorton(G, x, y)"
},
{
"ref":"vorts.vortons.Vorton.G",
"url":6,
"doc":"Alias for field number 0"
},
{
"ref":"vorts.vortons.Vorton.x",
"url":6,
"doc":"Alias for field number 1"
},
{
"ref":"vorts.vortons.Vorton.y",
"url":6,
"doc":"Alias for field number 2"
},
{
"ref":"vorts.vortons.Tracer",
"url":6,
"doc":"Tracer  a vorton with G=0 (no power)."
},
{
"ref":"vorts.vortons.Tracer.x",
"url":6,
"doc":"Alias for field number 0"
},
{
"ref":"vorts.vortons.Tracer.y",
"url":6,
"doc":"Alias for field number 1"
},
{
"ref":"vorts.vortons.Tracers",
"url":6,
"doc":"Collection of  Tracer s. Create tracer collection. Parameters      x, y : array_like (n_vortons,) tracer initial x and y positions"
},
{
"ref":"vorts.vortons.Tracers.n",
"url":6,
"doc":""
},
{
"ref":"vorts.vortons.Tracers.x",
"url":6,
"doc":""
},
{
"ref":"vorts.vortons.Tracers.y",
"url":6,
"doc":""
},
{
"ref":"vorts.vortons.Tracers.state_vec",
"url":6,
"doc":"",
"func":1
},
{
"ref":"vorts.vortons.Tracers.state_mat_full",
"url":6,
"doc":"Full state mat for tracers doesn't include G.",
"func":1
},
{
"ref":"vorts.vortons.Tracers.randu",
"url":6,
"doc":"",
"func":1
},
{
"ref":"vorts.vortons.Tracers.spiral",
"url":6,
"doc":"",
"func":1
},
{
"ref":"vorts.vortons.Tracers.plot",
"url":6,
"doc":"",
"func":1
},
{
"ref":"vorts.vortons.Vortons",
"url":6,
"doc":"Collection of  Vorton s. Create vorton collection. Parameters      G, x, y : array_like (n_vortons,) G: Gamma (strength of the circulation, with sign to indicate direction) In fluid dynamics, circulation $\\Gamma$ is the line integral of velocity or flux of vorticity vectors through a surface (here the xy-plane). x: x position y: y position"
},
{
"ref":"vorts.vortons.Vortons.state_vec",
"url":6,
"doc":"Return flattened state matrix (G not included). Needed to feed to  scipy.integrate.solve_ivp , which requires a 1-d array for the  y0 input.",
"func":1
},
{
"ref":"vorts.vortons.Vortons.state_mat_full",
"url":6,
"doc":"Return full state matrix: G and positions.",
"func":1
},
{
"ref":"vorts.vortons.Vortons.G_col",
"url":6,
"doc":"G as a column vector."
},
{
"ref":"vorts.vortons.Vortons.x",
"url":6,
"doc":""
},
{
"ref":"vorts.vortons.Vortons.y",
"url":6,
"doc":""
},
{
"ref":"vorts.vortons.Vortons.n",
"url":6,
"doc":"Number of vortons."
},
{
"ref":"vorts.vortons.Vortons.C",
"url":6,
"doc":"Calculate $C$.  C = \\sum_{\\alpha, \\beta = 1; \\alpha \\neq \\beta}^{N} \\Gamma_{\\alpha} \\Gamma_{\\beta} l_{\\alpha \\beta}^{2}  $C$ is supposed to be a conserved quantity in this system. - Chamecki (2005) eq. 15, which references Aref (1979)",
"func":1
},
{
"ref":"vorts.vortons.Vortons.H",
"url":6,
"doc":"Calculate $H$, the Hamiltonian of the system.  H = -\\frac{1}{4 \\pi} \\sum_{\\alpha, \\beta = 1; \\alpha \\neq \\beta}^{N} \\Gamma_{\\alpha} \\Gamma_{\\beta} \\ln | r_{\\alpha} - r_{\\beta} |  ",
"func":1
},
{
"ref":"vorts.vortons.Vortons.I",
"url":6,
"doc":"Calculate $I$, the angular impulse of the system.  I = \\sum_{\\alpha = 1}^{N} \\Gamma_{\\alpha} | r_{\\alpha} |^2  ",
"func":1
},
{
"ref":"vorts.vortons.Vortons.theta",
"url":6,
"doc":"Calculate $\\theta$, the action angles Chamecki eq. 19",
"func":1
},
{
"ref":"vorts.vortons.Vortons.plot",
"url":6,
"doc":"Plot the vortons. (Only their current positions, which are all this container knows about.)",
"func":1
},
{
"ref":"vorts.vortons.Vortons.mom",
"url":6,
"doc":"Compute  n -th moment. Parameters      n : int which moment https: en.wikipedia.org/wiki/Moment_(mathematics) abs_G : bool, optional (default False) whether to take the absolute value of G values center : bool, optional (default True) True: evaluate moment wrt. center-of-mass False: evaluate moment wrt. (0, 0)",
"func":1
},
{
"ref":"vorts.vortons.Vortons.cm",
"url":6,
"doc":"Compute center-of-mass using Gamma as mass.",
"func":1
},
{
"ref":"vorts.vortons.Vortons.center_coords",
"url":6,
"doc":"Make center-of-mass (0, 0).",
"func":1
},
{
"ref":"vorts.vortons.Vortons.regular_polygon",
"url":6,
"doc":"Create Vortons with positions corresponding to regular polygon. Parameters      n : int polygon order G : int, array-like, optional Gamma value(s) to use single value or array of values default: 1.0  kwargs are passed on to  vortons.regular_polygon_vertices . See signature there.",
"func":1
},
{
"ref":"vorts.vortons.Vortons.isos_triangle",
"url":6,
"doc":"Create Vortons with isosceles triangle vertices.  kwargs are passed on to  vortons.isos_triangle_vertices . See signature there.",
"func":1
},
{
"ref":"vorts.vortons.Vortons.maybe_with_tracers",
"url":6,
"doc":"Return new  Vortons with the tracers. (Temporary? hack to get full state_vec) If  Tracers is  None , just return  self .",
"func":1
},
{
"ref":"vorts.vortons.points_randn",
"url":6,
"doc":"",
"func":1
},
{
"ref":"vorts.vortons.points_randu",
"url":6,
"doc":"Sample from 2-d uniform.",
"func":1
},
{
"ref":"vorts.vortons.points_spiral",
"url":6,
"doc":"Create spiral of points.",
"func":1
},
{
"ref":"vorts.vortons.points_grid",
"url":6,
"doc":"",
"func":1
},
{
"ref":"vorts.vortons.rotmat_2d",
"url":6,
"doc":"Return rotation matrix for rotation  ang_deg in degrees. For left-multiplication of a column position vector. Note:  scipy.spatial.transform.Rotation can be used for 3-d rotations.",
"func":1
},
{
"ref":"vorts.vortons.rotate_2d",
"url":6,
"doc":"Rotate vector  x by  ang_deg degrees. Either  ang_deg or  rotmat must be provided. Parameters      x : array-like (1-d) the vector to be rotated ang_deg : int, float degrees to rotate  x about the origin positive -> counter-clockwise rotmat : array, shape (2, 2), optional rotation matrix  left-multiplies a column position vector to give rotated position Optionally can pass  rotmat instead to avoid computing it multiple times.",
"func":1
},
{
"ref":"vorts.vortons.regular_polygon_vertices",
"url":6,
"doc":"Regular polygon vertices. Parameters      n : int order (number of sides/vertices) c : 2-tuple / array-like center coordinate of the inscribing circle r_c : float, int radius of the inscribing circle",
"func":1
},
{
"ref":"vorts.vortons.isos_triangle_vertices",
"url":6,
"doc":"Isosceles triangle vertices. With fixed top point (0, 1) and fixed left & right y=-0.5. theta_deg : int, float the two angles between the base and connections to the top point at (0,1) 72 -> Lambda_c (1/sqrt(2 60 -> equi tri Lambda : float in (0, 1] 1 -> equi tri",
"func":1
}
]