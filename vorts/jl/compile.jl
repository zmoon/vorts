# Create a Julia sysimage for faster startup when calling from Py
# Can take a few minutes and be 100s of MB

using Pkg; Pkg.activate("./")

using PackageCompiler


create_sysimage(
    [:DifferentialEquations],#, :Plots],
    sysimage_path="sys_vorts.so",
    precompile_execution_file="vorts.jl",
)
