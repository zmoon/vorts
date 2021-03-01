# Create a Julia sysimage for faster startup
#
# The env should be activated when runinng:
# julia --project=. compile.jl
#
# Note that it can take a few minutes and the sysimage can be 100s of MB
#
# To try it out:
# julia --sysimage sys_vorts.so --project=. -e 'include(\"vorts.jl\"); main()'

using PackageCompiler


create_sysimage(
    [:DifferentialEquations, :Plots],
    sysimage_path="sys_vorts.so",
    precompile_execution_file="vorts.jl",
)
