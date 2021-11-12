project(
    'vortsf',
    'fortran',
    default_options: [
        'buildtype=debugoptimized',
        'default_library=shared',
        'warning_level=2',
    ]
)

srcs = [
    'src/m_params.f90',
    'src/m_io.f90',
    'src/m_vorts.f90',
    'src/capi.f90',
]

vorts_lib = library(
    meson.project_name(),
    sources: srcs,
    install: false,
)
