"""
Provides maps between MMSchema and OpenFF potential types.
TODO: make this more rigorous and expand the maps.
"""

_dihedrals_potentials_map = {
    "k*(1+cos(periodicity*theta-phase))": "CharmmMulti",
    # need to add all supported potentials in OpenFFTk
}

_dihedrals_improper_potentials_map = {
    "k*(1+cos(periodicity*theta-phase))": "CharmmMulti",
    # need to add all supported potentials in OpenFFTk
}
