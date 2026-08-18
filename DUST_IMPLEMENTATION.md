# Dust Implementation in Operational MPAS-Atmosphere 8.0

## Overview

This repository contains a dust-enabled MPAS-Atmosphere 8.0
implementation integrated into an operational numerical weather
prediction workflow.

The implementation includes dust initialization, static erodibility,
GOCART dust emissions, aerosol dry deposition, aerosol optical
properties, aerosol-radiation feedback and the associated MPAS
chemistry infrastructure.

## Scientific provenance

The dust capability builds scientifically on the MPAS dust
implementation introduced for MPAS-Atmosphere version 7.0 and
described in:

https://doi.org/10.1029/2023MS003636

That work provided an important scientific and software foundation
for dust modelling in MPAS.

The implementation in this repository is, however, not a direct copy
or simple port of the MPAS 7.0 implementation.

## Redevelopment for MPAS 8.0

The MPAS software structure changed substantially between versions
7.0 and 8.0.

Consequently, the dust and chemistry capability required extensive
redevelopment for integration into the MPAS 8.0 code base.

Selected elements of the MPAS 7.0 dust implementation were retained
where they remained applicable, while substantial portions of the
integration were reworked for MPAS 8.0.

The resulting changes span several parts of the model, including:

- atmosphere-core integration;
- chemistry infrastructure;
- atmospheric time integration;
- physics interfaces;
- microphysics interfaces;
- longwave and shortwave radiation interfaces;
- MPAS Registry definitions;
- init_atmosphere;
- static geographical-field processing;
- build-system integration;
- operational namelist and stream configuration.

A comparison with the official MPAS v8.0.0 source shows that the
operational branch adds an entire:

    src/core_atmosphere/chemistry/

subsystem that is not present in the pristine MPAS 8.0 source tree.

## Erodibility in init_atmosphere

A major addition in this MPAS 8.0 implementation is native handling
of the dust erodibility field:

    erod

inside `init_atmosphere`.

The static initialization code interpolates the geographical
erodibility dataset onto the MPAS mesh and stores the resulting
erodible-surface fractions in the initialized MPAS state.

The atmospheric model subsequently retrieves this native MPAS
`erod` field and supplies it directly to the GOCART dust-emission
calculation.

This removes the need to introduce the erodibility field externally
after MPAS initialization.

## Dust initialization

The initialization core explicitly handles the dust first-guess
fields:

    DUST_1
    DUST_2
    DUST_3
    DUST_4
    DUST_5

which are mapped into the MPAS chemistry species:

    dust_a01
    dust_a02
    dust_a03
    dust_a04
    dust_a05

The initialization logic vertically interpolates these fields to the
MPAS model levels and prevents negative initialized dust
concentrations.

## GOCART dust representation

The GOCART dust-emission implementation contains five dust size bins.

The source defines effective radii of approximately:

    0.73 micrometers
    1.40 micrometers
    2.40 micrometers
    4.50 micrometers
    8.00 micrometers

with particle densities of 2500 kg m-3 for the first bin and
2650 kg m-3 for the remaining bins.

The MPAS chemistry package selected operationally is:

    config_chem_opt = 'bin4aer_dust'

This package name should not be confused with the number of GOCART
dust-emission size bins; the GOCART implementation itself uses five
dust bins.

## Operational configuration

The current operational MPAS chemistry and dust configuration includes:

    config_use_aer_IN         = false
    config_gas_drydep_opt     = 1
    config_cu_scav_on         = true

    config_chem_opt           = 'bin4aer_dust'
    config_aer_op_opt         = 'volume_approx'

    config_dust_opt           = 'gocart'
    config_dust_factor        = 0.55

    config_drydep_on          = true
    config_wetscav_on         = true
    config_aer_drydep_opt     = 201
    config_kdepvel            = 1

    config_aer_ra_feedback    = true

The principal options have the following operational meaning:

- `config_chem_opt = 'bin4aer_dust'` selects the MPAS chemistry
  package used by the operational dust configuration.

- `config_dust_opt = 'gocart'` activates the GOCART dust-emission
  formulation.

- `config_dust_factor = 0.55` is the operational multiplicative
  tuning factor applied to GOCART dust emissions.

- `config_drydep_on = true` enables aerosol dry deposition.

- `config_aer_drydep_opt = 201` selects the operational aerosol
  dry-deposition treatment.

- `config_kdepvel = 1` selects the configured deposition-velocity
  treatment used by the aerosol dry-deposition pathway.

- `config_wetscav_on = true` enables wet scavenging of aerosol
  and dust.

- `config_cu_scav_on = true` enables aerosol scavenging associated
  with convective processes.

- `config_aer_op_opt = 'volume_approx'` selects the aerosol optical
  property treatment.

- `config_aer_ra_feedback = true` enables aerosol-radiation
  feedback.

- `config_use_aer_IN = false` disables aerosol-aware ice nucleation
  in the current operational configuration.

- `config_gas_drydep_opt = 1` selects the gas dry-deposition option
  used by the chemistry framework.

The value of `config_dust_factor` is an operational tuning factor
applied to GOCART dust emissions.

## Emission pathway

The operational pathway is conceptually:

    static erodibility data
             |
             v
      init_atmosphere
             |
             v
          erod
             |
             v
    MPAS atmospheric state
             |
             v
    chemistry emission driver
             |
             v
      GOCART dust emission
             |
             v
       five dust bins

The emission calculation combines erodible-surface fraction with
meteorological and land-surface information and subsequently applies
the configured dust-emission scaling factor.

## Deposition

Aerosol dry deposition is enabled operationally.

The MPAS chemistry framework includes aerosol deposition-velocity
calculations involving gravitational settling, aerosol diffusivity,
Schmidt number, aerodynamic resistance and friction velocity.

## Aerosol-radiation interaction

Aerosol-radiation feedback is enabled in the operational
configuration.

The MPAS Registry and atmospheric physics interfaces contain the
shortwave and longwave aerosol optical quantities required by the
dust/aerosol radiation pathway.

## Build

The model is compiled with chemistry support using:

    CPP_EXTRA_FLAGS="-D DO_CHEMISTRY"

for the relevant MPAS cores.

The supplied `comp.sh` documents the HPE Cray build environment used
for the operational implementation.

## Runtime configuration

The active MPAS namelist and stream configuration is provided under:

    config/mpas-8.0/

Site-specific paths and resource settings reflect the original
operational HPE Cray environment and must be adapted for other
systems.

## Attribution

The MPAS v7 dust work described in:

https://doi.org/10.1029/2023MS003636

should be cited when using or discussing the scientific foundation
of the dust capability.

The MPAS 8.0 source in this repository represents a substantial
redevelopment and operational integration of that capability for the
newer MPAS code structure, including additional initialization and
workflow functionality described above.
