# Operational MPAS-Atmosphere 8.0 with Dust

Operational MPAS-Atmosphere forecasting workflows, source-code
modifications, HPE Cray build configuration, dust/chemistry support,
initialization, post-processing and visualization utilities.

## Overview

This repository documents an operational MPAS-Atmosphere 8.0
forecasting system used in a high-performance computing environment.

The published system includes:

- MPAS-Atmosphere 8.0 source code
- dust-enabled chemistry support
- operational MPAS initialization
- MPAS forecast execution
- lateral-boundary and preprocessing workflows
- UPP post-processing
- MPAS-to-WRF conversion
- ASCII forecast products
- statistical utilities
- GeoJSON / polygon / isoband production
- vector-tile and MBTiles generation
- HPE Cray / Slurm operational configuration

## Operational Workflow

The operational forecast chain is approximately:

GFS / atmospheric input data
        |
        v
     ungrib
        |
        v
     metgrid
        |
        v
MPAS init_atmosphere
        |
        v
 MPAS-Atmosphere
   + dust/chemistry
        |
        +----------------------+
        |                      |
        v                      v
       UPP                 MPAS-to-WRF
        |                      |
        v                      v
 GRIB products            WRF-format fields
        |
        +----------------------+
        |
        v
 ASCII / statistics
        |
        v
 GeoJSON / isobands / polygons
        |
        v
 vector tiles / MBTiles

## Dust-Enabled MPAS

The dust capability in this MPAS-Atmosphere 8.0 branch builds on the
scientific implementation introduced for MPAS 7.0 and described in:

https://doi.org/10.1029/2023MS003636

Because the MPAS source-code structure changed substantially between
versions 7.0 and 8.0, the implementation required extensive
redevelopment rather than a direct source-code port.

Parts of the earlier dust implementation were retained where
appropriate, while the integration with the MPAS 8.0 atmosphere
core, chemistry infrastructure, physics, Registry, time integration
and initialization system was reworked.

A major addition in this implementation is native handling of the
erodibility field (`erod`) in `init_atmosphere`. Static erodibility is
interpolated onto the MPAS mesh during initialization and subsequently
used directly by the GOCART dust-emission pathway.

The operational implementation also includes dust initialization,
aerosol dry deposition and aerosol-radiation feedback.

See [DUST_IMPLEMENTATION.md](DUST_IMPLEMENTATION.md) for the detailed
technical description.

## HPE Cray Build Environment

The operational source has been built using an HPE Cray programming
environment including:

- PrgEnv-cray
- Cray MPICH
- Cray LibSci
- Cray NetCDF
- Cray Parallel-NetCDF
- Libfabric / OFI

The exact operational build script is provided with the source tree.

## Repository Layout

source/
  mpas-8.0-dust/
      MPAS 8.0 source tree with operational dust modifications

config/
  mpas-8.0/
      operational MPAS namelists and streams

operational/
      initialization, preprocessing, forecast, UPP,
      MPAS-to-WRF and ASCII/statistical workflows

postprocessing/
  json/
      GeoJSON, isoband, polygon and MBTiles workflows

docs/
      additional technical documentation

## Operational Configuration

The active operational MPAS configuration includes:

- `namelist.init_atmosphere`
- `streams.init_atmosphere`
- `namelist.atmosphere`
- `streams.atmosphere`
- atmospheric output stream lists
- diagnostic stream lists

Site-specific filesystem paths and Slurm settings reflect the
original HPE Cray environment and must be adapted for other systems.

## Runtime Physics Data

Large MPAS/WRF physics lookup datasets are intentionally excluded
from this source-code repository.

These include large Thompson microphysics, RRTMG, CAM and ozone
lookup datasets. They must be obtained from a compatible MPAS/WRF
runtime-data distribution.

## Operational Scripts

Important operational scripts include:

- `start.mpas`
- `restart.mpas`
- `modelenv.sh`
- `ungrib.mpas.sh`
- `metgrid.mpas.sh`
- `init.mpas.sh`
- `run.mpas.sh`
- `upp.mpas.sh`
- `mpas2wrf.sh`
- `ascii.mpas.sh`
- `prepare_ascii.mpas.sh`
- `stats.mpas.sh`
- `toctl.mpas.sh`

## Visualization Pipeline

The repository also includes the MPAS post-processing workflow used
to generate web-oriented meteorological products.

This includes:

- GeoJSON products
- filled contours / isobands
- polygon products
- vector tiles
- MBTiles

The scripts are organized under:

    postprocessing/json/

## Site-Specific Information

The scripts are published substantially as used operationally.

Private network addresses and internal infrastructure endpoints have
been removed from the public version.

Paths under `/scratch/lus/...`, HPE Cray module names and Slurm
settings are retained because they document the actual tested
operational environment.

## Author

**Radi Ajjaji**

Numerical Weather Prediction | Data Assimilation | HPC | AI Weather Models

LinkedIn:
https://www.linkedin.com/in/radi-ajjaji-10008071
