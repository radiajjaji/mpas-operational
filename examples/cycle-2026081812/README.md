# Operational MPAS Forecast Example — 2026081812

This directory contains representative logs from a real operational MPAS
forecast cycle.

The objective is to provide a compact reference showing the major stages of
the operational model workflow without storing every parallel post-processing
task log.

## Included logs

### `init.mpas.out`

Log from the MPAS initialization stage.

This stage prepares the native MPAS initial atmospheric state required by the
forecast model.

### `mpas.out`

Main operational MPAS-Atmosphere forecast log.

This is the principal model execution log and contains the model
initialization, configuration, physics initialization, time stepping,
diagnostics and normal forecast completion information.

### `mpas2wrf.189.out`

One representative log from the parallel MPAS-to-WRF conversion stage.

The operational workflow creates many such tasks corresponding to forecast
output times. Only one representative member is retained here because the
individual task logs have the same basic processing structure.

### `upp.mpas.204.out`

One representative Unified Post Processor execution log.

The operational workflow runs UPP in parallel for many forecast times. A
single representative log is retained as an example of the post-processing
stage.

## Operational chain

The example corresponds conceptually to:

    initial/boundary data
            |
            v
       MPAS initialization
            |
            | init.mpas.out
            v
       MPAS forecast
            |
            | mpas.out
            v
       native MPAS outputs
            |
            +----------------------+
            |                      |
            v                      v
        MPAS -> WRF              UPP
            |                      |
    mpas2wrf.189.out       upp.mpas.204.out
            |                      |
            +----------+-----------+
                       |
                       v
                downstream products

## Why only representative post-processing logs are stored

The operational forecast launches many parallel MPAS-to-WRF and UPP jobs,
usually one for each forecast valid time.

Those logs are highly repetitive.

To keep the repository compact while still providing reproducible operational
examples, this directory contains:

- the full initialization log;
- the full model forecast log;
- one representative MPAS-to-WRF log;
- one representative UPP log.

The complete operational scripts themselves are stored elsewhere in this
repository under `operational/`.

## Important note

The absolute paths, host names and environment information in these logs are
from the operational HPE Cray system on which the cycle was executed.

They are preserved intentionally because these are real operational logs.
Users reproducing the workflow on another platform must adapt those paths and
resource settings to their own environment.
