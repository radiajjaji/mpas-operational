# Operational MPAS 8.0 Configuration

This directory contains the active namelist and stream configuration
used by the operational dust-enabled MPAS-Atmosphere 8.0 workflow.

Included files:

- `namelist.init_atmosphere`
- `streams.init_atmosphere`
- `namelist.atmosphere`
- `streams.atmosphere`
- `stream_list.atmosphere.surface`
- `stream_list.atmosphere.diagnostics`
- `stream_list.atmosphere.output`

The configuration reflects the original operational HPE Cray
environment and may contain site-specific paths, stream names,
resource settings and model options that must be adapted on other
systems.

Historical or superseded namelist copies are intentionally omitted
from this repository.
