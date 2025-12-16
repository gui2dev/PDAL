(filters.mta)=

# filters.mta

The **MTA (Multiple-Time-Around) filter** performs zone analysis and correction for multi-echo LiDAR data
acquired during rapid scanner motion. When a scanner (rotating or oscillating) operates at high speed, 
echoes may be received after the scanner mechanism has returned to a similar orientation, creating range ambiguity.

```{eval-rst}
.. streamable::
```

```{note}
This filter expects:

- Point coordinates (`X`, `Y`, `Z`) in ECEF (Earth-Centered, Earth-Fixed) format;
- Beam geometry dimensions (`BeamOriginX/Y/Z`, `BeamDirectionX/Y/Z`) for `auto` mode;
- Shot metadata: `ShotTimestamp`, `EchoRange`, `UnambiguousRange`, `NumberOfReturns`;
- Round detection via `Segment` dimension (preferred) or `EdgeOfFlightLine` (fallback).

Use `filters.georeference` with `transform_beam: true` before this filter.
```

```{note}
**Point Ordering**: Points are expected to be pre-ordered by beam and timestamp.
The filter does not perform internal sorting; it processes points sequentially and groups them
by beam matching (origin/direction proximity) within buffered time windows.
Ensure your input pipeline preserves or establishes this ordering.
```

## Examples

### Basic pipeline: georeference + MTA correction (auto mode)

```json
[
    "input.rxp",
    {
        "type": "filters.georeference",
        "trajectory_file": "trajectory.out",
        "trajectory_options": {
            "type": "readers.sbet"
        },
        "scan2imu": "1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1",
        "coordinate_system": "NED",
        "time_offset": 0.0,
        "transform_beam": true
    },
    {
        "type": "filters.mta",
        "mode": "auto",
        "max_zones": 4
    },
    "output.las"
]
```



## Options

mode

: Zone selection strategy. One of: `auto`, `fixed`, `semi-auto`.
  - `auto`: Multi-round automatic correction using geometry and/or reflectance coherence.
  - `fixed`: Apply a user-specified zone to all points.
  - `semi-auto`: Select zone to place corrected range within `[gate_near, gate_far]` window.
  \[Default: `auto`\]

fixed_zone

: Zone number (1-based) for `fixed` mode. \[Default: 1\]

gate_near

: Near boundary (meters) for `semi-auto` mode. \[Default: 0.0\]

gate_far

: Far boundary (meters) for `semi-auto` mode. \[Default: 0.0\]

jump_tolerance

: Tolerance (meters) around multiples of `UnambiguousRange` for acceptable jumps between consecutive shots
  when running `auto` (Viterbi) mode. \[Default: 0.05\]

max_zones

: Maximum zone candidates to evaluate (typically 2–4 for most scanners). \[Default: 4\]

unambiguous_range

: Unambiguous range override in meters. If set, this value is used for all points,
  overriding the `UnambiguousRange` dimension. Otherwise, reads from dimension.
  \[Default: 0.0 (use dimension)\]

```{include} filter_opts.md
```

## See Also

- [filters.georeference](filters.georeference.md) — Prerequisite for ECEF and beam transformation
- [readers.rxp](readers.rxp.md) — RXP reader with shot metadata and beam geometry
