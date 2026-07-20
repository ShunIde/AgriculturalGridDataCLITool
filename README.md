# Farm Report CLI

Pulls live agricultural grid data from the [Open-Meteo Forecast API](https://open-meteo.com/en/docs)
(free, no API key) and matches it against a local `farm_layout.json` to produce either:

- a daily **Markdown operations report**, or
- a queue of **mock hardware commands** (irrigation valves, frost protection)

## Install

```bash
pip install -r requirements.txt
```

## Configure your farm

`farm_report.py` is a single CLI with two subcommands: `layout` (build a
config) and `run` (fetch live data and generate outputs).

**Option A — hand-edit real coordinates.** Copy `farm_layout.example.json` to
`farm_layout.json` and edit it with your actual fields.

**Option B — pull real field polygons from OpenStreetMap.** `farm_report.py layout`
queries the free Overpass API for real farmland/orchard/vineyard polygons in a
place or bounding box, computes each field's true centroid and area, and
writes a ready-to-use `farm_layout.json`:

```bash
# By place name (geocoded via Nominatim)
python farm_report.py layout --place "Story County, Iowa" --max-fields 15

# By explicit bounding box: south,west,north,east
python farm_report.py layout --bbox 41.90,-93.75,42.05,-93.55 --max-fields 15
```

OSM farmland tagging is crowd-sourced and uneven — some regions have every
field mapped, others have none, so try a few places if you get zero results.
Crop type is rarely tagged in OSM; when missing, the script falls back to
the landuse category (`farmland` / `orchard` / `vineyard`) as a placeholder
you should edit. Verify `test_osm_offline.py` passes if you modify the
geometry math:

```bash
python test_osm_offline.py
```

Either way, each field needs:

```jsonc
{
  "id": "field-01",              // unique slug
  "name": "North Wheat Field",
  "crop": "wheat",
  "area_hectares": 12.5,
  "latitude": 41.85,
  "longitude": -93.6,
  "irrigation_zone": "zone-A",   // groups fields under one valve/controller
  "thresholds": {
    "soil_moisture_min": 0.22,   // m3/m3 — below this, irrigation is triggered
    "soil_moisture_max": 0.42,   // above this, field is flagged "overwatered"
    "frost_temp_c": 2.0          // if tomorrow's forecast low <= this, frost alert fires
  }
}
```

All three `thresholds` are optional — sane defaults are used if omitted.

## Run it

```bash
# Markdown report only (default)
python farm_report.py run --config farm_layout.json --mode report

# Mock hardware commands only
python farm_report.py run --config farm_layout.json --mode commands

# Both, into a specific folder
python farm_report.py run --config farm_layout.json --mode both --output-dir out/
```

Output files are timestamped by date:
- `farm_report_YYYY-MM-DD.md`
- `hardware_commands_YYYY-MM-DD.json`

## How it works

1. **Batched fetch** — all field coordinates are sent to Open-Meteo in a *single*
   HTTP request (Open-Meteo accepts comma-separated lat/lon lists and returns
   one result object per pair), pulling `current` conditions (soil moisture,
   soil temperature, air temp, precipitation, wind) and a 2-day `daily`
   forecast (min/max temp, precipitation sum, FAO ET0).
2. **Analysis** — each field's live soil moisture and tomorrow's forecast low
   are compared against that field's thresholds to flag irrigation need,
   over-watering, and frost risk.
3. **Output** — the same analysis feeds either the Markdown report or the
   command generator, so both outputs are always consistent with each other.

## Mock hardware commands

The `CommandGenerator` is a simulation layer that turns flags into a JSON
command queue, e.g.:

```json
{
  "device": "valve_zone-A",
  "field_id": "field-01",
  "action": "OPEN",
  "duration_minutes": 70,
  "reason": "soil_moisture 0.150 below minimum 0.22"
}
```

Irrigation duration is a simple linear function of the moisture deficit
(`DEFICIT_TO_MINUTES = 1000`, clamped to 10–120 minutes) — tune the constants
in `CommandGenerator` or replace the formula with your controller's real
scheduling logic. Swap the dict shape in `_commands_for_field()` for your
actual device protocol (MQTT topic/payload, Modbus register writes, etc.)
when wiring this to real hardware.

## Variety recommendation (phase 1 breeding demo)

`farm_report.py variety` recommends maize lines per corn field using
multi-objective genomic selection — the same Pareto-frontier idea behind
[PyBrOpS](https://github.com/rzshrote/pybrops) / the paper this project was
inspired by ([Shrote & Beavis, *G3*, 2024](https://doi.org/10.1093/g3journal/jkae199)).

Two interchangeable engines implement this, selected with `--engine`:

- **`pymoo`** (`breeding.py`) — a from-scratch reimplementation on top of
  [`pymoo`](https://pymoo.org/) directly, parsing the VCF with the stdlib
  `gzip` module. No `cyvcf2` dependency, so it runs natively on Windows.
- **`pybrops`** (`breeding_pybrops.py`) — the real `pybrops` package,
  including its actual `OptimalContributionSubsetSelection` protocol and
  `cyvcf2`-backed VCF loading. Linux/macOS only, since `cyvcf2` ships no
  Windows wheel.

`--engine auto` (the default) picks `pybrops` if it's installed, otherwise
falls back to `pymoo` — so the same command works out of the box on either
platform.

```bash
# pymoo engine (Windows-safe, no cyvcf2)
pip install -r requirements-breeding.txt
python farm_report.py variety --config farm_layout.json --output-dir out/

# pybrops engine (Linux/macOS) — install into a venv, see caveats below
python3 -m venv .venv && source .venv/bin/activate
pip install -c <(echo "numpy<2") -r requirements-breeding-pybrops.txt
python farm_report.py variety --config farm_layout.json --output-dir out/ --engine pybrops
```

**pybrops install caveat:** `pybrops` 1.0.3 (latest on PyPI) is incompatible
with `numpy>=2.0` (it uses `np.float_`, removed in NumPy 2.0), and a plain
`pip install pybrops` pulls in a `cyvcf2` version that requires
`numpy>=2.0.0`. `requirements-breeding-pybrops.txt` pins `numpy<2` and
`cyvcf2<0.32` (which only needs `numpy>=1.16`) so the install resolves
cleanly — always install with that constraint applied, in a dedicated venv.

**What it does (both engines):**
1. Loads real genotype data — 942 lines x 2000 SNPs from the Wisconsin
   Diversity (WiDiv) maize panel, vendored from PyBrOpS's own examples
   (see `data/ATTRIBUTION.md`).
2. Synthesizes two negatively-correlated trait effects, `cold_tolerance`
   and `drought_tolerance`, over those markers — mirroring PyBrOpS's own
   demo almost exactly. **These are simulated trait effects, not measured
   phenotypes** — the panel has no public cold/drought phenotyping, so this
   demonstrates the selection methodology, not a real agronomic call.
3. Runs an NSGA-II subset-selection genetic algorithm to find the Pareto
   frontier trading off mean breeding value per trait against within-subset
   inbreeding (a genomic relationship / kinship matrix computed via VanRaden
   method 1) — conceptually the same problem Optimal Contribution Selection
   solves. With `--engine pymoo` this is a hand-written `pymoo` problem; with
   `--engine pybrops` it's the real `OptimalContributionSubsetSelection`
   protocol (which itself wraps `pymoo`'s NSGA-II under the hood). This is
   the expensive step, so the frontier is cached (`data/ocs_frontier_cache.npz`
   or `data/ocs_frontier_cache_pybrops.npz`, per engine) and only recomputed
   when its parameters change (or `--regen-frontier` is passed).
4. For each corn field, converts its `frost_temp_c` / `soil_moisture_min`
   thresholds into a `(cold_weight, drought_weight)` preference — a lower
   threshold means the grower leans more on the crop's innate tolerance
   before intervening — and picks the frontier point that best matches it,
   reporting the real WiDiv line IDs in that selection.

Output: `variety_recommendations_YYYY-MM-DD.md` in `--output-dir`, one
section per corn field with its trait weights and top candidate lines.

Only `crop: "corn"` fields are handled — other crops are logged and skipped,
since the maize panel is the only genotype data currently vendored.

## Multi-generation breeding simulation (phase 2, pybrops only)

`farm_report.py breeding-sim` goes a step further than `variety`'s
single-generation frontier: it actually simulates several generations of a
breeding program using real `pybrops` machinery.

```bash
python farm_report.py breeding-sim --config farm_layout.json --generations 5 --output-dir out/
```

Each generation: runs the same `OptimalContributionSubsetSelection` used by
`variety` (this time selecting biparental crosses, not just a flat parent
list), picks the best cross set for each corn field's trait preference,
mates them with pybrops's `TwoWayCross` mating protocol, and re-evaluates
GEBVs/kinship on the resulting progeny before repeating. Output
(`breeding_simulation_YYYY-MM-DD.md`) is a per-field table of mean breeding
value and kinship across generations — showing whether sustained selection
is actually improving the trait(s) that field cares about, and how fast it's
burning through genetic diversity to get there.

**Known approximation — no real genetic map:** simulating meiosis (which
markers recombine together vs. independently) requires a genetic map, and
no public genetic map exists for this specific 2000-SNP WiDiv panel. This
command fabricates one at a flat 1 cM/Mb genome-wide average (a commonly
cited maize-wide rule of thumb, not a measurement of this panel's actual
recombination landscape) via pybrops's `StandardGeneticMap` +
`HaldaneMapFunction`. Generation-to-generation trends here demonstrate the
simulation methodology, not a validated breeding forecast — see
`breeding_pybrops.py`'s module docstring for detail. This only affects
`breeding-sim`; `variety`'s single-generation frontier selects among
*existing* lines and never simulates recombination, so it's unaffected.

`breeding-sim` requires the `pybrops` engine (no `pymoo` equivalent exists,
since `pymoo` alone has no genotype/mating-simulation machinery) —
Linux/macOS only, see the install caveat above.

## Testing without network access

`test_offline.py` monkeypatches `GridDataClient.fetch` with synthetic data
covering all three branches (low moisture, over-watered, frost risk) so you
can validate the report/command logic without hitting the live API:

```bash
python test_offline.py
```

## Notes

- Open-Meteo's free tier is intended for non-commercial use; check their
  [terms](https://open-meteo.com/en/terms) if you plan to run this
  operationally at scale, and consider their commercial API tier if needed.
- `soil_moisture_0_to_1cm` is a shallow-layer estimate from Open-Meteo's
  underlying weather models (not a field sensor reading) — treat it as a
  directional signal, and calibrate thresholds against real soil probes if
  you have them.
