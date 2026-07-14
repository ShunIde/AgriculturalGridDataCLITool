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
