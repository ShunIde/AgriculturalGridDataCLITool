#!/usr/bin/env python3
"""
farm_report.py — A single CLI with three subcommands:

  layout   Generate a real farm_layout.json from OpenStreetMap farmland
           data (via Overpass + Nominatim), given a place name or bbox.

  run      Pull live agricultural grid data from Open-Meteo, match it
           against farm_layout.json, and produce a daily Markdown
           operations report and/or a set of mock hardware commands.

  variety  Recommend maize lines per corn field via climate-informed
           multi-objective genomic selection over a real SNP panel. Uses the
           real `pybrops` + `cyvcf2` engine when available (Linux/macOS —
           see breeding_pybrops.py), falling back to a from-scratch
           `pymoo` reimplementation otherwise (Windows-safe — see
           breeding.py). Override with `--engine {pymoo,pybrops}`.

  breeding-sim  Simulate several generations of OCS-driven selection +
           mating per corn field using the real pybrops engine (Linux/macOS
           only — see breeding_pybrops.py).

Usage:
    # 1) Generate a config from real field data
    python farm_report.py layout --place "Story County, Iowa" --max-fields 15

    # 2) Fetch live weather/soil data and generate outputs
    python farm_report.py run --config farm_layout.json --mode both --output-dir out/

    # 3) Recommend maize varieties per field
    python farm_report.py variety --config farm_layout.json --output-dir out/

    # 4) Simulate multiple generations of selection + mating (pybrops only)
    python farm_report.py breeding-sim --config farm_layout.json --generations 5

Data sources (all free, no API key required):
    - Nominatim (geocoding a place name -> bounding box)
    - Overpass API (real farmland/orchard/vineyard polygons)
    - Open-Meteo Forecast API (current + forecast soil/weather grid data)
    - Wisconsin Diversity (WiDiv) maize SNP panel, vendored from PyBrOpS's
      own examples (MIT-licensed) for the `variety` subcommand
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

import requests

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("farm_report")


# ==========================================================================
# Subcommand: layout  (OpenStreetMap -> farm_layout.json)
# ==========================================================================

NOMINATIM_URL = "https://nominatim.openstreetmap.org/search"
OVERPASS_URL = "https://overpass-api.de/api/interpreter"
OSM_USER_AGENT = "farm-report-cli/1.0 (educational/demo use)"

LANDUSE_TAGS = ["farmland", "orchard", "vineyard", "greenhouse_horticulture"]

# Rough default thresholds per landuse category, used when generating a
# config from OSM — edit these in the output file to match your agronomy.
DEFAULT_THRESHOLDS = {
    "farmland":  {"soil_moisture_min": 0.20, "soil_moisture_max": 0.42, "frost_temp_c": 2.0},
    "orchard":   {"soil_moisture_min": 0.18, "soil_moisture_max": 0.38, "frost_temp_c": 1.0},
    "vineyard":  {"soil_moisture_min": 0.15, "soil_moisture_max": 0.35, "frost_temp_c": 0.0},
    "greenhouse_horticulture": {"soil_moisture_min": 0.25, "soil_moisture_max": 0.45, "frost_temp_c": 3.0},
}


def geocode_place(place: str, session: requests.Session, timeout: int = 20) -> tuple[float, float, float, float]:
    """Returns (south, west, north, east) for a place name via Nominatim."""
    log.info("Geocoding %r via Nominatim...", place)
    resp = session.get(
        NOMINATIM_URL,
        params={"q": place, "format": "json", "limit": 1},
        headers={"User-Agent": OSM_USER_AGENT},
        timeout=timeout,
    )
    resp.raise_for_status()
    results = resp.json()
    if not results:
        raise ValueError(f"Nominatim found no results for place: {place!r}")

    bbox = results[0]["boundingbox"]  # [south, north, west, east] as strings
    south, north, west, east = (float(v) for v in bbox)
    log.info("Bounding box: south=%.4f west=%.4f north=%.4f east=%.4f", south, west, north, east)
    return south, west, north, east


def build_overpass_query(bbox: tuple[float, float, float, float]) -> str:
    south, west, north, east = bbox
    bbox_str = f"{south},{west},{north},{east}"
    tag_filter = "|".join(LANDUSE_TAGS)
    return f"""
[out:json][timeout:60];
(
  way["landuse"~"^({tag_filter})$"]({bbox_str});
);
out geom;
""".strip()


def fetch_overpass(query: str, session: requests.Session, timeout: int = 90) -> list[dict]:
    log.info("Querying Overpass API for farmland polygons...")
    resp = session.post(OVERPASS_URL, data={"data": query}, timeout=timeout)
    resp.raise_for_status()
    payload = resp.json()
    elements = payload.get("elements", [])
    log.info("Overpass returned %d element(s).", len(elements))
    return elements


def polygon_centroid_and_area_ha(points_latlon: list[tuple[float, float]]) -> tuple[float, float, float]:
    """
    points_latlon: list of (lat, lon) forming a closed or near-closed ring.
    Returns (centroid_lat, centroid_lon, area_hectares) using a local
    equirectangular projection centered on the ring's mean latitude —
    accurate enough for field-sized polygons.
    """
    if points_latlon[0] != points_latlon[-1]:
        points_latlon = points_latlon + [points_latlon[0]]  # close the ring

    lat0 = sum(p[0] for p in points_latlon) / len(points_latlon)
    lon0 = sum(p[1] for p in points_latlon) / len(points_latlon)
    lat0_rad = math.radians(lat0)

    m_per_deg_lat = 110_540.0
    m_per_deg_lon = 111_320.0 * math.cos(lat0_rad)

    xy = [
        ((lon - lon0) * m_per_deg_lon, (lat - lat0) * m_per_deg_lat)
        for lat, lon in points_latlon
    ]

    a_sum = 0.0
    cx_sum = 0.0
    cy_sum = 0.0
    for (x0, y0), (x1, y1) in zip(xy, xy[1:]):
        cross = x0 * y1 - x1 * y0
        a_sum += cross
        cx_sum += (x0 + x1) * cross
        cy_sum += (y0 + y1) * cross

    area_m2 = abs(a_sum) / 2.0
    if area_m2 < 1e-6:
        cx, cy = xy[0][0], xy[0][1]
    else:
        cx = cx_sum / (3.0 * a_sum)
        cy = cy_sum / (3.0 * a_sum)

    centroid_lat = lat0 + cy / m_per_deg_lat
    centroid_lon = lon0 + cx / m_per_deg_lon
    area_ha = area_m2 / 10_000.0

    return centroid_lat, centroid_lon, area_ha


def elements_to_fields(elements: list[dict], max_fields: int, min_area_ha: float) -> list[dict[str, Any]]:
    fields = []
    counter = 0
    for el in elements:
        if el.get("type") != "way":
            continue
        geometry = el.get("geometry")
        if not geometry or len(geometry) < 3:
            continue

        points = [(pt["lat"], pt["lon"]) for pt in geometry]
        try:
            lat, lon, area_ha = polygon_centroid_and_area_ha(points)
        except Exception as e:
            log.debug("Skipping way %s: geometry error %s", el.get("id"), e)
            continue

        if area_ha < min_area_ha:
            continue

        tags = el.get("tags", {})
        landuse = tags.get("landuse", "farmland")
        crop = tags.get("crop", landuse)
        name = tags.get("name") or f"{landuse.title()} Field {counter + 1:02d}"

        thresholds = DEFAULT_THRESHOLDS.get(landuse, DEFAULT_THRESHOLDS["farmland"])

        counter += 1
        fields.append(
            {
                "id": f"field-{counter:02d}",
                "name": name,
                "crop": crop,
                "area_hectares": round(area_ha, 2),
                "latitude": round(lat, 6),
                "longitude": round(lon, 6),
                "irrigation_zone": f"zone-{chr(ord('A') + (counter - 1) % 6)}",
                "thresholds": thresholds,
                "_osm_way_id": el.get("id"),  # traceability; harmless extra key
            }
        )

        if len(fields) >= max_fields:
            break

    return fields


def cmd_layout(args: argparse.Namespace) -> int:
    session = requests.Session()

    try:
        if args.place:
            bbox = geocode_place(args.place, session)
            farm_name = args.farm_name or args.place
        else:
            south, west, north, east = (float(v) for v in args.bbox.split(","))
            bbox = (south, west, north, east)
            farm_name = args.farm_name or "Custom Region"
    except (requests.RequestException, ValueError) as e:
        log.error("Failed to resolve location: %s", e)
        return 1

    query = build_overpass_query(bbox)
    try:
        elements = fetch_overpass(query, session)
    except requests.RequestException as e:
        log.error("Failed to query Overpass API: %s", e)
        return 1

    fields = elements_to_fields(elements, args.max_fields, args.min_area_ha)
    if not fields:
        log.warning(
            "No farmland polygons found (or all were below --min-area-ha). "
            "Try a different place, a larger bounding box, or lower --min-area-ha."
        )
        return 1

    output = {"farm_name": farm_name, "fields": fields}
    args.output.write_text(json.dumps(output, indent=2), encoding="utf-8")
    log.info("Wrote %d field(s) to %s", len(fields), args.output)
    return 0


# ==========================================================================
# Subcommand: run  (farm_layout.json + Open-Meteo -> report / commands)
# ==========================================================================

OPEN_METEO_URL = "https://api.open-meteo.com/v1/forecast"

CURRENT_VARS = [
    "temperature_2m",
    "relative_humidity_2m",
    "precipitation",
    "wind_speed_10m",
    "soil_temperature_0cm",
    "soil_moisture_0_to_1cm",
]

DAILY_VARS = [
    "temperature_2m_max",
    "temperature_2m_min",
    "precipitation_sum",
    "et0_fao_evapotranspiration",
]


@dataclass
class Field:
    id: str
    name: str
    crop: str
    area_hectares: float
    latitude: float
    longitude: float
    irrigation_zone: str
    thresholds: dict[str, float]

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Field":
        required = ["id", "name", "crop", "latitude", "longitude"]
        missing = [k for k in required if k not in d]
        if missing:
            raise ValueError(f"Field entry missing required keys {missing}: {d}")

        thresholds_in = d.get("thresholds", {})
        return cls(
            id=str(d["id"]),
            name=str(d["name"]),
            crop=str(d["crop"]),
            area_hectares=float(d.get("area_hectares", 0.0)),
            latitude=float(d["latitude"]),
            longitude=float(d["longitude"]),
            irrigation_zone=str(d.get("irrigation_zone", "default")),
            thresholds={
                "soil_moisture_min": float(thresholds_in.get("soil_moisture_min", 0.20)),
                "soil_moisture_max": float(thresholds_in.get("soil_moisture_max", 0.45)),
                "frost_temp_c": float(thresholds_in.get("frost_temp_c", 2.0)),
            },
        )


class FarmLayout:
    """Loads and validates farm_layout.json."""

    def __init__(self, path: Path):
        self.path = path
        raw = json.loads(path.read_text(encoding="utf-8"))
        self.farm_name: str = raw.get("farm_name", "Unnamed Farm")
        field_defs = raw.get("fields", [])
        if not field_defs:
            raise ValueError("farm_layout.json contains no 'fields' entries")
        self.fields: list[Field] = [Field.from_dict(f) for f in field_defs]


class GridDataClient:
    """
    Fetches current conditions + a 2-day forecast for every field in a
    single batched request. Open-Meteo accepts comma-separated lists of
    latitude/longitude and returns one result object per coordinate pair
    (in the same order), so N fields costs exactly one HTTP call.
    """

    def __init__(self, timeout: int = 20, session: requests.Session | None = None):
        self.timeout = timeout
        self.session = session or requests.Session()

    def fetch(self, fields: list[Field]) -> dict[str, dict]:
        params = {
            "latitude": ",".join(str(f.latitude) for f in fields),
            "longitude": ",".join(str(f.longitude) for f in fields),
            "current": ",".join(CURRENT_VARS),
            "daily": ",".join(DAILY_VARS),
            "forecast_days": 2,
            "timezone": "auto",
        }
        log.info("Requesting grid data for %d field(s) from Open-Meteo...", len(fields))
        resp = self.session.get(OPEN_METEO_URL, params=params, timeout=self.timeout)
        resp.raise_for_status()
        payload = resp.json()

        if isinstance(payload, dict):
            payload = [payload]

        if len(payload) != len(fields):
            raise RuntimeError(
                f"Expected {len(fields)} results from Open-Meteo, got {len(payload)}"
            )

        return {f.id: result for f, result in zip(fields, payload)}


@dataclass
class FieldAnalysis:
    field: Field
    air_temp_c: float | None
    humidity_pct: float | None
    wind_kmh: float | None
    precip_now_mm: float | None
    soil_temp_c: float | None
    soil_moisture: float | None
    temp_min_tomorrow_c: float | None
    temp_max_today_c: float | None
    precip_sum_today_mm: float | None
    et0_today_mm: float | None
    frost_risk: bool
    irrigation_needed: bool
    irrigation_deficit: float | None
    overwatered: bool


def _safe_get(lst: list | None, idx: int) -> Any:
    if not lst or idx >= len(lst):
        return None
    return lst[idx]


def analyze_field(field: Field, raw: dict) -> FieldAnalysis:
    current = raw.get("current", {})
    daily = raw.get("daily", {})

    soil_moisture = current.get("soil_moisture_0_to_1cm")
    soil_temp = current.get("soil_temperature_0cm")

    temp_min_tomorrow = _safe_get(daily.get("temperature_2m_min"), 1)
    temp_max_today = _safe_get(daily.get("temperature_2m_max"), 0)
    precip_sum_today = _safe_get(daily.get("precipitation_sum"), 0)
    et0_today = _safe_get(daily.get("et0_fao_evapotranspiration"), 0)

    frost_risk = (
        temp_min_tomorrow is not None
        and temp_min_tomorrow <= field.thresholds["frost_temp_c"]
    )

    irrigation_needed = False
    irrigation_deficit = None
    overwatered = False
    if soil_moisture is not None:
        if soil_moisture < field.thresholds["soil_moisture_min"]:
            irrigation_needed = True
            irrigation_deficit = round(field.thresholds["soil_moisture_min"] - soil_moisture, 4)
        elif soil_moisture > field.thresholds["soil_moisture_max"]:
            overwatered = True

    return FieldAnalysis(
        field=field,
        air_temp_c=current.get("temperature_2m"),
        humidity_pct=current.get("relative_humidity_2m"),
        wind_kmh=current.get("wind_speed_10m"),
        precip_now_mm=current.get("precipitation"),
        soil_temp_c=soil_temp,
        soil_moisture=soil_moisture,
        temp_min_tomorrow_c=temp_min_tomorrow,
        temp_max_today_c=temp_max_today,
        precip_sum_today_mm=precip_sum_today,
        et0_today_mm=et0_today,
        frost_risk=frost_risk,
        irrigation_needed=irrigation_needed,
        irrigation_deficit=irrigation_deficit,
        overwatered=overwatered,
    )


def analyze(fields: list[Field], raw_by_id: dict[str, dict]) -> list[FieldAnalysis]:
    return [analyze_field(f, raw_by_id[f.id]) for f in fields]


def _fmt(value: Any, unit: str = "", digits: int = 1) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.{digits}f}{unit}"
    return f"{value}{unit}"


class ReportGenerator:
    def __init__(self, layout: FarmLayout, results: list[FieldAnalysis]):
        self.layout = layout
        self.results = results

    def build(self) -> str:
        today = date.today().isoformat()
        lines: list[str] = []
        lines.append(f"# Daily Operations Report — {self.layout.farm_name}")
        lines.append(f"_Generated {today} from live Open-Meteo grid data_\n")

        alerts = self._alerts_section()
        if alerts:
            lines.append("## ⚠️ Alerts")
            lines.extend(alerts)
            lines.append("")

        lines.append("## Field Summary\n")
        lines.append("| Field | Crop | Soil Moisture | Soil Temp | Air Temp | Frost Risk (tmrw) | Action |")
        lines.append("|---|---|---|---|---|---|---|")
        for r in self.results:
            action = self._action_label(r)
            lines.append(
                f"| {r.field.name} | {r.field.crop} "
                f"| {_fmt(r.soil_moisture, digits=3)} "
                f"| {_fmt(r.soil_temp_c, ' °C')} "
                f"| {_fmt(r.air_temp_c, ' °C')} "
                f"| {'YES (' + _fmt(r.temp_min_tomorrow_c, ' °C') + ')' if r.frost_risk else 'no'} "
                f"| {action} |"
            )

        lines.append("\n## Field Detail\n")
        for r in self.results:
            f = r.field
            lines.append(f"### {f.name} (`{f.id}`)")
            lines.append(f"- Crop: **{f.crop}**  |  Area: {f.area_hectares} ha  |  Zone: {f.irrigation_zone}")
            lines.append(f"- Coordinates: {f.latitude}, {f.longitude}")
            lines.append(
                f"- Current conditions: {_fmt(r.air_temp_c, ' °C')} air, "
                f"{_fmt(r.humidity_pct, '%')} humidity, "
                f"{_fmt(r.wind_kmh, ' km/h')} wind, "
                f"{_fmt(r.precip_now_mm, ' mm')} precip now"
            )
            lines.append(
                f"- Soil: {_fmt(r.soil_moisture, digits=3)} moisture "
                f"(m³/m³, 0–1cm), {_fmt(r.soil_temp_c, ' °C')} surface temp"
            )
            lines.append(
                f"- Today's forecast: high {_fmt(r.temp_max_today_c, ' °C')}, "
                f"precip sum {_fmt(r.precip_sum_today_mm, ' mm')}, "
                f"ET0 {_fmt(r.et0_today_mm, ' mm')}"
            )
            lines.append(f"- Tomorrow's min temp: {_fmt(r.temp_min_tomorrow_c, ' °C')}")
            lines.append(f"- **Recommended action:** {self._action_label(r)}")
            lines.append("")

        return "\n".join(lines)

    def _action_label(self, r: FieldAnalysis) -> str:
        actions = []
        if r.irrigation_needed:
            actions.append(f"Irrigate (deficit {_fmt(r.irrigation_deficit, digits=3)})")
        if r.overwatered:
            actions.append("Hold irrigation / check drainage")
        if r.frost_risk:
            actions.append("Activate frost protection")
        return "; ".join(actions) if actions else "No action needed"

    def _alerts_section(self) -> list[str]:
        lines = []
        for r in self.results:
            if r.frost_risk:
                lines.append(
                    f"- **Frost risk** on `{r.field.name}` — forecast low of "
                    f"{_fmt(r.temp_min_tomorrow_c, ' °C')} tomorrow "
                    f"(threshold {r.field.thresholds['frost_temp_c']} °C)."
                )
            if r.irrigation_needed:
                lines.append(
                    f"- **Low soil moisture** on `{r.field.name}` — "
                    f"{_fmt(r.soil_moisture, digits=3)} vs. minimum "
                    f"{r.field.thresholds['soil_moisture_min']}."
                )
            if r.overwatered:
                lines.append(
                    f"- **Soil moisture above max** on `{r.field.name}` — "
                    f"{_fmt(r.soil_moisture, digits=3)} vs. maximum "
                    f"{r.field.thresholds['soil_moisture_max']}."
                )
        return lines


class CommandGenerator:
    """
    Translates field analysis into a queue of mock hardware commands.
    This is a simulation layer — swap `build()`'s output for real device
    payloads (MQTT topics, Modbus registers, etc.) when wiring to actual
    controllers.
    """

    MIN_IRRIGATION_MINUTES = 10
    MAX_IRRIGATION_MINUTES = 120
    DEFICIT_TO_MINUTES = 1000

    def __init__(self, layout: FarmLayout, results: list[FieldAnalysis]):
        self.layout = layout
        self.results = results

    def build(self) -> dict[str, Any]:
        commands = []
        for r in self.results:
            commands.extend(self._commands_for_field(r))

        return {
            "farm_name": self.layout.farm_name,
            "generated_at": date.today().isoformat(),
            "command_count": len(commands),
            "commands": commands,
        }

    def _commands_for_field(self, r: FieldAnalysis) -> list[dict[str, Any]]:
        f = r.field
        cmds: list[dict[str, Any]] = []

        if r.irrigation_needed and r.irrigation_deficit is not None:
            minutes = int(
                min(
                    self.MAX_IRRIGATION_MINUTES,
                    max(
                        self.MIN_IRRIGATION_MINUTES,
                        r.irrigation_deficit * self.DEFICIT_TO_MINUTES,
                    ),
                )
            )
            cmds.append(
                {
                    "device": f"valve_{f.irrigation_zone}",
                    "field_id": f.id,
                    "action": "OPEN",
                    "duration_minutes": minutes,
                    "reason": (
                        f"soil_moisture {r.soil_moisture:.3f} below minimum "
                        f"{f.thresholds['soil_moisture_min']}"
                    ),
                }
            )
        elif r.overwatered:
            cmds.append(
                {
                    "device": f"valve_{f.irrigation_zone}",
                    "field_id": f.id,
                    "action": "HOLD",
                    "duration_minutes": 0,
                    "reason": (
                        f"soil_moisture {r.soil_moisture:.3f} above maximum "
                        f"{f.thresholds['soil_moisture_max']}"
                    ),
                }
            )

        if r.frost_risk:
            cmds.append(
                {
                    "device": f"frost_protection_{f.id}",
                    "field_id": f.id,
                    "action": "ON",
                    "duration_minutes": None,
                    "reason": (
                        f"forecast low {r.temp_min_tomorrow_c:.1f}°C tomorrow "
                        f"<= threshold {f.thresholds['frost_temp_c']}°C"
                    ),
                }
            )

        return cmds


# ==========================================================================
# Subcommand: variety  (farm_layout.json + genomic panel -> variety report)
# ==========================================================================


def _resolve_variety_engine(requested: str) -> str:
    """
    'auto' tries the real pybrops engine first (Linux/macOS only) and falls
    back to the from-scratch pymoo engine if pybrops/cyvcf2 aren't
    installed, so the same command works out of the box on either platform.
    """
    if requested != "auto":
        return requested
    try:
        import pybrops  # noqa: F401
        import cyvcf2  # noqa: F401
    except ImportError:
        return "pymoo"
    return "pybrops"


def cmd_variety(args: argparse.Namespace) -> int:
    engine = _resolve_variety_engine(args.engine)

    try:
        layout = FarmLayout(args.config)
    except FileNotFoundError:
        log.error("Config file not found: %s", args.config)
        return 1
    except (ValueError, json.JSONDecodeError) as e:
        log.error("Invalid farm_layout.json: %s", e)
        return 1

    corn_fields = [f for f in layout.fields if f.crop.strip().lower() == "corn"]
    other_crops = sorted({f.crop for f in layout.fields if f.crop.strip().lower() != "corn"})
    if other_crops:
        log.info(
            "Skipping %d field(s) with crop(s) %s - only 'corn' is supported "
            "(the only genotype panel currently vendored).",
            len(layout.fields) - len(corn_fields), other_crops,
        )
    if not corn_fields:
        log.warning("No corn fields found in %s - nothing to recommend.", args.config)
        return 1

    if not args.vcf.exists():
        log.error("Genotype panel not found: %s", args.vcf)
        return 1

    cache_path = args.cache
    if cache_path is None:
        cache_path = args.vcf.parent / (
            "ocs_frontier_cache_pybrops.npz" if engine == "pybrops" else "ocs_frontier_cache.npz"
        )

    try:
        if engine == "pybrops":
            log.info("Using engine: pybrops (real cyvcf2-backed OptimalContributionSubsetSelection)")
            import breeding_pybrops

            pgmat = breeding_pybrops.load_pgmat(args.vcf)
            gpmod = breeding_pybrops.build_genomic_model(pgmat.nvrnt)
            frontier = breeding_pybrops.compute_or_load_ocs_frontier(
                pgmat, gpmod,
                vcf_path=args.vcf,
                cache_path=cache_path,
                k=args.k,
                pop_size=args.pop_size,
                n_gen=args.n_gen,
                seed=args.seed,
                force_regen=args.regen_frontier,
            )
            gebv = gpmod.gebv(pgmat).mat
            taxa = list(pgmat.taxa)
        else:
            log.info("Using engine: pymoo (from-scratch reimplementation, no cyvcf2 required)")
            import breeding

            panel = breeding.load_vcf_dosage(args.vcf)
            marker_effects = breeding.build_synthetic_marker_effects(panel.dosage.shape[1])
            gebv = breeding.compute_gebv(panel, marker_effects)
            grm = breeding.compute_grm(panel.dosage)

            frontier = breeding.compute_or_load_frontier(
                panel, gebv, grm,
                vcf_path=args.vcf,
                cache_path=cache_path,
                k=args.k,
                pop_size=args.pop_size,
                n_gen=args.n_gen,
                seed=args.seed,
                force_regen=args.regen_frontier,
            )
            taxa = panel.taxa
    except ImportError as e:
        log.error("%s", e)
        return 1

    import breeding  # report/recommendation code lives here regardless of engine

    recommendations = breeding.build_recommendations(corn_fields, frontier, gebv, taxa)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    today_str = date.today().isoformat()
    report_path = args.output_dir / f"variety_recommendations_{today_str}.md"
    report_path.write_text(
        breeding.VarietyReportGenerator(layout.farm_name, recommendations).build(),
        encoding="utf-8",
    )
    log.info("Variety report written to %s", report_path)
    return 0


def cmd_breeding_sim(args: argparse.Namespace) -> int:
    try:
        import breeding
        import breeding_pybrops
    except ImportError as e:
        log.error("%s", e)
        return 1

    try:
        layout = FarmLayout(args.config)
    except FileNotFoundError:
        log.error("Config file not found: %s", args.config)
        return 1
    except (ValueError, json.JSONDecodeError) as e:
        log.error("Invalid farm_layout.json: %s", e)
        return 1

    corn_fields = [f for f in layout.fields if f.crop.strip().lower() == "corn"]
    if not corn_fields:
        log.warning("No corn fields found in %s - nothing to simulate.", args.config)
        return 1

    if not args.vcf.exists():
        log.error("Genotype panel not found: %s", args.vcf)
        return 1

    pgmat = breeding_pybrops.load_pgmat(args.vcf)
    breeding_pybrops.assign_approximate_genetic_map(pgmat)
    gpmod = breeding_pybrops.build_genomic_model(pgmat.nvrnt)

    results = []
    for f in corn_fields:
        log.info("Simulating %d generations for field %s...", args.generations, f.id)
        results.append(
            breeding_pybrops.simulate_breeding_generations(
                pgmat, gpmod, f,
                n_generations=args.generations,
                n_crosses=args.n_crosses,
                nmating=args.nmating,
                nprogeny=args.nprogeny,
                pop_size=args.pop_size,
                n_gen_ga=args.n_gen,
                seed=args.seed,
            )
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    today_str = date.today().isoformat()
    report_path = args.output_dir / f"breeding_simulation_{today_str}.md"
    report_path.write_text(
        breeding_pybrops.BreedingSimReportGenerator(layout.farm_name, results).build(),
        encoding="utf-8",
    )
    log.info("Breeding simulation report written to %s", report_path)
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    try:
        layout = FarmLayout(args.config)
    except FileNotFoundError:
        log.error("Config file not found: %s", args.config)
        return 1
    except (ValueError, json.JSONDecodeError) as e:
        log.error("Invalid farm_layout.json: %s", e)
        return 1

    try:
        raw_by_id = GridDataClient().fetch(layout.fields)
    except requests.RequestException as e:
        log.error("Failed to fetch grid data from Open-Meteo: %s", e)
        return 1
    except RuntimeError as e:
        log.error("%s", e)
        return 1

    results = analyze(layout.fields, raw_by_id)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    today_str = date.today().isoformat()

    if args.mode in ("report", "both"):
        report_path = args.output_dir / f"farm_report_{today_str}.md"
        report_path.write_text(ReportGenerator(layout, results).build(), encoding="utf-8")
        log.info("Report written to %s", report_path)

    if args.mode in ("commands", "both"):
        cmds_path = args.output_dir / f"hardware_commands_{today_str}.json"
        cmds = CommandGenerator(layout, results).build()
        cmds_path.write_text(json.dumps(cmds, indent=2), encoding="utf-8")
        log.info("Commands written to %s", cmds_path)

    return 0


# ==========================================================================
# CLI entry point
# ==========================================================================

def build_argparser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--verbose", action="store_true", help="Enable debug logging")

    parser = argparse.ArgumentParser(
        prog="farm_report.py",
        description="Generate a farm_layout.json from real OSM data, and/or pull live "
        "Open-Meteo grid data against it to produce a report and mock hardware commands.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    layout_p = sub.add_parser(
        "layout",
        parents=[common],
        help="Generate farm_layout.json from real OpenStreetMap farmland data",
    )
    loc = layout_p.add_mutually_exclusive_group(required=True)
    loc.add_argument("--place", type=str, help='Place name to geocode, e.g. "Story County, Iowa"')
    loc.add_argument("--bbox", type=str, help="Explicit bounding box: south,west,north,east")
    layout_p.add_argument("--max-fields", type=int, default=15, help="Cap on number of fields to output")
    layout_p.add_argument("--min-area-ha", type=float, default=0.5, help="Skip polygons smaller than this (hectares)")
    layout_p.add_argument("--farm-name", type=str, default=None, help="Name to use for 'farm_name' in the output")
    layout_p.add_argument("--output", type=Path, default=Path("farm_layout.json"), help="Output file path")
    layout_p.set_defaults(func=cmd_layout)

    run_p = sub.add_parser(
        "run",
        parents=[common],
        help="Fetch live grid data and generate a report and/or mock hardware commands",
    )
    run_p.add_argument("--config", type=Path, default=Path("farm_layout.json"), help="Path to farm_layout.json")
    run_p.add_argument("--mode", choices=["report", "commands", "both"], default="report")
    run_p.add_argument("--output-dir", type=Path, default=Path("."), help="Directory to write output file(s) into")
    run_p.set_defaults(func=cmd_run)

    variety_p = sub.add_parser(
        "variety",
        parents=[common],
        help="Recommend maize lines per field via climate-informed multi-objective "
        "genomic selection (phase 1 demo - see breeding.py)",
    )
    variety_p.add_argument("--config", type=Path, default=Path("farm_layout.json"), help="Path to farm_layout.json")
    variety_p.add_argument("--output-dir", type=Path, default=Path("."), help="Directory to write the report into")
    variety_p.add_argument(
        "--vcf", type=Path, default=Path(__file__).parent / "data" / "widiv_2000SNPs.vcf.gz",
        help="Genotype panel (VCF, optionally gzipped)",
    )
    variety_p.add_argument(
        "--cache", type=Path, default=None,
        help="Pareto frontier cache file (recomputed automatically if params change). "
        "Defaults to data/ocs_frontier_cache.npz (pymoo engine) or "
        "data/ocs_frontier_cache_pybrops.npz (pybrops engine) if omitted.",
    )
    variety_p.add_argument("--k", type=int, default=20, help="Number of parent lines per candidate cross set")
    variety_p.add_argument("--pop-size", type=int, default=80, help="NSGA-II population size")
    variety_p.add_argument("--n-gen", type=int, default=120, help="NSGA-II generations")
    variety_p.add_argument("--seed", type=int, default=1, help="NSGA-II random seed")
    variety_p.add_argument("--regen-frontier", action="store_true", help="Force recomputation even if a valid cache exists")
    variety_p.add_argument(
        "--engine", choices=["auto", "pymoo", "pybrops"], default="auto",
        help="'pybrops' uses the real pybrops+cyvcf2 OptimalContributionSubsetSelection "
        "(Linux/macOS only, see breeding_pybrops.py); 'pymoo' uses the from-scratch "
        "Windows-compatible reimplementation (see breeding.py); 'auto' (default) picks "
        "pybrops if installed, else falls back to pymoo",
    )
    variety_p.set_defaults(func=cmd_variety)

    breeding_sim_p = sub.add_parser(
        "breeding-sim",
        parents=[common],
        help="Simulate several generations of OCS-driven selection and mating per corn "
        "field using the real pybrops engine (phase 2 - see breeding_pybrops.py)",
    )
    breeding_sim_p.add_argument("--config", type=Path, default=Path("farm_layout.json"), help="Path to farm_layout.json")
    breeding_sim_p.add_argument("--output-dir", type=Path, default=Path("."), help="Directory to write the report into")
    breeding_sim_p.add_argument(
        "--vcf", type=Path, default=Path(__file__).parent / "data" / "widiv_2000SNPs.vcf.gz",
        help="Genotype panel (VCF, optionally gzipped)",
    )
    breeding_sim_p.add_argument("--generations", type=int, default=5, help="Number of generations to simulate")
    breeding_sim_p.add_argument("--n-crosses", type=int, default=10, help="Number of biparental crosses selected per generation")
    breeding_sim_p.add_argument("--nmating", type=int, default=1, help="Matings per selected cross")
    breeding_sim_p.add_argument("--nprogeny", type=int, default=10, help="Progeny produced per mating")
    breeding_sim_p.add_argument("--pop-size", type=int, default=60, help="NSGA-II population size (per generation)")
    breeding_sim_p.add_argument("--n-gen", type=int, default=60, help="NSGA-II generations (per breeding generation)")
    breeding_sim_p.add_argument("--seed", type=int, default=1, help="Random seed")
    breeding_sim_p.set_defaults(func=cmd_breeding_sim)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_argparser().parse_args(argv)
    if args.verbose:
        log.setLevel(logging.DEBUG)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
