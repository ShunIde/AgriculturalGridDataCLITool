"""
Offline smoke test for osm_to_farm_layout.py — verifies polygon
centroid/area math and Overpass-element-to-field conversion using a
synthetic Overpass response, without hitting the live API.
"""
import json
from pathlib import Path

import farm_report as ofl


def make_square_way(way_id, center_lat, center_lon, half_side_m, tags):
    """Builds a roughly square way of known size for area-math verification."""
    import math
    dlat = half_side_m / 110_540.0
    dlon = half_side_m / (111_320.0 * math.cos(math.radians(center_lat)))
    corners = [
        (center_lat - dlat, center_lon - dlon),
        (center_lat - dlat, center_lon + dlon),
        (center_lat + dlat, center_lon + dlon),
        (center_lat + dlat, center_lon - dlon),
        (center_lat - dlat, center_lon - dlon),  # close ring
    ]
    return {
        "type": "way",
        "id": way_id,
        "geometry": [{"lat": lat, "lon": lon} for lat, lon in corners],
        "tags": tags,
    }


def run():
    # A 200m x 200m square == 4 hectares, centered at a known point.
    known_center = (41.90, -93.60)
    square = make_square_way(
        1001, known_center[0], known_center[1], half_side_m=100,
        tags={"landuse": "farmland", "crop": "soybean", "name": "Test Square Field"},
    )

    # A small sliver that should get filtered out by min_area_ha
    tiny = make_square_way(
        1002, 41.91, -93.61, half_side_m=5,  # ~0.01 ha
        tags={"landuse": "orchard"},
    )

    # A non-way element that should be ignored
    node = {"type": "node", "id": 2001, "lat": 41.92, "lon": -93.62, "tags": {}}

    elements = [square, tiny, node]

    fields = ofl.elements_to_fields(elements, max_fields=15, min_area_ha=0.5)

    assert len(fields) == 1, f"expected 1 field after filtering, got {len(fields)}"
    f = fields[0]

    # Area should be very close to 4 ha (200m x 200m)
    assert abs(f["area_hectares"] - 4.0) < 0.05, f"area off: {f['area_hectares']}"

    # Centroid should be very close to the known center
    assert abs(f["latitude"] - known_center[0]) < 0.0005, f"lat off: {f['latitude']}"
    assert abs(f["longitude"] - known_center[1]) < 0.0005, f"lon off: {f['longitude']}"

    assert f["crop"] == "soybean"
    assert f["name"] == "Test Square Field"
    assert f["id"] == "field-01"
    assert f["thresholds"]["frost_temp_c"] == ofl.DEFAULT_THRESHOLDS["farmland"]["frost_temp_c"]

    output = {"farm_name": "Test Region", "fields": fields}
    Path("test_out_farm_layout.json").write_text(json.dumps(output, indent=2), encoding="utf-8")

    print("All OSM conversion assertions passed.")
    print(json.dumps(fields, indent=2))


if __name__ == "__main__":
    run()
