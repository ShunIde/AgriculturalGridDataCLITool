"""
Offline smoke test: monkeypatches GridDataClient.fetch so the full
config -> analysis -> report/commands pipeline can be verified without
a live network call to Open-Meteo.
"""
from pathlib import Path
from unittest.mock import patch

import farm_report as fr


def fake_raw_for(field: fr.Field) -> dict:
    # Synthesize plausible-looking data per field so we can exercise
    # every branch: normal, low-moisture, over-watered, frost risk.
    if field.id == "field-01":
        # low soil moisture -> irrigation needed
        return {
            "current": {
                "temperature_2m": 18.4, "relative_humidity_2m": 55,
                "precipitation": 0.0, "wind_speed_10m": 12.1,
                "soil_temperature_0cm": 16.2, "soil_moisture_0_to_1cm": 0.15,
            },
            "daily": {
                "temperature_2m_max": [22.0, 21.0],
                "temperature_2m_min": [9.0, 8.5],
                "precipitation_sum": [0.0, 1.2],
                "et0_fao_evapotranspiration": [3.1, 2.8],
            },
        }
    if field.id == "field-02":
        # over-watered
        return {
            "current": {
                "temperature_2m": 19.0, "relative_humidity_2m": 70,
                "precipitation": 2.0, "wind_speed_10m": 8.0,
                "soil_temperature_0cm": 17.0, "soil_moisture_0_to_1cm": 0.50,
            },
            "daily": {
                "temperature_2m_max": [20.0, 19.5],
                "temperature_2m_min": [10.0, 9.0],
                "precipitation_sum": [5.0, 0.5],
                "et0_fao_evapotranspiration": [2.0, 2.2],
            },
        }
    # field-03: frost risk tomorrow, moisture fine
    return {
        "current": {
            "temperature_2m": 5.0, "relative_humidity_2m": 60,
            "precipitation": 0.0, "wind_speed_10m": 3.0,
            "soil_temperature_0cm": 4.0, "soil_moisture_0_to_1cm": 0.28,
        },
        "daily": {
            "temperature_2m_max": [10.0, 9.0],
            "temperature_2m_min": [3.0, 0.5],
            "precipitation_sum": [0.0, 0.0],
            "et0_fao_evapotranspiration": [1.0, 0.9],
        },
    }


def fake_fetch(self, fields):
    return {f.id: fake_raw_for(f) for f in fields}


def run():
    layout = fr.FarmLayout(Path("farm_layout.example.json"))

    with patch.object(fr.GridDataClient, "fetch", fake_fetch):
        raw_by_id = fr.GridDataClient().fetch(layout.fields)

    results = fr.analyze(layout.fields, raw_by_id)

    report_md = fr.ReportGenerator(layout, results).build()
    commands = fr.CommandGenerator(layout, results).build()

    Path("test_out_report.md").write_text(report_md, encoding="utf-8")
    Path("test_out_commands.json").write_text(
        __import__("json").dumps(commands, indent=2), encoding="utf-8"
    )

    # Basic assertions to sanity-check the logic branches
    r_by_id = {r.field.id: r for r in results}
    assert r_by_id["field-01"].irrigation_needed is True, "field-01 should need irrigation"
    assert r_by_id["field-02"].overwatered is True, "field-02 should be flagged overwatered"
    assert r_by_id["field-03"].frost_risk is True, "field-03 should show frost risk"

    cmd_devices = [c["device"] for c in commands["commands"]]
    assert any("valve_zone-A" in d for d in cmd_devices), "expected irrigation command for zone-A"
    assert any("valve_zone-B" in d for d in cmd_devices), "expected HOLD command for zone-B"
    assert any("frost_protection_field-03" in d for d in cmd_devices), "expected frost command"

    print("All offline assertions passed.")
    print("\n--- report preview ---\n")
    print(report_md[:800])
    print("\n--- commands preview ---\n")
    print(__import__("json").dumps(commands, indent=2)[:800])


if __name__ == "__main__":
    run()
