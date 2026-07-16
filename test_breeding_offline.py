"""
Offline smoke test for breeding.py's field-weighting and recommendation
logic. Builds a tiny synthetic Pareto frontier by hand (rather than running
the real NSGA-II search) so this runs instantly and does not require pymoo
to be installed just to check the recommendation math.
"""
from dataclasses import dataclass

import numpy as np

import breeding as br


@dataclass
class FakeField:
    id: str
    name: str
    thresholds: dict


def run():
    # weight derivation: lower thresholds should push weight toward that trait
    cold_w, drought_w = br.field_trait_weights(frost_temp_c=-8.0, soil_moisture_min=0.45)
    assert cold_w > drought_w, "very low frost threshold should weight cold_tolerance higher"

    cold_w2, drought_w2 = br.field_trait_weights(frost_temp_c=6.0, soil_moisture_min=0.02)
    assert drought_w2 > cold_w2, "very low soil-moisture threshold should weight drought_tolerance higher"

    # a tiny hand-built frontier: 4 taxa, 3 candidate 2-of-4 subsets spanning
    # the trade-off space (cold-favoring, drought-favoring, balanced)
    taxa = ["A", "B", "C", "D"]
    masks = np.array(
        [
            [True, True, False, False],   # A+B: high cold, low drought
            [False, False, True, True],   # C+D: low cold, high drought
            [True, False, False, True],   # A+D: balanced
        ]
    )
    objectives = np.array(
        [
            [-10.0, -1.0, 0.01],  # -meanBV_cold, -meanBV_drought, kinship
            [-1.0, -10.0, 0.01],
            [-5.0, -5.0, 0.005],
        ]
    )
    frontier = br.ParetoFrontier(taxa=taxa, masks=masks, objectives=objectives, k=2)

    cold_field = FakeField("f-cold", "Cold Field", {"frost_temp_c": -8.0, "soil_moisture_min": 0.45})
    drought_field = FakeField("f-drought", "Drought Field", {"frost_temp_c": 6.0, "soil_moisture_min": 0.02})

    cold_rec = br.recommend_for_field(cold_field, frontier)
    drought_rec = br.recommend_for_field(drought_field, frontier)

    assert set(cold_rec.selected_taxa) == {"A", "B"}, f"expected A+B for cold-leaning field, got {cold_rec.selected_taxa}"
    assert set(drought_rec.selected_taxa) == {"C", "D"}, f"expected C+D for drought-leaning field, got {drought_rec.selected_taxa}"

    # build_recommendations should rank each field's own top_lines by its own weights
    gebv = np.array(
        [
            [8.0, 2.0],   # A
            [12.0, 0.0],  # B
            [1.0, 9.0],   # C
            [-1.0, 11.0], # D
        ]
    )
    recs = br.build_recommendations([cold_field, drought_field], frontier, gebv, taxa)
    cold_top = dict((t, (c, d)) for t, c, d in recs[0].top_lines)
    assert list(cold_top.keys())[0] == "B", f"cold field's top line should be B (highest cold BV), got {recs[0].top_lines}"

    print("All breeding offline assertions passed.")
    print("cold weights:", (cold_w, drought_w), "-> selected", cold_rec.selected_taxa)
    print("drought weights:", (cold_w2, drought_w2), "-> selected", drought_rec.selected_taxa)


if __name__ == "__main__":
    run()
