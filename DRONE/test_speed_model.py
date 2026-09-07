#!/usr/bin/env python3
import math
import unittest

from speed_model import (
    A_SOUND,
    QuadConfig,
    analyze,
    catalogue,
    pack_voltage,
    pitch_speed_ms,
    recommend,
    rpm_no_load,
    solve_top_speed_ms,
    tip_speed_ms,
)


class Units(unittest.TestCase):
    def test_rpm_no_load(self):
        self.assertEqual(rpm_no_load(2550, 25.2), 64260)

    def test_pitch_speed(self):
        # 60k RPM, 5.1" pitch -> 60e3/60 * 0.12954 = 129.54 m/s
        v = pitch_speed_ms(60000, 5.1)
        self.assertAlmostEqual(v, 129.54, places=1)

    def test_tip_mach_formula(self):
        v = tip_speed_ms(37000, 5.1)
        self.assertAlmostEqual(v / A_SOUND, math.pi * 0.12954 * 37000 / 60 / A_SOUND, places=5)

    def test_6s_voltage(self):
        self.assertAlmostEqual(pack_voltage(6, False, 1.4), 6 * 4.2 - 1.4, places=6)
        self.assertAlmostEqual(pack_voltage(6, True, 1.8), 6 * 4.35 - 1.8, places=6)


class Ranking(unittest.TestCase):
    def setUp(self):
        self.results = [analyze(c) for c in catalogue()]
        self.by_name = {r.name: r for r in self.results}

    def test_reference_near_published_speed(self):
        ref = self.by_name["ref_mach_r5_class"]
        self.assertGreater(ref.estimated_top_kmh, 210)
        self.assertLess(ref.estimated_top_kmh, 280)

    def test_five_inch_beats_three_inch(self):
        self.assertGreater(
            self.by_name["speed5_6s_hdzero"].estimated_top_kmh,
            self.by_name["speed3_6s"].estimated_top_kmh + 25,
        )

    def test_six_s_beats_four_s_on_current(self):
        six = self.by_name["speed5_6s_hdzero"]
        four = self.by_name["speed5_4s_ultrahigh_kv"]
        self.assertLess(six.pack_current_a, four.pack_current_a)
        self.assertGreater(six.estimated_top_kmh, four.estimated_top_kmh - 15)

    def test_top_speed_below_pitch_speed(self):
        for r in self.results:
            self.assertLessEqual(r.estimated_top_kmh, r.pitch_speed_kmh + 0.1)

    def test_high_kv_6s_flags_tip_or_current(self):
        rec = self.by_name["speed5_6s_hdzero"]
        self.assertTrue(rec.tip_mach_warning or rec.tip_mach >= 0.85)

    def test_recommend_picks_a_five_inch_6s(self):
        rec = recommend(self.results)
        self.assertEqual(rec.cells, 6)
        self.assertGreaterEqual(rec.size_in, 5.0)

    def test_solver_zero_when_pitch_zero(self):
        cfg = QuadConfig(
            name="x",
            size_in=5.0,
            pitch_in=0.0,
            blades=2,
            cells=6,
            kv=2000,
            motor_stator="2207",
            auw_kg=0.5,
            cda_m2=0.006,
            pack_mah=1000,
            pack_c=100,
        )
        self.assertEqual(solve_top_speed_ms(cfg, 10000), 0.0)


if __name__ == "__main__":
    unittest.main()
