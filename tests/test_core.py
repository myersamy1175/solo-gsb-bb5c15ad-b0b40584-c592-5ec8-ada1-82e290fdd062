"""核心算法单元测试：python -m pytest（或 python -m unittest）。"""

import unittest

from app import seq_utils
from app.alignment import sw_align
from app.screening import screen


class TestSeqUtils(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(seq_utils.reverse_complement("ATGC"), "GCAT")
        self.assertEqual(seq_utils.reverse_complement("GCGATCCAGTCGGC"), "GCCGACTGGATCGC")

    def test_reverse_complement_involution(self):
        for s in ("ATGC", "GCGATCC", "NNN".replace("N", "A")):
            self.assertEqual(seq_utils.reverse_complement(seq_utils.reverse_complement(s)), s)

    def test_validate(self):
        seq_utils.validate("ACGTacgt\n")  # 合法
        with self.assertRaises(ValueError):
            seq_utils.validate("ACGTN")

    def test_duplex_tm_gc_rich_longer_higher(self):
        tm_at = seq_utils.duplex_tm("ATATATAT", 0.1, 250)
        tm_gc = seq_utils.duplex_tm("GCGCGCGC", 0.1, 250)
        self.assertGreater(tm_gc, tm_at)

    def test_hairpin_loop_penalty(self):
        tm_open = seq_utils.duplex_tm("GCGCGCGC", 0.1, 250)
        tm_loop = seq_utils.duplex_tm("GCGCGCGC", 0.1, 250, loop_nt=3)
        self.assertLess(tm_loop, tm_open)


class TestAlignment(unittest.TestCase):
    def test_perfect_3prime_dimer(self):
        # F 3' 尾 S, R 3' 尾 rc(S) => 经典 3'-3' 二聚体，14bp 完美
        stem = "GCGATCCAGTCGGC"
        f = "ATGCCTAGGCTAACG" + stem
        r = "TTAGCGTACGGTACG" + seq_utils.reverse_complement(stem)
        hits = sw_align(f, seq_utils.reverse_complement(r))
        self.assertTrue(hits)
        top = hits[0]
        self.assertEqual(top["mismatches"], 0)
        self.assertEqual(top["gaps"], 0)
        self.assertGreaterEqual(top["matched"], 14)
        # A 侧命中 3' 末端区间 16-29
        self.assertEqual((top["start1"] + 1, top["end1"] + 1), (16, 29))

    def test_no_alignment_for_unrelated_sequences(self):
        # rc 空间「相同碱基」才互补；两条无同源关系的短序列不应产生达阈命中
        hits = sw_align("ACGTCGTA", seq_utils.reverse_complement("CATGCATG"))
        self.assertEqual(hits, [])


class TestScreening(unittest.TestCase):
    def setUp(self):
        stem = "GCGATCCAGTCGGC"
        rc_stem = seq_utils.reverse_complement(stem)
        hstem = "GCGCGGCCGATGCT"
        rc_hstem = seq_utils.reverse_complement(hstem)
        self.primers = [
            {"name": "F1", "sequence": "ATGCCTAGGCTAACG" + stem},
            {"name": "R1", "sequence": "TTAGCGTACGGTACG" + rc_stem},
            {"name": "HP1", "sequence": hstem + "TTT" + rc_hstem + "GGTT"},
            {"name": "C1", "sequence": "GGTATACCGGATCGATTCG"},
        ]
        self.params = {
            "temperature": 50,
            "salt_mm": 100,
            "primer_nm": 250,
            "min_loop": 3,
        }

    def test_detects_3prime_dimer_and_rejects_pair(self):
        res = screen(self.primers, self.params)
        cross = [
            d for d in res["dimers"]
            if d["primer_a"] == "F1" and d["primer_b"] == "R1" and not d["self_dimer"]
        ]
        self.assertTrue(cross)
        best = cross[0]
        self.assertGreaterEqual(best["longest_match"], 14)
        self.assertTrue(best["anchored_3prime"])
        self.assertIn(best["risk"], ("high", "severe"))
        self.assertIn(["F1", "R1"], res["summary"]["rejected_pairs"])

    def test_detects_hairpin_with_localizable_arms(self):
        res = screen(self.primers, self.params)
        hp = [h for h in res["hairpins"] if h["primer"] == "HP1"]
        self.assertTrue(hp)
        h = hp[0]
        (a, b) = h["arms"][0], h["arms"][1]
        # 两臂区间 1-14 与 18-31，环 3 nt
        self.assertEqual((a["start"], a["end"]), (1, 14))
        self.assertEqual((b["start"], b["end"]), (18, 31))
        self.assertEqual(h["loop_nt"], 3)
        self.assertGreaterEqual(h["longest_match"], 14)
        # 3' 末端额外 4nt，发卡不锚定 3'
        self.assertFalse(h["anchored_3prime"])
        self.assertIn("HP1", res["summary"]["rejected_primers"])

    def test_higher_temperature_fewer_high_risk(self):
        cold = screen(self.primers, {**self.params, "temperature": 50})
        hot = screen(self.primers, {**self.params, "temperature": 75})
        self.assertGreaterEqual(
            len([h for h in cold["hairpins"] if h["risk"] in ("high", "severe")]),
            len([h for h in hot["hairpins"] if h["risk"] in ("high", "severe")]),
        )

    def test_arm_sequences_are_reverse_complementary(self):
        res = screen(self.primers, self.params)
        h = next(x for x in res["hairpins"] if x["primer"] == "HP1")
        seq = next(p["sequence"] for p in self.primers if p["name"] == "HP1")
        a = seq[h["arms"][0]["start"] - 1:h["arms"][0]["end"]]
        b = seq[h["arms"][1]["start"] - 1:h["arms"][1]["end"]]
        self.assertEqual(seq_utils.reverse_complement(a), b)


if __name__ == "__main__":
    unittest.main()
