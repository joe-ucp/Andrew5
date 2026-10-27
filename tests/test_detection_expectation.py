"""Enumeration tests use labelled molecule subsets, independent of the PMF."""
from itertools import combinations
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from revision_detection_expectation import detection_probability, expected_components


class DetectionExpectationTests(unittest.TestCase):
    def test_all_tiny_scalar_libraries_by_enumeration(self):
        for total in range(9):
            for count in range(total + 1):
                for depth in range(total + 1):
                    subsets = list(combinations(range(total), depth))
                    truth = sum(any(i < count for i in subset) for subset in subsets) / len(subsets)
                    actual = detection_probability(total, count, depth)
                    self.assertAlmostEqual(float(actual), truth, places=13,
                                           msg=f'N={total}, c={count}, d={depth}')

    def test_occupancy_with_unequal_blocks_and_outside_molecules(self):
        # Gene categories: two up genes (2,1 copies), one down gene (2),
        # and two outside-block molecules. Each molecule is distinguishable.
        labels = [0, 0, 1, 2, 2, 3, 3]
        for depth in range(len(labels) + 1):
            empirical = []
            for subset in combinations(range(len(labels)), depth):
                detected = {labels[i] for i in subset}
                kup = sum(g in detected for g in [0, 1])
                kdown = int(2 in detected)
                empirical.append([.5 * (kup / 2 - kdown), .5 * kup / 2,
                                  -.5 * kdown, kup, kdown])
            truth = np.mean(empirical, axis=0)
            actual = expected_components(np.array([[2], [1]]), np.array([[2]]), np.array([7]), depth)
            np.testing.assert_allclose([actual[k][0] for k in ['occupancy_component', 'occupancy_up',
                                                              'occupancy_down', 'k_up', 'k_down']],
                                       truth, atol=2e-14, rtol=0)

    def test_vector_boundary_support_and_rare_singleton(self):
        np.testing.assert_array_equal(detection_probability([0, 5, 5, 5], [0, 0, 1, 5], [0, 5, 5, 0]),
                                      [0., 0., 1., 0.])
        self.assertEqual(float(detection_probability(10**9, 1, 1)), 1e-9)
        self.assertEqual(float(detection_probability(11, 5, 1)), 5 / 11)
        self.assertEqual(float(detection_probability(11, 1, 5)), 5 / 11)
        empty = expected_components(np.empty((2, 0)), np.empty((3, 0)), np.array([], dtype=int), 10)
        self.assertEqual(empty['occupancy_component'].shape, (0,))

    def test_invalid_probabilities_and_block_shapes(self):
        for triple in [(2, 3, 1), (2, 1, 3), (-1, 0, 0), (2, -1, 1), (2, 1, -1),
                       (2.5, 1, 1), (2, .5, 1), (2, 1, .5), (np.nan, 0, 0),
                       (np.inf, 0, 0), (True, 1, 1), (2, True, 1), (2, 1, True)]:
            with self.subTest(triple=triple), self.assertRaises(ValueError):
                detection_probability(*triple)
        for up, down, total in [(np.array([1]), np.array([[1]]), np.array([2])),
                                (np.empty((0, 1)), np.array([[1]]), np.array([2])),
                                (np.array([[1, 0]]), np.array([[1]]), np.array([2])),
                                (np.array([[2]]), np.array([[2]]), np.array([3]))]:
            with self.assertRaises(ValueError):
                expected_components(up, down, total, 1)


if __name__ == '__main__':
    unittest.main()
