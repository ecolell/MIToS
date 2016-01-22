import unittest
from Bio import AlignIO
from mitos.MSA import coverage, columngappercentage, gapstrip


class TestMSA(unittest.TestCase):

    def setUp(self):
        self.msa = AlignIO.read('tests/samples/protein.sto', 'stockholm')

    def test_coverage(self):
        covers_expected = [1.0, 0.718, 0.674, 0.722, 0.727, 0.907, 0.907, 0.947,
                           0.952, 0.762 ]
        covers_obtain = map(lambda s: coverage(self.msa[0], s), self.msa)
        map(lambda c: self.assertAlmostEqual(*c, places=3),
            zip(covers_obtain, covers_expected))
        cover = coverage(self.msa[0], self.msa[1])
        self.assertAlmostEquals(cover, 0.718, 3)
        cover = coverage(self.msa[0], self.msa[2])
        self.assertAlmostEquals(cover, 0.674, 3)

    def test_columngappercentage(self):
        gapp_e = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.8, 0.6, 0.6, 0.6, 0.5, 0.6, 0.6]
        gapp_o = map(lambda i: columngappercentage(self.msa[:, i]),
                     range(self.msa.get_alignment_length()))
        map(lambda c: self.assertAlmostEqual(*c, places=3),
            zip(gapp_o, gapp_e))

    def test_gapstrip(self):
        result = gapstrip(self.msa, 'Whale',
                          coverage_limit=0.8,
                          gap_limit=0.000000001)
        self.assertEquals(len(result), 5)
        self.assertEquals(len(result[0]), 227)



if __name__ == '__main__':
    unittest.main()
