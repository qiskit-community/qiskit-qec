import unittest
from plerco.pec_python.qec_utils.indexer import Indexer


class TestIndexer(unittest.TestCase):
    def test_add_and_lookup(self):
        a = Indexer()
        lst = [4, 8, 2, 5, 1]
        idx = [0, 1, 2, 3, 4]
        # Add elements from l to the indexer
        # Confirm each assigned index
        for i in range(5):
            self.assertEqual(a[lst[i]], idx[i])
        # Lookup each element by its index
        # Confirm the returned element of l
        for i in range(5):
            self.assertEqual(a.rlookup(idx[i]), lst[i])

    def test_add_twice(self):
        a = Indexer()
        lst = [4, 8, 8]  # contains a duplicate
        idx = [0, 1, 1]
        for i in range(3):
            self.assertEqual(a[lst[i]], idx[i])
        for i in range(3):
            self.assertEqual(a.rlookup(idx[i]), lst[i])


if __name__ == "__main__":
    unittest.main()
