import unittest
from cistargetx.common.lrucache import LRUCache, cache

class LRUCacheTestCase(unittest.TestCase):
    def setUp(self):
        self.cache = LRUCache(3)

    def test_set_item(self):
        self.cache['1'] = 1
        self.assertTrue('1' in self.cache)

    def test_len(self):
        self.cache['1'] = 1
        self.assertEqual(len(self.cache), 1)
        self.cache['2'] = 2
        self.assertEqual(len(self.cache), 2)

    def test_get_item(self):
        self.cache['1'] = 1
        self.cache['2'] = 2
        self.assertEqual(self.cache['1'], 1)
        self.assertEqual(self.cache['2'], 2)

    def test_replace(self):
        self.cache['1'] = 1
        self.cache['2'] = 2
        self.assertEqual(self.cache['2'], 2)
        self.assertEqual(len(self.cache), 2)
        self.cache['2'] = 3
        self.assertEqual(self.cache['2'], 3)
        self.assertEqual(len(self.cache), 2)

    def test_lru(self):
        self.cache['1'] = 1
        self.cache['2'] = 2
        self.cache['3'] = 3
        self.assertEqual(len(self.cache), 3)
        self.cache['4'] = 4
        self.assertEqual(len(self.cache), 3)
        self.assertTrue('2' in self.cache)
        self.assertTrue('3' in self.cache)
        self.assertTrue('4' in self.cache)
        self.assertFalse('1' in self.cache)

    def test_function_cache(self):
        @cache(2)
        def function(x): return x * 2
        function(1)
        function(2)
        function(3)
        self.assertEqual(function(1), 2)
        self.assertEqual(function(2), 4)
        self.assertEqual(function(3), 6)


if __name__ == '__main__':
    unittest.main()
