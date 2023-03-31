import unittest

import ctxcore


class VersionTestCase(unittest.TestCase):
    """Version tests"""

    def test_version(self):
        """check ctxcore exposes a version attribute"""
        self.assertTrue(hasattr(ctxcore, "__version__"))
        self.assertIsInstance(ctxcore.__version__, str)


if __name__ == "__main__":
    unittest.main()
