import io
import unittest
from unittest.mock import patch
from urllib.error import HTTPError, URLError
import pubchem_lookup as pubchem

class PubChemTests(unittest.TestCase):
    def setUp(self):
        pubchem._cache.clear()

    def test_match_and_cache(self):
        with patch.object(pubchem, 'urlopen', return_value=io.BytesIO(b'{"IdentifierList":{"CID":[3034140]}}')) as request:
            self.assertEqual(pubchem.lookup('C=C')['cids'], [3034140])
            self.assertEqual(pubchem.lookup('C=C')['status'], 'found')
            self.assertEqual(request.call_count, 1)
            self.assertIn('same_stereo_isotope', request.call_args.args[0].full_url)

    def test_not_found_is_distinct_from_network_failure(self):
        error = HTTPError('url', 404, 'Not found', {}, io.BytesIO(b'{"Fault":{"Code":"PUGREST.NotFound"}}'))
        with patch.object(pubchem, 'urlopen', side_effect=error):
            self.assertEqual(pubchem.lookup('CC')['status'], 'not_found')
        with patch.object(pubchem, 'urlopen', side_effect=URLError('offline')):
            self.assertEqual(pubchem.lookup('CCC')['status'], 'unavailable')

if __name__ == '__main__':
    unittest.main()
