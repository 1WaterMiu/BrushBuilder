"""PubChem identity lookup, independent of surface construction and force fields."""
import json
import time
from threading import Lock
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen, Request

_cache = {}
_lock = Lock()

def lookup(smiles):
    if not isinstance(smiles, str) or not 1 <= len(smiles) <= 2000:
        raise ValueError('Invalid monomer SMILES')
    # Serialize requests and cache by structure, not density or polymer length.
    with _lock:
        cached = _cache.get(smiles)
        if cached and cached[0] > time.monotonic():
            return cached[1]
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastidentity/smiles/cids/JSON?' + urlencode(
            dict(smiles=smiles, identity_type='same_stereo_isotope', MaxRecords=10))
        result = dict(status='unavailable', cids=[])
        try:
            req = Request(url, headers={'Accept': 'application/json', 'User-Agent': 'BrushBuilder/1.0'})
            with urlopen(req, timeout=20) as response:
                payload = json.load(response)
            cids = payload['IdentifierList']['CID']
            if not isinstance(cids, list) or not all(type(cid) is int and cid > 0 for cid in cids):
                raise ValueError('Invalid PubChem response')
            result = dict(status='found' if cids else 'not_found', cids=cids)
        except HTTPError as error:
            if error.code == 404:
                try:
                    if json.load(error).get('Fault', {}).get('Code') == 'PUGREST.NotFound':
                        result = dict(status='not_found', cids=[])
                except (ValueError, OSError):
                    pass
        except (URLError, TimeoutError, OSError, ValueError, KeyError, TypeError):
            pass
        _cache[smiles] = (time.monotonic() + (30 if result['status'] == 'unavailable' else 86400), result)
        return result
