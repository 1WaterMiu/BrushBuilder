import json
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from urllib.parse import urlparse, parse_qs
from builder import ROOT, catalog, build, archive
from simulation import submit,job

def parameters(q):
    raw=q.get('repeats',20)
    if isinstance(raw,bool) or str(raw).strip()!=str(int(raw)): raise ValueError('Repeat units must be an integer')
    result=dict(key=q.get('template','CBAA'),density=q.get('density','high'),repeats=int(raw),methyl=q.get('methyl',True))
    if result['key'] not in ('CBAA','SBMA') or result['density'] not in ('high','low') or not 1<=result['repeats']<=60 or not isinstance(result['methyl'],bool):
        raise ValueError('Invalid construction parameters')
    return result

class Handler(SimpleHTTPRequestHandler):
    def __init__(self,*args,**kwargs): super().__init__(*args,directory=str(ROOT/'web'),**kwargs)
    def do_GET(self):
        url=urlparse(self.path)
        if not url.path.startswith('/api/'):
            return super().do_GET()
        try:
            if url.path=='/api/catalog': data=json.dumps(catalog()).encode(); mime='application/json'
            elif url.path=='/api/model':
                q=parse_qs(url.query)
                values={k:v[0] for k,v in q.items()}; values['methyl']=values.get('methyl','true')=='true'
                model=build(**parameters(values))
                data=json.dumps(model).encode(); mime='application/json'
            elif url.path in ('/api/job','/api/result','/api/download'):
                q=parse_qs(url.query); record=job(q.get('id',[''])[0])
                if url.path=='/api/job': data=json.dumps({k:v for k,v in record.items() if k!='directory'}).encode(); mime='application/json'
                else:
                    if record['status']!='completed': raise ValueError('Minimization not completed; MD-ready download is blocked')
                    from pathlib import Path
                    data=(Path(record['directory'])/('surface-md.zip' if url.path.endswith('download') else 'model.json')).read_bytes()
                    mime='application/zip' if url.path.endswith('download') else 'application/json'
            else: self.send_error(404); return
            self.send_response(200); self.send_header('Content-Type',mime)
            if mime=='application/zip': self.send_header('Content-Disposition','attachment; filename="brush-surface.zip"')
            self.send_header('Content-Length',str(len(data))); self.end_headers(); self.wfile.write(data)
        except (ValueError,AssertionError,RuntimeError) as e:
            self.send_error(400,str(e))

    def do_POST(self):
        if self.path!='/api/jobs': self.send_error(404); return
        try:
            # Local same-origin job creation only; avoid driving compute from foreign sites.
            if self.headers.get('Origin') not in (None,'http://127.0.0.1:8765','http://localhost:8765'):
                raise ValueError('Invalid origin')
            length=int(self.headers.get('Content-Length','0'))
            if not 0<length<4096: raise ValueError('Invalid request size')
            values=parameters(json.loads(self.rfile.read(length)))
            data=json.dumps(dict(id=submit(values))).encode()
            self.send_response(202); self.send_header('Content-Type','application/json'); self.send_header('Content-Length',str(len(data))); self.end_headers(); self.wfile.write(data)
        except (ValueError,TypeError) as e: self.send_error(400,str(e))

if __name__=='__main__':
    print('BrushBuilder: http://127.0.0.1:8765',flush=True)
    ThreadingHTTPServer(('127.0.0.1',8765),Handler).serve_forever()
