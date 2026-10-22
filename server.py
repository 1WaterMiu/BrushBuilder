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
    from chemistry import normalize
    spec=normalize(result['key'],q.get('headGroup'),q.get('linker'),q.get('spacer'))
    result.update(head_group=spec[0],linker=spec[1],spacer=spec[2])
    return result

class Handler(SimpleHTTPRequestHandler):
    def __init__(self,*args,**kwargs): super().__init__(*args,directory=str(ROOT/'web'),**kwargs)
    def do_GET(self):
        url=urlparse(self.path)
        if not url.path.startswith('/api/'):
            return super().do_GET()
        try:
            if url.path=='/api/catalog': data=json.dumps(catalog()).encode(); mime='application/json'
            elif url.path in ('/api/model','/api/structure'):
                q=parse_qs(url.query)
                values={k:v[0] for k,v in q.items()}
                if values.get('methyl','true') not in ('true','false'):raise ValueError('Methyl must be true or false')
                values['methyl']=values.get('methyl','true')=='true'
                model=build(**parameters(values))
                if url.path=='/api/structure':
                    from structure_export import archive_structure
                    data=archive_structure(model);mime='application/zip'
                else:data=json.dumps(model).encode(); mime='application/json'
            elif url.path in ('/api/job','/api/result','/api/download'):
                q=parse_qs(url.query); record=job(q.get('id',[''])[0])
                if url.path=='/api/job': data=json.dumps({k:v for k,v in record.items() if k!='directory'}).encode(); mime='application/json'
                else:
                    if record['status']!='completed': raise ValueError('Minimization not completed; MD-ready download is blocked')
                    from pathlib import Path
                    previous=json.loads((Path(record['directory'])/'model.json').read_text())
                    if 'parameterization' not in previous.get('chemistry',{}):
                        raise ValueError('此旧任务使用模板电荷；当前 CGenFF-only 流程禁止作为新结果下载')
                    data=(Path(record['directory'])/('surface-md.zip' if url.path.endswith('download') else 'model.json')).read_bytes()
                    mime='application/zip' if url.path.endswith('download') else 'application/json'
            else: self.send_error(404); return
            self.send_response(200); self.send_header('Content-Type',mime)
            if mime=='application/zip': self.send_header('Content-Disposition','attachment; filename="'+('brush-structure.zip' if url.path=='/api/structure' else 'brush-surface.zip')+'"')
            self.send_header('Content-Length',str(len(data))); self.end_headers(); self.wfile.write(data)
        except (ValueError,TypeError,AssertionError,RuntimeError) as e:
            self.send_error(400,'Invalid request',str(e))

    def do_POST(self):
        if self.path!='/api/jobs': self.send_error(404); return
        try:
            # Local same-origin job creation only; avoid driving compute from foreign sites.
            if self.headers.get('Origin') not in (None,'http://127.0.0.1:8765','http://localhost:8765'):
                raise ValueError('Invalid origin')
            length=int(self.headers.get('Content-Length','0'))
            if not 0<length<4096: raise ValueError('Invalid request size')
            values=parameters(json.loads(self.rfile.read(length)))
            from cgenff_library import status
            state=status((values['head_group'],values['linker'],values['spacer']),values['methyl'],values['repeats'])
            if not state['mdReady']:
                raise ValueError(state['reason'])
            data=json.dumps(dict(id=submit(values))).encode()
            self.send_response(202); self.send_header('Content-Type','application/json'); self.send_header('Content-Length',str(len(data))); self.end_headers(); self.wfile.write(data)
        except (ValueError,TypeError) as e: self.send_error(400,'Invalid request',str(e))

if __name__=='__main__':
    print('BrushBuilder: http://127.0.0.1:8765',flush=True)
    ThreadingHTTPServer(('127.0.0.1',8765),Handler).serve_forever()
