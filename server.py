import json
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from urllib.parse import urlparse, parse_qs
from builder import ROOT, catalog, build, archive

class Handler(SimpleHTTPRequestHandler):
    def __init__(self,*args,**kwargs): super().__init__(*args,directory=str(ROOT/'web'),**kwargs)
    def do_GET(self):
        url=urlparse(self.path)
        if not url.path.startswith('/api/'):
            return super().do_GET()
        try:
            if url.path=='/api/catalog': data=json.dumps(catalog()).encode(); mime='application/json'
            elif url.path in ('/api/model','/api/download'):
                q=parse_qs(url.query)
                model=build(q.get('template',['CBAA'])[0],q.get('density',['high'])[0],int(q.get('seed',['0'])[0]),int(q.get('repeats',['20'])[0]))
                data=archive(model) if url.path.endswith('download') else json.dumps(model).encode()
                mime='application/zip' if url.path.endswith('download') else 'application/json'
            else: self.send_error(404); return
            self.send_response(200); self.send_header('Content-Type',mime)
            if mime=='application/zip': self.send_header('Content-Disposition','attachment; filename="brush-surface.zip"')
            self.send_header('Content-Length',str(len(data))); self.end_headers(); self.wfile.write(data)
        except (ValueError,AssertionError) as e:
            self.send_error(400,str(e))

if __name__=='__main__':
    print('BrushBuilder: http://127.0.0.1:8765',flush=True)
    ThreadingHTTPServer(('127.0.0.1',8765),Handler).serve_forever()
