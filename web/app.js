const $=s=>document.querySelector(s);
let catalog,model,key='CBAA',density='high',table='repeatAtoms',view='side',angle=0,zoom=1,request=0,currentJob=null,poll=null;
function settings(){return {template:key,density,repeats:Number($('#repeats').value),methyl:$('#methyl').checked};}
function status(message,kind=''){const el=$('#status');el.textContent=message;el.className=kind;}
async function json(url,options){const r=await fetch(url,options);if(!r.ok)throw Error((await r.text()).replace(/<[^>]*>/g,' ').slice(-1800));return r.json();}
async function load(){
 const id=++request;currentJob=null;clearTimeout(poll);$('#download').disabled=true;$('#optimize').disabled=true;
 if(!$('#repeats').checkValidity()||!Number.isInteger(Number($('#repeats').value))){status('链长请输入 1–60 的整数','error');return;}
 status('正在独立生成坐标并重建拓扑…');
 try{const data=await json('/api/model?'+new URLSearchParams(settings()));if(id!==request)return;model=data;render();$('#optimize').disabled=false;status('初始坐标已生成。所有链共用独立几何模板；点击能量最小化后检查参数与 silica 位移。');}
 catch(e){if(id===request)status('生成失败：'+e.message,'error');}
}
function render(){
 const t=catalog.templates.find(t=>t.id===key);$('#linker').textContent=t.linker;
 $('#density').textContent=model.graftDensity.toFixed(4);$('#chains').textContent=model.chainCount;$('#area').textContent=model.area.toFixed(2);$('#charge').textContent=(Math.abs(model.charge)<1e-8?0:model.charge).toFixed(5);
 const done=model.status==='completed',top=model.topology;
 const fields={'几何模板':model.geometryVersion,'原子 / 固定':model.atoms.length.toLocaleString()+' / '+model.frozen.length.toLocaleString(),'链长 / 甲基':model.repeats+' / '+(model.methyl?'有':'无'),'最小化':done?'已收敛':'尚未完成','Silica 位移':done?'0 nm · 已核验':'待核验'};
 $('#source').innerHTML=Object.entries(fields).map(([k,v])=>'<dt>'+k+'</dt><dd>'+v+'</dd>').join('');
 const terms={Bond:top.bonds.length,Angle:top.angles.length,Dihedral:top.dihedrals.filter(r=>r[4]==='9').length,Improper:top.dihedrals.filter(r=>r[4]==='2').length,'1–4 pairs':top.pairs.length};
 $('#counts').innerHTML=Object.entries(terms).map(([k,v])=>'<dt>'+k+'</dt><dd>'+v.toLocaleString()+'</dd>').join('');
 $('#viewName').textContent=(view==='top'?'俯视':'侧视图')+' · '+(done?'优化后坐标':'初始坐标');renderTable();draw();
}
function renderTable(){
 if(!model)return;const top=model.topology;let rows,headers;
 if(table==='repeatAtoms'||table==='anchorAtoms'){
  headers=['Index','Type','Residue','Atom','Charge / e','Mass'];rows=top.atoms.filter(a=>table==='repeatAtoms'?a[2]==='3':Number(a[2])<=2).map(a=>[a[0],a[1],a[3],a[4],a[6],a[7]]);
 }else{rows=table==='proper'||table==='improper'?top.dihedrals.filter(r=>r[4]===(table==='proper'?'9':'2')):top[table];headers=table==='bonds'?['Atom i','Atom j','Function']:table==='angles'?['Atom i','Atom j','Atom k','Function']:['Atom i','Atom j','Atom k','Atom l','Function'];}
 $('#thead').innerHTML='<tr>'+headers.map(x=>'<th>'+x+'</th>').join('')+'</tr>';
 $('#tbody').innerHTML=rows.map(a=>'<tr>'+a.map(x=>'<td>'+x+'</td>').join('')+'</tr>').join('');
}
$('#optimize').onclick=async()=>{
 const version=request;$('#optimize').disabled=true;$('#download').disabled=true;status('正在提交 GROMACS 优化任务…');
 try{const result=await json('/api/jobs',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(settings())});if(version!==request)return;currentJob=result.id;watch(result.id,version);}
 catch(e){status(e.message,'error');$('#optimize').disabled=false;}
};
async function watch(id,version){
 try{const result=await json('/api/job?id='+id);if(version!==request)return;status(result.message,result.status==='failed'?'error':result.status==='completed'?'success':'running');
  if(result.status==='completed'){model=await json('/api/result?id='+id);if(version!==request)return;render();$('#download').disabled=false;$('#optimize').disabled=false;}
  else if(result.status==='failed'){$('#optimize').disabled=false;}
  else poll=setTimeout(()=>watch(id,version),1800);
 }catch(e){if(version===request){status(e.message,'error');$('#optimize').disabled=false;}}
}
$('#download').onclick=()=>{if(currentJob&&model.status==='completed')location.href='/api/download?id='+currentJob;};
const canvas=$('#canvas'),ctx=canvas.getContext('2d');
function draw(){
 if(!model)return;const dpr=devicePixelRatio||1,w=canvas.clientWidth,h=canvas.clientHeight;canvas.width=w*dpr;canvas.height=h*dpr;ctx.scale(dpr,dpr);ctx.clearRect(0,0,w,h);
 const box=model.box;let minz=Infinity,maxz=-Infinity;for(const a of model.atoms){minz=Math.min(minz,a.xyz[2]);maxz=Math.max(maxz,a.xyz[2]);}
 const scale=Math.min(w/(box[0]*1.35),h/((view==='top'?box[1]:maxz-minz)*1.25))*zoom;
 function project(a){let x=a[0]-box[0]/2,y=a[1]-box[1]/2,z=a[2]-(maxz+minz)/2;return [w/2+(x*Math.cos(angle)-y*Math.sin(angle))*scale,h/2-(view==='top'?(x*Math.sin(angle)+y*Math.cos(angle)):z)*scale,y];}
 const projected=model.atoms.map(a=>project(a.xyz));ctx.lineWidth=.65;
 for(const [i,j]of model.bonds){const a=model.atoms[i],b=model.atoms[j];if(a.name.startsWith('H')||b.name.startsWith('H'))continue;if(a.xyz.some((v,k)=>Math.abs(v-b.xyz[k])>box[k]/2))continue;ctx.strokeStyle=a.frozen?'#abb5a466':'#4b806766';ctx.beginPath();ctx.moveTo(...projected[i].slice(0,2));ctx.lineTo(...projected[j].slice(0,2));ctx.stroke();}
 const order=model.atoms.map((a,i)=>i).sort((i,j)=>projected[i][2]-projected[j][2]);
 for(const i of order){const a=model.atoms[i];if(a.name.startsWith('H'))continue;ctx.fillStyle=a.frozen?'#9ca893':a.name.startsWith('C')?'#477d64':'#cc996b';ctx.beginPath();ctx.arc(projected[i][0],projected[i][1],(a.frozen?1.3:1.7)*Math.sqrt(zoom),0,Math.PI*2);ctx.fill();}
}
let drag=null;canvas.onpointerdown=e=>{drag=e.clientX;canvas.setPointerCapture(e.pointerId);};canvas.onpointermove=e=>{if(drag!==null){angle+=(e.clientX-drag)*.01;drag=e.clientX;draw();}};canvas.onpointerup=canvas.onpointercancel=()=>drag=null;
canvas.addEventListener('wheel',e=>{e.preventDefault();zoom=Math.min(3,Math.max(.5,zoom*Math.exp(-e.deltaY*.001)));draw();},{passive:false});new ResizeObserver(draw).observe(canvas);
for(const v of ['side','top'])$('#'+v).onclick=()=>{view=v;$('#side').classList.toggle('active',v==='side');$('#top').classList.toggle('active',v==='top');render();};
$('#reset').onclick=()=>{angle=0;zoom=1;draw();};$('#repeats').onchange=load;$('#methyl').onchange=load;
document.querySelectorAll('[data-density]').forEach(b=>b.onclick=()=>{density=b.dataset.density;document.querySelectorAll('[data-density]').forEach(x=>x.classList.toggle('active',x===b));load();});
document.querySelectorAll('[data-table]').forEach(b=>b.onclick=()=>{table=b.dataset.table;document.querySelectorAll('[data-table]').forEach(x=>x.classList.toggle('active',x===b));renderTable();});
json('/api/catalog').then(data=>{catalog=data;$('#templates').innerHTML=data.templates.map(t=>'<button class="template '+(t.id===key?'active':'')+'" data-template="'+t.id+'"><strong>'+(t.id==='CBAA'?'CBAA3':'SBMA')+'</strong><small>'+t.linker.toUpperCase()+' / 独立几何模板</small></button>').join('');document.querySelectorAll('[data-template]').forEach(b=>b.onclick=()=>{key=b.dataset.template;document.querySelectorAll('[data-template]').forEach(x=>x.classList.toggle('active',x===b));load();});load();}).catch(e=>status(e.message,'error'));
