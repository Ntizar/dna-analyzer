/* DNA Analyzer — analisis de ADN crudo 100% en el navegador (privacidad total).
   Ascendencia continental (GrafAnc/GRAF-pop, Jin et al. HGG Advances 2025),
   rasgos, y sexo. Ningun dato sale de la maquina del usuario. */
'use strict';

/* ============================ CONSTANTES ============================ */
const MAGIC = 'DNAANC';
const SCORE_IDX1 = [9,9,9,9,3,3,6,18,18,24,16,17,23,20,20];
const SCORE_IDX2 = [11,12,10,13,5,4,8,20,19,25,15,14,22,21,7];
const SCORE_NAMES = ['EA1','EA2','EA3','EA4','AF1','AF2','AF3','EU1','EU2','EU3','SA1','SA2','IC1','IC2','IC3'];
const MIN_SNPS = 100;
const AFR_POS = [1.08, 1.10, 0.00];
const BASE = {'A':0,'C':1,'G':2,'T':3};
const COMP = {0:3, 3:0, 1:2, 2:1};
const BASE_CH = ['A','C','G','T'];

/* ============================ HELPERS VECTORES ============================ */
function subV(a,b){ return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]; }
function addV(a,b){ return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]; }
function rotX(p,t){ t=t*Math.PI/180; const y=p[1],z=p[2]; return [p[0], y*Math.cos(t)-z*Math.sin(t), y*Math.sin(t)+z*Math.cos(t)]; }
function rotY(p,t){ t=t*Math.PI/180; const x=p[0],z=p[2]; return [x*Math.cos(t)+z*Math.sin(t), p[1], -x*Math.sin(t)+z*Math.cos(t)]; }
function rotZ(p,t){ t=t*Math.PI/180; const x=p[0],y=p[1]; return [x*Math.cos(t)-y*Math.sin(t), x*Math.sin(t)+y*Math.cos(t), p[2]]; }
function setPoint(dist){ return [dist[2], dist[0], dist[1]]; } // dist=(e,f,a) -> (x=a,y=e,z=f)

function transformAll(eD,fD,aD,sD){
  let eP=setPoint(eD), fP=setPoint(fD), aP=setPoint(aD), sP=setPoint(sD);
  const dx=fP[0], dy=fP[1], dz=fP[2];
  const off=[dx,dy,dz];
  eP=subV(eP,off); fP=subV(fP,off); aP=subV(aP,off); sP=subV(sP,off);
  const t1=Math.atan2(aP[1],aP[0])*-180/Math.PI;
  eP=rotZ(eP,t1); fP=rotZ(fP,t1); aP=rotZ(aP,t1); sP=rotZ(sP,t1);
  const t2=Math.atan2(aP[2],aP[0])*180/Math.PI;
  eP=rotY(eP,t2); fP=rotY(fP,t2); aP=rotY(aP,t2); sP=rotY(sP,t2);
  const t3=Math.atan2(eP[1],eP[2])*180/Math.PI; const t3b=t3-90;
  eP=rotX(eP,t3b); fP=rotX(fP,t3b); aP=rotX(aP,t3b); sP=rotX(sP,t3b);
  eP=addV(eP,off); fP=addV(fP,off); aP=addV(aP,off); sP=addV(sP,off);
  const d2=[AFR_POS[0]-fP[0], AFR_POS[1]-fP[1], AFR_POS[2]-fP[2]];
  eP=addV(eP,d2); fP=addV(fP,d2); aP=addV(aP,d2); sP=addV(sP,d2);
  return {eP,fP,aP,sP};
}
function bary(eP,fP,aP,sP){
  const x1=eP[0],y1=eP[1],x2=fP[0],y2=fP[1],x3=aP[0],y3=aP[1],xp=sP[0],yp=sP[1];
  const det=(y2-y3)*(x1-x3)+(x3-x2)*(y1-y3);
  const eWt=((y2-y3)*(xp-x3)+(x3-x2)*(yp-y3))/det;
  const fWt=((y3-y1)*(xp-x3)+(x1-x3)*(yp-y3))/det;
  const aWt=1-eWt-fWt;
  return {eWt,fWt,aWt};
}

/* ============================ PANEL DE REFERENCIA ============================ */
let PANEL = null; // {n, chr, pos, rs, ref, alt, afs(Uint8Array n*26), rsToIdx, vtxExp, vtxExpGds, g0e,g0f,g0a}

async function loadPanel(){
  const resp = await fetch('data/panel.bin.gz');
  if(!resp.ok) throw new Error('No se pudo cargar el panel de referencia.');
  const stream = resp.body.pipeThrough(new DecompressionStream('gzip'));
  const buf = await new Response(stream).arrayBuffer();
  const dv = new DataView(buf);
  const magic = String.fromCharCode(dv.getUint8(0),dv.getUint8(1),dv.getUint8(2),dv.getUint8(3),dv.getUint8(4),dv.getUint8(5));
  if(magic!==MAGIC) throw new Error('Panel de referencia invalido.');
  const n = dv.getUint32(7, true);
  const chr = new Uint8Array(n), pos = new Uint32Array(n), rs = new Uint32Array(n);
  const ref = new Uint8Array(n), alt = new Uint8Array(n), afs = new Uint16Array(n*26);
  let o = 11;
  const rsToIdx = new Map();
  for(let i=0;i<n;i++){
    chr[i]=dv.getUint8(o); pos[i]=dv.getUint32(o+1,true); rs[i]=dv.getUint32(o+5,true);
    ref[i]=dv.getUint8(o+9); alt[i]=dv.getUint8(o+10);
    for(let j=0;j<26;j++) afs[i*26+j]=dv.getUint16(o+11+j*2,true);
    o += 63;
    rsToIdx.set(rs[i], i);
  }
  // precomputar vtxExpGenoDists y vtxExpGds
  const vtxExp = [];
  for(let v=0;v<3;v++){ const arr=[]; for(let r=0;r<3;r++) arr.push(new Float64Array(n)); vtxExp.push(arr); }
  const popExpPe=[0,0,0], popExpPf=[0,0,0], popExpPa=[0,0,0];
  for(let i=0;i<n;i++){
    const eur=afs[i*26]/65535, afr=afs[i*26+1]/65535, eas=afs[i*26+2]/65535;
    const pev=eur, pfv=afr, pav=eas, qev=1-pev, qfv=1-pfv, qav=1-pav;
    const aaPev=Math.log(pev)*2, bbPev=Math.log(qev)*2, abPev=Math.log(pev)+Math.log(qev)+Math.log(2);
    const aaPfv=Math.log(pfv)*2, bbPfv=Math.log(qfv)*2, abPfv=Math.log(pfv)+Math.log(qfv)+Math.log(2);
    const aaPav=Math.log(pav)*2, bbPav=Math.log(qav)*2, abPav=Math.log(pav)+Math.log(qav)+Math.log(2);
    for(let v=0;v<3;v++){
      const pv=(v===0?eur:v===1?afr:eas), qv=1-pv;
      const eGd=aaPev*pv*pv+bbPev*qv*qv+abPev*2*pv*qv;
      const fGd=aaPfv*pv*pv+bbPfv*qv*qv+abPfv*2*pv*qv;
      const aGd=aaPav*pv*pv+bbPav*qv*qv+abPav*2*pv*qv;
      vtxExp[v][0][i]=eGd; vtxExp[v][1][i]=fGd; vtxExp[v][2][i]=aGd;
      popExpPe[v]+=eGd; popExpPf[v]+=fGd; popExpPa[v]+=aGd;
    }
  }
  const vtxExpGds = [0,1,2].map(v=>[-popExpPe[v]/n, -popExpPf[v]/n, -popExpPa[v]/n]);
  const g0 = transformAll(vtxExpGds[0],vtxExpGds[1],vtxExpGds[2],vtxExpGds[0]);
  PANEL = {n, chr, pos, rs, ref, alt, afs, rsToIdx, vtxExp, vtxExpGds, g0e:g0.eP, g0f:g0.fP, g0a:g0.aP};
  return PANEL;
}

/* ============================ CODIFICACION GENOTIPO ============================ */
function baseCode(b){ return BASE[b]!==undefined ? BASE[b] : -1; }
// Devuelve el numero de alelos ALT (0,1,2) o -1 si no determinable
function countAlt(geno, refC, altC){
  let n=0;
  for(const c of geno){
    const b=baseCode(c);
    if(b===-1) return -1;
    if(b===altC) n++;
    else if(b===refC){ /* ref */ }
    else {
      const comp=COMP[b];
      if(comp===altC) n++;
      else if(comp===refC){ /* comp de ref */ }
      else return -1;
    }
  }
  return n;
}

/* ============================ PARSEO DE CSV ============================ */
// Detecta formato y devuelve {geno: Map(rs->[chr,pos,gt]), sexo, xHet}
function parseRaw(text){
  const lines = text.split(/\r?\n/);
  let header = null, colMap = null, start = 0;
  // busca cabecera de columnas (no comentarios #)
  for(let i=0;i<lines.length;i++){
    const l = lines[i].trim();
    if(!l || l.startsWith('#')) continue;
    const cols = l.split(/[\t,]/).map(x=>x.trim().toLowerCase());
    if(cols.some(c=>c==='rsid'||c==='rs_id'||c==='rs')){
      header = cols; start = i+1; break;
    }
  }
  // si no hay cabecera reconocible, intentar detectar por primeras lineas de datos
  const geno = new Map();
  let chrIdx=-1,posIdx=-1,gtIdx=-1,rsIdx=-1, a1Idx=-1,a2Idx=-1;
  if(header){
    rsIdx = header.findIndex(c=>c==='rsid'||c==='rs_id'||c==='rs');
    chrIdx = header.findIndex(c=>c==='chromosome'||c==='chrom'||c==='chr');
    posIdx = header.findIndex(c=>c==='position'||c==='pos');
    gtIdx = header.findIndex(c=>c==='result'||c==='genotype'||c==='gt');
    a1Idx = header.findIndex(c=>c==='allele1');
    a2Idx = header.findIndex(c=>c==='allele2');
  } else {
    // sin cabecera: asumir rsid,cromosoma,posicion,genotipo (MyHeritage)
    rsIdx=0; chrIdx=1; posIdx=2; gtIdx=3;
    start=0;
  }
  let xTotal=0, xHet=0, xChr = (chrIdx>=0 ? chrIdx : 1);
  for(let i=start;i<lines.length;i++){
    const l=lines[i]; if(!l || l.startsWith('#')) continue;
    const parts = l.split(/[\t,]/).map(x=>x.trim().replace(/"/g,''));
    if(parts.length<4) continue;
    const rs = parts[rsIdx];
    const gt = (gtIdx>=0 ? parts[gtIdx] : (parts[a1Idx]+parts[a2Idx])).toUpperCase();
    if(!/^[ACGT]{1,2}$/.test(gt)) continue;
    const chr = chrIdx>=0 ? parts[chrIdx] : '';
    const pos = posIdx>=0 ? parts[posIdx] : '';
    if((chr==='X'||chr==='x'||chr==='23'||chr==='XY')){ xTotal++; if(gt.length===2 && gt[0]!==gt[1]) xHet++; }
    geno.set(rs.replace(/^rs/,''), {chr, pos, gt});
  }
  const sex = xTotal>200 ? (xHet/xTotal < 0.1 ? 'M' : 'F') : '?';
  return {geno, sex, xHet, xTotal};
}

/* ============================ GRAFANC ============================ */
function computeAncestry(genoMap){
  const P = PANEL;
  const popPv=[0,0,0,0,0], refPopSnps=[0,0,0,0,0];
  const vtxPe=[0,0,0], vtxPf=[0,0,0], vtxPa=[0,0,0];
  const smpGd=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  const nomP1=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  const nomP2=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  const smpN=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  let numGeno=0, het=0;
  const af=(i,j)=>P.afs[i*26+j]/65535;
  for(const [rs, rec] of genoMap){
    const i = P.rsToIdx.get(Number(rs)); if(i===undefined) continue;
    const geno = countAlt(rec.gt, P.ref[i], P.alt[i]);
    if(geno===-1) continue;
    // 5 ref pops
    for(let pop=0;pop<5;pop++){
      const pv=af(i,pop); if(pv>0&&pv<1){ const qv=1-pv;
        if(geno===0) popPv[pop]+=Math.log(qv)*2;
        else if(geno===1) popPv[pop]+=Math.log(pv*qv*2);
        else if(geno===2) popPv[pop]+=Math.log(pv)*2;
        refPopSnps[pop]++;
      }
    }
    for(let v=0;v<3;v++){ vtxPe[v]+=P.vtxExp[v][0][i]; vtxPf[v]+=P.vtxExp[v][1][i]; vtxPa[v]+=P.vtxExp[v][2][i]; }
    for(let s=0;s<15;s++){
      const p1=af(i,SCORE_IDX1[s]), p2=af(i,SCORE_IDX2[s]);
      if(p1>0&&p1<1&&p2>0&&p2<1){
        const logp=Math.log(p2/p1), logq=Math.log((1-p2)/(1-p1));
        if(geno===2) smpGd[s]+=logp*2; else if(geno===1) smpGd[s]+=logp+logq; else if(geno===0) smpGd[s]+=logq*2;
        nomP1[s]+=(p1*logp+(1-p1)*logq)*2;
        nomP2[s]+=(p2*logp+(1-p2)*logq)*2;
        smpN[s]++;
      }
    }
    if(geno===1) het++;
    numGeno++;
  }
  if(numGeno<MIN_SNPS) return null;
  const popMean=[0,0,0,0,0];
  for(let p=0;p<5;p++) if(refPopSnps[p]>0) popMean[p]=-popPv[p]/refPopSnps[p];
  const subPopScores=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
  for(let s=0;s<15;s++){ if(smpN[s]===0) continue;
    const raw=smpGd[s]/smpN[s], n1=nomP1[s]/smpN[s], n2=nomP2[s]/smpN[s];
    subPopScores[s]=-1+2*(raw-n1)/(n2-n1);
  }
  const smpDist=[popMean[0],popMean[1],popMean[2]];
  const vtxExpDists=[0,1,2].map(v=>[-vtxPe[v]/numGeno, -vtxPf[v]/numGeno, -vtxPa[v]/numGeno]);
  const t = transformAll(vtxExpDists[0],vtxExpDists[1],vtxExpDists[2],smpDist);
  const b = bary(t.eP,t.fP,t.aP,t.sP);
  const gd1=b.eWt*P.g0e[0]+b.fWt*P.g0f[0]+b.aWt*P.g0a[0];
  const gd2=b.eWt*P.g0e[1]+b.fWt*P.g0f[1]+b.aWt*P.g0a[1];
  const gd3=t.sP[2];
  const ej=Math.max(b.eWt,0), fj=Math.max(b.fWt,0), aj=Math.max(b.aWt,0);
  const tot=ej+fj+aj;
  const ePct=ej*100/tot, fPct=fj*100/tot, aPct=aj*100/tot;
  const ancGroupId = assignGroup(aPct,fPct,gd1,gd2,gd3,subPopScores);
  return {numGeno, hetRate: het/numGeno, gd1, gd2, gd3, ePct, fPct, aPct,
          rawPe:b.eWt*100, rawPf:b.fWt*100, rawPa:b.aWt*100, ancGroupId, subPopScores};
}

function assignGroup(aPct,fPct,gd1,gd2,gd3,sp){
  const [ea1,ea2,ea4,af1,af2,af3,eu1,eu2,eu3,sa1,sa2,ic1,ic2,ic3]=
    [sp[0],sp[1],sp[3],sp[4],sp[5],sp[6],sp[7],sp[8],sp[9],sp[10],sp[11],sp[12],sp[13],sp[14]];
  let anc=0;
  if(aPct>50&&fPct>10&&ic1>0.4&&gd3>0.035) anc=700;
  else if(ic1>0.5&&fPct<15){ if(sa2<-0.1)anc=402; else if(sa2>0.7)anc=403; else if(sa1<-0.6)anc=404; else if(sa1>0.1)anc=405; else anc=401; }
  else if(aPct>15&&fPct>15) anc=800;
  else if(ic1>-0.3&&aPct>40){
    if(ea2>1.4)anc=501; else if(ea2>0.4)anc=502;
    else if(ea4>-0.7&&ea4<-0.2&&ea1<-0.3) anc=(ea4<-0.48?506:508);
    else if(ea2>-0.3&&ea1<-0.3)anc=503;
    else if(ea4>-0.2&&ea4>ea1+0.1)anc=504;
    else if(ea1<-0.6)anc=505; else if(ea1<-0.1)anc=507; else if(ea1<0.5)anc=509; else if(ea2>-1)anc=510; else anc=511;
  }
  else if(eu1<1.6&&ic2<-0.25){
    if(ic3<3.5-3.1*eu1){
      if(ic3>-0.3)anc=308; else if(eu2>0.3)anc=301;
      else if(eu3>-0.4-eu1){ if(eu3>eu1+1.1)anc=305; else if(eu3<1.3*eu1-1.3)anc=307; else anc=306; }
      else { if(eu3<eu1-1.0)anc=304; else if(eu2>-0.75&&eu1<-0.75)anc=302; else anc=303; }
    } else anc=(ic3>0.3?203:202);
  }
  else if(eu1>1.6&&ic2<0.4) anc=(ic3<-0.3?201:202);
  else if(eu1>2.2&&ic2<1.4) anc=106;
  else if(gd1>1.4758){ if(ic1<-0.2) anc=(gd3>0.05?603:602); else anc=308; }
  else {
    if(ic1<0.9*gd1-1.36) anc=601;
    else { if(af3<af1*4/3+0.4){ if(af3>0.5)anc=108; else { if(af1>0.6)anc=(af3>-0.1?104:105); else anc=103; } }
      else { if(af1>-0.8&&af2>-0.6&&af2<0.9*af1+0.8)anc=107; else { if(af3>0.5)anc=108; else if(af2<-0.1)anc=101; else anc=102; } } }
  }
  return anc;
}

/* ============================ RASGOS ============================ */
const TRAITS = [
  { rs:'4988235', name:'Lactosa (LCT/MCM6)', gene:'MCM6', effect:'T', desc:'Persistencia de lactasa',
    interp:g=>{const t=(g.match(/T/g)||[]).length; return t>=1?['Tolerante a la lactosa','Posee el alelo T de persistencia.']:['Probable intolerancia','No porta el alelo T. (Sitio multi-alelico, ver nota)'];}},
  { rs:'762551', name:'Cafeina (CYP1A2)', gene:'CYP1A2', effect:'A', desc:'Metabolismo de la cafeina',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return a===2?['Metabolizador rapido','AA']:a===1?['Metabolizador intermedio','AC']:['Metabolizador lento','CC'];}},
  { rs:'12913832', name:'Color de ojos (HERC2)', gene:'HERC2', effect:'A', desc:'Tendencia del color',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return a===2?['Tendencia a ojos claros','AA (azul/verde)']:a===1?['Color mixto','AG']:['Tendencia a ojos oscuros','GG'];}},
  { rs:'4680', name:'COMT (Val158Met)', gene:'COMT', effect:'A', desc:'Metabolismo de la dopamina',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return a===2?['Met/Met','Degradacion lenta de dopamina']:a===1?['Val/Met','Intermedio']:['Val/Val','Degradacion rapida'];}},
  { rs:'6265', name:'BDNF (Val66Met)', gene:'BDNF', effect:'T', desc:'Factor neurotrofico',
    interp:g=>{const t=(g.match(/T/g)||[]).length; return t===2?['Met/Met','']:t===1?['Val/Met','']:['Val/Val','Variante comun'];}},
  { rs:'9939609', name:'FTO (peso)', gene:'FTO', effect:'A', desc:'Tendencia al peso',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return a===2?['Dos alelos de riesgo','']:a===1?['Un alelo de riesgo','']:['Sin alelo de riesgo',''];}},
  { rs:'1229984', name:'Alcohol (ADH1B)', gene:'ADH1B', effect:'C', desc:'Metabolismo del alcohol',
    interp:g=>{const c=(g.match(/C/g)||[]).length; return c>=1?['Metabolismo alcohol rapido','Variante Arg48 (menor riesgo de alcoholismo)']:['Metabolismo alcohol tipico',''];}},
  { rs:'1801133', name:'MTHFR (folato)', gene:'MTHFR', effect:'T', desc:'Metabolismo del folato',
    interp:g=>{const t=(g.match(/T/g)||[]).length; return t===2?['TT (677TT)','Actividad reducida del folato']:t===1?['CT (677CT)','Intermedio']:['CC (677CC)','Normal'];}},
];

function computeTraits(genoMap){
  const out=[];
  for(const t of TRAITS){
    const rec=genoMap.get(t.rs);
    if(!rec || !/^[ACGT]{1,2}$/.test(rec.gt)) continue;
    const g=rec.gt.toUpperCase();
    const [label,detail]=t.interp(g);
    out.push({name:t.name, gene:t.gene, desc:t.desc, genotype:g, label, detail});
  }
  return out;
}

/* ============================ UI ============================ */
const $=s=>document.querySelector(s);
const state={};
async function handleFile(file){
  if(!file) return;
  setStatus('Descomprimiendo zip...', 5);
  let text=null;
  const lower=file.name.toLowerCase();
  if(lower.endsWith('.zip')){
    const zip=await JSZip.loadAsync(file);
    let csvName=null;
    zip.forEach((path,entry)=>{ if(!entry.dir && /\.(csv|txt)$/i.test(path) && !csvName) csvName=path; });
    if(!csvName) throw new Error('No se encontro un CSV dentro del zip.');
    text=await zip.file(csvName).async('string');
  } else if(/\.(csv|txt)$/i.test(lower)){
    text=await file.text();
  } else throw new Error('Formato no soportado. Sube un .zip o .csv.');
  setStatus('Procesando datos...', 30);
  const parsed=parseRaw(text);
  if(parsed.geno.size<1000) throw new Error('Se detectaron muy pocos SNP ('+parsed.geno.size+'). Revisa el fichero.');
  setStatus('Calculando ascendencia...', 55);
  await loadPanel();
  const anc=computeAncestry(parsed.geno);
  const traits=computeTraits(parsed.geno);
  setStatus('Renderizando...', 90);
  render(parsed, anc, traits);
  setStatus('Listo.', 100);
}
function setStatus(msg,pct){ $('#status').textContent=msg; if(pct!==undefined) $('#bar').style.width=pct+'%'; }
function render(parsed, anc, traits){
  $('#results').style.display='block';
  const sex=parsed.sex==='M'?'Hombre':parsed.sex==='F'?'Mujer':'No determinado';
  $('#sexo').textContent=sex;
  if(!anc){ $('#anc-box').innerHTML='<div class="warn">No hay suficientes SNP para calcular ascendencia ('+parsed.geno.size+' encontrados).</div>'; }
  else {
    $('#anc-box').innerHTML=
      '<div class="pct"><span style="color:#2f6fb2">Europea '+anc.ePct.toFixed(2)+'%</span>'+
      '<span style="color:#8a5a2b">Africana '+anc.fPct.toFixed(2)+'%</span>'+
      '<span style="color:#b23f3f">Este Asiatica '+anc.aPct.toFixed(2)+'%</span></div>'+
      '<div class="meta">GD1 '+anc.gd1.toFixed(4)+' · GD2 '+anc.gd2.toFixed(4)+' · GD3 '+anc.gd3.toFixed(4)+
      ' · SNP '+anc.numGeno+' · Heterocigosidad '+(anc.hetRate*100).toFixed(1)+'%</div>'+
      '<div class="group">Grupo subcontinental: <b>'+anc.ancGroupId+'</b> · '+groupName(anc.ancGroupId)+'</div>';
    drawEFAtriangle(anc);
  }
  const tEl=$('#traits'); tEl.innerHTML='';
  if(traits.length===0) tEl.innerHTML='<div class="warn">No se encontraron SNP de rasgos.</div>';
  traits.forEach(t=>{ const d=document.createElement('div'); d.className='card';
    d.innerHTML='<div class="t-name">'+t.name+' <span class="t-gene">'+t.gene+'</span></div>'+
      '<div class="t-geno">Genotipo: <b>'+t.genotype+'</b></div><div class="t-label">'+t.label+'</div>'+
      (t.detail?'<div class="t-detail">'+t.detail+'</div>':'');
    tEl.appendChild(d); });
  $('#resumen').textContent='SNP procesados: '+parsed.geno.size+' · Datos 100% en tu navegador, nada se guarda.';
}
function groupName(id){
  const m={300:'Europeo',303:'Europa Occidental',304:'Europa Meridional',302:'Norte de Europa',301:'Finlandia',305:'Nordeste de Europa',306:'Sudeste de Europa',307:'Balcanes',308:'Otra Europa',100:'Africano',101:'Nigeria',102:'Africa Occidental',103:'Africa Central',104:'Kenya',105:'Africa Meridional',106:'Africa Nororiental',107:'Afroamericano',108:'Otra Africa',200:'MEN',201:'Norte de Africa',202:'Oriente Medio 2',203:'Oriente Medio 1',400:'Sur de Asia',401:'India',402:'Gujarati',403:'Pakistan',404:'Sri Lanka',405:'Bangladesh',500:'Este de Asia',501:'Ryukyu',502:'Japon',503:'Corea',504:'Norte de Asia',505:'Norte China 1',506:'Norte China 2',507:'Sur China 1',508:'Sur China 2',509:'Sudeste Asiatico',510:'Tailandia',511:'Otra Asia Oriental',600:'Americano',601:'Latinoam. 1',602:'Latinoam. 2',603:'Nativo Americano',700:'Oceania',800:'Multiascendencia'};
  return m[id]||'';
}
function drawEFAtriangle(anc){
  const svg=$('#efa'); svg.innerHTML='';
  const E=[1.4758,1.4370], F=[1.0800,1.1000], A=[1.7045,1.1000], S=[anc.gd1,anc.gd2];
  const xmin=0.95,xmax=1.85,ymin=0.95,ymax=1.62;
  const W=460,H=400;
  const sx=x=>(x-xmin)/(xmax-xmin)*W, sy=y=>H-(y-ymin)/(ymax-ymin)*H;
  const ns='http://www.w3.org/2000/svg';
  function pt(p){ return sx(p[0]).toFixed(1)+','+sy(p[1]).toFixed(1); }
  const tri=document.createElementNS(ns,'polygon');
  tri.setAttribute('points',pt(E)+' '+pt(F)+' '+pt(A));
  tri.setAttribute('fill','#eef3f8'); tri.setAttribute('stroke','#5b6b7a');
  svg.appendChild(tri);
  function dot(p,color,r,label,dx,dy){
    const c=document.createElementNS(ns,'circle'); c.setAttribute('cx',sx(p[0])); c.setAttribute('cy',sy(p[1])); c.setAttribute('r',r); c.setAttribute('fill',color); svg.appendChild(c);
    if(label){ const tx=document.createElementNS(ns,'text'); tx.setAttribute('x',sx(p[0])+dx); tx.setAttribute('y',sy(p[1])+dy); tx.setAttribute('font-size',11); tx.setAttribute('fill',color); tx.textContent=label; svg.appendChild(tx); }
  }
  dot(E,'#2f6fb2',8,'Europeo',6,16); dot(F,'#8a5a2b',8,'Africano',6,16); dot(A,'#b23f3f',8,'E.Asiatico',6,16);
  dot(S,'#0f9d58',12,'Tu muestra',-30,26);
}

/* ============================ INIT ============================ */
window.addEventListener('DOMContentLoaded',()=>{
  const input=$('#file'); const zone=$('#drop');
  input.addEventListener('change',e=>handleFile(e.target.files[0]).catch(err=>{ setStatus('Error: '+err.message,0); }));
  ['dragover','dragenter'].forEach(ev=>zone.addEventListener(ev,e=>{e.preventDefault();zone.classList.add('over');}));
  ['dragleave','drop'].forEach(ev=>zone.addEventListener(ev,e=>{e.preventDefault();zone.classList.remove('over');}));
  zone.addEventListener('drop',e=>{ const f=e.dataTransfer.files[0]; handleFile(f).catch(err=>{ setStatus('Error: '+err.message,0); }); });
});
