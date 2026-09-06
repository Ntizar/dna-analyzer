// Test del render con DOM simulado (incluye mapa, verificacion, rasgos enriquecidos).
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const dataSrc = fs.readFileSync(path.join(__dirname, '..', 'js', 'data.js'), 'utf8');
const appSrc  = fs.readFileSync(path.join(__dirname, '..', 'js', 'app.js'), 'utf8');
const panelBytes = fs.readFileSync(path.join(__dirname, '..', 'data', 'panel.bin.gz'));

function el(tag){ return { tag, children:[], textContent:'', innerHTML:'', style:{}, attrs:{},
  setAttribute(k,v){ this.attrs[k]=v; }, appendChild(c){ this.children.push(c); return c; } }; }
const nodes = {};
const documentMock = {
  querySelector(sel){ if(!nodes[sel]) nodes[sel]=el(sel); return nodes[sel]; },
  createElement(t){ return el(t); }, createElementNS(ns,t){ return el(t); },
};
const sandbox = {
  console, Math, Map, Uint8Array, Uint32Array, Float64Array, String, Object, Array,
  TextDecoder, Promise, Response, DecompressionStream,
  window:{ addEventListener(){} }, document: documentMock, JSZip:null, L:undefined,
  fetch: async (u)=>{ if(String(u).includes('panel.bin.gz')) return new Response(panelBytes); throw new Error('no route'); },
};
vm.createContext(sandbox);
vm.runInContext(dataSrc + '\n' + appSrc, sandbox);

(async () => {
  await sandbox.loadPanel();
  const text = fs.readFileSync('C:/Users/d_ant/dna_work/MyHeritage_raw_dna_data.csv','utf8');
  const parsed = sandbox.parseRaw(text);
  const anc = sandbox.computeAncestry(parsed.geno);
  const traits = sandbox.computeTraits(parsed.geno);
  const verif = sandbox.computeVerification(parsed.geno);
  sandbox.render(parsed, anc, traits, verif);   // no debe lanzar
  console.log('render OK, sin errores');
  console.log('results display:', nodes['#results'] && nodes['#results'].style.display);
  console.log('sexo:', nodes['#sexo'] && nodes['#sexo'].textContent);
  console.log('anc-box len>0:', (nodes['#anc-box'] && nodes['#anc-box'].innerHTML.length>0));
  console.log('verif-box tiene HTML:', (nodes['#verif-box'] && nodes['#verif-box'].innerHTML.length>0));
  console.log('traits cards:', nodes['#traits'] && nodes['#traits'].children.length);
  console.log('efa children:', nodes['#efa'] && nodes['#efa'].children.length);
  console.log('map div existe (render no fallo):', !!nodes['#map']);
  console.log('resumen:', nodes['#resumen'] && nodes['#resumen'].textContent);
})();
