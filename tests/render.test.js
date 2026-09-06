// Verifica el render de la UI con un DOM simulado (sin navegador).
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const appSrc = fs.readFileSync(path.join(__dirname, 'js', 'app.js'), 'utf8');
const panelBytes = fs.readFileSync(path.join(__dirname, 'data', 'panel.bin.gz'));

// --- DOM simulado minimal ---
function el(tag){ return {
  tag, children:[], textContent:'', innerHTML:'', style:{}, attrs:{},
  setAttribute(k,v){ this.attrs[k]=v; },
  appendChild(c){ this.children.push(c); return c; },
}; }
const nodes = {};
const documentMock = {
  querySelector(sel){ if(!nodes[sel]) nodes[sel]=el(sel); return nodes[sel]; },
  createElement(tag){ return el(tag); },
  createElementNS(ns, tag){ return el(tag); },
};
const sandbox = {
  console, Math, Map, Uint8Array, Uint32Array, Float64Array, String, Object, Array,
  TextDecoder, Promise, Response, DecompressionStream,
  window:{ addEventListener(){} },
  document: documentMock,
  JSZip:null,
  fetch: async (u)=>{ if(String(u).includes('panel.bin.gz')) return new Response(panelBytes); throw new Error('no route'); },
};
vm.createContext(sandbox);
vm.runInContext(appSrc, sandbox);

(async () => {
  await sandbox.loadPanel();
  const text = fs.readFileSync('C:/Users/d_ant/dna_work/MyHeritage_raw_dna_data.csv','utf8');
  const parsed = sandbox.parseRaw(text);
  const anc = sandbox.computeAncestry(parsed.geno);
  const traits = sandbox.computeTraits(parsed.geno);
  // Ejecutar el render (mock DOM) — debe no lanzar
  sandbox.render(parsed, anc, traits);
  console.log('render() OK — sin errores');
  console.log('resultados display:', nodes['#results'] && nodes['#results'].style.display);
  console.log('sexo:', nodes['#sexo'] && nodes['#sexo'].textContent);
  console.log('anc-box tiene contenido:', (nodes['#anc-box'] && nodes['#anc-box'].innerHTML.length>0));
  console.log('traits cards:', nodes['#traits'] && nodes['#traits'].children.length);
  console.log('efa svg children (vertices+muestra):', nodes['#efa'] && nodes['#efa'].children.length);
  console.log('resumen:', nodes['#resumen'] && nodes['#resumen'].textContent);
})();
