// Harness: ejecuta app.js en un contexto vm y llama a las funciones reales.
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const appSrc = fs.readFileSync(path.join(__dirname, 'js', 'app.js'), 'utf8');
const panelBytes = fs.readFileSync(path.join(__dirname, 'data', 'panel.bin.gz'));

const sandbox = {
  console, Math, Map, Uint8Array, Uint32Array, Float64Array, String, Object, Array,
  TextDecoder, Promise,
  Response, DecompressionStream,
  window: { addEventListener(){} },
  document: { querySelector(){ return null; } },
  JSZip: null,
  fetch: async (url) => { if(String(url).includes('panel.bin.gz')) return new Response(panelBytes); throw new Error('no route '+url); },
};
vm.createContext(sandbox);
vm.runInContext(appSrc, sandbox);

(async () => {
  await sandbox.loadPanel();
  const text = fs.readFileSync('C:/Users/d_ant/dna_work/MyHeritage_raw_dna_data.csv', 'utf8');
  const parsed = sandbox.parseRaw(text);
  console.log('SNP parseados:', parsed.geno.size, '| sexo:', parsed.sex, '| xHet:', parsed.xHet);
  const anc = sandbox.computeAncestry(parsed.geno);
  console.log('=== RESULTADO JS (GrafAnc) ===');
  console.log('numGeno:', anc.numGeno);
  console.log('GD1/GD2/GD3:', anc.gd1.toFixed(6), anc.gd2.toFixed(6), anc.gd3.toFixed(6));
  console.log('EUR/AFR/EAS:', anc.ePct.toFixed(2)+'%', anc.fPct.toFixed(2)+'%', anc.aPct.toFixed(2)+'%');
  console.log('AncGroupID:', anc.ancGroupId);
  console.log('subpop:', anc.subPopScores.map(x=>x.toFixed(3)).join(', '));
  const traits = sandbox.computeTraits(parsed.geno);
  console.log('=== RASGOS ===');
  traits.forEach(t=>console.log(' -', t.name, '|', t.genotype, '|', t.label, '|', t.detail||''));
})();
