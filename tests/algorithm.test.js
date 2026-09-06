// Test: algoritmo GrafAnc + verificacion, ejecutando data.js + app.js en un contexto vm.
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const dataSrc = fs.readFileSync(path.join(__dirname, '..', 'js', 'data.js'), 'utf8');
const appSrc  = fs.readFileSync(path.join(__dirname, '..', 'js', 'app.js'), 'utf8');
const panelBytes = fs.readFileSync(path.join(__dirname, '..', 'data', 'panel.bin.gz'));

const sandbox = {
  console, Math, Map, Uint8Array, Uint32Array, Float64Array, String, Object, Array,
  TextDecoder, Promise, Response, DecompressionStream,
  window: { addEventListener(){} },
  document: { querySelector(){ return null; } },
  JSZip: null, L: undefined,
  fetch: async (u)=>{ if(String(u).includes('panel.bin.gz')) return new Response(panelBytes); throw new Error('no route'); },
};
vm.createContext(sandbox);
vm.runInContext(dataSrc + '\n' + appSrc, sandbox);

(async () => {
  await sandbox.loadPanel();
  const text = fs.readFileSync('C:/Users/d_ant/dna_work/MyHeritage_raw_dna_data.csv','utf8');
  const parsed = sandbox.parseRaw(text);
  const anc = sandbox.computeAncestry(parsed.geno);
  const verif = sandbox.computeVerification(parsed.geno);
  const traits = sandbox.computeTraits(parsed.geno);
  console.log('SNP:', parsed.geno.size, '| sexo:', parsed.sex);
  console.log('=== ASCENDENCIA ===');
  console.log('GD1/GD2/GD3:', anc.gd1.toFixed(4), anc.gd2.toFixed(4), anc.gd3.toFixed(4));
  console.log('EUR/AFR/EAS:', anc.ePct.toFixed(2)+'%', anc.fPct.toFixed(2)+'%', anc.aPct.toFixed(2)+'%', '| grupo', anc.ancGroupId);
  console.log('=== VERIFICACION ===');
  console.log('corr EUR/AFR/EAS:', verif.cEur.toFixed(4), verif.cAfr.toFixed(4), verif.cEas.toFixed(4));
  console.log('sigma GD1/GD2/GD3:', verif.sigmaGD1.toFixed(4), verif.sigmaGD2.toFixed(4), verif.sigmaGD3.toFixed(4));
  console.log('=== RASGOS ('+traits.length+') ===');
  traits.forEach(t=>console.log(' -', t.name, '|', t.genotype, '|', t.label, '|', (t.cat||'')));
})();
