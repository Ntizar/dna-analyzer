/* DNA Analyzer — datos: rasgos explicados en lenguaje humano + poblaciones de referencia para el mapa. */
'use strict';

/* Cada rasgo tiene: bio (que es, en cristiano), interp (que significa PARA TI) y
   opcionalmente consejo (accion practica). Todo en lenguaje llano. */
const TRAITS = [
  { rs:'4988235', name:'Tolerancia a la lactosa', gene:'MCM6 (LCT)', effect:'T', cat:'Nutrición',
    bio:'La lactasa es la enzima que digiere el azúcar de la leche (lactosa). Muchas personas dejan de producirla al crecer; el alelo T la mantiene activa toda la vida.',
    interp:g=>{const t=(g.match(/T/g)||[]).length;
      return t>=1?['✓ Tu cuerpo digiere bien la leche','Portas el alelo T: mantienes la enzima lactasa activa de por vida, así que la leche no te da problemas.']
                 :['⚠️ Tu cuerpo probablemente digiere mal la leche','No portas el alelo T. Es probable que tengas intolerancia a la lactosa: al tomar lácteos puedes notar hinchazón o molestias. (Este sitio es complejo, no es 100% concluyente.)'];}},
  { rs:'762551', name:'Metabolismo de la cafeína', gene:'CYP1A2', effect:'A', cat:'Metabolismo',
    bio:'CYP1A2 es la enzima que elimina la cafeína de tu cuerpo. Hay versiones más rápidas y más lentas.',
    interp:g=>{const a=(g.match(/A/g)||[]).length;
      return a===2?['⚡ Eliminas la cafeína muy rápido','A/A: el café apenas te afecta y no te quita el sueño. Puedes tomarlo tarde sin problema.']
                  :a===1?['Normal: la cafeína te afecta de forma media','A/C: ni muy rápido ni muy lento. El café de tarde puede quitarte el sueño, pero moderado.']
                        :['🐢 La cafeína dura más en tu cuerpo','C/C: eres más sensible. Mejor no tomar café por la tarde, te costará dormir.'];}},
  { rs:'12913832', name:'Color de ojos', gene:'HERC2', effect:'A', cat:'Rasgo físico',
    bio:'Este gen controla cuánta melanina (pigmento) hay en el iris. Menos melanina = ojos más claros.',
    interp:g=>{const a=(g.match(/A/g)||[]).length;
      return a===2?['👁️ Tendencia a ojos claros (azul/verde)','A/A: genéticamente muy probable que tengas ojos azules o verdes. (El color exacto lo deciden varios genes, pero este es de los que más mandan.)']
                  :a===1?['👁️ Color de ojos mixto','A/G.']:['👁️ Tendencia a ojos oscuros','G/G.'];}},
  { rs:'4680', name:'COMT (Val158Met)', gene:'COMT', effect:'A', cat:'Neurobiología',
    bio:'COMT es la enzima que descompone la dopamina, una sustancia del cerebro ligada al ánimo, la motivación y el placer.',
    interp:g=>{const a=(g.match(/A/g)||[]).length;
      return a===2?['🧠 Versión Met/Met','Tu cerebro mantiene la dopamina más tiempo (degradación lenta).']
                  :a===1?['🧠 Versión intermedia (Val/Met)','Equilibrio de dopamina: ni muy alta ni muy baja.']:['🧠 Versión Val/Val','Degradación rápida de dopamina.'];}},
  { rs:'6265', name:'BDNF (Val66Met)', gene:'BDNF', effect:'T', cat:'Neurobiología',
    bio:'BDNF es una proteína que protege y nutre las neuronas. Ciertas variantes producen algo menos.',
    interp:g=>{const t=(g.match(/T/g)||[]).length;
      return t===2?['🧠 Met/Met','Versión menos común, algo menos de BDNF.']
                  :t===1?['🧠 Val/Met','Intermedio.']:['🧠 Val/Val (la más común)','La variante habitual de esta proteína protectora.'];}},
  { rs:'9939609', name:'Tendencia al peso (FTO)', gene:'FTO', effect:'A', cat:'Metabolismo',
    bio:'FTO influye en el apetito y en cuánta energía gastas. Ciertas variantes dan algo más de tendencia a ganar peso.',
    interp:g=>{const a=(g.match(/A/g)||[]).length;
      return a===2?['⚖️ Dos alelos de riesgo','Algo más de tendencia a subir de peso. La dieta y el ejercicio mandan mucho más que esto.']
                  :a===1?['⚖️ Un alelo de riesgo','Tendencia moderada. Tu peso depende sobre todo de tus hábitos.']:['⚖️ Sin alelo de riesgo','Genotipo de menor tendencia.'];}},
  { rs:'1229984', name:'Metabolismo del alcohol', gene:'ADH1B', effect:'C', cat:'Metabolismo',
    bio:'ADH1B es la primera enzima que procesa el alcohol. La variante Arg48 es más activa: trabaja más rápido.',
    interp:g=>{const c=(g.match(/C/g)||[]).length;
      return c>=1?['🍺 Procesas el alcohol más rápido que la mayoría','Variante Arg48: tu cuerpo descompone el alcohol deprisa, por eso tienes menor riesgo de alcoholismo. Ojo: no significa que beber sea sano, el alcohol sigue dañando.']
                 :['🍺 Metabolismo típico','.'];}},
  { rs:'671', name:'Enrojecimiento por alcohol', gene:'ALDH2', effect:'A', cat:'Metabolismo',
    bio:'ALDH2 es la segunda enzima del alcohol: elimina el acetaldehído, la sustancia que da la cara roja y la "resaca".',
    interp:g=>{const a=(g.match(/A/g)||[]).length;
      return a===0?['✅ No te pones rojo al beber','G/G: tu cuerpo elimina el acetaldehído con normalidad (lo típico en europeos).']
                  :['🔴 Puedes ponerte rojo al beber','Portas la variante ALDH2*2 (más común en asiáticos).'];}},
  { rs:'1801133', name:'MTHFR (folato)', gene:'MTHFR', effect:'A', cat:'Metabolismo',
    bio:'MTHFR ayuda a activar el folato (vitamina B9). Una variante reduce un poco su eficiencia.',
    interp:g=>{const a=(g.match(/A/g)||[]).length;
      return a===2?['677TT','Menos eficiencia con el folato: conviene asegurarlo en la dieta (verduras de hoja, legumbres).']
                  :a===1?['677CT (intermedio)','Una copia de la variante: efecto leve, nada que vigilar de especial.']:['677CC','Eficiencia normal del folato.'];}},
  { rs:'1544410', name:'Vitamina D (VDR)', gene:'VDR', effect:'T', cat:'Metabolismo',
    bio:'VDR es el receptor que "recibe" la vitamina D. Ciertas variantes afectan un poco a cómo la usa tu cuerpo.',
    interp:g=>{const t=(g.match(/T/g)||[]).length;
      return t===2?['BsmI TT','Algo de variación en el uso de vitamina D.']:t===1?['Heterocigoto (C/T)','Efecto leve. Mantener la vitamina D normal es suficiente.']:['BsmI CC','.'];}},
  { rs:'1042713', name:'Receptor β2 (ADRB2)', gene:'ADRB2', effect:'A', cat:'Fisiología',
    bio:'ADRB2 es un receptor que se encuentra en los pulmones y el corazón. Ciertas variantes cambian cómo responde a los medicamentos para el asma.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===2?'Arg16/Arg16':a===1?'Arg16/Gly16':'Gly16/Gly16','Variante común, sin relevancia práctica para ti.'];}},
  { rs:'3827760', name:'Pelo grueso / incisivos', gene:'EDAR', effect:'A', cat:'Rasgo físico',
    bio:'EDAR influye en el grosor del pelo y la forma de los incisivos. La variante A es más típica de personas con ascendencia asiática.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===2?'Tienes la variante A/A (poco habitual en europeos)':'Curioso: esta variante se asocia a pelo más grueso e incisivos "en pala". Es inusual en europeos, así que es un detalle llamativo de tu ADN.'];}},
  { rs:'713598', name:'Gusto amargo (PTC)', gene:'TAS2R38', effect:'C', cat:'Sentido',
    bio:'TAS2R38 detecta el amargor. Según tus alelos, ciertos alimentos amargos (brócoli, café, coles) te saben muy amargos o casi nada.',
    interp:g=>{return ['Genotipo '+g,'Heterocigoto: percibes el amargor de forma intermedia (ni super sensible ni insensible).'];}},
  { rs:'1805008', name:'Pelo rojo', gene:'MC1R', effect:'T', cat:'Rasgo físico',
    bio:'MC1R decide entre producir pigmento rojo (pelirrojo) o oscuro. Ciertas variantes dan más tonos rojizos.',
    interp:g=>{const t=(g.match(/T/g)||[]).length; return [t===0?'No tienes la variante de pelo rojo':'C/C: tu pelo no tiende al rojo por este gen.'];}},
  { rs:'1426654', name:'Piel clara', gene:'SLC24A5', effect:'A', cat:'Rasgo físico',
    bio:'SLC24A5 es uno de los genes que decide lo clara u oscura que es tu piel. La variante A es la típica europea de piel clara.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===2?'Tienes la variante europea de piel clara':'A/A: el gen de la piel clara, el más frecuente en poblaciones de Europa.'];}},
  { rs:'1800414', name:'Pigmentación (SLC45A2)', gene:'SLC45A2', effect:'A', cat:'Rasgo físico',
    bio:'Otro gen de la pigmentación. Junto con el anterior, ayuda a dar tonos de piel más claros.',
    interp:g=>{return ['Genotipo '+g,'Variante de pigmentación, junto con SLC24A5.'];}},
  { rs:'12203592', name:'Piel / pecas (IRF4)', gene:'IRF4', effect:'T', cat:'Rasgo físico',
    bio:'IRF4 influye en el color del pelo, la piel y en la tendencia a tener pecas.',
    interp:g=>{const t=(g.match(/T/g)||[]).length; return [t>=1?'Tienes el alelo T (más pecas / piel clara)':'Se asocia a piel clara y a más facilidad para las pecas.'];}},
  { rs:'1799945', name:'Hemocromatosis (HFE H63D) — portador', gene:'HFE', effect:'G', cat:'Salud (portador)',
    bio:'HFE controla cuánto hierro absorbes. La variante H63D es común. Ser portador (una sola copia) NO da la enfermedad; solo dos copias (o una + C282Y) pueden acumular demasiado hierro.',
    interp:g=>{const gg=(g.match(/G/g)||[]).length;
      return gg>=1?['⚠️ Eres portador de H63D (no estás enfermo)','Muy común (~25% de europeos). Llevas una copia: no te afecta. Solo importaría si tu pareja también fuera portadora. No es un diagnóstico.']:['Sin la variante H63D','.'];}},
  { rs:'3892097', name:'CYP2D6 (metabolismo de fármacos)', gene:'CYP2D6', effect:'T', cat:'Farmacogenética',
    bio:'CYP2D6 metaboliza muchos medicamentos. La variante *4 hace que algunos fármacos se eliminen más lento.',
    interp:g=>{const t=(g.match(/T/g)||[]).length; return [t===0?'No portas la variante *4 (metabolismo normal en este marcador)':'Solo es un marcador entre varios de CYP2D6; no es una evaluación completa.'];}},
  { rs:'4986893', name:'CYP2C19 (metabolismo de fármacos)', gene:'CYP2C19', effect:'A', cat:'Farmacogenética',
    bio:'CYP2C19 metaboliza fármacos como el clopidogrel. La variante *3 lo hace más lento.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===0?'No portas la variante *3 (metabolismo normal en este marcador)':'Solo un marcador; no es una evaluación completa.'];}},
  { rs:'17822931', name:'Cerumen seco/húmedo', gene:'ABCC11', effect:'A', cat:'Rasgo físico',
    bio:'ABCC11 decide si tu cera del oído es seca o húmeda, y también influye en el olor corporal. La variante ancestral da cera seca y menos olor.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===0?'Probablemente tienes cera seca (y menos olor corporal)':'C/C: variante ancestral, lo típico en europeos y asiáticos orientales.'];}},
  { rs:'1051730', name:'Tabaco (cantidad de cigarrillos)', gene:'CHRNA3', effect:'A', cat:'Hábito',
    bio:'CHRNA3 influye en cuántos cigarrillos fuma una persona y en lo difícil que es dejarlo.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a>=1?'Portas el alelo A (algo más de tendencia a fumar más)':'Si fumas, este genotipo puede asociarse a más cigarrillos al día. No fumar es lo mejor.'];}},
  { rs:'1800497', name:'Receptor de dopamina (DRD2)', gene:'DRD2/ANKK1', effect:'A', cat:'Neurobiología',
    bio:'DRD2 es el receptor de la recompensa. El alelo A1 se asocia a menos receptores y a más búsqueda de recompensa (comida, tabaco, etc.).',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===0?'No portas el alelo A1 (A2/A2)':'G/G: más receptores D2, respuesta de recompensa normal.'];}},
  { rs:'1799971', name:'Sensibilidad a opioides', gene:'OPRM1', effect:'G', cat:'Farmacogenética',
    bio:'OPRM1 participa en cómo sientes el dolor y cómo respondes a los analgésicos opioides. La variante G aumenta la sensibilidad.',
    interp:g=>{const g2=(g.match(/G/g)||[]).length; return [g2===0?'No portas el alelo G (A118G)':'A/A: sensibilidad al dolor y respuesta a opioides típicas. Solo es un marcador.'];}},
  { rs:'1805009', name:'Pelo rojo (MC1R)', gene:'MC1R', effect:'A', cat:'Rasgo físico',
    bio:'Otra variante de MC1R (R163Q) ligada a tonos de pelo más claros o rojizos.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===0?'No portas la variante R163Q':'G/G: sin la variante de tonos rojizos.'];}},
  { rs:'1805006', name:'Pelo rojo (MC1R)', gene:'MC1R', effect:'T', cat:'Rasgo físico',
    bio:'Variante D294H de MC1R, que reduce la pigmentación oscura.',
    interp:g=>{const t=(g.match(/T/g)||[]).length; return [t===0?'No portas la variante D294H':'C/C: sin la variante de pelo rojo.'];}},
  { rs:'1726866', name:'Gusto amargo (TAS2R38)', gene:'TAS2R38', effect:'A', cat:'Sentido',
    bio:'Segundo marcador del receptor del gusto amargo.',
    interp:g=>{return ['Genotipo '+g,'Percepción intermedia del amargor.'];}},
  { rs:'10246939', name:'Gusto amargo (TAS2R38)', gene:'TAS2R38', effect:'T', cat:'Sentido',
    bio:'Tercer marcador del gusto amargo.',
    interp:g=>{return ['Genotipo '+g,'Percepción intermedia del amargor.'];}},
  { rs:'731236', name:'Vitamina D (VDR TaqI)', gene:'VDR', effect:'A', cat:'Metabolismo',
    bio:'Otra variante del receptor de vitamina D.',
    interp:g=>{return ['Genotipo '+g,'Variación leve en el uso de vitamina D.'];}},
  { rs:'2073618', name:'Vitamina D (VDR)', gene:'VDR', effect:'A', cat:'Metabolismo',
    bio:'Variante FokI del receptor de vitamina D.',
    interp:g=>{return ['Genotipo '+g,'Variación leve en el uso de vitamina D.'];}},
  { rs:'7294', name:'Vitamina D (VDR)', gene:'VDR', effect:'T', cat:'Metabolismo',
    bio:'Marcador adicional del receptor de vitamina D.',
    interp:g=>{return ['Genotipo '+g,'Variación leve en el uso de vitamina D.'];}},
  { rs:'1799963', name:'Protrombina (F2 G20210A)', gene:'F2', effect:'A', cat:'Salud (trombofilia)',
    bio:'La variante G20210A de la protrombina aumenta algo el riesgo de coágulos (trombosis). Ser portador da un riesgo ligeramente mayor.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===0?'No portas la variante G20210A':'G/G: genotipo común, sin la variante. No es un diagnóstico.'];}},
  { rs:'6025', name:'Factor V Leiden (trombofilia)', gene:'F5', effect:'A', cat:'Salud (trombofilia)',
    bio:'El Factor V Leiden es la causa genética de trombosis más frecuente en europeos. Es un marcador de riesgo, NO un diagnóstico.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===0?'No se detecta el alelo de riesgo Leiden':'C/C. Este marcador necesita interpretación clínica; consúltalo con un profesional.'];}},
  { rs:'2108622', name:'CYP4F2 (metabolismo de fármacos)', gene:'CYP4F2', effect:'T', cat:'Farmacogenética',
    bio:'CYP4F2 procesa vitamina K y la warfarina (un anticoagulante). La variante T lo hace más lento.',
    interp:g=>{return ['Genotipo '+g,'Relevante solo para ajustar la dosis de warfarina. Un marcador más.'];}},
  { rs:'1042714', name:'Receptor β2 (Gln27Glu)', gene:'ADRB2', effect:'G', cat:'Fisiología',
    bio:'Otra variante del receptor beta-2, relacionada con la función pulmonar y la respuesta a ciertos fármacos.',
    interp:g=>{return ['Genotipo '+g,'Variante común de ADRB2, sin relevancia práctica para ti.'];}},
  { rs:'1801394', name:'MTRR (vitamina B12)', gene:'MTRR', effect:'A', cat:'Metabolismo',
    bio:'MTRR participa en reciclar la vitamina B12 y el folato.',
    interp:g=>{return ['Genotipo '+g,'Variante del metabolismo de B12/folato. Mantén una dieta variada.'];}},
];

/* Ubicaciones geograficas de poblaciones de referencia (lat, lon, grupo continental, color) para el mapa. */
const REF_POPS = [
  { name:'España / Portugal (IBS)', lat:40.4, lon:-3.7, grp:'EUR', color:'#2f6fb2' },
  { name:'Italia (TSI)', lat:41.9, lon:12.5, grp:'EUR', color:'#2f6fb2' },
  { name:'Francia', lat:46.6, lon:2.4, grp:'EUR', color:'#2f6fb2' },
  { name:'Irlanda / UK (CEU/GBR)', lat:53.5, lon:-7.0, grp:'EUR', color:'#2f6fb2' },
  { name:'Alemania', lat:51.2, lon:10.4, grp:'EUR', color:'#2f6fb2' },
  { name:'Finlandia (FIN)', lat:62.0, lon:25.0, grp:'EUR', color:'#2f6fb2' },
  { name:'Polonia', lat:52.1, lon:19.4, grp:'EUR', color:'#2f6fb2' },
  { name:'Nigeria (YRI)', lat:7.6, lon:4.5, grp:'AFR', color:'#8a5a2b' },
  { name:'Ghana', lat:5.6, lon:-1.2, grp:'AFR', color:'#8a5a2b' },
  { name:'Irán', lat:32.4, lon:53.7, grp:'MEN', color:'#a07c2b' },
  { name:'Marruecos / N. África', lat:31.8, lon:-7.1, grp:'MEN', color:'#a07c2b' },
  { name:'China (CHB)', lat:35.0, lon:103.0, grp:'EAS', color:'#b23f3f' },
  { name:'Japón (JPT)', lat:36.2, lon:138.0, grp:'EAS', color:'#b23f3f' },
  { name:'India (GIH)', lat:22.0, lon:79.0, grp:'SAS', color:'#7a4fb2' },
  { name:'Pakistán', lat:30.0, lon:70.0, grp:'SAS', color:'#7a4fb2' },
  { name:'México / LATAM (MXL)', lat:19.4, lon:-99.1, grp:'AMR', color:'#2f8f6b' },
];

/* Zonas / regiones geograficas para el mapa (nombre, centro, radio en grados, color, grupo). */
const REGIONS = [
  { name:'Península Ibérica', lat:40.4, lon:-3.7, r:3.6, color:'#2f6fb2', grp:'EUR' },
  { name:'Europa Occidental', lat:47.5, lon:2.5, r:5.0, color:'#3a7fc0', grp:'EUR' },
  { name:'Italia / Sur de Europa', lat:42.0, lon:12.5, r:3.0, color:'#2f6fb2', grp:'EUR' },
  { name:'Norte de Europa', lat:55.5, lon:12.0, r:6.0, color:'#5a8fca', grp:'EUR' },
  { name:'Norte de África', lat:32.0, lon:-5.0, r:4.0, color:'#a07c2b', grp:'MEN' },
  { name:'África Occidental', lat:8.0, lon:-5.0, r:5.0, color:'#8a5a2b', grp:'AFR' },
  { name:'Este de Asia', lat:35.0, lon:110.0, r:9.0, color:'#b23f3f', grp:'EAS' },
  { name:'Sur de Asia', lat:22.0, lon:78.0, r:7.0, color:'#7a4fb2', grp:'SAS' },
  { name:'Latinoamérica', lat:19.0, lon:-99.0, r:6.0, color:'#2f8f6b', grp:'AMR' },
];
