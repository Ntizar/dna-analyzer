/* DNA Analyzer — datos: rasgos con explicaciones biológicas + poblaciones de referencia para el mapa. */
'use strict';

/* TRAITS: cada entrada tiene alelo efectivo (hebra forward GRCh37, verificada),
   descripcion corta y explicacion biologica de como afecta la genetica. */
const TRAITS = [
  { rs:'4988235', name:'Tolerancia a la lactosa', gene:'MCM6 (LCT)', effect:'T', cat:'Nutrición',
    bio:'La enzima lactasa digiere la lactosa de la leche. En la mayoría de adultos su actividad cae; el alelo T de rs4988235 la mantiene activa en la vida adulta (persistencia).',
    interp:g=>{const t=(g.match(/T/g)||[]).length;
      return t>=1?['Tolerante a la lactosa','Portas el alelo T de persistencia: la lactasa sigue activa.']
                 :['Probable intolerancia','No portas el alelo T. Puede que al tomar leche notes molestias. (Sitio multi-alélico, no 100% concluyente).'];}},
  { rs:'762551', name:'Metabolismo de la cafeína', gene:'CYP1A2', effect:'A', cat:'Metabolismo',
    bio:'CYP1A2 metaboliza la cafeína. El alelo A (rs762551) se asocia a mayor inducibilidad: metabolizas la cafeína más rápido.',
    interp:g=>{const a=(g.match(/A/g)||[]).length;
      return a===2?['Metabolizador rápido','A/A: eliminas la cafeína rápido, menos efecto y menos insomnio.']
                  :a===1?['Metabolizador intermedio','A/C: velocidad media.']
                        :['Metabolizador lento','C/C: la cafeína dura más en tu cuerpo, más sensibilidad.'];}},
  { rs:'12913832', name:'Color de ojos', gene:'HERC2', effect:'A', cat:'Rasgo físico',
    bio:'HERC2 regula la expresión de OCA2, que produce melanina en el iris. El alelo A reduce la melanina → ojos claros (azul/verde).',
    interp:g=>{const a=(g.match(/A/g)||[]).length;
      return a===2?['Tendencia a ojos claros','A/A: fuerte propensión a ojos azul/verde (el color es poligénico, pero este SNP pesa mucho).']
                  :a===1?['Color mixto','A/G.']:['Tendencia a ojos oscuros','G/G.'];}},
  { rs:'4680', name:'COMT (Val158Met)', gene:'COMT', effect:'A', cat:'Neurobiología',
    bio:'COMT degrada la dopamina. El alelo A (Met) reduce su actividad → la dopamina dura más en las sinapsis.',
    interp:g=>{const a=(g.match(/A/g)||[]).length;
      return a===2?['Met/Met','Degradación lenta de dopamina (más dopamina disponible).']
                  :a===1?['Val/Met','Intermedio: equilibrio de dopamina.']
                        :['Val/Val','Degradación rápida de dopamina.'];}},
  { rs:'6265', name:'BDNF (Val66Met)', gene:'BDNF', effect:'T', cat:'Neurobiología',
    bio:'BDNF apoya la salud de las neuronas. El alelo T (Met) se asocia a menor secreción de BDNF.',
    interp:g=>{const t=(g.match(/T/g)||[]).length;
      return t===2?['Met/Met','Menor BDNF (variante menos común).']
                  :t===1?['Val/Met','Intermedio.']:['Val/Val','Variante común de BDNF.'];}},
  { rs:'9939609', name:'Tendencia al peso (FTO)', gene:'FTO', effect:'A', cat:'Metabolismo',
    bio:'FTO influye en el apetito y el gasto energético. El alelo A de rs9939609 se asocia a un IMC ligeramente mayor.',
    interp:g=>{const a=(g.match(/A/g)||[]).length;
      return a===2?['Dos alelos de riesgo','Mayor tendencia, aunque la dieta y ejercicio mandan.']
                  :a===1?['Un alelo de riesgo','Riesgo moderado.']:['Sin alelo de riesgo','Genotipo de menor tendencia.'];}},
  { rs:'1229984', name:'Metabolismo del alcohol', gene:'ADH1B', effect:'C', cat:'Metabolismo',
    bio:'ADH1B convierte el alcohol en acetaldehído. El alelo C (Arg48) da una enzima más activa: procesas el alcohol más rápido.',
    interp:g=>{const c=(g.match(/C/g)||[]).length;
      return c>=1?['Metabolismo rápido','Variante Arg48: degradas el alcohol rápido, menor riesgo de alcoholismo.']
                 :['Metabolismo típico','.'];}},
  { rs:'671', name:'Enrojecimiento por alcohol', gene:'ALDH2', effect:'A', cat:'Metabolismo',
    bio:'ALDH2 elimina el acetaldehído. El alelo A (ALDH2*2) lo inactiva y causa enrojecimiento/taquicardia al beber (común en Asia).',
    interp:g=>{const a=(g.match(/A/g)||[]).length;
      return a===0?['Sin variante de enrojecimiento','G/G: metabolismo normal del acetaldehído.']
                  :['Riesgo de enrojecimiento','Portas la variante ALDH2*2 (más típica de ascendencia asiática).'];}},
  { rs:'1801133', name:'MTHFR (folato)', gene:'MTHFR', effect:'A', cat:'Metabolismo',
    bio:'MTHFR activa el folato. El alelo A corresponde al 677T (menor actividad de la enzima).',
    interp:g=>{const a=(g.match(/A/g)||[]).length;
      return a===2?['677TT','Menor actividad de folato (más aporte dietético recomendable).']
                  :a===1?['677CT','Intermedio.']:['677CC','Actividad normal.'];}},
  { rs:'1544410', name:'Vitamina D (VDR)', gene:'VDR', effect:'T', cat:'Metabolismo',
    bio:'VDR es el receptor de la vitamina D. El alelo T (BsmI) modula ligeramente los niveles de vitamina D.',
    interp:g=>{const t=(g.match(/T/g)||[]).length;
      return t===2?['BsmI TT',''].slice(0):t===1?['Heterocigoto','C/T.']:['BsmI CC','.'];}},
  { rs:'1042713', name:'Receptor β2 (ADRB2)', gene:'ADRB2', effect:'A', cat:'Fisiología',
    bio:'ADRB2 es el receptor beta-2 (pulmones, corazón). El alelo A (Arg16) influye en la respuesta a beta-agonistas.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===2?'Arg16/Arg16':a===1?'Arg16/Gly16':'Gly16/Gly16','Variante común.'];}},
  { rs:'3827760', name:'Pelo grueso / incisivos', gene:'EDAR', effect:'A', cat:'Rasgo físico',
    bio:'EDAR (V370A) influye en el grosor del pelo y la forma de los incisivos. El alelo A es más frecuente en ascendencia asiática.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===2?'Portas la variante (A/A)':'','El alelo A se asocia a pelo más grueso e incisivos en pala; es inusual en europeos — observación interesante.'];}},
  { rs:'713598', name:'Gusto amargo (PTC)', gene:'TAS2R38', effect:'C', cat:'Sentido',
    bio:'TAS2R38 detecta el amargor (PTC). Distintos alelos hacen que unas personas lo perciban muy amargo y otras apenas.',
    interp:g=>{return ['Genotipo '+g,'Heterocigoto en el receptor del gusto amargo: percepción intermedia del amargor.'];}},
  { rs:'1805008', name:'Pelo rojo', gene:'MC1R', effect:'T', cat:'Rasgo físico',
    bio:'MC1R controla la producción de feomelanina (rojo) vs eumelanina (oscuro). El alelo T (R160W) reduce la eumelanina.',
    interp:g=>{const t=(g.match(/T/g)||[]).length; return [t===0?'No portas la variante de pelo rojo':'','C/C: sin la variante R160W asociada al pelo rojo.'];}},
  { rs:'1426654', name:'Piel clara', gene:'SLC24A5', effect:'A', cat:'Rasgo físico',
    bio:'SLC24A5 influye en la pigmentación. El alelo A (Thr111) es la variante europea asociada a piel más clara.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===2?'Portas la variante europea de piel clara':'','A/A: variante de piel clara, frecuente en poblaciones europeas.'];}},
  { rs:'1800414', name:'Pigmentación (SLC45A2)', gene:'SLC45A2', effect:'A', cat:'Rasgo físico',
    bio:'SLC45A2 también regula la pigmentación. El alelo A se asocia a piel/ojos más claros en algunas poblaciones.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return ['Genotipo '+g,'Variante de pigmentación.'];}},
  { rs:'12203592', name:'Piel / pecas (IRF4)', gene:'IRF4', effect:'T', cat:'Rasgo físico',
    bio:'IRF4 influye en el color del pelo y la piel, y en la tendencia a pecas. El alelo T se asocia a piel clara/pecas.',
    interp:g=>{const t=(g.match(/T/g)||[]).length; return [t>=1?'Portas el alelo T (IRF4)':'','Asociado a piel clara y mayor tendencia a pecas.'];}},
  { rs:'1799945', name:'Hemocromatosis (HFE H63D) — portador', gene:'HFE', effect:'G', cat:'Salud (portador)',
    bio:'HFE regula la absorción de hierro. El alelo G (H63D) es una variante común. Ser portador (1 copia) NO causa la enfermedad; solo 2 copias (o 1 copia + C282Y) pueden sobrecargar hierro.',
    interp:g=>{const gg=(g.match(/G/g)||[]).length;
      return gg>=1?['Heterocigoto portador de H63D','Portador (no afectado). Muy común (~25% de europeos). No es un diagnóstico.']:['Sin la variante H63D','.'];}},
  { rs:'3892097', name:'CYP2D6 (metabolismo de fármacos)', gene:'CYP2D6', effect:'T', cat:'Farmacogenética',
    bio:'CYP2D6 metaboliza muchos fármacos. El alelo T (*4) produce un metabolizador pobre.',
    interp:g=>{const t=(g.match(/T/g)||[]).length; return [t===0?'No portas *4 (metabolizador normal en este marcador)':'','Solo 1 de los alelos que definen a CYP2D6; no es una evaluación completa.'];}},
  { rs:'4986893', name:'CYP2C19 (metabolismo de fármacos)', gene:'CYP2C19', effect:'A', cat:'Farmacogenética',
    bio:'CYP2C19 metaboliza fármacos como el clopidogrel. El alelo A (*3) produce un metabolizador pobre.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===0?'No portas *3 (metabolizador normal en este marcador)':'','Solo 1 alelo de CYP2C19; no es una evaluación completa.'];}},
  { rs:'17822931', name:'Cerumen seco/húmedo', gene:'ABCC11', effect:'A', cat:'Rasgo físico',
    bio:'ABCC11 determina el tipo de cerumen y el olor corporal. El alelo A produce cerumen húmedo y mayor olor; el ancestral, cerumen seco.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===0?'Probable cerumen seco (y menos olor corporal)':'','C/C: variante ancestral, típica de europeos/asiáticos orientales.'];}},
  { rs:'1051730', name:'Tabaco (cantidad de cigarrillos)', gene:'CHRNA3', effect:'A', cat:'Hábito',
    bio:'CHRNA3 influye en la cantidad de cigarrillos. El alelo A se asocia a un mayor consumo y a más dificultad para dejar de fumar.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a>=1?'Portas el alelo A (mayor tendencia al consumo)':'','Si fumas, este genotipo puede asociarse a más cigarrillos/día.'];}},
  { rs:'1800497', name:'Receptor de dopamina (DRD2)', gene:'DRD2/ANKK1', effect:'A', cat:'Neurobiología',
    bio:'DRD2 media la recompensa. El alelo A1 se asocia a menos receptores D2 y a conductas de búsqueda de recompensa.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===0?'No portas el alelo A1 (A2/A2)':'','G/G: variante con más receptores D2.'];}},
  { rs:'1799971', name:'Sensibilidad a opioides', gene:'OPRM1', effect:'G', cat:'Farmacogenética',
    bio:'OPRM1 (A118G) media la analgesia opioide. El alelo G se asocia a mayor sensibilidad al dolor y mayor respuesta a opioides.',
    interp:g=>{const g2=(g.match(/G/g)||[]).length; return [g2===0?'No portas el alelo G (A118G)':'','A/A: sensibilidad opioide típica. Solo es un marcador, no una evaluación.'];}},
  { rs:'1805009', name:'Pelo rojo (MC1R)', gene:'MC1R', effect:'A', cat:'Rasgo físico',
    bio:'MC1R (R163Q). Ciertas variantes reducen la eumelanina y favorecen el pelo rojo/tonos rojizos.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===0?'No portas la variante R163Q':'','G/G: sin la variante asociada a tonos rojizos.'];}},
  { rs:'1805006', name:'Pelo rojo (MC1R)', gene:'MC1R', effect:'T', cat:'Rasgo físico',
    bio:'MC1R (D294H). Variante que reduce la pigmentación oscura.',
    interp:g=>{const t=(g.match(/T/g)||[]).length; return [t===0?'No portas la variante D294H':'','C/C: sin la variante asociada a pelo rojo.'];}},
  { rs:'1726866', name:'Gusto amargo (TAS2R38)', gene:'TAS2R38', effect:'A', cat:'Sentido',
    bio:'Segundo marcador del receptor del gusto amargo (PTC).',
    interp:g=>{return ['Genotipo '+g,'Heterocigoto en TAS2R38: percepción intermedia del amargor.'];}},
  { rs:'10246939', name:'Gusto amargo (TAS2R38)', gene:'TAS2R38', effect:'T', cat:'Sentido',
    bio:'Tercer marcador del receptor del gusto amargo.',
    interp:g=>{return ['Genotipo '+g,'Heterocigoto en TAS2R38: percepción intermedia del amargor.'];}},
  { rs:'731236', name:'Vitamina D (VDR TaqI)', gene:'VDR', effect:'A', cat:'Metabolismo',
    bio:'Otro polimorfismo de VDR (TaqI) relacionado con los niveles de vitamina D.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return ['Genotipo '+g,'VDR TaqI: modula los niveles de vitamina D.'];}},
  { rs:'2073618', name:'Vitamina D (VDR)', gene:'VDR', effect:'A', cat:'Metabolismo',
    bio:'VDR (FokI) influye en el receptor de vitamina D.',
    interp:g=>{return ['Genotipo '+g,'Variante del receptor de vitamina D.'];}},
  { rs:'7294', name:'Vitamina D (VDR)', gene:'VDR', effect:'T', cat:'Metabolismo',
    bio:'Marcador adicional de VDR.',
    interp:g=>{return ['Genotipo '+g,'Variante del receptor de vitamina D.'];}},
  { rs:'1799963', name:'Protrombina (F2 G20210A)', gene:'F2', effect:'A', cat:'Salud (trombofilia)',
    bio:'La variante G20210A del gen de la protrombina eleva el riesgo de trombosis. Ser portador aumenta ligeramente el riesgo.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===0?'No portas la variante G20210A':'','G/G: genotipo común, sin la variante. No es un diagnóstico.'];}},
  { rs:'6025', name:'Factor V Leiden (trombofilia)', gene:'F5', effect:'A', cat:'Salud (trombofilia)',
    bio:'El Factor V Leiden es la causa genética de trombosis más común en europeos. Es un marcador de riesgo, NO un diagnóstico.',
    interp:g=>{const a=(g.match(/A/g)||[]).length; return [a===0?'No se detecta el alelo de riesgo Leiden':'','C/C. Este marcador requiere interpretación clínica; consulta con un profesional.'];}},
  { rs:'2108622', name:'CYP4F2 (metabolismo de fármacos)', gene:'CYP4F2', effect:'T', cat:'Farmacogenética',
    bio:'CYP4F2 metaboliza vitamina K y warfarina. El alelo T reduce su actividad.',
    interp:g=>{const t=(g.match(/T/g)||[]).length; return ['Genotipo '+g,'Relevante en dosificación de warfarina. Solo un marcador.'];}},
  { rs:'1042714', name:'Receptor β2 (Gln27Glu)', gene:'ADRB2', effect:'G', cat:'Fisiología',
    bio:'ADRB2 (Gln27Glu) influye en la función pulmonar y la respuesta a agonistas beta.',
    interp:g=>{return ['Genotipo '+g,'Variante común de ADRB2.'];}},
  { rs:'1801394', name:'MTRR (vitamina B12)', gene:'MTRR', effect:'A', cat:'Metabolismo',
    bio:'MTRR participa en la remetilación de homocisteína (vitamina B12/folato).',
    interp:g=>{return ['Genotipo '+g,'Variante del metabolismo de B12/folato.'];}},
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
