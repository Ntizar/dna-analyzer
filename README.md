# 🧬 DNA Analyzer

Analiza tu ADN crudo (MyHeritage / 23andMe / Ancestry) **100% en tu navegador**, sin que los datos salgan de tu máquina.

> 🔒 **Privacidad total.** Tu fichero se descomprime y procesa localmente en tu navegador. No hay servidor, ni subida, ni almacenamiento. Solo se descarga un panel genético de referencia (anónimo, no contiene datos de personas).

## Qué calcula

- **Ascendencia continental** (método **GrafAnc / GRAF-pop**, NCBI — *Jin et al., HGG Advances 2025*): porcentaje Europeo / Africano / Este-asiático, gráfico del triángulo E-F-A.
- **Grupo subcontinental** (p. ej. Europa Occidental / Meridional).
- **Sexo** (por homocigosidad del cromosoma X).
- **Rasgos genéticos** (curiosidad, no médico): lactosa, cafeína, color de ojos, COMT, BDNF, FTO, alcohol, MTHFR…

## Cómo usar

1. Descarga tu ADN crudo desde **MyHeritage**, **23andMe** o **Ancestry** (formato `.zip` o `.csv`).
2. Arrastra el fichero sobre la página o selecciónalo.
3. El análisis se ejecuta en tu navegador y se muestran los resultados al instante.

## Método

El núcleo es el algoritmo **GRAF-pop/GrafAnc**: compara los genotipos de ~71.000 SNP de ascendencia con las frecuencias alélicas de **26 poblaciones de referencia** (UK Biobank + 1000 Genomas + HGDP) y calcula distancias genéticas y proporciones de mezcla mediante coordenadas baricéntricas sobre el triángulo Europeo–Africano–Este-asiático. Se reimplementó en JavaScript puro (sin dependencias) y se validó contra la implementación de referencia (mismo resultado: GD1≈1.472, GD2≈1.434, EUR 98.98% para un europeo).

## Estructura

```
dna-analyzer/
├── index.html          # interfaz
├── css/style.css       # estilos
├── js/app.js           # lógica (parseo, GrafAnc, rasgos, UI)
├── js/jszip.min.js     # descompresión de ZIP (vendored)
├── data/panel.bin.gz   # panel de referencia compacto (~16 MB)
└── tests/              # tests de Node (algoritmo + render)
```

## Despliegue

GitHub Pages. El panel de referencia se sirve como estático y se descarga una vez (se cachea).

---

Hecho con ❤️ por **David Antizar**. Método: GrafAnc/GRAF-pop (NCBI). No es un consejo médico.
