### Importar datos a QIIME2
qiime tools import 
		--type 'SampleData[PairedEndSequencesWithQuality]' 
		--input-path manifest.csv 
		--output-path demux.qza 
		--input-format PairedEndFastqManifestPhred33V2

### Limpieza de primers y filtro de longitud 
	qiime cutadapt trim-paired 
		--i-demultiplexed-sequences demux.qza 
		--p-front-f GCATCGATGAAGAACGCAGC 
		--p-front-r TCCTCCGCTTATTGATATGC 
		--p-minimum-length 100 
		--o-trimmed-sequences trimmed.qza 
		--verbose

### Denoising con dada2
	qiime dada2 denoise-paired
		--i-demultiplexed-seqs trimmed.qza 
		--p-trunc-len-f 220 
		--p-trunc-len-r 160 
		--p-trim-left-f 10 
		--p-trim-left-r 10 
		--p-max-ee-f 3 
		--p-max-ee-r 3 
		--o-table table_test.qza 
		--o-representative-sequences rep-seqs_test.qza 
		--o-denoising-stats denoising-stats_test.qza

### Eliminacion de siingletons
	qiime feature-table filter-features 
		--i-table ../3.denoising_dada2/table_test.qza 
		--p-min-frequency 2 
		--o-filtered-table table_no_singletons.qza

### Clasificacion taxonomica
	qiime feature-classifier classify-sklearn 
		--i-classifier unite_ver10_99_19.02.2025-Q2-2024.10.qza 
		--i-reads ../3.denoising_dada2/rep-seqs_test.qza 
		--o-classification taxonomy.qza

	qiime taxa barplot 
		--i-table ../4.eliminacion_singletons/table_no_singletons.qza 
		--i-taxonomy taxonomy.qza 
		--m-metadata-file ../0.datos_base/metadata.tsv 
		--o-visualization taxa-bar-plots.qzv

### Metricas de diversidad

## Rarefaccion de datos
	qiime feature-table summarize 
		--i-table ../4.eliminacion_singletons/table_no_singletons.qza 
		--o-visualization table2_rarefaccion.qzv 
		--m-sample-metadata-file ../0.datos_base/metadata.tsv

## Obtencion de metricas de diversidad
	qiime diversity core-metrics 
		--i-table ../4.eliminacion_singletons/table_no_singletons.qza 
		--p-sampling-depth 87821 
		--m-metadata-file ../0.datos_base/metadata.tsv 
		--output-dir diversity-metrics-results

## Visualizacion de metricas de alfa diversidad
	qiime diversity alpha-group-significance 
		--i-alpha-diversity observed_features_vector.qza 
		--m-metadata-file ../../0.datos_base/metadata.tsv 
		--o-visualization observed_features_group_significance.qzv

	qiime diversity alpha-group-significance 
		--i-alpha-diversity shannon_vector.qza 
		--m-metadata-file ../../0.datos_base/metadata.tsv 
		--o-visualization shannon-group-significance.qzv

	qiime diversity alpha-group-significance 
		--i-alpha-diversity evenness_vector.qza 
		--m-metadata-file ../../0.datos_base/metadata.tsv 
		--o-visualization evenness_group_significance.qzv

### Transformar resultados a archivos compatibles con R
	qiime tools export 
		--input-path ../5.clasificacion_taxonomica/taxonomy.qza 
		--output-path exported-taxonomy

	qiime tools export 
		--input-path ../3.denoising_dada2/rep-seqs_test.qza 
		--output-path exported-repseqs	
