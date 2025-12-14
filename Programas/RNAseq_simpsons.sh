## Este script sacará análisis de la calidad de secuencias FASTQ, realizará un paso de 'Trimming' 
## (filtrado de posiviones y secuencias de baja calidad), y cuantificará niveles de expresión de 
## los genes dados en el archivo 'Referencia.fasta'.

## El script está diseñado para cortar y pegar los comandos y ejecutarlos en un terminal de Linux
## Lo único que hay que tener cuidado es con el directorio de trabajo (definir ruta en la variable
## projectDir), donde tienen que estar los directorios 'Fastqs' y 'Genes', así como los archivos
## 'Referencia.fasta'


## Crear sesión de conda  con herramientas para QC, alineamiento y conteo de lecturas mapeadas

conda create -n simpsons -c bioconda -c conda-forge -c defaults fastqc multiqc fastp salmon

conda activate simpsons

## Generar path del proyecto

#projectDir="~/TallerGrupal_Ficheros"
projectDir="/mnt/c/Users/Ruben/Documents/Ruben/Cursos/UNIR_MasterBioinfo/Secuenciacion_y_omicas/actividad2/mubio03_act2/TallerGrupal_Ficheros"

## Generar árbol de directorios en nuestro directorio

mkdir $projectDir/Resultados/

cd $projectDir/Resultados

mkdir FastQC_raw Fastp FastQC_trimm Salmon

## Generar una lista con las lecturas que nos interesen (Obeso1 y normopeso)

cd $projectDir/Fastqs

ls *Simpson* | cut -d$'_' -f1 | sort | uniq | grep -v 'MargeSimpson' > lista.txt

mv lista.txt $projectDir

## Comprobar la calidad de las lecturas crudas con FastQC y MultiQC

cd $projectDir

cat lista.txt | while read line; do
  mkdir $projectDir/Resultados/FastQC_raw/$line\_R1.fastq
  fastqc -o $projectDir/Resultados/FastQC_raw/$line\_R1.fastq $projectDir/Fastqs/$line\_R1.fastq.gz 
  mkdir $projectDir/Resultados/FastQC_raw/$line\_R2.fastq
  fastqc -o $projectDir/Resultados/FastQC_raw/$line\_R2.fastq $projectDir/Fastqs/$line\_R2.fastq.gz
done

cd $projectDir/Resultados/FastQC_raw

multiqc .

## Eliminar posiciones de baja calidad y eliminar adaptadores con fastp

cd $projectDir

cat lista.txt | while read line; do
  fastp --in1 $projectDir/Fastqs/$line\_R1.fastq.gz \
  --in2 $projectDir/Fastqs/$line\_R2.fastq.gz \
  --out1 $projectDir/Resultados/Fastp/$line\_R1.trimm.fastq.gz \
  --out2 $projectDir/Resultados/Fastp/$line\_R2.trimm.fastq.gz \
  --detect_adapter_for_pe \
  --cut_front \
  --cut_tail \
  --cut_window_size 12 \
  --cut_mean_quality 20 \
  --length_required 35 \
  --json $projectDir/Resultados/Fastp/$line.json \
  --html $projectDir/Resultados/Fastp/$line.html
done

## Comprobar de nuevo la calidad de las lecturas filtradas con FastQC y MultiQC

cd $projectDir

cat lista.txt | while read line; do
  mkdir $projectDir/Resultados/FastQC_trimm/$line\_R1.trimm.fastq
  fastqc -o $projectDir/Resultados/FastQC_trimm/$line\_R1.trimm.fastq $projectDir/Resultados/Fastp/$line\_R1.trimm.fastq.gz 
  mkdir $projectDir/Resultados/FastQC_trimm/$line\_R2.trimm.fastq
  fastqc -o $projectDir/Resultados/FastQC_trimm/$line\_R2.trimm.fastq $projectDir/Resultados/Fastp/$line\_R2.trimm.fastq.gz 
done

cd $projectDir/Resultados/FastQC_trimm

multiqc .

## Indexar transcriptoma de referencia con salmon

cd $projectDir

salmon index -t $projectDir/Referencia.fasta -i Referencias_index -k 31

## Sumarizar los conteos con salmon

cat lista.txt | while read line; do
  salmon quant -i referencias_index \
  -l IU \
  -1 $projectDir/Resultados/Fastp/$line\_R1.trimm.fastq.gz \
  -2 $projectDir/Resultados/Fastp/$line\_R2.trimm.fastq.gz \
  --validateMappings \
  -o $projectDir/Resultados/Salmon/$line\_quant
done

## Generar matriz de expresión

cat lista.txt | while read line; do
  cut -f1 $projectDir/Resultados/Salmon/$line\_quant/quant.sf > $projectDir/Resultados/Salmon/$line\_referencia.sf
  cut -f5 $projectDir/Resultados/Salmon/$line\_quant/quant.sf > $projectDir/Resultados/Salmon/$line\_quant.sf
  paste $projectDir/Resultados/Salmon/$line\_referencia.sf $projectDir/Resultados/Salmon/$line\_quant.sf > $projectDir/Resultados/Salmon/$line\_quant.csv
  rm $projectDir/Resultados/Salmon/$line\_referencia.sf $projectDir/Resultados/Salmon/$line\_quant.sf
done

touch rnaseqHeader 
echo "Gene AbrahamSimpson HomerSimpson BartSimpson LisaSimpson MaggieSimpson" | tr ' ' '\t' > rnaseqHeader
paste *.csv | cut -f1,2,4,6,8,10 | tail -n +2 | cat rnaseqHeader - > SimpsonsExpression.csv
rm rnaseqHeader *_quant.csv

## Usaremos el archivo SimpsonsExpression.csv para realizar el análisis de expresión diferencial

mkdir $projectDir/Resultados/R

ln -s $projectDir/Resultados/Salmon/SimpsonsExpression.csv $projectDir/Resultados/R/SimpsonsExpression.csv
