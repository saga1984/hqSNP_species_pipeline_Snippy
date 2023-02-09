#!/bin/bash

############################################################################################################
# correr snippy en modo de una sola secuencia sobre todos los ensambles en un directorio (ASSEMBLY)        #
# este script selecciona ensambles en base a su nombre buscando la especie como palabra clave,             #
# esta palabra clave se guarda en la variable "ref_name", la cual puede ser necesario modificar su codigo  #
# para identificar el nombre de la especie                                                                 #
############################################################################################################

run_snippy() {
   #################
   # correr Snippy #
   #################

   # asignar nombre clave de genoma de referencia
   ref_name=$(basename ${ref} .fa)

   # for loop para todos los ensambles dentro de la carpeta de ensambles ASSEMBLY
   for ensamble in ASSEMBLY/*${especie}*.fa; do
      # nombre corto de ensambles
      ensamble_name="$(basename ${ensamble} .fa)"
      # fila en blanco (espacio)
      echo -e "\n"
      # imprime nombre de ensamble
      echo "##########################################################"
      echo "Ensamble: ${ensamble_name}; Referencia: ${ref_name}"
      echo "##########################################################"
      echo ""

      # asignar directorio de salida para guardar resultados de Snippy
      dir="SNIPPY_${ref_name}"
      # ejecuta Snippy
      snippy --cpus $(nproc) --force --ref ${ref} --outdir ${dir}/coreSNP_${ensamble_name} --ctgs ${ensamble} \
      --ram $(grep "MemTotal" /proc/meminfo | awk '{print $2/(1024 * 1024)}' | cut -d "." -f "1") # info de la ram, se divide para Gigas y se eliminan decimales
   done

#####################################################################################
# correr snippy-core, snippy-clean. gubbins, snp-sites, fasttree, raxml y snp-dists #
#####################################################################################

   # moverse al directorio correspondiente
   cd ${dir}

   echo "#########################"
   echo " ejecutando snippy core "
   echo -e "#########################\n"

   # ejecuta snippy-core
   snippy-core --ref ../${ref} --prefix core $(echo coreSNP*)

   echo "############################################################"
   echo "  limpiando alineamiento para generar SNPs de alta calidad  "
   echo -e "############################################################\n"

   # ejecuta snippy-clean, para limpiar alineamiento
   snippy-clean_full_aln core.full.aln > clean.full.aln
   # correr gubbins con arbol rapido, para limpiar alineamiento
   run_gubbins --verbose --threads $(nproc) --tree_builder fasttree --prefix gubbins clean.full.aln
   # si gubbins falla con error (porque son muy pocas secuencias)
   if [[ $? != 0 ]]; then
      # guardar el resultado de porcentajes de missing data
      run_gubbins --verbose --threads $(nproc) --tree_builder fasttree --prefix gubbins clean.full.aln > tmp_porcentajes.txt
      # filtrar para obtener el valor maximo de porcentajes de missing data encontrados
      maximo=$(cat tmp_porcentajes.txt | grep 'percentage' | awk '{print $NF}' | awk '/[0-9]/' | sort -r | sed -n '1p')
      # volver a correr gubbins con nuevo criterio porcentaje de missing data
      run_gubbins --verbose --threads $(nproc) --filter_percentage $(echo ${maximo} + 1 | bc | cut -d '.' -f '1') --tree_builder fasttree --prefix gubbins clean.full.aln
      # eliminar archivos temporales
      rm tmp*
   fi
   # correr snp-sites, para limpiar alineamiento
   snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln

   echo "#############################################"
   echo "  obteniendo reconstrucciones filogeneticas  "
   echo -e "#############################################\n"

   # correr fasttree, para hacer reconstruccion filogenetica por maximum likelihood
   FastTree -gtr -nt clean.core.aln > FastTree_clean.core.tree
   # correr raxml, para hacer segunda reconstruccion filogenetica por maximum likelihood (con 100 bootstraps)
   raxmlHPC -f a -p 1234567890 -s clean.core.aln -x 1234567890 -# 100 -m GTRGAMMA -n clean.core.newick
   # correr snp-dists, para obtener matriz de distancias de SNPs
   snp-dists -j $(nproc) clean.core.aln > Genero_SNP_matrix.tsv
}

# ---------------------------------------
#  uso correcto del script: script_mapeo
# ---------------------------------------

# variable que indica la forma correcta de usar el script
USO() {
echo -e "  uso:\t$(basename $0) -a <reference assembly> -e <species>\n"
echo -e "  -a)\t assembly used to map for SNPs (mandatory argument)"
echo -e "  -e)\t species of interest to filter assemblies (mandatory argument)\n"
echo -e "  -h)\t shows this help menu"
}

# si el script no se ejecuta correctamente (sin argumentos/opciones), entonces indicar uso correcto y salir con error
if [[ $# -eq 0 ]]; then
   USO
   exit 2
fi

# ----------------------------------------------------------------------
#  parsear de argumentos y ejecutar el script script_mapeo_ref_virus.sh
# ----------------------------------------------------------------------

# parseo de opciones/argumentos
while getopts ":a:e:h"  opciones; do
   case "${opciones}" in
      a)
         # si no existe la opcion -e manda mensaje de uso y sal con error
         if [[ $# -lt 3 ]]; then
            echo -e "\n\tSpecies option/argument needed \n"
            USO
         else
            ref="${OPTARG}" # introducir ensamble de referencia
         fi
         ;;
      e)
         especie="${OPTARG}" # introducir especie de interes
         # si no existe la opcion -a manda mensaje de uso y sal con error
         if [[ ! -z ${ref} ]]; then
            # mensaje de inicio de analisis
            echo "##################################################################################################"
            echo -e "Comenzando llamado de variantes, usando como genoma de referencia: ${ref} y especie: ${especie}"
            echo "##################################################################################################"
            echo ""

            # correr Snippy
            run_snippy

            # mensaje de finalizacion de analisis
            echo "#######################################"
            echo -e "  Llamado de variantes terminado  "
            echo "#######################################"
            echo ""
         else
            echo -e "\n\tReference assembly needed"
            USO
         fi
         ;;
      h)
         USO
         ;;
      ?)
         echo "invalid option -${OPTARG}"
         USO
         ;;
   esac
done
