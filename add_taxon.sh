#!/bin/bash

############################################################################################################################################
# Written by Kirsty Kwok
# Last updated on 26/03/2024
# Description: this script will process diamond m8 / blastn file and annotate each row with superkingdom, kingdom, family, genus, and species info
# Usage: bash add_taxon.sh -i <inputFile> -f <inputFormat> m8/blastn
# Optional options: -d <path to acc2taxid.tsv.gz> -k <key> -n <headerName> -c <column number for acc/staxid> -e <if header exists in input file>
# Output file: tab-delimited m8annotated txt
# Dependent packages include: taxonkit, csvtk and access to default acc2taxid.tsv.gz file
# Method for annotating m8 file was adapted from https://bioinf.shenwei.me/taxonkit/tutorial/.
# superkingdom phylum class order family genus species subspecies/rank {k};{p};{c};{o};{f};{g};{s};{t}
# subspecies {S} strain {T} subspecies/strain {t}
############################################################################################################################################

usage() {
    echo "Usage: "
    echo "$0 -i <inputFile> -f <inputFormat> <m8/blastn>"
    echo "Optional parameters:"
    echo "-e input file contains header, it is a boolean flag, default is 0 (no header)"
    echo "-d <pathToAcc2taxid.tsv.gz>";
    echo "   default acc2taxid for m8: /home3/2893911k/kirsty_db/accession2taxid/protacc2taxid.tsv.gz"
    echo "   default acc2taxid for blastn: /home3/2893911k/kirsty_db/accession2taxid/nucl_gb_acc2taxid.tsv.gz"
    echo "   custom db format: tab-delimited gz text file with two columns (1st: accession 2nd: taxid) with header. Accession number should be without version number."
    echo "-k <acc/staxid>, taxon will be added based on this provided key, required for blastn file. Accession number will be used for m8."
    echo "   if key is staxid, make sure blastn column is staxid not staxids (semicolon separated, multiple staxids), using staxids instead of staxid will cause error."
    echo "-n <headerName>, separate by comma"
    echo "   default headers for m8: 'qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore', ncol:12"
    echo "   default headers for blastn: 'qseqid,sseqid,sscinames,staxid,staxids,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore', ncol:17"
    echo "-c <columnNumberOfAcc/staxid>, 1-based index"
    echo "   default column number for m8: 2"
    echo "   default column number for blastn: 4"
    1>&2
    exit 1;
}

if [ $# -eq 0 ]; then usage; exit 1; fi

acc2taxid=""; key=""; header=""; col_key=0, headerExists="no"; headerExists=0;

while getopts ":i:f:d:k:n:c:e" opt
do
        case "${opt}" in
                i) inputFile=${OPTARG} ;;
                f) inputFormat=${OPTARG} ;;
                d) acc2taxid=${OPTARG} ;;
                k) key=${OPTARG} ;;
                n) header=${OPTARG} ;;
                c) col_key=${OPTARG} ;;
                e) headerExists=1 ;;
                ?) echo "Invalid option" >&2
                exit 1
                ;;
        esac
done

shift $((OPTIND-1))

#Check if all files and parameters are  correctly provided.
if [[ "$inputFormat" != "m8" && "$inputFormat" != "blastn" ]]; then echo "INVALID input format."; exit 1; fi
if [[ "$inputFormat" == "m8" ]]; then key="acc"; echo "***Accession number will be used for m8 file.***"; fi
if [[ "$key" != "staxid" && "$key" != "acc" ]]; then echo "INVALID key."; exit 1; fi

#Check input File and acc2taxid file and file id
if [[ ! -f "$inputFile" ]];  then  echo "INVALID input file."; exit 1; fi
if [[ "$inputFormat" == "m8" &&  "$acc2taxid" ==  "" ]]; then acc2taxid="/home3/2893911k/kirsty_db/accession2taxid/protacc2taxid.tsv.gz"; fi
if [[ "$inputFormat" == "blastn" && "$acc2taxid" ==  "" && "$key" == "acc" ]]; then acc2taxid="/home3/2893911k/kirsty_db/accession2taxid/nucl_gb_acc2taxid.tsv.gz"; fi
if [[ "$inputFormat" == "blastn" && "$key" == "staxid" ]]
then
    acc2taxid="taxdump.tar.gz in .taxonkit directory will be used."
else
    if [[ ! -f "$acc2taxid" ]];  then  echo "INVALID acc2taxid path."; exit 1; fi
fi
FileID=$(basename $inputFile | cut -f 1 -d '.')
if [[ $FileID == "" ]]; then echo "INVALID file id."; exit 1; fi

#Set default header for m8 and blastn and check if number of header matches with number of column in the input files.
if [[ "$header" == "" && "$inputFormat" == "m8" ]]; then header="qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore"; col_key=2; fi
if [[ "$header" == "" && "$inputFormat" == "blastn" && "$key" == "staxid" ]]; then header="qseqid,sseqid,sscinames,staxid,staxids,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore"; col_key=4; fi
if [[ "$header" == "" && "$inputFormat" == "blastn" && "$key" == "acc" ]]; then header="qseqid,sseqid,sscinames,staxid,staxids,pident,length,qlen,slen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore"; col_key=2; fi
if [[ "$col_key" == 0 ]]; then echo "Please provide a valid column number for key."; exit 1; fi
comma_count=$(awk -F',' '{print NF-1}' <<< "${header}")
nheader=$((comma_count+1))
ncol=$(csvtk ncol -Ht $inputFile)
if [[ $nheader != $ncol ]]; then echo "Number of header does not match with number of columns in the input file, please provide headers and try again."; exit 1; fi

#Print out the parameters
echo "input file: ""$inputFile"
echo "input format: ""$inputFormat"
echo "file id: ""$FileID"
echo "key: ""$key"
file_prefix=${FileID}.${inputFormat}.${key}
echo "file_prefix: ""$file_prefix"
echo "column number for key: ""$col_key"
echo "header: ""$header"
echo "header exists: ""$headerExists"
echo "acc2taxid file: ""$acc2taxid"

if [[ "$headerExists" == 1 ]]; then sed '1d' $inputFile > ${file_prefix}.tmp.input; else cat $inputFile > ${file_prefix}.tmp.input; fi

# This module is used if provided file is m8 format
# Extract accession and deduplicate and remove versions
if [[ "$inputFormat" ==  "m8" ]]
then
    #Extract accession from m8 file and remove duplicates
    cut -f $col_key ${file_prefix}.tmp.input | cut -f 1 -d '.' | csvtk uniq -Ht > ${file_prefix}.acc
    #Find taxid for each accession
    zcat $acc2taxid | grep -w -f ${file_prefix}.acc > ${file_prefix}.hits
    #Extract taxid, remove duplicates, and prepare taxid2name.tsv, superkingdom, kingdom, family, genus, and species name are retrived for the taxids.
    cut -f 2 ${file_prefix}.hits | csvtk uniq -Ht | taxonkit reformat -f '{k}@{K}@{f}@{g}@{s}' -I 1 > ${file_prefix}.taxid2name
    sed -i 's/\"/*/g' ${file_prefix}.taxid2name
    sed -i "s/\'/*/g" ${file_prefix}.taxid2name
    #Join taxid2name.tsv with hits file
    csvtk join --left-join --na NA -Ht -f "2;1" ${file_prefix}.hits ${file_prefix}.taxid2name -To ${file_prefix}.acc_taxid2name

    #Add accesion number to the m8 file
    rm -f ${file_prefix}.tmpblastn
    while IFS= read -r line
    do
        accNum=$(echo "$line" | cut -f $col_key -d$'\t' |  cut -f 1 -d '.')
        echo -e "$line\t$accNum" >> ${file_prefix}.tmpblastn
    done < ${file_prefix}.tmp.input
    ncol=$(csvtk ncol -Ht ${file_prefix}.tmpblastn)

    #Annotate m8 file with taxon info
    csvtk join --left-join --na NA -Ht -f "${ncol};1" ${file_prefix}.tmpblastn ${file_prefix}.acc_taxid2name  \
      | csvtk add-header -t --names $header,accNum,taxid,taxon \
      | csvtk sep -t -f taxon -n superkingdom,kingdom,family,genus,species -s '@' --remove --drop > ${file_prefix}.tmp

# This module is used if provided file is blastn format and staxid is used as taxon retreival key
elif [[ "$inputFormat" ==  "blastn" && "$key" == "staxid" ]]
then
    #Extract taxid from blastn file and remove duplicates, prepare taxid2name.tsv, superkingdom, kingdom, family, genus, and species name are retrived using uniq taxids originally listed on blastn file.
    cut -f $col_key ${file_prefix}.tmp.input | csvtk uniq -Ht  | taxonkit reformat -f '{k}@{K}@{f}@{g}@{s}' -I 1 > ${file_prefix}.taxid2name

    #Annotate blastn file with taxon info
    csvtk join --left-join --na NA -Ht -f "$col_key;1"  ${file_prefix}.tmp.input ${file_prefix}.taxid2name\
    | csvtk add-header -t --names $header,taxon \
    | csvtk sep -t -f taxon -n superkingdom,kingdom,family,genus,species -s '@' --remove --drop > ${file_prefix}.tmp

# This module is used if provided file is blastn format and accession number is used as taxon retreival key
elif [[ "$inputFormat" ==  "blastn" && "$key" == "acc" ]]
then
    # Extract accession from blastn file and remove duplicates
    cut -f $col_key ${file_prefix}.tmp.input | cut -f 4 -d '|' |  cut -f 1 -d '.' | csvtk uniq -Ht > ${file_prefix}.acc
    zcat $acc2taxid | grep -w -f ${file_prefix}.acc > ${file_prefix}.hits
    cut -f 2 ${file_prefix}.hits | csvtk uniq -Ht | taxonkit reformat -f '{k}@{K}@{f}@{g}@{s}' -I 1 > ${file_prefix}.taxid2name
    csvtk join --left-join --na NA -Ht -f "2;1" ${file_prefix}.hits ${file_prefix}.taxid2name -To ${file_prefix}.acc_taxid2name

    # Extract accession number from blastn file and append to the file
    rm -f ${file_prefix}.tmpblastn
    while IFS= read -r line
    do
        accNum=$(echo "$line" | cut -f $col_key -d$'\t' | cut -f 4 -d '|' |  cut -f 1 -d '.')
        echo -e "$line\t$accNum" >> ${file_prefix}.tmpblastn
    done < ${file_prefix}.tmp.input
    ncol=$(csvtk ncol -Ht ${file_prefix}.tmpblastn)

    csvtk join --left-join --na NA -Ht -f "${ncol};1" ${file_prefix}.tmpblastn ${file_prefix}.acc_taxid2name \
    | csvtk add-header -t --names $header,accNum,taxid,taxon \
    | csvtk sep -t -f taxon -n superkingdom,kingdom,family,genus,species -s '@' --remove --drop > ${file_prefix}.tmp
fi

# This step adds a FileID column to the final annotated m8 file for subsequent analysis
# Need double quote "' '"for using bash variable
ID=${FileID} && cat ${file_prefix}.tmp | csvtk mutate2 -t -e " '$ID' " --at 1 -n file_id  > ${FileID}.${inputFormat}annotated

if [ -s "${FileID}.${inputFormat}annotated" ]
then
    original_nrow=$(csvtk nrow -Ht ${file_prefix}.tmp.input)
    annotated_nrow=$(csvtk nrow -t ${FileID}.${inputFormat}annotated)
    if [[ $original_nrow == $annotated_nrow ]]
    then
        echo "All data rows should be annotated, removing intermediate files..."
        rm -f ${file_prefix}.acc ${file_prefix}.hits ${file_prefix}*taxid2name ${file_prefix}.tmp*
    else
        if [[ $original_nrow > $annotated_nrow ]]; then echo "Some data rows are missing after annotation, please check the output and intermediate files."; fi
        if [[ $original_nrow < $annotated_nrow ]]; then echo "Some data rows are potentially duplicated after annotation, please check the output and intermediate files."; fi
    fi
    echo "Output file: ${FileID}.${inputFormat}annotated"
fi