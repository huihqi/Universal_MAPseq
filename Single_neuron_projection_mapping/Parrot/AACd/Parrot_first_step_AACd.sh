#!/bin/bash

set -euo pipefail
shopt -s nullglob

Prefix="Parrot_AACd"
CSI_file="CSI_parrot.txt"

samples=(
    "P4LAACd"
    "P4LAACdv"
    "P4LNon_AACd"
    "P4LAACv"
    "P4LNLC"
    "P4LNAO"
    "P4LMO"
    "P4LLPO"
    "P4LDMm"
    "P4LDM"
    "P4LnXIIts"
    "P4LPAM-RAM"
    "P4RAACd"
    "P4RAACdv"
    "P4RNon_AACd"
    "P4RAACv"
    "P4RNLC"
    "P4RNAO"
    "P4RMO"
    "P4RLPO"
    "P4RDMm"
    "P4RDM"
    "P4RnXIIts"
    "P4RPAM-RAM"
    "P4RCerebellum"
    "P5LAACd"
    "P5LAACdv"
    "P5LNon_AACd"
    "P5LAACv"
    "P5LNLC"
    "P5LNAO"
    "P5LMO"
    "P5LLPO"
    "P5LDMm"
    "P5LDM"
    "P5LnXIIts"
    "P5LPAM-RAM"
    "P5RAACd"
    "P5RAACdv"
    "P5RNon_AACd"
    "P5RAACv"
    "P5RNLC"
    "P5RNAO"
    "P5RMO"
    "P5RLPO"
    "P5RDMm"
    "P5RDM"
    "P5RnXIIts"
    "P5RPAM-RAM"
    "P5RCerebellum"
)

for index in "${!samples[@]}"; do

    # Bash array begins at 0, so add 1 to generate BC1–BC50
    i=$((index + 1))
    sample="${samples[$index]}"

    echo "=========================================="
    echo "Processing BC${i}: ${sample}"

    # Find paired-end FASTQ files.
    # This pattern works regardless of the exact S number.
    R1_files=( "${sample}"_S*_L003_R1_001.fastq.gz )
    R2_files=( "${sample}"_S*_L003_R2_001.fastq.gz )

    if (( ${#R1_files[@]} == 0 )); then
        echo "ERROR: R1 file not found for ${sample}" >&2
        continue
    fi

    if (( ${#R2_files[@]} == 0 )); then
        echo "ERROR: R2 file not found for ${sample}" >&2
        continue
    fi

    if (( ${#R1_files[@]} > 1 )); then
        echo "ERROR: More than one R1 file found for ${sample}:" >&2
        printf '  %s\n' "${R1_files[@]}" >&2
        continue
    fi

    if (( ${#R2_files[@]} > 1 )); then
        echo "ERROR: More than one R2 file found for ${sample}:" >&2
        printf '  %s\n' "${R2_files[@]}" >&2
        continue
    fi

    R1="${R1_files[0]}"
    R2="${R2_files[0]}"

    echo "R1: ${R1}"
    echo "R2: ${R2}"

    R1_stripped="${Prefix}_BC${i}_R1_stripped.txt"
    PE_file="${Prefix}_BC${i}_PE.txt"
    BU_file="${Prefix}_BC${i}_BU.txt"
    counts_file="${Prefix}_BC${i}_BU_counts.txt"
    seq_file="${Prefix}_BC${i}_BU_seq.txt"

    # Extract the first 30 bases from R1 sequence lines
    zcat "${R1}" |
        awk 'NR % 4 == 2' |
        cut -b 1-30 \
        > "${R1_stripped}"

    # Extract the first 20 bases from R2 and append them to R1
    zcat "${R2}" |
        awk 'NR % 4 == 2' |
        cut -b 1-20 |
        paste -d '' "${R1_stripped}" - \
        > "${PE_file}"

    # Read the first 8 bases of the corresponding CSI sequence
    C=$(awk -v row="${i}" 'NR == row {print substr($0,1,8); exit}' "${CSI_file}")

    if [[ -z "${C}" ]]; then
        echo "ERROR: CSI sequence was not found on row ${i}" >&2
        continue
    fi

    echo "CSI for BC${i}: ${C}"

    # Keep sequences ending with the corresponding CSI,
    # remove sequences containing N,
    # retain the first 42 bases,
    # and count identical sequences
    awk -v csi="${C}" '
        index($0, "N") == 0 && substr($0, length($0)-length(csi)+1) == csi {
            print substr($0, 1, 42)
        }
    ' "${PE_file}" |
        sort |
        uniq -c |
        sort -k1,1nr \
        > "${BU_file}"

    # Separate read counts and barcode sequences
    awk '{print $1}' "${BU_file}" > "${counts_file}"
    awk '{print $2}' "${BU_file}" > "${seq_file}"

    echo "Finished BC${i}: ${sample}"

done

echo "=========================================="
echo "Finished processing BC1–BC50."