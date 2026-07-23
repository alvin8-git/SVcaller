// MELT_MERGE collects the 4 per-type VCFs from the MELT_CALL fan-out and produces
// the single ${meta.id}.melt.vcf.gz the ensemble consumes. This is the exact merge
// that used to live at the tail of the old serial MELT_CALL: concatenate, keep
// canonical chr + PASS, strip INFO to SVTYPE/MEITYPE/SVLEN/END, inject the MEITYPE
// header, sort, gzip. Content is identical to the old single-task output.
process MELT_MERGE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("${meta.id}.melt.vcf.gz"), emit: vcf

    script:
    """
    # The 4 fan-out tasks each emit ${meta.id}.<TYPE>.final_comp.vcf. Sorting the names
    # gives ALU, HERVK, LINE1, SVA -- the same order the old glob produced, so the header
    # is taken from ALU exactly as before. The final coordinate sort makes record order
    # independent of concat order regardless.
    first=\$(ls *.final_comp.vcf 2>/dev/null | sort | head -1 || true)
    if [ -z "\$first" ]; then
        echo "ERROR: MELT_MERGE received no *.final_comp.vcf from the fan-out tasks." >&2
        echo "Refusing to emit an empty VCF that would be read as 'no MEIs found'." >&2
        ls -l >&2 || true
        exit 1
    fi

    (
        grep "^#" "\$first"
        for vcf in \$(ls *.final_comp.vcf | sort); do
            grep -v "^#" "\$vcf"
        done
    ) | awk '
        BEGIN{OFS="\\t"}
        /^##fileformat/{
            print
            print "##INFO=<ID=MEITYPE,Number=1,Type=String,Description=\\"Mobile element insertion type (ALU/LINE1/SVA/HERVK)\\">"
            next
        }
        /^#/{print;next}
        \$1!~/^chr([0-9]+|X|Y|M)\$/{next}
        \$7!="PASS" && \$7!="."{next}
        {
            # Strip INFO to essentials; normalise SVTYPE → INS; store ME type as MEITYPE
            svtype="INS"; meitype=""; svlen=""; end_=""
            n=split(\$8,info,";")
            for(i=1;i<=n;i++){
                f=info[i]
                if(f~/^SVTYPE=/){
                    t=substr(f,8)
                    if(t~/^(ALU|LINE1|LINE|SVA|HERVK)\$/) { meitype=t; svtype="INS" }
                    else svtype=t
                }
                else if(f~/^SVLEN=/) svlen=f
                else if(f~/^END=/)   end_=f
            }
            new_info="SVTYPE=" svtype
            if(meitype!="") new_info=new_info ";MEITYPE=" meitype
            if(svlen!="")   new_info=new_info ";" svlen
            if(end_!="")    new_info=new_info ";" end_
            \$8=new_info
            # Add FORMAT/GT if absent
            if(\$9=="" || \$9=="."){
                \$9="GT"; \$10="0/1"
            }
            print
        }' | sort -k1,1 -k2,2n | gzip > ${meta.id}.melt.vcf.gz
    """
}
