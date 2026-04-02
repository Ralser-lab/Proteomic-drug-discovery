rule run_preprocessingdia:
    input:
        script = "preprocessing_dia/src/preprocessingdevice.py",
        metadata = "preprocessing_dia/input/20240314_AF_50-0121_metadata.xlsx",
        prmatrix = "preprocessing_dia/input/SB_PROTAC_prmatrix_240314a.tsv"
    log: 
        "preprocessing_dia/logs/preprocessing_dia.log"
    output:
        "preprocessing_dia/logs/preprocessing_dia.done"
    shell:
        """
        python3 {input.script} > {log} 2>&1
        touch {output}
        """