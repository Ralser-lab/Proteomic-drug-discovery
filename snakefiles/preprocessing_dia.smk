rule run_preprocessingdia:
    input:
        script = "preprocessing-dia/src/preprocessingdevice.py",
        metadata = "preprocessing-dia/input/20240314_AF_50-0121_metadata.xlsx",
        prmatrix = "preprocessing-dia/input/SB_PROTAC_prmatrix_240314a.tsv"
    log: 
        "preprocessing-dia/logs/preprocessing-dia.log"
    output:
        "preprocessing-dia/logs/preprocessing-dia.done"
    shell:
        """
        python3 {input.script} > {log} 2>&1
        touch {output}
        """