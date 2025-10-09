rule run_deepsearch:
    input:
        script = "code/protacs_22_gbdt_deepsearch.py",
        logs = "logs/all_manuscript_steps.done"
    params:
        hyper = "configs/HYPER.json"
    shell:
        """
        python3 {input.script} {params.hyper}
        """