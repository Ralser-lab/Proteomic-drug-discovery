rule run_deepsearch:
    input:
        script = "code/protacs_22_gbdt_train_narrow.py",
        logs = "logs/manuscript_steps.done"
    params:
        hyper = "configs/cfg_optuna_narrow.json"
    shell:
        """
        python3 {input.script} {params.hyper}
        """