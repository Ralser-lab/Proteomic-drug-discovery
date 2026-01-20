rule run_deepsearch:
    input:
        script = "code/protacs_23_gbdt_train_wide.py",
        logs = "logs/manuscript_steps.done"
    params:
        hyper = "configs/cfg_optuna_wide.json"
    shell:
        """
        python3 {input.script} {params.hyper}
        """