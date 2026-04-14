rule train_gbdt_wide:
    input:
        script = "code/protacs_16_gbdt_train_retrain_wide.py",
        logs = "logs/manuscript_steps.done"
    shell:
        """
        python3 {input.script}
        """
