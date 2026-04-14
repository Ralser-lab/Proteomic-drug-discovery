rule train_gbdt:
    input:
        script = "code/protacs_17_gbdt_train_retrain_narrow.py",
        logs = "logs/manuscript_steps.done"
    shell:
        """
        python3 {input.script} 
        """