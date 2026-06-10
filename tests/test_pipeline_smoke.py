from __future__ import annotations

import subprocess
import sys


def test_experiment_help_commands_do_not_require_raw_data():
    modules = [
        "src.experiments.exp1.run_classification",
        "src.experiments.exp2.run_multiscale",
        "src.experiments.exp3.run_variant_effect",
    ]
    for module in modules:
        result = subprocess.run(
            [sys.executable, "-m", module, "--help"],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0, result.stderr
        assert "--config" in result.stdout
