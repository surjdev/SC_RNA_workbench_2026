import yaml

from workbench_utils import load_config


def test_load_default_config():
    cfg = load_config()
    assert "project" in cfg
    assert "qc" in cfg
    assert cfg["project"]["random_seed"] == 42
    assert cfg["qc"]["min_genes"] == 200


def test_load_custom_config(tmp_path):
    custom_yaml = tmp_path / "custom.yaml"
    custom_data = {
        "project": {"name": "Custom_Experiment", "random_seed": 99},
        "qc": {"min_genes": 350},
    }
    with open(custom_yaml, "w") as f:
        yaml.dump(custom_data, f)

    cfg = load_config(custom_yaml)
    assert cfg["project"]["name"] == "Custom_Experiment"
    assert cfg["project"]["random_seed"] == 99
    assert cfg["qc"]["min_genes"] == 350
    # Retains defaults for unspecified values
    assert cfg["qc"]["min_counts"] == 500
