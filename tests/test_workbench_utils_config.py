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


def test_config_errors_and_save(tmp_path):
    import pytest

    from workbench_utils.config import print_config_summary, save_config, validate_config

    # Missing file
    with pytest.raises(FileNotFoundError):
        load_config(tmp_path / "non_existent.yaml")

    # Missing section
    invalid_cfg = {"project": {"random_seed": 42}}
    with pytest.raises(KeyError, match="Missing required configuration section"):
        validate_config(invalid_cfg)

    # Missing random seed
    invalid_seed_cfg = {
        "project": {},
        "qc": {},
        "normalization": {},
        "reduction": {},
        "clustering": {},
    }
    with pytest.raises(KeyError, match="Missing 'random_seed'"):
        validate_config(invalid_seed_cfg)

    # Save config and load back
    cfg = load_config()
    out_yaml = tmp_path / "saved.yaml"
    save_config(cfg, out_yaml)
    assert out_yaml.exists()

    # Print summary
    print_config_summary(cfg)
