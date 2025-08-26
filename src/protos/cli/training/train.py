import argparse
import os
from pathlib import Path
from models.config.config_manager import get_config, Config
from models.model_trainer import ModelTrainer


def train_model(run_id):
    """
    Train a model with the given run ID.
    
    Args:
        run_id (str): Run identifier for the configuration and output logging
    
    Returns:
        bool: True if training completed successfully
    """
    # Construct the configuration file path using the run_id
    config_path = Path(f"configs/config_{run_id.zfill(6)}.ini")
    
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    config_dict = get_config(config_path)
    config = Config(config_dict)

    # Initialize the model trainer with the loaded configuration
    trainer = ModelTrainer(config)
    trainer.train()
    
    return True


def main():
    """Command-line entry point for model training."""
    parser = argparse.ArgumentParser(description="Train a model with given configuration")
    parser.add_argument('--run_id', type=str, required=True, 
                      help='Run identifier for the configuration and output logging')
    args = parser.parse_args()

    try:
        train_model(args.run_id)
        print(f"Training completed successfully for run ID: {args.run_id}")
    except Exception as e:
        print(f"Error during training: {str(e)}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())