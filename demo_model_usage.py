"""Demo of how the model system would work in practice."""

from pathlib import Path
import sys

from protos.models import (
    StandardModel, ModelRegistry, ModelDownloader,
    list_available_models, get_model_definition,
    SequenceFormat, validate_model_compatibility
)


def demo_model_discovery():
    """Show how to discover available models."""
    print("=== Model Discovery ===\n")
    
    # List all available models
    models = list_available_models()
    print(f"Available models: {models}\n")
    
    # Get details about ESM-2
    esm2_def = get_model_definition("esm2")
    print("ESM-2 Details:")
    print(f"  Full name: {esm2_def.full_name}")
    print(f"  Description: {esm2_def.description}")
    print(f"  Input: {[f.value for f in esm2_def.input_formats]}")
    print(f"  Output: {esm2_def.output_format.value}")
    print(f"  Variants: {list(esm2_def.sources.keys())}")
    print(f"  Paper: {esm2_def.paper_url}")
    
    # Get details about Lambda
    lambda_def = get_model_definition("lambda")
    print("\nLambda Details:")
    print(f"  Full name: {lambda_def.full_name}")
    print(f"  Description: {lambda_def.description}")
    print(f"  Input: {[f.value for f in lambda_def.input_formats]}")
    print(f"  Output: {lambda_def.output_format.value}")
    print(f"  Architecture: {lambda_def.model_config.get('architecture')}")


def demo_model_usage():
    """Show how to use a model (without actually loading weights)."""
    print("\n\n=== Model Usage Demo ===\n")
    
    # Create a model instance
    print("Creating ESM-2 model instance...")
    model = StandardModel("esm2", "esm2_t6_8M", device="cpu")
    
    # Get model info
    info = model.get_model_info()
    print(f"\nModel info:")
    print(f"  Name: {info['name']}")
    print(f"  Version: {info['version']}")
    print(f"  Device: {info['device']}")
    print(f"  Input formats: {info['input_formats']}")
    print(f"  Output format: {info['output_format']}")
    
    print("\nIn practice, you would:")
    print("1. Load the model: model.load_model()")
    print("2. Make predictions: prediction = model.predict('protein_name')")
    print("3. Access results: embeddings = prediction.prediction['embeddings']")


def demo_model_workflow():
    """Show a complete workflow."""
    print("\n\n=== Complete Workflow Demo ===\n")
    
    # Step 1: Check what formats we have for an entity
    available_formats = ["structure", "sequence", "grn"]
    print(f"Entity 'BACR_HALSA' has formats: {available_formats}")
    
    # Step 2: Find compatible models
    print("\nChecking model compatibility:")
    
    models_to_check = {
        "esm2": ["sequence"],
        "lambda": ["structure", "sequence", "grn"],
        "esmfold": ["sequence"],
        "boltz1": ["sequence", "msa"]
    }
    
    compatible_models = []
    for model_name, required in models_to_check.items():
        if validate_model_compatibility(available_formats, required):
            compatible_models.append(model_name)
            print(f"  ✓ {model_name}: Compatible")
        else:
            missing = set(required) - set(available_formats)
            print(f"  ✗ {model_name}: Missing {missing}")
    
    # Step 3: Use a compatible model
    print(f"\nUsing compatible model: {compatible_models[0]}")
    
    # Create sequence format
    seq = SequenceFormat(
        sequence="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK",
        sequence_id="BACR_HALSA"
    )
    
    print(f"Sequence length: {len(seq.sequence)}")
    
    # Would create and use model
    print("\nWorkflow steps:")
    print("1. model = StandardModel('esm2', 'esm2_t33_650M')")
    print("2. model.load_model()")
    print("3. prediction = model.predict('BACR_HALSA')")
    print("4. embeddings = prediction.prediction['embeddings']")
    print("5. Save/analyze embeddings...")


def demo_model_registry():
    """Show model registry usage."""
    print("\n\n=== Model Registry Demo ===\n")
    
    registry = ModelRegistry()
    
    # List models (from registry file)
    print("Registry operations:")
    print("1. Discover models: registry.discover_models()")
    print("2. List models: registry.list_models()")
    print("3. Get model info: registry.get_model_info('standard/esm2_650M')")
    print("4. Load model: model = registry.load_model('standard/esm2_650M')")
    print("5. Find compatible: registry.get_models_for_entity('protein_name')")


def demo_downloader():
    """Show model downloader usage."""
    print("\n\n=== Model Downloader Demo ===\n")
    
    downloader = ModelDownloader()
    
    # Check download sizes
    print("Download sizes:")
    for model in ["esm2", "ankh", "lambda"]:
        try:
            definition = get_model_definition(model)
            for variant, source in list(definition.sources.items())[:2]:
                size = source.size_mb or 0
                print(f"  {model}/{variant}: {size:.1f} MB")
        except:
            pass
    
    print("\nDownloader operations:")
    print("1. Check space: downloader.check_disk_space(2500)")
    print("2. Download: path = downloader.download_model('esm2', 'esm2_t33_650M', target_dir)")
    print("3. List downloaded: downloader.list_downloaded_models()")


def main():
    """Run all demos."""
    print("=" * 70)
    print("PROTOS MODEL SYSTEM DEMO")
    print("=" * 70)
    print("\nThis demo shows how the model system works without")
    print("actually loading models or downloading weights.\n")
    
    demo_model_discovery()
    demo_model_usage()
    demo_model_workflow()
    demo_model_registry()
    demo_downloader()
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("\nThe Protos model system provides:")
    print("✓ Unified interface for all AI models")
    print("✓ Automatic format validation and conversion")
    print("✓ Model discovery and compatibility checking")
    print("✓ Automated downloading and installation")
    print("✓ Integration with Protos entity system")
    print("\nReady for production use!")


if __name__ == "__main__":
    main()