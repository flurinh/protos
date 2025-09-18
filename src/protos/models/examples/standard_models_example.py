"""Examples of using standard models in Protos.

This example demonstrates:
1. Using pre-defined standard models
2. Automatic dependency installation
3. Model downloading and loading
4. Making predictions with different model types
"""

from pathlib import Path
from protos.models import (
    StandardModel, ModelRegistry, ModelDownloader, ModelInstaller,
    list_available_models, get_model_definition, get_models_by_input,
    InputFormat, OutputFormat
)


def example_list_available_models():
    """Show all available standard models."""
    print("=== Available Standard Models ===\n")
    
    models = list_available_models()
    for model_name in models:
        definition = get_model_definition(model_name)
        print(f"{model_name}:")
        print(f"  Full name: {definition.full_name}")
        print(f"  Description: {definition.description}")
        print(f"  Input formats: {[f.value for f in definition.input_formats]}")
        print(f"  Output format: {definition.output_format.value}")
        print(f"  Framework: {definition.framework.value}")
        print()


def example_find_models_by_capability():
    """Find models based on input/output requirements."""
    print("=== Finding Models by Capability ===\n")
    
    # Find all models that can process sequences
    sequence_models = get_models_by_input(InputFormat.SEQUENCE)
    print(f"Models that accept sequences: {sequence_models}")
    
    # Find all models that output embeddings
    embedding_models = get_models_by_output(OutputFormat.EMBEDDING)
    print(f"Models that output embeddings: {embedding_models}")
    
    # Find all structure prediction models
    structure_models = get_models_by_output(OutputFormat.STRUCTURE)
    print(f"Structure prediction models: {structure_models}")
    print()


def example_install_and_download_model():
    """Install dependencies and download a model."""
    print("=== Installing and Downloading ESM-2 ===\n")
    
    # Initialize installer and downloader
    installer = ModelInstaller()
    downloader = ModelDownloader()
    
    # Check requirements for ESM-2
    print("Checking requirements for ESM-2...")
    status = installer.check_requirements("esm2")
    for req, met in status.items():
        print(f"  {'✓' if met else '✗'} {req}")
    
    # Install dependencies (in practice, this would actually install)
    print("\nInstalling dependencies...")
    # installer.install_model("esm2")
    print("(Skipping actual installation in example)")
    
    # Download smallest ESM-2 variant
    print("\nDownloading ESM-2 (8M parameters)...")
    # weights_path = downloader.download_model("esm2", "esm2_t6_8M", Path("models/esm2/weights"))
    print("(Skipping actual download in example)")
    print()


def example_use_esm2_model():
    """Use ESM-2 for sequence embeddings."""
    print("=== Using ESM-2 Model ===\n")
    
    # Create ESM-2 model
    model = StandardModel(
        model_name="esm2",
        model_variant="esm2_t6_8M",
        device="cpu"  # Use CPU for example
    )
    
    # Show model info
    info = model.get_model_info()
    print(f"Model: {info['full_name']} ({info['variant']})")
    print(f"Framework: {info['framework']}")
    print(f"Device: {info['device']}")
    
    # Load model (would load weights in practice)
    print("\nLoading model weights...")
    # model.load_model()
    print("(Skipping actual loading in example)")
    
    # Make prediction on a sequence
    sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"
    
    print(f"\nProcessing sequence (length: {len(sequence)})")
    # prediction = model.predict("my_protein")
    # embeddings = prediction.prediction['embeddings']
    print("(Would return embeddings of shape [seq_len, 320])")
    print()


def example_use_ankh_model():
    """Use Ankh for sequence embeddings."""
    print("=== Using Ankh Model ===\n")
    
    # Create Ankh model
    model = StandardModel(
        model_name="ankh",
        model_variant="base"
    )
    
    print("Ankh is optimized for efficiency while maintaining high quality embeddings")
    print("It uses a T5-based architecture and is particularly good for long sequences")
    
    # Would use similar to ESM-2
    # model.load_model()
    # prediction = model.predict("protein_name")
    print()


def example_use_boltz_model():
    """Use Boltz for structure prediction."""
    print("=== Using Boltz-1 Model ===\n")
    
    # Create Boltz model
    model = StandardModel(
        model_name="boltz1",
        model_variant="model"
    )
    
    definition = get_model_definition("boltz1")
    print(f"Boltz-1 is a structure prediction model")
    print(f"Input formats: {[f.value for f in definition.input_formats]}")
    print(f"Output: {definition.output_format.value}")
    print(f"GPU Required: {not definition.requirements.supports_cpu}")
    
    # Structure prediction workflow
    print("\nStructure prediction workflow:")
    print("1. Load sequence (and optionally MSA)")
    print("2. Run Boltz-1 prediction")
    print("3. Save predicted structure as PDB/CIF")
    print("4. Structure automatically registered in Protos")
    
    # Example usage (would require GPU)
    # model.load_model()
    # prediction = model.predict("target_protein")
    # structure = prediction.prediction  # Would be structure coordinates
    print()


def example_model_registry_usage():
    """Use the model registry for management."""
    print("=== Model Registry ===\n")
    
    registry = ModelRegistry()
    
    # Discover models
    print("Discovering available models...")
    models = registry.discover_models()
    print(f"Found {len(models)} model configurations")
    
    # Create a new model instance
    print("\nCreating ESM-2 model through registry...")
    # model = registry.create_model(
    #     model_type="standard",
    #     model_name="esm2_small",
    #     config=ModelConfig(...)
    # )
    
    # Load existing model
    print("Loading model from registry...")
    # model = registry.load_model("standard/esm2_small")
    
    # Find compatible models for an entity
    print("\nFinding models compatible with a protein entity...")
    # compatible = registry.get_models_for_entity("BACR_HALSA", output_format="embedding")
    # print(f"Compatible models: {compatible}")
    print()


def example_custom_model_with_standard_interface():
    """Create a custom model using the standard interface."""
    print("=== Custom Model with Standard Interface ===\n")
    
    from protos.models import ModelDefinition, ModelSource, ModelRequirements
    from protos.models.model_definitions import STANDARD_MODELS, ModelFramework
    
    # Define a custom model
    custom_definition = ModelDefinition(
        name="my_custom_model",
        full_name="My Custom Model",
        version="1.0",
        description="Custom property predictor",
        framework=ModelFramework.PYTORCH,
        input_formats=[InputFormat.SEQUENCE, InputFormat.STRUCTURE],
        output_format=OutputFormat.PROPERTY,
        sources={
            "default": ModelSource(
                url="https://myserver.com/model.pth",
                size_mb=100
            )
        },
        pip_dependencies=["torch", "numpy"],
        requirements=ModelRequirements(
            min_gpu_memory_gb=4,
            supports_cpu=True
        ),
        preprocessing_config={
            "max_length": 512,
            "use_plddt": True
        }
    )
    
    # Register the custom model
    STANDARD_MODELS["my_custom_model"] = custom_definition
    
    # Now can use it like any standard model
    model = StandardModel("my_custom_model")
    print(f"Created custom model: {model.name}")
    print(f"Input formats: {model.config.input_formats}")
    print(f"Output format: {model.config.output_format}")
    print()


def example_batch_processing():
    """Process multiple proteins with a model."""
    print("=== Batch Processing ===\n")
    
    # Create model
    model = StandardModel("esm2", "esm2_t6_8M")
    
    # Create dataset of proteins to process
    print("Creating dataset for batch processing...")
    proteins = ["BACR_HALSA", "ChR2", "NpHR", "C1C2"]
    
    model.create_dataset(
        "opsin_embeddings",
        proteins,
        metadata={"purpose": "embedding generation", "model": "esm2_t6_8M"}
    )
    
    # Process entire dataset
    print("Processing dataset...")
    # model.load_model()
    # predictions = model.predict_dataset("opsin_embeddings", batch_size=2)
    
    print(f"Would process {len(proteins)} proteins")
    print("Results automatically saved and tracked in entity registry")
    print()


def main():
    """Run all examples."""
    print("=== Protos Standard Models Examples ===\n")
    
    example_list_available_models()
    example_find_models_by_capability()
    example_install_and_download_model()
    example_use_esm2_model()
    example_use_ankh_model()
    example_use_boltz_model()
    example_model_registry_usage()
    example_custom_model_with_standard_interface()
    example_batch_processing()
    
    print("=== Summary ===\n")
    print("The standard model system provides:")
    print("- Pre-configured definitions for popular models")
    print("- Automatic dependency management")
    print("- Model weight downloading from official sources")
    print("- Unified interface for all models")
    print("- Automatic format conversion")
    print("- Integration with Protos entity system")
    print("\nNo need to write model-specific code!")


if __name__ == "__main__":
    main()