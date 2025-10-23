"""Model installer for managing models dependencies.

This module handles installing required dependencies for different models,
including pip packages, conda packages, and system requirements.
"""

import subprocess
import sys
import json
from pathlib import Path
from typing import List, Optional, Dict, Tuple
import platform

from protos.io.paths import ProtosPaths
from protos.models.model_definitions import get_model_definition, list_available_models


class ModelInstaller:
    """
    Manages installation of models dependencies.
    
    Features:
    - Checks for existing installations
    - Installs pip dependencies
    - Provides conda environment setup instructions
    - Validates system requirements
    - Tracks installed models
    """
    
    def __init__(self, paths: Optional[ProtosPaths] = None):
        """Initialize installer."""
        self.paths = paths or ProtosPaths()
        self.install_registry = self.paths.get_processor_path("models") / ".install_registry.json"
        self._load_registry()
        
        # Detect environment
        self.in_conda = self._detect_conda()
        self.has_gpu = self._detect_gpu()
        self.platform = platform.system().lower()
    
    def _load_registry(self):
        """Load installation registry."""
        if self.install_registry.exists():
            with open(self.install_registry, 'r') as f:
                self.registry = json.load(f)
        else:
            self.registry = {"installed": {}, "environments": {}}
    
    def _save_registry(self):
        """Save installation registry."""
        self.install_registry.parent.mkdir(parents=True, exist_ok=True)
        with open(self.install_registry, 'w') as f:
            json.dump(self.registry, f, indent=2)
    
    def _detect_conda(self) -> bool:
        """Detect if running in conda environment."""
        return 'CONDA_DEFAULT_ENV' in os.environ
    
    def _detect_gpu(self) -> bool:
        """Detect if GPU is available."""
        try:
            import torch
            return torch.cuda.is_available()
        except ImportError:
            # Check nvidia-smi
            try:
                subprocess.run(['nvidia-smi'], capture_output=True, check=True)
                return True
            except:
                return False
    
    def check_requirements(self, model_name: str) -> Dict[str, bool]:
        """
        Check if models requirements are met.
        
        Returns:
            Dict with requirement status
        """
        definition = get_model_definition(model_name)
        status = {
            "pip_dependencies": True,
            "conda_dependencies": True,
            "gpu_available": True,
            "disk_space": True,
            "memory": True
        }
        
        # Check pip dependencies
        for dep in definition.pip_dependencies:
            if not self._check_pip_package(dep):
                status["pip_dependencies"] = False
                break
        
        # Check GPU requirements
        if not definition.requirements.supports_cpu and not self.has_gpu:
            status["gpu_available"] = False
        
        # Check memory requirements
        if definition.requirements.min_gpu_memory_gb and self.has_gpu:
            gpu_memory = self._get_gpu_memory()
            if gpu_memory and gpu_memory < definition.requirements.min_gpu_memory_gb:
                status["memory"] = False
        
        return status
    
    def _check_pip_package(self, package: str) -> bool:
        """Check if pip package is installed."""
        # Handle package specs like "torch>=1.9.0"
        package_name = package.split('[')[0].split('>=')[0].split('==')[0]
        
        try:
            __import__(package_name.replace('-', '_'))
            return True
        except ImportError:
            return False
    
    def _get_gpu_memory(self) -> Optional[float]:
        """Get available GPU memory in GB."""
        try:
            import torch
            if torch.cuda.is_available():
                return torch.cuda.get_device_properties(0).total_memory / (1024**3)
        except:
            pass
        return None
    
    def install_model(self, model_name: str, force: bool = False) -> bool:
        """
        Install dependencies for a models.
        
        Args:
            model_name: Name of the models
            force: Force reinstall even if already installed
            
        Returns:
            True if successful
        """
        definition = get_model_definition(model_name)
        
        # Check if already installed
        if model_name in self.registry["installed"] and not force:
            print(f"{model_name} is already installed")
            return True
        
        print(f"Installing dependencies for {model_name}...")
        
        # Check requirements first
        status = self.check_requirements(model_name)
        
        # Warn about unmet requirements
        if not status["gpu_available"] and not definition.requirements.supports_cpu:
            print(f"WARNING: {model_name} requires GPU but none detected")
            response = input("Continue anyway? [y/N]: ")
            if response.lower() != 'y':
                return False
        
        # Install pip dependencies
        if definition.pip_dependencies:
            success = self._install_pip_dependencies(definition.pip_dependencies)
            if not success:
                return False
        
        # Handle conda dependencies
        if definition.conda_dependencies and self.in_conda:
            self._show_conda_instructions(definition.conda_dependencies)
        
        # Run setup commands
        if definition.setup_commands:
            success = self._run_setup_commands(definition.setup_commands)
            if not success:
                return False
        
        # Update registry
        self.registry["installed"][model_name] = {
            "version": definition.version,
            "installed_at": str(Path.cwd()),
            "pip_deps": definition.pip_dependencies,
            "conda_deps": definition.conda_dependencies
        }
        self._save_registry()
        
        print(f"Successfully installed {model_name}")
        return True
    
    def _install_pip_dependencies(self, dependencies: List[str]) -> bool:
        """Install pip dependencies."""
        print(f"Installing pip packages: {', '.join(dependencies)}")
        
        for dep in dependencies:
            # Special handling for torch with CUDA
            if dep.startswith("torch") and self.has_gpu:
                dep = self._get_torch_cuda_package(dep)
            
            try:
                subprocess.check_call([
                    sys.executable, "-m", "pip", "install", dep
                ])
            except subprocess.CalledProcessError as e:
                print(f"Failed to install {dep}: {e}")
                return False
        
        return True
    
    def _get_torch_cuda_package(self, package: str) -> str:
        """Get appropriate PyTorch package for CUDA."""
        # This is simplified - in practice would detect CUDA version
        if "torchvision" in package:
            return package
        elif "torch-" in package:  # torch-geometric, etc.
            return package
        else:
            # Default to CUDA 11.8
            return f"{package} --index-url https://download.pytorch.org/whl/cu118"
    
    def _show_conda_instructions(self, dependencies: List[str]):
        """Show conda installation instructions."""
        print("\nConda dependencies required:")
        print("Run the following commands in your conda environment:")
        for dep in dependencies:
            print(f"  conda install -c conda-forge {dep}")
        print()
    
    def _run_setup_commands(self, commands: List[str]) -> bool:
        """Run setup commands."""
        print("Running setup commands...")
        
        for cmd in commands:
            print(f"  {cmd}")
            try:
                subprocess.check_call(cmd, shell=True)
            except subprocess.CalledProcessError as e:
                print(f"Command failed: {e}")
                return False
        
        return True
    
    def create_environment(self, model_name: str, env_name: Optional[str] = None) -> str:
        """
        Create conda environment for a models.
        
        Returns:
            Environment setup script
        """
        definition = get_model_definition(model_name)
        
        if not env_name:
            env_name = f"protos_{model_name}"
        
        # Generate environment.yml
        env_config = {
            "name": env_name,
            "channels": ["conda-forge", "pytorch", "nvidia"],
            "dependencies": []
        }
        
        # Add Python
        env_config["dependencies"].append("python=3.9")
        
        # Add conda dependencies
        env_config["dependencies"].extend(definition.conda_dependencies)
        
        # Add pip section
        pip_deps = definition.pip_dependencies.copy()
        if pip_deps:
            env_config["dependencies"].append({
                "pip": pip_deps
            })
        
        # Save environment file
        env_file = self.paths.get_processor_path("models") / f"{model_name}_environment.yml"
        
        import yaml
        with open(env_file, 'w') as f:
            yaml.dump(env_config, f, default_flow_style=False)
        
        # Generate setup script
        script = f"""#!/bin/bash
# Setup script for {model_name}

# Create conda environment
conda env create -f {env_file}

# Activate environment
conda activate {env_name}

# Install Protos in development mode
pip install -e .

echo "Environment '{env_name}' created successfully!"
echo "Activate with: conda activate {env_name}"
"""
        
        script_file = env_file.with_suffix('.sh')
        with open(script_file, 'w') as f:
            f.write(script)
        
        script_file.chmod(0o755)
        
        print(f"Environment configuration saved to: {env_file}")
        print(f"Setup script saved to: {script_file}")
        
        return str(script_file)
    
    def list_installed_models(self) -> List[str]:
        """List installed models."""
        return list(self.registry["installed"].keys())
    
    def uninstall_model(self, model_name: str):
        """Remove models from registry (doesn't uninstall packages)."""
        if model_name in self.registry["installed"]:
            del self.registry["installed"][model_name]
            self._save_registry()
            print(f"Removed {model_name} from registry")
        else:
            print(f"{model_name} not found in registry")
    
    def install_all_models(self, skip_large: bool = True):
        """Install all available models."""
        models = list_available_models()
        
        for model_name in models:
            definition = get_model_definition(model_name)
            
            # Skip large models if requested
            if skip_large:
                total_size = sum(s.size_mb or 0 for s in definition.sources.values())
                if total_size > 5000:  # 5GB
                    print(f"Skipping {model_name} (too large: {total_size}MB)")
                    continue
            
            try:
                self.install_model(model_name)
            except Exception as e:
                print(f"Failed to install {model_name}: {e}")
                continue
    
    def check_all_models(self) -> Dict[str, Dict[str, bool]]:
        """Check requirements for all models."""
        results = {}
        
        for model_name in list_available_models():
            results[model_name] = self.check_requirements(model_name)
        
        return results


def install_model_cli(model_name: str, force: bool = False):
    """CLI function for installing models."""
    installer = ModelInstaller()
    
    # Show models info
    definition = get_model_definition(model_name)
    print(f"\nModel: {definition.full_name} ({definition.version})")
    print(f"Description: {definition.description}")
    
    # Check requirements
    print("\nChecking requirements...")
    status = installer.check_requirements(model_name)
    
    for req, met in status.items():
        symbol = "✓" if met else "✗"
        print(f"  {symbol} {req}")
    
    # Install
    if all(status.values()) or force:
        installer.install_model(model_name, force)
    else:
        print("\nSome requirements are not met.")
        response = input("Install anyway? [y/N]: ")
        if response.lower() == 'y':
            installer.install_model(model_name, force)


if __name__ == "__main__":
    import os
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python model_installer.py <model_name> [--force]")
        print("\nAvailable models:")
        for model in list_available_models():
            print(f"  - {model}")
        sys.exit(1)
    
    model_name = sys.argv[1]
    force = "--force" in sys.argv
    
    install_model_cli(model_name, force)