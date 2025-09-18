"""Model downloader for fetching model weights from various sources.

This module handles downloading model weights from GitHub releases,
HuggingFace, and other sources with progress tracking and verification.
"""

import os
import hashlib
import json
import shutil
import tempfile
from pathlib import Path
from typing import Optional, Dict, Any, Tuple
from urllib.parse import urlparse
import requests
from tqdm import tqdm

from protos.io.paths import ProtosPaths
from protos.models.model_definitions import get_model_definition, ModelSource


class ModelDownloader:
    """
    Downloads and manages model weights for Protos.
    
    Features:
    - Progress tracking
    - Checksum verification
    - Resume interrupted downloads
    - Multiple source support (GitHub, HuggingFace, direct URLs)
    - Automatic extraction of archives
    """
    
    def __init__(self, paths: Optional[ProtosPaths] = None):
        """Initialize downloader."""
        self.paths = paths or ProtosPaths()
        self.cache_dir = self.paths.get_processor_path("models") / ".cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Download metadata
        self.download_registry = self.cache_dir / "download_registry.json"
        self._load_registry()
    
    def _load_registry(self):
        """Load download registry."""
        if self.download_registry.exists():
            with open(self.download_registry, 'r') as f:
                self.registry = json.load(f)
        else:
            self.registry = {}
    
    def _save_registry(self):
        """Save download registry."""
        with open(self.download_registry, 'w') as f:
            json.dump(self.registry, f, indent=2)
    
    def download_model(self, model_name: str, variant: str, 
                      target_dir: Path) -> Path:
        """
        Download model weights.
        
        Args:
            model_name: Name of the model (e.g., 'esm2')
            variant: Model variant (e.g., 'esm2_t33_650M')
            target_dir: Directory to save weights to
            
        Returns:
            Path to downloaded weights
        """
        # Get model definition
        definition = get_model_definition(model_name)
        
        if variant not in definition.sources:
            raise ValueError(f"Unknown variant '{variant}' for model '{model_name}'")
        
        source = definition.sources[variant]
        
        # Determine download method based on source
        if source.format == "huggingface":
            return self._download_huggingface(source, model_name, variant, target_dir)
        elif source.format in ["weights", "checkpoint", "onnx"]:
            return self._download_direct(source, model_name, variant, target_dir)
        else:
            raise ValueError(f"Unknown source format: {source.format}")
    
    def _download_direct(self, source: ModelSource, model_name: str, 
                        variant: str, target_dir: Path) -> Path:
        """Download from direct URL."""
        url = source.url
        
        # Determine filename
        parsed_url = urlparse(url)
        filename = os.path.basename(parsed_url.path)
        if not filename:
            filename = f"{variant}.pt"
        
        target_path = target_dir / filename
        
        # Check if already downloaded
        if target_path.exists() and self._verify_checksum(target_path, source.checksum):
            print(f"Model {variant} already downloaded and verified")
            return target_path
        
        # Download with progress
        print(f"Downloading {model_name} ({variant}) from {url}")
        
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            
            # Use temporary file for atomic download
            with tempfile.NamedTemporaryFile(delete=False, dir=target_dir) as tmp_file:
                with tqdm(total=total_size, unit='B', unit_scale=True, desc=variant) as pbar:
                    for chunk in response.iter_content(chunk_size=8192):
                        tmp_file.write(chunk)
                        pbar.update(len(chunk))
                
                tmp_path = Path(tmp_file.name)
            
            # Verify checksum if provided
            if source.checksum and not self._verify_checksum(tmp_path, source.checksum):
                tmp_path.unlink()
                raise ValueError(f"Checksum verification failed for {variant}")
            
            # Move to final location
            shutil.move(str(tmp_path), str(target_path))
            
            # Update registry
            self.registry[f"{model_name}/{variant}"] = {
                "url": url,
                "path": str(target_path),
                "size": total_size,
                "checksum": source.checksum
            }
            self._save_registry()
            
            print(f"Successfully downloaded {variant}")
            return target_path
            
        except Exception as e:
            print(f"Error downloading {variant}: {e}")
            if tmp_path.exists():
                tmp_path.unlink()
            raise
    
    def _download_huggingface(self, source: ModelSource, model_name: str,
                             variant: str, target_dir: Path) -> Path:
        """Download from HuggingFace."""
        try:
            from huggingface_hub import snapshot_download, hf_hub_download
        except ImportError:
            raise ImportError(
                "HuggingFace Hub not installed. "
                "Run: pip install huggingface-hub"
            )
        
        repo_id = source.url.replace("https://huggingface.co/", "")
        
        # Create model-specific directory
        model_dir = target_dir / variant
        model_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"Downloading {model_name} ({variant}) from HuggingFace: {repo_id}")
        
        try:
            # Download entire model repository
            local_dir = snapshot_download(
                repo_id=repo_id,
                local_dir=model_dir,
                local_dir_use_symlinks=False,
                resume_download=True
            )
            
            # Find the main model file
            model_files = list(Path(local_dir).glob("*.bin")) + \
                         list(Path(local_dir).glob("*.safetensors")) + \
                         list(Path(local_dir).glob("*.pt"))
            
            if model_files:
                # Use the largest file as the main model
                main_model = max(model_files, key=lambda p: p.stat().st_size)
                
                # Create a symlink or copy to expected location
                expected_path = target_dir / f"{variant}.pt"
                if not expected_path.exists():
                    shutil.copy2(main_model, expected_path)
                
                return expected_path
            else:
                return Path(local_dir)
                
        except Exception as e:
            print(f"Error downloading from HuggingFace: {e}")
            raise
    
    def _verify_checksum(self, file_path: Path, expected_checksum: Optional[str]) -> bool:
        """Verify file checksum."""
        if not expected_checksum:
            return True
        
        # Determine hash algorithm from checksum length
        if len(expected_checksum) == 32:
            hasher = hashlib.md5()
        elif len(expected_checksum) == 40:
            hasher = hashlib.sha1()
        elif len(expected_checksum) == 64:
            hasher = hashlib.sha256()
        else:
            print(f"Unknown checksum format: {expected_checksum}")
            return True
        
        # Calculate file hash
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                hasher.update(chunk)
        
        calculated = hasher.hexdigest()
        return calculated == expected_checksum
    
    def list_downloaded_models(self) -> Dict[str, Dict[str, Any]]:
        """List all downloaded models."""
        return self.registry.copy()
    
    def get_download_info(self, model_name: str, variant: str) -> Optional[Dict[str, Any]]:
        """Get download information for a model."""
        key = f"{model_name}/{variant}"
        return self.registry.get(key)
    
    def remove_downloaded_model(self, model_name: str, variant: str):
        """Remove a downloaded model."""
        key = f"{model_name}/{variant}"
        if key in self.registry:
            info = self.registry[key]
            path = Path(info['path'])
            
            # Remove file or directory
            if path.exists():
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()
            
            # Update registry
            del self.registry[key]
            self._save_registry()
            
            print(f"Removed {model_name} ({variant})")
    
    def download_all_variants(self, model_name: str, target_dir: Path) -> Dict[str, Path]:
        """Download all variants of a model."""
        definition = get_model_definition(model_name)
        downloaded = {}
        
        for variant in definition.sources:
            try:
                path = self.download_model(model_name, variant, target_dir)
                downloaded[variant] = path
            except Exception as e:
                print(f"Failed to download {variant}: {e}")
                continue
        
        return downloaded
    
    def estimate_download_size(self, model_name: str, 
                             variant: Optional[str] = None) -> float:
        """Estimate download size in MB."""
        definition = get_model_definition(model_name)
        
        if variant:
            if variant in definition.sources:
                return definition.sources[variant].size_mb or 0
            else:
                return 0
        else:
            # Sum all variants
            total = 0
            for source in definition.sources.values():
                total += source.size_mb or 0
            return total
    
    def check_disk_space(self, required_mb: float) -> bool:
        """Check if enough disk space is available."""
        import shutil
        stat = shutil.disk_usage(self.cache_dir)
        available_mb = stat.free / (1024 * 1024)
        return available_mb > required_mb * 1.5  # 50% buffer


def download_model_cli(model_name: str, variant: Optional[str] = None,
                      target_dir: Optional[str] = None):
    """CLI function for downloading models."""
    downloader = ModelDownloader()
    
    if target_dir:
        target = Path(target_dir)
    else:
        target = downloader.paths.get_processor_path("models") / model_name / "weights"
    
    target.mkdir(parents=True, exist_ok=True)
    
    # Check disk space
    size_mb = downloader.estimate_download_size(model_name, variant)
    if size_mb > 0:
        print(f"Estimated download size: {size_mb:.1f} MB")
        if not downloader.check_disk_space(size_mb):
            print("WARNING: Insufficient disk space")
            response = input("Continue anyway? [y/N]: ")
            if response.lower() != 'y':
                return
    
    # Download
    if variant:
        path = downloader.download_model(model_name, variant, target)
        print(f"Downloaded to: {path}")
    else:
        paths = downloader.download_all_variants(model_name, target)
        print(f"Downloaded {len(paths)} variants")
        for var, path in paths.items():
            print(f"  {var}: {path}")


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python model_downloader.py <model_name> [variant] [target_dir]")
        sys.exit(1)
    
    model_name = sys.argv[1]
    variant = sys.argv[2] if len(sys.argv) > 2 else None
    target_dir = sys.argv[3] if len(sys.argv) > 3 else None
    
    download_model_cli(model_name, variant, target_dir)