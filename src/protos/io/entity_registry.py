"""
Entity Registry for tracking biological objects across multiple formats.

This module provides centralized tracking of all biological entities in the system,
maintaining relationships between different data formats for the same entity.
"""

import json
import hashlib
import os
import time
import tempfile
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, asdict

from .paths import ProtosPaths


@dataclass
class EntityInfo:
    """Information about an entity in a specific format."""
    hash_id: str
    original_id: str
    format_type: str
    file_path: str
    metadata: Dict[str, Any]
    created: str
    modified: str


@dataclass
class Entity:
    """Complete entity with all formats."""
    hash_id: str
    original_id: str
    aliases: List[str]
    formats: Dict[str, Dict[str, Any]]
    created: str
    modified: str


class EntityRegistry:
    """
    Central registry for tracking entities across all formats.
    
    Key principles:
    - Hash IDs are internal only - never exposed to users
    - All public methods work with human-readable names
    - File paths always use human-readable names
    """
    
    def __init__(self, paths: Optional[ProtosPaths] = None):
        """
        Initialize EntityRegistry with ProtosPaths.
        
        Args:
            paths: ProtosPaths instance for path management
        """
        self.paths = paths or ProtosPaths()
        self.registry_file = Path(self.paths.get_global_registry_path())
        self._registry = self._load_registry()
        
    def _load_registry(self) -> Dict[str, Dict]:
        """Load registry from JSON file."""
        if self.registry_file.exists():
            with open(self.registry_file, 'r') as f:
                data = json.load(f)
                return data.get('entities', {})
        return {}
    
    def _save_registry(self):
        """Save registry to JSON file with atomic write and retry logic."""
        # Ensure parent directory exists
        self.registry_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Atomic write: write to temp file then rename
        temp_fd, temp_path = tempfile.mkstemp(
            dir=self.registry_file.parent,
            prefix='.tmp_registry_',
            suffix='.json'
        )
        
        try:
            # Write to temp file
            with os.fdopen(temp_fd, 'w') as f:
                json.dump({'entities': self._registry}, f, indent=2)
            
            # Atomic rename with retry logic for Windows
            max_retries = 5
            retry_delay = 0.1
            
            for attempt in range(max_retries):
                try:
                    # On Windows, need to remove target if it exists
                    if os.name == 'nt' and self.registry_file.exists():
                        self.registry_file.unlink()
                    
                    # Atomic rename
                    shutil.move(temp_path, str(self.registry_file))
                    break
                    
                except (OSError, PermissionError) as e:
                    if attempt < max_retries - 1:
                        time.sleep(retry_delay * (attempt + 1))
                    else:
                        raise RuntimeError(
                            f"Failed to save registry after {max_retries} attempts: {e}"
                        )
        except:
            # Clean up temp file on error
            if os.path.exists(temp_path):
                os.unlink(temp_path)
            raise
    
    def refresh(self):
        """Reload registry from disk to see updates from other processes."""
        self._registry = self._load_registry()
    
    def _generate_hash_id(self, name: str) -> str:
        """
        Generate hash ID for internal use only.
        
        Args:
            name: Human-readable name
            
        Returns:
            10-character hash ID
        """
        # Use SHA256 and take first 10 characters
        hash_obj = hashlib.sha256(name.encode('utf-8'))
        return hash_obj.hexdigest()[:10]
    
    def _resolve_to_hash(self, name: str) -> Optional[str]:
        """
        Convert human name to hash ID for internal lookups.
        
        Args:
            name: Human-readable name or alias
            
        Returns:
            Hash ID if found, None otherwise
        """
        # First check if it's already a hash ID (shouldn't happen in public API)
        if name in self._registry:
            return name
            
        # Search by original_id
        for hash_id, entity in self._registry.items():
            if entity.get('original_id') == name:
                return hash_id
                
        # Search by aliases
        for hash_id, entity in self._registry.items():
            if name in entity.get('aliases', []):
                return hash_id
                
        return None
    
    def register_entity(self, name: str, format_type: str, file_path: str, 
                       metadata: Optional[Dict[str, Any]] = None) -> str:
        """
        Register an entity - users provide human-readable name.
        
        Args:
            name: Human-readable entity name
            format_type: Format type (e.g., 'structure', 'sequence')
            file_path: Path to file (should use human-readable name)
            metadata: Optional metadata dict
            
        Returns:
            The human-readable name (not hash ID)
        """
        metadata = metadata or {}
        
        # Check if entity already exists in memory
        hash_id = self._resolve_to_hash(name)
        
        if hash_id is None:
            # Not in memory - refresh and check again
            self.refresh()
            hash_id = self._resolve_to_hash(name)
            
            if hash_id is None:
                # Still not found - new entity
                hash_id = self._generate_hash_id(name)
                self._registry[hash_id] = {
                    'original_id': name,
                    'aliases': [],
                    'formats': {},
                    'created': datetime.now().isoformat(),
                    'modified': datetime.now().isoformat()
                }
        
        # Update format information
        self._registry[hash_id]['formats'][format_type] = {
            'file_path': file_path,
            'metadata': metadata,
            'created': datetime.now().isoformat()
        }
        self._registry[hash_id]['modified'] = datetime.now().isoformat()
        
        self._save_registry()
        return name  # Always return human-readable name
    
    def find_entity(self, name: str, format_type: Optional[str] = None) -> Optional[EntityInfo]:
        """
        Find entity by human-readable name or alias.
        
        Args:
            name: Human-readable name or alias
            format_type: Optional format filter
            
        Returns:
            EntityInfo if found, None otherwise
        """
        hash_id = self._resolve_to_hash(name)
        if hash_id is None:
            return None
            
        entity = self._registry[hash_id]
        
        # If format specified, check if it exists
        if format_type:
            if format_type not in entity['formats']:
                return None
            format_data = entity['formats'][format_type]
            return EntityInfo(
                hash_id=hash_id,
                original_id=entity['original_id'],
                format_type=format_type,
                file_path=format_data['file_path'],
                metadata=format_data.get('metadata', {}),
                created=format_data.get('created', entity['created']),
                modified=entity['modified']
            )
        
        # Return info for first available format
        if entity['formats']:
            format_type = list(entity['formats'].keys())[0]
            format_data = entity['formats'][format_type]
            return EntityInfo(
                hash_id=hash_id,
                original_id=entity['original_id'],
                format_type=format_type,
                file_path=format_data['file_path'],
                metadata=format_data.get('metadata', {}),
                created=format_data.get('created', entity['created']),
                modified=entity['modified']
            )
            
        return None
    
    def list_entities(self, format_type: Optional[str] = None) -> List[str]:
        """
        List all entities - returns human-readable names.
        
        Args:
            format_type: Optional filter by format type
            
        Returns:
            List of human-readable entity names
        """
        names = []
        for entity in self._registry.values():
            if format_type is None or format_type in entity['formats']:
                names.append(entity['original_id'])
        return sorted(names)
    
    def get_entity_formats(self, name: str) -> List[str]:
        """
        Get all available formats for an entity.
        
        Args:
            name: Human-readable entity name
            
        Returns:
            List of format types
        """
        hash_id = self._resolve_to_hash(name)
        if hash_id is None:
            return []
            
        return list(self._registry[hash_id]['formats'].keys())
    
    def add_alias(self, name: str, alias: str):
        """
        Add an alias for an entity.
        
        Args:
            name: Original entity name
            alias: Alternative name to add
        """
        hash_id = self._resolve_to_hash(name)
        if hash_id is None:
            raise ValueError(f"Entity '{name}' not found")
            
        aliases = self._registry[hash_id].get('aliases', [])
        if alias not in aliases and alias != self._registry[hash_id]['original_id']:
            aliases.append(alias)
            self._registry[hash_id]['aliases'] = aliases
            self._registry[hash_id]['modified'] = datetime.now().isoformat()
            self._save_registry()
    
    def remove_format(self, name: str, format_type: str):
        """
        Remove a format from an entity.
        
        Args:
            name: Human-readable entity name
            format_type: Format to remove
        """
        hash_id = self._resolve_to_hash(name)
        if hash_id is None:
            return
            
        if format_type in self._registry[hash_id]['formats']:
            del self._registry[hash_id]['formats'][format_type]
            self._registry[hash_id]['modified'] = datetime.now().isoformat()
            
            # If no formats left, remove entity entirely
            if not self._registry[hash_id]['formats']:
                del self._registry[hash_id]
                
            self._save_registry()
    
    def entity_exists(self, name: str, format_type: Optional[str] = None) -> bool:
        """
        Check if entity exists.
        
        Args:
            name: Human-readable entity name
            format_type: Optional format type to check
            
        Returns:
            True if entity exists (in specified format if given)
        """
        hash_id = self._resolve_to_hash(name)
        if hash_id is None:
            return False
            
        if format_type:
            return format_type in self._registry[hash_id]['formats']
            
        return True
    
    def get_entity_metadata(self, name: str, format_type: str) -> Dict[str, Any]:
        """
        Get metadata for a specific entity format.
        
        Args:
            name: Human-readable entity name
            format_type: Format type
            
        Returns:
            Metadata dictionary
        """
        info = self.find_entity(name, format_type)
        if info:
            return info.metadata
        return {}
    
    def update_metadata(self, name: str, format_type: str, metadata: Dict[str, Any]):
        """
        Update metadata for an entity format.
        
        Args:
            name: Human-readable entity name
            format_type: Format type
            metadata: New metadata to merge
        """
        hash_id = self._resolve_to_hash(name)
        if hash_id and format_type in self._registry[hash_id]['formats']:
            current = self._registry[hash_id]['formats'][format_type].get('metadata', {})
            current.update(metadata)
            self._registry[hash_id]['formats'][format_type]['metadata'] = current
            self._registry[hash_id]['modified'] = datetime.now().isoformat()
            self._save_registry()