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
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, asdict

from .paths import ProtosPaths


# Define relationship types with their directionality
RELATIONSHIP_TYPES = {
    "derived_from": {"inverse": "derives_to", "symmetric": False},
    "subset_of": {"inverse": "contains", "symmetric": False},
    "merged_from": {"inverse": "merged_into", "symmetric": False},
    "version_of": {"inverse": "has_version", "symmetric": False},
    "aligned_to": {"inverse": "aligned_to", "symmetric": True},
    "annotated_by": {"inverse": "annotates", "symmetric": False},
}


@dataclass
class Relationship:
    """A directed relationship between two entities."""
    type: str  # Relationship type (e.g., 'derived_from')
    source: str  # Source entity UUID
    target: str  # Target entity UUID (usually self)
    metadata: Dict[str, Any]  # Additional relationship data
    created: str  # When relationship was created


@dataclass
class EntityInfo:
    """Information about an entity in a specific format."""
    entity_id: str  # UUID (internal use only)
    original_id: str  # Human-readable name
    format_type: str
    file_path: str
    metadata: Dict[str, Any]
    created: str
    modified: str


@dataclass
class Entity:
    """Complete entity with all formats."""
    entity_id: str  # UUID (internal use only)
    original_id: str  # Human-readable name
    aliases: List[str]
    formats: Dict[str, Dict[str, Any]]
    relationships: List[Dict[str, Any]]  # List of relationship objects
    created: str
    modified: str


class EntityRegistry:
    """
    Central registry for tracking entities across all formats.
    
    Key principles:
    - Entity IDs (UUIDs) are internal only - never exposed to users
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
        self._registry, self._name_index = self._load_registry()
        self._build_name_index()
        
    def _load_registry(self) -> Tuple[Dict[str, Dict], Dict[str, str]]:
        """Load registry from JSON file.
        
        Returns:
            Tuple of (entities dict, name_index dict)
        """
        if self.registry_file.exists():
            with open(self.registry_file, 'r') as f:
                data = json.load(f)
                entities = data.get('entities', {})
                name_index = data.get('name_index', {})
                return entities, name_index
        return {}, {}
    
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
                json.dump({
                    'entities': self._registry,
                    'name_index': self._name_index
                }, f, indent=2)
            
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
        self._registry, self._name_index = self._load_registry()
        self._build_name_index()
    
    def _build_name_index(self):
        """Build name-to-ID index for fast lookups."""
        # Clear existing index
        self._name_index = {}
        
        # Build index from entities
        for entity_id, entity in self._registry.items():
            # Add original ID
            original_id = entity.get('original_id')
            if original_id:
                self._name_index[original_id] = entity_id
                
            # Add aliases
            for alias in entity.get('aliases', []):
                self._name_index[alias] = entity_id
    
    def _generate_entity_id(self, name: str = None) -> str:
        """
        Generate UUID for internal use only.
        
        Args:
            name: Human-readable name (no longer used, kept for compatibility)
            
        Returns:
            UUID string (36 characters)
        """
        # Generate a new UUID independent of the name
        return str(uuid.uuid4())
    
    def _resolve_to_id(self, name: str) -> Optional[str]:
        """
        Convert human name to entity ID (UUID) for internal lookups.
        
        Args:
            name: Human-readable name or alias
            
        Returns:
            Entity ID (UUID) if found, None otherwise
        """
        # First check if it's already an entity ID (shouldn't happen in public API)
        if name in self._registry:
            return name
            
        # Use name index for O(1) lookup
        return self._name_index.get(name)
    
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
            The human-readable name (not entity ID)
        """
        metadata = metadata or {}
        
        # Check if entity already exists in memory
        entity_id = self._resolve_to_id(name)
        
        if entity_id is None:
            # Not in memory - refresh and check again
            self.refresh()
            entity_id = self._resolve_to_id(name)
            
            if entity_id is None:
                # Still not found - new entity
                entity_id = self._generate_entity_id()  # UUID generation, no name needed
                self._registry[entity_id] = {
                    'original_id': name,
                    'aliases': [],
                    'formats': {},
                    'relationships': [],  # Initialize empty relationships list
                    'created': datetime.now().isoformat(),
                    'modified': datetime.now().isoformat()
                }
                # Update name index
                self._name_index[name] = entity_id
        
        # Update format information
        self._registry[entity_id]['formats'][format_type] = {
            'file_path': file_path,
            'metadata': metadata,
            'created': datetime.now().isoformat()
        }
        self._registry[entity_id]['modified'] = datetime.now().isoformat()
        
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
        entity_id = self._resolve_to_id(name)
        if entity_id is None:
            return None
            
        entity = self._registry[entity_id]
        
        # If format specified, check if it exists
        if format_type:
            if format_type not in entity['formats']:
                return None
            format_data = entity['formats'][format_type]
            return EntityInfo(
                entity_id=entity_id,
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
                entity_id=entity_id,
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
        entity_id = self._resolve_to_id(name)
        if entity_id is None:
            return []
            
        return list(self._registry[entity_id]['formats'].keys())
    
    def add_alias(self, name: str, alias: str):
        """
        Add an alias for an entity.
        
        Args:
            name: Original entity name
            alias: Alternative name to add
        """
        entity_id = self._resolve_to_id(name)
        if entity_id is None:
            raise ValueError(f"Entity '{name}' not found")
            
        aliases = self._registry[entity_id].get('aliases', [])
        if alias not in aliases and alias != self._registry[entity_id]['original_id']:
            aliases.append(alias)
            self._registry[entity_id]['aliases'] = aliases
            self._registry[entity_id]['modified'] = datetime.now().isoformat()
            # Update name index
            self._name_index[alias] = entity_id
            self._save_registry()
    
    def remove_format(self, name: str, format_type: str):
        """
        Remove a format from an entity.
        
        Args:
            name: Human-readable entity name
            format_type: Format to remove
        """
        entity_id = self._resolve_to_id(name)
        if entity_id is None:
            return
            
        if format_type in self._registry[entity_id]['formats']:
            del self._registry[entity_id]['formats'][format_type]
            self._registry[entity_id]['modified'] = datetime.now().isoformat()
            
            # If no formats left, remove entity entirely
            if not self._registry[entity_id]['formats']:
                # Remove from name index
                original_id = self._registry[entity_id].get('original_id')
                if original_id and original_id in self._name_index:
                    del self._name_index[original_id]
                for alias in self._registry[entity_id].get('aliases', []):
                    if alias in self._name_index:
                        del self._name_index[alias]
                # Remove entity
                del self._registry[entity_id]
                
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
        entity_id = self._resolve_to_id(name)
        if entity_id is None:
            return False
            
        if format_type:
            return format_type in self._registry[entity_id]['formats']
            
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
        entity_id = self._resolve_to_id(name)
        if entity_id and format_type in self._registry[entity_id]['formats']:
            current = self._registry[entity_id]['formats'][format_type].get('metadata', {})
            current.update(metadata)
            self._registry[entity_id]['formats'][format_type]['metadata'] = current
            self._registry[entity_id]['modified'] = datetime.now().isoformat()
            self._save_registry()
    
    # Relationship management methods
    
    def add_relationship(self, source_name: str, target_name: str, 
                        rel_type: str, metadata: Optional[Dict[str, Any]] = None):
        """
        Add a directed relationship between two entities.
        
        Args:
            source_name: Human-readable name of source entity
            target_name: Human-readable name of target entity
            rel_type: Relationship type (must be in RELATIONSHIP_TYPES)
            metadata: Optional metadata for the relationship
            
        Raises:
            ValueError: If relationship type is invalid or entities don't exist
        """
        if rel_type not in RELATIONSHIP_TYPES:
            raise ValueError(f"Invalid relationship type: {rel_type}")
            
        # Resolve names to IDs
        source_id = self._resolve_to_id(source_name)
        target_id = self._resolve_to_id(target_name)
        
        if not source_id:
            raise ValueError(f"Source entity '{source_name}' not found")
        if not target_id:
            raise ValueError(f"Target entity '{target_name}' not found")
            
        # Create relationship object
        relationship = {
            'type': rel_type,
            'source': source_id,
            'target': target_id,
            'metadata': metadata or {},
            'created': datetime.now().isoformat()
        }
        
        # Add to target entity's relationships
        if 'relationships' not in self._registry[target_id]:
            self._registry[target_id]['relationships'] = []
            
        # Check if relationship already exists
        existing = [r for r in self._registry[target_id]['relationships']
                   if r['type'] == rel_type and r['source'] == source_id]
        
        if not existing:
            self._registry[target_id]['relationships'].append(relationship)
            self._registry[target_id]['modified'] = datetime.now().isoformat()
            self._save_registry()
    
    def get_relationships(self, name: str, rel_type: Optional[str] = None,
                         direction: str = "both") -> List[Dict[str, Any]]:
        """
        Get relationships for an entity, optionally filtered by type and direction.
        
        Args:
            name: Human-readable entity name
            rel_type: Optional relationship type filter
            direction: "incoming", "outgoing", or "both"
            
        Returns:
            List of relationship dictionaries with resolved names
        """
        entity_id = self._resolve_to_id(name)
        if not entity_id:
            return []
            
        relationships = []
        
        # Get incoming relationships (where this entity is the target)
        if direction in ("incoming", "both"):
            for rel in self._registry[entity_id].get('relationships', []):
                if rel_type and rel['type'] != rel_type:
                    continue
                    
                # Resolve source name
                source_name = self._registry.get(rel['source'], {}).get('original_id', 'Unknown')
                
                relationships.append({
                    'type': rel['type'],
                    'source_name': source_name,
                    'target_name': name,
                    'direction': 'incoming',
                    'metadata': rel.get('metadata', {}),
                    'created': rel.get('created')
                })
        
        # Get outgoing relationships (where this entity is the source)
        if direction in ("outgoing", "both"):
            for other_id, other_entity in self._registry.items():
                for rel in other_entity.get('relationships', []):
                    if rel['source'] == entity_id:
                        if rel_type and rel['type'] != rel_type:
                            continue
                            
                        # For outgoing, we need the inverse type
                        rel_info = RELATIONSHIP_TYPES.get(rel['type'], {})
                        display_type = rel_info.get('inverse', rel['type'])
                        
                        relationships.append({
                            'type': display_type,
                            'source_name': name,
                            'target_name': other_entity['original_id'],
                            'direction': 'outgoing',
                            'metadata': rel.get('metadata', {}),
                            'created': rel.get('created')
                        })
        
        return relationships
    
    def get_related_entities(self, name: str, rel_type: Optional[str] = None,
                           direction: str = "both") -> List[str]:
        """
        Return related entity names (resolved from IDs).
        
        Args:
            name: Human-readable entity name
            rel_type: Optional relationship type filter
            direction: "incoming", "outgoing", or "both"
            
        Returns:
            List of human-readable entity names
        """
        relationships = self.get_relationships(name, rel_type, direction)
        
        related = []
        for rel in relationships:
            if rel['direction'] == 'incoming':
                related.append(rel['source_name'])
            else:
                related.append(rel['target_name'])
                
        return list(set(related))  # Remove duplicates
    
    def remove_relationship(self, source_name: str, target_name: str, rel_type: str):
        """
        Remove a specific relationship between two entities.
        
        Args:
            source_name: Human-readable name of source entity
            target_name: Human-readable name of target entity  
            rel_type: Relationship type to remove
        """
        # Resolve names to IDs
        source_id = self._resolve_to_id(source_name)
        target_id = self._resolve_to_id(target_name)
        
        if not source_id or not target_id:
            return  # Silently ignore if entities don't exist
            
        # Remove from target entity's relationships
        if target_id in self._registry:
            relationships = self._registry[target_id].get('relationships', [])
            updated = [r for r in relationships 
                      if not (r['type'] == rel_type and r['source'] == source_id)]
            
            if len(updated) < len(relationships):
                self._registry[target_id]['relationships'] = updated
                self._registry[target_id]['modified'] = datetime.now().isoformat()
                self._save_registry()
    
    def rename_entity(self, old_name: str, new_name: str):
        """
        Rename an entity while preserving its ID and relationships.
        
        Args:
            old_name: Current entity name
            new_name: New entity name
            
        Raises:
            ValueError: If old entity doesn't exist or new name conflicts
        """
        # Find the entity
        entity_id = self._resolve_to_id(old_name)
        if not entity_id:
            raise ValueError(f"Entity '{old_name}' not found")
            
        # Check if new name is already taken
        existing_id = self._resolve_to_id(new_name)
        if existing_id and existing_id != entity_id:
            raise ValueError(f"Entity with name '{new_name}' already exists")
            
        # Update the entity's original_id
        self._registry[entity_id]['original_id'] = new_name
        self._registry[entity_id]['modified'] = datetime.now().isoformat()
        
        # Update name index
        # Remove old name from index
        if old_name in self._name_index:
            del self._name_index[old_name]
        # Add new name to index
        self._name_index[new_name] = entity_id
        
        # Save changes
        self._save_registry()
    
    def find_by_content_hash(self, content_hash: str) -> List[EntityInfo]:
        """
        Find entities with matching content hash.
        
        Args:
            content_hash: SHA256 hash of file content
            
        Returns:
            List of EntityInfo objects with matching content
        """
        matches = []
        
        for entity_id, entity in self._registry.items():
            for format_type, format_data in entity['formats'].items():
                # Check if this format has a content hash
                if format_data.get('metadata', {}).get('content_hash') == content_hash:
                    matches.append(EntityInfo(
                        entity_id=entity_id,
                        original_id=entity['original_id'],
                        format_type=format_type,
                        file_path=format_data['file_path'],
                        metadata=format_data.get('metadata', {}),
                        created=format_data.get('created', entity['created']),
                        modified=entity['modified']
                    ))
        
        return matches