"""
Conflict resolution for entity registration.

This module provides strategies for handling conflicts when registering
new entities that may duplicate existing content or names.
"""

import hashlib
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any, List, Tuple
from enum import Enum
from dataclasses import dataclass
import logging

from protos.io.core.entity_registry import EntityRegistry, EntityInfo


class ConflictResolutionStrategy(Enum):
    """Strategies for handling naming conflicts during registration."""
    SKIP = "skip"              # Skip registration, return existing
    VERSION = "version"        # Create versioned name (e.g., protein_v2)
    REPLACE = "replace"        # Replace existing (with backup)
    MERGE = "merge"            # Merge with existing (for tables)
    ASK = "ask"                # Prompt user for decision


class ConflictAction(Enum):
    """Actions taken to resolve conflicts."""
    SKIP = "skip"
    VERSION = "version"
    REPLACE = "replace"
    MERGE = "merge"
    ADD_ALIAS = "add_alias"
    ERROR = "error"


@dataclass
class ResolutionResult:
    """Result of conflict resolution."""
    action: ConflictAction
    entity_name: Optional[str] = None  # Final entity name to use
    entity: Optional[EntityInfo] = None  # Existing entity if applicable
    backup_path: Optional[Path] = None  # Backup location if replaced
    message: Optional[str] = None  # Explanation message


class ConflictResolver:
    """Handles conflicts during entity registration."""
    
    def __init__(self, entity_registry: EntityRegistry, backup_dir: Optional[Path] = None):
        """
        Initialize conflict resolver.
        
        Args:
            entity_registry: The entity registry to check against
            backup_dir: Directory for backups when replacing
        """
        self.registry = entity_registry
        self.backup_dir = backup_dir
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def resolve_duplicate_content(self, 
                                existing_entities: List[EntityInfo],
                                new_name: str,
                                format_type: str) -> ResolutionResult:
        """
        Handle case where file content already exists.
        
        Args:
            existing_entities: List of entities with same content
            new_name: Proposed name for new entity
            format_type: Format type being registered
            
        Returns:
            ResolutionResult with action to take
        """
        # Check if any existing entity has the same name
        for entity in existing_entities:
            if entity.original_id == new_name:
                # Exact duplicate - safe to skip
                return ResolutionResult(
                    action=ConflictAction.SKIP,
                    entity=entity,
                    message=f"Identical file already registered as '{new_name}'"
                )
        
        # Same content, different name
        # Could add as alias, but for now just skip with info
        primary_entity = existing_entities[0]
        return ResolutionResult(
            action=ConflictAction.SKIP,
            entity=primary_entity,
            message=f"Content already registered as '{primary_entity.original_id}'"
        )
    
    def resolve_name_conflict(self,
                            existing_name: str,
                            new_file: Path,
                            format_type: str,
                            strategy: ConflictResolutionStrategy = ConflictResolutionStrategy.SKIP) -> ResolutionResult:
        """
        Handle case where name already exists with different content.
        
        Args:
            existing_name: The conflicting entity name
            new_file: Path to new file being registered
            format_type: Format type
            strategy: Conflict resolution strategy
            
        Returns:
            ResolutionResult with action to take
        """
        # Get existing entity info
        existing_entity = self.registry.find_entity(existing_name, format_type)
        
        if strategy == ConflictResolutionStrategy.SKIP:
            return ResolutionResult(
                action=ConflictAction.SKIP,
                entity=existing_entity,
                message=f"Entity '{existing_name}' already exists"
            )
        
        elif strategy == ConflictResolutionStrategy.VERSION:
            # Generate versioned name
            versioned_name = self._generate_version_name(existing_name, format_type)
            return ResolutionResult(
                action=ConflictAction.VERSION,
                entity_name=versioned_name,
                message=f"Creating versioned entity '{versioned_name}'"
            )
        
        elif strategy == ConflictResolutionStrategy.REPLACE:
            # Create backup of existing
            if existing_entity and self.backup_dir:
                backup_path = self._backup_existing(existing_entity)
                return ResolutionResult(
                    action=ConflictAction.REPLACE,
                    entity_name=existing_name,
                    backup_path=backup_path,
                    message=f"Replacing '{existing_name}' (backup at {backup_path})"
                )
            else:
                return ResolutionResult(
                    action=ConflictAction.REPLACE,
                    entity_name=existing_name,
                    message=f"Replacing '{existing_name}' (no backup created)"
                )
        
        else:
            # Unsupported strategy for file conflicts
            return ResolutionResult(
                action=ConflictAction.ERROR,
                message=f"Strategy '{strategy}' not supported for file conflicts"
            )
    
    def resolve_table_conflict(self,
                             existing_table: str,
                             new_table_path: Path,
                             strategy: ConflictResolutionStrategy = ConflictResolutionStrategy.MERGE) -> ResolutionResult:
        """
        Handle conflicts for table-based data (property, GRN).
        
        Args:
            existing_table: Name of existing table
            new_table_path: Path to new table file
            strategy: Conflict resolution strategy
            
        Returns:
            ResolutionResult with action to take
        """
        if strategy == ConflictResolutionStrategy.MERGE:
            # Tables can be merged
            return ResolutionResult(
                action=ConflictAction.MERGE,
                entity_name=existing_table,
                message=f"Will merge with existing table '{existing_table}'"
            )
        elif strategy == ConflictResolutionStrategy.VERSION:
            versioned_name = self._generate_version_name(existing_table, 'table')
            return ResolutionResult(
                action=ConflictAction.VERSION,
                entity_name=versioned_name,
                message=f"Creating new table '{versioned_name}'"
            )
        else:
            # Fall back to standard name conflict resolution
            return self.resolve_name_conflict(
                existing_table, new_table_path, 'table', strategy
            )
    
    def _generate_version_name(self, base_name: str, format_type: str) -> str:
        """
        Generate a versioned name that doesn't conflict.
        
        Args:
            base_name: Original name
            format_type: Format type to check
            
        Returns:
            Available versioned name
        """
        version = 2
        while True:
            versioned_name = f"{base_name}_v{version}"
            if not self.registry.entity_exists(versioned_name, format_type):
                return versioned_name
            version += 1
    
    def _backup_existing(self, entity_info: EntityInfo) -> Path:
        """
        Create backup of existing entity file.
        
        Args:
            entity_info: Entity to backup
            
        Returns:
            Path to backup file
        """
        if not self.backup_dir:
            raise ValueError("No backup directory configured")
        
        self.backup_dir.mkdir(parents=True, exist_ok=True)
        
        # Create timestamped backup name
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        source_path = Path(entity_info.file_path)
        backup_name = f"{timestamp}_{entity_info.original_id}{source_path.suffix}"
        backup_path = self.backup_dir / backup_name
        
        # Copy file
        if source_path.exists():
            import shutil
            shutil.copy2(source_path, backup_path)
            self.logger.info(f"Created backup: {backup_path}")
        
        return backup_path
    
    def check_conflicts(self,
                       file_path: Path,
                       proposed_name: str,
                       format_type: str,
                       content_hash: Optional[str] = None) -> Tuple[bool, Optional[ResolutionResult]]:
        """
        Check for all types of conflicts before registration.
        
        Args:
            file_path: Path to file being registered
            proposed_name: Proposed entity name
            format_type: Format type
            content_hash: Optional pre-computed content hash
            
        Returns:
            Tuple of (has_conflict, resolution_result)
        """
        # Check content duplicate if hash provided
        if content_hash:
            existing_by_content = self.registry.find_by_content_hash(content_hash)
            if existing_by_content:
                return True, self.resolve_duplicate_content(
                    existing_by_content, proposed_name, format_type
                )
        
        # Check name conflict
        if self.registry.entity_exists(proposed_name, format_type):
            return True, self.resolve_name_conflict(
                proposed_name, file_path, format_type,
                ConflictResolutionStrategy.SKIP
            )
        
        # No conflicts
        return False, None