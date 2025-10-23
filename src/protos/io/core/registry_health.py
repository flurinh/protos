"""
Registry health checks and dataset integrity validation.

This module provides tools to ensure the entity registry and datasets
remain consistent and all referenced files are accessible.
"""

import os
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional, Set
from dataclasses import dataclass, field

from protos.io.paths import ProtosPaths
from protos.io.core.entity_registry import EntityRegistry
from protos.io.core.dataset_manager import DatasetManager


@dataclass
class DatasetIntegrityReport:
    """Report on dataset integrity check results."""
    dataset_name: str
    total_entities: int = 0
    valid_entities: List[str] = field(default_factory=list)
    missing_registrations: List[str] = field(default_factory=list)
    missing_files: List[Dict[str, str]] = field(default_factory=list)
    inaccessible_files: List[Dict[str, str]] = field(default_factory=list)
    
    @property
    def is_healthy(self) -> bool:
        """Check if dataset is completely healthy."""
        return (len(self.missing_registrations) == 0 and
                len(self.missing_files) == 0 and
                len(self.inaccessible_files) == 0)
    
    @property
    def health_percentage(self) -> float:
        """Calculate health percentage."""
        if self.total_entities == 0:
            return 100.0
        return (len(self.valid_entities) / self.total_entities) * 100
    
    def add_valid(self, entity_name: str):
        """Add a valid entity."""
        self.valid_entities.append(entity_name)
    
    def add_missing_registration(self, entity_name: str):
        """Add entity with missing registration."""
        self.missing_registrations.append(entity_name)
    
    def add_missing_file(self, entity_name: str, file_path: str):
        """Add entity with missing file."""
        self.missing_files.append({
            'entity': entity_name,
            'path': file_path
        })
    
    def add_inaccessible(self, entity_name: str, file_path: str):
        """Add entity with inaccessible file."""
        self.inaccessible_files.append({
            'entity': entity_name,
            'path': file_path
        })
    
    def display_summary(self):
        """Display a summary of the integrity check."""
        print(f"\nDataset: {self.dataset_name}")
        print(f"Total entities: {self.total_entities}")
        print(f"Valid entities: {len(self.valid_entities)} ({self.health_percentage:.1f}%)")
        
        if self.missing_registrations:
            print(f"Missing registrations: {len(self.missing_registrations)}")
            for name in self.missing_registrations[:5]:
                print(f"  - {name}")
            if len(self.missing_registrations) > 5:
                print(f"  ... and {len(self.missing_registrations) - 5} more")
        
        if self.missing_files:
            print(f"Missing files: {len(self.missing_files)}")
            for item in self.missing_files[:5]:
                print(f"  - {item['entity']}: {item['path']}")
            if len(self.missing_files) > 5:
                print(f"  ... and {len(self.missing_files) - 5} more")
        
        if self.inaccessible_files:
            print(f"Inaccessible files: {len(self.inaccessible_files)}")
            for item in self.inaccessible_files[:5]:
                print(f"  - {item['entity']}: {item['path']}")
            if len(self.inaccessible_files) > 5:
                print(f"  ... and {len(self.inaccessible_files) - 5} more")
        
        if self.is_healthy:
            print("✓ Dataset is healthy")
        else:
            print("✗ Dataset has issues")


class RegistryHealthChecker:
    """Performs health checks on entity registry and datasets."""
    
    def __init__(self, entity_registry: EntityRegistry, paths: ProtosPaths):
        """
        Initialize health checker.
        
        Args:
            entity_registry: The entity registry to check
            paths: ProtosPaths instance for data directories
        """
        self.registry = entity_registry
        self.paths = paths
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def find_unregistered_files(self, processor_type: str) -> List[Path]:
        """
        Find files in data directories not in registry.
        
        Args:
            processor_type: Type of processor (e.g., 'structure', 'sequence')
            
        Returns:
            List of unregistered file paths
        """
        unregistered = []
        
        # Get processor directory
        processor_path = Path(self.paths.get_processor_path(processor_type))
        if not processor_path.exists():
            return unregistered
        
        # Define subdirectories to check based on processor type
        subdirs_to_check = self._get_data_subdirs(processor_type)
        
        # Scan each subdirectory
        for subdir in subdirs_to_check:
            data_dir = processor_path / subdir
            if data_dir.exists():
                # Get appropriate file patterns
                patterns = self._get_file_patterns(processor_type)
                
                for pattern in patterns:
                    for file_path in data_dir.glob(pattern):
                        if file_path.is_file():
                            # Check if file is registered
                            relative_path = file_path.relative_to(Path(self.paths.data_root))
                            
                            # Search for this file in registry
                            found = False
                            for entity_data in self.registry._registry.values():
                                for format_data in entity_data['formats'].values():
                                    if format_data['file_path'] == str(relative_path):
                                        found = True
                                        break
                                if found:
                                    break
                            
                            if not found:
                                unregistered.append(file_path)
        
        return unregistered
    
    def check_dataset_integrity(self, dataset_name: str, 
                              processor_type: str) -> DatasetIntegrityReport:
        """
        Check if all dataset entities exist and are accessible.
        
        Args:
            dataset_name: Name of dataset to check
            processor_type: Processor type
            
        Returns:
            DatasetIntegrityReport with results
        """
        report = DatasetIntegrityReport(dataset_name)
        
        # Load dataset
        dataset_manager = DatasetManager(processor_type=processor_type, 
                                       paths=self.paths,
                                       entity_registry=self.registry)
        
        if not dataset_manager.dataset_exists(dataset_name):
            self.logger.warning(f"Dataset '{dataset_name}' not found")
            return report
        
        # Get dataset entities
        try:
            entity_names = dataset_manager.get_dataset_entities(dataset_name)
            report.total_entities = len(entity_names)
        except Exception as e:
            self.logger.error(f"Error loading dataset: {e}")
            return report
        
        # Check each entity
        for entity_name in entity_names:
            # Check if entity is registered
            if not self.registry.entity_exists(entity_name, processor_type):
                report.add_missing_registration(entity_name)
                continue
            
            # Get entity info
            entity_info = self.registry.find_entity(entity_name, processor_type)
            if not entity_info:
                report.add_missing_registration(entity_name)
                continue
            
            # Check if file exists
            file_path = Path(entity_info.file_path)
            if not file_path.is_absolute():
                file_path = Path(self.paths.data_root) / file_path
            
            if not file_path.exists():
                report.add_missing_file(entity_name, str(file_path))
            elif not file_path.is_file():
                report.add_inaccessible(entity_name, str(file_path) + " (not a file)")
            elif not os.access(file_path, os.R_OK):
                report.add_inaccessible(entity_name, str(file_path) + " (not readable)")
            else:
                report.add_valid(entity_name)
        
        return report
    
    def check_orphaned_files(self) -> Dict[str, List[Path]]:
        """
        Find files that are registered but don't exist.
        
        Returns:
            Dict mapping entity names to missing file paths
        """
        orphaned = {}
        
        for entity_id, entity_data in self.registry._registry.items():
            entity_name = entity_data['original_id']
            
            for format_type, format_data in entity_data['formats'].items():
                file_path = Path(format_data['file_path'])
                if not file_path.is_absolute():
                    file_path = Path(self.paths.data_root) / file_path
                
                if not file_path.exists():
                    if entity_name not in orphaned:
                        orphaned[entity_name] = []
                    orphaned[entity_name].append(file_path)
        
        return orphaned
    
    def check_registry_consistency(self) -> Dict[str, Any]:
        """
        Check internal consistency of registry.
        
        Returns:
            Dict with consistency check results
        """
        results = {
            'total_entities': len(self.registry._registry),
            'total_names': len(self.registry._name_index),
            'issues': []
        }
        
        # Check name index consistency
        for name, entity_id in self.registry._name_index.items():
            if entity_id not in self.registry._registry:
                results['issues'].append({
                    'type': 'orphaned_name_index',
                    'name': name,
                    'entity_id': entity_id
                })
        
        # Check entity names match index
        for entity_id, entity_data in self.registry._registry.items():
            name = entity_data['original_id']
            if name not in self.registry._name_index:
                results['issues'].append({
                    'type': 'missing_name_index',
                    'name': name,
                    'entity_id': entity_id
                })
            elif self.registry._name_index[name] != entity_id:
                results['issues'].append({
                    'type': 'name_index_mismatch',
                    'name': name,
                    'entity_id': entity_id,
                    'indexed_id': self.registry._name_index[name]
                })
        
        # Check relationships reference valid entities
        for entity_id, entity_data in self.registry._registry.items():
            for rel in entity_data.get('relationships', []):
                if rel['source'] not in self.registry._registry:
                    results['issues'].append({
                        'type': 'invalid_relationship_source',
                        'entity': entity_data['original_id'],
                        'rel_source': rel['source']
                    })
                if rel['target'] not in self.registry._registry:
                    results['issues'].append({
                        'type': 'invalid_relationship_target',
                        'entity': entity_data['original_id'],
                        'rel_target': rel['target']
                    })
        
        results['is_consistent'] = len(results['issues']) == 0
        return results
    
    def generate_health_report(self) -> Dict[str, Any]:
        """
        Generate comprehensive health report.
        
        Returns:
            Dict with full health check results
        """
        report = {
            'timestamp': datetime.now().isoformat(),
            'registry_size': len(self.registry._registry),
            'consistency': self.check_registry_consistency(),
            'orphaned_files': self.check_orphaned_files(),
            'processor_health': {}
        }
        
        # Check each processor type
        for processor_type in ['structure', 'sequence', 'grn', 'property', 'embedding']:
            try:
                unregistered = self.find_unregistered_files(processor_type)
                report['processor_health'][processor_type] = {
                    'unregistered_files': len(unregistered),
                    'examples': [str(f) for f in unregistered[:5]]
                }
            except Exception as e:
                report['processor_health'][processor_type] = {
                    'error': str(e)
                }
        
        return report
    
    def _get_data_subdirs(self, processor_type: str) -> List[str]:
        """Get data subdirectories to check for processor type."""
        subdirs = {
            'structure': ['mmcif', 'cache', 'pdb'],
            'sequence': ['fasta', 'alignments'],
            'grn': ['tables', 'assignments', 'ref'],
            'property': ['tables', 'cache'],
            'embedding': ['embeddings', 'models']
        }
        return subdirs.get(processor_type, [])
    
    def _get_file_patterns(self, processor_type: str) -> List[str]:
        """Get file patterns to search for processor type."""
        from protos.io.format_registry import format_registry, ProcessorType
        
        # Map string processor type to enum
        type_map = {
            'structure': ProcessorType.STRUCTURE,
            'sequence': ProcessorType.SEQUENCE,
            'grn': ProcessorType.GRN,
            'property': ProcessorType.PROPERTY,
            'embedding': ProcessorType.EMBEDDING,
            'ligand': ProcessorType.MOLECULE,
            'graph': ProcessorType.GRAPH
        }
        
        proc_type = type_map.get(processor_type)
        if proc_type:
            extensions = format_registry.get_extensions_for_processor(proc_type)
            return [f'*{ext}' for ext in extensions]
        else:
            return ['*']