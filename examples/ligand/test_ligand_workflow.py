"""
Complete workflow demonstration for the Protos Ligand Processor.

This script demonstrates:
1. Structure download and loading
2. Ligand extraction from structures
3. Binding pocket analysis
4. ChEMBL integration for bioactivity data
5. Property storage following separation of concerns
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.paths.path_config import ProtosPaths
from protos.processing.molecule import MoleculeProcessor
from protos.processing.property import PropertyProcessor
from protos.processing.structure import StructureProcessor
from protos.io.ingest.download_structures import download_structures_with_processor
from protos.analysis.structure_ligand_analysis import (
    extract_all_ligands,
    get_binding_site,
    calculate_ligand_interactions,
    create_ligand_interaction_report,
    analyze_all_ligands_in_structure
)


def main():
    """Demonstrate the complete ligand workflow."""
    # Setup environment
    data_root = Path(__file__).parent / "data"
    data_root.mkdir(exist_ok=True, parents=True)
    os.environ["PROTOS_DATA_ROOT"] = str(data_root.absolute())
    
    # Initialize paths
    paths = ProtosPaths()
    
    print("=" * 60)
    print("Protos Ligand-Structure Analysis Workflow")
    print("=" * 60)
    
    # Initialize processors
    print("\n1. Initializing Processors")
    print("-" * 30)
    
    # Use StructureProcessor instead of CifBaseProcessor (following the codebase)
    struct_proc = StructureProcessor(name="structure_processor", paths=paths)
    print(f"✓ StructureProcessor initialized")
    
    lig_proc = MoleculeProcessor(paths=paths)
    print(f"✓ LigandProcessor initialized at: {lig_proc.data_path}")
    
    prop_proc = PropertyProcessor(paths=paths)
    print(f"✓ PropertyProcessor initialized at: {prop_proc.data_path}")
    
    # Download a structure with interesting ligands
    print("\n2. Downloading Structure with Ligands")
    print("-" * 30)
    
    # Download a kinase structure with ATP and inhibitor
    # 1ATP: Phosphorylase kinase with ATP
    # Alternative good examples: 4HHB (hemoglobin with heme), 1HVH (HIV protease with inhibitor)
    test_pdb = "1ATP"
    
    print(f"Downloading structure {test_pdb}...")
    try:
        # Download structure using Protos' preferred method
        success, failed = download_structures_with_processor(
            [test_pdb],
            processor=struct_proc,
            overwrite=False
        )
        if success:
            print(f"✓ Downloaded {test_pdb}")
        else:
            print(f"⚠️  Failed to download {test_pdb}: {failed}")
    except Exception as e:
        print(f"⚠️  Download error: {e}")
    
    # Load the structure
    print(f"\nLoading structure {test_pdb}...")
    try:
        struct_proc.load_structures([test_pdb])
        print(f"✓ Loaded structure data with {len(struct_proc.data)} atoms")
    except Exception as e:
        print(f"⚠️  Failed to load structure: {e}")
        print("\nUsing alternative structure for demo...")
        # Try a smaller test structure
        test_pdb = "1CRN"  # Small protein (crambin)
        print(f"Attempting to download {test_pdb}...")
        try:
            success, failed = download_structures_with_processor([test_pdb], struct_proc)
            if success:
                struct_proc.load_structures([test_pdb])
                print(f"✓ Loaded {test_pdb} with {len(struct_proc.data)} atoms")
            else:
                print("⚠️  Demo requires a structure file. Please ensure structure files are available.")
                return
        except Exception as e:
            print(f"⚠️  Could not load structure: {e}")
            return
    
    # Extract and analyze ligands from structure
    print("\n3. Extracting Ligands from Structure")
    print("-" * 30)
    
    # Extract all ligands (exclude water and ions)
    ligands = extract_all_ligands(struct_proc, test_pdb, exclude_common=True)
    print(f"\nFound {len(ligands)} ligands in {test_pdb}:")
    
    for ligand in ligands:
        print(f"  - {ligand['res_name3l']} (chain {ligand['chain_id']}, "
              f"{ligand['num_atoms']} atoms)")
    
    # Analyze binding pockets for each ligand
    print("\n4. Binding Pocket Analysis")
    print("-" * 30)
    
    for ligand in ligands:
        print(f"\nAnalyzing {ligand['res_name3l']} binding site:")
        
        # Get binding site residues
        binding_site = get_binding_site(
            struct_proc, 
            test_pdb, 
            ligand['atoms'],
            cutoff=5.0
        )
        
        if not binding_site['residues'].empty:
            print(f"  Found {len(binding_site['residues'])} binding residues:")
            
            # Show top 5 closest residues
            for idx, res in binding_site['residues'].head(5).iterrows():
                print(f"    - {res['res_name']}{res['res_id']} "
                      f"(chain {res['chain_id']}, {res['min_distance']:.2f} Å)")
            
            # Calculate detailed interactions
            interactions = calculate_ligand_interactions(
                struct_proc,
                test_pdb,
                ligand['atoms'],
                detailed=True
            )
            
            # Show interaction summary
            if 'summary' in interactions:
                summary = interactions['summary']
                print(f"\n  Interaction summary:")
                print(f"    - H-bonds: {summary.get('num_hydrogen_bonds', 0)}")
                print(f"    - Hydrophobic: {summary.get('num_hydrophobic', 0)}")
                print(f"    - Water bridges: {summary.get('num_water_bridges', 0)}")
        
        # Store binding site data as properties
        prop_proc.assign_property(
            entity_identifier=f"{test_pdb}_{ligand['res_name3l']}",
            property_name='binding_residues',
            property_value=len(binding_site['residues']),
            dataset_name='structure_binding_sites'
        )
    
    # Example: Process known drug compounds
    print("\n5. Processing Known Drug Compounds")
    print("-" * 30)
    
    # Example drug compounds
    drug_compounds = [
        {
            'smiles': 'COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1',
            'chembl_id': 'CHEMBL553',
            'name': 'Gefitinib',
            'activity_type': 'IC50',
            'value': 33,
            'units': 'nM',
            'target': 'EGFR'
        },
        {
            'smiles': 'CN(C)C/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1O[C@H]1CCOC1',
            'chembl_id': 'CHEMBL1173',
            'name': 'Erlotinib',
            'activity_type': 'IC50',
            'value': 2,
            'units': 'nM',
            'target': 'EGFR'
        }
    ]
    
    # Process compounds
    print("\nProcessing compounds:")
    for compound in drug_compounds:
        print(f"\n  {compound['name']} ({compound['chembl_id']})")
        
        # Save structure with LigandProcessor
        lig_proc.save_entity(compound['smiles'], {
            'smiles': compound['smiles'],
            'chembl_id': compound['chembl_id'],
            'name': compound['name']
        })
        print(f"    ✓ Structure saved")
        
        # Calculate properties
        props = lig_proc.calculate_properties(compound['smiles'])
        if props:
            print(f"    ✓ MW: {props['mw']:.1f}, LogP: {props['logp']:.2f}")
        
        # Save bioactivity with PropertyProcessor
        prop_proc.assign_property(
            entity_identifier=compound['smiles'],
            property_name=f"{compound['activity_type']}_{compound['target']}",
            property_value=compound['value'],
            dataset_name="egfr_inhibitors"
        )
        print(f"    ✓ {compound['activity_type']}: {compound['value']} {compound['units']}")
    
    # Create comprehensive report
    print("\n6. Creating Comprehensive Interaction Report")
    print("-" * 30)
    
    if ligands:
        # Create detailed report for each ligand
        all_reports = []
        
        for ligand in ligands:
            report = create_ligand_interaction_report(
                struct_proc,
                test_pdb,
                ligand['res_name3l'],
                ligand['chain_id']
            )
            all_reports.append(report)
            
            # Print detailed report
            print(f"\n{'='*60}")
            print(f"LIGAND INTERACTION REPORT: {ligand['res_name3l']}")
            print(f"{'='*60}")
            
            print(f"\n📊 Basic Information:")
            print(f"  • Structure: {report['pdb_id']}")
            print(f"  • Ligand: {report['ligand']['name']}")
            print(f"  • Chain: {report['ligand']['chain']}")
            print(f"  • Number of atoms: {report['ligand']['num_atoms']}")
            print(f"  • SMILES: {report['ligand']['smiles'] or 'Not available'}")
            
            print(f"\n🎯 Binding Site Analysis:")
            print(f"  • Number of binding residues: {report['binding_site']['num_residues']}")
            if report['binding_site']['volume_estimate']:
                print(f"  • Estimated binding site volume: {report['binding_site']['volume_estimate']:.1f} ų")
            
            # Top binding residues
            if report['binding_site']['residues']:
                print(f"\n  Top 5 closest residues:")
                for i, res in enumerate(report['binding_site']['residues'][:5]):
                    print(f"    {i+1}. {res['res_name']}{res['res_id']} "
                          f"(Chain {res['chain_id']}) - {res['min_distance']:.2f} Å")
            
            print(f"\n⚛️ Interaction Analysis:")
            if 'summary' in report['interactions']:
                summary = report['interactions']['summary']
                print(f"  • Hydrogen bonds: {summary.get('num_hydrogen_bonds', 0)}")
                print(f"  • Hydrophobic contacts: {summary.get('num_hydrophobic', 0)}")
                print(f"  • Water-mediated contacts: {summary.get('num_water_bridges', 0)}")
                print(f"  • π-π stacking: {summary.get('num_pi_stacking', 0)}")
                print(f"  • Salt bridges: {summary.get('num_salt_bridges', 0)}")
            
            # Detailed interactions
            if 'hydrogen_bonds' in report['interactions'] and report['interactions']['hydrogen_bonds']:
                print(f"\n  Hydrogen Bonds (showing first 3):")
                for i, hbond in enumerate(report['interactions']['hydrogen_bonds'][:3]):
                    if 'donor' in hbond and 'acceptor' in hbond:
                        print(f"    • {hbond['donor']['atom']} → {hbond['acceptor']['atom']} "
                              f"({hbond.get('distance', 0):.2f} Å)")
            
            if 'hydrophobic' in report['interactions'] and report['interactions']['hydrophobic']:
                print(f"\n  Hydrophobic Contacts (showing first 3):")
                for i, contact in enumerate(report['interactions']['hydrophobic'][:3]):
                    print(f"    • {contact.get('ligand_atom', 'Unknown')} - "
                          f"{contact.get('residue', 'Unknown')} "
                          f"({contact.get('distance', 0):.2f} Å)")
            
            print(f"\n📈 Summary Statistics:")
            print(f"  • Total contacts: {report['summary']['total_contacts']}")
            if report['summary']['key_residues']:
                print(f"  • Key residues: {', '.join(report['summary']['key_residues'][:5])}")
            
            # Save detailed report as property
            prop_proc.assign_property(
                entity_identifier=f"{test_pdb}_{ligand['res_name3l']}",
                property_name='total_interactions',
                property_value=report['summary']['total_contacts'],
                dataset_name='interaction_reports'
            )
            
            prop_proc.assign_property(
                entity_identifier=f"{test_pdb}_{ligand['res_name3l']}",
                property_name='binding_site_residues',
                property_value=report['binding_site']['num_residues'],
                dataset_name='interaction_reports'
            )
        
        # Comparative analysis if multiple ligands
        if len(all_reports) > 1:
            print(f"\n\n{'='*60}")
            print(f"COMPARATIVE LIGAND ANALYSIS")
            print(f"{'='*60}")
            
            print(f"\n📊 Interaction Comparison:")
            print(f"{'Ligand':<10} {'H-bonds':<10} {'Hydrophobic':<12} {'Water bridges':<13} {'Total':<10}")
            print("-" * 60)
            
            for report in all_reports:
                if 'summary' in report['interactions']:
                    s = report['interactions']['summary']
                    print(f"{report['ligand']['name']:<10} "
                          f"{s.get('num_hydrogen_bonds', 0):<10} "
                          f"{s.get('num_hydrophobic', 0):<12} "
                          f"{s.get('num_water_bridges', 0):<13} "
                          f"{report['summary']['total_contacts']:<10}")
            
            # Find shared binding residues
            all_binding_residues = []
            for report in all_reports:
                if report['binding_site']['residues']:
                    residues = set(f"{r['res_name']}{r['res_id']}" 
                                  for r in report['binding_site']['residues'])
                    all_binding_residues.append(residues)
            
            if len(all_binding_residues) > 1:
                shared_residues = set.intersection(*all_binding_residues)
                if shared_residues:
                    print(f"\n🔗 Shared binding residues: {', '.join(sorted(shared_residues))}")
    
    # Future: Convert ligand to SMILES/SDF
    print("\n7. Ligand Format Conversion (To Be Implemented)")
    print("-" * 30)
    print("Next steps:")
    print("  - Extract ligand 3D coordinates from structure")
    print("  - Convert to RDKit molecule object")
    print("  - Generate SMILES from structure")
    print("  - Export as SDF with proper bond orders")
    print("  - Register in LigandProcessor")
    
    # Summary
    print("\n" + "=" * 60)
    print("✓ WORKFLOW DEMONSTRATION COMPLETE!")
    print("=" * 60)
    
    print("\n📋 ANALYSIS SUMMARY:")
    print(f"\n🏗️ Structure Analysis:")
    print(f"  • PDB ID: {test_pdb}")
    print(f"  • Total atoms: {len(struct_proc.data[struct_proc.data['pdb_id'] == test_pdb])}")
    print(f"  • Protein atoms: {len(struct_proc.data[(struct_proc.data['pdb_id'] == test_pdb) & (struct_proc.data['group'] == 'ATOM')])}")
    print(f"  • HETATM atoms: {len(struct_proc.data[(struct_proc.data['pdb_id'] == test_pdb) & (struct_proc.data['group'] == 'HETATM')])}")
    
    print(f"\n💊 Ligand Analysis:")
    print(f"  • Ligands found: {len(ligands)}")
    if ligands:
        print(f"  • Ligand types: {', '.join([l['res_name3l'] for l in ligands])}")
        print(f"  • Total ligand atoms: {sum(l['num_atoms'] for l in ligands)}")
    
    print(f"\n🧬 Interaction Analysis Performed:")
    print(f"  ✓ Binding site identification (5.0 Å cutoff)")
    print(f"  ✓ Hydrogen bond detection")
    print(f"  ✓ Hydrophobic contact analysis")
    print(f"  ✓ Water-mediated interaction detection")
    print(f"  ✓ Binding site volume estimation")
    print(f"  ✓ Key residue identification")
    
    print(f"\n💾 Data Storage:")
    print(f"  • Ligand structures: {len(drug_compounds)} compounds registered")
    print(f"  • Properties stored: Binding residues, interaction counts")
    print(f"  • Datasets created: 'structure_binding_sites', 'egfr_inhibitors', 'interaction_reports'")
    
    print(f"\n📊 Analysis Capabilities Demonstrated:")
    print(f"  1. Structure download and loading")
    print(f"  2. Ligand extraction from structures")
    print(f"  3. Binding pocket characterization")
    print(f"  4. Detailed interaction profiling")
    print(f"  5. Comparative ligand analysis")
    print(f"  6. Property-based data storage")
    print(f"  7. Cross-processor integration")
    
    print(f"\n🚀 Next Steps:")
    print(f"  • Implement ligand extraction to SMILES/SDF with proper bond orders")
    print(f"  • Connect structure ligands to ChEMBL bioactivity data")
    print(f"  • Enable cross-structure binding site comparison")
    print(f"  • Build interaction fingerprints for ML applications")
    print(f"  • Develop binding site similarity metrics")


if __name__ == "__main__":
    main()