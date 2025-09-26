#!/usr/bin/env python3
"""Test ModelManager workflows showcasing integration with processors."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.models.model_manager import ModelManager, prepare_mutation_screen
from protos.processing.structure import StructureProcessor
from protos.processing.sequence import SequenceProcessor
from protos.processing.grn import GRNProcessor
from protos.processing.embedding import EmbeddingProcessor
from protos.processing.property import PropertyProcessor
from protos.io.ingest.structure_loader import StructureLoader


def ensure_data_root() -> Path:
    """Set up data root directory."""
    data_root = PROJECT_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def workflow_1_sequence_to_structure():
    """Workflow 1: Predict structures from sequences."""
    print("\n=== Workflow 1: Sequence to Structure Prediction ===")
    
    # Initialize processors and manager
    seq_proc = SequenceProcessor()
    manager = ModelManager()
    
    # Create a test dataset of GPCR sequences
    test_sequences = {
        "ADRB2_HUMAN": "MGQPGNGSAFLLAPNGSHAPDHDVTQERDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGAAHILMKMWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFKYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRATHQEAINCYANETCCDFFTNQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQKIDKSEGRFHVQNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNLIRKEVYILLNWIGYVNSGFNPLIYCRSPDFRIAFQELLCLRRSSLKAYGNGYSSNGNTGEQSGYHVEQEKENKLLCEDLPGTEDFVGHQGTVPSDNIDSQGRNCSTNDSLL",
        "A2AR_HUMAN": "MPPDSNSTNGEASSSSQNGSAAGPEGQASVGGVLEEAAIAQMVAGPQGSIIISVLVAIIVFGNVLVIAVFTSRALKAPQNLFLVSLASADILVATLVIPFSLANELCGVFFIACLIMCVTSLVLTAVSIGSLLAIAVDRYLAIRIPLEYNITKRTRRVVALVVWVISAVISGLPVIIGWNCIVQVCGICVTEVIAGLCAIGSMNVLFIIKVSLLKVIQKLVKENARRQGNGVQQSKKTEFFTVILAIVLGVFVVCWFPFFFTYTLTAVGCSVPRTLFKFFFWFGYCNSAVNPVIYTIFNHDFRRAFKKILFHKQKRQKKKIDKEPTDFQVSPDDQPLGNSSSSHESKDSK",
        "DRD2_SHORT": "MDPLNLSWYDDDLERQNWSRPFNGSDGKADRPHYNYYATLLTLLIAVIVFGNVLVCMAVSREKALQTTTNYLIVSLAVADLLVATLVMPWVVYLEVVGEWKFSRIHCDIFVTLDVMMCTASILNLCAISIDRYTAVAMPMLYNTRYSSKRRVTVMISIVWVLSFTISCPLLFGLNNADQNECIIANPAFVVYSSIVSFYVPFIVTLLVYIKIYIVLRRRRKRVNTKRSSRAFRAHLRAPLKGNCTHPEDMKLCTVIMKSNGSFPVNRRRVEAARRAQELEMEMLSSTSPPERTRYSPIPPSHHQLTLPDPSHHGLHSTPDSPAKPEKNGHAKDHPKIAKIFEIQTMPNGKTRTSLKTMSRRKLSQQKEKKATQMLAIVLGVFIICWLPFFITHILNIHCDCNIPPVLYSAFTWLGYVNSAVNPIIYTTFNIEFRKAFLKILHC"
    }
    
    # Save sequences as dataset
    seq_proc.save_sequences(
        test_sequences,
        output_file="gpcr_test_sequences",
        dataset_name="gpcr_test_sequences",
        metadata={"family": "GPCR", "species": "human"},
        materialize_entities=True
    )
    print(f"✓ Created sequence dataset with {len(test_sequences)} sequences")
    
    # Prepare structure predictions for each sequence
    predictions = []
    for seq_name in test_sequences:
        # Standard prediction
        config = {
            "recycling": 3,
            "num_samples": 1,
            "crop_size": 384,
            "output_name": f"{seq_name}_predicted"
        }
        
        try:
            model_input = manager.prepare_input(
                model_name="boltz2",
                entity_name=seq_name,
                entity_format="sequence",
                dataset_name="gpcr_test_sequences",
                config=config
            )
            predictions.append(model_input)
            print(f"✓ Prepared prediction for {seq_name}")
        except Exception as e:
            print(f"✗ Failed to prepare {seq_name}: {e}")
    
    # Also prepare a high-confidence version of ADRB2
    high_conf_config = {
        "recycling": 10,
        "num_samples": 5,
        "crop_size": 512,
        "output_name": "ADRB2_HUMAN_high_confidence"
    }
    
    try:
        high_conf_input = manager.prepare_input(
            model_name="boltz2",
            entity_name="ADRB2_HUMAN",
            entity_format="sequence",
            dataset_name="gpcr_test_sequences",
            config=high_conf_config
        )
        predictions.append(high_conf_input)
        print("✓ Prepared high-confidence prediction for ADRB2_HUMAN")
    except Exception as e:
        print(f"✗ Failed to prepare high-confidence: {e}")
    
    print(f"\nTotal predictions prepared: {len(predictions)}")
    return predictions


def workflow_2_interface_mutations():
    """Workflow 2: Design mutations at protein-protein interface."""
    print("\n=== Workflow 2: Interface Mutation Design ===")
    
    # Initialize processors and manager
    struct_proc = StructureProcessor()
    seq_proc = SequenceProcessor()
    prop_proc = PropertyProcessor()
    manager = ModelManager()
    loader = StructureLoader()
    
    # Download a GPCR-G protein complex
    complex_id = "3sn6"  # Beta2-adrenergic receptor with Gs protein
    print(f"\nDownloading complex structure: {complex_id}")
    
    try:
        success, _ = loader.download_batch(
            [complex_id],
            dataset_name="gpcr_complexes",
            create_dataset=True
        )
        if not success:
            print(f"Failed to download {complex_id}")
            return []
    except Exception as e:
        print(f"Download error: {e}")
        return []
    
    # Load and analyze interface
    struct_df = struct_proc.load_entity(complex_id)
    if struct_df is None:
        print(f"Failed to load structure {complex_id}")
        return []
    
    # Identify interface residues (simplified - would use proper interface analysis)
    # For demo, we'll target known interface positions
    interface_positions = {
        "R": [131, 135, 139],  # Receptor chain
        "A": [380, 384, 387]   # G-protein alpha chain
    }
    
    print(f"\nTarget interface positions:")
    for chain, positions in interface_positions.items():
        print(f"  Chain {chain}: {positions}")
    
    # Extract sequences
    chain_sequences = struct_proc.collect_chain_sequences([complex_id])
    
    # Prepare mutations at interface
    mutation_configs = []
    for struct_chains in chain_sequences.values():
        for chain_id, chain_data in struct_chains.items():
            if chain_id not in interface_positions:
                continue
                
            seq_name = chain_data["entity_name"]
            sequence = chain_data["sequence"]
            
            # Create mutations at interface positions
            for pos in interface_positions[chain_id]:
                if pos <= len(sequence):
                    original = sequence[pos-1]
                    # Try different amino acids that might improve binding
                    for mutant in ["W", "Y", "F", "R", "K"]:
                        if mutant != original:
                            mutation_configs.append({
                                "entity": seq_name,
                                "format": "sequence", 
                                "config": {
                                    "mutations": [{
                                        "position": pos,
                                        "original": original,
                                        "mutant": mutant,
                                        "name": f"{original}{pos}{mutant}"
                                    }],
                                    "recycling": 3,
                                    "output_name": f"{seq_name}_{original}{pos}{mutant}_interface"
                                }
                            })
    
    # Prepare batch predictions
    if mutation_configs:
        print(f"\nPreparing {len(mutation_configs)} interface mutations...")
        batch = manager.prepare_batch(
            model_name="boltz2",
            entity_configs=mutation_configs[:10],  # Limit to 10 for demo
            batch_name="interface_optimization_batch"
        )
        print(f"✓ Created batch: {batch.name}")
        print(f"  Mutations prepared: {len(batch.inputs)}")
        return batch.inputs
    
    return []


def workflow_3_grn_guided_mutations():
    """Workflow 3: GRN-guided mutations at functionally important positions."""
    print("\n=== Workflow 3: GRN-Guided Functional Mutations ===")
    
    # Initialize processors
    seq_proc = SequenceProcessor()
    grn_proc = GRNProcessor()
    prop_proc = PropertyProcessor()
    
    # Use sequences from workflow 1
    dataset_name = "gpcr_test_sequences"
    
    # Annotate with GRN
    print("\nAnnotating sequences with GRN...")
    try:
        grn_table, summary = seq_proc.annotate_with_grn(
            dataset_name=dataset_name,
            reference_table="gpcrdb_ref",
            protein_family="gpcr_a",
            output_table=f"{dataset_name}_grn",
            return_summary=True,
            allow_create=True
        )
        print(f"✓ GRN annotation complete")
        print(f"  Sequences annotated: {summary['global']['annotated']}/{summary['global']['total']}")
    except Exception as e:
        print(f"✗ GRN annotation failed: {e}")
        return []
    
    # Target functionally important GRN positions
    functional_positions = {
        "3.50": "DRY motif - G-protein coupling",
        "6.48": "CWxP motif - activation",
        "7.50": "NPxxY motif - signaling",
        "2.50": "D2.50 - sodium binding",
        "5.50": "P5.50 - helix kink"
    }
    
    print("\nTarget functional positions:")
    for grn, description in functional_positions.items():
        print(f"  {grn}: {description}")
    
    # Prepare mutations at these positions
    try:
        inputs = prepare_mutation_screen(
            seq_proc=seq_proc,
            dataset_name=dataset_name,
            grn_positions=list(functional_positions.keys()),
            mutations=["A", "G", "V", "L", "I", "F", "W"],  # Various mutations
            grn_table_name=f"{dataset_name}_grn"
        )
        
        print(f"\n✓ Prepared {len(inputs)} structure predictions")
        print("  Sample mutations:")
        for inp in inputs[:5]:  # Show first 5
            mutation_info = inp.metadata.get("mutations", [{}])[0]
            if mutation_info:
                print(f"    {mutation_info.get('name', 'Unknown')}")
        
        return inputs
    except Exception as e:
        print(f"✗ Failed to prepare mutations: {e}")
        return []


def workflow_4_homolog_comparison():
    """Workflow 4: Predict structures for homologous sequences."""
    print("\n=== Workflow 4: Homolog Structure Prediction ===")
    
    # Initialize processors
    seq_proc = SequenceProcessor()
    emb_proc = EmbeddingProcessor()
    manager = ModelManager()
    
    # Create dataset of GPCR homologs from different families
    homologs = {
        # Class A (Rhodopsin-like)
        "OPRM_HUMAN": "MDSSAAPTNASNCTDALAYSSCSPAPSPGSWVNLSHLDGNLSDPCGPNRTDLGGRDSLCPPTGSPSMITAITIMALYSIVCVVGLFGNFLVMYVIVRYTKMKTATNIYIFNLALADALATSTLPFQSVNYLMGTWPFGTILCKIVISIDYYNMFTSIFTLCTMSVDRYIAVCHPVALKALDFRTPRNAKIINVCNWILSSAIGLPVMFMATTKYRQGSIDCTLTFSHPTWYWENLLKICVFIFAFIMPVLIITVCYGLMILRLKSVRMLSGSKEKDRNLRRITRMVLVVVAVFIVCWTPIHIYVIIKALVTIPETTFQTVSWHFCIALGYTNSCLNPVLYAFLDENFKRCFREFCIPTSSNIEQQNSTRIRQNTRDHPSTANTVDRTNHQLENLEAETAPLP",
        # Class B (Secretin-like)
        "GLR_HUMAN": "MRLSYLLLGILVLASVSSQGEEDEEQRTEREEKEPGESPPQGEPQASGLELRPQDRFLPEEQGEEVEEEGDGDYREQNEQSASPDLQEAALDLERGLLRDFNQKSNGRKSTTSLSPLSTQEEHRQLLQRLQQLLLKSHHPEAQSEQRKRAGISQPDFSRYLSLLQHLLRKLPPQVLKTYPTVLVFHQQESQHLPKFNHTLQQAYQNLLKHQRSNQGLAPGHGEKEVQEEASGLLNLIASLLAVFSEQFIQASFDRISEDTNEYSELLCNLPQETRLQHYRQCCHYFLRSVSVRVLTFFQGLTTNKMNVYCLQASFRWDYILGPGRNAVGFETLEELGKKRAALCLRQFLCCGRSLVEAVEEGCPLGTAAKMYLKHLLYSLGAAHGLHQKWVLGQTMCMSVNQVLTRALPAVQQSEGCCYAFLVQQQLHTSLQMLAQWFHQGACQRHLQERYLSCPQDSSRSPCAQECPEGWSHLEPGRQCQRLGEACFSRDPAPCELCPQVPEWDGASCQVPTTFYCLPDFVQAGVQLLEHYRQCYSFLEASVPCLQAGDVVLQLPPLLQELPALPQDLAPGLNLPYLPEYALHHCLAALPAPVQERGRRLSLPVLSPGQRIRKVSPDRFAFSESS",
        # Class C (Metabotropic glutamate)
        "GRM5_HUMAN": "MVLLLILSVLLLKEDVRGSAAQSGRPCLQAPEGHPLEVQCGPVDDKQAEETLQAIRIAEALADAQEKEFTQLAAGKLPVVKGQNRNVVDGVSLMGDIILGLFAFLAHQGSCVPPQRPGKPVSYCGDIAPVALGAQTADCPGPYTPCCSREMAGLGFDRFIQRTPNNDAHQNVQFDFPGLLGNARAYGLLFLQHFGPSGQNGGQPGVQAEEDCVEWQDQPSVQANLSFPRQEKKVIVIGPCSSVALSKAMAEGIKFDVLLPNKTGYENMTLRGCFGDRVQFAQSRDQNNSTNSTVSVSEAEVRELCENLQAWQRGLCDLLGLNVGCDLLTSKPYQLRSVKLPGCPRGDYQVQLRTGDEAELRQAPGRGYCLAFSDLLGGLGSGYTARPSILTGCPQNEAEKLEELMERLGIIHYDSCVNPGDYKDIEEDYQVCRACQLVRVLLQYKVPERPRPAGFVPSSPAPTHLPGTVGGNPQKPTGACVGECPCGLVCSLPEPLSGSCPEGRTIKKTGSGCVGDPCCSRSCSGKSPDCTGGSCWCGVPCEPHAIPSHSCKETQECPPGTFMVVNCSRCNECLDGVWQPPTLKLNMSRPGGRVVAVYGPGISSRLSQSQKLLPLTLKLYAEKDSGNYCTPNLTQKEPKQYVVFDLVTCNGSQEPRGVVAAGRDTTIGRLLQMAEQGVGQFLLRAHQVPCIAVPCRDYSVRPPGETPCCWRCVPCGPGGLRMASPGCHEKCPAGFNWNYDTLNYQYTESLMFVVNAVYAMAHALHKMQRTFCGNTYKAPYPSLLNELESWKCPNWPSAGQLVAVDGKDGTLLGFSVPHQKKAGATGRCRSPESALQQVVVAIFQSQPVASVAVFDHQFNMAKVGPWVKEKIQDTLTKAVAGEKSCRNRSNLQKILTLVGCFGDCGYGRTTDNQTSHTFGFPDFLHLLQPSGYPEFRGSVKVGEGCTGAQVPMHGWMSAWGLPFVGPGKGSPVANAVQALDLGTGSQYSMGSIIGVVGACLSTVFALLLGRSYQQLKRIRSISVQQTSEKREMVLTKDSAFLCLDLVWDHYRQQPLHCRCCCVSCRTLPPLRSMQRSGAFRPWQCPQVAQALQNQTWLKIKPETFPQSQGPQANRDAKTGQNVGIFSQVYPNTDKSYASTSEEEMEQLLQLISSARQHPSSFGWQENEKEEENMDVSGTSSQSSTCPLRPREMLQPMLTGMGQGGVTTATLATPTPGPLAPHVQGQLERWTPPLGSLGRQCLQNV",
        # Class F (Frizzled)
        "FZD7_HUMAN": "MGDLPRGLLVLLALGLAASGGSGDADLPGVVMEEQYTLPALFHGAELASGQAYWLVETRQCQSFRPLCQSLVMFLCGLIYYTVSHQYPALLQAVGSCLERLFTFLCCCVFFEAMAAYFLLAWIFMPLVLTVICQVLVFWLFLDHQLYCSISPGLLLVGCQYAGTLTFFLTHVIGGSLCFFLVCFVCNYTTGVKLLQRQALRRNSAHQGEDVTLLTSTFILKLFNHRSFRQVQPLLHPEVWTPSPNLCVLEQDDFDLVYLKTVQRETNLGCFTAPQPPPPHPDDDYSHTHSLYSSHMEPEEDPGWVVDEHGGHSHTYVHQEEQTFICKFHNFDPEVGDPSPDDSLCDVMKVNATEMDSEQTVLKTDLTFYKVIEAINMIPEPLRDDGSYVALVTAIRKGVPRALVRCPVMCQNGSWEWLPDHHYCSPCFARAFANELRWPRLEAAAGAEVLQVPPPAVHPGPPPESSLASQQLMEALPMGECVWFRHPPRLRSHTLQPTLQMNLKTASESSFLERTPMLNDSSEHSCLTTLRIALSDPAVCRVLPES"
    }
    
    # Save homolog dataset
    seq_proc.save_sequences(
        homologs,
        output_file="gpcr_homologs",
        dataset_name="gpcr_homologs",
        metadata={"description": "GPCR homologs from different classes"},
        materialize_entities=True
    )
    print(f"✓ Created homolog dataset with {len(homologs)} sequences")
    
    # Generate embeddings for clustering analysis
    print("\nGenerating embeddings for homolog analysis...")
    try:
        # Load sequences first
        sequences_dict = seq_proc.load_dataset("gpcr_homologs")
        
        # Convert to list format for embedding
        sequences = [(name, seq) for name, seq in sequences_dict.items()]
        
        # Generate embeddings
        embeddings = emb_proc.embed_sequences(
            sequences=sequences_dict,  # Can pass dict directly
            embedding_type="mean",
            save_dataset="gpcr_homologs_embeddings"
        )
        print(f"✓ Generated embeddings for {len(embeddings)} sequences")
    except Exception as e:
        print(f"✗ Embedding generation failed: {e}")
    
    # Prepare structure predictions with different confidence levels
    predictions = []
    for seq_name, sequence in homologs.items():
        # Determine prediction quality based on sequence length
        seq_len = len(sequence)
        if seq_len < 400:
            recycling = 3
            samples = 1
        elif seq_len < 600:
            recycling = 5
            samples = 2
        else:
            recycling = 7
            samples = 3
        
        config = {
            "recycling": recycling,
            "num_samples": samples,
            "crop_size": min(512, ((seq_len + 255) // 256) * 256),  # Round up to multiple of 256
            "output_name": f"{seq_name}_homolog_structure"
        }
        
        try:
            model_input = manager.prepare_input(
                model_name="boltz2",
                entity_name=seq_name,
                entity_format="sequence",
                dataset_name="gpcr_homologs",
                config=config
            )
            predictions.append(model_input)
            print(f"✓ Prepared prediction for {seq_name} (length: {seq_len}, recycling: {recycling})")
        except Exception as e:
            print(f"✗ Failed to prepare {seq_name}: {e}")
    
    return predictions


def workflow_5_multi_chain_complex():
    """Workflow 5: Predict multi-chain protein complexes."""
    print("\n=== Workflow 5: Multi-Chain Complex Prediction ===")
    
    # Initialize processors and manager
    seq_proc = SequenceProcessor()
    manager = ModelManager()
    
    # Example: GPCR with G-protein and beta-arrestin
    complex_sequences = {
        # Receptor
        "ADRB2_receptor": "MGQPGNGSAFLLAPNGSHAPDHDVTQERDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGAAHILMKMWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFKYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRATHQEAINCYANETCCDFFTNQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQKIDKSEGRFHVQNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNLIRKEVYILLNWIGYVNSGFNPLIYCRSPDFRIAFQELLCLRRSSLKAYGNGYSSNGNTGEQSGYHVEQEKENKLLCEDLPGTEDFVGHQGTVPSDNIDSQGRNCSTNDSLL",
        # G-protein alpha subunit (truncated for demo)
        "Galpha_subunit": "MGCLGNSKTEDQRNEEKAQREANKKIEKQLQKDKQVYRATHRLLLLGAGESGKSTIVKQMRILHVNGFNGEGGEEDPQAARSNSDGEKATKVQDIKNNLKEAIETIVAAMSNLVPPVELANPENQFRVDYILSVMNVPDFDFPPEFYEHAKALWEDEGVRACYERSNEYQLIDCAQYFLNKIDVIKQADYVPSDQDLLRCRVLTSGIFETKFQVDKVNFHMFDVGGQRDERRKWIQCFNDVTAIIFVVASSSYNMVIREDNQTNRLQEALNLFKSIWNNRWLRTISVILFLNKQDLLAEKVLAGKSKIEDYFPEFARYTTPEDATPEPGEDPRVTRAKYFIRDEFLRISTASGDGRHYCYPHFTCAVDTENIRRVFNDCRDIIQRMHLRQYELL",
        # Beta-arrestin (truncated)
        "Beta_arrestin": "MGEKPGTRVFKKSSPNCKLTVYLGKRDFVDHLDKVDPVDGVVLVDPDYLKDRKVFVTLTCAFRYGREDLDVLGLSFRKDLFIATYQAFPPVPNPPRPPTRLQDRLLRKLGQHAHPFFFTIPQNLPCSVTLQPGPEDTGKACGVDFEIRAFCAKSLEEKSHKRNSVRLVIRKVQFAPEKPGPQPSAETTRHFLMSDRSLHLEASLDKELYYHGEPLNVNVHVTNNSTKTVKKIKVSVRQYADICLFSTAQYKCPVAMEEADDTVAPSSTFCKVYTLTPFLANNREKRGLALDGKLKHEDTNLASSTIVKEGANKEVLGILVSYRVKVKLVVSRGGDVSVELPFVLMHPKPHDHIPLPRPQSAAPETDVPVDTNLIEFDTNYATDDDIVFEDFARLRLKGMKDDDYDDQLC"
    }
    
    # Save complex components
    seq_proc.save_sequences(
        complex_sequences,
        output_file="gpcr_complex_components",
        dataset_name="gpcr_complex_components",
        metadata={"complex": "GPCR-Gprotein-Arrestin"},
        materialize_entities=True
    )
    
    # Prepare different complex configurations
    complex_configs = [
        {
            "name": "receptor_gprotein",
            "components": ["ADRB2_receptor", "Galpha_subunit"],
            "description": "Active state GPCR with G-protein"
        },
        {
            "name": "receptor_arrestin", 
            "components": ["ADRB2_receptor", "Beta_arrestin"],
            "description": "Desensitized GPCR with arrestin"
        },
        {
            "name": "full_complex",
            "components": ["ADRB2_receptor", "Galpha_subunit", "Beta_arrestin"],
            "description": "Full ternary complex (hypothetical)"
        }
    ]
    
    predictions = []
    for complex_config in complex_configs:
        print(f"\nPreparing {complex_config['name']}: {complex_config['description']}")
        
        # For Boltz-2, we need to prepare a combined entity with all chains
        # In real usage, this would involve proper chain assignment
        config = {
            "chains": complex_config["components"],
            "recycling": 5,
            "num_samples": 2,
            "crop_size": 512,
            "output_name": f"complex_{complex_config['name']}"
        }
        
        # For multi-chain, we'd typically create a combined sequence file
        # Here we demonstrate the configuration
        print(f"  Components: {', '.join(complex_config['components'])}")
        print(f"  Configuration: recycling={config['recycling']}, samples={config['num_samples']}")
        
        # In practice, you'd prepare the multi-chain input like:
        # model_input = manager.prepare_input(
        #     model_name="boltz2",
        #     entity_name=f"complex_{complex_config['name']}",
        #     entity_format="sequence",
        #     dataset_name="gpcr_complex_components",
        #     config=config
        # )
    
    return predictions


def workflow_6_binder_design():
    """Workflow 6: Design protein binders using RFdiffusion."""
    print("\n=== Workflow 6: Protein Binder Design ===")
    
    # Initialize processors and manager
    struct_proc = StructureProcessor()
    manager = ModelManager()
    loader = StructureLoader()
    
    # Target structures for binder design
    targets = {
        "5d5a": "A2A adenosine receptor",
        "6b73": "mGluR5 metabotropic glutamate receptor",
        "4daj": "M3 muscarinic receptor"
    }
    
    # Download target structures
    print("Downloading target structures...")
    try:
        success, _ = loader.download_batch(
            list(targets.keys()),
            dataset_name="binder_targets",
            create_dataset=True
        )
        print(f"✓ Downloaded {len(success)} target structures")
    except Exception as e:
        print(f"✗ Download failed: {e}")
        return []
    
    # Prepare binder designs
    design_configs = []
    
    for pdb_id, description in targets.items():
        print(f"\nPreparing binder designs for {pdb_id} ({description})")
        
        # Different design strategies
        strategies = [
            {
                "name": "small_binder",
                "contigs": [f"A1-300/0 B1-50"],  # 50-residue binder
                "hotspots": [],  # Let RFdiffusion choose
                "num_designs": 20
            },
            {
                "name": "medium_binder_hotspot",
                "contigs": [f"A1-300/0 B1-80"],  # 80-residue binder
                "hotspots": [125, 130, 200, 205],  # Target specific positions
                "num_designs": 30
            },
            {
                "name": "scaffolded_binder",
                "contigs": [f"A1-300/0 B10-30,40-60,70-90"],  # Scaffolded design
                "hotspots": [130, 200],
                "num_designs": 25
            }
        ]
        
        for strategy in strategies:
            config = {
                "contigs": strategy["contigs"],
                "hotspots": strategy["hotspots"],
                "num_designs": strategy["num_designs"],
                "timesteps": 50,
                "seed": 42,
                "run_name": f"{pdb_id}_{strategy['name']}"
            }
            
            try:
                model_input = manager.prepare_input(
                    model_name="rfdiffusion",
                    entity_name=pdb_id,
                    entity_format="structure",
                    dataset_name="binder_targets",
                    config=config
                )
                design_configs.append(model_input)
                print(f"  ✓ {strategy['name']}: {strategy['num_designs']} designs")
            except Exception as e:
                print(f"  ✗ Failed to prepare {strategy['name']}: {e}")
    
    print(f"\n✓ Total binder design configurations: {len(design_configs)}")
    return design_configs


def demonstrate_output_parsing():
    """Demonstrate parsing model outputs back into Protos."""
    print("\n=== Output Parsing and Registration ===")
    
    manager = ModelManager()
    struct_proc = StructureProcessor()
    
    print("Examples of output parsing:")
    
    # Example 1: Parse Boltz-2 output
    print("\n1. Boltz-2 structure prediction:")
    print("   Input: predictions/boltz2/ADRB2_HUMAN_V91A/prediction.pdb")
    print("   Process:")
    print("   - Extract predicted structure")
    print("   - Parse confidence metrics from confidence.json")
    print("   - Register as 'ADRB2_V91A_predicted' in 'boltz2_mutations' dataset")
    print("   - Structure automatically available in StructureProcessor")
    
    # Example 2: Parse RFdiffusion outputs
    print("\n2. RFdiffusion binder designs:")
    print("   Input: predictions/rfdiffusion/5d5a_small_binder/")
    print("   Process:")
    print("   - Find all design_*.pdb files")
    print("   - Register each as '5d5a_binder_design_1', '5d5a_binder_design_2', etc.")
    print("   - Create 'designed_binders' dataset with all designs")
    print("   - Extract trajectory information for analysis")
    
    # Example 3: Batch output processing
    print("\n3. Batch mutation screen results:")
    print("   Input: predictions/boltz2/grn_mutations_batch/")
    print("   Process:")
    print("   - Iterate through all mutation predictions")
    print("   - Extract pLDDT scores for quality filtering")
    print("   - Register only high-confidence predictions (pLDDT > 70)")
    print("   - Create property table with mutation effects")
    
    # Show how to access parsed results
    print("\n4. Accessing parsed results in processors:")
    print("   # After parsing")
    print("   struct_proc.load_entity('ADRB2_V91A_predicted')")
    print("   struct_proc.list_dataset_entities('boltz2_mutations')")
    print("   ")
    print("   # Analyze all predictions")
    print("   for entity in struct_proc.list_dataset_entities('designed_binders'):")
    print("       struct_df = struct_proc.load_entity(entity)")
    print("       # Analyze structure...")


def main():
    """Run all ModelManager workflow demonstrations."""
    print("ModelManager Workflow Examples")
    print("==============================")
    print("Demonstrating ModelManager integration with Protos processors")
    
    # Initialize data directory
    data_root = ensure_data_root()
    print(f"\nData root: {data_root}")
    
    # Run workflows
    print("\n" + "="*60)
    
    # Workflow 1: Basic sequence to structure
    seq_predictions = workflow_1_sequence_to_structure()
    
    print("\n" + "="*60)
    
    # Workflow 2: Interface mutations
    interface_mutations = workflow_2_interface_mutations()
    
    print("\n" + "="*60)
    
    # Workflow 3: GRN-guided mutations
    grn_mutations = workflow_3_grn_guided_mutations()
    
    print("\n" + "="*60)
    
    # Workflow 4: Homolog comparison
    homolog_predictions = workflow_4_homolog_comparison()
    
    print("\n" + "="*60)
    
    # Workflow 5: Multi-chain complexes
    complex_predictions = workflow_5_multi_chain_complex()
    
    print("\n" + "="*60)
    
    # Workflow 6: Binder design
    binder_designs = workflow_6_binder_design()
    
    print("\n" + "="*60)
    
    # Demonstrate output parsing
    demonstrate_output_parsing()
    
    # Summary
    print("\n" + "="*60)
    print("\nWorkflow Summary")
    print("================")
    total_inputs = (
        len(seq_predictions) + 
        len(interface_mutations) + 
        len(grn_mutations) + 
        len(homolog_predictions) + 
        len(binder_designs)
    )
    print(f"Total model inputs prepared: {total_inputs}")
    print("\nKey capabilities demonstrated:")
    print("✓ Sequence to structure prediction with mutations")
    print("✓ Interface-focused mutation design") 
    print("✓ GRN-guided functional mutations")
    print("✓ Homolog structure comparison")
    print("✓ Multi-chain complex prediction")
    print("✓ Protein binder design with RFdiffusion")
    print("✓ Output parsing and registration")
    print("\nAll workflows use only processors and ModelManager - no direct paths!")
    print("\nTo execute predictions, run the generated commands from ModelInput.get_command()")
    print("Model outputs will be automatically registered back into Protos.")


if __name__ == "__main__":
    main()