
from protos.io.paths import ProtosPaths
from protos.processing.grn.grn_utils import GRNConfigManager, get_grn_interval, sort_grns_str


paths = ProtosPaths()
print(paths.data_root)

config = GRNConfigManager(paths=paths)

# Get microbial opsin config
grn_config_str = config.get_config(protein_family='microbial_opsins', strict=True)

grns_str_str = []
for region_name, (start_grn, end_grn) in grn_config_str.items():
    print(f"GRN interval: {region_name}: {start_grn} - {end_grn}")
    region_grns = get_grn_interval(start_grn, end_grn, grns_str=None)
    grns_str_str.extend(region_grns)
grns_str_str = list(set(grns_str_str))

grns_str_str = sort_grns_str(grns_str_str)
print(f"Using {len(grns_str_str)} GRN positions for annotation")