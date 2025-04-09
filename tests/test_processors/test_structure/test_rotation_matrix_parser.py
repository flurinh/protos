import pytest
import os
import re

def test_rotation_matrix_extraction():
    """Test that we can extract rotation matrix and translation vector from GTalign output."""
    
    # Sample output from GTalign with rotation matrix and translation vector
    sample_output = """gtalign 0.15.00                            Thu Mar  6 12:41:40 2025

 Command line:
gtalign --qrs=data\\temp\\cif\\5ahz_alignment_1741261299.cif --rfs=data\\temp\\cif\\5awz_alignment_1741261299.cif -o data\\temp\\gtalign\\alignment_output -s 0.0 --speed=9 --depth=2 --sort=0 --nhits=100 -v

 Query (245 residues):
data\\temp\\cif\\5ahz_alignment_1741261299.cif Chn:A (M:1)

 Searched:
data\\temp\\cif\\5awz_alignment_1741261299.cif
             1 structure(s)
             236 total residues

 Legend:
TM-score (Refn./Query), Reference/Query length-normalized TM-score
2TM-score, secondary TM-score excluding unmatched helices
d0 (Refn./Query), Normalizing inter-residue distance d0 for Reference/Query
RMSD, Root-mean-square deviation (A); Chn, Chain; (M), Model
+, pairs of aligned residues within a distance of 5 A
Secondary structure: h, Helix; e, Strand; t, Turn


                             Query_length-normalized_TM-score| RMSD|Reference_alignment_boundaries|
                  Reference_length-normalized_TM-score|    Query_alignment_boundaries|
   No.|                   Reference_description|        #Aligned_residues|             Reference_length|

     1 ..._alignment_1741261299.cif Chn:A (M:1) 0.8764 0.8463  2.44   229     2-243       3-236     236



1.           
>data\\temp\\cif\\5awz_alignment_1741261299.cif Chn:A (M:1)
  Length: Refn. = 236, Query = 245

 TM-score (Refn./Query) = 0.87641 / 0.84628, d0 (Refn./Query) = 5.70 / 5.80,  RMSD = 2.44 A
 Identities = 53/247 (21%), Matched = 214/247 (86%), Gaps = 18/247 (7%)

struct        hh hhhhhhhhhhhhhhhhhhhhhhhhhhhhhh    thhhhhhhhhhhhhhhhhhhhhhhht   t eeeee t  t 
Query:     2 AAV-RENALLSSSLWVNVALAGIAILVFVYMGRTIRPGRPRLIWGATLMIPLVSISSYLGLLSGL-TVGMIEMPAGHAL 79   
                   ++L+++++W V++A++++A++VF++++++++++ +RL+++++++I+++++++Y++++++    +++++      
Refn.:     3 PN-PFQTTLGTDAQW-VVFAVMALAAIVFSIAVQFRPLP-LRLTYYVNIAICTIAATAYYAMAVNGGDNKPTAGT----- 74   
struct            e thhhhhhh hhhhhhhhhhhhhhhhhh   th h hhhhhhhhhhhhhhhhhhhhhht ht   e        

struct          eee   hhhhhhhhhhhhhhhhhhhhht   thhhhhhhhhhhhhhhhhhhhh     thhhhhhhhhhhhhhhhhh
Query:    80 AGEMVRSQWGRYLTTWALSTPMILLALGLLADVDLGSLFTVIAADIGMCVTGLAAAMTTSALLFRWAFYAISCAFFVVVL 159  
              +++++++++RY ++W+++TP++LL+L+LL+++++++++++++ADI+M+++G+++A+T+  ++++W++++++C++++V++
Refn.:    75 GADERQVIYARY-IDWVFTTPLLLLDLVLLTNMPATMIAWIMGADIAMIAFGIIGAFTV--GSYKWFYFVVGCIMLAVLA 151  
struct       th  eeee  hh hhhhhhhhhhhhhhhht    thhhhhhhhhhhhhhhhhhhhh     tthhhhhhhhhhhhhhhhh

struct       hhhhhhhhh  hhhhtt  hhhhhhhhhhhhhhhhhhhhhhtttht      thhh hhhhhhhhhhhhhhhhhhhhhhh
Query:   160 SALVTDWAA--SASSAGTAEIFDTLRVLTVVLWLGYPIVWAVGVEGLALVQSVGVT-SWAYSVLDVFAKYVFAFILLRWV 236  
             ++++++  +  ++ + +++++++TL+++++VLW++YPIVW++G +G++++ +V+V+ ++A+++LD++AK+++A+++L++V
Refn.:   152 WGMINPIFKEELQKHKEYTGAYTTLLIYLIVLWVIYPIVWGLG-AGGHII-GVDVEEIIAMGILDLLAKPLYAIGVLITV 229  
struct       hhhhhhht  h    th  hhhhhhhhhhhhhhhhhhhhhhhh hh      thhhhhhhhhhhhhhhhhhhhhhhhhhh

struct       hh  hh 
Query:   237 ANNERTV 243  
             +++++  
Refn.:   230 EVVYGKL 236  
struct       hhhhh  

 Rotation [3,3] and translation [3,1] for Query:
   -0.865071   0.444565   0.232410        -13.257886
   -0.218879   0.082370  -0.972269         73.631096
   -0.451381  -0.891951   0.026050         21.836479

Query length: 245
Total length of reference structures: 236
Search space: 57820
Time elapsed from process initiation: 0.537817 sec
Query batch execution time: 0.020580 sec
Query batch size: 1
Device: "NVIDIA GeForce RTX 3050 Ti Laptop GPU"
"""

    # Pattern to extract the rotation matrix and translation vector
    rotation_pattern = r"\s+Rotation\s+\[3,3\]\s+and\s+translation\s+\[3,1\]\s+for\s+Query:\s*\n((?:\s*[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n){3})"
    
    # Use regex to extract the matrix text
    match = re.search(rotation_pattern, sample_output)
    assert match is not None, "Pattern not found in output"
    
    matrix_text = match.group(1)
    lines = matrix_text.strip().split('\n')
    
    print("\nExtracted matrix lines:")
    for line in lines:
        print(f"  '{line}'")
    
    # We should have 3 lines for a 3x3 rotation matrix
    assert len(lines) == 3
    
    # Parse the rotation matrix and translation vector
    rotation_matrix = []
    translation_vector = []
    
    for line in lines:
        # Extract all floating point numbers from the line
        nums = re.findall(r'[-]?\d+\.\d+', line)
        assert len(nums) >= 4, f"Expected at least 4 numbers, got {len(nums)}"
        
        # First 3 are rotation matrix, 4th is translation vector
        rotation_matrix.append([float(nums[0]), float(nums[1]), float(nums[2])])
        translation_vector.append(float(nums[3]))
    
    # Check the parsed values
    print("\nParsed rotation matrix:")
    for row in rotation_matrix:
        print(f"  {row}")
    
    print("\nParsed translation vector:")
    print(f"  {translation_vector}")
    
    # Verify the extracted values
    assert len(rotation_matrix) == 3
    assert len(translation_vector) == 3
    
    # Check a specific value for correctness
    assert abs(rotation_matrix[0][0] - (-0.865071)) < 0.000001
    assert abs(translation_vector[2] - 21.836479) < 0.000001