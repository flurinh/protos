# Define max bond length dictionary outside the function
BOND_LENGTH_DICT = {
    ('C', 'C'): 1.65,  # Increased to account for some double bonds
    ('C', 'N'): 1.60,  # Increased to account for various bond types
    ('C', 'O'): 1.55,  # Increased for carbonyl bonds
    ('C', 'S'): 1.95,  # Increased slightly
    ('C', 'H'): 1.20,  # Increased slightly
    ('N', 'N'): 1.55,  # Increased for various types of N-N bonds
    ('N', 'O'): 1.50,  # Increased slightly
    ('N', 'S'): 1.90,  # N-S bonds can vary
    ('N', 'H'): 1.10,  # Increased slightly
    ('O', 'O'): 1.60,  # Increased for peroxide bonds
    ('O', 'S'): 1.80,  # O-S bonds in sulfoxides and sulfones
    ('O', 'H'): 1.05,  # Increased slightly
    ('S', 'S'): 2.10,  # S-S bonds can be quite long
    ('S', 'H'): 1.45,  # Increased slightly
    ('H', 'H'): 0.80,  # H-H bond (rare, but possible in some cases)
    ('BR', 'C'): 2.00,  # C-Br bond
    ('BR', 'N'): 2.10,  # N-Br bond (rare)
    ('BR', 'O'): 2.00,  # O-Br bond (rare)
    ('BR', 'S'): 2.30,  # S-Br bond
    ('BR', 'H'): 1.50,  # H-Br bond
    ('BR', 'BR'): 2.40,  # Br-Br bond (rare)
    ('CL', 'C'): 1.90,  # C-Cl bond
    ('CL', 'N'): 2.00,  # N-Cl bond (rare)
    ('CL', 'O'): 1.90,  # O-Cl bond (rare)
    ('CL', 'S'): 2.20,  # S-Cl bond
    ('CL', 'H'): 1.40,  # H-Cl bond
    ('CL', 'CL'): 2.20,  # Cl-Cl bond (rare)
    ('BR', 'CL'): 2.30,  # Br-Cl bond (very rare)
}

