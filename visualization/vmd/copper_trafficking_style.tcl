# Generic VMD representation for copper-trafficking protein complexes.
# Load your PSF/PDB/DCD before sourcing this file. Review the selections below
# against the system's residue and segment naming; historical archives used
# several different names for copper atoms and protein segments.

# User-editable VMD selection expressions.
set protein_selection "protein"
set copper_selection "resname CU1 CU CU1P COP"
set water_selection "water"

# Clear existing representations for the top molecule only.
set molid [molinfo top]
for {set rep [expr {[molinfo $molid get numreps] - 1}]} {$rep >= 0} {incr rep -1} {
    mol delrep $rep $molid
}

# Protein cartoon.
set protein [atomselect $molid $protein_selection]
if {[$protein num] > 0} {
    mol representation NewCartoon 0.3 10.0 4.1 0
    mol color Chain
    mol selection $protein_selection
    mol material Goodsell
    mol addrep $molid
} else {
    puts stderr "WARNING: protein selection matched no atoms: $protein_selection"
}
$protein delete

# Copper/metal spheres. Validate that the expression selects only the intended metal atoms.
set copper [atomselect $molid $copper_selection]
if {[$copper num] > 0} {
    mol representation VDW 0.9 20.0
    mol color Element
    mol selection $copper_selection
    mol material AOShiny
    mol addrep $molid
} else {
    puts stderr "WARNING: copper selection matched no atoms: $copper_selection"
}
$copper delete

# Water is hidden by default to keep protein-level images uncluttered.
# Uncomment and adapt the next block when solvent visualization is required.
# mol representation Lines 1.0
# mol color Name
# mol selection $water_selection
# mol material Transparent
# mol addrep $molid

axes location Off
display projection Orthographic
puts "Applied generic copper-trafficking representation to molecule $molid."
