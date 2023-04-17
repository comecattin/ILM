# # Load trajectory and reference structure
# mol addfile <path/to/trajectory/and/reference> type {xtc gro} first 0 last -1 step 1 waitfor all

# # Create new molecule representation for protein
# mol new 1
# mol modselect 1 0 protein
# mol modstyle 1 0 Licorice
# mol modcolor 1 0 ColorID 0

# Disable automatic centering and orientation
# mol modauto 0 1

# Select atoms to align to z-axis
set sel [atomselect top "residue 32 to 53"]

# Calculate RMSD to reference structure and rotate protein
set ref [atomselect 1 "residue 32 to 53"]
set rmsd [measure rmsd $sel $ref]
set angle [expr {(-acos(1-$rmsd/2.0))*180.0/3.14159}]
rotate x $angle axis {1 0 0}
rotate y $angle axis {0 1 0}

# Save current orientation as reference structure
molinfo top set reference [atomselect top "protein and name CA"]

# Define function to align protein to reference structure
proc align_protein {} {
    # Load next frame of trajectory
    animate goto next

    # Align protein to reference structure
    rmsd_traj fit [atomselect top "protein"] [atomselect top "protein" frame 0] [atomselect top "reference"]

    # Rotate protein to align selected atoms with z-axis
    set rmsd [measure rmsd $sel [atomselect top "protein"]]
    set angle [expr {(-acos(1-$rmsd/2.0))*180.0/3.14159}]
    rotate x $angle axis {1 0 0}
    rotate y $angle axis {0 1 0}
}

# Call align_protein function every time a new trajectory frame is loaded
molinfo top set molinfoframecallback align_protein