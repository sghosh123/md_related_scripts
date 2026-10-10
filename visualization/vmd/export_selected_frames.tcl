# Export selected frames from the top VMD molecule as individual PDBs.
# Load the topology/trajectory first, then edit these settings and source this file.

set output_dir "selected_frames"
set first_frame 0
set last_frame [expr {[molinfo top get numframes] - 1}]
set stride 100
set selection_text "all"

set nframes [molinfo top get numframes]
if {$nframes < 1} { error "Top molecule has no frames." }
if {$first_frame < 0 || $last_frame >= $nframes || $last_frame < $first_frame} {
    error "Invalid frame range: $first_frame..$last_frame; trajectory has $nframes frames."
}
if {$stride < 1} { error "stride must be a positive integer." }
if {![file exists $output_dir]} {
    file mkdir $output_dir
}
set sel [atomselect top $selection_text]
if {[$sel num] == 0} {
    $sel delete
    error "Selection matched no atoms: $selection_text"
}
for {set frame $first_frame} {$frame <= $last_frame} {incr frame $stride} {
    $sel frame $frame
    set filename [format "%s/frame_%06d.pdb" $output_dir $frame]
    $sel writepdb $filename
}
$sel delete
puts "Exported frames $first_frame..$last_frame (stride $stride) to $output_dir"
