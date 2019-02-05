# fibril_detection
Find raw micrographs that contain amyloid fibrils

Run with ./find_fibrils.py <mrc files>

Takes a samples of the power spectrum between 40-90 px radiating out for the center, 
averages and looks for peaks corresponding to merdianol layer lines from fibular structures.
Also has code to sample the same way in the 4.8 Angstrom reigon to look for corresponding peaks 
90 degress to the meridianol lines, but the signal to noise is generally too low for this to be effective.

Sometimes gets tricked by very strong carbon edges, but there doesn't seem to be any way around that for the time being.

Currently uses a set threshold for determining if an image has fibrils, later will make this adjustable and scale it to defocus.

