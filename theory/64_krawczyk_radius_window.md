# Krawczyk radius-window autopsy

Module 135 certifies the coupled five-channel Hodge fixed point at 256 and 512 bits but narrowly misses inclusion in three coordinates at 160 bits.  This does not by itself distinguish a genuine failure of contraction from a poor enclosure radius at the lowest precision.

Module 136 keeps the physical map, candidate branch, native Arb fixed-Hodge solves, point preconditioner and strict Krawczyk criterion unchanged.  It varies only the radius of the interval neighborhood around the same candidate and records the coordinatewise image/box width ratios.

The purpose is diagnostic and rigorous: if any radius yields `K(X) subset interior(X)` in all five physical coordinates, a fixed point is certified at that precision.  If no radius does, the next repair must change the observer chart or linear enclosure, not the tolerance and not the Navier-Stokes model.
