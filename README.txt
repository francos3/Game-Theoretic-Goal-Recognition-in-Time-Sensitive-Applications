## Authors
This code is authored by Sara Bernardini, Fabio Fagnini and Santiago Franco for the AAMAS 2025 paper titled: "Game-Theoretic Goal Recognition in Time-Sensitive Applications".

## License

This project is licensed under the MIT License. See the LICENSE file for details.

##Instructions

To configure, build and run any of the experiments in the paper.

1. To compile (Works in Ubuntu 22.04 and should work in most linux configurations), run from this folder: 
cmake .;make -j4

3. To prebuild the yen paths for the 30x30 orthogonal grid, populate SEED (1-199 in the paper) usingf the maxium number of DESTINATIONS (5 in the paper):  
./yen AG filename_SaT="Grid-30.csv" S=$SEED R=5 stocastic_fictitious_play Budget=10 acyclic save_paths >& log_createpaths.txt

4. To run Best Observer Response on the 30x30 orthogonal grid, populate the saturating factor "$t" (1-4 in the paper), SEED (1-199 in the paper) and number of DESTINATIONS (2-5 in the paper); the code will load the pregenerated paths up to the number of DESTINATIONS selected:   
./yen AG filename_SaT="Grid-30.csv" S=$SEED R=$DESTINATIONS stocastic_fictitious_play Stg=0 Budget=10 acyclic i=10000 sat=$t type=3 load_paths >& log_BestObserverResponse.txt

5. To run Fixed Response on the 30x30 orthogonal grid, populate perimeter radious $RAD (perimeter=2) in the paper, $t (1-4 in the paper), SEED (1-199 in the paper) and number of DESTINATIONS (2-5 in the paper):
./yen AG filename_SaT="Grid-30.csv" S=$SEED R=$DESTINATIONS stocastic_fictitious_play Stg=0 Budget=10 acyclic i=1 sat=$t type=3 load_paths fix_rad=$RAD fixed_strategy >& log_FixedResponse.txt

6. Repeat step 1 to prebuild the yen paths for the second set of experiments with the MAPF benchmark of Shanghai. You can reuse steps 3-5 to rerun the observer strategies. You need to replace the filename by using filename_SaT="Shanghai_0_256.csv" instead of filename_SaT="Grid-30.csv". Also set i=1000 iterations instead of i=10000 as we limit the number of iterations to 1000. Note that the Best Observer Response can take up to 3 or 4 days in some instances to run due to how large this grids are (256x256) are comparted to the first set of experiments in a 30x30 grid.

7. Note that the A Priori Observer Response will always obtain a Nash Equilibira reward of 1/$DESTINATIONS because it always chooses the farthest destination, and it makes the selection from the origin, there is no need to simulate it.

