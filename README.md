# membranePRE_backcalculation
This is Python Codes used to back calculate membrane PRE value for protein-membrane complex PDB files derived from a MD trajectory.

(Step One)

Build data base

    BEFORE RUN open PRE-db-calculation.py:

    1)   make sure PDB directory is set correctly in line # 5
    2)   r2.dat is in the same directory with PRE-db-calculation.py


    TO RUN:

    python PRE-db-calculation.py > The_file_Name_You_want.dat

Copy data base to the folder of ensemble selection :


(Step Two)


Ensemble selection:

    Before RUN:
        weight factor and experimental error are required in a Tab delimited Text file, 
        which can be generated and exported by MS EXCEL.

        Example:
        column 1    column 2                    column 3
        res#        experimental error          weight factor
    TO RUN example:
        python asteroids_ensemble_selection.py pre-database PRE-measured.dat error_and_weightFactor.txt 32 > run_log
    
    
