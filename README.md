# epidemic-game-coevolution
1. Figures 2 and 7 consist of four subplots each, and the source code corresponds to the folders *_baby, *_child, *_adult, and *_elderly.

2. main1.py to main10.py correspond to experimental results under different relative vaccination costs. To expedite the code execution, we utilized "G_C_percentage_list_all=np.linspace(0,1,20)" to divide the relative vaccination cost into 20 values. main1.py to main10.py sequentially correspond to two of these relative vaccination cost values.

3. The file "get_contact_matrix.py" constructs the contact matrix based on age structure.

3. The file "get_seed.py" generates initial infections.
