# inapps_sim
This is a character matrix simulator based on Lewis (2001) Mk model with the inclusion of logically dependent characters and multiple levels of hierarchy.

For the use of the following simulator, use the function simulate hierarchical(). 
The following inputs are mandatory:
1. A tree in Newick format (my_tree)
2. A list of the amount of characters that is desired to simulate (char_list). In this list, each number corresponds to the number of characters that will be generated in a certain level of hierarchy. For example, a char_list = [100,30] simulates 130 characters, from which 100 are standard characetrs and 30 are logically dependent characetrs to one of the first 100.
3. Degrees of freedom for the Chi-squaed distribution (df). The degrees of freedom of the distribution determine how concentrated will be the distribution of parent characters. A smaller df will generate dataset in chich most characters arise from a smaller number of parent characters.
