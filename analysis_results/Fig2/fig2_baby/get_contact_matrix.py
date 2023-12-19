import os
import numpy as np
thisdir = os.path.abspath('')
datadir=os.path.join(thisdir,'data')
# path to the main directory of the repository
#maindir = os.path.split(os.path.split(thisdir)[0])[0]

# path to the analysis_results subdirectory
#analysisdir = os.path.split(thisdir)[0]

# path to the data subdirectory
#datadir = os.path.join(os.path.split(os.path.split(thisdir)[0])[0], 'data')

def read_contact_matrix(location, country, level, setting, num_agebrackets=85): #Beijing, China, subnational, household/overall
    """
    Read in the contact for each setting.

    Args:
        location (str)        : name of the location
        country (str)         : name of the country
        level (str)           : name of level (country or subnational)
        setting (str)         : name of the contact setting
        num_agebrackets (int) : the number of age brackets for the matrix

    Returns:
        A numpy matrix of contact.
    """
    setting_type, setting_suffix = 'F', 'setting'
    if setting == 'overall':
        setting_type, setting_suffix = 'M', 'contact_matrix'


    if level == 'country':
        file_name = country + '_' + level + '_level_' + setting_type + '_' + setting + '_' + setting_suffix + '_' + '%i' % num_agebrackets + '.csv'
    else: ##我们关注这儿
        file_name = country + '_' + level + '_' + location + '_' + setting_type + '_' + setting + '_' + setting_suffix + '_' + '%i' % num_agebrackets + '.csv'

    file_path = os.path.join(datadir, 'contact_matrices', file_name)
    M = np.loadtxt(file_path, delimiter=',')
    return M