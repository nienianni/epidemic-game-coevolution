import os
import pandas as pd
thisdir = os.path.abspath('')
datadir=os.path.join(thisdir,'data')


def get_ages(location, country, level, num_agebrackets=85):
    """
    Get the age count for the synthetic population of the location.

    Args:
        location (str)        : name of the location
        country (str)         : name of the country
        level (str)           : name of level (country or subnational)
        num_agebrackets (int) : the number of age brackets

    Returns:
        dict: A dictionary of the age count.
    """

    if country == 'Europe':
        country = location
        level = 'country'

    if level == 'country':
        file_name = country + '_' + level + '_level_age_distribution_' + '%i' % num_agebrackets + '.csv'
    else:
        file_name = country + '_' + level + '_' + location + '_age_distribution_' + '%i' % num_agebrackets + '.csv'
    file_path = os.path.join(datadir, 'age_population', file_name)
    df = pd.read_csv(file_path, delimiter=',', header=None)
    df.columns = ['age', 'age_count']
    ages = dict(zip(df.age.values.astype(int), df.age_count.values))
    return ages

def get_seeds(location, country, level,percent_of_initial_infected_seeds, initial_infected_age):
    ages=get_ages(location, country, level, num_agebrackets=85)
    total_population = sum(ages.values())  # 总人口数
    initial_infected_number = min(total_population * percent_of_initial_infected_seeds,ages[initial_infected_age]) #
    return initial_infected_number,ages

