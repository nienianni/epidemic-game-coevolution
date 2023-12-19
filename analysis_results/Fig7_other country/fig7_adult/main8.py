from theory import *
import numpy as np
from get_contact_matrix import *
from get_seed import *
from multiprocessing import Pool
from tools import *
from tqdm import tqdm

num_agebrackets=85

###流行病传播参数
initial_infected_age = 20  # some initial age to seed infections within the population
percent_of_initial_infected_seeds = 1e-5 #初始感染人口占比

lambda1=0.3
lambda2=0.4
r=0.8#恢复率

###疫苗接种+年龄博弈参数
x0_ABCD_list=[0.5,0.5,0.5,0.5]#所有年龄初始疫苗接种
G_A_percentage=0.2 #疫苗接种相对成本
G_B_percentage=0.5
G_C_percentage_list_all=np.linspace(0,1,20)
G_C_percentage_list=G_C_percentage_list_all[14:16]
G_D_percentage=0.2

beta_list=[10]#选择强度

#n_simulation=1#模拟总次数
t_max=100 #SIR传播步数
game_t_max=2800#季节数


f_R='result_R_a8.txt'
f_Rv='result_Rv_a8.txt'

if __name__ == '__main__':
    M=read_contact_matrix('New_York', 'United_States', 'subnational', 'overall',num_agebrackets=85)
    H=read_contact_matrix('New_York', 'United_States', 'subnational', 'household',num_agebrackets=85)
    seed_num,ages=get_seeds('New_York', 'United_States', 'subnational', percent_of_initial_infected_seeds, initial_infected_age)
    for beta in beta_list:
          args = []
          args.append([num_agebrackets, ages, seed_num, initial_infected_age, lambda1, lambda2, r, M, H, t_max, game_t_max, x0_ABCD_list, G_A_percentage, G_B_percentage, G_C_percentage_list, G_D_percentage, beta])
          R_list,R_A_list, R_B_list, R_C_list, R_D_list, Rv_list, Rv_A_list, Rv_B_list, Rv_C_list,Rv_D_list=theory(*args)
          print_file_1(f_R, R_list)
          print_file_1(f_R, R_A_list)
          print_file_1(f_R, R_B_list)
          print_file_1(f_R, R_C_list)
          print_file_1(f_R, R_D_list)

          print_file_1(f_Rv, Rv_list)
          print_file_1(f_Rv, Rv_A_list)
          print_file_1(f_Rv, Rv_B_list)
          print_file_1(f_Rv, Rv_C_list)
          print_file_1(f_Rv, Rv_D_list)

