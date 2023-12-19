import math
import numpy as np
import random
from matplotlib import pyplot as plt
from tqdm import tqdm
import copy

class game_theory():
    def __init__(self, num_agebrackets, t_max, ages, x0_ABCD):
        print(sum(ages.values()))
        self.t_max=t_max
        self.ages=ages
        self.num_agebrackets=num_agebrackets

        self.x_A=x0_ABCD[0]
        self.x_B = x0_ABCD[1]
        self.x_C = x0_ABCD[2]
        self.x_D=x0_ABCD[3]

        self.states = np.zeros((4, num_agebrackets, self.t_max + 1))
        self.susceptible = 0
        self.infected = 1
        self.recovered = 2
        self.vaccinated = 3

        self.uA = np.sum(list(ages.values())[0:3]) / sum(self.ages.values())
        self.uB = np.sum(list(ages.values())[3:18]) / sum(self.ages.values())
        self.uC = np.sum(list(ages.values())[18:60]) / sum(self.ages.values())
        self.uD = np.sum(list(ages.values())[60:65]) / sum(self.ages.values())

    def get_x_ABCD(self): #四个年龄段的疫苗接种人数占比
        ages=list(self.ages.values())
        x_A=np.sum(self.states[self.vaccinated, 0:3, 0]) / np.sum(ages[0:3])
        x_B=np.sum(self.states[self.vaccinated, 3:18, 0]) / np.sum(ages[3:18])
        x_C = np.sum(self.states[self.vaccinated, 18:60, 0]) / np.sum(ages[18:60])
        x_D = np.sum(self.states[self.vaccinated, 60:85, 0]) / np.sum(ages[60:85])
        return x_A, x_B, x_C, x_D


    def find_steady_state(self, seed_num, initial_infected_age, M, H, lambda1, lambda2, r):
        #### SIRRv皆为个体数量
        #print(M[84][84])
        #初始I
        self.states[self.infected][initial_infected_age][0] = seed_num
        #print(self.states[self.infected][initial_infected_age][0])
        #初始Rv
        for a in range(self.num_agebrackets):
            if a<3:
                self.states[self.vaccinated][a][0]=self.x_A*self.ages[a]#int(+0.5)
                #print(self.states[self.vaccinated][a][0])
            elif a>=3 and a<18:
                self.states[self.vaccinated][a][0] = self.x_B * self.ages[a]#int( + 0.5)
            elif a>=18 and a<60:
                self.states[self.vaccinated][a][0] = self.x_C * self.ages[a]#int( + 0.5)
            else:
                self.states[self.vaccinated][a][0] = self.x_D * self.ages[a]#int( + 0.5)
        #初始S
        for a in range(self.num_agebrackets):  # 设置t=0时，各个年龄下S态、Rv态人数
            self.states[self.susceptible][a][0] = copy.deepcopy(self.ages[a]) - self.states[self.infected][a][0] - self.states[self.vaccinated][a][0] # 年龄为a的总人数-t=0时，年龄为a的感染态人数

        # print(self.states[self.vaccinated, 3:18, 0])
        # print(self.states[self.susceptible, 3:18, 0])
        # print(self.states[self.infected, 3:18, 0])
        # print(self.states[self.recovered, 3:18, 0])

        t = 0
        while t < self.t_max:
            # print("t:",t)
            # print("S", self.states[self.susceptible, :,t].sum())
            # print("I", self.states[self.infected, :, t].sum())
            # print("R", self.states[self.recovered, :, t].sum())
            # print("Rv",self.states[self.vaccinated, :, 0].sum())
            # print(self.states[self.infected, :, -1])
            # print(self.states[self.infected, :, -1].sum())

            for i in range(self.num_agebrackets):
                A=0
                B=0
                for j in range(self.num_agebrackets):
                    A=A+M[i][j]*(self.states[self.infected][j][t]/ self.ages[j])
                    B1=0
                    for l in range(self.num_agebrackets):
                        # if(t==19):
                        #     print(i,j,l,t)
                        #     print("1",self.states[self.infected][j][t])
                        #     print("2",self.ages[j])
                        #     print("3",self.states[self.infected][l][t])
                        #     print("4",self.ages[l] )
                        B1=B1+4.11*H[i][j]*H[i][l]*(self.states[self.infected][j][t]/self.ages[j])*(self.states[self.infected][l][t]/self.ages[l])
                    B=B+B1

                self.states[self.susceptible][i][t + 1] = self.states[self.susceptible][i][t]-lambda1*self.states[self.susceptible][i][t]*A -lambda2*self.states[self.susceptible][i][t]*B
                self.states[self.infected][i][t + 1] = self.states[self.infected][i][t]-r*self.states[self.infected][i][t]+lambda1*self.states[self.susceptible][i][t]*A +lambda2*self.states[self.susceptible][i][t]*B
                self.states[self.recovered][i][t + 1] = self.states[self.recovered][i][t] + r*self.states[self.infected][i][t]

            if self.states[self.infected, :, t].sum()<1: #如果所有年龄的感染态人数之和<1,结束时间演化
                self.final_t = t
                t=self.t_max-1
            t=t+1

        ages = list(self.ages.values())
        #print(self.states[self.recovered, :, self.final_t])
        self.result_R_A = np.sum(self.states[self.recovered, 0:3, self.final_t]) / np.sum(ages[0:3])
        self.result_R_B = np.sum(self.states[self.recovered, 3:18, self.final_t]) / np.sum(ages[3:18])
        self.result_R_C = np.sum(self.states[self.recovered, 18:60, self.final_t]) / np.sum(ages[18:60])
        self.result_R_D = np.sum(self.states[self.recovered, 60:85, self.final_t]) / np.sum(ages[60:85])

        # print("R",self.result_R_B )

        self.result_Rv_A = np.sum(self.states[self.vaccinated, 0:3, 0]) / np.sum(ages[0:3])
        self.result_Rv_B = np.sum(self.states[self.vaccinated, 3:18, 0]) / np.sum(ages[3:18])
        self.result_Rv_C = np.sum(self.states[self.vaccinated, 18:60, 0]) / np.sum(ages[18:60])
        self.result_Rv_D = np.sum(self.states[self.vaccinated, 60:85, 0]) / np.sum(ages[60:85])

        # print("Rv",self.result_Rv_B )
        # print("\n")

        self.result_R = self.states[self.recovered, :, self.final_t].sum()/sum(self.ages.values()) #sum(self.ages.values())城市总人数
        self.result_Rv = self.states[self.vaccinated, :, 0].sum() / sum(self.ages.values())  # sum(self.ages.values())城市总人数

        #self.time_plot()


    def get_dx_help(self, beta, x,G,result_R,u):
        if (1-x)==0: #说明该年龄层全部人都接种了
            w_x=0
        else:
            w_x=result_R/(1-x)
        dx=x*u*(1-x)*((1-w_x)*math.tanh((beta/2)*(-G))+w_x*math.tanh((beta/2)*(1-G)))
        x=x+dx
        return x

    def get_dx(self,beta,G_A,G_B,G_C,G_D):
        x_A, x_B, x_C, x_D = self.get_x_ABCD()  # 四个年龄段的接种人口占比
        self.x_A=self.get_dx_help(beta, x_A, G_A, self.result_R_A, self.uA)
        self.x_B=self.get_dx_help(beta, x_B, G_B, self.result_R_B, self.uB)
        self.x_C=self.get_dx_help(beta, x_C, G_C, self.result_R_C, self.uC)
        self.x_D=self.get_dx_help(beta, x_D, G_D, self.result_R_D, self.uD)


    def time_plot(self):
        plt.title("Matplotlib demo")
        plt.xlabel("t")
        plt.ylabel("size")
        x=self.final_t
        S = np.sum(self.states[self.susceptible, :, 0:self.final_t], axis=0)/sum(self.ages.values())
        I = np.sum(self.states[self.infected, :, 0:self.final_t], axis=0)/sum(self.ages.values())
        R=np.sum(self.states[self.recovered, :, 0:self.final_t], axis=0)/sum(self.ages.values())
        Rv=np.sum(self.states[self.vaccinated, :, 0], axis=0)/sum(self.ages.values())
        l1 = plt.plot(range(x), S, label='S', linewidth=0.5)
        l2 = plt.plot(range(x), I, label='I', linewidth=0.5)
        l3 = plt.plot(range(x), R, label='R', linewidth=0.5)
        l4 = plt.plot(range(x), np.ones(x)*Rv, label='Rv', linewidth=0.5)
        plt.legend()
        plt.show()

def theory(args):
    num_agebrackets, ages, seed_num, initial_infected_age, lambda1, lambda2, r, M, H, t_max, game_t_max, x0_ABCD_list, G_A, G_B, G_C_list, G_D, beta = args
    R_GD_list = []
    R_A_GD_list = []
    R_B_GD_list = []
    R_C_GD_list = []
    R_D_GD_list = []
    Rv_GD_list = []
    Rv_A_GD_list = []
    Rv_B_GD_list = []
    Rv_C_GD_list = []
    Rv_D_GD_list = []


    for G_C in tqdm(G_C_list):
        game_t = 1
        i = 0
        label = 0
        # 初始化--设置初始接种
        game_theory_model = game_theory(num_agebrackets, t_max, ages, x0_ABCD_list)
        R_t_list = []
        R_A_t_list = []
        R_B_t_list = []
        R_C_t_list = []
        R_D_t_list = []
        Rv_t_list = []
        Rv_A_t_list = []
        Rv_B_t_list = []
        Rv_C_t_list = []
        Rv_D_t_list = []

        # 季节数
        while game_t < game_t_max:
            print(game_t)
            ### 感染传播SIR
            game_theory_model.find_steady_state(seed_num, initial_infected_age, M, H, lambda1, lambda2, r)
            R_A_t_list.append(game_theory_model.result_R_A)
            R_B_t_list.append(game_theory_model.result_R_B)
            R_C_t_list.append(game_theory_model.result_R_C)
            R_D_t_list.append(game_theory_model.result_R_D)
            R_t_list.append(game_theory_model.result_R)
            # print(game_theory_model.result_R,game_theory_model.result_R_A,game_theory_model.result_R_B,game_theory_model.result_R_C,game_theory_model.result_R_D)
            # print("\n")
            Rv_A_t_list.append(game_theory_model.result_Rv_A)
            Rv_B_t_list.append(game_theory_model.result_Rv_B)
            Rv_C_t_list.append(game_theory_model.result_Rv_C)
            Rv_D_t_list.append(game_theory_model.result_Rv_D)
            Rv_t_list.append(game_theory_model.result_Rv)
            # print(game_theory_model.result_Rv, game_theory_model.result_Rv_A, game_theory_model.result_Rv_B,
            #       game_theory_model.result_Rv_C, game_theory_model.result_Rv_D)
            # print("\n")
            ### 判断博弈是否停止
            if game_t > 1:
                if abs(Rv_A_t_list[i] - Rv_A_t_list[i - 1]) < 0.0000001 and abs(Rv_B_t_list[i] - Rv_B_t_list[i - 1]) < 0.000001 and abs(Rv_C_t_list[i] - Rv_C_t_list[i - 1]) < 0.000001 and abs(Rv_D_t_list[i] - Rv_D_t_list[i - 1]) < 0.000001 and label == 0:
                    game_t = game_t_max - 2
                    label = 1

            ### 博弈+疫苗接种
            game_theory_model.get_dx(beta,G_A,G_B,G_C,G_D)
            game_t += 1
            i+=1
        R_GD_list.append(np.array(R_t_list[-2:]).mean())
        R_A_GD_list.append(np.array(R_A_t_list[-2:]).mean())
        R_B_GD_list.append(np.array(R_B_t_list[-2:]).mean())
        R_C_GD_list.append(np.array(R_C_t_list[-2:]).mean())
        R_D_GD_list.append(np.array(R_D_t_list[-2:]).mean())

        Rv_GD_list.append(np.array(Rv_t_list[-2:]).mean())
        Rv_A_GD_list.append(np.array(Rv_A_t_list[-2:]).mean())
        Rv_B_GD_list.append(np.array(Rv_B_t_list[-2:]).mean())
        Rv_C_GD_list.append(np.array(Rv_C_t_list[-2:]).mean())
        Rv_D_GD_list.append(np.array(Rv_D_t_list[-2:]).mean())
    return R_GD_list, R_A_GD_list, R_B_GD_list, R_C_GD_list, R_D_GD_list, Rv_GD_list, Rv_A_GD_list, Rv_B_GD_list, Rv_C_GD_list, Rv_D_GD_list



    # for q in q_list:
    #     R_c_list = []
    #     Rv_c_list = []
    #     for c_percentage in tqdm(c_percentage_list):
    #         game_t = 1
    #         i = 0
    #         label = 0
    #         R_t_list = []
    #         Rv_t_list = []
    #         game_theory_model = game_theory(hyperedges_list, x0_percentage)  # 初始化--设置初始接种
    #         while game_t < game_t_max:  # 季节数
    #             #print("season", game_t)
    #             game_theory_model.find_steady_state(t_max,I0_percentage,g,max_size,r_list)  # 感染传播SIR
    #             R_t_list.append(game_theory_model.r_R)
    #             Rv_t_list.append(game_theory_model.x)
    #             if game_t > 1:
    #                 if abs(Rv_t_list[i] - Rv_t_list[i - 1]) < 0.000001 and label == 0:
    #                     game_t = game_t_max - 10
    #                     label = 1
    #             game_theory_model.get_dx(beta, c_percentage, q)  # 接种疫苗
    #             i = i + 1
    #             game_t += 1
    #         R_c = np.array(R_t_list[-10:]).mean()
    #         Rv_c = np.array(Rv_t_list[-10:]).mean()
    #         if R_c>=1:
    #             R_c=1
    #         if Rv_c>=1:
    #             Rv_c=1
    #         R_c_list.append(R_c)
    #         Rv_c_list.append(Rv_c)
    #     beta_R_c_list.append(R_c_list)
    #     beta_Rv_c_list.append(Rv_c_list)
    # return beta_R_c_list, beta_Rv_c_list #二维数组[q[c,c,c,c],q[c,c,c,c,c]]






