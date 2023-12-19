######结果输出
def print_file(file_txt,list):
    with open(file_txt,'a') as f:
        for i in list:
            for j in i:
                f.write("%f" %j)
                f.write(' ')
            f.write('\n')
        f.write('\n')
        f.close()

def print_file_1(file_txt,list):
    with open(file_txt,'a') as f:
        for i in list:
            f.write("%f" %i)
            f.write(' ')
        f.write('\n')
        f.close()

def print_file_p(file_txt):
    with open(file_txt,'a') as f:
        f.write('\n\n\n')
        f.close()
