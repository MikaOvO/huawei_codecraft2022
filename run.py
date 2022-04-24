import os
from unittest import result

data_dir = os.getcwd() + '\\data'
output_dir = os.getcwd() + '\\output'

main_file = os.getcwd() + '\SDK_C++\CodeCraft-2022\src\CodeCraft-2022.cpp'

os.system('g++ -o codecraft.exe %s' % (main_file))

for path, dir_list, _ in os.walk(data_dir):
    for dir in dir_list:
        ## only run sample
        ##if dir != '30_10_50_c20000':
        ##    continue
        data_child_dir = os.path.join(path, dir)
        solution_dir = output_dir + '\\' + dir
        if not os.path.exists(solution_dir):
            os.makedirs(solution_dir)
        ## modify here
        
        result_file = output_dir + '\\' +  'result.txt'
        info = 'x!' ## 控制是否打印到result里,debug的时候设置为! 否则设置为任意其他字符
        params = ' info=%s data_dir=%s debug_file=%s output_dir=%s result_file=%s' % \
                 (info, data_child_dir, solution_dir + '\\' + 'debug.txt', solution_dir ,  result_file)
        ##
        
        os.system('codecraft.exe%s' % (params))