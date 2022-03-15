import os

from pytest import param

data_dir = os.getcwd() + '\\data'
output_dir = os.getcwd() + '\\output'

main_file = os.getcwd() + '\SDK_C++\CodeCraft-2022\src\CodeCraft-2022.cpp'

for path, dir_list, _ in os.walk(data_dir):
    for dir in dir_list:
        ## only run sample
        ## if dir != 'sample':
        ##   continue
        data_child_dir = os.path.join(path, dir)
        solution_dir = output_dir + '\\' + dir
        if not os.path.exists(solution_dir):
            os.makedirs(solution_dir)
        ## modify here
        params = ' data_dir=%s debug_file=%s output_dir=%s result_file=%s' % \
                 (data_child_dir, output_dir + '\\' + dir + '_debug.txt', solution_dir ,  output_dir + '\\' + dir + '_result.txt')
        ##
        
        os.system('g++ -o codecraft.exe %s' % (main_file))
        os.system('codecraft.exe%s' % (params))