import os

from pytest import param

data_dir = os.getcwd() + '\\data'
output_dir = os.getcwd() + '\\output'

main_file = os.getcwd() + '\SDK_C++\CodeCraft-2022\src\CodeCraft-2022.cpp'

for path, dir_list, _ in os.walk(data_dir):
    for dir in dir_list:
        ## only run sample
        ## if dir != 'sample':
        ##    continue
        data_child_dir = os.path.join(path, dir)
        
        ## modify here
        params = ' data_dir=%s debug_file=%s output_dir=%s' % \
                 (data_child_dir, output_dir + '\\' + dir + '_debug.txt', output_dir)
        ##
        
        os.system('g++ -o codecraft.exe %s' % (main_file))
        os.system('codecraft.exe%s' % (params))