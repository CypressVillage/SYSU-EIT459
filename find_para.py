import os
import itertools

patch_data = {
    'scStr.modulation.subcarrierSpacing': [15e3, 30e3, 60e3],
    'scStr.modulation.nGuardSymbols': [1, 7],
    'scStr.modulation.samplingRate': [15e3*2048, 30e3*2048, 60e3*2048]
}

def Patch_LTEA(patchPack):
    path = os.getcwd() + '\\Scenarios\\LTEAcompliant.m'
    with open(path, 'r') as file:
        data = file.readlines()
    for key, value in patchPack.items():
        for i in range(len(data)):
            if key in data[i]:
                data[i] = key + ' = ' + str(value) + ';\n'
                print('Patched: ' + key + ' = ' + str(value))
    with open(path, 'w') as file:
        file.writelines(data)

# 生成所有组合
keys = patch_data.keys()
values = patch_data.values()

combinations = list(itertools.product(*values))

# 将组合转换为字典
result_dicts = [dict(zip(keys, combination)) for combination in combinations]

# 打印结果
for result in result_dicts:
    Patch_LTEA(result)
    os.system('matlab -batch "main" > NUL')
    print('-------------------')