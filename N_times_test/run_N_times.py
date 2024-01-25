import os
import numpy as np


substep_number = 6

sweep_parameter = {
    "name": "t_1",
    "expt_number": 32,
    "description": f"new substep {substep_number} simulation for different t",
    "parameter_list": [60, 70, 80, 90, 100, 110, 120]
}



N = 100

# for _ in range(0, N):
#     os.system(f"python reorder_vector_and_new_protocol.py")
# # for _ in range(0, N):
# #     

print( sweep_parameter["parameter_list"])
for index, t in enumerate(sweep_parameter["parameter_list"]):
    print('\x1b[6;30;42m' + f"{index}/{len(sweep_parameter['parameter_list'])}" + '\x1b[0m')
    print("current parameter is: ", t)
    os.system(f"python substep_test_script.py {t} \"{sweep_parameter['description']}\" {substep_number}")




