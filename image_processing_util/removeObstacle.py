def removeObstacle(I_list):
    I_obstacle = I_list[0] > 0
    for i in range(1, len(I_list)):
        I_obstacle = I_obstacle * (I_list[i] > 0)
    for i in range(len(I_list)):
        I_list[i][I_obstacle > 0] = 0
    return I_list