import numpy as np
import sys


# constants
# on, off, hydrolysis, bending for gtp, bending for gdp, straightening for gtp, straightening for gdp
K_H = 0.1 # GTP hydrolysis rate constant, 1/s
K_PLUS = 0.6 # 1/(um*s), association const
C_TUB = 10.0 # BtubAB concentration, uM
G_LONG_D = -12.8
G_LONG_T = -13.2
K_OFF_GTP = K_PLUS * (10 ** 6) / np.exp(-G_LONG_T) # off const for gtp, 1/s
K_OFF_GDP = K_PLUS * (10 ** 6) / np.exp(-G_LONG_D) # off const of gdp, 1/s

K_STR_T = 300.0 # GTP-tubulin dimer straightening rate constant with zero or one protofilament neighbor, 1/s
K_STR_D = 50.0 # GDP-tubulin dimer straightening rate constant with zero or one protofilament neighbor, 1/s
G_BEND_T = 6.0 # deformation energy of a straight tubulin dimer, connected to a given GTP-tubulin, kT
G_BEND_D = 8.0 # deformation energy of a straight tubulin dimer, connected to a given GDP-tubulin, kT
G_LAT_T = -7.5 # free energy of a lateral bond that a GTP tubulin forms with an adjacent dimer (T/D + T), kT
G_LAT_D = -5.6 # free energy of a lateral bond that a GDP tubulin forms with an adjacent dimer (D+D), kT
LAMBD = 1.0 # barrier coef

MAX_IT = 500 * 900 # max number of iterations
NPASS = 800
MAX_LEN = 70 # max lenght of microtubule
ARR_SIZE = 7000
PROTOFIL = 4 # number of protofilaments
INIT_LENGHT = 4 # length of microtubule in the beginning of simulation
INF = 10 ** 9


def print_microtubule(lenn = 50): # the function prints part of microtubule of lenght = lenn
    arr = [[0 for i in range(lenn)] for j in range(PROTOFIL)]
    for p in range(PROTOFIL):
        for j in range(lenn):
            if h[i][j] == 2 and c[i][j] == 2:
                arr[i][j] = 4
            elif h[i][j] == 2 and c[i][j] == 1:
                arr[i][j] = 3
            elif h[i][j] == 1 and c[i][j] == 2:
                arr[i][j] = 2
            elif h[i][j] == 1 and c[i][j] == 1:
                arr[i][j] = 1
            elif h[i][j] == 3:
                arr[i][j] = 5
    for p in range(PROTOFIL):
        print(*h[p][:lenn])
    print("__")
    for p in range(PROTOFIL):
        print(*c[p][:lenn])


def save_for_graphs_len_time(a): # saves to a file average length of a microtubule, time
    file = open("len_time_graph_1", "w")
    s = ""
    for i in range(len(a)):
        s += str(a[i][0]) + " " + str(a[i][1]) + "\n"
    file.write(s)
    file.close()


def each_cur_length():
    for p in range(PROTOFIL):
        if first_to_str[p] == -1:
            len_each[p].append(last_added[p] + 1)
        else:
            len_each[p].append(first_to_str[p])


def save_len_each_time():
    file = open("len_each_time_graph_1", "w")
    s = ""
    for i in range(len(len_each[0])):
        s += str(len_each[0][i]) + " " + str(len_each[1][i]) + " " + str(len_each[2][i]) + " " + str(len_each[3][i]) + " " + str(len_time[i][1]) + "\n"
    file.write(s)
    file.close()


def save_len_curved():
    file = open("len_curved_1", "w")
    s = ""
    for i in range(len(len_curved)):
        s += str(len_curved[i]) + "\n"
    file.write(s)
    file.close()


def save_matrix(h, it): # saves information about condition of tubulins in microtubule
    if it == 1:
        file = open("matrix_for_animation_1", "w") # cleans everything that was in the file
    else:
        file = open("matrix_for_animation_1", "a") # addes text to the end of the file
    s = ""
    lst = max(last_added)
    for i in range(PROTOFIL): # initialization
        new_h = h[i][max(0, lst - MAX_LEN + 1):max(0, lst - MAX_LEN + 1) + MAX_LEN] # if last_added[i] > MAX_LEN, we need to take only the last MAX_LEN elements
        new_c = c[i][max(0, lst - MAX_LEN + 1):max(0, lst - MAX_LEN + 1) + MAX_LEN] # if last_added[i] > MAX_LEN, we need to take only the last MAX_LEN elements
        arr = [0 for i in range(len(new_h))]
        for i in range(len(new_h)):
            if new_h[i] == 2 and new_c[i] == 2:
                arr[i] = 4
            elif new_h[i] == 2 and new_c[i] == 1:
                arr[i] = 3
            elif new_h[i] == 1 and new_c[i] == 2:
                arr[i] = 2
            elif new_h[i] == 1 and new_c[i] == 1:
                arr[i] = 1
            elif new_h[i] == 3:
                arr[i] = 5
        s += " ".join(map(str, arr)) + " "
    s += "\n"
    file.write(s)
    file.close()

def cnt_cur_length():
    lens = []
    curved = []
    for i in range(PROTOFIL):
        if first_to_str[i] == -1:
            lens.append(last_added[i] + 1)
            curved.append(0)
        else:
            lens.append(first_to_str[i])
            curved.append(last_added[i] + 1 - first_to_str[i])
    len_now = np.median(lens)
    std_now = np.std(lens)
    std_len_str.append(std_now)
    len_time.append((len_now, cur_time))
    len_curved.append(np.mean(curved))


def check_error_arrays_hc(): # checks whether the arrays work together
    for p in range(PROTOFIL):
        for j in range(1, len(h[0])):
            if c[p][j - 1] == 0 and c[p][j] != 0:
                print("ERROR C", p, j, it)
            if h[p][j - 1] == 0 and h[p][j] != 0:
                print("ERROR H", p, j, it)
            if h[p][j - 1] == 0 and c[p][j - 1] != 0:
                print("BOTH ERROR1", it)
            if h[p][j - 1] != 0 and c[p][j - 1] == 0:
                print("BOTH ERROR2", p, j, it)


def check_error_2(): # checks whether the array first_to_str works correctly
    for p in range(PROTOFIL):
        if first_to_str[p] != -1 and (not 2 in c[p]):
            print("oh no first_to_str doesn't work right", it, p)


def check_hydr_error(): # checks whether the array can_be_hydr works correctly
    for p, j in can_be_hydr:
        if h[p][j] != 2 and h[p][j + 1] != 0:
            print(it, "HYDR ERR", p, j)


def check_set():
    ct = set(curved_t)
    cd = set(curved_d)
    cbh = set(can_be_hydr)
    ct = list(ct)
    if sorted(ct) != sorted(curved_t):
        print("Oh no ct!=curved t")
    cd = list(cd)
    if sorted(cd) != sorted(curved_d):
        print("Oh no cd!=curved d")
    cbh = list(cbh)
    if sorted(cbh) != sorted(can_be_hydr):
        print("Oh no cbh!=can_be_hydr")
        print(cbh)
        print("____")
        print(can_be_hydr)
        print("=========")


def check_errors():
    check_error_arrays_hc()
    check_error_2()
    check_hydr_error()


def create_microtub(lenn=20): # creates microtubule of straight gtp tubulins of length = lenn
    for p in range(PROTOFIL):
        first_to_str[p] = INIT_LENGHT
        for j in range(INIT_LENGHT, INIT_LENGHT + lenn):
            h[p][j] = 2
            c[p][j] = 2
            last_added[p] = j
            curved_t.append((p, j))


def get_lateral_energy(a, b):
    if h[a[0]][a[1]] == 1 and h[b[0]][b[1]] == 1:
        return G_LAT_D
    return G_LAT_T

print(sys.argv[1])
seed_num = int(sys.argv[1])
np.random.seed(seed_num)

# initialisation

c = [[0 for i in range(ARR_SIZE)] for j in range(PROTOFIL)] # array that describes conformation of tubulins
# 0 -- no tubulin
# 1 -- straight
# 2 -- curved

h = [[0 for i in range(ARR_SIZE)] for j in range(PROTOFIL)] # array that describes gtp/gdp form of tubulins
# 0 -- no tubulin
# 1 -- gdp
# 2 -- gtp
# 3 -- special type of tubulins, they cannot be destroyed

for i in range(PROTOFIL): # initialization
    for j in range(INIT_LENGHT):
        c[i][j] = 1 # straight
        h[i][j] = 3 # cannot be destroyed

last_added = [INIT_LENGHT - 1 for i in range(4)] # index of last tubulins in each protofilament
first_to_str = [-1 for i in range(4)] # last straight in each protofilament
curved_t = [] # array of pairs of indices of curved gtp tubulins
curved_d = [] # array of pairs of indices of curved gdp tubulins
can_be_hydr = [] # array of pairs of indices of tubulins that can be hydrolyzed
len_time = [] # array that contains pairs of (current average microtubule length, current time)
control_t = [0 for i in range(6)] # array of frequency of every event
len_each = []
std_len_str = []
len_curved = []
for i in range(PROTOFIL):
    len_each.append([])

# simulation

cur_time = 0.0
it = 0
while it < MAX_IT:
    it += 1
    ind_on, t_on, ind_off_t, t_off_t, ind_off_d, t_off_d, ind_hydr, t_hydr, ind_str, t_str, ind_bend, t_bend = INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF

    # on
    k_on = C_TUB * K_PLUS
    rand = np.random.rand(PROTOFIL) # PROTOFIL рандомных чисел
    ind_on = np.argmin(-np.log(rand) / k_on)
    t_on = -np.log(rand[ind_on]) / k_on

    # off
    cnt_off_gtp = len(curved_t)
    if cnt_off_gtp > 0:
        rand = np.random.rand(cnt_off_gtp) # cnt_off рандомных чисел
        ind_off_t = np.argmin(-np.log(rand) / K_OFF_GTP)
        t_off_t = -np.log(rand[ind_off_t]) / K_OFF_GTP

    cnt_off_gdp = len(curved_d)
    if cnt_off_gdp > 0:
        rand = np.random.rand(cnt_off_gdp) # cnt_off рандомных чисел
        ind_off_d = np.argmin(-np.log(rand) / K_OFF_GDP)
        t_off_d = -np.log(rand[ind_off_d]) / K_OFF_GDP

    # straghtening 
    t_str = INF
    for p in range(PROTOFIL):
        if first_to_str[p] == -1: # no curved in this protofilament
            continue

        if h[p][first_to_str[p]] == 2: # gtp
            rand = np.random.rand(1)
            t_curr = float(-np.log(rand) / K_STR_T)
            if t_curr < t_str:
                t_str = t_curr
                ind_str = p
        elif h[p][first_to_str[p]] == 1:
            rand = np.random.rand(1)
            t_curr = float(-np.log(rand) / K_STR_D)
            if t_curr < t_str:
                t_str = t_curr
                ind_str = p

    # bending
    k_bend = INF
    for p in range(PROTOFIL):
        if h[p][last_added[p]] == 3: # or h[p][last_added[p]] == 0:
            # print("-")
            continue
        to_bend = -1
        if first_to_str[p] == -1: # no curved in this protofilament
            to_bend = last_added[p]
        elif first_to_str[p] - 1 >= 0 and h[p][first_to_str[p] - 1] != 3:
            to_bend = first_to_str[p] - 1
        if to_bend == -1:
            # print("-")
            continue
        
        energy = 0.0
        if p == 0: # seam
            if c[PROTOFIL - 1][to_bend - 2] == 1:
                energy += 0.5 * get_lateral_energy((0, to_bend), (PROTOFIL - 1, to_bend - 2))
            if c[PROTOFIL - 1][to_bend - 1] == 1:
                energy += 0.5 * get_lateral_energy((0, to_bend), (PROTOFIL - 1, to_bend - 1))
            if c[1][to_bend] == 1:
                energy += get_lateral_energy((0, to_bend), (1, to_bend))
        elif p == PROTOFIL - 1: # seam
            if c[0][to_bend + 1] == 1:
                energy += 0.5 * get_lateral_energy((0, to_bend + 1), (PROTOFIL - 1, to_bend))
            if c[0][to_bend + 2] == 1:
                energy += 0.5 * get_lateral_energy((0, to_bend + 2), (PROTOFIL - 1, to_bend))
            if c[PROTOFIL - 2][to_bend] == 1:
                energy += get_lateral_energy((PROTOFIL - 1, to_bend), (PROTOFIL - 2, to_bend))
        else: # not seam
            if c[(p + 1) % PROTOFIL][to_bend] == 1:
                energy += get_lateral_energy(((p + 1) % PROTOFIL, to_bend), (p, to_bend))
            if c[(p - 1) % PROTOFIL][to_bend] == 1:
                energy += get_lateral_energy(((p - 1) % PROTOFIL, to_bend), (p, to_bend))

        if h[p][to_bend] == 2:
            energy += G_BEND_T
            k_bend = K_STR_T * np.exp(energy)
        elif h[p][to_bend] == 1:
            energy += G_BEND_D
            k_bend = K_STR_D * np.exp(energy)

        rand = np.random.rand(1)
        t_curr = float(-np.log(rand) / k_bend)
        if t_curr < t_bend:
            t_bend = t_curr
            ind_bend = p

    cnt_hydr = len(can_be_hydr)
    if cnt_hydr > 0:
        rand = np.random.rand(cnt_hydr) # cnt_hydr рандомных чисел
        ind_hydr = np.argmin(-np.log(rand) / K_H)
        t_hydr = -np.log(rand[ind_hydr]) / K_H
        

    # EVENTS ________________________________________________________________________________________________________________________________________________________________
    t_all = [t_on, t_off_t, t_off_d, t_hydr, t_str, t_bend]
    t_min_ind = t_all.index(min(t_all))
    t_min = -1

    if t_min_ind == 0: # event -- on
        t_min = t_on

        last_added[ind_on] += 1
        h[ind_on][last_added[ind_on]] = 2 # new gtp tubulin added
        c[ind_on][last_added[ind_on]] = 2 # new tubulin is curved
        curved_t.append((ind_on, last_added[ind_on]))
        if first_to_str[ind_on] == -1:
            first_to_str[ind_on] = last_added[ind_on]

    elif t_min_ind == 1 or t_min_ind == 2: # event -- off

        p, ind_in_protofil = -1, -1
        if t_min_ind == 1:
            t_min = t_off_t
            p, ind_in_protofil = curved_t[ind_off_t]
        else:
            t_min = t_off_d
            p, ind_in_protofil = curved_d[ind_off_d]

        if first_to_str[p] == ind_in_protofil:
            first_to_str[p] = -1

        for j in range(ind_in_protofil, last_added[p] + 1): # remove everything above this tubulin
            if h[p][j] == 2:
                curved_t.remove((p, j))
                if (p, j) in can_be_hydr:
                    can_be_hydr.remove((p, j))
            elif h[p][j] == 1:
                curved_d.remove((p, j))

            h[p][j] = 0 # tubulin removed
            c[p][j] = 0

        last_added[p] = ind_in_protofil - 1

        if h[p][last_added[p]] == 2 and (p, last_added[p]) in can_be_hydr:
            can_be_hydr.remove((p, last_added[p]))

    elif t_min_ind == 3: # event -- hydrolysis
        t_min = t_hydr
        p, j = can_be_hydr[ind_hydr]
        can_be_hydr.remove((p, j)) # removing from the list of hydr.
        h[p][j] = 1 # gtp -> gdp form
        if (p, j) in curved_t:
            curved_t.remove((p, j))
            curved_d.append((p, j))
    
    elif t_min_ind == 4: # event -- straightening
        t_min = t_str
        c[ind_str][first_to_str[ind_str]] = 1

        if h[ind_str][first_to_str[ind_str]] == 2:
            curved_t.remove((ind_str, first_to_str[ind_str]))
        elif h[ind_str][first_to_str[ind_str]] == 1:
            curved_d.remove((ind_str, first_to_str[ind_str]))

        if h[ind_str][first_to_str[ind_str] - 1] == 2 and h[ind_str][first_to_str[ind_str]] == 2:
            can_be_hydr.append((ind_str, first_to_str[ind_str] - 1))

        if last_added[ind_str] >= first_to_str[ind_str] + 1:
            first_to_str[ind_str] += 1
        else:
            first_to_str[ind_str] = -1

    elif t_min_ind == 5: # event -- bending
        
        t_min = t_bend

        if first_to_str[ind_bend] == -1: # no curved in this protofilament
            to_bend = last_added[ind_bend]
        elif first_to_str[ind_bend] - 1 >= 0 and h[ind_bend][first_to_str[ind_bend] - 1] != 3:
            to_bend = first_to_str[ind_bend] - 1

        c[ind_bend][to_bend] = 2

        if h[ind_bend][to_bend] == 2:
            curved_t.append((ind_bend, to_bend))
        elif h[ind_bend][to_bend] == 1:
            curved_d.append((ind_bend, to_bend))

        if (ind_bend, to_bend - 1) in can_be_hydr:
            can_be_hydr.remove((ind_bend, to_bend - 1))

        first_to_str[ind_bend] = to_bend

    cur_time += t_min

    if it % NPASS == 0 or it == 1:
        cnt_cur_length()
        each_cur_length()
    if it % NPASS == 0 or it == 1:
        save_matrix(h, it)

    if it % 100000 == 0:
        print(it)
    control_t[t_min_ind] += 1
    
print(cur_time, "cur_time")

save_for_graphs_len_time(len_time)
save_len_each_time()
save_len_curved()
print(control_t)
percent_each_event = [0 for i in range(len(control_t))]
for i in range(len(percent_each_event)):
    percent_each_event[i] = control_t[i] / sum(control_t) * 100
print(percent_each_event)


velocities = []
for i in range(1, len(len_time)):
    velocities.append(((len_time[i][0] - INIT_LENGHT) / len_time[i][1]) * 8) # len_now / cur_time * 8

# print("Average velocity: ", np.mean(velocities))
vel = (len_time[-1][0] - INIT_LENGHT) / len_time[-1][1] * 8
# print(len_time[-10:-1])
print(vel, "velocity")

file = open("velocity_info.txt", "w")
file.write("Average velocity: " + str(np.mean(velocities)) + "\n")
file.write(" ".join(map(str, velocities)))
file.close()



