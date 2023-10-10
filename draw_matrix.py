import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from PIL import Image, ImageDraw


MAX_LEN = 70 # max lenght of microtubule
PROTOFIL = 4 # number of protofilaments
SCREEN_WIDTH = 1200
SCREEN_HEIGHT = 800


file = open("matrix_for_animation_1", "r")
file_time = open("len_time_graph_1", "r")
times = [0.0]
for line in file_time:
    l, t = map(float, line.split())
    times.append(t)
duration = [int((times[i + 1] - times[i]) * 100) for i in range(len(times) - 1)]
matrices = file.read().split('\n')
file.close()
file_time.close()

images = []
for i in range(1000):
    im = Image.new('RGB', (SCREEN_WIDTH, SCREEN_HEIGHT), "white")
    draw = ImageDraw.Draw(im)

    arr = np.array(list(map(int, matrices[i].split())))
    arr = np.reshape(arr, (PROTOFIL, MAX_LEN))

    for p in range(PROTOFIL):
        for j in range(MAX_LEN):
            x = j * 11 + 5
            y = p * 11 + 395
            r = 4
            if arr[p][j] == 0:
                # print("kek")
                draw.ellipse((x - r, y - r, x + r, y + r), fill="white")

            elif arr[p][j] == 1: # straight, gdp
                draw.ellipse((x - r, y - r, x + r, y + r), fill="yellowgreen")

            elif arr[p][j] == 2: # curved, gdp
                draw.ellipse((x - r, y - r, x + r, y + r), fill="darkgreen")

            elif arr[p][j] == 3: # straight, gtp
                draw.ellipse((x - r, y - r, x + r, y + r), fill="red")

            elif arr[p][j] == 4: # curved, gtp
                draw.ellipse((x - r, y - r, x + r, y + r), fill="darkred")

            elif arr[p][j] == 5: # special
                draw.ellipse((x - r, y - r, x + r, y + r), fill="lightpink")
    images.append(im)
print(len(images), len(duration))

images[0].save(
    'animation_microtub_1.gif',
    save_all=True,
    append_images=images[1:],
    optimize=True,
    duration=1,
    loop=0
)