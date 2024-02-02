# bMT_mc4statemodel

## Описание

Программа для моделирования работы микротрубочки. Основана на алгоритме Гиллеспи. Каждый тубулин имеет 2 конформации: curved, straight; 2 химических состояния: T, D (гидролизованный и негидролизованный). 

## Структура проекта

* `gillespie.py` - основной файл программы, алгоритм Гиллеспи
* 'draw_graphs.py' - программа, рисующая график зависимости длины микротрубочки от времени
* 'draw_curved.py' - программа, рисующая график зависимости длины изогнутного конца микротрубочки от времени
* 'draw_matrix.py' - программа, рисующая микротрубочку
* 'get_velocities_from_data.py' - программа, по графику определяющая скорости роста микротрубочки
