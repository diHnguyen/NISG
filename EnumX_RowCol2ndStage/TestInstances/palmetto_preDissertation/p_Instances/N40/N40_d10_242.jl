global arcs = [1 4; 1 7; 1 9; 1 11; 1 12; 1 26; 1 31; 1 33; 2 9; 2 32; 2 34; 2 38; 3 14; 3 17; 3 18; 3 30; 3 39; 4 23; 4 24; 4 31; 5 7; 5 16; 5 17; 5 25; 5 31; 5 40; 6 13; 6 17; 6 19; 6 21; 6 30; 6 40; 7 15; 7 17; 7 28; 7 33; 8 10; 8 15; 8 27; 8 30; 8 36; 8 38; 9 18; 9 19; 10 3; 10 6; 10 14; 11 20; 11 23; 12 22; 12 33; 12 38; 13 37; 13 40; 14 22; 14 32; 15 11; 15 17; 15 21; 15 28; 15 39; 15 40; 16 8; 16 15; 16 23; 16 25; 16 35; 17 8; 17 27; 17 38; 18 4; 18 6; 19 2; 19 8; 19 10; 19 14; 19 26; 19 28; 20 4; 20 7; 20 8; 20 9; 20 19; 20 32; 21 10; 21 19; 21 26; 22 7; 22 17; 23 6; 23 21; 23 22; 23 24; 23 31; 24 7; 24 10; 24 23; 24 32; 24 35; 25 3; 25 5; 25 38; 26 13; 26 19; 26 24; 26 29; 27 3; 27 16; 27 20; 27 38; 28 37; 29 7; 29 9; 30 26; 30 31; 30 34; 30 39; 31 4; 31 8; 31 30; 31 34; 32 7; 32 9; 32 12; 32 15; 32 19; 32 27; 32 30; 32 37; 33 8; 33 13; 33 14; 33 22; 33 24; 33 35; 34 5; 34 6; 34 9; 34 36; 35 13; 35 28; 35 34; 36 18; 37 16; 37 22; 37 30; 38 20; 38 40; 39 14; 39 20; 39 33]
global d_x = [3.0, 10.0, 7.0, 6.0, 3.0, 6.0, 7.0, 9.0, 6.0, 8.0, 8.0, 9.0, 8.0, 7.0, 6.0, 3.0, 9.0, 5.0, 3.0, 10.0, 1.0, 10.0, 4.0, 6.0, 2.0, 5.0, 8.0, 7.0, 6.0, 3.0, 4.0, 6.0, 1.0, 9.0, 1.0, 4.0, 8.0, 7.0, 6.0, 7.0, 7.0, 3.0, 5.0, 3.0, 5.0, 7.0, 5.0, 6.0, 6.0, 7.0, 1.0, 6.0, 7.0, 7.0, 3.0, 10.0, 9.0, 9.0, 10.0, 6.0, 3.0, 3.0, 10.0, 5.0, 2.0, 6.0, 5.0, 1.0, 7.0, 4.0, 4.0, 8.0, 1.0, 2.0, 6.0, 6.0, 10.0, 6.0, 5.0, 6.0, 7.0, 8.0, 9.0, 6.0, 8.0, 5.0, 2.0, 9.0, 2.0, 1.0, 10.0, 2.0, 3.0, 7.0, 7.0, 4.0, 9.0, 5.0, 1.0, 2.0, 6.0, 9.0, 1.0, 6.0, 5.0, 7.0, 4.0, 4.0, 1.0, 8.0, 7.0, 2.0, 4.0, 10.0, 1.0, 9.0, 4.0, 10.0, 2.0, 5.0, 8.0, 3.0, 1.0, 4.0, 4.0, 5.0, 10.0, 1.0, 1.0, 8.0, 4.0, 1.0, 7.0, 1.0, 7.0, 7.0, 9.0, 10.0, 5.0, 9.0, 1.0, 1.0, 1.0, 5.0, 3.0, 10.0, 9.0, 8.0, 9.0, 10.0, 1.0]
global b_x = 5
global d_y = [3.0, 1.0, 3.0, 1.0, 7.0, 8.0, 10.0, 8.0, 9.0, 9.0, 5.0, 3.0, 1.0, 8.0, 2.0, 3.0, 5.0, 8.0, 2.0, 4.0, 1.0, 7.0, 8.0, 7.0, 8.0, 1.0, 10.0, 4.0, 3.0, 8.0, 10.0, 5.0, 1.0, 7.0, 3.0, 10.0, 4.0, 6.0, 7.0, 1.0, 8.0, 8.0, 3.0, 6.0, 10.0, 1.0, 9.0, 9.0, 2.0, 4.0, 4.0, 5.0, 1.0, 6.0, 10.0, 9.0, 2.0, 1.0, 6.0, 6.0, 9.0, 4.0, 4.0, 8.0, 5.0, 8.0, 10.0, 5.0, 8.0, 3.0, 7.0, 10.0, 10.0, 9.0, 3.0, 3.0, 8.0, 3.0, 10.0, 3.0, 3.0, 10.0, 7.0, 1.0, 10.0, 7.0, 8.0, 9.0, 4.0, 5.0, 8.0, 4.0, 10.0, 4.0, 6.0, 4.0, 6.0, 8.0, 5.0, 7.0, 1.0, 4.0, 1.0, 8.0, 10.0, 9.0, 10.0, 3.0, 9.0, 3.0, 7.0, 9.0, 7.0, 9.0, 2.0, 5.0, 3.0, 2.0, 4.0, 3.0, 2.0, 1.0, 5.0, 5.0, 5.0, 8.0, 4.0, 9.0, 8.0, 5.0, 6.0, 3.0, 4.0, 8.0, 2.0, 6.0, 10.0, 2.0, 4.0, 9.0, 6.0, 1.0, 8.0, 8.0, 8.0, 3.0, 2.0, 1.0, 1.0, 1.0, 2.0]
global b_y = 10
global p = [0.386, 0.136, 0.552, 0.022, 0.105, 0.27, 0.777, 0.506, 0.82, 0.455, 0.79, 0.848, 0.989, 0.733, 0.439, 0.256, 0.311, 0.818, 0.53, 0.252, 0.921, 0.127, 0.568, 0.485, 0.835, 0.346, 0.207, 0.345, 0.922, 0.755, 0.518, 0.98, 0.122, 0.498, 0.2, 0.488, 0.283, 0.642, 0.102, 0.838, 0.244, 0.374, 0.149, 0.13, 0.625, 0.311, 0.653, 0.217, 0.212, 0.117, 0.461, 0.211, 0.415, 0.985, 0.634, 0.114, 0.347, 0.301, 0.723, 0.352, 0.453, 0.633, 0.728, 0.701, 0.055, 0.135, 0.538, 0.777, 0.875, 0.841, 0.275, 0.341, 0.607, 0.874, 0.641, 0.207, 0.369, 0.994, 0.169, 0.354, 0.369, 0.509, 0.811, 0.305, 0.286, 0.959, 0.304, 0.236, 0.453, 0.76, 0.908, 0.575, 0.301, 0.496, 0.126, 0.279, 0.047, 0.078, 0.327, 0.203, 0.024, 0.189, 0.531, 0.936, 0.904, 0.438, 0.029, 0.465, 0.777, 0.667, 0.259, 0.515, 0.797, 0.35, 0.078, 0.992, 0.496, 0.133, 0.444, 0.834, 0.978, 0.004, 0.477, 0.851, 0.794, 0.213, 0.827, 0.588, 0.574, 0.722, 0.008, 0.742, 0.346, 0.455, 0.291, 0.075, 0.286, 0.612, 0.543, 0.528, 0.312, 0.933, 0.518, 0.969, 0.236, 0.877, 0.055, 0.617, 0.709, 0.049, 0.133]
global q = [0.598, 0.629, 0.848, 0.251, 0.351, 0.619, 0.941, 0.876, 0.936, 0.686, 0.915, 0.947, 0.993, 0.787, 0.562, 0.281, 0.89, 0.829, 0.945, 0.756, 0.998, 0.674, 0.639, 0.874, 0.98, 0.652, 0.7, 0.674, 0.99, 0.855, 0.581, 0.996, 0.494, 0.507, 0.822, 0.714, 0.847, 0.713, 0.758, 0.881, 0.436, 0.479, 0.617, 0.29, 0.926, 0.976, 0.66, 0.683, 0.87, 0.166, 0.632, 0.367, 0.704, 0.985, 0.788, 0.447, 0.933, 0.613, 0.909, 0.437, 0.738, 0.87, 0.823, 0.987, 0.415, 0.57, 0.85, 0.956, 0.899, 0.873, 0.624, 0.723, 0.658, 0.937, 0.689, 0.968, 0.858, 0.997, 0.258, 0.869, 0.633, 0.726, 0.982, 0.694, 0.786, 0.983, 0.965, 0.353, 0.739, 0.822, 0.964, 0.968, 0.842, 0.572, 0.345, 0.318, 0.866, 0.775, 0.542, 0.509, 0.915, 0.456, 0.888, 0.964, 0.936, 0.565, 0.451, 0.517, 0.927, 0.889, 0.53, 0.675, 0.821, 0.952, 0.226, 0.996, 0.921, 0.531, 0.95, 0.938, 0.987, 0.948, 0.989, 0.987, 0.882, 0.467, 0.895, 0.675, 0.685, 0.994, 0.043, 0.896, 0.494, 0.986, 0.43, 0.111, 0.795, 0.67, 0.609, 0.772, 0.379, 0.996, 0.748, 0.972, 0.395, 0.991, 0.284, 0.704, 0.8, 0.723, 0.434]
global origin = 1
global destination = 40