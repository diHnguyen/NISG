global arcs = [1 2; 1 14; 1 16; 1 17; 1 22; 1 36; 2 4; 2 8; 2 14; 2 29; 2 31; 3 4; 3 19; 3 30; 4 16; 4 23; 5 22; 5 25; 5 27; 5 30; 6 11; 6 14; 6 20; 6 25; 6 33; 6 35; 7 16; 7 34; 7 35; 7 36; 8 16; 8 30; 8 39; 9 15; 10 5; 10 34; 10 35; 10 38; 11 3; 11 12; 11 27; 11 34; 12 2; 12 15; 12 33; 13 10; 13 24; 13 25; 13 34; 13 37; 13 38; 14 17; 14 18; 15 4; 15 11; 15 14; 15 29; 16 10; 16 17; 16 21; 16 25; 16 28; 16 34; 17 4; 17 10; 17 13; 17 18; 17 22; 17 39; 18 12; 18 19; 18 23; 18 25; 18 40; 19 3; 19 34; 19 35; 20 10; 21 6; 21 18; 21 24; 21 33; 21 35; 22 2; 22 8; 22 21; 22 26; 22 38; 23 5; 23 18; 23 22; 23 30; 23 31; 23 40; 24 13; 24 27; 25 9; 25 10; 25 23; 26 17; 26 25; 26 32; 26 37; 27 6; 27 30; 27 31; 27 37; 27 40; 28 3; 28 9; 28 22; 28 36; 29 7; 29 16; 29 36; 29 37; 30 13; 30 23; 30 34; 31 36; 32 22; 32 31; 33 3; 33 19; 33 25; 33 37; 34 2; 34 19; 34 27; 34 29; 35 12; 35 13; 35 16; 35 17; 35 27; 35 30; 35 31; 35 40; 36 5; 36 7; 36 11; 36 14; 36 22; 36 24; 36 30; 36 40; 37 14; 37 27; 37 28; 38 10; 38 26; 38 40; 39 22]
global d_x = [8.0, 1.0, 10.0, 5.0, 10.0, 3.0, 9.0, 1.0, 1.0, 5.0, 2.0, 8.0, 8.0, 1.0, 4.0, 7.0, 10.0, 6.0, 10.0, 1.0, 4.0, 8.0, 9.0, 10.0, 6.0, 9.0, 2.0, 7.0, 6.0, 3.0, 3.0, 7.0, 10.0, 9.0, 8.0, 6.0, 1.0, 1.0, 6.0, 4.0, 5.0, 2.0, 3.0, 8.0, 3.0, 2.0, 5.0, 8.0, 1.0, 6.0, 10.0, 2.0, 1.0, 4.0, 8.0, 6.0, 2.0, 5.0, 9.0, 8.0, 2.0, 5.0, 2.0, 5.0, 4.0, 5.0, 3.0, 9.0, 5.0, 10.0, 9.0, 7.0, 4.0, 1.0, 9.0, 10.0, 4.0, 7.0, 10.0, 9.0, 6.0, 8.0, 3.0, 6.0, 3.0, 7.0, 9.0, 10.0, 3.0, 5.0, 4.0, 4.0, 8.0, 5.0, 4.0, 9.0, 6.0, 6.0, 6.0, 3.0, 4.0, 3.0, 6.0, 8.0, 5.0, 7.0, 6.0, 3.0, 3.0, 8.0, 5.0, 4.0, 10.0, 10.0, 6.0, 3.0, 2.0, 9.0, 8.0, 9.0, 9.0, 1.0, 4.0, 4.0, 9.0, 9.0, 3.0, 10.0, 1.0, 1.0, 3.0, 3.0, 1.0, 5.0, 5.0, 4.0, 5.0, 4.0, 6.0, 8.0, 2.0, 6.0, 3.0, 10.0, 1.0, 2.0, 9.0, 1.0, 3.0, 3.0, 7.0, 3.0, 3.0]
global b_x = 5
global d_y = [4.0, 2.0, 1.0, 6.0, 9.0, 9.0, 8.0, 1.0, 1.0, 9.0, 1.0, 5.0, 10.0, 2.0, 10.0, 6.0, 8.0, 6.0, 2.0, 1.0, 2.0, 5.0, 1.0, 5.0, 1.0, 2.0, 8.0, 9.0, 2.0, 5.0, 2.0, 9.0, 7.0, 4.0, 7.0, 9.0, 2.0, 2.0, 3.0, 4.0, 6.0, 1.0, 4.0, 4.0, 10.0, 8.0, 4.0, 10.0, 3.0, 9.0, 7.0, 10.0, 5.0, 8.0, 1.0, 2.0, 6.0, 5.0, 7.0, 2.0, 6.0, 7.0, 8.0, 4.0, 4.0, 4.0, 6.0, 4.0, 6.0, 2.0, 4.0, 1.0, 3.0, 10.0, 9.0, 6.0, 8.0, 10.0, 6.0, 8.0, 2.0, 1.0, 7.0, 5.0, 3.0, 8.0, 2.0, 10.0, 10.0, 8.0, 4.0, 10.0, 3.0, 7.0, 3.0, 2.0, 2.0, 6.0, 1.0, 7.0, 8.0, 10.0, 3.0, 10.0, 6.0, 3.0, 5.0, 5.0, 10.0, 8.0, 4.0, 2.0, 10.0, 6.0, 2.0, 7.0, 8.0, 3.0, 1.0, 5.0, 7.0, 10.0, 9.0, 7.0, 5.0, 10.0, 9.0, 3.0, 7.0, 4.0, 1.0, 10.0, 4.0, 4.0, 1.0, 5.0, 2.0, 8.0, 7.0, 7.0, 8.0, 9.0, 3.0, 10.0, 9.0, 7.0, 4.0, 4.0, 10.0, 6.0, 6.0, 10.0, 7.0]
global b_y = 10
global p = [0.776, 0.444, 0.004, 0.069, 0.69, 0.115, 0.627, 0.748, 0.134, 0.495, 0.807, 0.439, 0.535, 0.791, 0.272, 0.502, 0.009, 0.751, 0.594, 0.539, 0.451, 0.055, 0.299, 0.66, 0.147, 0.23, 0.765, 0.582, 0.969, 0.078, 0.75, 0.428, 0.912, 0.129, 0.321, 0.805, 0.93, 0.519, 0.489, 0.581, 0.667, 0.848, 0.508, 0.037, 0.932, 0.315, 0.596, 0.504, 0.715, 0.52, 0.399, 0.574, 0.544, 0.98, 0.541, 0.471, 0.712, 0.387, 0.3, 0.891, 0.102, 0.578, 0.011, 0.1, 0.491, 0.615, 0.866, 0.835, 0.612, 0.603, 0.059, 0.176, 0.062, 0.087, 0.719, 0.741, 0.127, 0.58, 0.025, 0.92, 0.027, 0.097, 0.331, 0.892, 0.729, 0.45, 0.484, 0.769, 0.247, 0.988, 0.383, 0.389, 0.013, 0.645, 0.122, 0.997, 0.774, 0.75, 0.791, 0.638, 0.215, 0.909, 0.771, 0.319, 0.583, 0.675, 0.612, 0.123, 0.417, 0.287, 0.666, 0.89, 0.952, 0.75, 0.291, 0.735, 0.448, 0.158, 0.228, 0.705, 0.114, 0.223, 0.242, 0.217, 0.488, 0.055, 0.734, 0.026, 0.705, 0.305, 0.56, 0.643, 0.128, 0.483, 0.485, 0.395, 0.979, 0.727, 0.494, 0.51, 0.893, 0.302, 0.461, 0.917, 0.167, 0.862, 0.327, 0.293, 0.65, 0.306, 0.039, 0.356, 0.521]
global q = [0.814, 0.648, 0.708, 0.941, 0.848, 0.505, 0.881, 0.758, 0.291, 0.518, 0.88, 0.875, 0.784, 0.973, 0.73, 0.968, 0.441, 0.802, 0.877, 0.965, 0.984, 0.444, 0.707, 0.746, 0.992, 0.425, 0.988, 0.896, 0.979, 0.337, 0.902, 0.98, 0.916, 0.52, 0.46, 0.836, 0.993, 0.946, 0.637, 0.86, 0.747, 0.976, 0.578, 0.497, 0.998, 0.681, 0.654, 0.689, 0.99, 0.777, 0.497, 0.865, 0.811, 0.991, 0.589, 0.893, 0.967, 0.588, 0.848, 0.943, 0.939, 0.578, 0.615, 0.538, 0.948, 0.786, 0.975, 0.986, 0.982, 0.615, 0.7, 0.767, 0.784, 0.481, 0.992, 0.744, 0.407, 0.654, 0.961, 0.969, 0.61, 0.787, 0.335, 0.926, 0.989, 0.565, 0.859, 0.9, 0.663, 0.995, 0.87, 0.785, 0.245, 0.662, 0.761, 0.999, 0.854, 0.925, 0.964, 0.769, 0.428, 0.938, 0.775, 0.4, 0.943, 0.948, 0.828, 0.283, 0.626, 0.692, 0.729, 0.91, 0.97, 0.876, 0.984, 0.973, 0.641, 0.933, 0.658, 0.858, 0.347, 0.41, 0.306, 0.743, 0.704, 0.289, 0.981, 0.615, 0.99, 0.959, 0.904, 0.94, 0.263, 0.78, 0.573, 0.587, 0.996, 0.898, 0.662, 0.768, 0.999, 0.825, 0.879, 0.999, 0.518, 0.899, 0.579, 0.527, 0.954, 0.875, 0.339, 0.455, 0.661]
global origin = 1
global destination = 40