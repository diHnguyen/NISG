global arcs = [1 6; 1 19; 1 25; 1 31; 1 41; 2 4; 2 5; 2 18; 2 20; 2 24; 2 25; 2 27; 2 28; 2 43; 3 6; 3 10; 3 12; 3 13; 3 18; 3 27; 3 31; 4 32; 4 47; 5 7; 5 18; 5 27; 5 34; 6 4; 6 12; 7 12; 7 15; 7 23; 7 46; 7 48; 8 11; 8 26; 8 36; 9 2; 9 39; 9 40; 9 41; 9 47; 10 3; 10 4; 10 6; 10 8; 10 14; 10 15; 10 26; 10 40; 10 45; 11 15; 11 18; 11 33; 11 38; 11 50; 12 2; 12 19; 12 20; 12 39; 13 10; 13 25; 13 28; 13 40; 14 2; 15 12; 15 31; 15 37; 15 41; 15 45; 16 4; 16 28; 16 38; 17 22; 17 24; 17 29; 17 44; 18 5; 18 25; 18 34; 18 35; 19 2; 19 10; 19 28; 19 32; 19 36; 19 37; 19 40; 20 2; 20 3; 20 15; 20 33; 20 34; 20 39; 20 40; 20 43; 20 44; 21 5; 21 12; 21 25; 21 28; 21 45; 22 6; 22 28; 22 36; 22 50; 23 11; 23 21; 23 24; 23 30; 23 38; 23 49; 24 12; 24 21; 24 23; 24 38; 25 19; 25 33; 25 46; 26 10; 26 21; 26 37; 26 39; 27 24; 27 26; 27 33; 27 39; 27 50; 28 12; 28 16; 28 42; 29 12; 29 13; 29 15; 30 25; 30 36; 31 25; 31 26; 31 30; 31 33; 31 47; 32 2; 32 15; 32 26; 32 29; 32 31; 33 3; 33 13; 33 29; 33 34; 33 48; 33 49; 34 2; 34 9; 34 17; 34 18; 34 20; 34 24; 34 27; 34 36; 34 48; 35 6; 35 13; 35 18; 35 28; 35 36; 35 37; 36 5; 36 14; 36 48; 37 3; 37 19; 37 30; 37 38; 37 41; 37 50; 38 9; 39 25; 39 38; 39 48; 40 8; 41 11; 41 12; 41 13; 41 16; 41 43; 41 46; 41 48; 41 50; 42 12; 42 15; 42 28; 42 39; 42 40; 43 2; 43 13; 43 14; 43 18; 43 22; 43 42; 43 48; 44 2; 44 17; 44 32; 44 35; 44 39; 45 2; 45 19; 45 37; 46 7; 46 22; 47 2; 47 3; 47 10; 47 13; 47 19; 47 31; 47 37; 47 39; 47 45; 47 50; 48 10; 48 27; 48 29; 48 32; 48 44; 48 49; 49 4; 49 9; 49 11; 49 46]
global d_x = [10.0, 1.0, 3.0, 4.0, 6.0, 10.0, 4.0, 4.0, 1.0, 2.0, 1.0, 3.0, 9.0, 3.0, 2.0, 1.0, 7.0, 4.0, 7.0, 1.0, 8.0, 7.0, 2.0, 6.0, 10.0, 6.0, 9.0, 3.0, 3.0, 4.0, 4.0, 3.0, 3.0, 2.0, 6.0, 4.0, 3.0, 3.0, 1.0, 5.0, 5.0, 4.0, 3.0, 10.0, 2.0, 1.0, 3.0, 10.0, 5.0, 4.0, 9.0, 4.0, 5.0, 10.0, 9.0, 2.0, 5.0, 7.0, 10.0, 1.0, 8.0, 9.0, 10.0, 5.0, 4.0, 8.0, 8.0, 2.0, 8.0, 9.0, 6.0, 5.0, 10.0, 4.0, 4.0, 4.0, 1.0, 9.0, 2.0, 10.0, 5.0, 2.0, 9.0, 7.0, 9.0, 2.0, 6.0, 7.0, 4.0, 7.0, 6.0, 7.0, 8.0, 6.0, 9.0, 2.0, 1.0, 5.0, 2.0, 5.0, 4.0, 1.0, 10.0, 10.0, 7.0, 2.0, 8.0, 10.0, 5.0, 1.0, 2.0, 8.0, 4.0, 4.0, 8.0, 10.0, 5.0, 8.0, 6.0, 1.0, 4.0, 10.0, 1.0, 9.0, 3.0, 7.0, 2.0, 3.0, 7.0, 2.0, 1.0, 6.0, 2.0, 5.0, 1.0, 3.0, 7.0, 4.0, 10.0, 8.0, 5.0, 10.0, 6.0, 8.0, 1.0, 6.0, 8.0, 4.0, 2.0, 3.0, 9.0, 6.0, 4.0, 7.0, 5.0, 8.0, 3.0, 9.0, 10.0, 3.0, 7.0, 2.0, 4.0, 9.0, 8.0, 2.0, 2.0, 6.0, 10.0, 6.0, 5.0, 7.0, 10.0, 4.0, 1.0, 10.0, 1.0, 8.0, 3.0, 5.0, 7.0, 2.0, 7.0, 9.0, 6.0, 4.0, 2.0, 5.0, 7.0, 10.0, 4.0, 6.0, 9.0, 8.0, 5.0, 4.0, 5.0, 2.0, 3.0, 9.0, 9.0, 5.0, 2.0, 1.0, 5.0, 10.0, 5.0, 9.0, 5.0, 9.0, 3.0, 7.0, 2.0, 10.0, 3.0, 5.0, 8.0, 6.0, 5.0, 10.0, 1.0, 5.0, 8.0, 5.0, 5.0, 2.0, 9.0, 4.0, 9.0, 10.0, 10.0]
global b_x = 5
global d_y = [3.0, 9.0, 3.0, 5.0, 6.0, 5.0, 7.0, 5.0, 9.0, 6.0, 8.0, 4.0, 7.0, 7.0, 9.0, 7.0, 2.0, 9.0, 1.0, 5.0, 7.0, 7.0, 10.0, 1.0, 1.0, 5.0, 2.0, 8.0, 8.0, 1.0, 1.0, 2.0, 9.0, 2.0, 1.0, 6.0, 9.0, 2.0, 8.0, 10.0, 7.0, 3.0, 1.0, 6.0, 8.0, 8.0, 7.0, 10.0, 10.0, 1.0, 1.0, 10.0, 8.0, 2.0, 2.0, 9.0, 1.0, 3.0, 2.0, 2.0, 2.0, 8.0, 1.0, 3.0, 10.0, 3.0, 10.0, 5.0, 3.0, 4.0, 8.0, 2.0, 4.0, 6.0, 3.0, 8.0, 6.0, 6.0, 1.0, 8.0, 6.0, 1.0, 4.0, 5.0, 8.0, 7.0, 6.0, 2.0, 4.0, 3.0, 3.0, 3.0, 9.0, 3.0, 3.0, 7.0, 3.0, 8.0, 7.0, 6.0, 5.0, 9.0, 9.0, 9.0, 2.0, 5.0, 8.0, 3.0, 2.0, 6.0, 5.0, 10.0, 2.0, 8.0, 1.0, 10.0, 3.0, 8.0, 3.0, 10.0, 3.0, 8.0, 6.0, 10.0, 1.0, 4.0, 7.0, 7.0, 3.0, 2.0, 7.0, 6.0, 1.0, 5.0, 7.0, 7.0, 6.0, 3.0, 6.0, 6.0, 10.0, 8.0, 3.0, 6.0, 5.0, 8.0, 9.0, 2.0, 3.0, 8.0, 2.0, 2.0, 2.0, 9.0, 8.0, 9.0, 7.0, 10.0, 1.0, 10.0, 6.0, 4.0, 9.0, 7.0, 1.0, 1.0, 3.0, 9.0, 5.0, 6.0, 10.0, 7.0, 7.0, 4.0, 3.0, 2.0, 1.0, 4.0, 1.0, 4.0, 3.0, 4.0, 7.0, 5.0, 1.0, 8.0, 9.0, 4.0, 8.0, 8.0, 6.0, 2.0, 10.0, 1.0, 2.0, 7.0, 10.0, 1.0, 6.0, 10.0, 6.0, 10.0, 9.0, 1.0, 2.0, 7.0, 4.0, 9.0, 8.0, 6.0, 7.0, 8.0, 5.0, 7.0, 1.0, 7.0, 9.0, 6.0, 1.0, 2.0, 9.0, 7.0, 9.0, 6.0, 5.0, 6.0, 2.0, 4.0, 4.0, 9.0, 2.0]
global b_y = 10
global p = [0.06, 0.809, 0.625, 0.339, 0.51, 0.367, 0.007, 0.297, 0.534, 0.967, 0.306, 0.434, 0.301, 0.17, 0.365, 0.964, 0.874, 0.761, 0.389, 0.829, 0.954, 0.104, 0.813, 0.742, 0.298, 0.656, 0.448, 0.423, 0.87, 0.237, 0.266, 0.834, 0.663, 0.874, 0.353, 0.45, 0.852, 0.306, 0.075, 0.37, 0.104, 0.022, 0.453, 0.295, 0.659, 0.191, 0.885, 0.683, 0.558, 0.447, 0.64, 0.889, 0.593, 0.245, 0.753, 0.894, 0.387, 0.375, 0.19, 0.37, 0.314, 0.318, 0.191, 0.638, 0.537, 0.593, 0.861, 0.248, 0.13, 0.516, 0.754, 0.942, 0.273, 0.373, 0.916, 0.664, 0.336, 0.117, 0.907, 0.206, 0.831, 0.38, 0.298, 0.376, 0.295, 0.194, 0.382, 0.929, 0.326, 0.675, 0.489, 0.397, 0.891, 0.321, 0.931, 0.284, 0.917, 0.137, 0.602, 0.103, 0.023, 0.613, 0.789, 0.132, 0.941, 0.722, 0.086, 0.097, 0.99, 0.232, 0.675, 0.833, 0.692, 0.498, 0.261, 0.491, 0.207, 0.14, 0.297, 0.001, 0.28, 0.685, 0.278, 0.767, 0.576, 0.175, 0.49, 0.716, 0.618, 0.51, 0.22, 0.821, 0.977, 0.757, 0.547, 0.624, 0.352, 0.784, 0.881, 0.755, 0.46, 0.303, 0.199, 0.034, 0.116, 0.846, 0.865, 0.81, 0.688, 0.308, 0.243, 0.225, 0.703, 0.96, 0.432, 0.8, 0.312, 0.041, 0.538, 0.457, 0.599, 0.613, 0.437, 0.995, 0.986, 0.432, 0.627, 0.044, 0.257, 0.269, 0.146, 0.051, 0.651, 0.133, 0.79, 0.42, 0.102, 0.669, 0.146, 0.452, 0.982, 0.362, 0.074, 0.472, 0.78, 0.819, 0.213, 0.382, 0.577, 0.435, 0.939, 0.776, 0.133, 0.786, 0.414, 0.679, 0.764, 0.831, 0.651, 0.373, 0.208, 0.306, 0.95, 0.801, 0.374, 0.884, 0.397, 0.357, 0.455, 0.045, 0.954, 0.924, 0.126, 0.173, 0.924, 0.959, 0.545, 0.793, 0.072, 0.851, 0.16, 0.29, 0.523, 0.067, 0.493, 0.664, 0.456, 0.525, 0.19, 0.018, 0.932]
global q = [0.801, 0.964, 0.949, 0.984, 0.554, 0.919, 0.315, 0.82, 0.807, 0.974, 0.553, 0.499, 0.348, 0.188, 0.932, 0.974, 0.999, 0.971, 0.506, 0.899, 0.957, 0.269, 0.945, 0.851, 0.651, 0.835, 0.994, 0.626, 0.927, 0.748, 0.878, 0.937, 0.982, 0.895, 0.389, 0.902, 0.915, 0.761, 0.431, 0.871, 0.604, 0.585, 0.581, 0.489, 0.906, 0.483, 0.97, 0.732, 0.837, 0.736, 0.869, 0.986, 0.774, 0.443, 0.801, 0.983, 0.944, 0.425, 0.394, 0.476, 0.453, 0.973, 0.47, 0.676, 0.831, 0.828, 0.927, 0.712, 0.404, 0.983, 0.837, 0.953, 0.658, 0.8, 0.918, 0.689, 0.835, 0.445, 0.91, 0.756, 0.924, 0.752, 0.552, 0.703, 0.441, 0.622, 0.688, 0.961, 0.597, 0.84, 0.51, 0.649, 0.916, 0.568, 0.961, 0.874, 0.925, 0.754, 0.623, 0.185, 0.823, 0.81, 0.9, 0.237, 0.983, 0.977, 0.455, 0.391, 0.999, 0.659, 0.696, 0.967, 0.914, 0.63, 0.615, 0.519, 0.379, 0.736, 0.755, 0.449, 0.282, 0.888, 0.569, 0.847, 0.626, 0.995, 0.85, 0.778, 0.681, 0.564, 0.478, 0.922, 0.993, 0.796, 0.623, 0.678, 0.563, 0.979, 0.961, 0.811, 0.781, 0.37, 0.617, 0.3, 0.753, 0.982, 0.908, 0.945, 0.869, 0.897, 0.938, 0.402, 0.873, 0.993, 0.771, 0.83, 0.388, 0.72, 0.544, 0.698, 0.739, 0.669, 0.552, 0.995, 0.991, 0.613, 0.631, 0.103, 0.648, 0.328, 0.396, 0.426, 0.989, 0.536, 0.983, 0.835, 0.166, 0.863, 0.833, 0.945, 0.983, 0.998, 0.927, 0.578, 0.986, 0.907, 0.488, 0.977, 0.678, 0.806, 0.984, 0.797, 0.661, 0.869, 0.863, 0.988, 0.909, 0.967, 0.888, 0.966, 0.954, 0.36, 0.984, 0.922, 0.786, 0.905, 0.568, 0.375, 0.62, 0.879, 0.978, 0.967, 0.432, 0.48, 0.974, 0.967, 0.73, 0.833, 0.181, 0.956, 0.615, 0.797, 0.889, 0.63, 0.839, 0.736, 0.827, 0.675, 0.795, 0.456, 0.994]
global origin = 1
global destination = 50