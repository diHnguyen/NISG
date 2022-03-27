global arcs = [1 23; 1 28; 1 33; 1 38; 2 11; 2 13; 2 30; 3 32; 4 2; 4 19; 4 20; 4 24; 4 28; 4 36; 5 2; 5 3; 5 14; 5 23; 5 27; 5 28; 5 30; 5 34; 5 36; 5 38; 5 40; 6 9; 6 16; 6 23; 6 24; 6 31; 7 9; 7 14; 7 27; 7 28; 8 12; 8 17; 8 20; 8 27; 9 3; 10 5; 10 18; 10 25; 10 27; 10 35; 11 4; 11 7; 11 19; 11 21; 12 24; 13 5; 13 7; 14 10; 14 15; 14 16; 14 19; 14 27; 14 28; 14 30; 14 31; 15 9; 15 11; 15 26; 15 28; 16 9; 16 10; 16 14; 16 17; 16 24; 16 34; 17 7; 17 8; 17 10; 17 22; 17 25; 18 4; 18 8; 18 24; 18 31; 18 33; 19 9; 19 14; 19 26; 19 31; 19 39; 20 8; 20 11; 20 21; 20 23; 20 31; 20 39; 21 17; 21 18; 21 28; 21 32; 22 9; 22 12; 22 16; 22 38; 23 8; 23 19; 23 28; 23 29; 23 34; 23 35; 23 36; 24 16; 24 29; 24 30; 24 31; 24 37; 25 14; 25 21; 25 24; 25 28; 26 13; 26 22; 26 33; 26 37; 26 38; 27 8; 27 23; 27 25; 27 26; 27 34; 28 8; 28 12; 28 13; 28 33; 28 37; 28 39; 29 7; 29 27; 29 28; 29 39; 30 4; 30 11; 30 14; 30 23; 31 2; 31 4; 31 6; 31 18; 31 19; 32 22; 32 28; 32 34; 32 40; 33 6; 34 13; 34 14; 34 31; 34 35; 35 3; 35 11; 35 17; 35 22; 36 13; 36 26; 36 28; 37 2; 37 20; 37 31; 38 10; 38 25; 39 3; 39 8; 39 17; 39 37; 39 38]
global d_x = [1.0, 3.0, 10.0, 2.0, 5.0, 1.0, 4.0, 5.0, 8.0, 10.0, 10.0, 6.0, 8.0, 2.0, 5.0, 1.0, 1.0, 8.0, 3.0, 3.0, 3.0, 5.0, 10.0, 8.0, 10.0, 2.0, 8.0, 2.0, 3.0, 8.0, 6.0, 5.0, 8.0, 2.0, 7.0, 2.0, 9.0, 1.0, 3.0, 7.0, 4.0, 1.0, 10.0, 2.0, 10.0, 7.0, 9.0, 10.0, 8.0, 1.0, 4.0, 5.0, 1.0, 9.0, 9.0, 7.0, 6.0, 6.0, 7.0, 2.0, 2.0, 4.0, 5.0, 5.0, 3.0, 9.0, 6.0, 8.0, 6.0, 3.0, 9.0, 1.0, 10.0, 10.0, 5.0, 4.0, 9.0, 5.0, 3.0, 1.0, 8.0, 9.0, 1.0, 9.0, 10.0, 4.0, 1.0, 10.0, 10.0, 2.0, 5.0, 1.0, 8.0, 7.0, 7.0, 4.0, 9.0, 6.0, 4.0, 5.0, 5.0, 7.0, 8.0, 7.0, 6.0, 4.0, 7.0, 6.0, 4.0, 2.0, 6.0, 8.0, 5.0, 3.0, 4.0, 9.0, 2.0, 7.0, 9.0, 4.0, 9.0, 6.0, 5.0, 8.0, 4.0, 5.0, 4.0, 4.0, 1.0, 6.0, 7.0, 3.0, 1.0, 10.0, 4.0, 8.0, 5.0, 4.0, 2.0, 7.0, 8.0, 7.0, 2.0, 7.0, 8.0, 2.0, 9.0, 7.0, 6.0, 5.0, 2.0, 4.0, 4.0, 1.0, 5.0, 5.0, 10.0, 6.0, 8.0, 4.0, 3.0, 4.0, 5.0, 7.0, 1.0, 3.0, 9.0, 9.0, 4.0]
global b_x = 5
global d_y = [8.0, 6.0, 10.0, 4.0, 6.0, 8.0, 5.0, 1.0, 7.0, 3.0, 4.0, 10.0, 4.0, 10.0, 10.0, 2.0, 7.0, 4.0, 6.0, 9.0, 4.0, 2.0, 7.0, 7.0, 4.0, 1.0, 9.0, 10.0, 1.0, 6.0, 1.0, 7.0, 1.0, 4.0, 2.0, 7.0, 4.0, 10.0, 5.0, 6.0, 6.0, 7.0, 8.0, 7.0, 8.0, 9.0, 9.0, 6.0, 4.0, 8.0, 9.0, 6.0, 8.0, 8.0, 9.0, 2.0, 3.0, 5.0, 5.0, 2.0, 10.0, 3.0, 10.0, 1.0, 8.0, 6.0, 10.0, 6.0, 5.0, 4.0, 4.0, 7.0, 1.0, 9.0, 7.0, 1.0, 1.0, 10.0, 10.0, 5.0, 4.0, 1.0, 7.0, 4.0, 10.0, 7.0, 6.0, 5.0, 2.0, 3.0, 6.0, 8.0, 2.0, 10.0, 10.0, 8.0, 7.0, 8.0, 7.0, 10.0, 5.0, 3.0, 2.0, 5.0, 2.0, 7.0, 10.0, 4.0, 3.0, 7.0, 8.0, 6.0, 10.0, 5.0, 4.0, 1.0, 3.0, 4.0, 2.0, 8.0, 10.0, 10.0, 4.0, 5.0, 1.0, 10.0, 8.0, 6.0, 7.0, 1.0, 1.0, 1.0, 2.0, 3.0, 3.0, 8.0, 2.0, 3.0, 6.0, 5.0, 1.0, 10.0, 3.0, 7.0, 1.0, 1.0, 7.0, 8.0, 8.0, 5.0, 8.0, 10.0, 1.0, 9.0, 1.0, 5.0, 1.0, 8.0, 3.0, 7.0, 5.0, 6.0, 5.0, 6.0, 6.0, 2.0, 9.0, 8.0, 3.0]
global b_y = 10
global p = [0.732, 0.62, 0.335, 0.597, 0.809, 0.141, 0.956, 0.488, 0.91, 0.805, 0.331, 0.99, 0.383, 0.494, 0.8, 0.296, 0.023, 0.961, 0.054, 0.403, 0.991, 0.434, 0.124, 0.569, 0.987, 0.092, 0.042, 0.749, 0.048, 0.713, 0.466, 0.405, 0.506, 0.748, 0.984, 0.181, 0.477, 0.007, 0.942, 0.109, 0.758, 0.58, 0.515, 0.902, 0.901, 0.48, 0.026, 0.375, 0.648, 0.901, 0.609, 0.387, 0.01, 0.324, 0.257, 0.585, 0.786, 0.14, 0.528, 0.979, 0.901, 0.27, 0.903, 0.499, 0.238, 0.316, 0.212, 0.741, 0.863, 0.36, 0.825, 0.372, 0.577, 0.377, 0.056, 0.376, 0.83, 0.153, 0.789, 0.44, 0.493, 0.23, 0.314, 0.48, 0.004, 0.568, 0.888, 0.932, 0.476, 0.199, 0.819, 0.087, 0.337, 0.338, 0.098, 0.062, 0.372, 0.751, 0.53, 0.383, 0.735, 0.813, 0.228, 0.036, 0.618, 0.448, 0.721, 0.806, 0.805, 0.84, 0.423, 0.754, 0.832, 0.301, 0.567, 0.878, 0.597, 0.67, 0.406, 0.367, 0.212, 0.614, 0.301, 0.828, 0.254, 0.233, 0.001, 0.938, 0.539, 0.177, 0.814, 0.136, 0.567, 0.574, 0.948, 0.377, 0.076, 0.044, 0.345, 0.197, 0.059, 0.796, 0.627, 0.644, 0.175, 0.075, 0.034, 0.743, 0.281, 0.144, 0.798, 0.516, 0.313, 0.392, 0.397, 0.355, 0.936, 0.262, 0.912, 0.896, 0.052, 0.716, 0.504, 0.77, 0.431, 0.712, 0.959, 0.864, 0.653]
global q = [0.818, 0.868, 0.658, 0.958, 0.904, 0.348, 0.996, 0.752, 0.934, 0.942, 0.77, 0.998, 0.434, 0.628, 0.817, 0.455, 0.871, 0.978, 0.368, 0.692, 0.999, 0.935, 0.131, 0.571, 0.989, 0.598, 0.709, 0.85, 0.896, 0.888, 0.852, 0.481, 0.775, 0.77, 0.998, 0.214, 0.684, 0.641, 0.979, 0.813, 0.986, 0.744, 0.999, 0.972, 0.922, 0.754, 0.184, 0.489, 0.747, 0.919, 0.616, 0.54, 0.576, 0.834, 0.451, 0.619, 0.893, 0.997, 0.657, 0.981, 0.926, 0.82, 0.987, 0.855, 0.743, 0.641, 0.279, 0.999, 0.963, 0.639, 0.94, 0.461, 0.752, 0.69, 0.401, 0.576, 0.998, 0.639, 0.927, 0.539, 0.642, 0.882, 0.73, 0.756, 0.112, 0.788, 0.977, 0.943, 0.974, 0.962, 0.968, 0.996, 0.561, 0.451, 0.489, 0.103, 0.759, 0.867, 0.739, 0.738, 0.954, 0.846, 0.592, 0.882, 0.629, 0.857, 0.761, 0.944, 0.974, 0.935, 0.585, 0.755, 0.943, 0.301, 0.747, 0.884, 0.851, 0.712, 0.935, 0.943, 0.486, 0.856, 0.811, 0.951, 0.665, 0.241, 0.392, 0.957, 0.738, 0.416, 0.941, 0.716, 0.72, 0.935, 0.973, 0.671, 0.466, 0.946, 0.511, 0.838, 0.226, 0.882, 0.671, 0.841, 0.942, 0.956, 0.094, 0.87, 0.775, 0.52, 0.873, 0.769, 0.791, 0.533, 0.59, 0.627, 0.948, 0.51, 0.956, 0.934, 0.428, 0.988, 0.984, 0.932, 0.812, 0.879, 0.978, 0.873, 0.96]
global origin = 1
global destination = 40