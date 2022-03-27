global arcs = [1 6; 1 12; 1 18; 1 19; 1 30; 1 37; 2 14; 2 15; 2 38; 2 46; 2 50; 3 13; 3 17; 3 22; 3 35; 3 37; 3 49; 3 50; 4 3; 4 22; 4 35; 5 8; 5 10; 5 18; 5 43; 5 46; 5 49; 6 3; 6 7; 6 18; 7 9; 8 29; 9 12; 9 41; 9 42; 10 15; 10 33; 10 34; 10 42; 10 46; 10 49; 10 50; 11 15; 11 40; 12 4; 12 27; 12 31; 12 34; 12 42; 12 49; 13 2; 13 17; 13 21; 13 29; 13 39; 13 43; 14 12; 14 19; 14 20; 15 13; 15 17; 16 39; 16 43; 17 5; 17 6; 17 13; 17 16; 17 25; 17 31; 17 50; 18 5; 18 7; 18 12; 18 23; 18 31; 19 11; 19 16; 19 29; 19 34; 20 7; 20 10; 20 14; 20 33; 21 2; 21 5; 21 28; 22 7; 22 11; 22 32; 23 24; 23 30; 23 39; 23 42; 24 6; 24 10; 24 42; 24 45; 24 48; 25 6; 25 9; 25 14; 25 22; 25 24; 25 34; 25 36; 25 40; 25 44; 25 45; 26 28; 26 30; 26 50; 27 3; 27 29; 27 38; 27 49; 27 50; 28 34; 29 8; 29 12; 29 18; 29 44; 30 29; 30 42; 31 4; 31 9; 31 39; 31 46; 32 30; 32 31; 32 39; 33 8; 33 17; 33 29; 33 45; 34 25; 34 27; 34 35; 34 40; 35 6; 35 7; 35 11; 36 9; 36 27; 36 44; 37 9; 37 12; 37 14; 37 30; 37 44; 38 19; 38 22; 38 25; 38 28; 38 34; 38 40; 38 49; 39 3; 39 6; 39 20; 39 21; 39 22; 39 34; 39 49; 40 39; 40 48; 41 10; 41 11; 41 22; 41 29; 41 48; 42 2; 42 15; 42 20; 42 46; 43 25; 43 35; 43 50; 44 11; 44 16; 44 18; 44 20; 44 27; 44 28; 44 31; 44 33; 44 45; 44 46; 44 47; 44 48; 45 3; 45 19; 45 37; 45 50; 46 12; 46 19; 46 29; 46 30; 46 42; 47 5; 47 16; 47 35; 48 12; 48 13; 48 15; 48 22; 48 26; 48 27; 48 34; 48 44; 48 45; 48 50; 49 11; 49 17; 49 19; 49 21; 49 24]
global d_x = [3.0, 8.0, 3.0, 10.0, 7.0, 9.0, 3.0, 3.0, 6.0, 6.0, 7.0, 9.0, 7.0, 6.0, 2.0, 7.0, 6.0, 5.0, 7.0, 1.0, 2.0, 10.0, 6.0, 6.0, 4.0, 6.0, 6.0, 10.0, 6.0, 1.0, 4.0, 9.0, 5.0, 5.0, 1.0, 4.0, 7.0, 9.0, 4.0, 5.0, 1.0, 10.0, 1.0, 9.0, 9.0, 9.0, 9.0, 6.0, 6.0, 5.0, 8.0, 10.0, 1.0, 8.0, 6.0, 4.0, 1.0, 7.0, 2.0, 2.0, 5.0, 5.0, 6.0, 9.0, 9.0, 6.0, 1.0, 3.0, 4.0, 5.0, 5.0, 3.0, 8.0, 7.0, 5.0, 9.0, 4.0, 8.0, 9.0, 3.0, 1.0, 2.0, 7.0, 6.0, 5.0, 5.0, 5.0, 5.0, 4.0, 2.0, 5.0, 7.0, 7.0, 8.0, 4.0, 2.0, 8.0, 1.0, 6.0, 5.0, 4.0, 5.0, 1.0, 6.0, 7.0, 5.0, 10.0, 5.0, 3.0, 2.0, 2.0, 7.0, 8.0, 3.0, 8.0, 1.0, 8.0, 3.0, 9.0, 1.0, 8.0, 8.0, 3.0, 8.0, 5.0, 8.0, 6.0, 10.0, 5.0, 4.0, 8.0, 1.0, 5.0, 4.0, 7.0, 6.0, 3.0, 5.0, 8.0, 2.0, 7.0, 4.0, 1.0, 7.0, 1.0, 8.0, 1.0, 4.0, 4.0, 7.0, 2.0, 6.0, 10.0, 4.0, 1.0, 4.0, 9.0, 1.0, 2.0, 10.0, 9.0, 7.0, 2.0, 5.0, 1.0, 6.0, 4.0, 7.0, 3.0, 9.0, 5.0, 4.0, 6.0, 2.0, 1.0, 10.0, 9.0, 8.0, 1.0, 10.0, 10.0, 9.0, 5.0, 9.0, 3.0, 3.0, 8.0, 4.0, 3.0, 8.0, 1.0, 5.0, 2.0, 8.0, 9.0, 10.0, 4.0, 1.0, 7.0, 8.0, 5.0, 4.0, 1.0, 4.0, 10.0, 6.0, 1.0, 5.0, 9.0, 5.0, 9.0, 4.0, 8.0, 10.0, 6.0, 5.0]
global b_x = 5
global d_y = [4.0, 8.0, 6.0, 6.0, 6.0, 9.0, 9.0, 4.0, 8.0, 10.0, 8.0, 9.0, 4.0, 2.0, 7.0, 2.0, 6.0, 8.0, 8.0, 1.0, 7.0, 1.0, 2.0, 2.0, 1.0, 10.0, 7.0, 10.0, 7.0, 8.0, 2.0, 8.0, 9.0, 9.0, 6.0, 10.0, 5.0, 2.0, 9.0, 3.0, 1.0, 7.0, 7.0, 2.0, 6.0, 5.0, 10.0, 2.0, 9.0, 9.0, 7.0, 2.0, 9.0, 7.0, 1.0, 4.0, 2.0, 8.0, 3.0, 8.0, 5.0, 8.0, 3.0, 5.0, 3.0, 5.0, 8.0, 3.0, 4.0, 2.0, 3.0, 5.0, 3.0, 3.0, 5.0, 7.0, 3.0, 5.0, 2.0, 8.0, 5.0, 10.0, 3.0, 6.0, 8.0, 1.0, 5.0, 2.0, 7.0, 1.0, 8.0, 2.0, 6.0, 8.0, 2.0, 8.0, 9.0, 2.0, 7.0, 1.0, 7.0, 2.0, 2.0, 2.0, 4.0, 1.0, 10.0, 2.0, 2.0, 4.0, 8.0, 9.0, 8.0, 5.0, 10.0, 9.0, 4.0, 9.0, 4.0, 8.0, 7.0, 2.0, 3.0, 8.0, 5.0, 2.0, 6.0, 10.0, 7.0, 2.0, 3.0, 8.0, 3.0, 9.0, 5.0, 2.0, 4.0, 5.0, 5.0, 4.0, 5.0, 2.0, 2.0, 2.0, 1.0, 5.0, 8.0, 3.0, 5.0, 2.0, 6.0, 4.0, 6.0, 1.0, 2.0, 1.0, 3.0, 4.0, 6.0, 5.0, 8.0, 1.0, 10.0, 9.0, 5.0, 8.0, 2.0, 6.0, 5.0, 2.0, 7.0, 3.0, 4.0, 2.0, 7.0, 2.0, 3.0, 2.0, 3.0, 7.0, 3.0, 7.0, 4.0, 10.0, 3.0, 2.0, 8.0, 8.0, 2.0, 7.0, 5.0, 1.0, 3.0, 2.0, 8.0, 4.0, 8.0, 2.0, 3.0, 2.0, 9.0, 4.0, 1.0, 10.0, 7.0, 9.0, 4.0, 3.0, 6.0, 8.0, 8.0, 10.0, 8.0, 4.0, 10.0, 2.0]
global b_y = 10
global p = [0.126, 0.698, 0.686, 0.98, 0.685, 0.211, 0.695, 0.302, 0.202, 0.086, 0.413, 0.242, 0.023, 0.896, 0.421, 0.113, 0.491, 0.026, 0.506, 0.552, 0.171, 0.878, 0.966, 0.878, 0.867, 0.795, 0.531, 0.372, 0.695, 0.957, 0.756, 0.669, 0.896, 0.47, 0.931, 0.351, 0.28, 0.839, 0.356, 0.982, 0.108, 0.644, 0.712, 0.076, 0.411, 0.157, 0.745, 0.447, 0.149, 0.823, 0.682, 0.417, 0.416, 0.261, 0.723, 0.855, 0.452, 0.025, 0.437, 0.429, 0.544, 0.684, 0.282, 0.424, 0.202, 0.492, 0.355, 0.611, 0.163, 0.986, 0.051, 0.788, 0.377, 0.29, 0.218, 0.878, 0.769, 0.235, 0.956, 0.365, 0.093, 0.833, 0.05, 0.734, 0.659, 0.855, 0.44, 0.314, 0.065, 0.394, 0.441, 0.547, 0.02, 0.605, 0.612, 0.771, 0.491, 0.62, 0.471, 0.519, 0.579, 0.27, 0.17, 0.904, 0.1, 0.19, 0.958, 0.414, 0.589, 0.659, 0.237, 0.888, 0.795, 0.72, 0.732, 0.45, 0.144, 0.528, 0.237, 0.527, 0.928, 0.693, 0.476, 0.012, 0.972, 0.405, 0.476, 0.72, 0.299, 0.611, 0.532, 0.206, 0.856, 0.598, 0.719, 0.189, 0.397, 0.294, 0.074, 0.001, 0.344, 0.232, 0.713, 0.302, 0.227, 0.702, 0.604, 0.676, 0.417, 0.586, 0.642, 0.673, 0.442, 0.083, 0.1, 0.672, 0.07, 0.251, 0.304, 0.928, 0.708, 0.133, 0.398, 0.128, 0.448, 0.359, 0.873, 0.811, 0.524, 0.248, 0.108, 0.747, 0.293, 0.698, 0.568, 0.619, 0.985, 0.744, 0.772, 0.099, 0.305, 0.369, 0.561, 0.333, 0.544, 0.895, 0.573, 0.172, 0.708, 0.438, 0.245, 0.582, 0.84, 0.105, 0.811, 0.803, 0.674, 0.676, 0.448, 0.179, 0.658, 0.199, 0.958, 0.512, 0.641, 0.41, 0.772, 0.943, 0.502, 0.386, 0.204, 0.5, 0.235, 0.45, 0.372, 0.163]
global q = [0.37, 0.775, 0.983, 0.99, 0.956, 0.769, 0.93, 0.522, 0.696, 0.181, 0.962, 0.518, 0.709, 0.911, 0.956, 0.157, 0.719, 0.318, 0.633, 0.552, 0.782, 0.882, 0.976, 0.997, 0.91, 0.968, 0.909, 0.786, 0.901, 0.995, 0.896, 0.926, 0.986, 0.783, 0.939, 0.869, 0.66, 0.983, 0.683, 0.988, 0.437, 0.735, 0.911, 0.439, 0.462, 0.205, 0.775, 0.592, 0.828, 0.867, 0.974, 0.848, 0.97, 0.927, 0.75, 0.928, 0.965, 0.427, 0.779, 0.53, 0.831, 0.962, 0.739, 0.735, 0.827, 0.767, 0.756, 0.914, 0.611, 0.997, 0.338, 0.924, 0.735, 0.474, 0.878, 0.909, 0.843, 0.984, 0.987, 0.516, 0.999, 0.875, 0.949, 0.766, 0.904, 0.895, 0.807, 0.481, 0.516, 0.696, 0.494, 0.637, 0.93, 0.819, 0.915, 0.999, 0.631, 0.632, 0.859, 0.647, 0.771, 0.325, 0.975, 0.952, 0.161, 0.396, 0.963, 0.518, 0.854, 0.798, 0.544, 0.978, 0.995, 0.901, 0.907, 0.542, 0.688, 0.696, 0.393, 0.732, 0.994, 0.905, 0.703, 0.457, 0.975, 0.839, 0.813, 0.829, 0.472, 0.664, 0.723, 0.485, 0.896, 0.91, 0.9, 0.902, 0.656, 0.968, 0.998, 0.147, 0.403, 0.786, 0.819, 0.344, 0.462, 0.954, 0.991, 0.718, 0.823, 0.592, 0.815, 0.853, 0.491, 0.949, 0.131, 0.891, 0.783, 0.274, 0.503, 0.938, 0.848, 0.975, 0.87, 0.534, 0.822, 0.746, 0.982, 0.82, 0.782, 0.777, 0.545, 0.949, 0.711, 0.781, 0.938, 0.693, 0.999, 0.935, 0.886, 0.787, 0.563, 0.43, 0.721, 0.8, 0.849, 0.926, 0.693, 0.475, 0.993, 0.803, 0.804, 0.869, 0.914, 0.983, 0.993, 0.903, 0.911, 0.735, 0.835, 0.612, 0.968, 0.419, 0.979, 0.798, 0.717, 0.608, 0.995, 0.953, 0.829, 0.572, 0.28, 0.718, 0.774, 0.679, 0.91, 0.944]
global origin = 1
global destination = 50