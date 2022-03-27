global arcs = [1 9; 1 18; 1 23; 1 29; 1 30; 1 46; 2 4; 2 16; 2 25; 2 32; 3 6; 3 32; 3 37; 3 41; 3 43; 4 8; 4 26; 4 37; 5 10; 5 13; 5 19; 5 38; 5 39; 5 41; 6 8; 6 9; 6 35; 6 38; 6 46; 7 6; 7 11; 7 27; 7 28; 8 5; 8 26; 8 32; 8 43; 9 3; 9 20; 9 25; 9 37; 9 39; 10 5; 10 29; 10 36; 10 37; 10 46; 10 50; 11 9; 11 19; 11 23; 11 32; 11 42; 12 7; 12 24; 12 26; 12 37; 13 9; 13 31; 13 42; 13 47; 13 48; 14 21; 14 43; 15 11; 15 17; 15 19; 15 23; 16 27; 16 30; 16 46; 17 3; 17 6; 17 13; 17 20; 17 27; 17 30; 17 35; 17 43; 18 16; 18 17; 18 19; 18 21; 18 30; 18 39; 18 41; 19 17; 19 46; 20 18; 20 38; 20 41; 21 6; 21 14; 21 33; 21 38; 21 44; 21 45; 22 16; 22 24; 22 50; 23 8; 23 9; 23 15; 23 17; 23 20; 23 22; 23 30; 23 33; 24 4; 24 5; 24 18; 25 22; 25 42; 25 49; 26 16; 26 35; 26 38; 26 48; 27 10; 27 16; 27 21; 27 33; 27 35; 27 44; 27 48; 28 9; 28 14; 28 20; 29 18; 29 22; 29 23; 29 27; 29 37; 29 49; 30 8; 30 24; 30 40; 31 8; 31 12; 31 35; 31 47; 31 49; 32 33; 32 34; 32 49; 33 2; 33 4; 33 12; 33 13; 33 14; 33 22; 33 35; 33 36; 33 47; 33 48; 34 8; 34 14; 34 18; 34 39; 35 7; 35 39; 35 41; 35 48; 35 49; 36 16; 36 32; 37 5; 37 13; 37 24; 37 31; 37 38; 37 47; 37 50; 38 2; 38 4; 38 5; 38 40; 38 47; 39 5; 39 19; 39 23; 39 43; 39 47; 39 48; 40 28; 40 30; 40 37; 41 4; 41 21; 41 35; 41 44; 42 5; 42 9; 42 25; 42 33; 42 35; 42 36; 42 39; 42 40; 43 9; 43 25; 43 36; 44 7; 44 20; 44 25; 44 37; 44 40; 44 42; 44 47; 45 9; 45 23; 45 24; 45 25; 45 27; 45 43; 45 50; 46 10; 46 22; 47 2; 47 6; 47 7; 47 18; 47 20; 47 48; 48 9; 48 12; 48 15; 48 23; 48 42; 48 50; 49 17; 49 33; 49 34; 49 40; 49 45; 49 50]
global d_x = [8.0, 7.0, 5.0, 8.0, 5.0, 1.0, 9.0, 10.0, 6.0, 3.0, 3.0, 5.0, 3.0, 4.0, 3.0, 9.0, 8.0, 2.0, 1.0, 1.0, 10.0, 3.0, 5.0, 5.0, 10.0, 6.0, 4.0, 3.0, 9.0, 10.0, 6.0, 7.0, 5.0, 3.0, 2.0, 3.0, 7.0, 7.0, 1.0, 5.0, 2.0, 1.0, 9.0, 2.0, 1.0, 7.0, 6.0, 8.0, 8.0, 3.0, 10.0, 8.0, 5.0, 1.0, 6.0, 10.0, 8.0, 9.0, 6.0, 9.0, 4.0, 6.0, 5.0, 2.0, 1.0, 8.0, 8.0, 6.0, 5.0, 5.0, 10.0, 5.0, 2.0, 10.0, 3.0, 6.0, 2.0, 3.0, 2.0, 7.0, 9.0, 4.0, 7.0, 10.0, 2.0, 9.0, 7.0, 2.0, 8.0, 1.0, 4.0, 4.0, 1.0, 1.0, 2.0, 6.0, 9.0, 7.0, 3.0, 8.0, 8.0, 3.0, 7.0, 3.0, 7.0, 8.0, 1.0, 6.0, 9.0, 8.0, 10.0, 3.0, 8.0, 4.0, 7.0, 9.0, 9.0, 8.0, 5.0, 7.0, 6.0, 4.0, 7.0, 10.0, 8.0, 7.0, 10.0, 8.0, 2.0, 7.0, 3.0, 10.0, 5.0, 7.0, 6.0, 1.0, 3.0, 4.0, 5.0, 10.0, 3.0, 9.0, 8.0, 10.0, 2.0, 9.0, 3.0, 10.0, 10.0, 5.0, 8.0, 8.0, 10.0, 7.0, 9.0, 4.0, 7.0, 5.0, 8.0, 4.0, 3.0, 5.0, 9.0, 6.0, 2.0, 1.0, 5.0, 1.0, 9.0, 5.0, 1.0, 3.0, 4.0, 10.0, 8.0, 10.0, 5.0, 2.0, 4.0, 7.0, 10.0, 7.0, 3.0, 7.0, 2.0, 8.0, 2.0, 7.0, 9.0, 5.0, 7.0, 1.0, 5.0, 8.0, 5.0, 2.0, 9.0, 10.0, 6.0, 3.0, 8.0, 2.0, 2.0, 9.0, 1.0, 10.0, 2.0, 3.0, 6.0, 8.0, 2.0, 4.0, 7.0, 9.0, 5.0, 4.0, 3.0, 10.0, 5.0, 2.0, 2.0, 10.0, 10.0, 4.0, 3.0, 2.0, 3.0, 6.0, 4.0, 1.0, 10.0, 6.0, 6.0, 6.0, 10.0, 8.0]
global b_x = 5
global d_y = [6.0, 1.0, 7.0, 6.0, 8.0, 2.0, 2.0, 4.0, 9.0, 6.0, 3.0, 2.0, 7.0, 3.0, 4.0, 7.0, 3.0, 3.0, 6.0, 4.0, 5.0, 3.0, 3.0, 5.0, 8.0, 1.0, 8.0, 6.0, 4.0, 5.0, 5.0, 7.0, 5.0, 5.0, 7.0, 6.0, 3.0, 1.0, 6.0, 6.0, 8.0, 2.0, 4.0, 2.0, 7.0, 3.0, 3.0, 10.0, 1.0, 5.0, 7.0, 2.0, 2.0, 7.0, 1.0, 7.0, 8.0, 5.0, 7.0, 9.0, 2.0, 1.0, 7.0, 7.0, 6.0, 9.0, 6.0, 6.0, 7.0, 6.0, 4.0, 7.0, 3.0, 7.0, 5.0, 1.0, 1.0, 9.0, 2.0, 5.0, 6.0, 5.0, 1.0, 6.0, 6.0, 8.0, 6.0, 8.0, 5.0, 7.0, 4.0, 3.0, 6.0, 6.0, 3.0, 1.0, 1.0, 10.0, 6.0, 6.0, 5.0, 6.0, 7.0, 9.0, 6.0, 8.0, 1.0, 2.0, 10.0, 1.0, 10.0, 9.0, 7.0, 3.0, 2.0, 8.0, 4.0, 1.0, 6.0, 4.0, 6.0, 5.0, 10.0, 9.0, 6.0, 3.0, 10.0, 2.0, 2.0, 8.0, 3.0, 8.0, 3.0, 2.0, 10.0, 5.0, 9.0, 1.0, 8.0, 4.0, 5.0, 10.0, 5.0, 9.0, 1.0, 5.0, 10.0, 10.0, 6.0, 8.0, 3.0, 8.0, 5.0, 4.0, 6.0, 4.0, 9.0, 3.0, 10.0, 2.0, 4.0, 3.0, 6.0, 1.0, 7.0, 3.0, 10.0, 7.0, 4.0, 6.0, 2.0, 9.0, 4.0, 7.0, 3.0, 1.0, 3.0, 9.0, 10.0, 8.0, 1.0, 10.0, 5.0, 6.0, 4.0, 4.0, 10.0, 5.0, 2.0, 5.0, 7.0, 4.0, 7.0, 4.0, 2.0, 8.0, 2.0, 1.0, 8.0, 10.0, 8.0, 6.0, 1.0, 6.0, 8.0, 5.0, 5.0, 10.0, 9.0, 10.0, 8.0, 2.0, 9.0, 6.0, 3.0, 8.0, 1.0, 7.0, 8.0, 4.0, 6.0, 1.0, 9.0, 9.0, 9.0, 1.0, 10.0, 9.0, 7.0, 1.0, 1.0, 9.0, 4.0, 9.0, 4.0, 5.0]
global b_y = 10
global p = [0.234, 0.835, 0.672, 0.521, 0.089, 0.292, 0.582, 0.916, 0.544, 0.974, 0.23, 0.387, 0.624, 0.897, 0.867, 0.574, 0.858, 0.729, 0.895, 0.631, 0.429, 0.169, 0.67, 0.332, 0.566, 0.75, 0.706, 0.425, 0.559, 0.856, 0.049, 0.308, 0.185, 0.915, 0.239, 0.04, 0.915, 0.576, 0.361, 0.642, 0.034, 0.004, 0.643, 0.297, 0.706, 0.943, 0.676, 0.805, 0.078, 0.83, 0.746, 0.495, 0.724, 0.752, 0.51, 0.079, 0.853, 0.496, 0.538, 0.914, 0.906, 0.42, 0.493, 0.769, 0.564, 0.337, 0.253, 0.395, 0.004, 0.358, 0.036, 0.345, 0.711, 0.594, 0.642, 0.438, 0.664, 0.737, 0.338, 0.699, 0.679, 0.756, 0.318, 0.77, 0.354, 0.524, 0.115, 0.923, 0.641, 0.634, 0.242, 0.879, 0.112, 0.28, 0.513, 0.297, 0.557, 0.974, 0.604, 0.338, 0.371, 0.381, 0.962, 0.919, 0.481, 0.968, 0.77, 0.317, 0.953, 0.563, 0.825, 0.77, 0.856, 0.58, 0.773, 0.293, 0.248, 0.048, 0.116, 0.9, 0.688, 0.591, 0.42, 0.345, 0.698, 0.391, 0.107, 0.317, 0.294, 0.707, 0.382, 0.952, 0.419, 0.465, 0.186, 0.778, 0.033, 0.459, 0.653, 0.698, 0.506, 0.171, 0.768, 0.369, 0.024, 0.121, 0.163, 0.169, 0.929, 0.392, 0.126, 0.212, 0.998, 0.456, 0.824, 0.711, 0.939, 0.844, 0.092, 0.932, 0.133, 0.675, 0.014, 0.883, 0.407, 0.074, 0.747, 0.749, 0.866, 0.966, 0.348, 0.567, 0.608, 0.781, 0.926, 0.545, 0.287, 0.055, 0.369, 0.905, 0.798, 0.192, 0.469, 0.2, 0.574, 0.471, 0.043, 0.85, 0.987, 0.474, 0.442, 0.194, 0.53, 0.84, 0.113, 0.343, 0.882, 0.85, 0.816, 0.828, 0.899, 0.464, 0.931, 0.04, 0.06, 0.997, 0.958, 0.674, 0.734, 0.67, 0.061, 0.881, 0.562, 0.089, 0.269, 0.387, 0.172, 0.777, 0.551, 0.498, 0.635, 0.305, 0.97, 0.12, 0.698, 0.51, 0.594, 0.329, 0.127, 0.151, 0.631, 0.229, 0.546, 0.362, 0.638, 0.608]
global q = [0.27, 0.91, 0.845, 0.659, 0.657, 0.838, 0.994, 0.944, 0.609, 0.992, 0.34, 0.53, 0.637, 0.906, 0.926, 0.646, 0.885, 0.972, 0.978, 0.718, 0.911, 0.306, 0.982, 0.53, 0.824, 0.852, 0.802, 0.815, 0.903, 0.887, 0.386, 0.682, 0.364, 0.919, 0.28, 0.341, 0.947, 0.576, 0.956, 0.98, 0.742, 0.526, 0.697, 0.991, 0.912, 0.951, 0.832, 0.909, 0.94, 0.946, 0.8, 0.899, 0.739, 0.902, 0.981, 0.445, 0.884, 0.623, 0.7, 0.932, 0.932, 0.758, 0.839, 0.912, 0.746, 0.697, 0.357, 0.407, 0.06, 0.775, 0.716, 0.4, 0.714, 0.916, 0.91, 0.846, 0.854, 0.884, 0.759, 0.984, 0.774, 0.898, 0.884, 0.946, 0.394, 0.861, 0.173, 0.979, 0.661, 0.907, 0.412, 0.976, 0.44, 0.894, 0.664, 0.856, 0.963, 0.997, 0.639, 0.493, 0.541, 0.427, 0.987, 0.93, 0.836, 0.982, 0.995, 0.958, 0.969, 0.704, 0.941, 0.928, 0.978, 0.754, 0.82, 0.769, 0.284, 0.186, 0.124, 0.913, 0.957, 0.793, 0.948, 0.489, 0.749, 0.826, 0.995, 0.668, 0.975, 0.922, 0.554, 0.965, 0.512, 0.648, 0.444, 0.904, 0.218, 0.473, 0.81, 0.878, 0.696, 0.578, 0.936, 0.895, 0.393, 0.122, 0.682, 0.823, 0.95, 0.95, 0.739, 0.296, 0.999, 0.692, 0.932, 0.952, 0.953, 0.96, 0.189, 0.957, 0.436, 0.884, 0.138, 0.94, 0.687, 0.479, 0.975, 0.91, 0.896, 0.974, 0.558, 0.956, 0.626, 0.819, 0.969, 0.696, 0.725, 0.604, 0.369, 0.987, 0.808, 0.222, 0.922, 0.448, 0.801, 0.89, 0.129, 0.99, 0.994, 0.518, 0.532, 0.229, 0.616, 0.84, 0.902, 0.76, 0.911, 0.986, 0.978, 0.867, 0.953, 0.535, 0.949, 0.3, 0.31, 0.998, 0.984, 0.869, 0.767, 0.854, 0.535, 0.891, 0.692, 0.308, 0.752, 0.71, 0.93, 0.962, 0.554, 0.626, 0.963, 0.847, 0.981, 0.284, 0.863, 0.615, 0.911, 0.345, 0.234, 0.292, 0.915, 0.547, 0.901, 0.964, 0.817, 0.664]
global origin = 1
global destination = 50