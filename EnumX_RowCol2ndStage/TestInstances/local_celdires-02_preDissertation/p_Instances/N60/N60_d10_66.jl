global arcs = [1 10; 1 15; 1 35; 1 41; 1 46; 2 11; 2 47; 2 51; 2 56; 2 58; 3 15; 3 22; 3 34; 3 53; 4 3; 4 8; 4 10; 4 20; 4 25; 4 27; 4 33; 4 37; 4 38; 4 40; 4 48; 4 52; 4 58; 5 17; 5 22; 5 48; 6 3; 6 19; 6 20; 6 22; 6 41; 6 50; 7 3; 7 4; 7 6; 7 21; 7 28; 7 30; 7 38; 7 42; 7 60; 8 24; 8 28; 8 36; 8 44; 8 57; 9 5; 9 6; 9 11; 9 15; 9 16; 9 31; 9 35; 9 41; 10 2; 10 6; 10 12; 10 26; 10 28; 10 46; 11 20; 11 28; 11 41; 11 50; 12 4; 12 16; 12 34; 12 45; 13 2; 13 5; 13 23; 13 51; 13 53; 13 55; 14 2; 14 8; 14 9; 14 25; 14 46; 15 3; 15 28; 15 53; 16 4; 16 23; 16 25; 16 31; 16 33; 16 49; 17 4; 17 11; 17 16; 17 25; 17 54; 17 60; 18 8; 18 15; 18 17; 18 28; 18 35; 18 50; 18 54; 19 7; 19 9; 19 10; 19 34; 19 39; 19 48; 19 59; 20 10; 20 21; 20 27; 20 31; 20 33; 20 36; 20 40; 20 53; 21 35; 21 39; 21 42; 21 52; 21 53; 21 59; 22 8; 22 15; 22 51; 22 56; 23 6; 23 10; 23 21; 23 49; 24 3; 24 15; 25 15; 25 31; 25 40; 25 43; 26 20; 26 21; 26 49; 26 51; 26 56; 26 58; 27 12; 27 26; 28 18; 28 29; 28 31; 28 51; 29 4; 29 16; 29 40; 29 55; 30 3; 30 10; 30 18; 30 46; 31 8; 31 9; 31 29; 31 44; 31 59; 31 60; 32 20; 32 27; 32 46; 32 48; 32 54; 32 56; 33 12; 33 16; 33 24; 33 37; 33 40; 33 48; 33 52; 34 15; 34 19; 34 22; 34 24; 34 26; 34 44; 34 47; 34 54; 35 17; 35 22; 35 51; 36 12; 36 21; 36 22; 36 25; 36 47; 37 2; 37 6; 37 28; 37 46; 37 49; 38 23; 38 28; 38 51; 39 3; 39 9; 39 11; 39 26; 39 31; 39 38; 39 42; 39 52; 39 58; 40 4; 40 13; 40 17; 40 25; 40 28; 40 33; 40 36; 40 53; 40 57; 40 59; 41 6; 41 7; 41 56; 42 14; 42 23; 42 24; 42 28; 42 40; 42 44; 42 45; 42 58; 43 25; 43 29; 43 31; 43 52; 43 56; 44 29; 44 32; 44 38; 44 47; 44 52; 44 55; 45 20; 45 23; 45 27; 45 59; 46 11; 46 12; 46 16; 46 20; 47 5; 47 14; 47 24; 47 33; 47 38; 47 52; 48 12; 48 19; 48 33; 48 41; 48 50; 48 58; 49 4; 49 13; 49 19; 49 22; 49 29; 49 35; 49 36; 49 38; 49 40; 49 50; 49 54; 50 7; 50 23; 50 25; 50 26; 50 39; 50 41; 50 42; 50 59; 51 11; 51 12; 51 16; 51 22; 51 36; 52 4; 52 6; 52 22; 52 53; 52 54; 52 60; 53 3; 53 7; 53 10; 53 28; 53 33; 53 37; 53 38; 53 40; 53 42; 53 50; 53 57; 54 15; 54 29; 54 39; 54 43; 54 49; 55 2; 55 6; 55 8; 55 21; 55 32; 55 37; 55 56; 56 15; 56 17; 56 20; 56 23; 56 27; 56 28; 56 41; 56 42; 56 55; 57 26; 57 30; 57 46; 57 51; 58 6; 58 11; 58 29; 58 30; 58 34; 58 37; 58 49; 58 54; 58 55; 58 57; 59 18; 59 31; 59 34; 59 36; 59 41; 59 44; 59 60]
global d_x = [4.0, 2.0, 6.0, 2.0, 10.0, 4.0, 7.0, 6.0, 8.0, 6.0, 4.0, 2.0, 4.0, 5.0, 6.0, 7.0, 1.0, 7.0, 6.0, 5.0, 5.0, 1.0, 10.0, 1.0, 5.0, 2.0, 7.0, 1.0, 3.0, 4.0, 4.0, 3.0, 3.0, 3.0, 6.0, 7.0, 4.0, 1.0, 5.0, 1.0, 5.0, 6.0, 8.0, 7.0, 1.0, 5.0, 9.0, 8.0, 7.0, 6.0, 6.0, 2.0, 10.0, 9.0, 2.0, 8.0, 1.0, 2.0, 10.0, 9.0, 3.0, 5.0, 9.0, 2.0, 6.0, 10.0, 2.0, 3.0, 6.0, 2.0, 4.0, 9.0, 10.0, 4.0, 5.0, 3.0, 1.0, 6.0, 7.0, 4.0, 2.0, 4.0, 9.0, 9.0, 2.0, 5.0, 10.0, 7.0, 7.0, 7.0, 9.0, 6.0, 9.0, 9.0, 3.0, 6.0, 5.0, 10.0, 9.0, 3.0, 3.0, 4.0, 6.0, 3.0, 1.0, 8.0, 6.0, 3.0, 8.0, 7.0, 1.0, 4.0, 1.0, 9.0, 1.0, 4.0, 7.0, 5.0, 1.0, 6.0, 1.0, 4.0, 7.0, 10.0, 2.0, 4.0, 7.0, 7.0, 1.0, 7.0, 7.0, 1.0, 7.0, 5.0, 10.0, 7.0, 2.0, 5.0, 6.0, 5.0, 10.0, 8.0, 1.0, 4.0, 2.0, 9.0, 5.0, 7.0, 1.0, 4.0, 5.0, 6.0, 6.0, 10.0, 1.0, 7.0, 7.0, 8.0, 9.0, 1.0, 10.0, 6.0, 10.0, 8.0, 2.0, 3.0, 3.0, 3.0, 7.0, 1.0, 6.0, 7.0, 6.0, 1.0, 3.0, 8.0, 5.0, 5.0, 1.0, 2.0, 6.0, 2.0, 1.0, 1.0, 8.0, 6.0, 2.0, 4.0, 9.0, 5.0, 8.0, 5.0, 3.0, 5.0, 7.0, 10.0, 1.0, 2.0, 3.0, 8.0, 8.0, 3.0, 8.0, 8.0, 8.0, 10.0, 1.0, 3.0, 7.0, 3.0, 6.0, 4.0, 8.0, 6.0, 7.0, 3.0, 1.0, 2.0, 2.0, 9.0, 8.0, 7.0, 3.0, 9.0, 6.0, 5.0, 7.0, 8.0, 3.0, 7.0, 6.0, 8.0, 8.0, 5.0, 2.0, 10.0, 9.0, 4.0, 2.0, 8.0, 2.0, 4.0, 4.0, 5.0, 4.0, 6.0, 5.0, 9.0, 9.0, 10.0, 9.0, 3.0, 7.0, 5.0, 4.0, 9.0, 2.0, 8.0, 2.0, 4.0, 4.0, 8.0, 3.0, 9.0, 9.0, 1.0, 2.0, 2.0, 5.0, 1.0, 1.0, 6.0, 10.0, 7.0, 4.0, 6.0, 7.0, 6.0, 9.0, 4.0, 9.0, 5.0, 9.0, 4.0, 10.0, 8.0, 6.0, 8.0, 5.0, 7.0, 5.0, 10.0, 10.0, 10.0, 1.0, 9.0, 8.0, 8.0, 1.0, 5.0, 7.0, 4.0, 7.0, 7.0, 4.0, 4.0, 4.0, 3.0, 5.0, 8.0, 5.0, 3.0, 4.0, 3.0, 7.0, 9.0, 5.0, 4.0, 10.0, 2.0, 9.0, 9.0, 9.0, 1.0, 2.0, 6.0, 7.0, 10.0, 9.0, 4.0, 2.0, 9.0, 1.0, 3.0, 5.0, 2.0, 3.0, 6.0, 10.0, 4.0, 7.0, 10.0, 3.0, 7.0, 1.0, 10.0, 10.0]
global b_x = 5
global d_y = [8.0, 8.0, 5.0, 5.0, 8.0, 9.0, 2.0, 7.0, 1.0, 2.0, 2.0, 8.0, 10.0, 5.0, 7.0, 4.0, 2.0, 7.0, 1.0, 1.0, 6.0, 9.0, 2.0, 7.0, 8.0, 8.0, 7.0, 2.0, 4.0, 10.0, 4.0, 1.0, 1.0, 2.0, 3.0, 8.0, 4.0, 5.0, 5.0, 5.0, 9.0, 3.0, 1.0, 1.0, 3.0, 10.0, 8.0, 2.0, 6.0, 9.0, 8.0, 8.0, 7.0, 5.0, 8.0, 1.0, 2.0, 7.0, 5.0, 4.0, 7.0, 1.0, 9.0, 7.0, 9.0, 10.0, 10.0, 10.0, 10.0, 10.0, 8.0, 8.0, 5.0, 5.0, 2.0, 10.0, 8.0, 4.0, 9.0, 1.0, 5.0, 1.0, 3.0, 1.0, 6.0, 10.0, 10.0, 6.0, 2.0, 9.0, 2.0, 4.0, 3.0, 3.0, 3.0, 5.0, 3.0, 9.0, 10.0, 5.0, 9.0, 4.0, 3.0, 1.0, 7.0, 8.0, 7.0, 9.0, 8.0, 8.0, 7.0, 9.0, 6.0, 2.0, 7.0, 4.0, 6.0, 2.0, 6.0, 7.0, 1.0, 2.0, 8.0, 6.0, 8.0, 8.0, 2.0, 2.0, 1.0, 1.0, 3.0, 5.0, 7.0, 10.0, 4.0, 10.0, 6.0, 4.0, 2.0, 5.0, 1.0, 1.0, 2.0, 8.0, 2.0, 2.0, 1.0, 1.0, 3.0, 5.0, 7.0, 5.0, 4.0, 3.0, 9.0, 2.0, 8.0, 10.0, 9.0, 10.0, 9.0, 1.0, 2.0, 8.0, 4.0, 1.0, 9.0, 10.0, 5.0, 10.0, 8.0, 7.0, 3.0, 5.0, 2.0, 7.0, 3.0, 10.0, 10.0, 5.0, 3.0, 2.0, 4.0, 5.0, 7.0, 3.0, 3.0, 8.0, 6.0, 6.0, 3.0, 6.0, 10.0, 5.0, 9.0, 10.0, 6.0, 4.0, 4.0, 1.0, 6.0, 7.0, 2.0, 9.0, 8.0, 3.0, 8.0, 8.0, 1.0, 6.0, 5.0, 8.0, 7.0, 10.0, 9.0, 8.0, 10.0, 1.0, 2.0, 1.0, 8.0, 9.0, 8.0, 8.0, 7.0, 4.0, 1.0, 5.0, 6.0, 2.0, 9.0, 7.0, 7.0, 9.0, 8.0, 1.0, 7.0, 8.0, 9.0, 3.0, 4.0, 9.0, 10.0, 5.0, 1.0, 4.0, 3.0, 2.0, 6.0, 1.0, 4.0, 1.0, 6.0, 8.0, 5.0, 1.0, 7.0, 6.0, 5.0, 10.0, 9.0, 3.0, 7.0, 4.0, 8.0, 7.0, 7.0, 5.0, 5.0, 9.0, 7.0, 5.0, 6.0, 3.0, 4.0, 10.0, 7.0, 3.0, 7.0, 4.0, 8.0, 4.0, 4.0, 3.0, 3.0, 10.0, 6.0, 7.0, 5.0, 10.0, 8.0, 7.0, 6.0, 10.0, 5.0, 1.0, 8.0, 9.0, 7.0, 3.0, 7.0, 9.0, 8.0, 6.0, 2.0, 7.0, 10.0, 2.0, 3.0, 8.0, 6.0, 2.0, 6.0, 3.0, 1.0, 4.0, 5.0, 2.0, 6.0, 8.0, 10.0, 3.0, 7.0, 2.0, 7.0, 6.0, 9.0, 7.0, 3.0, 4.0, 3.0, 10.0, 6.0, 2.0, 10.0, 7.0, 5.0, 6.0, 9.0, 1.0, 10.0, 3.0, 3.0, 2.0, 6.0, 6.0, 5.0]
global b_y = 10
global p = [0.907, 0.393, 0.839, 0.758, 0.805, 0.743, 0.884, 0.016, 0.487, 0.098, 0.19, 0.605, 0.015, 0.911, 0.867, 0.145, 0.845, 0.057, 0.865, 0.338, 0.327, 0.701, 0.849, 0.9, 0.137, 0.583, 0.852, 0.407, 0.713, 0.447, 0.382, 0.345, 0.832, 0.169, 0.166, 0.65, 0.444, 0.971, 0.718, 0.151, 0.472, 0.915, 0.608, 0.87, 0.593, 0.105, 0.818, 0.241, 0.69, 0.513, 0.744, 0.73, 0.433, 0.484, 0.429, 0.566, 0.915, 0.742, 0.957, 0.014, 0.053, 0.401, 0.734, 0.161, 0.572, 0.151, 0.115, 0.556, 0.615, 0.933, 0.2, 0.805, 0.978, 0.829, 0.104, 0.822, 0.515, 0.916, 0.353, 0.121, 0.571, 0.826, 0.693, 0.568, 0.514, 0.685, 0.486, 0.801, 0.689, 0.187, 0.527, 0.724, 0.57, 0.013, 0.576, 0.699, 0.513, 0.265, 0.62, 0.532, 0.369, 0.554, 0.174, 0.983, 0.861, 0.582, 0.078, 0.083, 0.554, 0.945, 0.789, 0.595, 0.594, 0.737, 0.343, 0.24, 0.806, 0.618, 0.266, 0.469, 0.24, 0.369, 0.022, 0.228, 0.349, 0.582, 0.024, 0.219, 0.012, 0.722, 0.806, 0.299, 0.452, 0.772, 0.692, 0.276, 0.613, 0.68, 0.796, 0.909, 0.341, 0.385, 0.137, 0.695, 0.55, 0.115, 0.705, 0.849, 0.495, 0.47, 0.935, 0.2, 0.237, 0.356, 0.497, 0.93, 0.815, 0.637, 0.184, 0.764, 0.257, 0.285, 0.44, 0.21, 0.661, 0.192, 0.645, 0.372, 0.36, 0.073, 0.88, 0.417, 0.537, 0.61, 0.723, 0.139, 0.922, 0.838, 0.39, 0.034, 0.636, 0.92, 0.234, 0.169, 0.352, 0.962, 0.387, 0.266, 0.022, 0.914, 0.643, 0.55, 0.443, 0.337, 0.359, 0.396, 0.622, 0.786, 0.205, 0.647, 0.009, 0.911, 0.953, 0.659, 0.468, 0.499, 0.978, 0.599, 0.268, 0.505, 0.865, 0.056, 0.849, 0.477, 0.438, 0.519, 0.164, 0.196, 0.562, 0.847, 0.137, 0.985, 0.827, 0.966, 0.705, 0.013, 0.015, 0.279, 0.803, 0.512, 0.771, 0.703, 0.665, 0.493, 0.968, 0.069, 0.374, 0.889, 0.139, 0.371, 0.163, 0.084, 0.373, 0.533, 0.73, 0.353, 0.996, 0.94, 0.141, 0.159, 0.984, 0.916, 0.362, 0.461, 0.455, 0.972, 0.122, 0.672, 0.171, 0.738, 0.791, 0.474, 0.748, 0.289, 0.73, 0.453, 0.52, 0.032, 0.728, 0.423, 0.089, 0.635, 0.346, 0.055, 0.107, 0.116, 0.238, 0.995, 0.262, 0.04, 0.472, 0.175, 0.844, 0.215, 0.497, 0.29, 0.66, 0.43, 0.538, 0.973, 0.777, 0.98, 0.572, 0.347, 0.536, 0.847, 0.873, 0.883, 0.626, 0.021, 0.978, 0.162, 0.452, 0.007, 0.477, 0.703, 0.093, 0.221, 0.112, 0.842, 0.391, 0.444, 0.315, 0.237, 0.052, 0.437, 0.555, 0.287, 0.695, 0.089, 0.263, 0.169, 0.585, 0.447, 0.456, 0.9, 0.675, 0.717, 0.06, 0.057, 0.769, 0.545, 0.36, 0.004, 0.342, 0.972, 0.088, 0.452, 0.174, 0.909, 0.787, 0.661, 0.633, 0.073, 0.169, 0.014, 0.79]
global q = [0.935, 0.662, 0.844, 0.877, 0.974, 0.889, 0.934, 0.503, 0.53, 0.911, 0.881, 0.862, 0.891, 0.93, 0.989, 0.33, 0.922, 0.123, 0.872, 0.614, 0.839, 0.751, 0.945, 0.95, 0.29, 0.597, 0.933, 0.457, 0.807, 0.807, 0.634, 0.423, 0.896, 0.383, 0.582, 0.934, 0.786, 0.981, 0.958, 0.256, 0.903, 0.939, 0.708, 0.9, 0.905, 0.816, 0.909, 0.306, 0.784, 0.711, 0.958, 0.761, 0.553, 0.695, 0.775, 0.67, 0.963, 0.967, 0.976, 0.112, 0.227, 0.635, 0.798, 0.535, 0.641, 0.167, 0.265, 0.6, 0.993, 0.971, 0.325, 0.935, 0.981, 0.917, 0.951, 0.925, 0.583, 0.939, 0.895, 0.517, 0.61, 0.866, 0.966, 0.584, 0.559, 0.909, 0.991, 0.812, 0.792, 0.533, 0.71, 0.857, 0.697, 0.929, 0.838, 0.817, 0.993, 0.986, 0.731, 0.632, 0.729, 0.907, 0.652, 0.989, 0.899, 0.706, 0.763, 0.201, 0.65, 0.976, 0.8, 0.742, 0.888, 0.765, 0.525, 0.321, 0.982, 0.954, 0.833, 0.598, 0.613, 0.832, 0.159, 0.373, 0.965, 0.944, 0.476, 0.357, 0.851, 0.834, 0.937, 0.815, 0.874, 0.822, 0.968, 0.279, 0.667, 0.94, 0.909, 0.931, 0.479, 0.934, 0.635, 0.722, 0.904, 0.161, 0.882, 0.949, 0.624, 0.811, 0.991, 0.485, 0.362, 0.899, 0.716, 0.943, 0.838, 0.82, 0.653, 0.869, 0.865, 0.288, 0.973, 0.27, 0.669, 0.577, 0.848, 0.903, 0.577, 0.155, 0.941, 0.89, 0.605, 0.612, 0.876, 0.163, 0.95, 0.957, 0.521, 0.748, 0.927, 0.955, 0.764, 0.175, 0.795, 0.968, 0.999, 0.327, 0.765, 0.949, 0.95, 0.862, 0.57, 0.732, 0.742, 0.82, 0.801, 0.982, 0.324, 0.985, 0.85, 0.917, 0.967, 0.848, 0.86, 0.638, 0.989, 0.953, 0.607, 0.902, 0.901, 0.36, 0.934, 0.665, 0.701, 0.945, 0.61, 0.907, 0.609, 0.996, 0.536, 0.985, 0.997, 0.985, 0.896, 0.851, 0.567, 0.747, 0.932, 0.833, 0.875, 0.926, 0.951, 0.571, 0.999, 0.883, 0.882, 0.924, 0.822, 0.954, 0.171, 0.299, 0.929, 0.578, 0.79, 0.649, 0.997, 0.986, 0.307, 0.878, 0.985, 0.945, 0.906, 0.519, 0.635, 0.983, 0.252, 0.797, 0.959, 0.826, 0.817, 0.724, 0.785, 0.494, 0.989, 0.957, 0.93, 0.814, 0.854, 0.545, 0.116, 0.819, 0.792, 0.124, 0.2, 0.882, 0.81, 0.999, 0.634, 0.251, 0.929, 0.601, 0.903, 0.985, 0.686, 0.518, 0.737, 0.559, 0.561, 0.991, 0.911, 0.984, 0.65, 0.853, 0.819, 0.977, 0.947, 0.901, 0.886, 0.076, 0.98, 0.369, 0.668, 0.728, 0.585, 0.703, 0.538, 0.431, 0.851, 0.998, 0.669, 0.686, 0.969, 0.36, 0.745, 0.457, 0.998, 0.771, 0.699, 0.769, 0.861, 0.244, 0.731, 0.77, 0.867, 0.975, 0.784, 0.982, 0.398, 0.088, 0.785, 0.588, 0.589, 0.327, 0.595, 0.986, 0.443, 0.718, 0.547, 0.959, 0.902, 0.766, 0.708, 0.831, 0.967, 0.052, 0.959]
global origin = 1
global destination = 60