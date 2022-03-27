global arcs = [1 15; 1 22; 1 33; 1 51; 2 5; 2 6; 2 10; 2 11; 2 32; 2 38; 2 40; 3 4; 3 13; 3 42; 4 7; 4 12; 4 19; 4 25; 4 31; 4 34; 4 36; 4 37; 4 39; 4 47; 4 50; 4 53; 5 6; 5 10; 5 11; 5 35; 5 38; 5 47; 5 48; 5 50; 5 60; 6 19; 6 26; 6 39; 6 49; 7 4; 7 14; 7 15; 7 27; 7 43; 7 54; 8 2; 8 5; 8 27; 8 49; 9 20; 9 26; 9 33; 9 34; 9 49; 9 50; 9 57; 9 58; 10 20; 10 26; 10 37; 10 52; 10 55; 10 58; 10 59; 11 3; 11 6; 11 41; 11 42; 11 52; 12 3; 12 4; 12 11; 12 28; 12 29; 12 31; 12 45; 13 5; 13 7; 13 10; 13 25; 13 30; 13 34; 13 41; 13 42; 13 50; 13 56; 13 57; 14 30; 14 33; 14 40; 15 31; 15 58; 16 10; 16 39; 16 43; 16 47; 17 2; 17 3; 17 6; 17 11; 17 18; 17 31; 17 37; 17 44; 17 46; 17 47; 17 56; 18 2; 18 13; 18 26; 18 41; 19 3; 19 21; 19 24; 19 30; 19 32; 19 39; 19 50; 20 6; 20 33; 20 49; 20 60; 21 25; 21 32; 21 53; 22 2; 22 9; 22 17; 22 32; 22 35; 22 36; 22 37; 22 44; 22 48; 22 54; 23 14; 23 19; 23 25; 23 26; 23 32; 23 53; 24 8; 24 13; 24 14; 24 19; 24 41; 24 56; 24 57; 25 17; 25 23; 25 41; 26 20; 26 34; 26 36; 26 39; 27 19; 27 20; 27 31; 27 39; 27 41; 28 7; 28 26; 28 32; 28 37; 28 44; 28 54; 28 56; 28 59; 29 8; 29 51; 30 6; 30 8; 30 35; 30 40; 30 48; 30 54; 30 59; 31 11; 31 40; 31 50; 32 36; 32 45; 32 49; 32 55; 33 4; 33 7; 33 13; 33 15; 33 48; 34 5; 34 9; 34 10; 34 23; 34 27; 34 30; 34 31; 34 42; 34 58; 35 15; 35 17; 35 25; 35 29; 35 40; 35 45; 35 50; 35 58; 35 59; 36 6; 36 10; 36 21; 36 24; 36 31; 36 38; 36 55; 36 59; 37 22; 37 25; 37 27; 37 39; 37 42; 37 43; 37 50; 38 4; 38 24; 38 29; 38 42; 38 47; 39 4; 39 9; 39 10; 39 21; 39 45; 39 48; 39 56; 40 4; 40 6; 40 15; 40 55; 41 7; 41 13; 41 21; 41 25; 41 33; 41 46; 42 2; 42 8; 42 19; 42 31; 42 35; 42 38; 42 50; 42 51; 43 4; 43 27; 43 49; 43 60; 44 4; 44 15; 44 16; 44 27; 44 33; 44 50; 44 59; 45 12; 45 23; 45 24; 45 27; 45 51; 45 55; 46 21; 46 27; 46 37; 46 38; 46 48; 46 49; 47 21; 47 22; 47 33; 47 34; 47 35; 48 4; 48 12; 48 17; 48 28; 48 34; 48 50; 49 15; 49 38; 49 40; 49 55; 50 6; 50 12; 50 15; 50 21; 50 38; 51 2; 51 15; 51 17; 51 43; 51 46; 52 2; 52 14; 52 16; 52 44; 52 46; 52 56; 53 2; 53 4; 53 12; 53 22; 53 30; 53 40; 53 43; 54 39; 54 46; 54 53; 55 8; 55 21; 55 25; 55 31; 56 3; 56 17; 56 18; 56 24; 56 25; 57 9; 57 18; 57 45; 58 19; 58 34; 58 46; 58 53; 58 55; 58 59; 59 2; 59 4; 59 16; 59 25; 59 33; 59 35; 59 36]
global d_x = [10.0, 10.0, 2.0, 8.0, 6.0, 5.0, 10.0, 7.0, 8.0, 7.0, 6.0, 9.0, 4.0, 2.0, 8.0, 9.0, 1.0, 9.0, 5.0, 4.0, 9.0, 8.0, 4.0, 7.0, 8.0, 6.0, 3.0, 5.0, 10.0, 7.0, 3.0, 1.0, 4.0, 8.0, 6.0, 10.0, 9.0, 1.0, 2.0, 2.0, 6.0, 2.0, 1.0, 4.0, 5.0, 1.0, 1.0, 2.0, 2.0, 7.0, 7.0, 8.0, 6.0, 2.0, 1.0, 8.0, 2.0, 9.0, 1.0, 10.0, 4.0, 10.0, 9.0, 8.0, 8.0, 4.0, 5.0, 7.0, 8.0, 5.0, 1.0, 2.0, 10.0, 9.0, 10.0, 1.0, 10.0, 7.0, 5.0, 6.0, 7.0, 2.0, 7.0, 2.0, 3.0, 1.0, 4.0, 4.0, 9.0, 6.0, 9.0, 3.0, 1.0, 5.0, 6.0, 8.0, 3.0, 5.0, 5.0, 4.0, 3.0, 9.0, 4.0, 7.0, 2.0, 7.0, 9.0, 8.0, 1.0, 10.0, 5.0, 3.0, 8.0, 3.0, 8.0, 2.0, 7.0, 6.0, 4.0, 1.0, 9.0, 2.0, 8.0, 8.0, 7.0, 10.0, 10.0, 7.0, 1.0, 7.0, 4.0, 9.0, 5.0, 5.0, 3.0, 9.0, 2.0, 5.0, 8.0, 6.0, 10.0, 2.0, 3.0, 4.0, 7.0, 8.0, 7.0, 4.0, 9.0, 3.0, 8.0, 5.0, 1.0, 1.0, 2.0, 7.0, 1.0, 6.0, 1.0, 5.0, 9.0, 7.0, 8.0, 5.0, 6.0, 5.0, 4.0, 8.0, 1.0, 6.0, 1.0, 3.0, 5.0, 8.0, 1.0, 7.0, 1.0, 1.0, 7.0, 9.0, 6.0, 8.0, 4.0, 1.0, 4.0, 4.0, 1.0, 10.0, 5.0, 8.0, 3.0, 6.0, 2.0, 6.0, 6.0, 5.0, 3.0, 3.0, 5.0, 8.0, 7.0, 10.0, 3.0, 5.0, 3.0, 3.0, 1.0, 8.0, 6.0, 3.0, 8.0, 3.0, 10.0, 6.0, 2.0, 4.0, 3.0, 2.0, 2.0, 4.0, 10.0, 3.0, 10.0, 5.0, 7.0, 3.0, 8.0, 10.0, 8.0, 9.0, 1.0, 1.0, 8.0, 2.0, 6.0, 2.0, 8.0, 5.0, 8.0, 4.0, 6.0, 6.0, 8.0, 6.0, 5.0, 9.0, 6.0, 1.0, 6.0, 6.0, 10.0, 3.0, 10.0, 5.0, 9.0, 9.0, 9.0, 6.0, 7.0, 3.0, 5.0, 1.0, 1.0, 5.0, 1.0, 6.0, 9.0, 6.0, 2.0, 1.0, 2.0, 6.0, 2.0, 3.0, 1.0, 6.0, 2.0, 7.0, 9.0, 4.0, 7.0, 3.0, 6.0, 10.0, 3.0, 1.0, 2.0, 5.0, 8.0, 4.0, 10.0, 7.0, 6.0, 10.0, 3.0, 9.0, 3.0, 8.0, 1.0, 5.0, 4.0, 5.0, 7.0, 2.0, 9.0, 4.0, 4.0, 3.0, 6.0, 8.0, 4.0, 9.0, 5.0, 8.0, 6.0, 6.0, 4.0, 6.0, 5.0, 9.0, 6.0, 3.0, 9.0, 6.0, 7.0, 6.0, 8.0, 3.0, 4.0, 5.0, 10.0, 7.0, 7.0, 6.0, 10.0, 3.0, 9.0, 1.0, 8.0, 5.0, 10.0]
global b_x = 5
global d_y = [8.0, 10.0, 5.0, 7.0, 9.0, 5.0, 1.0, 8.0, 5.0, 2.0, 5.0, 10.0, 7.0, 8.0, 8.0, 7.0, 5.0, 3.0, 3.0, 4.0, 9.0, 3.0, 4.0, 9.0, 2.0, 7.0, 8.0, 6.0, 3.0, 5.0, 6.0, 10.0, 3.0, 3.0, 10.0, 9.0, 3.0, 8.0, 1.0, 9.0, 6.0, 10.0, 7.0, 7.0, 10.0, 8.0, 4.0, 4.0, 9.0, 5.0, 10.0, 7.0, 3.0, 5.0, 9.0, 4.0, 6.0, 9.0, 4.0, 1.0, 9.0, 2.0, 1.0, 10.0, 4.0, 1.0, 9.0, 5.0, 7.0, 5.0, 7.0, 4.0, 8.0, 4.0, 10.0, 9.0, 9.0, 10.0, 3.0, 5.0, 3.0, 10.0, 8.0, 4.0, 5.0, 5.0, 9.0, 5.0, 4.0, 10.0, 8.0, 9.0, 3.0, 5.0, 6.0, 1.0, 1.0, 8.0, 2.0, 3.0, 5.0, 8.0, 4.0, 3.0, 5.0, 6.0, 9.0, 3.0, 8.0, 3.0, 1.0, 2.0, 5.0, 1.0, 5.0, 6.0, 5.0, 7.0, 2.0, 4.0, 8.0, 3.0, 4.0, 8.0, 8.0, 3.0, 6.0, 8.0, 4.0, 7.0, 2.0, 1.0, 3.0, 5.0, 2.0, 7.0, 7.0, 8.0, 5.0, 4.0, 1.0, 6.0, 2.0, 2.0, 8.0, 10.0, 8.0, 2.0, 6.0, 8.0, 3.0, 1.0, 10.0, 3.0, 6.0, 2.0, 2.0, 5.0, 3.0, 4.0, 9.0, 1.0, 4.0, 2.0, 6.0, 7.0, 4.0, 8.0, 1.0, 10.0, 8.0, 6.0, 9.0, 4.0, 8.0, 9.0, 8.0, 1.0, 10.0, 5.0, 3.0, 8.0, 1.0, 8.0, 2.0, 9.0, 7.0, 9.0, 2.0, 10.0, 2.0, 2.0, 8.0, 10.0, 1.0, 1.0, 4.0, 10.0, 8.0, 6.0, 7.0, 7.0, 3.0, 2.0, 9.0, 3.0, 3.0, 5.0, 8.0, 6.0, 7.0, 5.0, 7.0, 3.0, 7.0, 10.0, 10.0, 9.0, 1.0, 4.0, 2.0, 4.0, 1.0, 5.0, 6.0, 7.0, 4.0, 1.0, 7.0, 3.0, 7.0, 10.0, 8.0, 9.0, 3.0, 6.0, 3.0, 6.0, 1.0, 1.0, 2.0, 8.0, 1.0, 4.0, 9.0, 8.0, 5.0, 8.0, 2.0, 2.0, 6.0, 6.0, 4.0, 3.0, 9.0, 3.0, 10.0, 2.0, 10.0, 5.0, 3.0, 4.0, 7.0, 9.0, 9.0, 8.0, 5.0, 1.0, 4.0, 10.0, 7.0, 8.0, 9.0, 7.0, 10.0, 6.0, 8.0, 4.0, 7.0, 9.0, 10.0, 5.0, 6.0, 8.0, 8.0, 9.0, 10.0, 5.0, 9.0, 5.0, 9.0, 5.0, 2.0, 3.0, 6.0, 4.0, 6.0, 10.0, 10.0, 5.0, 9.0, 2.0, 2.0, 3.0, 9.0, 9.0, 1.0, 6.0, 7.0, 7.0, 3.0, 1.0, 1.0, 3.0, 9.0, 8.0, 9.0, 2.0, 4.0, 9.0, 6.0, 1.0, 6.0, 2.0, 2.0, 1.0, 10.0, 3.0, 4.0, 1.0, 1.0, 4.0, 10.0, 3.0, 3.0, 6.0, 4.0, 8.0, 6.0, 5.0, 10.0]
global b_y = 10
global p = [0.752, 0.003, 0.19, 0.748, 0.534, 0.377, 0.18, 0.024, 0.07, 0.675, 0.602, 0.515, 0.372, 0.292, 0.817, 0.94, 0.98, 0.79, 0.327, 0.141, 0.001, 0.478, 0.382, 0.167, 0.358, 0.779, 0.943, 0.008, 0.628, 0.609, 0.594, 0.576, 0.745, 0.792, 0.707, 0.547, 0.324, 0.544, 0.91, 0.051, 0.817, 0.819, 0.761, 0.394, 0.691, 0.347, 0.502, 0.816, 0.357, 0.689, 0.191, 0.406, 0.204, 0.753, 0.136, 0.998, 0.402, 0.015, 0.33, 0.915, 0.975, 0.84, 0.515, 0.593, 0.806, 0.121, 0.163, 0.028, 0.4, 0.772, 0.381, 0.404, 0.806, 0.975, 0.52, 0.202, 0.637, 0.475, 0.522, 0.425, 0.19, 0.102, 0.68, 0.249, 0.751, 0.142, 0.428, 0.86, 0.119, 0.447, 0.665, 0.291, 0.377, 0.597, 0.874, 0.326, 0.752, 0.604, 0.248, 0.619, 0.291, 0.915, 0.024, 0.035, 0.116, 0.555, 0.097, 0.375, 0.763, 0.136, 0.57, 0.783, 0.463, 0.713, 0.899, 0.91, 0.55, 0.703, 0.017, 0.873, 0.467, 0.863, 0.978, 0.875, 0.279, 0.181, 0.54, 0.028, 0.211, 0.295, 0.49, 0.49, 0.724, 0.89, 0.363, 0.553, 0.126, 0.974, 0.769, 0.662, 0.32, 0.377, 0.155, 0.235, 0.704, 0.907, 0.855, 0.611, 0.437, 0.78, 0.715, 0.884, 0.377, 0.258, 0.516, 0.623, 0.917, 0.712, 0.314, 0.874, 0.322, 0.459, 0.981, 0.002, 0.238, 0.074, 0.989, 0.549, 0.205, 0.231, 0.117, 0.369, 0.898, 0.384, 0.732, 0.809, 0.548, 0.162, 0.601, 0.91, 0.447, 0.099, 0.162, 0.285, 0.917, 0.73, 0.087, 0.458, 0.2, 0.419, 0.461, 0.515, 0.996, 0.782, 0.976, 0.029, 0.562, 0.437, 0.988, 0.653, 0.986, 0.902, 0.445, 0.231, 0.394, 0.866, 0.997, 0.88, 0.713, 0.055, 0.02, 0.061, 0.342, 0.384, 0.357, 0.784, 0.795, 0.077, 0.199, 0.67, 0.146, 0.563, 0.244, 0.549, 0.376, 0.37, 0.106, 0.951, 0.445, 0.533, 0.095, 0.905, 0.778, 0.803, 0.031, 0.679, 0.39, 0.207, 0.004, 0.47, 0.453, 0.402, 0.707, 0.723, 0.336, 0.724, 0.768, 0.589, 0.237, 0.222, 0.238, 0.865, 0.308, 0.918, 0.143, 0.457, 0.888, 0.047, 0.615, 0.501, 0.768, 0.841, 0.896, 0.034, 0.2, 0.378, 0.66, 0.917, 0.234, 0.258, 0.616, 0.292, 0.024, 0.279, 0.696, 0.111, 0.126, 0.205, 0.65, 0.413, 0.136, 0.509, 0.739, 0.987, 0.446, 0.738, 0.562, 0.551, 0.894, 0.709, 0.349, 0.383, 0.895, 0.895, 0.009, 0.77, 0.737, 0.176, 0.082, 0.994, 0.057, 0.289, 0.679, 0.292, 0.465, 0.496, 0.436, 0.712, 0.001, 0.284, 0.801, 0.776, 0.978, 0.871, 0.846, 0.264, 0.343, 0.398, 0.542, 0.307, 0.416, 0.717, 0.241, 0.683, 0.292, 0.679, 0.846, 0.801, 0.743, 0.081, 0.374, 0.81, 0.184, 0.382, 0.068, 0.708, 0.845, 0.36, 0.152, 0.319, 0.35]
global q = [0.91, 0.09, 0.35, 0.881, 0.538, 0.763, 0.761, 0.039, 0.202, 0.762, 0.703, 0.981, 0.57, 0.473, 0.847, 0.998, 0.983, 0.943, 0.893, 0.718, 0.715, 0.964, 0.961, 0.581, 0.844, 0.851, 0.96, 0.144, 0.887, 0.753, 0.621, 0.895, 0.849, 0.996, 0.776, 0.657, 0.697, 0.815, 0.99, 0.352, 0.99, 0.978, 0.854, 0.864, 0.873, 0.52, 0.625, 0.822, 0.927, 0.834, 0.638, 0.83, 0.244, 0.788, 0.606, 0.999, 0.731, 0.303, 0.583, 0.971, 0.98, 0.968, 0.864, 0.936, 0.829, 0.287, 0.348, 0.303, 0.872, 0.974, 0.523, 0.795, 0.998, 0.982, 0.535, 0.314, 0.919, 0.536, 0.625, 0.985, 0.44, 0.648, 0.826, 0.275, 0.82, 0.512, 0.813, 0.867, 0.182, 0.746, 0.92, 0.54, 0.434, 0.956, 0.929, 0.615, 0.95, 0.609, 0.424, 0.858, 0.84, 0.943, 0.095, 0.678, 0.795, 0.706, 0.181, 0.806, 0.968, 0.682, 0.73, 0.831, 0.911, 0.876, 0.996, 0.989, 0.861, 0.939, 0.086, 0.906, 0.778, 0.93, 0.994, 0.879, 0.765, 0.264, 0.819, 0.321, 0.818, 0.48, 0.537, 0.78, 0.998, 0.893, 0.976, 0.967, 0.488, 0.994, 0.802, 0.849, 0.844, 0.578, 0.962, 0.516, 0.916, 0.983, 0.873, 0.69, 0.61, 0.861, 0.978, 0.904, 0.721, 0.763, 0.709, 0.631, 0.958, 0.971, 0.568, 0.974, 0.843, 0.598, 0.997, 0.946, 0.687, 0.851, 0.997, 0.58, 0.374, 0.95, 0.615, 0.41, 0.982, 0.798, 0.787, 0.888, 0.548, 0.746, 0.667, 0.988, 0.487, 0.282, 0.797, 0.428, 0.999, 0.868, 0.7, 0.896, 0.201, 0.454, 0.647, 0.638, 0.998, 0.917, 0.989, 0.337, 0.875, 0.602, 0.999, 0.828, 0.987, 0.957, 0.976, 0.842, 0.711, 0.943, 0.999, 0.991, 0.899, 0.603, 0.262, 0.756, 0.849, 0.644, 0.78, 0.895, 0.831, 0.387, 0.564, 0.69, 0.61, 0.897, 0.523, 0.68, 0.406, 0.668, 0.508, 0.967, 0.814, 0.543, 0.969, 0.973, 0.952, 0.811, 0.79, 0.92, 0.749, 0.513, 0.684, 0.818, 0.516, 0.52, 0.823, 0.759, 0.979, 0.993, 0.972, 0.811, 0.529, 0.348, 0.245, 0.981, 0.76, 0.961, 0.528, 0.878, 0.928, 0.927, 0.671, 0.62, 0.82, 0.931, 0.913, 0.292, 0.385, 0.53, 0.765, 0.92, 0.358, 0.876, 0.68, 0.779, 0.741, 0.771, 0.843, 0.349, 0.137, 0.758, 0.848, 0.909, 0.702, 0.909, 0.764, 0.988, 0.989, 0.866, 0.596, 0.953, 0.954, 0.821, 0.586, 0.806, 0.985, 0.918, 0.168, 0.818, 0.892, 0.694, 0.542, 0.998, 0.951, 0.885, 0.737, 0.938, 0.67, 0.859, 0.55, 0.913, 0.235, 0.972, 0.894, 0.854, 0.978, 0.876, 0.872, 0.663, 0.71, 0.589, 0.985, 0.879, 0.776, 0.777, 0.362, 0.742, 0.959, 0.766, 0.919, 0.916, 0.76, 0.37, 0.468, 0.865, 0.971, 0.806, 0.62, 0.936, 0.857, 0.881, 0.407, 0.608, 0.617]
global origin = 1
global destination = 60