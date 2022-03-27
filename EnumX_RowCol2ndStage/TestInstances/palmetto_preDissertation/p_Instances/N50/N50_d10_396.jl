global arcs = [1 36; 1 41; 1 43; 1 49; 2 16; 2 18; 2 19; 2 36; 2 45; 3 11; 3 34; 3 36; 3 47; 4 30; 4 44; 4 48; 5 4; 5 21; 5 39; 5 40; 5 44; 6 2; 6 3; 6 12; 6 18; 6 34; 6 39; 6 41; 7 13; 7 15; 7 16; 7 27; 7 29; 7 33; 7 36; 7 46; 7 47; 8 4; 8 10; 8 15; 8 18; 8 21; 8 32; 8 41; 8 48; 9 8; 9 10; 9 26; 9 30; 9 32; 10 7; 10 15; 10 21; 10 26; 10 38; 11 3; 11 10; 11 17; 11 24; 11 30; 11 38; 11 44; 11 45; 11 48; 12 26; 12 34; 12 43; 12 50; 13 4; 13 7; 13 11; 13 24; 13 43; 14 9; 14 12; 14 22; 14 30; 15 8; 15 21; 15 31; 15 34; 15 35; 15 39; 16 19; 16 34; 16 38; 16 46; 16 47; 17 9; 17 14; 17 26; 17 40; 18 14; 18 35; 18 44; 18 45; 19 20; 19 27; 19 39; 20 6; 20 7; 20 12; 20 16; 20 19; 20 22; 20 31; 20 38; 20 50; 21 7; 21 13; 21 25; 21 50; 22 4; 22 13; 22 16; 22 30; 23 5; 24 7; 24 18; 24 21; 24 29; 24 32; 25 7; 25 10; 25 16; 25 17; 25 21; 25 28; 25 37; 25 40; 25 47; 26 24; 26 32; 26 34; 26 39; 27 12; 27 23; 27 31; 27 40; 27 43; 27 45; 27 49; 28 7; 28 9; 28 39; 28 40; 28 46; 29 3; 29 8; 29 25; 30 7; 30 8; 30 12; 30 27; 30 34; 31 2; 31 16; 31 32; 32 17; 32 19; 32 34; 33 12; 33 40; 33 41; 33 42; 34 12; 34 15; 34 43; 35 7; 35 16; 35 19; 35 32; 35 40; 36 3; 36 7; 36 12; 36 27; 37 16; 37 36; 37 39; 38 7; 38 9; 38 12; 38 20; 38 43; 39 6; 39 13; 39 19; 39 28; 40 15; 40 19; 40 26; 40 32; 40 33; 41 9; 42 8; 42 20; 42 33; 42 36; 42 39; 42 43; 42 44; 42 50; 43 8; 43 17; 44 6; 44 45; 44 49; 45 7; 45 10; 45 11; 45 30; 45 33; 45 36; 45 43; 46 6; 46 15; 46 24; 46 28; 46 38; 46 39; 47 3; 47 30; 48 8; 48 30; 48 35; 48 36; 48 39; 48 45; 49 13]
global d_x = [10.0, 5.0, 6.0, 5.0, 3.0, 8.0, 10.0, 5.0, 1.0, 1.0, 5.0, 7.0, 3.0, 1.0, 5.0, 1.0, 10.0, 8.0, 4.0, 5.0, 6.0, 6.0, 1.0, 5.0, 6.0, 5.0, 1.0, 9.0, 9.0, 3.0, 5.0, 9.0, 3.0, 9.0, 4.0, 9.0, 3.0, 10.0, 9.0, 8.0, 3.0, 10.0, 3.0, 3.0, 7.0, 4.0, 6.0, 5.0, 6.0, 10.0, 4.0, 2.0, 6.0, 6.0, 4.0, 8.0, 7.0, 2.0, 8.0, 7.0, 2.0, 5.0, 10.0, 10.0, 4.0, 1.0, 10.0, 8.0, 4.0, 10.0, 2.0, 5.0, 6.0, 3.0, 3.0, 9.0, 6.0, 9.0, 5.0, 9.0, 9.0, 4.0, 9.0, 7.0, 7.0, 5.0, 4.0, 9.0, 9.0, 5.0, 1.0, 6.0, 8.0, 9.0, 2.0, 10.0, 9.0, 6.0, 1.0, 9.0, 2.0, 1.0, 3.0, 1.0, 6.0, 5.0, 6.0, 3.0, 6.0, 8.0, 2.0, 8.0, 9.0, 9.0, 8.0, 9.0, 8.0, 6.0, 1.0, 5.0, 5.0, 4.0, 7.0, 10.0, 6.0, 8.0, 6.0, 8.0, 8.0, 8.0, 6.0, 4.0, 8.0, 10.0, 1.0, 10.0, 3.0, 7.0, 4.0, 3.0, 2.0, 4.0, 8.0, 6.0, 8.0, 4.0, 5.0, 5.0, 6.0, 1.0, 10.0, 6.0, 5.0, 2.0, 3.0, 2.0, 5.0, 4.0, 3.0, 3.0, 9.0, 2.0, 2.0, 2.0, 10.0, 4.0, 5.0, 2.0, 9.0, 8.0, 8.0, 1.0, 8.0, 10.0, 3.0, 1.0, 8.0, 4.0, 4.0, 7.0, 9.0, 6.0, 2.0, 3.0, 5.0, 4.0, 2.0, 9.0, 1.0, 1.0, 9.0, 5.0, 4.0, 8.0, 6.0, 7.0, 4.0, 6.0, 1.0, 1.0, 4.0, 10.0, 8.0, 7.0, 7.0, 10.0, 3.0, 6.0, 5.0, 5.0, 1.0, 5.0, 3.0, 10.0, 1.0, 1.0, 3.0, 7.0, 4.0, 4.0, 1.0, 9.0, 7.0, 8.0, 9.0, 10.0, 1.0, 10.0, 5.0, 5.0]
global b_x = 5
global d_y = [7.0, 9.0, 9.0, 1.0, 2.0, 4.0, 5.0, 9.0, 4.0, 3.0, 3.0, 7.0, 3.0, 7.0, 8.0, 3.0, 4.0, 7.0, 4.0, 7.0, 10.0, 1.0, 2.0, 8.0, 10.0, 10.0, 2.0, 7.0, 2.0, 3.0, 4.0, 9.0, 6.0, 5.0, 3.0, 8.0, 4.0, 3.0, 6.0, 3.0, 2.0, 2.0, 7.0, 4.0, 4.0, 7.0, 6.0, 9.0, 8.0, 9.0, 8.0, 8.0, 7.0, 10.0, 2.0, 10.0, 7.0, 4.0, 6.0, 5.0, 2.0, 4.0, 10.0, 9.0, 10.0, 6.0, 1.0, 9.0, 4.0, 9.0, 5.0, 7.0, 4.0, 1.0, 7.0, 7.0, 5.0, 6.0, 4.0, 10.0, 3.0, 8.0, 6.0, 7.0, 6.0, 2.0, 6.0, 3.0, 6.0, 6.0, 2.0, 2.0, 2.0, 8.0, 5.0, 5.0, 4.0, 4.0, 10.0, 10.0, 7.0, 6.0, 5.0, 10.0, 5.0, 4.0, 8.0, 4.0, 1.0, 10.0, 2.0, 10.0, 2.0, 1.0, 3.0, 2.0, 9.0, 4.0, 1.0, 4.0, 5.0, 3.0, 6.0, 4.0, 8.0, 4.0, 4.0, 5.0, 8.0, 3.0, 5.0, 5.0, 2.0, 2.0, 10.0, 3.0, 6.0, 10.0, 2.0, 6.0, 1.0, 3.0, 8.0, 8.0, 4.0, 4.0, 1.0, 7.0, 1.0, 10.0, 10.0, 9.0, 6.0, 2.0, 3.0, 10.0, 7.0, 2.0, 6.0, 3.0, 10.0, 4.0, 2.0, 5.0, 5.0, 7.0, 10.0, 4.0, 6.0, 9.0, 9.0, 4.0, 6.0, 5.0, 9.0, 3.0, 9.0, 5.0, 4.0, 5.0, 8.0, 2.0, 2.0, 6.0, 2.0, 5.0, 9.0, 8.0, 5.0, 5.0, 7.0, 2.0, 10.0, 5.0, 10.0, 1.0, 1.0, 8.0, 3.0, 9.0, 1.0, 3.0, 5.0, 8.0, 2.0, 6.0, 3.0, 1.0, 2.0, 3.0, 1.0, 9.0, 2.0, 1.0, 9.0, 4.0, 1.0, 2.0, 9.0, 5.0, 6.0, 9.0, 2.0, 2.0, 7.0, 7.0, 2.0, 6.0, 9.0, 6.0]
global b_y = 10
global p = [0.833, 0.612, 0.138, 0.505, 0.972, 0.23, 0.314, 0.153, 0.755, 0.059, 0.788, 0.161, 0.059, 0.431, 0.886, 0.317, 0.244, 0.327, 0.042, 0.257, 0.518, 0.007, 0.96, 0.756, 0.238, 0.043, 0.651, 0.042, 0.243, 0.524, 0.802, 0.511, 0.109, 0.765, 0.528, 0.185, 0.998, 0.215, 0.676, 0.993, 0.336, 0.635, 0.666, 0.159, 0.438, 0.699, 0.397, 0.588, 0.63, 0.125, 0.744, 0.884, 0.478, 0.87, 0.839, 0.518, 0.04, 0.192, 0.179, 0.594, 0.228, 0.451, 0.384, 0.968, 0.949, 0.002, 0.826, 0.226, 0.361, 0.036, 0.729, 0.307, 0.81, 0.428, 0.957, 0.708, 0.723, 0.59, 0.84, 0.85, 0.185, 0.286, 0.757, 0.572, 0.776, 0.905, 0.999, 0.212, 0.133, 0.338, 0.173, 0.118, 0.21, 0.776, 0.227, 0.635, 0.737, 0.891, 0.993, 0.906, 0.68, 0.919, 0.779, 0.92, 0.272, 0.218, 0.17, 0.323, 0.646, 0.329, 0.234, 0.595, 0.653, 0.243, 0.829, 0.663, 0.522, 0.879, 0.584, 0.337, 0.864, 0.018, 0.719, 0.787, 0.075, 0.455, 0.085, 0.501, 0.703, 0.792, 0.407, 0.53, 0.421, 0.79, 0.939, 0.624, 0.58, 0.97, 0.647, 0.415, 0.012, 0.43, 0.608, 0.747, 0.166, 0.991, 0.656, 0.661, 0.188, 0.309, 0.303, 0.998, 0.478, 0.172, 0.287, 0.909, 0.872, 0.736, 0.936, 0.803, 0.144, 0.902, 0.213, 0.433, 0.356, 0.851, 0.664, 0.71, 0.328, 0.635, 0.018, 0.976, 0.523, 0.641, 0.114, 0.201, 0.982, 0.222, 0.818, 0.306, 0.14, 0.933, 0.203, 0.333, 0.624, 0.156, 0.785, 0.725, 0.698, 0.585, 0.615, 0.277, 0.376, 0.232, 0.032, 0.559, 0.39, 0.988, 0.605, 0.192, 0.574, 0.411, 0.832, 0.063, 0.384, 0.281, 0.725, 0.058, 0.366, 0.704, 0.862, 0.08, 0.126, 0.322, 0.857, 0.131, 0.733, 0.469, 0.365, 0.283, 0.021, 0.627, 0.431, 0.091, 0.221, 0.937, 0.303, 0.429, 0.861, 0.209]
global q = [0.909, 0.672, 0.551, 0.982, 0.974, 0.587, 0.447, 0.717, 0.964, 0.842, 0.928, 0.659, 0.481, 0.738, 0.968, 0.642, 0.332, 0.718, 0.741, 0.634, 0.613, 0.307, 0.974, 0.795, 0.512, 0.494, 0.682, 0.516, 0.561, 0.999, 0.851, 0.589, 0.806, 0.821, 0.681, 0.579, 0.998, 0.719, 0.76, 0.999, 0.966, 0.717, 0.887, 0.919, 0.927, 0.833, 0.875, 0.964, 0.931, 0.378, 0.845, 0.973, 0.619, 0.993, 0.885, 0.976, 0.332, 0.557, 0.995, 0.806, 0.236, 0.655, 0.707, 0.978, 0.951, 0.335, 0.956, 0.387, 0.445, 0.058, 0.859, 0.877, 0.938, 0.834, 0.973, 0.931, 0.831, 0.612, 0.98, 0.966, 0.955, 0.741, 0.927, 0.77, 0.777, 0.973, 0.999, 0.514, 0.946, 0.568, 0.932, 0.619, 0.603, 0.779, 0.693, 0.646, 0.949, 0.966, 0.996, 0.964, 0.992, 0.948, 0.799, 0.989, 0.719, 0.223, 0.402, 0.434, 0.826, 0.696, 0.719, 0.83, 0.893, 0.745, 0.938, 0.717, 0.565, 0.982, 0.631, 0.651, 0.969, 0.848, 0.906, 0.955, 0.136, 0.804, 0.567, 0.903, 0.885, 0.922, 0.785, 0.578, 0.927, 0.961, 0.966, 0.81, 0.863, 0.975, 0.98, 0.955, 0.967, 0.79, 0.695, 0.754, 0.747, 0.999, 0.674, 0.986, 0.726, 0.435, 0.702, 0.999, 0.971, 0.885, 0.604, 0.93, 0.879, 0.894, 0.993, 0.861, 0.444, 0.942, 0.925, 0.982, 0.747, 0.927, 0.683, 0.852, 0.36, 0.905, 0.076, 0.976, 0.816, 0.808, 0.419, 0.881, 0.996, 0.97, 0.951, 0.565, 0.423, 0.978, 0.62, 0.941, 0.915, 0.946, 0.893, 0.947, 0.743, 0.894, 0.784, 0.347, 0.658, 0.289, 0.637, 0.712, 0.655, 0.992, 0.778, 0.719, 0.744, 0.66, 0.934, 0.247, 0.802, 0.29, 0.952, 0.074, 0.439, 0.759, 0.863, 0.831, 0.459, 0.71, 0.998, 0.326, 0.98, 0.758, 0.863, 0.731, 0.713, 0.866, 0.521, 0.977, 0.25, 0.969, 0.97, 0.704, 0.907, 0.885]
global origin = 1
global destination = 50