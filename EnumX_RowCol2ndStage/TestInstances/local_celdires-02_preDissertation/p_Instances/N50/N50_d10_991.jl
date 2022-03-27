global arcs = [1 4; 1 8; 1 10; 1 27; 1 38; 1 39; 1 42; 1 48; 2 4; 2 14; 2 17; 2 35; 2 37; 2 38; 2 42; 3 4; 3 5; 3 24; 3 27; 3 29; 3 41; 4 8; 4 11; 4 31; 5 7; 5 15; 5 16; 5 29; 6 14; 6 18; 6 20; 6 23; 6 26; 6 27; 6 38; 7 2; 7 13; 7 18; 7 20; 8 6; 8 20; 8 35; 9 7; 9 8; 9 12; 9 14; 9 33; 9 44; 9 50; 10 21; 10 25; 10 29; 11 6; 11 14; 11 15; 11 29; 11 31; 11 45; 12 6; 12 21; 12 41; 13 7; 13 29; 13 31; 13 37; 13 41; 13 44; 13 47; 14 2; 14 8; 14 19; 14 20; 15 18; 15 20; 15 31; 15 41; 15 42; 16 6; 16 17; 16 18; 16 25; 16 36; 17 2; 17 20; 17 21; 17 24; 17 38; 18 12; 18 22; 18 33; 18 35; 19 10; 19 23; 19 38; 20 3; 20 29; 20 35; 20 40; 21 8; 22 3; 22 10; 22 18; 22 41; 23 3; 23 6; 23 7; 23 26; 23 44; 23 47; 24 10; 24 15; 24 22; 24 23; 24 25; 24 44; 24 48; 25 17; 25 27; 25 40; 25 46; 25 47; 25 49; 26 17; 27 30; 28 5; 28 21; 28 46; 28 50; 29 4; 29 10; 29 36; 29 40; 29 48; 30 4; 30 5; 30 13; 30 50; 31 3; 31 18; 31 23; 31 38; 32 14; 32 19; 32 23; 32 24; 32 37; 32 41; 33 5; 33 9; 33 14; 33 18; 33 29; 33 31; 33 38; 34 37; 34 38; 34 41; 34 48; 35 2; 35 6; 35 18; 35 48; 36 5; 36 22; 36 28; 36 34; 36 38; 36 44; 37 3; 37 4; 37 12; 37 16; 37 36; 37 40; 37 48; 38 4; 38 12; 38 13; 38 14; 38 23; 38 36; 39 2; 39 3; 39 22; 39 35; 39 43; 39 49; 40 2; 40 7; 40 8; 40 9; 40 18; 40 29; 40 32; 40 39; 40 50; 41 2; 41 10; 41 18; 41 32; 41 39; 41 42; 42 33; 42 38; 42 41; 42 47; 43 19; 43 25; 44 8; 44 11; 44 29; 44 33; 45 11; 45 20; 45 23; 45 26; 45 48; 45 50; 46 21; 46 32; 46 39; 47 3; 47 9; 47 32; 48 18; 48 26; 49 15; 49 17; 49 19; 49 23; 49 43; 49 47]
global d_x = [1.0, 10.0, 6.0, 5.0, 2.0, 9.0, 10.0, 3.0, 6.0, 10.0, 8.0, 5.0, 6.0, 6.0, 7.0, 7.0, 7.0, 7.0, 7.0, 3.0, 10.0, 8.0, 5.0, 6.0, 8.0, 2.0, 6.0, 4.0, 5.0, 4.0, 9.0, 7.0, 2.0, 4.0, 5.0, 4.0, 2.0, 1.0, 9.0, 1.0, 5.0, 6.0, 5.0, 7.0, 3.0, 9.0, 1.0, 8.0, 9.0, 1.0, 5.0, 9.0, 8.0, 10.0, 3.0, 9.0, 10.0, 6.0, 7.0, 4.0, 8.0, 7.0, 2.0, 3.0, 1.0, 10.0, 7.0, 8.0, 4.0, 3.0, 7.0, 3.0, 2.0, 7.0, 5.0, 10.0, 2.0, 6.0, 2.0, 9.0, 4.0, 10.0, 1.0, 3.0, 5.0, 4.0, 5.0, 2.0, 6.0, 1.0, 1.0, 3.0, 10.0, 2.0, 9.0, 3.0, 4.0, 7.0, 6.0, 5.0, 7.0, 6.0, 10.0, 6.0, 5.0, 8.0, 7.0, 4.0, 2.0, 8.0, 6.0, 3.0, 7.0, 1.0, 2.0, 4.0, 9.0, 3.0, 7.0, 10.0, 3.0, 1.0, 1.0, 5.0, 8.0, 2.0, 3.0, 6.0, 3.0, 9.0, 1.0, 1.0, 9.0, 5.0, 1.0, 6.0, 8.0, 8.0, 7.0, 2.0, 2.0, 7.0, 4.0, 8.0, 3.0, 9.0, 7.0, 3.0, 1.0, 6.0, 4.0, 2.0, 3.0, 10.0, 3.0, 7.0, 5.0, 3.0, 7.0, 5.0, 6.0, 1.0, 3.0, 9.0, 6.0, 10.0, 10.0, 4.0, 2.0, 3.0, 6.0, 7.0, 9.0, 1.0, 7.0, 8.0, 4.0, 3.0, 2.0, 10.0, 4.0, 4.0, 2.0, 9.0, 8.0, 7.0, 8.0, 8.0, 8.0, 2.0, 8.0, 8.0, 4.0, 10.0, 5.0, 1.0, 1.0, 9.0, 9.0, 4.0, 3.0, 8.0, 10.0, 8.0, 8.0, 6.0, 2.0, 4.0, 7.0, 8.0, 6.0, 7.0, 7.0, 6.0, 7.0, 9.0, 4.0, 5.0, 4.0, 9.0, 8.0, 10.0, 1.0, 3.0, 8.0, 3.0, 5.0, 9.0, 1.0, 7.0, 3.0, 7.0]
global b_x = 5
global d_y = [4.0, 9.0, 2.0, 10.0, 9.0, 8.0, 1.0, 4.0, 8.0, 9.0, 8.0, 6.0, 1.0, 3.0, 6.0, 4.0, 4.0, 3.0, 9.0, 4.0, 4.0, 1.0, 3.0, 6.0, 7.0, 4.0, 5.0, 6.0, 4.0, 6.0, 4.0, 8.0, 6.0, 7.0, 4.0, 8.0, 8.0, 3.0, 4.0, 3.0, 9.0, 9.0, 4.0, 9.0, 8.0, 10.0, 7.0, 10.0, 7.0, 4.0, 2.0, 8.0, 3.0, 7.0, 8.0, 5.0, 2.0, 6.0, 10.0, 5.0, 10.0, 8.0, 5.0, 6.0, 3.0, 3.0, 5.0, 1.0, 6.0, 1.0, 8.0, 8.0, 7.0, 2.0, 5.0, 8.0, 9.0, 6.0, 9.0, 1.0, 1.0, 1.0, 4.0, 2.0, 10.0, 7.0, 4.0, 2.0, 7.0, 3.0, 2.0, 2.0, 3.0, 10.0, 5.0, 1.0, 8.0, 3.0, 1.0, 7.0, 1.0, 8.0, 4.0, 6.0, 8.0, 4.0, 9.0, 8.0, 8.0, 5.0, 8.0, 8.0, 3.0, 5.0, 5.0, 2.0, 6.0, 8.0, 3.0, 7.0, 8.0, 2.0, 9.0, 5.0, 2.0, 5.0, 7.0, 10.0, 6.0, 5.0, 8.0, 10.0, 1.0, 4.0, 2.0, 7.0, 8.0, 10.0, 7.0, 10.0, 3.0, 6.0, 9.0, 8.0, 10.0, 5.0, 2.0, 9.0, 10.0, 9.0, 1.0, 7.0, 6.0, 6.0, 8.0, 5.0, 7.0, 10.0, 2.0, 1.0, 10.0, 8.0, 4.0, 1.0, 4.0, 2.0, 9.0, 2.0, 7.0, 5.0, 6.0, 7.0, 5.0, 2.0, 1.0, 10.0, 5.0, 2.0, 4.0, 4.0, 8.0, 1.0, 2.0, 3.0, 7.0, 9.0, 2.0, 8.0, 5.0, 4.0, 1.0, 6.0, 7.0, 4.0, 5.0, 6.0, 1.0, 7.0, 3.0, 5.0, 3.0, 8.0, 3.0, 7.0, 8.0, 6.0, 2.0, 5.0, 8.0, 1.0, 3.0, 9.0, 6.0, 7.0, 2.0, 7.0, 2.0, 8.0, 9.0, 4.0, 8.0, 3.0, 8.0, 2.0, 3.0, 5.0, 5.0, 5.0, 9.0, 7.0, 5.0, 8.0]
global b_y = 10
global p = [0.818, 0.651, 0.481, 0.518, 0.106, 0.694, 0.961, 0.942, 0.437, 0.776, 0.46, 0.871, 0.995, 0.141, 0.167, 0.107, 0.956, 0.196, 0.488, 0.605, 0.266, 0.849, 0.399, 0.811, 0.744, 0.433, 0.724, 0.549, 0.175, 0.194, 0.068, 0.769, 0.717, 0.153, 0.144, 0.126, 0.531, 0.1, 0.938, 0.281, 0.096, 0.465, 0.547, 0.327, 0.982, 0.936, 0.42, 0.156, 0.214, 0.397, 0.935, 0.408, 0.071, 0.325, 0.583, 0.186, 0.98, 0.482, 0.508, 0.227, 0.151, 0.996, 0.275, 0.016, 0.813, 0.696, 0.153, 0.158, 0.893, 0.671, 0.596, 0.613, 0.609, 0.202, 0.626, 0.92, 0.036, 0.553, 0.668, 0.266, 0.836, 0.452, 0.974, 0.031, 0.451, 0.779, 0.994, 0.037, 0.597, 0.356, 0.009, 0.108, 0.61, 0.441, 0.887, 0.085, 0.725, 0.656, 0.229, 0.811, 0.018, 0.768, 0.302, 0.749, 0.743, 0.307, 0.03, 0.585, 0.816, 0.425, 0.064, 0.433, 0.956, 0.57, 0.998, 0.736, 0.567, 0.973, 0.773, 0.258, 0.138, 0.514, 0.522, 0.098, 0.633, 0.05, 0.04, 0.315, 0.122, 0.024, 0.8, 0.521, 0.194, 0.466, 0.829, 0.179, 0.865, 0.029, 0.4, 0.569, 0.666, 0.637, 0.235, 0.918, 0.912, 0.85, 0.919, 0.684, 0.311, 0.685, 0.05, 0.196, 0.823, 0.854, 0.666, 0.714, 0.922, 0.586, 0.254, 0.558, 0.531, 0.752, 0.954, 0.934, 0.204, 0.627, 0.403, 0.878, 0.629, 0.441, 0.767, 0.969, 0.347, 0.787, 0.484, 0.476, 0.202, 0.264, 0.009, 0.953, 0.076, 0.687, 0.958, 0.152, 0.508, 0.255, 0.166, 0.136, 0.174, 0.009, 0.339, 0.049, 0.102, 0.629, 0.732, 0.755, 0.663, 0.029, 0.855, 0.575, 0.789, 0.599, 0.879, 0.716, 0.092, 0.691, 0.758, 0.847, 0.774, 0.464, 0.506, 0.66, 0.353, 0.204, 0.569, 0.764, 0.002, 0.785, 0.931, 0.794, 0.908, 0.891, 0.432, 0.617, 0.047, 0.583, 0.898, 0.955, 0.61, 0.265, 0.408, 0.012]
global q = [0.944, 0.724, 0.658, 0.647, 0.224, 0.695, 0.966, 0.944, 0.627, 0.915, 0.872, 0.993, 0.998, 0.844, 0.185, 0.722, 0.966, 0.97, 0.744, 0.661, 0.886, 0.956, 0.666, 0.878, 0.76, 0.685, 0.839, 0.59, 0.55, 0.636, 0.592, 0.945, 0.808, 0.69, 0.57, 0.739, 0.797, 0.133, 0.949, 0.942, 0.117, 0.888, 0.802, 0.391, 0.996, 0.944, 0.934, 0.582, 0.266, 0.607, 0.973, 0.565, 0.184, 0.454, 0.604, 0.635, 0.983, 0.514, 0.934, 0.923, 0.578, 0.997, 0.993, 0.437, 0.855, 0.734, 0.175, 0.2, 0.969, 0.718, 0.791, 0.868, 0.944, 0.511, 0.86, 0.941, 0.509, 0.902, 0.724, 0.935, 0.96, 0.547, 0.985, 0.031, 0.621, 0.939, 0.999, 0.632, 0.69, 0.553, 0.607, 0.611, 0.807, 0.569, 0.994, 0.121, 0.909, 0.767, 0.461, 0.814, 0.585, 0.795, 0.79, 0.872, 0.768, 0.595, 0.062, 0.632, 0.838, 0.748, 0.084, 0.559, 0.979, 0.978, 0.999, 0.981, 0.634, 0.992, 0.949, 0.829, 0.529, 0.52, 0.618, 0.425, 0.831, 0.053, 0.944, 0.915, 0.765, 0.232, 0.926, 0.836, 0.384, 0.523, 0.924, 0.179, 0.874, 0.917, 0.939, 0.941, 0.978, 0.765, 0.453, 0.984, 0.981, 0.932, 0.94, 0.937, 0.534, 0.858, 0.987, 0.531, 0.983, 0.998, 0.974, 0.932, 0.973, 0.813, 0.285, 0.815, 0.604, 0.91, 0.982, 0.977, 0.255, 0.987, 0.697, 0.994, 0.804, 0.477, 0.893, 0.977, 0.518, 0.93, 0.916, 0.953, 0.78, 0.96, 0.978, 0.972, 0.682, 0.77, 0.997, 0.459, 0.579, 0.339, 0.8, 0.762, 0.277, 0.778, 0.653, 0.425, 0.169, 0.718, 0.787, 0.823, 0.86, 0.06, 0.879, 0.876, 0.845, 0.699, 0.975, 0.85, 0.259, 0.911, 0.878, 0.912, 0.914, 0.712, 0.796, 0.737, 0.599, 0.962, 0.97, 0.863, 0.484, 0.93, 0.962, 0.996, 0.987, 0.968, 0.618, 0.757, 0.868, 0.804, 0.907, 0.975, 0.941, 0.66, 0.845, 0.286]
global origin = 1
global destination = 50