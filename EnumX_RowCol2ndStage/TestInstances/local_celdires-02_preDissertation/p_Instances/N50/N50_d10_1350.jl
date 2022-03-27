global arcs = [1 2; 1 7; 1 20; 1 33; 1 37; 2 6; 2 12; 2 18; 2 22; 2 36; 3 13; 3 22; 3 29; 3 33; 3 39; 3 44; 3 45; 4 28; 4 39; 5 2; 5 12; 5 42; 6 2; 6 30; 6 38; 6 42; 6 43; 7 11; 7 16; 7 19; 7 26; 7 33; 7 37; 8 18; 8 36; 8 45; 8 48; 9 18; 9 20; 9 25; 9 46; 9 48; 10 16; 10 34; 10 41; 10 49; 11 10; 11 24; 11 27; 11 28; 11 43; 12 5; 12 18; 12 26; 12 32; 12 34; 12 35; 13 2; 13 6; 13 8; 13 36; 13 47; 14 7; 14 12; 14 20; 14 28; 14 43; 14 48; 15 49; 16 7; 16 9; 16 19; 16 32; 16 38; 17 4; 17 7; 17 23; 17 25; 17 26; 17 38; 18 2; 18 15; 18 30; 18 48; 19 12; 19 13; 19 24; 19 43; 19 50; 20 21; 20 36; 20 42; 21 7; 21 14; 21 15; 21 32; 21 45; 21 48; 22 6; 22 29; 23 20; 23 30; 23 36; 23 39; 23 47; 24 2; 24 9; 24 13; 24 17; 24 21; 24 41; 25 39; 25 40; 25 42; 25 43; 26 11; 26 33; 26 45; 26 46; 27 28; 27 43; 28 22; 28 30; 28 34; 28 37; 28 39; 28 45; 29 8; 29 12; 29 18; 29 25; 29 33; 29 36; 30 3; 30 9; 30 13; 30 25; 30 28; 30 34; 30 37; 31 3; 31 12; 31 14; 31 15; 31 24; 31 25; 31 27; 31 49; 32 7; 32 8; 32 10; 32 16; 32 23; 33 6; 33 18; 33 24; 33 26; 33 35; 33 36; 34 6; 34 7; 34 13; 34 33; 35 13; 35 26; 35 27; 35 31; 35 36; 36 3; 36 5; 36 6; 36 9; 36 11; 36 24; 36 30; 36 40; 37 15; 37 34; 38 35; 38 50; 39 12; 39 42; 39 49; 40 7; 40 43; 41 7; 41 23; 41 25; 41 26; 41 37; 41 38; 41 40; 42 27; 42 28; 42 34; 42 41; 43 11; 43 14; 43 16; 43 17; 43 33; 43 36; 43 41; 44 22; 44 25; 44 49; 44 50; 45 18; 45 30; 45 44; 46 5; 46 6; 46 7; 46 30; 46 38; 46 41; 46 43; 47 4; 47 12; 47 16; 47 37; 47 38; 47 45; 48 10; 48 21; 48 33; 49 10; 49 12; 49 14; 49 22; 49 38; 49 41; 49 44]
global d_x = [7.0, 7.0, 10.0, 8.0, 6.0, 8.0, 4.0, 3.0, 8.0, 3.0, 3.0, 8.0, 10.0, 5.0, 2.0, 1.0, 3.0, 6.0, 9.0, 5.0, 9.0, 2.0, 6.0, 6.0, 8.0, 8.0, 6.0, 2.0, 1.0, 10.0, 5.0, 4.0, 4.0, 5.0, 10.0, 3.0, 1.0, 5.0, 9.0, 10.0, 6.0, 9.0, 6.0, 6.0, 4.0, 7.0, 7.0, 10.0, 7.0, 4.0, 4.0, 9.0, 3.0, 8.0, 4.0, 2.0, 5.0, 5.0, 5.0, 2.0, 9.0, 7.0, 2.0, 2.0, 3.0, 3.0, 6.0, 2.0, 10.0, 3.0, 7.0, 7.0, 4.0, 10.0, 6.0, 3.0, 8.0, 5.0, 2.0, 3.0, 3.0, 3.0, 7.0, 5.0, 9.0, 4.0, 6.0, 7.0, 2.0, 9.0, 10.0, 1.0, 2.0, 9.0, 9.0, 3.0, 6.0, 3.0, 10.0, 6.0, 1.0, 5.0, 7.0, 10.0, 5.0, 7.0, 1.0, 4.0, 9.0, 2.0, 4.0, 6.0, 2.0, 3.0, 2.0, 1.0, 4.0, 6.0, 6.0, 7.0, 8.0, 1.0, 7.0, 3.0, 4.0, 2.0, 2.0, 7.0, 9.0, 2.0, 7.0, 5.0, 7.0, 8.0, 8.0, 10.0, 1.0, 1.0, 5.0, 6.0, 6.0, 6.0, 2.0, 1.0, 7.0, 1.0, 3.0, 5.0, 9.0, 5.0, 6.0, 10.0, 10.0, 7.0, 7.0, 7.0, 10.0, 2.0, 8.0, 3.0, 7.0, 6.0, 10.0, 6.0, 4.0, 7.0, 4.0, 4.0, 3.0, 3.0, 8.0, 3.0, 2.0, 9.0, 4.0, 2.0, 10.0, 6.0, 1.0, 5.0, 2.0, 10.0, 3.0, 6.0, 6.0, 9.0, 4.0, 2.0, 8.0, 5.0, 3.0, 8.0, 9.0, 9.0, 5.0, 1.0, 6.0, 3.0, 6.0, 2.0, 7.0, 9.0, 1.0, 4.0, 6.0, 8.0, 7.0, 7.0, 6.0, 6.0, 10.0, 6.0, 9.0, 2.0, 3.0, 7.0, 4.0, 9.0, 10.0, 1.0, 2.0, 3.0, 7.0, 9.0, 9.0, 5.0, 6.0, 1.0, 3.0, 9.0, 9.0, 3.0, 5.0]
global b_x = 5
global d_y = [1.0, 4.0, 10.0, 8.0, 10.0, 4.0, 5.0, 7.0, 2.0, 2.0, 5.0, 6.0, 2.0, 8.0, 9.0, 8.0, 7.0, 7.0, 6.0, 5.0, 9.0, 7.0, 1.0, 9.0, 10.0, 6.0, 1.0, 1.0, 2.0, 8.0, 10.0, 6.0, 2.0, 4.0, 8.0, 1.0, 5.0, 3.0, 10.0, 8.0, 7.0, 1.0, 7.0, 4.0, 1.0, 2.0, 2.0, 4.0, 7.0, 4.0, 9.0, 3.0, 8.0, 2.0, 7.0, 5.0, 8.0, 1.0, 3.0, 10.0, 9.0, 2.0, 2.0, 10.0, 2.0, 9.0, 9.0, 4.0, 6.0, 8.0, 2.0, 4.0, 4.0, 5.0, 4.0, 4.0, 2.0, 4.0, 3.0, 4.0, 5.0, 8.0, 2.0, 9.0, 4.0, 5.0, 5.0, 5.0, 2.0, 7.0, 10.0, 4.0, 10.0, 4.0, 8.0, 3.0, 7.0, 3.0, 6.0, 1.0, 5.0, 7.0, 9.0, 7.0, 3.0, 9.0, 10.0, 6.0, 4.0, 5.0, 7.0, 5.0, 1.0, 8.0, 1.0, 4.0, 7.0, 8.0, 5.0, 2.0, 3.0, 2.0, 6.0, 5.0, 4.0, 6.0, 8.0, 2.0, 9.0, 5.0, 7.0, 5.0, 7.0, 8.0, 3.0, 2.0, 10.0, 9.0, 7.0, 2.0, 10.0, 8.0, 2.0, 5.0, 7.0, 10.0, 10.0, 9.0, 7.0, 2.0, 3.0, 10.0, 2.0, 8.0, 10.0, 2.0, 7.0, 3.0, 10.0, 9.0, 10.0, 1.0, 5.0, 6.0, 6.0, 4.0, 4.0, 9.0, 4.0, 2.0, 9.0, 2.0, 9.0, 7.0, 4.0, 7.0, 6.0, 1.0, 3.0, 3.0, 1.0, 4.0, 8.0, 1.0, 8.0, 10.0, 8.0, 2.0, 3.0, 8.0, 10.0, 1.0, 1.0, 5.0, 5.0, 3.0, 9.0, 2.0, 7.0, 7.0, 3.0, 8.0, 6.0, 4.0, 8.0, 5.0, 8.0, 9.0, 2.0, 8.0, 3.0, 4.0, 7.0, 7.0, 4.0, 3.0, 3.0, 3.0, 4.0, 9.0, 2.0, 1.0, 7.0, 9.0, 9.0, 5.0, 6.0, 6.0, 6.0, 5.0, 4.0, 6.0, 2.0]
global b_y = 10
global p = [0.639, 0.622, 0.535, 0.122, 0.443, 0.095, 0.79, 0.007, 0.865, 0.616, 0.006, 0.945, 0.237, 0.121, 0.885, 0.071, 0.476, 0.486, 0.475, 0.401, 0.38, 0.682, 0.897, 0.111, 0.543, 0.175, 0.562, 0.874, 0.215, 0.397, 0.174, 0.205, 0.699, 0.872, 0.913, 0.321, 0.679, 0.408, 0.186, 0.392, 0.394, 0.687, 0.608, 0.494, 0.108, 0.987, 0.7, 0.355, 0.447, 0.202, 0.043, 0.16, 0.965, 0.531, 0.969, 0.612, 0.488, 0.718, 0.945, 0.502, 0.294, 0.146, 0.91, 0.721, 0.326, 0.416, 0.957, 0.588, 0.544, 0.943, 0.081, 0.042, 0.204, 0.19, 0.887, 0.278, 0.353, 0.502, 0.33, 0.191, 0.761, 0.401, 0.781, 0.856, 0.666, 0.984, 0.72, 0.208, 0.192, 0.322, 0.003, 0.011, 0.548, 0.113, 0.111, 0.953, 0.266, 0.227, 0.78, 0.356, 0.294, 0.056, 0.89, 0.469, 0.715, 0.644, 0.177, 0.199, 0.603, 0.502, 0.018, 0.335, 0.946, 0.664, 0.477, 0.882, 0.41, 0.573, 0.997, 0.869, 0.026, 0.686, 0.285, 0.14, 0.226, 0.126, 0.206, 0.542, 0.619, 0.155, 0.634, 0.097, 0.02, 0.431, 0.31, 0.285, 0.941, 0.221, 0.781, 0.676, 0.785, 0.594, 0.597, 0.123, 0.135, 0.118, 0.724, 0.565, 0.81, 0.996, 0.128, 0.54, 0.228, 0.653, 0.136, 0.778, 0.543, 0.24, 0.33, 0.792, 0.593, 0.29, 0.201, 0.636, 0.777, 0.88, 0.298, 0.615, 0.302, 0.625, 0.565, 0.842, 0.119, 0.335, 0.372, 0.87, 0.805, 0.732, 0.597, 0.276, 0.912, 0.414, 0.749, 0.622, 0.248, 0.343, 0.981, 0.214, 0.4, 0.982, 0.017, 0.951, 0.805, 0.566, 0.549, 0.901, 0.352, 0.856, 0.987, 0.343, 0.871, 0.606, 0.805, 0.973, 0.658, 0.731, 0.735, 0.4, 0.126, 0.207, 0.125, 0.323, 0.841, 0.373, 0.361, 0.516, 0.585, 0.267, 0.103, 0.127, 0.825, 0.838, 0.59, 0.659, 0.759, 0.443, 0.718, 0.525, 0.719, 0.815, 0.216, 0.288, 0.032]
global q = [0.896, 0.676, 0.951, 0.431, 0.757, 0.866, 0.918, 0.795, 0.904, 0.641, 0.26, 0.973, 0.55, 0.402, 0.929, 0.285, 0.483, 0.97, 0.97, 0.714, 0.908, 0.951, 0.999, 0.801, 0.931, 0.345, 0.863, 0.917, 0.236, 0.633, 0.977, 0.215, 0.788, 0.874, 0.979, 0.502, 0.871, 0.637, 0.423, 0.561, 0.979, 0.935, 0.941, 0.621, 0.129, 0.991, 0.756, 0.495, 0.615, 0.704, 0.917, 0.545, 0.969, 0.982, 0.994, 0.706, 0.712, 0.721, 0.986, 0.669, 0.8, 0.392, 0.957, 0.721, 0.604, 0.938, 0.986, 0.725, 0.615, 0.948, 0.616, 0.922, 0.802, 0.702, 0.96, 0.434, 0.626, 0.911, 0.958, 0.855, 0.786, 0.534, 0.883, 0.969, 0.892, 0.989, 0.843, 0.456, 0.784, 0.444, 0.58, 0.879, 0.794, 0.137, 0.203, 0.992, 0.592, 0.613, 0.877, 0.834, 0.351, 0.164, 0.928, 0.754, 0.88, 0.652, 0.538, 0.774, 0.691, 0.64, 0.612, 0.768, 0.965, 0.717, 0.71, 0.901, 0.66, 0.836, 0.997, 0.913, 0.12, 0.951, 0.874, 0.546, 0.89, 0.867, 0.368, 0.939, 0.724, 0.934, 0.666, 0.821, 0.958, 0.44, 0.555, 0.905, 0.958, 0.543, 0.978, 0.949, 0.919, 0.613, 0.933, 0.607, 0.783, 0.802, 0.965, 0.914, 0.834, 0.996, 0.889, 0.747, 0.386, 0.816, 0.775, 0.895, 0.978, 0.446, 0.504, 0.913, 0.935, 0.297, 0.696, 0.651, 0.968, 0.956, 0.628, 0.919, 0.999, 0.664, 0.59, 0.982, 0.445, 0.469, 0.522, 0.919, 0.865, 0.769, 0.885, 0.564, 0.951, 0.548, 0.752, 0.666, 0.384, 0.487, 0.997, 0.718, 0.709, 0.998, 0.934, 0.957, 0.981, 0.884, 0.904, 0.919, 0.503, 0.904, 0.988, 0.431, 0.933, 0.696, 0.821, 0.974, 0.687, 0.781, 0.988, 0.956, 0.655, 0.696, 0.687, 0.708, 0.868, 0.418, 0.514, 0.533, 0.908, 0.634, 0.166, 0.642, 0.893, 0.913, 0.915, 0.702, 0.951, 0.911, 0.84, 0.788, 0.773, 0.966, 0.482, 0.475, 0.965]
global origin = 1
global destination = 50