global arcs = [1 13; 1 30; 2 36; 3 13; 3 20; 3 21; 3 33; 3 38; 4 9; 4 12; 4 21; 4 25; 4 37; 5 2; 5 16; 5 31; 6 10; 6 30; 7 5; 7 11; 7 19; 8 2; 8 20; 8 30; 9 25; 9 26; 10 2; 10 3; 10 15; 10 18; 10 35; 10 37; 10 40; 11 7; 11 17; 11 19; 11 29; 11 37; 12 6; 12 15; 12 30; 12 31; 12 32; 12 33; 12 37; 12 40; 13 7; 13 10; 13 18; 13 20; 13 22; 13 23; 13 25; 13 31; 14 16; 14 28; 15 5; 16 9; 16 30; 16 31; 17 8; 17 15; 17 23; 17 26; 17 28; 17 36; 18 2; 18 23; 18 32; 19 6; 19 18; 19 30; 19 34; 19 40; 20 6; 20 9; 20 27; 21 7; 21 11; 21 22; 21 36; 22 9; 22 16; 22 25; 22 32; 22 34; 22 38; 23 18; 23 19; 23 37; 23 40; 24 37; 24 38; 25 23; 25 33; 25 38; 26 10; 26 22; 26 25; 26 28; 26 29; 27 2; 27 6; 27 8; 27 9; 27 13; 27 14; 27 24; 27 29; 27 30; 27 37; 28 5; 28 9; 28 20; 28 22; 28 32; 28 37; 28 38; 29 21; 29 23; 29 40; 30 11; 30 12; 30 18; 30 20; 30 25; 30 36; 30 39; 31 2; 31 5; 31 7; 31 10; 31 17; 32 4; 32 35; 32 38; 32 40; 33 16; 33 22; 33 26; 33 35; 34 19; 34 31; 34 37; 35 12; 35 26; 35 33; 36 9; 36 20; 36 29; 37 5; 37 6; 37 16; 37 22; 37 29; 38 2; 38 6; 38 12; 38 18; 38 32; 38 36; 38 37; 39 10; 39 23; 39 27]
global d_x = [3.0, 3.0, 2.0, 8.0, 6.0, 2.0, 3.0, 3.0, 3.0, 2.0, 1.0, 4.0, 1.0, 7.0, 9.0, 9.0, 3.0, 3.0, 3.0, 8.0, 4.0, 8.0, 1.0, 1.0, 10.0, 3.0, 3.0, 10.0, 6.0, 3.0, 7.0, 10.0, 4.0, 1.0, 2.0, 2.0, 10.0, 6.0, 2.0, 9.0, 1.0, 3.0, 4.0, 2.0, 7.0, 8.0, 9.0, 6.0, 7.0, 2.0, 10.0, 10.0, 9.0, 5.0, 6.0, 8.0, 1.0, 8.0, 10.0, 9.0, 6.0, 7.0, 9.0, 1.0, 3.0, 5.0, 3.0, 2.0, 3.0, 2.0, 1.0, 6.0, 3.0, 1.0, 4.0, 9.0, 10.0, 4.0, 7.0, 3.0, 9.0, 9.0, 4.0, 5.0, 4.0, 2.0, 5.0, 6.0, 7.0, 6.0, 9.0, 4.0, 2.0, 1.0, 9.0, 6.0, 10.0, 5.0, 2.0, 10.0, 3.0, 8.0, 10.0, 10.0, 10.0, 9.0, 4.0, 1.0, 6.0, 1.0, 5.0, 6.0, 2.0, 5.0, 4.0, 9.0, 3.0, 9.0, 6.0, 3.0, 9.0, 1.0, 1.0, 8.0, 4.0, 7.0, 6.0, 2.0, 5.0, 9.0, 7.0, 10.0, 8.0, 5.0, 10.0, 4.0, 4.0, 8.0, 1.0, 5.0, 10.0, 5.0, 7.0, 8.0, 9.0, 5.0, 1.0, 8.0, 7.0, 3.0, 5.0, 8.0, 1.0, 6.0, 4.0, 3.0, 5.0, 8.0, 8.0, 2.0, 8.0, 2.0, 4.0, 3.0, 5.0]
global b_x = 5
global d_y = [6.0, 6.0, 6.0, 6.0, 10.0, 7.0, 10.0, 6.0, 2.0, 7.0, 10.0, 4.0, 3.0, 10.0, 4.0, 9.0, 9.0, 4.0, 4.0, 6.0, 7.0, 1.0, 10.0, 6.0, 10.0, 3.0, 10.0, 10.0, 8.0, 9.0, 1.0, 1.0, 9.0, 7.0, 8.0, 10.0, 10.0, 5.0, 1.0, 8.0, 6.0, 8.0, 1.0, 5.0, 1.0, 4.0, 3.0, 8.0, 9.0, 4.0, 5.0, 6.0, 8.0, 2.0, 9.0, 1.0, 9.0, 9.0, 8.0, 6.0, 2.0, 6.0, 1.0, 9.0, 1.0, 1.0, 3.0, 6.0, 8.0, 2.0, 6.0, 9.0, 8.0, 10.0, 6.0, 6.0, 6.0, 1.0, 6.0, 8.0, 1.0, 1.0, 10.0, 2.0, 2.0, 9.0, 4.0, 6.0, 5.0, 7.0, 9.0, 5.0, 2.0, 5.0, 10.0, 8.0, 3.0, 7.0, 4.0, 2.0, 9.0, 8.0, 9.0, 4.0, 2.0, 9.0, 8.0, 3.0, 10.0, 10.0, 3.0, 5.0, 6.0, 3.0, 3.0, 4.0, 9.0, 10.0, 6.0, 5.0, 5.0, 5.0, 4.0, 2.0, 2.0, 1.0, 1.0, 10.0, 2.0, 10.0, 10.0, 7.0, 3.0, 8.0, 6.0, 6.0, 4.0, 3.0, 7.0, 4.0, 10.0, 5.0, 4.0, 2.0, 9.0, 4.0, 6.0, 10.0, 4.0, 2.0, 1.0, 1.0, 4.0, 8.0, 3.0, 10.0, 7.0, 4.0, 6.0, 2.0, 6.0, 1.0, 7.0, 4.0, 9.0]
global b_y = 10
global p = [0.391, 0.136, 0.793, 0.022, 0.106, 0.834, 0.452, 0.849, 0.152, 0.048, 0.87, 0.486, 0.459, 0.78, 0.059, 0.452, 0.07, 0.387, 0.821, 0.829, 0.222, 0.573, 0.567, 0.667, 0.079, 0.763, 0.255, 0.771, 0.747, 0.75, 0.151, 0.242, 0.839, 0.465, 0.312, 0.868, 0.029, 0.201, 0.167, 0.133, 0.933, 0.652, 0.138, 0.269, 0.545, 0.211, 0.409, 0.809, 0.37, 0.334, 0.035, 0.648, 0.26, 0.069, 0.124, 0.626, 0.981, 0.429, 0.383, 0.282, 0.014, 0.394, 0.567, 0.749, 0.877, 0.227, 0.791, 0.415, 0.625, 0.159, 0.491, 0.268, 0.798, 0.165, 0.252, 0.597, 0.712, 0.908, 0.888, 0.936, 0.327, 0.155, 0.749, 0.243, 0.796, 0.351, 0.057, 0.528, 0.68, 0.029, 0.859, 0.457, 0.065, 0.134, 0.823, 0.906, 0.929, 0.639, 0.011, 0.561, 0.872, 0.225, 0.937, 0.459, 0.516, 0.712, 0.713, 0.507, 0.457, 0.523, 0.855, 0.84, 0.591, 0.251, 0.515, 0.162, 0.804, 0.952, 0.621, 0.525, 0.785, 0.885, 0.929, 0.857, 0.552, 0.701, 0.273, 0.98, 0.04, 0.418, 0.063, 0.083, 0.886, 0.831, 0.496, 0.663, 0.776, 0.227, 0.526, 0.83, 0.22, 0.473, 0.13, 0.222, 0.746, 0.066, 0.459, 0.165, 0.204, 0.76, 0.229, 0.435, 0.764, 0.302, 0.125, 0.405, 0.953, 0.168, 0.157, 0.927, 0.715, 0.863, 0.272, 0.339, 0.644]
global q = [0.518, 0.429, 0.822, 0.423, 0.186, 0.998, 0.999, 0.919, 0.235, 0.977, 0.998, 0.636, 0.931, 0.942, 0.795, 0.758, 0.715, 0.468, 0.845, 0.853, 0.336, 0.798, 0.682, 0.702, 0.7, 0.778, 0.341, 0.795, 0.897, 0.788, 0.986, 0.864, 0.882, 0.825, 0.774, 0.969, 0.31, 0.515, 0.249, 0.618, 0.984, 0.8, 0.729, 0.845, 0.741, 0.4, 0.585, 0.954, 0.816, 0.89, 0.944, 0.96, 0.369, 0.646, 0.579, 0.737, 0.981, 0.474, 0.706, 0.923, 0.769, 0.963, 0.812, 0.867, 0.982, 0.793, 0.866, 0.82, 0.678, 0.922, 0.876, 0.51, 0.999, 0.737, 0.635, 0.861, 0.938, 0.958, 0.947, 0.948, 0.746, 0.386, 0.798, 0.558, 0.99, 0.513, 0.286, 0.563, 0.935, 0.216, 0.977, 0.862, 0.921, 0.439, 0.905, 0.945, 0.947, 0.971, 0.465, 0.896, 0.907, 0.53, 0.937, 0.847, 0.719, 0.975, 0.856, 0.984, 0.981, 0.717, 0.908, 0.949, 0.769, 0.557, 0.696, 0.521, 0.975, 0.996, 0.686, 0.864, 0.925, 0.911, 0.982, 0.959, 0.803, 0.97, 0.704, 0.987, 0.589, 0.95, 0.082, 0.183, 0.991, 0.858, 0.976, 0.729, 0.889, 0.801, 0.752, 0.853, 0.776, 0.983, 0.393, 0.995, 0.856, 0.141, 0.976, 0.776, 0.742, 0.764, 0.255, 0.754, 0.891, 0.904, 0.28, 0.779, 0.97, 0.82, 0.76, 0.941, 0.893, 0.917, 0.866, 0.612, 0.787]
global origin = 1
global destination = 40