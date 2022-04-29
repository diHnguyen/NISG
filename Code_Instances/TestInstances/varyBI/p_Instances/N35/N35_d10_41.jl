global arcs = [1 23; 1 25; 2 10; 2 15; 2 16; 2 18; 2 27; 3 2; 3 5; 3 8; 4 3; 4 8; 4 24; 5 9; 5 12; 5 26; 6 5; 6 13; 6 25; 7 19; 7 20; 7 21; 7 32; 7 34; 8 7; 8 9; 8 12; 8 18; 8 27; 8 33; 9 19; 9 31; 10 11; 10 16; 10 26; 10 31; 10 33; 10 34; 11 3; 12 5; 12 9; 12 10; 12 13; 13 8; 14 3; 14 13; 14 23; 14 30; 15 8; 15 28; 15 33; 16 5; 16 24; 16 26; 16 27; 16 30; 16 31; 16 35; 17 5; 17 7; 17 22; 17 25; 17 27; 17 33; 18 22; 18 24; 18 25; 18 28; 18 32; 18 34; 19 2; 19 11; 19 16; 19 23; 20 7; 20 13; 20 14; 20 24; 21 2; 21 9; 21 16; 21 30; 22 13; 22 20; 22 24; 22 28; 23 33; 23 35; 24 8; 24 16; 25 13; 25 28; 25 35; 26 6; 26 16; 26 21; 26 29; 27 9; 27 19; 27 21; 28 21; 28 24; 28 27; 29 11; 29 17; 29 22; 29 23; 30 7; 31 7; 31 9; 31 26; 31 33; 32 18; 32 21; 33 4; 33 14; 33 18; 33 21; 33 24; 34 3; 34 6; 34 8; 34 12; 34 29; 34 33]
global d_x = [1.0, 5.0, 9.0, 1.0, 1.0, 5.0, 8.0, 6.0, 1.0, 6.0, 7.0, 10.0, 9.0, 2.0, 6.0, 5.0, 4.0, 8.0, 10.0, 1.0, 5.0, 9.0, 4.0, 10.0, 3.0, 4.0, 3.0, 5.0, 7.0, 10.0, 3.0, 5.0, 7.0, 6.0, 6.0, 8.0, 10.0, 5.0, 4.0, 9.0, 1.0, 2.0, 1.0, 9.0, 10.0, 10.0, 2.0, 4.0, 10.0, 6.0, 10.0, 6.0, 1.0, 8.0, 3.0, 2.0, 8.0, 1.0, 7.0, 9.0, 9.0, 3.0, 2.0, 1.0, 6.0, 8.0, 8.0, 4.0, 2.0, 2.0, 8.0, 1.0, 5.0, 10.0, 9.0, 1.0, 2.0, 5.0, 5.0, 6.0, 3.0, 6.0, 8.0, 4.0, 7.0, 9.0, 5.0, 4.0, 8.0, 2.0, 3.0, 4.0, 6.0, 4.0, 7.0, 2.0, 6.0, 6.0, 4.0, 7.0, 1.0, 6.0, 8.0, 4.0, 9.0, 5.0, 9.0, 5.0, 6.0, 4.0, 6.0, 9.0, 4.0, 7.0, 2.0, 3.0, 10.0, 10.0, 10.0, 8.0, 7.0, 4.0, 2.0, 6.0, 8.0]
global b_x = 5
global d_y = [9.0, 9.0, 4.0, 8.0, 4.0, 5.0, 9.0, 4.0, 7.0, 2.0, 4.0, 9.0, 9.0, 7.0, 8.0, 1.0, 10.0, 1.0, 9.0, 2.0, 7.0, 5.0, 3.0, 10.0, 8.0, 6.0, 8.0, 5.0, 2.0, 9.0, 3.0, 10.0, 10.0, 2.0, 5.0, 3.0, 10.0, 6.0, 8.0, 8.0, 6.0, 1.0, 7.0, 8.0, 4.0, 7.0, 7.0, 5.0, 8.0, 3.0, 9.0, 1.0, 4.0, 7.0, 3.0, 2.0, 6.0, 9.0, 10.0, 7.0, 2.0, 3.0, 1.0, 8.0, 4.0, 7.0, 2.0, 1.0, 8.0, 1.0, 1.0, 10.0, 6.0, 2.0, 7.0, 8.0, 7.0, 8.0, 9.0, 3.0, 1.0, 5.0, 9.0, 3.0, 2.0, 5.0, 7.0, 10.0, 8.0, 9.0, 1.0, 9.0, 9.0, 3.0, 8.0, 7.0, 6.0, 3.0, 1.0, 10.0, 2.0, 9.0, 3.0, 8.0, 9.0, 9.0, 5.0, 3.0, 6.0, 5.0, 4.0, 6.0, 4.0, 9.0, 2.0, 4.0, 6.0, 7.0, 1.0, 5.0, 4.0, 2.0, 8.0, 2.0, 1.0]
global b_y = 10
global p = [0.18, 0.722, 0.87, 0.565, 0.15, 0.149, 0.013, 0.797, 0.086, 0.3, 0.547, 0.187, 0.535, 0.81, 0.024, 0.308, 0.119, 0.467, 0.763, 0.83, 0.416, 0.679, 0.928, 0.062, 0.396, 0.12, 0.604, 0.488, 0.393, 0.766, 0.874, 0.627, 0.428, 0.107, 0.166, 0.76, 0.415, 0.716, 0.27, 0.768, 0.111, 0.849, 0.153, 0.086, 0.41, 0.999, 0.055, 0.265, 0.669, 0.894, 0.914, 0.246, 0.012, 0.03, 0.985, 0.677, 0.954, 0.588, 0.401, 0.334, 0.54, 0.835, 0.872, 0.897, 0.042, 0.664, 0.202, 0.3, 0.677, 0.559, 0.74, 0.699, 0.165, 0.68, 0.341, 0.792, 0.158, 0.979, 0.106, 0.927, 0.365, 0.117, 0.961, 0.116, 0.164, 0.025, 0.009, 0.923, 0.142, 0.024, 0.417, 0.017, 0.561, 0.057, 0.488, 0.093, 0.451, 0.153, 0.812, 0.463, 0.855, 0.277, 0.364, 0.751, 0.888, 0.127, 0.574, 0.311, 0.728, 0.283, 0.46, 0.032, 0.688, 0.88, 0.285, 0.887, 0.494, 0.829, 0.987, 0.842, 0.321, 0.097, 0.75, 0.245, 0.736]
global q = [0.551, 0.869, 0.893, 0.937, 0.307, 0.384, 0.482, 0.933, 0.337, 0.482, 0.808, 0.308, 0.763, 0.968, 0.541, 0.378, 0.51, 0.97, 0.775, 0.922, 0.513, 0.846, 0.99, 0.263, 0.627, 0.608, 0.672, 0.559, 0.683, 0.976, 0.958, 0.897, 0.938, 0.835, 0.169, 0.924, 0.749, 0.826, 0.279, 0.918, 0.685, 0.905, 0.221, 0.803, 0.632, 0.999, 0.51, 0.936, 0.89, 0.94, 0.919, 0.289, 0.056, 0.33, 0.987, 0.891, 0.961, 0.662, 0.88, 0.436, 0.628, 0.869, 0.96, 0.93, 0.107, 0.82, 0.654, 0.825, 0.913, 0.742, 0.967, 0.889, 0.623, 0.878, 0.773, 0.905, 0.339, 0.995, 0.353, 0.938, 0.999, 0.745, 0.99, 0.856, 0.813, 0.697, 0.299, 0.983, 0.804, 0.501, 0.879, 0.763, 0.72, 0.326, 0.797, 0.295, 0.864, 0.876, 0.834, 0.87, 0.93, 0.891, 0.659, 0.818, 0.93, 0.839, 0.589, 0.823, 0.951, 0.941, 0.724, 0.92, 0.85, 0.996, 0.38, 0.924, 0.507, 0.96, 0.999, 0.939, 0.807, 0.847, 0.875, 0.261, 0.856]
global origin = 1
global destination = 35