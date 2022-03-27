global arcs = [1 15; 1 16; 1 24; 2 10; 2 11; 2 15; 2 17; 2 25; 2 36; 3 6; 3 7; 3 20; 4 12; 4 21; 4 23; 5 7; 5 18; 5 22; 5 23; 5 25; 6 10; 6 11; 6 13; 6 20; 6 24; 6 33; 7 21; 7 25; 8 3; 8 9; 8 19; 8 20; 8 21; 8 23; 8 26; 8 30; 8 33; 9 2; 9 14; 9 15; 9 18; 9 19; 9 22; 10 27; 11 8; 11 14; 11 19; 11 28; 11 29; 12 34; 12 36; 12 37; 13 19; 13 26; 14 8; 14 15; 14 37; 15 18; 15 23; 16 2; 16 13; 16 30; 16 36; 17 23; 18 7; 18 19; 19 7; 19 22; 19 27; 20 8; 20 29; 20 35; 21 22; 21 30; 21 32; 22 8; 22 12; 22 21; 22 28; 22 32; 22 38; 23 8; 23 18; 23 19; 24 3; 24 30; 25 4; 25 19; 25 34; 26 14; 26 17; 27 3; 27 18; 27 23; 28 10; 28 21; 28 22; 28 40; 29 2; 29 4; 29 5; 29 10; 29 12; 29 26; 29 36; 30 9; 30 16; 30 29; 31 10; 31 16; 31 19; 31 21; 31 32; 32 17; 32 18; 32 25; 33 2; 33 13; 33 18; 34 17; 34 21; 34 23; 35 8; 35 30; 35 34; 35 39; 36 15; 37 5; 37 14; 37 25; 37 31; 38 4; 38 11; 39 7; 39 11; 39 16]
global d_x = [9.0, 9.0, 2.0, 10.0, 3.0, 7.0, 6.0, 9.0, 2.0, 5.0, 7.0, 9.0, 10.0, 10.0, 4.0, 4.0, 2.0, 6.0, 1.0, 2.0, 3.0, 10.0, 4.0, 8.0, 8.0, 2.0, 9.0, 3.0, 9.0, 1.0, 2.0, 5.0, 9.0, 2.0, 5.0, 7.0, 9.0, 4.0, 7.0, 9.0, 10.0, 1.0, 8.0, 8.0, 5.0, 1.0, 2.0, 4.0, 5.0, 10.0, 3.0, 5.0, 3.0, 3.0, 9.0, 4.0, 7.0, 5.0, 7.0, 10.0, 9.0, 9.0, 2.0, 2.0, 2.0, 7.0, 3.0, 3.0, 10.0, 4.0, 1.0, 8.0, 4.0, 3.0, 2.0, 6.0, 2.0, 8.0, 7.0, 4.0, 7.0, 10.0, 10.0, 10.0, 2.0, 6.0, 9.0, 2.0, 8.0, 1.0, 9.0, 3.0, 7.0, 2.0, 9.0, 2.0, 3.0, 9.0, 4.0, 9.0, 9.0, 6.0, 6.0, 9.0, 6.0, 5.0, 9.0, 6.0, 4.0, 5.0, 7.0, 8.0, 6.0, 1.0, 10.0, 5.0, 4.0, 7.0, 9.0, 2.0, 8.0, 9.0, 2.0, 8.0, 10.0, 6.0, 7.0, 8.0, 2.0, 1.0, 1.0, 10.0, 4.0, 10.0, 1.0, 3.0]
global b_x = 5
global d_y = [1.0, 2.0, 9.0, 5.0, 1.0, 3.0, 3.0, 5.0, 3.0, 6.0, 2.0, 10.0, 7.0, 6.0, 1.0, 4.0, 7.0, 5.0, 5.0, 1.0, 5.0, 10.0, 2.0, 2.0, 4.0, 8.0, 3.0, 5.0, 10.0, 7.0, 4.0, 7.0, 9.0, 5.0, 10.0, 7.0, 1.0, 6.0, 6.0, 6.0, 10.0, 2.0, 6.0, 7.0, 4.0, 9.0, 4.0, 8.0, 3.0, 4.0, 4.0, 10.0, 10.0, 9.0, 6.0, 2.0, 3.0, 4.0, 6.0, 7.0, 1.0, 2.0, 4.0, 8.0, 10.0, 6.0, 9.0, 4.0, 8.0, 9.0, 6.0, 2.0, 1.0, 2.0, 10.0, 3.0, 9.0, 2.0, 7.0, 4.0, 8.0, 7.0, 9.0, 9.0, 6.0, 10.0, 1.0, 3.0, 10.0, 9.0, 9.0, 10.0, 10.0, 7.0, 8.0, 4.0, 6.0, 2.0, 10.0, 7.0, 7.0, 2.0, 4.0, 5.0, 4.0, 2.0, 4.0, 3.0, 3.0, 1.0, 5.0, 9.0, 2.0, 9.0, 1.0, 8.0, 6.0, 6.0, 3.0, 5.0, 2.0, 8.0, 3.0, 3.0, 1.0, 8.0, 6.0, 6.0, 6.0, 3.0, 4.0, 7.0, 1.0, 6.0, 2.0, 3.0]
global b_y = 10
global p = [0.235, 0.417, 0.868, 0.885, 0.15, 0.732, 0.495, 0.51, 0.702, 0.249, 0.193, 0.679, 0.911, 0.79, 0.768, 0.273, 0.095, 0.707, 0.243, 0.123, 0.032, 0.263, 0.041, 0.777, 0.58, 0.64, 0.91, 0.319, 0.559, 0.554, 0.809, 0.4, 0.559, 0.775, 0.131, 0.489, 0.908, 0.929, 0.184, 0.164, 0.935, 0.689, 0.816, 0.456, 0.498, 0.3, 0.611, 0.294, 0.076, 0.194, 0.396, 0.125, 0.078, 0.766, 0.708, 0.658, 0.875, 0.622, 0.96, 0.315, 0.585, 0.1, 0.109, 0.928, 0.775, 0.155, 0.219, 0.828, 0.001, 0.499, 0.242, 0.409, 0.281, 0.882, 0.711, 0.134, 0.247, 0.502, 0.991, 0.043, 0.75, 0.479, 0.263, 0.044, 0.696, 0.533, 0.593, 0.382, 0.195, 0.162, 0.33, 0.599, 0.442, 0.393, 0.283, 0.146, 0.896, 0.966, 0.929, 0.979, 0.854, 0.735, 0.197, 0.985, 0.138, 0.167, 0.917, 0.808, 0.261, 0.958, 0.931, 0.913, 0.81, 0.175, 0.963, 0.7, 0.91, 0.575, 0.772, 0.475, 0.559, 0.489, 0.457, 0.855, 0.665, 0.418, 0.892, 0.993, 0.522, 0.82, 0.789, 0.045, 0.784, 0.35, 0.106, 0.587]
global q = [0.445, 0.859, 0.904, 0.944, 0.827, 0.974, 0.715, 0.886, 0.813, 0.523, 0.458, 0.81, 0.956, 0.838, 0.899, 0.312, 0.916, 0.955, 0.533, 0.305, 0.488, 0.998, 0.621, 0.988, 0.955, 0.664, 0.948, 0.667, 0.994, 0.671, 0.828, 0.565, 0.937, 0.942, 0.811, 0.542, 0.988, 0.95, 0.34, 0.852, 0.973, 0.792, 0.925, 0.663, 0.799, 0.346, 0.998, 0.987, 0.125, 0.932, 0.522, 0.253, 0.596, 0.911, 0.75, 0.995, 0.93, 0.644, 0.963, 0.658, 0.685, 0.136, 0.295, 0.96, 0.813, 0.257, 0.284, 0.965, 0.38, 0.874, 0.914, 0.805, 0.838, 0.979, 0.903, 0.331, 0.935, 0.639, 0.999, 0.122, 0.855, 0.561, 0.924, 0.162, 0.905, 0.941, 0.652, 0.651, 0.523, 0.536, 0.521, 0.764, 0.512, 0.845, 0.856, 0.596, 0.958, 0.986, 0.991, 0.994, 0.968, 0.959, 0.952, 0.996, 0.535, 0.481, 0.955, 0.937, 0.334, 0.959, 0.947, 0.945, 0.878, 0.32, 0.993, 0.775, 0.934, 0.847, 0.832, 0.74, 0.88, 0.752, 0.568, 0.996, 0.815, 0.971, 0.997, 0.993, 0.839, 0.988, 0.843, 0.061, 0.893, 0.902, 0.776, 0.964]
global origin = 1
global destination = 40