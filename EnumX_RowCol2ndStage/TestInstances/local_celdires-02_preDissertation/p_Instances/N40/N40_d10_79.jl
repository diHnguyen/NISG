global arcs = [1 39; 2 7; 2 20; 2 28; 3 27; 3 33; 3 36; 3 40; 4 15; 4 19; 4 24; 4 40; 5 17; 5 29; 5 40; 6 4; 6 7; 6 9; 6 14; 6 36; 7 11; 7 18; 7 30; 8 16; 8 17; 8 19; 8 25; 8 29; 8 30; 8 40; 9 5; 9 22; 9 35; 9 39; 10 6; 10 15; 10 28; 10 32; 11 10; 11 35; 11 39; 12 16; 12 21; 12 36; 12 37; 12 40; 13 15; 13 20; 13 25; 13 30; 14 3; 14 7; 14 17; 14 19; 14 30; 14 32; 14 40; 15 8; 15 10; 15 17; 15 35; 15 38; 16 5; 16 19; 16 30; 17 9; 17 33; 17 35; 18 4; 18 32; 19 4; 19 12; 19 22; 19 31; 20 30; 21 16; 21 35; 21 37; 21 38; 22 16; 22 23; 23 2; 23 25; 23 29; 23 32; 24 2; 24 26; 24 37; 25 3; 25 8; 25 12; 25 39; 26 6; 26 19; 26 30; 26 35; 27 19; 27 23; 27 28; 27 31; 28 35; 29 11; 29 32; 29 34; 29 39; 30 3; 30 13; 30 14; 30 32; 30 34; 30 35; 30 37; 31 13; 32 7; 32 14; 32 29; 33 4; 33 15; 33 17; 34 5; 34 8; 34 33; 35 11; 35 20; 35 22; 35 29; 35 34; 36 33; 37 4; 37 22; 37 30; 38 2; 38 3; 38 6; 38 12; 38 30; 39 4; 39 7; 39 12; 39 17; 39 37; 39 38]
global d_x = [6.0, 1.0, 3.0, 3.0, 7.0, 1.0, 3.0, 5.0, 9.0, 4.0, 7.0, 2.0, 7.0, 7.0, 4.0, 8.0, 8.0, 7.0, 10.0, 7.0, 7.0, 8.0, 7.0, 10.0, 5.0, 1.0, 6.0, 1.0, 5.0, 10.0, 6.0, 2.0, 7.0, 2.0, 8.0, 3.0, 9.0, 2.0, 8.0, 7.0, 4.0, 7.0, 9.0, 1.0, 8.0, 6.0, 8.0, 6.0, 6.0, 8.0, 4.0, 8.0, 9.0, 6.0, 9.0, 2.0, 10.0, 7.0, 2.0, 3.0, 4.0, 7.0, 1.0, 6.0, 6.0, 7.0, 5.0, 9.0, 4.0, 2.0, 3.0, 10.0, 5.0, 9.0, 7.0, 2.0, 2.0, 1.0, 6.0, 10.0, 8.0, 10.0, 1.0, 8.0, 10.0, 4.0, 7.0, 3.0, 7.0, 7.0, 9.0, 6.0, 8.0, 6.0, 4.0, 2.0, 3.0, 4.0, 10.0, 2.0, 6.0, 7.0, 2.0, 6.0, 2.0, 10.0, 4.0, 10.0, 3.0, 6.0, 9.0, 8.0, 1.0, 6.0, 5.0, 7.0, 6.0, 4.0, 4.0, 7.0, 3.0, 8.0, 2.0, 6.0, 3.0, 7.0, 6.0, 1.0, 7.0, 1.0, 2.0, 4.0, 2.0, 9.0, 3.0, 10.0, 10.0, 7.0, 4.0, 8.0, 10.0, 6.0]
global b_x = 5
global d_y = [4.0, 7.0, 9.0, 10.0, 5.0, 5.0, 3.0, 3.0, 10.0, 9.0, 4.0, 4.0, 4.0, 5.0, 1.0, 1.0, 3.0, 8.0, 7.0, 6.0, 6.0, 8.0, 4.0, 10.0, 10.0, 8.0, 2.0, 7.0, 6.0, 4.0, 5.0, 6.0, 3.0, 7.0, 5.0, 6.0, 10.0, 3.0, 6.0, 8.0, 5.0, 10.0, 6.0, 1.0, 8.0, 2.0, 8.0, 6.0, 10.0, 1.0, 9.0, 1.0, 3.0, 5.0, 3.0, 5.0, 9.0, 7.0, 9.0, 8.0, 2.0, 5.0, 3.0, 9.0, 6.0, 8.0, 2.0, 8.0, 8.0, 7.0, 10.0, 3.0, 3.0, 2.0, 9.0, 4.0, 1.0, 3.0, 4.0, 3.0, 3.0, 9.0, 1.0, 6.0, 7.0, 2.0, 6.0, 3.0, 1.0, 5.0, 10.0, 8.0, 1.0, 5.0, 6.0, 4.0, 8.0, 7.0, 1.0, 7.0, 5.0, 1.0, 10.0, 10.0, 7.0, 3.0, 7.0, 10.0, 7.0, 2.0, 4.0, 3.0, 4.0, 3.0, 6.0, 4.0, 10.0, 1.0, 6.0, 4.0, 6.0, 10.0, 5.0, 9.0, 1.0, 7.0, 2.0, 8.0, 7.0, 10.0, 5.0, 7.0, 7.0, 7.0, 2.0, 7.0, 9.0, 10.0, 2.0, 9.0, 6.0, 6.0]
global b_y = 10
global p = [0.664, 0.209, 0.104, 0.887, 0.496, 0.721, 0.585, 0.112, 0.813, 0.992, 0.072, 0.404, 0.946, 0.108, 0.957, 0.066, 0.364, 0.735, 0.084, 0.249, 0.126, 0.239, 0.12, 0.888, 0.388, 0.937, 0.889, 0.545, 0.288, 0.868, 0.563, 0.3, 0.132, 0.847, 0.571, 0.48, 0.581, 0.086, 0.306, 0.108, 0.014, 0.724, 0.87, 0.926, 0.078, 0.513, 0.481, 0.027, 0.695, 0.727, 0.03, 0.751, 0.859, 0.729, 0.328, 0.348, 0.222, 0.661, 0.308, 0.111, 0.323, 0.255, 0.287, 0.746, 0.171, 0.753, 0.506, 0.495, 0.07, 0.009, 0.776, 0.94, 0.456, 0.727, 0.912, 0.31, 0.034, 0.174, 0.747, 0.882, 0.901, 0.859, 0.051, 0.861, 0.956, 0.868, 0.273, 0.12, 0.94, 0.255, 0.791, 0.031, 0.684, 0.342, 0.501, 0.464, 0.101, 0.576, 0.638, 0.437, 0.456, 0.694, 0.279, 0.005, 0.758, 0.763, 0.05, 0.507, 0.396, 0.887, 0.192, 0.914, 0.041, 0.093, 0.846, 0.83, 0.571, 0.188, 0.373, 0.732, 0.6, 0.311, 0.745, 0.647, 0.771, 0.992, 0.14, 0.115, 0.873, 0.487, 0.767, 0.779, 0.858, 0.338, 0.761, 0.845, 0.151, 0.126, 0.702, 0.719, 0.977, 0.004]
global q = [0.751, 0.268, 0.613, 0.956, 0.997, 0.753, 0.625, 0.999, 0.939, 0.997, 0.807, 0.794, 0.99, 0.534, 0.99, 0.986, 0.865, 0.957, 0.297, 0.918, 0.207, 0.382, 0.39, 0.999, 0.416, 0.96, 0.957, 0.872, 0.731, 0.968, 0.902, 0.394, 0.987, 0.994, 0.85, 0.864, 0.959, 0.46, 0.925, 0.682, 0.941, 0.986, 0.979, 0.999, 0.58, 0.906, 0.573, 0.865, 0.92, 0.917, 0.122, 0.901, 0.951, 0.835, 0.822, 0.967, 0.739, 0.885, 0.875, 0.908, 0.487, 0.353, 0.694, 0.746, 0.837, 0.959, 0.684, 0.993, 0.696, 0.825, 0.975, 0.959, 0.625, 0.985, 0.973, 0.881, 0.91, 0.369, 0.896, 0.93, 0.914, 0.997, 0.132, 0.918, 0.957, 0.967, 0.805, 0.3, 0.98, 0.443, 0.831, 0.653, 0.705, 0.638, 0.615, 0.641, 0.794, 0.666, 0.927, 0.788, 0.763, 0.869, 0.751, 0.943, 0.782, 0.963, 0.459, 0.666, 0.955, 0.979, 0.71, 0.947, 0.182, 0.642, 0.85, 0.852, 0.747, 0.346, 0.673, 0.961, 0.879, 0.529, 0.812, 0.656, 0.824, 0.992, 0.906, 0.275, 0.967, 0.575, 0.997, 0.876, 0.925, 0.414, 0.922, 0.902, 0.522, 0.449, 0.918, 0.738, 0.98, 0.215]
global origin = 1
global destination = 40