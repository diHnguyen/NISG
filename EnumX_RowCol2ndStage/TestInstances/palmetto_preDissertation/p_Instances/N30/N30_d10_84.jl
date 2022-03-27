global arcs = [1 11; 1 13; 1 16; 2 3; 2 12; 2 19; 2 22; 3 2; 4 3; 4 8; 4 10; 4 14; 5 2; 5 10; 5 15; 5 19; 5 27; 6 10; 6 16; 6 17; 6 24; 6 25; 6 27; 7 5; 7 28; 8 4; 8 26; 9 7; 9 18; 9 22; 9 23; 10 4; 10 12; 10 15; 10 19; 11 17; 11 23; 11 30; 12 4; 12 9; 12 10; 13 18; 14 11; 14 12; 14 15; 14 24; 14 30; 15 18; 16 2; 16 13; 16 21; 16 29; 17 3; 17 19; 17 21; 18 10; 19 29; 20 17; 20 22; 20 25; 20 29; 21 8; 21 20; 21 23; 22 12; 22 24; 23 14; 24 8; 24 11; 24 23; 25 2; 25 14; 26 4; 26 19; 26 23; 26 29; 27 4; 27 6; 27 20; 27 21; 28 18; 29 9; 29 14; 29 23; 29 27]
global d_x = [10.0, 9.0, 2.0, 1.0, 3.0, 2.0, 3.0, 10.0, 10.0, 8.0, 5.0, 3.0, 3.0, 5.0, 6.0, 7.0, 10.0, 10.0, 2.0, 9.0, 1.0, 1.0, 8.0, 3.0, 7.0, 6.0, 9.0, 2.0, 5.0, 1.0, 1.0, 5.0, 2.0, 3.0, 5.0, 8.0, 5.0, 8.0, 5.0, 8.0, 2.0, 4.0, 9.0, 1.0, 7.0, 5.0, 6.0, 9.0, 4.0, 5.0, 10.0, 5.0, 6.0, 8.0, 6.0, 5.0, 8.0, 3.0, 5.0, 7.0, 4.0, 1.0, 10.0, 10.0, 5.0, 2.0, 10.0, 4.0, 10.0, 2.0, 10.0, 7.0, 5.0, 2.0, 8.0, 2.0, 2.0, 2.0, 1.0, 2.0, 3.0, 5.0, 10.0, 3.0, 9.0]
global b_x = 5
global d_y = [7.0, 6.0, 7.0, 6.0, 6.0, 5.0, 6.0, 5.0, 10.0, 5.0, 3.0, 2.0, 10.0, 6.0, 2.0, 7.0, 7.0, 3.0, 6.0, 6.0, 6.0, 9.0, 5.0, 5.0, 7.0, 9.0, 1.0, 1.0, 7.0, 5.0, 9.0, 8.0, 3.0, 3.0, 4.0, 9.0, 4.0, 2.0, 3.0, 7.0, 2.0, 4.0, 6.0, 9.0, 1.0, 3.0, 1.0, 1.0, 1.0, 10.0, 6.0, 8.0, 10.0, 6.0, 6.0, 6.0, 10.0, 3.0, 10.0, 10.0, 10.0, 8.0, 1.0, 4.0, 2.0, 6.0, 8.0, 2.0, 1.0, 8.0, 9.0, 3.0, 1.0, 5.0, 10.0, 9.0, 1.0, 4.0, 2.0, 2.0, 4.0, 3.0, 2.0, 6.0, 2.0]
global b_y = 10
global p = [0.987, 0.434, 0.42, 0.695, 0.893, 0.343, 0.465, 0.745, 0.352, 0.959, 0.748, 0.365, 0.747, 0.671, 0.879, 0.313, 0.769, 0.138, 0.914, 0.318, 0.356, 0.435, 0.979, 0.1, 0.02, 0.321, 0.126, 0.551, 0.947, 0.953, 0.693, 0.252, 0.791, 0.413, 0.865, 0.571, 0.676, 0.108, 0.545, 0.002, 0.834, 0.241, 0.802, 0.253, 0.887, 0.532, 0.778, 0.925, 0.785, 0.148, 0.645, 0.248, 0.573, 0.518, 0.661, 0.574, 0.287, 0.912, 0.915, 0.954, 0.7, 0.92, 0.747, 0.847, 0.574, 0.162, 0.519, 0.617, 0.934, 0.311, 0.094, 0.695, 0.416, 0.844, 0.197, 0.33, 0.767, 0.289, 0.082, 0.088, 0.394, 0.425, 0.606, 0.886, 0.223]
global q = [0.991, 0.675, 0.524, 0.742, 0.994, 0.553, 0.967, 0.98, 0.889, 0.997, 0.777, 0.749, 0.885, 0.827, 0.945, 0.52, 0.96, 0.514, 0.963, 0.382, 0.778, 0.771, 0.991, 0.765, 0.416, 0.334, 0.957, 0.723, 0.985, 0.961, 0.762, 0.608, 0.924, 0.976, 0.993, 0.791, 0.729, 0.51, 0.848, 0.504, 0.908, 0.998, 0.907, 0.911, 0.942, 0.774, 0.846, 0.973, 0.831, 0.756, 0.689, 0.482, 0.759, 0.857, 0.932, 0.948, 0.3, 0.964, 0.952, 0.967, 0.902, 0.992, 0.874, 0.984, 0.784, 0.336, 0.852, 0.628, 0.992, 0.325, 0.484, 0.865, 0.77, 0.97, 0.584, 0.91, 0.776, 0.647, 0.613, 0.429, 0.629, 0.557, 0.957, 0.996, 0.328]
global origin = 1
global destination = 30